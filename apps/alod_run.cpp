#include "alod/problems.hpp"
#include "alod/mesh_state.hpp"
#include "helmholtz/boundary.h"
#include "helmholtz/operators.h"
#ifdef ALOD_LEGACY_ORACLE
#include "helmholtz/model.h"
#include "helmholtz/adaptive/estimator.h"
#else
#include "alod/afem.hpp"
#include "alod/slod.hpp"
#endif
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <omp.h>
#ifndef ALOD_LEGACY_ORACLE
#include "adaptive_run.hpp"
#endif

int main(int argc,char** argv) {
#ifndef ALOD_LEGACY_ORACLE
 if(argc>1&&std::string(argv[1])=="adaptive")return adaptive_main(argc-1,argv+1);
#endif
 try {
    if(argc==2 && std::string(argv[1])=="--help") {
        std::cout<<"Usage: alod_run E1|E2 AFEM|UFEM|SLOD [--initial-level=N] [--states=N] [--theta=X] [--target=X] [--maximum-nodes=N] [--threads=N] [--emit-solution=0|1]\n";
        return 0;
    }
    if(argc<3)throw std::invalid_argument("expected problem and baseline method; see --help");
    const std::string id=argv[1],method=argv[2];
    if(method!="AFEM" && method!="UFEM" && method!="SLOD")throw std::invalid_argument("method must be AFEM, UFEM or SLOD");
    std::map<std::string,std::string> opts{{"initial-level","2"},{"states","3"},{"theta","0.15"},{"target","0"},{"maximum-nodes","20000"},{"threads","1"},{"emit-solution","0"}};
    std::set<std::string> seen;
    for(int i=3;i<argc;++i){std::string arg=argv[i];auto equal=arg.find('=');
        if(!arg.starts_with("--")||equal==std::string::npos)throw std::invalid_argument("options require --name=value");
        auto key=arg.substr(2,equal-2);
        if(!opts.contains(key)||!seen.insert(key).second)throw std::invalid_argument("unknown or duplicate option");
        opts[key]=arg.substr(equal+1);
    }
    auto integer=[&](const std::string& key){std::size_t n=0;int v=std::stoi(opts.at(key),&n);if(n!=opts.at(key).size())throw std::invalid_argument("invalid integer");return v;};
    auto real=[&](const std::string& key){std::size_t n=0;double v=std::stod(opts.at(key),&n);if(n!=opts.at(key).size()||!std::isfinite(v))throw std::invalid_argument("invalid finite number");return v;};
    int level=integer("initial-level"),states=integer("states"),cap=integer("maximum-nodes"),threads=integer("threads"),emit=integer("emit-solution");
    double theta=real("theta"),target=real("target");
    if(level<0||level>24||states<1||states>10000||cap<4||threads<1||threads>64||(emit!=0&&emit!=1)||theta<=0||theta>1||target<0)
        throw std::invalid_argument("baseline option outside supported range");
    omp_set_num_threads(threads);
    auto problem=alod::make_problem(id);auto quad=alod::paper_quadrature(id);
    using namespace lod2d;using namespace lod2d::helmholtz;
    // Bound the initial uniform allocation before constructing prolongations.
    if(std::ldexp(static_cast<double>(problem.initial_mesh.elems.size()),level)>2.0*cap)
        throw std::runtime_error("initial mesh exceeds maximum-nodes budget");
    alod::MeshState state(problem.initial_mesh);
    for(int i=0;i<level;++i){std::vector<int> marks(state.mesh.elems.size());std::iota(marks.begin(),marks.end(),0);state.refine(marks);}
    std::cout<<std::setprecision(17);
    for(int step=0;step<states;++step) {
        const auto& mesh=state.mesh;validate_boundary_tags(mesh);
        if(mesh.nodes.size()>static_cast<std::size_t>(cap))throw std::runtime_error("mesh exceeds maximum-nodes budget");
        TriMesh evaluation_mesh;ComplexVector values;
        double residual=0,patch_residual=0,constraint_residual=0;
        std::vector<double> indicators;std::vector<int> marks;
        if(method=="SLOD") {
            if(mesh.elems.size()*8+mesh.nodes.size()*4>static_cast<std::size_t>(cap))
                throw std::runtime_error("SLOD reference mesh exceeds maximum-nodes budget");
#ifdef ALOD_LEGACY_ORACLE
            HelmholtzProblemConfig cfg;cfg.H=level+step;cfg.h=cfg.H+4;cfg.ell=3;cfg.wavenumber=16;
            cfg.initial_mesh=problem.initial_mesh;cfg.quadrature=quad;cfg.quadrature_context=problem.quadrature_context;
            cfg.patch_solver.kind=HelmholtzPatchSolverKind::DirectSchur;
            cfg.patch_solver.symbolic_cache_slots=8;
            auto model=HelmholtzLodModel::build(cfg);
            auto sol=model.solve_source(problem.source);evaluation_mesh=model.problem().fine;values=sol.fine_values;
            residual=sol.petrov_residual;patch_residual=std::max(model.correctors().diagnostics.max_primal_residual,model.correctors().diagnostics.max_adjoint_residual);
            constraint_residual=model.correctors().diagnostics.max_constraint_residual;
#else
            auto sol=alod::solve_slod(mesh,problem,quad);
            evaluation_mesh=std::move(sol.fine);values=std::move(sol.values);
            residual=sol.residual;patch_residual=sol.patch_residual;constraint_residual=sol.constraint_residual;
#endif
        } else {
            evaluation_mesh=mesh;
            auto op=assemble_helmholtz_operators(mesh,16);
            auto load=assemble_helmholtz_load(mesh,problem.source,quad,problem.quadrature_context);
            values=solve_helmholtz_fem(op,load,HelmholtzFemSolverKind::Umfpack);
            ComplexVector r=op.system*values-load,b=load;
            for(int i:op.dirichlet_nodes){r[i]=0;b[i]=0;}
            residual=r.norm()/std::max(1e-12,b.norm());
            if(method=="AFEM") {
                auto estimate=adaptive::diagnostics::estimate_conforming_p1_residual(mesh,op,values,load,problem.source,quad,problem.quadrature_context);
                if(estimate.algebraic_relative_difference>1e-8)throw std::runtime_error("strong residual reconstruction failed");
                indicators=estimate.element_squared;
                marks=adaptive::mark_doerfler(indicators,theta);
            }
        }
        if(method!="AFEM"){marks.resize(mesh.elems.size());std::iota(marks.begin(),marks.end(),0);}
        auto err=compute_helmholtz_error(evaluation_mesh,values,16,problem.exact,problem.exact_gradient,quad,problem.quadrature_context);
        auto norm=compute_helmholtz_error(evaluation_mesh,ComplexVector::Zero(values.size()),16,problem.exact,problem.exact_gradient,quad,problem.quadrature_context);
        double relative=err.energy/norm.energy;
        if(!values.allFinite()||!std::isfinite(relative)||!std::isfinite(residual)||residual>1e-9)
            throw std::runtime_error("baseline finite-value/residual gate failed");
        bool reached=target>0 && relative<=target;
        std::cout<<"{\"problem\":\""<<id<<"\",\"method\":\""<<method<<"\",\"state\":"<<step
            <<",\"nodes\":"<<mesh.nodes.size()<<",\"elements\":"<<mesh.elems.size()
            <<",\"free_dof\":"<<mesh.nodes.size()-dirichlet_nodes(mesh).size()
            <<",\"reference_nodes\":"<<evaluation_mesh.nodes.size()<<",\"reference_elements\":"<<evaluation_mesh.elems.size()
            <<",\"mesh_fingerprint\":\""<<alod::mesh_fingerprint(mesh)<<"\",\"reference_fingerprint\":\""<<alod::mesh_fingerprint(evaluation_mesh)
            <<"\",\"energy_error\":"<<err.energy<<",\"l2_error\":"<<err.l2<<",\"exact_energy_norm\":"<<norm.energy
            <<",\"relative_energy_error\":"<<relative<<",\"residual\":"<<residual<<",\"patch_residual\":"<<patch_residual
            <<",\"constraint_residual\":"<<constraint_residual<<",\"target_reached\":"<<(reached?"true":"false")
            <<",\"stop_reason\":\""<<(reached?"exact_error_target":step+1==states?"state_limit":"continue")<<"\",\"marked_elements\":[";
        for(std::size_t i=0;i<marks.size();++i){if(i)std::cout<<',';std::cout<<marks[i];}
        std::cout<<"],\"marked_ids\":[";
        for(std::size_t i=0;i<marks.size();++i){if(i)std::cout<<',';std::cout<<state.elements[marks[i]].id;}
        std::cout<<"],\"element_ids\":[";
        for(std::size_t i=0;i<state.elements.size();++i){if(i)std::cout<<',';std::cout<<state.elements[i].id;}
        std::cout<<"],\"indicators_squared\":[";
        for(std::size_t i=0;i<indicators.size();++i){if(i)std::cout<<',';std::cout<<indicators[i];}
        std::cout<<"]";
        if(emit){std::cout<<",\"solution\":[";for(int i=0;i<values.size();++i){if(i)std::cout<<',';std::cout<<'['<<values[i].real()<<','<<values[i].imag()<<']';}std::cout<<']';}
        std::cout<<"}\n"<<std::flush;
        if(reached||step+1==states)break;
        state.refine(marks);
    }
 } catch(const std::exception& e) {std::cerr<<"alod_run: "<<e.what()<<'\n';return 1;}
}
