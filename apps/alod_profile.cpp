#include "fixed_support.hpp"
#include "alod/batch.hpp"
#include "alod/regional.hpp"
#include <chrono>
using namespace alod;
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point start){return std::chrono::duration<double>(Clock::now()-start).count();}
int main(int argc,char** argv){try{
    auto input=fixed::parse(argc,argv);auto H=fixed::coarse_mesh(input);auto h=lod2d::refine_mesh_nvb(H,input.gap);
    LodLimits limits;limits.threads=input.threads;limits.maximum_reference_nodes=input.cap;
    LodSpace space(H,h,16,input.ell,input.policy=="area"?InterpolationPolicy::ManuscriptAreaWeighted:InterpolationPolicy::ArchivedArithmetic,limits);
    AdditiveKernelRieszContext riesz(space);AdjointTestCache training_cache;
    std::vector<Problem> problems;for(auto& m:input.members)problems.push_back(m.problem);
    auto q=paper_quadrature(input.problem);auto loads=assemble_load_batch(h.mesh,problems,q,input.threads);
    auto trained=train_regional(space,riesz,training_cache,loads,{.9,6,false});
    if(trained.phi.cols()==0)throw std::runtime_error("profile requires a nonempty dictionary");
    auto mask=riesz.regional_mask(.9);double maximum_difference=0;
    auto compare=[&](const ComplexMatrix& a,const ComplexMatrix& b){double d=(a-b).norm()/std::max(1.,a.norm());maximum_difference=std::max(maximum_difference,d);if(d>1e-9)throw std::runtime_error("profile numerical equivalence failed");};
    std::map<std::string,std::vector<double>> times;
    std::map<std::string,std::vector<std::size_t>> columns;
    std::map<std::string,std::size_t> dense_bytes;
    auto measure=[&](const std::string& key,auto f){auto t=Clock::now();f();times[key].push_back(elapsed(t));};
    ReferenceFemContext reference(space.operators());auto ref=reference.solve(loads);
    for(int repeat=0;repeat<3;++repeat){
        for(bool reuse:{false,true}){
            AdjointTestCache aot;RegionalEvaluator evaluator(space,riesz,aot,reuse);
            std::string mode=reuse?"cached":"uncached";
            for(auto temperature:{"cold","hot"})measure(mode+"_"+temperature,[&]{auto result=evaluator.evaluate(loads,trained.phi,mask);compare(trained.accepted.values,result.values);});
            columns[mode].push_back(aot.solved_columns());
            dense_bytes[mode+"_persistent_cache"]=evaluator.dense_cache_bytes();
        }
        ErrorBatch baseline,fused;
        measure("error_three_pass",[&]{baseline=integrate_audit_batch(h.mesh,space.energy(),trained.accepted.values,ref,problems,q,input.threads,false);});
        measure("error_fused",[&]{fused=integrate_audit_batch(h.mesh,space.energy(),trained.accepted.values,ref,problems,q,input.threads,true);});
        for(auto pair:{std::pair{baseline.exact_error,fused.exact_error},std::pair{baseline.exact_norm,fused.exact_norm},std::pair{baseline.reference_error,fused.reference_error}}){
            double d=(pair.first-pair.second).norm()/std::max(1.,pair.first.norm());maximum_difference=std::max(maximum_difference,d);if(d>1e-12)throw std::runtime_error("fused integration changed numerical results");}
        // Candidate only: compare raw triangular and actual projected PG blocks.
        measure("raw_triangular_candidate",[&]{
            const auto& raw=trained.accepted.raw_tests;const auto& phi=trained.phi;
            ComplexMatrix lower=raw.adjoint()*space.operators().system*space.trial();
            if(lower.norm()/std::max(1.,space.trial().norm())>1e-8)throw std::runtime_error("raw triangular gate");
            ComplexMatrix d=(raw.adjoint()*space.operators().system*phi).fullPivLu().solve(raw.adjoint()*loads);
            auto adjusted=space.solve(loads-space.operators().system*phi*d).values;
            compare(trained.accepted.values,(adjusted+phi*d).eval());
        });
        ComplexMatrix rhs=space.energy().cast<Complex>()*trained.phi,z;
        measure("separate_forward_adjoint_factors",[&]{ReferenceFemContext f(space.operators());AdjointTestCache a;compare(ref,f.solve(loads));z=a.solve(space.operators(),rhs);});
        measure("unified_symmetric_factor_candidate",[&]{ReferenceFemContext f(space.operators());compare(ref,f.solve(loads));ComplexMatrix candidate=f.solve(rhs.conjugate()).conjugate();compare(z,candidate);});
        measure("ell_rebuild",[&]{LodSpace rebuilt(H,h,16,input.ell+1,input.policy=="area"?InterpolationPolicy::ManuscriptAreaWeighted:InterpolationPolicy::ArchivedArithmetic,limits);});
        measure("ell_reference_reuse",[&]{LodSpace reused(space,input.ell+1);if(reused.reference_identity()!=space.reference_identity())throw std::runtime_error("reference reuse identity");});
    }
    // Batched load assembly must agree with the existing sequential RHS path.
    for(int j=0;j<loads.cols();++j)compare(loads.col(j),lod2d::helmholtz::assemble_helmholtz_load(h.mesh,problems[j].source,q,problems[j].quadrature_context));
    std::cout<<std::setprecision(17)<<"{\"schema\":1,\"threads\":"<<input.threads<<",\"rhs\":"<<loads.cols()<<",\"reference_nodes\":"<<loads.rows()<<",\"rank\":"<<trained.phi.cols()<<",\"space_identity\":\""<<space.identity()<<"\",\"maximum_relative_difference\":"<<maximum_difference<<",\"seconds\":{";
    bool first=true;for(auto& [name,values]:times){if(!first)std::cout<<',';first=false;std::cout<<'"'<<name<<"\":";fixed::vector(values);}std::cout<<"},\"aot_solved_columns\":{";first=true;
    for(auto& [name,values]:columns){if(!first)std::cout<<',';first=false;std::cout<<'"'<<name<<"\":";fixed::vector(values);}
    dense_bytes["load_block"]=loads.size()*sizeof(Complex);dense_bytes["reference_block"]=ref.size()*sizeof(Complex);dense_bytes["accepted_solution"]=trained.accepted.values.size()*sizeof(Complex);
    std::cout<<"},\"accounted_dense_bytes\":{";first=true;for(auto& [name,bytes]:dense_bytes){if(!first)std::cout<<',';first=false;std::cout<<'"'<<name<<"\":"<<bytes;}std::cout<<"}}\n";
    return 0;
}catch(const std::exception& error){std::cerr<<"alod_profile: "<<error.what()<<'\n';return 1;}}
