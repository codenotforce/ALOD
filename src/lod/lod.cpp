#include "alod/timing.hpp"
#include "alod/lod.hpp"
#include "alod/patch_cache.hpp"
#include "fingerprint.hpp"
#include "alod/mesh_state.hpp"
#include "helmholtz/boundary.h"
#include "helmholtz/corrector.h"
#include "helmholtz/patch_system.h"
#include "helmholtz/patch_solver.h"
#include "lod/patches.h"
#include <Eigen/SparseLU>
#include <algorithm>
#include <bit>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace alod {
using namespace lod2d;
using namespace lod2d::helmholtz;
struct LodSpace::Impl {
    std::shared_ptr<LodHierarchyData> hierarchy;
    TriMesh& coarse;
    RefineOutput& reference;
    Sparse &interpolation,&coarse_basis,&energy;
    HelmholtzOperators& operators;
    explicit Impl(std::shared_ptr<LodHierarchyData> h):hierarchy(std::move(h)),
        coarse(hierarchy->coarse),reference(hierarchy->reference),interpolation(hierarchy->interpolation),
        coarse_basis(hierarchy->coarse_basis),energy(hierarchy->energy),operators(hierarchy->operators),
        coarse_nodes(hierarchy->coarse_nodes),reference_identity(hierarchy->reference_identity){}
    ComplexSparseMatrix trial, test, reduced;
    Eigen::SparseLU<ComplexSparseMatrix> solver;
    std::vector<int>& coarse_nodes;
    std::string& reference_identity;
    std::string identity;
    LodLimits limits;
    InterpolationPolicy policy;
    double patch_residual = 0, constraint_residual = 0;
};

LodSpace::LodSpace(TriMesh coarse, RefineOutput reference, double k, int ell,
                   InterpolationPolicy policy, LodLimits limits)
    :LodSpace(std::move(coarse),std::move(reference),k,ell,policy,limits,nullptr){}
LodSpace::LodSpace(const LodSpace& previous,int ell)
    :LodSpace({}, {},previous.operators().wavenumber,ell,previous.impl_->policy,previous.limits(),&previous){}
LodSpace::LodSpace(TriMesh coarse,RefineOutput reference,double k,int ell,
                   InterpolationPolicy policy,LodLimits limits,const ComplexSparseMatrix& trial,const ComplexSparseMatrix* reduced)
    :LodSpace(std::move(coarse),std::move(reference),k,ell,policy,limits,nullptr,&trial,reduced){}
LodSpace::LodSpace(std::shared_ptr<const LodHierarchyData> hierarchy,int ell,InterpolationPolicy policy,
                   LodLimits limits,const ComplexSparseMatrix& trial,const ComplexSparseMatrix* reduced)
    :LodSpace({}, {},hierarchy?hierarchy->operators.wavenumber:0,ell,policy,limits,nullptr,&trial,reduced,hierarchy){}
LodSpace::LodSpace(TriMesh coarse,RefineOutput reference,double k,int ell,
                   InterpolationPolicy policy,LodLimits limits,const LodSpace* reused,const ComplexSparseMatrix* accepted_trial,const ComplexSparseMatrix* accepted_reduced,std::shared_ptr<const LodHierarchyData> shared) : impl_(std::make_unique<Impl>(reused?reused->impl_->hierarchy:shared?std::const_pointer_cast<LodHierarchyData>(shared):std::make_shared<LodHierarchyData>())) {
    auto& p = *impl_;
    const auto& source_reference=(reused||shared)?p.reference:reference;
    const auto& source_coarse=(reused||shared)?p.coarse:coarse;
    const int nh = source_reference.mesh.nodes.size(), nH = source_coarse.nodes.size();
    if (!std::isfinite(k) || k <= 0 || ell < 1 || ell > 16 || limits.threads < 1
        || limits.maximum_reference_nodes < 1 || nh > limits.maximum_reference_nodes
        || nH == 0 || nh == 0 || source_coarse.elems.empty() || source_reference.mesh.elems.empty()
        || source_reference.P_node.rows() != nh || source_reference.P_node.cols() != nH
        || source_reference.P_elem.rows() != static_cast<int>(source_reference.mesh.elems.size())
        || source_reference.P_elem.cols() != static_cast<int>(source_coarse.elems.size())
        || (policy != InterpolationPolicy::ArchivedArithmetic && policy != InterpolationPolicy::ManuscriptAreaWeighted))
        throw std::invalid_argument("invalid LOD hierarchy, policy, parameters or reference resource limit");
    if(!reused&&!shared){p.coarse=std::move(coarse);p.reference=std::move(reference);}
    if(shared&&shared->policy!=policy)throw std::invalid_argument("shared hierarchy interpolation policy mismatch");
    if(!reused&&!shared)p.hierarchy->policy=policy;
    p.limits=limits;p.policy=policy;
    validate_boundary_tags(p.coarse);validate_boundary_tags(p.reference.mesh);
    if(static_cast<std::size_t>(nh)*nH>limits.maximum_patch_entries)
        throw std::runtime_error("LOD constraint workspace bound exceeds resource limit");
    int column=0;
    {PhaseTimer timing("lod_hierarchy_operators",-1);
    if(reused||shared){
        column=p.coarse_nodes.size();
    }else {
    std::vector<Eigen::Triplet<double>> entries;
    for (int e=0; e<static_cast<int>(p.reference.mesh.elems.size()); ++e)
        for (int j=0; j<3; ++j) entries.emplace_back(3*e+j,p.reference.mesh.elems[e][j],1);
    Sparse cg(3*p.reference.mesh.elems.size(),nh); cg.setFromTriplets(entries.begin(),entries.end());
    p.interpolation = build_quasi_interp(p.coarse,p.reference.mesh,p.reference.P_dg,cg,nh,nH,policy);
    std::vector<bool> boundary(nH,false);
    for (int i:dirichlet_nodes(p.coarse)) boundary[i]=true;
    entries.clear();
    std::vector<Eigen::Triplet<double>> expected_entries;
    column=0;
    for (int i=0;i<nH;++i) if (!boundary[i]) {
        p.coarse_nodes.push_back(i); expected_entries.emplace_back(i,i,1);
        for (Sparse::InnerIterator it(p.reference.P_node,i);it;++it) entries.emplace_back(it.row(),column,it.value());
        ++column;
    }
    if (!column) throw std::invalid_argument("LOD needs free coarse nodes");
    Sparse expected(nH,nH); expected.setFromTriplets(expected_entries.begin(),expected_entries.end());
    if (Sparse(p.interpolation*p.reference.P_node-expected).norm()>1e-9)
        throw std::invalid_argument("reference embedding does not reproduce coarse P1");
    p.coarse_basis.resize(nh,column);p.coarse_basis.setFromTriplets(entries.begin(),entries.end());
    p.operators=assemble_helmholtz_operators(p.reference.mesh,k);
    p.energy=p.operators.stiffness+k*k*p.operators.mass;
    if (ComplexSparseMatrix(p.operators.system-ComplexSparseMatrix(p.operators.system.transpose())).norm()>1e-12*p.operators.system.norm())
        throw std::runtime_error("conjugate two-sided LOD requires complex symmetry");
    // Immutable generated embeddings and boundary/operator parameters define the state.
    FingerprintBuilder embeddings;
    for(const Sparse* matrix:{&p.reference.P_node,&p.reference.P_elem,&p.reference.P_dg,&p.interpolation}) {
        embeddings.add_i64(matrix->rows());embeddings.add_i64(matrix->cols());
        for(int col=0;col<matrix->outerSize();++col)for(Sparse::InnerIterator it(*matrix,col);it;++it) {
            embeddings.add_i64(it.row());embeddings.add_i64(it.col());embeddings.add_double(it.value());
        }
    }
    p.reference_identity=mesh_fingerprint(p.coarse)+":"+mesh_fingerprint(p.reference.mesh)+":"+embeddings.finish()
        +":"+std::to_string(std::bit_cast<std::uint64_t>(k))+":"+std::to_string(static_cast<int>(policy));
    }
    }
    if(accepted_trial){
        PhaseTimer timing("lod_basis_restore",-1);
        if(accepted_trial->rows()!=nh||accepted_trial->cols()!=column)
            throw std::invalid_argument("accepted LOD basis dimensions mismatch");
        for(int c=0;c<accepted_trial->outerSize();++c)for(ComplexSparseMatrix::InnerIterator it(*accepted_trial,c);it;++it)
            if(!std::isfinite(it.value().real())||!std::isfinite(it.value().imag()))throw std::invalid_argument("nonfinite accepted LOD basis");
        p.trial=*accepted_trial;
        const ComplexSparseMatrix defect=p.interpolation.cast<Complex>()*(p.trial-p.coarse_basis.cast<Complex>());
        if(defect.norm()>1e-8*std::max(1.,p.trial.norm()))throw std::invalid_argument("accepted LOD interpolation invariant failed");
    }else{
    if(limits.patch_cache)limits.patch_cache->begin_generation();
    auto patches=build_patches(p.coarse,ell);
    const std::vector<TriMesh> no_meshes;
    const std::vector<Sparse> no_prolongations;
    HelmholtzPatchAssembler assembler(p.coarse,p.reference.mesh,p.reference.P_elem,p.reference.P_dg,
        p.interpolation,patches,no_meshes,no_prolongations,no_prolongations,p.operators);
    if(limits.patch_cache)assembler.prepare_local_dependencies();
    std::vector<HelmholtzElementCorrector> correctors(p.coarse.elems.size());
    // Assembly is shared and immutable; independent target factors are thread-local.
    std::vector<double> residuals(correctors.size()), constraints(correctors.size());
    std::vector<std::string> errors(correctors.size());
{PhaseTimer timing("lod_local_correctors",-1);
#pragma omp parallel num_threads(limits.threads)
    {
#pragma omp for schedule(dynamic,1)
    for (int target=0;target<static_cast<int>(correctors.size());++target) {
        try {
            std::shared_ptr<const HelmholtzPatchSolveResult> saved;std::string cache_key;
            const bool cacheable=static_cast<bool>(limits.patch_cache);
            auto system=assembler.geometry(target);
            if(cacheable)cache_key=assembler.local_dependency_key(system);
            if(cacheable)saved=limits.patch_cache->find_shared(cache_key);
            const bool hit=static_cast<bool>(saved);
            // A hit never assembles a local matrix or performs constraint QR.
            if(!hit)system=assembler.assemble(target,false);
            std::vector<int> cache_rows;
            if(cacheable){
                cache_rows.resize(system.local_vertices.size());std::iota(cache_rows.begin(),cache_rows.end(),0);
                std::sort(cache_rows.begin(),cache_rows.end(),[&](int a,int b){
                    const auto& x=p.reference.mesh.nodes[system.local_vertices[a]];const auto& y=p.reference.mesh.nodes[system.local_vertices[b]];
                    return x.x()!=y.x()?x.x()<y.x():x.y()<y.y();
                });
            }
            const std::size_t local_size=system.local_vertices.size();
            const std::size_t constraint_rows=hit?saved->multipliers.rows():system.constraints.rows();
            if (local_size*constraint_rows>limits.maximum_patch_entries
                || local_size*(constraint_rows+3)>limits.maximum_dense_entries)
                throw std::runtime_error("LOD patch resource limit exceeded");
            if(!hit){HelmholtzPatchSolverConfig config;config.symbolic_cache_slots=1;
                auto result=std::make_shared<HelmholtzPatchSolveResult>(solve_helmholtz_patch(system,config));
                if(cacheable){
                    auto original=result->corrector;
                    for(int i=0;i<static_cast<int>(cache_rows.size());++i)result->corrector.row(i)=original.row(cache_rows[i]);
                }
                saved=std::move(result);}
            if(cacheable)limits.patch_cache->insert_shared(std::move(cache_key),saved);
            const auto& solved=*saved;
            residuals[target]=std::max(solved.diagnostics.primal_residual,solved.diagnostics.adjoint_residual);
            constraints[target]=solved.diagnostics.constraint_residual;
            for(int row=0;row<solved.corrector.rows();++row)
                for(int col=0;col<solved.corrector.cols();++col)
                    if(std::abs(solved.corrector(row,col))>1e-14)
                        correctors[target].push_back({system.local_vertices[cacheable?cache_rows[row]:row],col,solved.corrector(row,col)});
        } catch(const std::exception& e) { errors[target]=e.what(); }
    }
    release_helmholtz_patch_cache();
    }
    }
    for(const auto& error:errors) if(!error.empty())throw std::runtime_error(error);
    p.patch_residual=*std::max_element(residuals.begin(),residuals.end());
    p.constraint_residual=*std::max_element(constraints.begin(),constraints.end());
    if(p.patch_residual>1e-8 || p.constraint_residual>1e-8)throw std::runtime_error("LOD local residual gate failed");
    {
        PhaseTimer timing("lod_basis_assembly",-1);
        auto full=build_helmholtz_corrected_basis(p.reference.P_node,p.coarse,nh,correctors);
        // Element correctors are no longer needed. Release their potentially
        // large storage before allocating the extracted trial/test bases.
        std::vector<HelmholtzElementCorrector>().swap(correctors);
        std::vector<ComplexTriplet> values;
        for(int c=0;c<column;++c)
            for(ComplexSparseMatrix::InnerIterator it(full,p.coarse_nodes[c]);it;++it)values.emplace_back(it.row(),c,it.value());
        p.trial.resize(nh,column);p.trial.setFromTriplets(values.begin(),values.end());
    } // Full basis and triplets must not overlap the reduced factorization.
    }
    p.test=p.trial.conjugate();
    if(accepted_reduced&&accepted_reduced->cols()){
        PhaseTimer timing("lod_reduced_restore",-1);
        if(!accepted_trial||accepted_reduced->rows()!=column||accepted_reduced->cols()!=column)
            throw std::invalid_argument("accepted reduced matrix dimensions invalid");
        p.reduced=*accepted_reduced;
        for(int c=0;c<p.reduced.outerSize();++c)for(ComplexSparseMatrix::InnerIterator it(p.reduced,c);it;++it)
            if(!std::isfinite(it.value().real())||!std::isfinite(it.value().imag()))throw std::invalid_argument("nonfinite reduced matrix");
        ComplexMatrix probe(column,3);
        for(int r=0;r<column;++r){probe(r,0)=1.;probe(r,1)=double((r*17)%29)/29.;probe(r,2)=Complex(double((r*7)%31)/31.,double((r*11)%23)/23.);}
        const ComplexMatrix expected=p.test.adjoint()*(p.operators.system*(p.trial*probe));
        if((p.reduced*probe-expected).norm()>1e-9*std::max(1.,expected.norm()))throw std::invalid_argument("accepted reduced operator invariant failed");
    }else {PhaseTimer timing("lod_reduced_assembly",-1);
    p.reduced=p.test.adjoint()*p.operators.system*p.trial;
    p.reduced.prune(Complex(0,0),1e-14);p.reduced.makeCompressed();}
    {PhaseTimer timing("lod_reduced_factor",-1);p.solver.compute(p.reduced);}
    if(p.solver.info()!=Eigen::Success)throw std::runtime_error("LOD coarse factorization failed");
    p.identity=p.reference_identity+":ell="+std::to_string(ell);
}
LodSpace::~LodSpace()=default;
LodSolution LodSpace::solve(const ComplexMatrix& loads) {
    auto& p=*impl_;
    if(loads.rows()!=p.trial.rows() || loads.cols()<1 || !loads.allFinite()
        || static_cast<std::size_t>(loads.size())>p.limits.maximum_dense_entries)
        throw std::invalid_argument("LOD load dimensions, values or dense resource limit");
    ComplexMatrix rhs=p.test.adjoint()*loads;
    LodSolution result;result.coefficients=solve_reduced(rhs);
    if(p.solver.info()!=Eigen::Success || !result.coefficients.allFinite())throw std::runtime_error("LOD coarse solve failed");
    result.values=p.trial*result.coefficients;
    result.pg_relative_residual=(p.test.adjoint()*(p.operators.system*result.values-loads)).norm()/std::max(1.0,rhs.norm());
    if(result.pg_relative_residual>1e-9)throw std::runtime_error("LOD PG residual gate failed");
    return result;
}
ComplexMatrix LodSpace::solve_reduced(const ComplexMatrix& rhs) const {
    auto& p=*impl_;
    if(rhs.rows()!=p.reduced.rows()||rhs.cols()<1||!rhs.allFinite()
        ||static_cast<std::size_t>(rhs.size())>p.limits.maximum_dense_entries)
        throw std::invalid_argument("reduced RHS dimensions, values or resource limit");
    ComplexMatrix result=p.solver.solve(rhs);
    if(p.solver.info()!=Eigen::Success||!result.allFinite())throw std::runtime_error("reduced solve failed");
    return result;
}
std::shared_ptr<const LodHierarchyData> LodSpace::hierarchy()const{return impl_->hierarchy;}
const TriMesh& LodSpace::coarse()const{return impl_->coarse;}
const TriMesh& LodSpace::fine()const{return impl_->reference.mesh;}
const Sparse& LodSpace::prolongation()const{return impl_->reference.P_node;}
const Sparse& LodSpace::element_prolongation()const{return impl_->reference.P_elem;}
const Sparse& LodSpace::interpolation()const{return impl_->interpolation;}
const Sparse& LodSpace::coarse_basis()const{return impl_->coarse_basis;}
const Sparse& LodSpace::energy()const{return impl_->energy;}
const HelmholtzOperators& LodSpace::operators()const{return impl_->operators;}
const ComplexSparseMatrix& LodSpace::trial()const{return impl_->trial;}
const ComplexSparseMatrix& LodSpace::test()const{return impl_->test;}
const ComplexSparseMatrix& LodSpace::reduced()const{return impl_->reduced;}
const std::vector<int>& LodSpace::coarse_nodes()const{return impl_->coarse_nodes;}
const std::string& LodSpace::reference_identity()const{return impl_->reference_identity;}
const std::string& LodSpace::identity()const{return impl_->identity;}
const LodLimits& LodSpace::limits()const{return impl_->limits;}
double LodSpace::patch_residual()const{return impl_->patch_residual;}
double LodSpace::constraint_residual()const{return impl_->constraint_residual;}
}
