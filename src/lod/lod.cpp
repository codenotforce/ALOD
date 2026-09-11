#include "alod/lod.hpp"
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
#include <stdexcept>

namespace alod {
using namespace lod2d;
using namespace lod2d::helmholtz;
struct LodSpace::Impl {
    TriMesh coarse;
    RefineOutput reference;
    Sparse interpolation, coarse_basis, energy;
    HelmholtzOperators operators;
    ComplexSparseMatrix trial, test, reduced;
    Eigen::SparseLU<ComplexSparseMatrix> solver;
    std::vector<int> coarse_nodes;
    std::string reference_identity, identity;
    LodLimits limits;
    double patch_residual = 0, constraint_residual = 0;
};

LodSpace::LodSpace(TriMesh coarse, RefineOutput reference, double k, int ell,
                   InterpolationPolicy policy, LodLimits limits) : impl_(std::make_unique<Impl>()) {
    auto& p = *impl_;
    const int nh = reference.mesh.nodes.size(), nH = coarse.nodes.size();
    if (!std::isfinite(k) || k <= 0 || ell < 1 || ell > 16 || limits.threads < 1
        || limits.maximum_reference_nodes < 1 || nh > limits.maximum_reference_nodes
        || nH == 0 || nh == 0 || coarse.elems.empty() || reference.mesh.elems.empty()
        || reference.P_node.rows() != nh || reference.P_node.cols() != nH
        || reference.P_elem.rows() != static_cast<int>(reference.mesh.elems.size())
        || reference.P_elem.cols() != static_cast<int>(coarse.elems.size())
        || (policy != InterpolationPolicy::ArchivedArithmetic && policy != InterpolationPolicy::ManuscriptAreaWeighted))
        throw std::invalid_argument("invalid LOD hierarchy, policy, parameters or reference resource limit");
    p.coarse = std::move(coarse); p.reference = std::move(reference); p.limits = limits;
    validate_boundary_tags(p.coarse);validate_boundary_tags(p.reference.mesh);
    if(static_cast<std::size_t>(nh)*nH>limits.maximum_patch_entries)
        throw std::runtime_error("LOD constraint workspace bound exceeds resource limit");
    std::vector<Eigen::Triplet<double>> entries;
    for (int e=0; e<static_cast<int>(p.reference.mesh.elems.size()); ++e)
        for (int j=0; j<3; ++j) entries.emplace_back(3*e+j,p.reference.mesh.elems[e][j],1);
    Sparse cg(3*p.reference.mesh.elems.size(),nh); cg.setFromTriplets(entries.begin(),entries.end());
    p.interpolation = build_quasi_interp(p.coarse,p.reference.mesh,p.reference.P_dg,cg,nh,nH,policy);
    std::vector<bool> boundary(nH,false);
    for (int i:dirichlet_nodes(p.coarse)) boundary[i]=true;
    entries.clear();
    std::vector<Eigen::Triplet<double>> expected_entries;
    int column=0;
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
    auto patches=build_patches(p.coarse,ell);
    const std::vector<TriMesh> no_meshes;
    const std::vector<Sparse> no_prolongations;
    HelmholtzPatchAssembler assembler(p.coarse,p.reference.mesh,p.reference.P_elem,p.reference.P_dg,
        p.interpolation,patches,no_meshes,no_prolongations,no_prolongations,p.operators);
    std::vector<HelmholtzElementCorrector> correctors(p.coarse.elems.size());
    // Assembly is shared and immutable; independent target factors are thread-local.
    std::vector<double> residuals(correctors.size()), constraints(correctors.size());
    std::vector<std::string> errors(correctors.size());
#pragma omp parallel for schedule(dynamic,1) num_threads(limits.threads)
    for (int target=0;target<static_cast<int>(correctors.size());++target) {
        try {
            auto system=assembler.assemble(target);
            if (static_cast<std::size_t>(system.constraints.size())>limits.maximum_patch_entries
                || static_cast<std::size_t>(system.helmholtz.rows())*(system.constraints.rows()+3)>limits.maximum_dense_entries)
                throw std::runtime_error("LOD patch resource limit exceeded");
            HelmholtzPatchSolverConfig config; config.symbolic_cache_slots=1;
            auto solved=solve_helmholtz_patch(system,config);
            residuals[target]=std::max(solved.diagnostics.primal_residual,solved.diagnostics.adjoint_residual);
            constraints[target]=solved.diagnostics.constraint_residual;
            for(int row=0;row<solved.corrector.rows();++row)
                for(int col=0;col<solved.corrector.cols();++col)
                    if(std::abs(solved.corrector(row,col))>1e-14)
                        correctors[target].push_back({system.local_vertices[row],col,solved.corrector(row,col)});
        } catch(const std::exception& e) { errors[target]=e.what(); }
    }
    for(const auto& error:errors) if(!error.empty())throw std::runtime_error(error);
    p.patch_residual=*std::max_element(residuals.begin(),residuals.end());
    p.constraint_residual=*std::max_element(constraints.begin(),constraints.end());
    if(p.patch_residual>1e-8 || p.constraint_residual>1e-8)throw std::runtime_error("LOD local residual gate failed");
    auto full=build_helmholtz_corrected_basis(p.reference.P_node,p.coarse,nh,correctors);
    std::vector<ComplexTriplet> values;
    for(int c=0;c<column;++c)
        for(ComplexSparseMatrix::InnerIterator it(full,p.coarse_nodes[c]);it;++it)values.emplace_back(it.row(),c,it.value());
    p.trial.resize(nh,column);p.trial.setFromTriplets(values.begin(),values.end());
    p.test=p.trial.conjugate();
    p.reduced=p.test.adjoint()*p.operators.system*p.trial;
    p.reduced.prune(Complex(0,0),1e-14);p.reduced.makeCompressed();p.solver.compute(p.reduced);
    if(p.solver.info()!=Eigen::Success)throw std::runtime_error("LOD coarse factorization failed");
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
    p.identity=p.reference_identity+":ell="+std::to_string(ell);
}
LodSpace::~LodSpace()=default;
LodSolution LodSpace::solve(const ComplexMatrix& loads) {
    auto& p=*impl_;
    if(loads.rows()!=p.trial.rows() || loads.cols()<1 || !loads.allFinite()
        || static_cast<std::size_t>(loads.size())>p.limits.maximum_dense_entries)
        throw std::invalid_argument("LOD load dimensions, values or dense resource limit");
    ComplexMatrix rhs=p.test.adjoint()*loads;
    LodSolution result;result.coefficients=p.solver.solve(rhs);
    if(p.solver.info()!=Eigen::Success || !result.coefficients.allFinite())throw std::runtime_error("LOD coarse solve failed");
    result.values=p.trial*result.coefficients;
    result.pg_relative_residual=(p.test.adjoint()*(p.operators.system*result.values-loads)).norm()/std::max(1.0,rhs.norm());
    if(result.pg_relative_residual>1e-9)throw std::runtime_error("LOD PG residual gate failed");
    return result;
}
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
const std::vector<int>& LodSpace::coarse_nodes()const{return impl_->coarse_nodes;}
const std::string& LodSpace::reference_identity()const{return impl_->reference_identity;}
const std::string& LodSpace::identity()const{return impl_->identity;}
const LodLimits& LodSpace::limits()const{return impl_->limits;}
double LodSpace::patch_residual()const{return impl_->patch_residual;}
double LodSpace::constraint_residual()const{return impl_->constraint_residual;}
}
