#include "alod/slod.hpp"
#include "helmholtz/corrector.h"
#include "helmholtz/boundary.h"
#include "lod/quasi_interp.h"
#include "lod/patches.h"
#include "mesh/refine.h"
#include <Eigen/SparseLU>
#include <algorithm>
#include <stdexcept>

namespace alod {
SlodSolution solve_slod(const lod2d::TriMesh& coarse,const Problem& problem,
                       const lod2d::helmholtz::QuadraturePolicy& quad,int ell,int reference_gap) {
    using namespace lod2d;using namespace lod2d::helmholtz;
    if(ell!=3 || reference_gap!=4) throw std::invalid_argument("P1 SLOD requires ell=3 and reference gap=4");
    auto areas=compute_area(coarse);
    if(*std::max_element(areas.begin(),areas.end())>*std::min_element(areas.begin(),areas.end())*(1+1e-12))
        throw std::invalid_argument("P1 SLOD requires a uniform coarse mesh");
    auto fine=refine_mesh_nvb(coarse,reference_gap);
    std::vector<Eigen::Triplet<double>> entries;
    for(int e=0;e<static_cast<int>(fine.mesh.elems.size());++e)
        for(int j=0;j<3;++j) entries.emplace_back(3*e+j,fine.mesh.elems[e][j],1);
    Eigen::SparseMatrix<double> cg(3*fine.mesh.elems.size(),fine.mesh.nodes.size());
    cg.setFromTriplets(entries.begin(),entries.end());
    auto I=build_quasi_interp(coarse,fine.mesh,fine.P_dg,cg,fine.mesh.nodes.size(),coarse.nodes.size());
    std::vector<Eigen::Triplet<double>> expected_entries;
    std::vector<bool> boundary(coarse.nodes.size(),false);
    for(int i:dirichlet_nodes(coarse))boundary[i]=true;
    for(int i=0;i<static_cast<int>(boundary.size());++i)if(!boundary[i])expected_entries.emplace_back(i,i,1);
    Eigen::SparseMatrix<double> expected(coarse.nodes.size(),coarse.nodes.size());
    expected.setFromTriplets(expected_entries.begin(),expected_entries.end());
    if(Eigen::SparseMatrix<double>(I*fine.P_node-expected).norm()>1e-9)
        throw std::runtime_error("SLOD quasi-interpolation does not reproduce coarse P1");
    auto patches=build_patches(coarse,ell);
    auto operators=assemble_helmholtz_operators(fine.mesh,problem.wavenumber);
    // Direct Schur does not need the geometric V-cycle hierarchy.
    const std::vector<TriMesh> no_meshes;
    const std::vector<Eigen::SparseMatrix<double>> no_prolongations;
    HelmholtzPatchAssembler assembler(coarse,fine.mesh,fine.P_elem,fine.P_dg,I,patches,
        no_meshes,no_prolongations,no_prolongations,operators);
    std::vector<HelmholtzElementCorrector> correctors(coarse.elems.size());
    SlodSolution result;
    for(int target=0;target<static_cast<int>(correctors.size());++target) {
        auto system=assembler.assemble(target);
        auto solved=solve_helmholtz_patch(system);
        result.patch_residual=std::max({result.patch_residual,solved.diagnostics.primal_residual,solved.diagnostics.adjoint_residual});
        result.constraint_residual=std::max(result.constraint_residual,solved.diagnostics.constraint_residual);
        for(int row=0;row<solved.corrector.rows();++row)
            for(int col=0;col<solved.corrector.cols();++col)
                if(std::abs(solved.corrector(row,col))>1e-14)
                    correctors[target].push_back({system.local_vertices[row],col,solved.corrector(row,col)});
    }
    if(result.patch_residual>1e-8 || result.constraint_residual>1e-8)
        throw std::runtime_error("SLOD local equation/constraint residual gate failed");
    auto full=build_helmholtz_corrected_basis(fine.P_node,coarse,fine.mesh.nodes.size(),correctors);
    std::vector<bool> is_dirichlet(coarse.nodes.size(),false);
    for(int node:dirichlet_nodes(coarse)) is_dirichlet[node]=true;
    std::vector<ComplexTriplet> triplets;
    int col=0;
    for(int i=0;i<full.outerSize();++i) if(!is_dirichlet[i]) {
        for(ComplexSparseMatrix::InnerIterator it(full,i);it;++it)triplets.emplace_back(it.row(),col,it.value());
        ++col;
    }
    if(col==0)throw std::invalid_argument("SLOD mesh has no free coarse nodes");
    ComplexSparseMatrix trial(full.rows(),col);trial.setFromTriplets(triplets.begin(),triplets.end());
    ComplexSparseMatrix test=trial.conjugate();
    ComplexSparseMatrix A=test.adjoint()*operators.system*trial;
    A.prune(Complex(0,0),1e-14);A.makeCompressed();
    Eigen::SparseLU<ComplexSparseMatrix> solver;solver.compute(A);
    if(solver.info()!=Eigen::Success)throw std::runtime_error("SLOD coarse factorization failed");
    auto load=assemble_helmholtz_load(fine.mesh,problem.source,quad,problem.quadrature_context);
    ComplexVector rhs=test.adjoint()*load, coefficients=solver.solve(rhs);
    if(solver.info()!=Eigen::Success || !coefficients.allFinite())throw std::runtime_error("SLOD coarse solve failed");
    result.values=trial*coefficients;
    result.residual=(test.adjoint()*(operators.system*result.values-load)).norm()/std::max(1.0,rhs.norm());
    if(result.residual>1e-9)throw std::runtime_error("SLOD PG residual gate failed");
    result.fine=std::move(fine.mesh);return result;
}
}
