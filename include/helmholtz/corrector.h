#pragma once
#include "helmholtz/patch_solver.h"
namespace lod2d::helmholtz {
struct HelmholtzCorrectorEntry {
    int row = -1;
    int local_coarse_vertex = -1;
    Complex value = 0.0;
};
using HelmholtzElementCorrector = std::vector<HelmholtzCorrectorEntry>;
ComplexSparseMatrix build_helmholtz_corrector_matrix(
    const TriMesh &coarse,
    int fine_node_count,
    const std::vector<HelmholtzElementCorrector> &correctors);

ComplexSparseMatrix build_helmholtz_corrected_basis(
    const Eigen::SparseMatrix<double> &coarse_to_fine,
    const TriMesh &coarse,
    int fine_node_count,
    const std::vector<HelmholtzElementCorrector> &correctors);

} // namespace lod2d::helmholtz
