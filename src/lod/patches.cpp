#include "lod/patches.h"
#include <Eigen/Sparse>
#include <stdexcept>
#include <vector>

namespace lod2d {
namespace {

Eigen::SparseMatrix<double> element_adjacency(const TriMesh &coarse) {
    const int NH = static_cast<int>(coarse.nodes.size());
    const int NTH = static_cast<int>(coarse.elems.size());
    std::vector<Eigen::Triplet<double>> ivt_t;
    for (int t = 0; t < NTH; ++t) {
        for (int i = 0; i < 3; ++i) {
            int v = coarse.elems[t][i];
            ivt_t.emplace_back(v, t, 1.0);
        }
    }
    Eigen::SparseMatrix<double> Ivt(NH, NTH);
    Ivt.setFromTriplets(ivt_t.begin(), ivt_t.end());

    Eigen::SparseMatrix<double> Itt = (Ivt.transpose() * Ivt).pruned();
    for (int k = 0; k < Itt.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Itt, k); it; ++it) {
            it.valueRef() = 1.0;
        }
    }
    Itt.makeCompressed();
    return Itt;
}

Eigen::SparseMatrix<double> binary_with_identity(
    const Eigen::SparseMatrix<double> &patch) {
    if (patch.rows() != patch.cols())
        throw std::invalid_argument("patch matrix must be square");
    std::vector<Eigen::Triplet<double>> bin_t;
    for (int k = 0; k < patch.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(patch, k); it; ++it) {
            bin_t.emplace_back(it.row(), it.col(), 1.0);
        }
    }
    // Diagonal (MATLAB: patch = spones(patch) + speye(NTH))
    for (int i = 0; i < patch.rows(); ++i) {
        bin_t.emplace_back(i, i, 1.0);
    }

    Eigen::SparseMatrix<double> result(patch.rows(), patch.cols());
    result.setFromTriplets(bin_t.begin(), bin_t.end());
    result.makeCompressed();
    return result;
}

} // namespace

Eigen::SparseMatrix<double> expand_patches(
    const TriMesh &coarse,
    const Eigen::SparseMatrix<double> &base,
    const int additional_layers) {
    const int element_count = static_cast<int>(coarse.elems.size());
    if (additional_layers < 0)
        throw std::invalid_argument(
            "additional patch layers must be nonnegative");
    if (base.rows() != element_count || base.cols() != element_count)
        throw std::invalid_argument(
            "base patch matrix does not match the coarse mesh");
    Eigen::SparseMatrix<double> patch = base;
    if (additional_layers > 0) {
        const Eigen::SparseMatrix<double> adjacency =
            element_adjacency(coarse);
        for (int layer = 0; layer < additional_layers; ++layer)
            patch = adjacency * patch;
    }
    return binary_with_identity(patch);
}

Eigen::SparseMatrix<double> build_patches(const TriMesh &coarse, int ell) {
    if (ell < 0)
        throw std::invalid_argument(
            "patch oversampling level must be nonnegative");
    const int element_count = static_cast<int>(coarse.elems.size());
    Eigen::SparseMatrix<double> identity(element_count, element_count);
    identity.setIdentity();
    return expand_patches(coarse, identity, ell);
}

} // namespace lod2d
