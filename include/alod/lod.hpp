#pragma once
#include "helmholtz/operators.h"
#include "lod/quasi_interp.h"
#include "mesh/refine.h"
#include <memory>
#include <string>

namespace alod {
using lod2d::helmholtz::Complex;
using lod2d::helmholtz::ComplexMatrix;
using lod2d::helmholtz::ComplexVector;
using lod2d::helmholtz::ComplexSparseMatrix;
using Sparse = Eigen::SparseMatrix<double>;
using InterpolationPolicy = lod2d::QuasiInterpolationPolicy;

// Fixed-state resource ceilings, checked before patch factors and dense RHS blocks.
struct LodLimits {
    int maximum_reference_nodes = 20000;
    std::size_t maximum_patch_entries = 8000000;
    std::size_t maximum_dense_entries = 8000000;
    int threads = 1;
};
struct LodSolution {
    ComplexMatrix values, coefficients;
    double pg_relative_residual = 0;
};
class LodSpace {
public:
    LodSpace(lod2d::TriMesh coarse, lod2d::RefineOutput reference, double wavenumber,
             int ell, InterpolationPolicy policy = InterpolationPolicy::ManuscriptAreaWeighted,
             LodLimits limits = {});
    ~LodSpace();
    LodSpace(const LodSpace&) = delete;
    LodSpace& operator=(const LodSpace&) = delete;
    const lod2d::TriMesh& coarse() const;
    const lod2d::TriMesh& fine() const;
    const Sparse& prolongation() const;
    const Sparse& element_prolongation() const;
    const Sparse& interpolation() const;
    const Sparse& coarse_basis() const;
    const Sparse& energy() const;
    const lod2d::helmholtz::HelmholtzOperators& operators() const;
    const ComplexSparseMatrix& trial() const;
    const ComplexSparseMatrix& test() const;
    const std::vector<int>& coarse_nodes() const;
    const std::string& reference_identity() const;
    const std::string& identity() const;
    const LodLimits& limits() const;
    double patch_residual() const;
    double constraint_residual() const;
    // The factor belongs to this context. Calls on one context are sequential;
    // individual RHS columns are solved as a block, never concurrently on one LU.
    LodSolution solve(const ComplexMatrix& loads);
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}
