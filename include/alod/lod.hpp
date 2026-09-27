#pragma once
#include <omp.h>
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

class LodPatchCache;
class LocalFactorCache;

// Fixed-state resource ceilings, checked before patch factors and dense RHS blocks.
struct LodLimits {
    int maximum_reference_nodes = 20000;
    std::size_t maximum_patch_entries = 8000000;
    std::size_t maximum_dense_entries = 8000000;
    int threads = omp_get_max_threads();
    bool parallel_riesz_gather = true;
    std::shared_ptr<LodPatchCache> patch_cache;
    std::shared_ptr<LocalFactorCache> riesz_cache;
};
struct LodSolution {
    ComplexMatrix values, coefficients;
    double pg_relative_residual = 0;
};
// Immutable after construction; ell-only spaces and estimator contexts share it.
struct LodHierarchyData {
    lod2d::TriMesh coarse;
    lod2d::RefineOutput reference;
    Sparse interpolation,coarse_basis,energy;
    lod2d::helmholtz::HelmholtzOperators operators;
    std::vector<int> coarse_nodes;
    std::string reference_identity;
    InterpolationPolicy policy=InterpolationPolicy::ManuscriptAreaWeighted;
};
class LodSpace {
public:
    LodSpace(lod2d::TriMesh coarse, lod2d::RefineOutput reference, double wavenumber,
             int ell, InterpolationPolicy policy = InterpolationPolicy::ManuscriptAreaWeighted,
             LodLimits limits = {});
    // Rebuild only ell-dependent correctors; reuse validated hierarchy operators.
    LodSpace(const LodSpace& reference_space, int ell);
    // Restore an accepted, checksummed trial basis without local patch solves.
    LodSpace(lod2d::TriMesh, lod2d::RefineOutput, double, int,
             InterpolationPolicy, LodLimits, const ComplexSparseMatrix& accepted_trial,
             const ComplexSparseMatrix* accepted_reduced=nullptr);
    // Restore another accepted basis over an immutable, already assembled hierarchy.
    LodSpace(std::shared_ptr<const LodHierarchyData>, int, InterpolationPolicy,
             LodLimits, const ComplexSparseMatrix&, const ComplexSparseMatrix* = nullptr);
    ~LodSpace();
    LodSpace(const LodSpace&) = delete;
    LodSpace& operator=(const LodSpace&) = delete;
    std::shared_ptr<const LodHierarchyData> hierarchy() const;
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
    const ComplexSparseMatrix& reduced() const;
    const std::vector<int>& coarse_nodes() const;
    const std::string& reference_identity() const;
    const std::string& identity() const;
    const LodLimits& limits() const;
    double patch_residual() const;
    double constraint_residual() const;
    // The factor belongs to this context. Calls on one context are sequential;
    // individual RHS columns are solved as a block, never concurrently on one LU.
    LodSolution solve(const ComplexMatrix& loads);
    // Sequential access to the state-owned factor; no additional factorization.
    ComplexMatrix solve_reduced(const ComplexMatrix& rhs) const;
private:
    LodSpace(lod2d::TriMesh, lod2d::RefineOutput, double, int, InterpolationPolicy, LodLimits, const LodSpace*, const ComplexSparseMatrix* = nullptr,const ComplexSparseMatrix* = nullptr, std::shared_ptr<const LodHierarchyData> = {});
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}
