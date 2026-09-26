#pragma once
#include "alod/problems.hpp"
#include "alod/lod.hpp"
namespace alod {
// Parallelize over RHS while preserving every element/reduction order inside
// a member. The quadrature rule cache is thread-local; no shared factor is
// solved concurrently. Outer batches provide the dense-memory bound.
ComplexMatrix assemble_load_batch(const lod2d::TriMesh&,const std::vector<Problem>&,
    const lod2d::helmholtz::QuadraturePolicy&,int threads,
    std::vector<lod2d::helmholtz::SourceMomentData>* moments=nullptr);
struct ErrorBatch { Eigen::VectorXd exact_norm,exact_error,energy,reference_error; };
ErrorBatch integrate_error_batch(const lod2d::TriMesh&,const Sparse& energy,const ComplexMatrix&,
    const std::vector<Problem>&,const lod2d::helmholtz::QuadraturePolicy&,int threads);
ErrorBatch integrate_audit_batch(const lod2d::TriMesh&,const Sparse& energy,
    const ComplexMatrix& values,const ComplexMatrix& reference,const std::vector<Problem>&,
    const lod2d::helmholtz::QuadraturePolicy&,int threads,bool fused=true);
class ReferenceFemContext {
public:
    explicit ReferenceFemContext(const lod2d::helmholtz::HelmholtzOperators&);
    ~ReferenceFemContext();
    ComplexMatrix solve(const ComplexMatrix& loads);
    double factor_seconds() const;
    double solve_seconds() const;
    double relative_residual() const;
private:
    struct Impl;std::unique_ptr<Impl> impl_;
};
}
