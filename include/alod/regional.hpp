#pragma once
#include "alod/estimator.hpp"
#include <memory>
namespace alod {
// A single reference-operator factor, independent of RHS, ell and dictionary.
class AdjointTestCache {
public:
    AdjointTestCache();
    ~AdjointTestCache();
    ComplexMatrix solve(const lod2d::helmholtz::HelmholtzOperators&,const ComplexMatrix& rhs);
    std::string identity() const;
    std::size_t factorizations() const;
    std::size_t solved_columns() const;
    double relative_residual() const;
private:
    struct Impl;std::unique_ptr<Impl> impl_;
};
enum class EnrichmentTests { KernelLift, ArchivedAdjoint };
struct RegionalConfig { double radius=.6; int rank_cap=24; bool inherit=true; EnrichmentTests tests=EnrichmentTests::KernelLift; };
struct RegionalEvaluation {
    ComplexMatrix values, seeds, raw_tests, tests;
    Eigen::VectorXd eta;
    double kernel_lift_residual=0;
    double pg_residual=0, aot_residual=0, raw_base_block=0, raw_dictionary_block=0;
};
class RegionalEvaluator {
public:
    RegionalEvaluator(const LodSpace&,AdditiveKernelRieszContext&,AdjointTestCache&,bool reuse=true,EnrichmentTests=EnrichmentTests::KernelLift);
    RegionalEvaluator(const LodSpace&,AdjointTestCache&,bool reuse=true,EnrichmentTests=EnrichmentTests::KernelLift);
    ~RegionalEvaluator();
    RegionalEvaluation evaluate_frozen(const ComplexMatrix& loads,const ComplexMatrix& phi);
    RegionalEvaluation evaluate(const ComplexMatrix& loads,const ComplexMatrix& phi,const std::vector<int>& mask);
    // Accounted persistent dense Eigen buffers; excludes solver internals,
    // sparse storage and temporary expression allocations.
    std::size_t dense_cache_bytes() const;
private:
    struct Impl;std::unique_ptr<Impl> impl_;
};
struct RegionalResult {
    ComplexMatrix raw_kernel, phi, full_phi, full_values;
    Eigen::VectorXd targets, full_eta;
    RegionalEvaluation accepted;
    std::string stop;
    int working_rank=0, inherited_rank=0, evaluations=0, compression_trials=0, selected_patches=0;
    bool training_met=false, compression_accepted=false;
    double kernel_residual=0, base_orthogonality=0, gram_residual=0;
};
// Incoming representatives must already be injected to the current reference.
// A mesh change repairs them by M_D E; an ell-only change only reprojects them.
RegionalResult train_regional(LodSpace&,AdditiveKernelRieszContext&,AdjointTestCache&,
    const ComplexMatrix& training_loads,const RegionalConfig& = {},
    const ComplexMatrix& incoming_raw = {},bool mesh_changed=false);
RegionalEvaluation evaluate_regional(const LodSpace&,AdditiveKernelRieszContext&,AdjointTestCache&,
    const ComplexMatrix& loads,const ComplexMatrix& phi,const std::vector<int>& selected_mask,EnrichmentTests=EnrichmentTests::KernelLift);
}
