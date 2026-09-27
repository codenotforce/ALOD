#pragma once
#include "alod/lod.hpp"
#include "alod/afem.hpp"

namespace alod {
class KernelDefectOperator;
enum class RieszPatchPolicy { ManuscriptN2, ArchivedSupportExpanded };
struct KernelRieszPatch {
    int coarse_node = -1;
    std::vector<int> coarse_hat_support, coarse_elements, discrete_dofs;
    std::vector<int> active_constraint_rows, independent_constraint_rows;
    Eigen::MatrixXd constraints;
};
struct ResidualRieszBatch {
    ComplexMatrix global_residuals;
    Eigen::MatrixXd node_eta_squared, element_eta_squared;
    Eigen::VectorXd eta;
    std::vector<std::vector<int>> marked_elements;
    int parallel_threads = 1;
    double patch_solve_seconds = 0;
};
class AdditiveKernelRieszContext {
public:
    struct Result {
        ComplexMatrix values, selected_values;
        Eigen::MatrixXd node_eta_squared;
        Eigen::VectorXd eta, selected_eta;
        double identity_relative_error = 0, selected_identity_relative_error = 0;
        double constraint_relative_residual = 0, local_relative_residual = 0, solve_seconds = 0;
    };
    explicit AdditiveKernelRieszContext(const LodSpace&,RieszPatchPolicy = RieszPatchPolicy::ManuscriptN2);
    ~AdditiveKernelRieszContext();
    AdditiveKernelRieszContext(const AdditiveKernelRieszContext&) = delete;
    AdditiveKernelRieszContext& operator=(const AdditiveKernelRieszContext&) = delete;
    // Concurrent applications share immutable factors; workspaces are call-local
    // and usage counters atomic. Construction/destruction still require ownership.
    Result apply(const ComplexMatrix& residual);
    // Same additive kernel inverse, without estimator/diagnostic side products.
    ComplexMatrix apply_action(const ComplexMatrix& residual);
    Result apply_selected(const ComplexMatrix& residual,const std::vector<int>& mask,bool full_estimator=false);
    ResidualRieszBatch estimate(const ComplexMatrix& loads,const ComplexMatrix& values,double theta);
    std::vector<int> regional_mask(double radius) const;
    const std::vector<KernelRieszPatch>& patches() const;
    const std::string& reference_identity() const;
    const std::string& patch_policy_name() const;
    double factorization_seconds() const;
    std::size_t factorizations() const;
    std::size_t applications() const;
    std::size_t applied_columns() const;
private:
    friend class KernelDefectOperator;
    Result apply_impl(const ComplexMatrix&,const std::vector<int>&,bool,bool);
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

enum class MemberRole { Train, Test, Shift, Pure };
struct FamilyMember { int id; MemberRole role; };
struct FamilyMarking {
    std::vector<int> member_ids, marked_elements;
    Eigen::VectorXd frozen_scales, mean_squared, worst_squared;
    int worst_member_id = -1;
    double mean_bulk = 0, worst_bulk = 0;
};
// Explicit training IDs select columns by identity, independent of input ordering.
// Only computed solution energies are accepted; audit members cannot drive marking.
FamilyMarking mark_family(const Eigen::MatrixXd& element_squared,
    const Eigen::VectorXd& solution_energy, const std::vector<FamilyMember>& members,
    const std::vector<int>& training_ids, double theta);

// This fine-grid strong residual controls reference refinement. It is distinct
// from the kernel estimator on coarse elements and from AFEM's own solution.
struct ReferenceStrongResidual {
    lod2d::helmholtz::adaptive::diagnostics::HelmholtzP1ResidualEstimate fine;
};
ReferenceStrongResidual reference_strong_residual(const LodSpace&,
    const ComplexVector& embedded_solution, const ComplexVector& load,
    const lod2d::helmholtz::ComplexFunction& source,
    const lod2d::helmholtz::QuadraturePolicy& quadrature = {},
    const lod2d::helmholtz::QuadratureContext& context = {});
}
