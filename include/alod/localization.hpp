#pragma once
#include "alod/estimator.hpp"
#include <optional>
namespace alod {
struct LocalizationWarmStart { std::string identity; ComplexMatrix block; };
struct LocalizationEigenConfig {
    int maximum_iterations=750;
    double relative_tolerance=1e-9;
    bool eigenvalue_relative_residual=false; // archived mode used max(1,lambda)
    int dense_cross_check_max_dimension=64;
    int dense_fallback_max_dimension=0;
    std::optional<LocalizationWarmStart> warm_start;
};
struct LocalizationSpectrum {
    double lambda_max=0, relative_residual=0, dense_lambda_max=0, dense_relative_difference=0;
    int iterations=0;
    bool converged=false, dense_cross_checked=false, used_dense_fallback=false;
    bool used_sparse_generalized_solver=false, used_warm_start=false;
    ComplexVector dominant_vector;
    ComplexMatrix dominant_subspace;
};
struct LocalizationResult {
    double theta=0;
    LocalizationSpectrum spectrum;
    LocalizationWarmStart warm_start;
    std::string identity;
};
// Adjoint defect A* B_test and the uncorrected coarse energy denominator.
// A converged Ritz pair is an approximation, not a certified upper bound.
LocalizationResult localization_theta(const LodSpace&,AdditiveKernelRieszContext&,
    const LocalizationEigenConfig& = {},bool adjoint=true);
struct LocalizationRatio { std::string status; std::optional<double> value; };
LocalizationRatio localization_ratio(double theta,double nominal_eta);
}
