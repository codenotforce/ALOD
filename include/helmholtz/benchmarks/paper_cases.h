#pragma once

#include "alod/legacy_case_id.hpp"
#include "helmholtz/quadrature.h"
#include "helmholtz/manufactured.h"
#include "mesh/types.h"

#include <optional>
#include <string>
#include <string_view>

namespace lod2d::helmholtz::benchmarks {

struct PaperCaseData {
    experiments::PaperCase id = experiments::PaperCase::R1;
    double wavenumber = 0.0;
    TriMesh initial_mesh;
    ComplexFunction source;
    ComplexFunction exact;
    ComplexGradientFunction exact_gradient;
    ComplexFunction exact_laplacian;
    QuadratureContext quadrature_context;
    std::optional<double> gaussian_sigma;
    std::optional<double> gaussian_normalization;
    std::optional<double> singular_oscillatory_fraction;
    std::optional<double> singular_cutoff_outer_radius;
    bool singular_quintic_cutoff = false;
    std::optional<double> smooth_wave_amplitude;
    std::optional<std::string> singular_solution_profile;
};

PaperCaseData make_paper_case(
    experiments::PaperCase id,
    double wavenumber);

// R1 with the localized wave packet translated to a user-supplied interior
// center.  The polynomial/sine envelope continues to enforce the same mixed
// homogeneous boundary conditions, so only the right-hand side changes.
PaperCaseData make_shifted_r1_paper_case(
    double wavenumber,
    const Point2 &center);

// E2c mixed corner/Gaussian manufactured solution with only the oscillatory
// packet center translated.  The singular component, operator, domain, and
// homogeneous boundary conditions are unchanged.
PaperCaseData make_shifted_boundary_gaussian_s_paper_case(
    double wavenumber,
    double smooth_wave_amplitude,
    double gaussian_alpha,
    const Point2 &center);

// Parameterized E2c manufactured family used by the fixed-space multi-RHS
// transfer experiment.  The differential operator and meshes stay fixed;
// only the manufactured right-hand side changes.  wave_phase is a constant
// complex phase offset in radians.
PaperCaseData make_parameterized_boundary_gaussian_s_paper_case(
    double wavenumber,
    double singular_coefficient,
    double smooth_wave_amplitude,
    double wave_phase,
    double gaussian_alpha,
    const Point2 &center);

PaperCaseData make_paper_case(
    experiments::PaperCase id,
    double wavenumber,
    double singular_oscillatory_fraction,
    double singular_cutoff_outer_radius = 0.5,
    bool singular_quintic_cutoff = false,
    double smooth_wave_amplitude = 0.0,
    std::string_view singular_solution_profile = "radial-cutoff");

// Case S keeps the legacy parameter hooks for controlled comparisons.  The
// two-argument factory selects the revised paper benchmark for S: the corner
// singularity is multiplied by a smooth outer-boundary weight and a Gaussian
// oscillatory component of amplitude 0.25 is centered at (-1/2,1/2).  The
// explicit overload retains the historical radial-cutoff profile.  Diagnostic
// profile suffixes alpha80/alpha144 select narrower Gaussian envelopes;
// wave-only variants suppress the corner component for fixed-space amplitude
// calibration without changing the production default.

double normalized_gaussian_constant(
    double sigma,
    const Point2 &center = Point2(0.35, 0.55));

// Frozen C-infinity cut-off and its first two radial derivatives.
double singular_cutoff(double radius);
double singular_cutoff_prime(double radius);
double singular_cutoff_second(double radius);

} // namespace lod2d::helmholtz::benchmarks
