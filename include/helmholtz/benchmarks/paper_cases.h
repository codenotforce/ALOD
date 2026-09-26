#pragma once
#include "alod/legacy_case_id.hpp"
#include "helmholtz/quadrature.h"
#include "helmholtz/manufactured.h"
namespace lod2d::helmholtz::benchmarks {
struct PaperCaseData {
    experiments::PaperCase id = experiments::PaperCase::R1;
    double wavenumber = 0;
    TriMesh initial_mesh;
    ComplexFunction source, exact, exact_laplacian;
    ComplexGradientFunction exact_gradient;
    QuadratureContext quadrature_context;
    // Optional joint evaluation avoids repeating the analytic envelope and phase.
    std::function<std::pair<Complex, Eigen::Vector2cd>(const Point2&)> exact_jet;
};
PaperCaseData make_shifted_r1_paper_case(double, const Point2&);
PaperCaseData make_parameterized_boundary_gaussian_s_paper_case(double, double, double, double, double, const Point2&);
}
