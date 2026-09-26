#include "helmholtz/benchmarks/paper_cases.h"
#include "helmholtz/boundary.h"
#include "alod/meshes.hpp"
#include "mesh/refine.h"
#include <cmath>
#include <stdexcept>
namespace lod2d::helmholtz::benchmarks {
namespace {
constexpr double kAlpha=2.0/3.0;
constexpr double kR1Localization=80.0;
struct SmoothStepValue {
    double value = 0.0;
    double first = 0.0;
    double second = 0.0;
};
struct SingularAmplitude {
    double value = 0.0;
    Eigen::Vector2d gradient = Eigen::Vector2d::Zero();
    double laplacian = 0.0;
};
struct SmoothWaveEnvelope {
    double value = 0.0;
    Eigen::Vector2d gradient = Eigen::Vector2d::Zero();
    double laplacian = 0.0;
};
TriMesh make_r1_mixed_boundary_mesh() {
    TriMesh mesh = make_helmholtz_unit_square_mesh();
    const auto [edges, boundary] = compute_edges(mesh);
    mesh.boundary_edges.clear();
    for (std::size_t index = 0; index < edges.size(); ++index) {
        if (!boundary[index]) continue;
        const Edge edge = edges[index];
        const Point2 midpoint = 0.5 * (
            mesh.nodes[edge[0]] + mesh.nodes[edge[1]]);
        BoundaryTag tag = BoundaryTag::Robin;
        if (std::abs(midpoint.y()) < 1e-14
            || std::abs(midpoint.y() - 1.0) < 1e-14) {
            tag = BoundaryTag::Dirichlet;
        } else if (std::abs(midpoint.x()) < 1e-14) {
            tag = BoundaryTag::Neumann;
        }
        mesh.boundary_edges.push_back({edge, tag});
    }
    synchronize_dirichlet_nodes(mesh);
    validate_boundary_tags(mesh);
    return mesh;
}
SmoothWaveEnvelope localized_r1_amplitude(
    const Point2 &point,
    const Point2 &center) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    const double x = point.x();
    const double y = point.y();
    const double one_minus_x = 1.0 - x;
    const double polynomial = x * x * one_minus_x * one_minus_x;
    const double polynomial_first = 2.0 * x - 6.0 * x * x
        + 4.0 * x * x * x;
    const double polynomial_second = 2.0 - 12.0 * x + 12.0 * x * x;
    const double sine = std::sin(pi * y);
    const double sine_first = pi * std::cos(pi * y);
    const double sine_second = -pi * pi * sine;
    const double dx = x - center.x();
    const double dy = y - center.y();
    const double gaussian = std::exp(
        -kR1Localization * (dx * dx + dy * dy));
    const double gaussian_x = -2.0 * kR1Localization * dx * gaussian;
    const double gaussian_y = -2.0 * kR1Localization * dy * gaussian;
    const double gaussian_laplacian = (
        4.0 * kR1Localization * kR1Localization * (dx * dx + dy * dy)
        - 4.0 * kR1Localization) * gaussian;

    SmoothWaveEnvelope result;
    result.value = polynomial * sine * gaussian;
    result.gradient.x() = sine * (
        polynomial_first * gaussian + polynomial * gaussian_x);
    result.gradient.y() = polynomial * (
        sine_first * gaussian + sine * gaussian_y);
    result.laplacian = polynomial_second * sine * gaussian
        + polynomial * sine_second * gaussian
        + polynomial * sine * gaussian_laplacian
        + 2.0 * polynomial_first * sine * gaussian_x
        + 2.0 * polynomial * sine_first * gaussian_y;
    return result;
}
SmoothWaveEnvelope boundary_weight(const Point2 &point) {
    const auto factor = [](const double coordinate) {
        const double one_minus_square = 1.0 - coordinate * coordinate;
        return SmoothStepValue{
            one_minus_square * one_minus_square,
            -4.0 * coordinate * one_minus_square,
            -4.0 + 12.0 * coordinate * coordinate};
    };
    const SmoothStepValue x = factor(point.x());
    const SmoothStepValue y = factor(point.y());
    SmoothWaveEnvelope result;
    result.value = x.value * y.value;
    result.gradient = Eigen::Vector2d(
        x.first * y.value, x.value * y.first);
    result.laplacian = x.second * y.value + x.value * y.second;
    return result;
}
SmoothWaveEnvelope boundary_gaussian_wave_envelope(
    const Point2 &point, const double alpha, const Point2 &center,
    const double normalization) {
    const SmoothWaveEnvelope weight = boundary_weight(point);
    const double x = point.x();
    const double y = point.y();
    const double polynomial = x * y;
    const Eigen::Vector2d polynomial_gradient(y, x);
    const double weighted_polynomial = polynomial * weight.value;
    const Eigen::Vector2d weighted_polynomial_gradient =
        polynomial_gradient * weight.value + polynomial * weight.gradient;
    const double weighted_polynomial_laplacian =
        polynomial * weight.laplacian
        + 2.0 * polynomial_gradient.dot(weight.gradient);

    const Eigen::Vector2d offset(
        x - center.x(), y - center.y());
    const double gaussian = std::exp(
        -alpha * offset.squaredNorm());
    const Eigen::Vector2d gaussian_gradient =
        -2.0 * alpha * gaussian * offset;
    const double gaussian_laplacian =
        (4.0 * alpha * alpha
             * offset.squaredNorm()
         - 4.0 * alpha) * gaussian;
    SmoothWaveEnvelope result;
    result.value = normalization * weighted_polynomial * gaussian;
    result.gradient = normalization * (
        gaussian * weighted_polynomial_gradient
        + weighted_polynomial * gaussian_gradient);
    result.laplacian = normalization * (
        gaussian * weighted_polynomial_laplacian
        + weighted_polynomial * gaussian_laplacian
        + 2.0 * weighted_polynomial_gradient.dot(gaussian_gradient));
    return result;
}
SingularAmplitude corner_singularity(const Point2 &point) {
    const double radius = point.norm();
    if (radius == 0.0) return {};
    double angle = std::atan2(point.y(), point.x());
    if (angle < 0.0) angle += 2.0 * std::acos(-1.0);
    if (angle > 1.5 * std::acos(-1.0) + 1e-12)
        throw std::invalid_argument(
            "singular manufactured solution evaluated outside the L-shaped domain");
    const double sine = std::sin(kAlpha * angle);
    const double cosine = std::cos(kAlpha * angle);
    const double radial_power = std::pow(radius, kAlpha);
    const double radial_derivative =
        kAlpha * std::pow(radius, kAlpha - 1.0);
    const Eigen::Vector2d radial(point.x() / radius, point.y() / radius);
    const Eigen::Vector2d angular(-radial.y(), radial.x());
    SingularAmplitude result;
    result.value = radial_power * sine;
    result.gradient = radial_derivative * (
        sine * radial + cosine * angular);
    return result;
}
PaperCaseData make_r1(double wavenumber, const Point2 &center) {
    PaperCaseData result;
    result.id = experiments::PaperCase::R1;
    result.wavenumber = wavenumber;
    result.initial_mesh = make_r1_mixed_boundary_mesh();
    result.quadrature_context.integrand_class = QuadratureClass::LocalizedGaussian;
    result.quadrature_context.feature_point = center;
    result.quadrature_context.feature_scale = 1.0 / std::sqrt(160.0);
    result.exact = [=](const Point2 &point) {
        const SmoothWaveEnvelope amplitude = localized_r1_amplitude(point, center);
        return amplitude.value
            * std::exp(Complex(0.0, wavenumber * point.x()));
    };
    result.exact_gradient = [=](const Point2 &point) {
        const SmoothWaveEnvelope amplitude = localized_r1_amplitude(point, center);
        const Complex phase = std::exp(
            Complex(0.0, wavenumber * point.x()));
        Eigen::Vector2cd gradient = phase
            * amplitude.gradient.cast<Complex>();
        gradient.x() += phase * Complex(0.0, wavenumber * amplitude.value);
        return gradient;
    };
    result.exact_jet = [=](const Point2 &point) {
        const SmoothWaveEnvelope amplitude = localized_r1_amplitude(point, center);
        const Complex phase = std::exp(Complex(0.0, wavenumber * point.x()));
        Eigen::Vector2cd gradient = phase * amplitude.gradient.cast<Complex>();
        gradient.x() += phase * Complex(0.0, wavenumber * amplitude.value);
        return std::make_pair(amplitude.value * phase, gradient);
    };
    result.exact_laplacian = [=](const Point2 &point) {
        const SmoothWaveEnvelope amplitude = localized_r1_amplitude(point, center);
        const Complex phase = std::exp(
            Complex(0.0, wavenumber * point.x()));
        return phase * (
            amplitude.laplacian
            + Complex(0.0, 2.0 * wavenumber * amplitude.gradient.x())
            - wavenumber * wavenumber * amplitude.value);
    };
    result.source = [=](const Point2 &point) {
        const SmoothWaveEnvelope amplitude = localized_r1_amplitude(point, center);
        const Complex phase = std::exp(
            Complex(0.0, wavenumber * point.x()));
        return -phase * (
            amplitude.laplacian
            + Complex(0.0, 2.0 * wavenumber * amplitude.gradient.x()));
    };
    return result;
}
}
PaperCaseData make_shifted_r1_paper_case(
    const double wavenumber,
    const Point2 &center) {
    if (!std::isfinite(wavenumber) || !(wavenumber > 0.0))
        throw std::invalid_argument("R1 wavenumber must be positive");
    if (!(center.x() > 0.0 && center.x() < 1.0
          && center.y() > 0.0 && center.y() < 1.0)) {
        throw std::invalid_argument("shifted R1 center must lie in the unit square");
    }
    return make_r1(wavenumber, center);
}
PaperCaseData make_parameterized_boundary_gaussian_s_paper_case(
    const double wavenumber,
    const double singular_coefficient,
    const double smooth_wave_amplitude,
    const double wave_phase,
    const double gaussian_alpha,
    const Point2 &center) {
    if (!std::isfinite(wavenumber) || !std::isfinite(singular_coefficient) || !(wavenumber > 0.0) || !(singular_coefficient >= 0.0)
        || !(smooth_wave_amplitude >= 0.0)
        || !(smooth_wave_amplitude <= 1.0)
        || !std::isfinite(wave_phase)
        || !(gaussian_alpha > 0.0)
        || !(center.x() < 0.0 && center.x() > -1.0)
        || !(center.y() > 0.0 && center.y() < 1.0)) {
        throw std::invalid_argument(
            "invalid parameterized boundary-Gaussian case-S parameters");
    }
    if (gaussian_alpha != 80.0) throw std::invalid_argument("E2 production alpha must be 80");
    PaperCaseData result;
    result.id = experiments::PaperCase::S;
    result.wavenumber = wavenumber;
    result.initial_mesh = make_helmholtz_l_shape_mesh();
    result.quadrature_context.integrand_class = QuadratureClass::ReentrantSingular;
    result.quadrature_context.feature_point = Point2::Zero();
    result.quadrature_context.feature_scale = 0.25;
    const double wave_alpha = gaussian_alpha;
    const double singular_scale = singular_coefficient;
    const Complex wave_coefficient = smooth_wave_amplitude * std::exp(Complex(0.0, wave_phase));
    const Point2 wave_center = center;
    // Parameter-only normalization is shared by every quadrature point.
    const double wave_normalization = 1.0 /
        (wave_center.x() * wave_center.y() * boundary_weight(wave_center).value);
    result.exact = [=](const Point2 &point) {
        const SingularAmplitude corner = corner_singularity(point);
        const SmoothWaveEnvelope weight = boundary_weight(point);
        const SmoothWaveEnvelope wave =
            boundary_gaussian_wave_envelope(point, wave_alpha, wave_center, wave_normalization);
        const Complex phase = std::exp(Complex(
            0.0, wavenumber * (point.x() - wave_center.x())));
        return singular_scale * weight.value * corner.value
            + wave_coefficient * wave.value * phase;
    };
    result.exact_gradient = [=](const Point2 &point) {
        const SingularAmplitude corner = corner_singularity(point);
        const SmoothWaveEnvelope weight = boundary_weight(point);
        const SmoothWaveEnvelope wave =
            boundary_gaussian_wave_envelope(point, wave_alpha, wave_center, wave_normalization);
        const Complex phase = std::exp(Complex(
            0.0, wavenumber * (point.x() - wave_center.x())));
        Eigen::Vector2cd gradient =
            singular_scale * (weight.value * corner.gradient
             + corner.value * weight.gradient).cast<Complex>();
        gradient += wave_coefficient * phase
            * wave.gradient.cast<Complex>();
        gradient.x() += wave_coefficient * phase
            * Complex(0.0, wavenumber * wave.value);
        return gradient;
    };
    result.exact_jet = [=](const Point2 &point) {
        const SingularAmplitude corner = corner_singularity(point);
        const SmoothWaveEnvelope weight = boundary_weight(point);
        const SmoothWaveEnvelope wave =
            boundary_gaussian_wave_envelope(point, wave_alpha, wave_center, wave_normalization);
        const Complex phase = std::exp(Complex(
            0.0, wavenumber * (point.x() - wave_center.x())));
        const Complex value = singular_scale * weight.value * corner.value
            + wave_coefficient * wave.value * phase;
        Eigen::Vector2cd gradient = singular_scale * (weight.value * corner.gradient
            + corner.value * weight.gradient).cast<Complex>();
        gradient += wave_coefficient * phase * wave.gradient.cast<Complex>();
        gradient.x() += wave_coefficient * phase
            * Complex(0.0, wavenumber * wave.value);
        return std::make_pair(value, gradient);
    };
    result.exact_laplacian = [=](const Point2 &point) {
        const SingularAmplitude corner = corner_singularity(point);
        const SmoothWaveEnvelope weight = boundary_weight(point);
        const SmoothWaveEnvelope wave =
            boundary_gaussian_wave_envelope(point, wave_alpha, wave_center, wave_normalization);
        const Complex phase = std::exp(Complex(
            0.0, wavenumber * (point.x() - wave_center.x())));
        const double singular_laplacian =
            corner.value * weight.laplacian
            + 2.0 * weight.gradient.dot(corner.gradient);
        return singular_scale * singular_laplacian
            + wave_coefficient * phase
                * (wave.laplacian
                   + Complex(0.0, 2.0 * wavenumber * wave.gradient.x())
                   - wavenumber * wavenumber * wave.value);
    };
    result.source = [=](const Point2 &point) {
        const SingularAmplitude corner = corner_singularity(point);
        const SmoothWaveEnvelope weight = boundary_weight(point);
        const SmoothWaveEnvelope wave =
            boundary_gaussian_wave_envelope(point, wave_alpha, wave_center, wave_normalization);
        const Complex phase = std::exp(Complex(
            0.0, wavenumber * (point.x() - wave_center.x())));
        const double singular_laplacian =
            corner.value * weight.laplacian
            + 2.0 * weight.gradient.dot(corner.gradient);
        return -singular_scale * singular_laplacian
            - singular_scale * wavenumber * wavenumber
                * weight.value * corner.value
            - wave_coefficient * phase
                * (wave.laplacian
                   + Complex(0.0, 2.0 * wavenumber * wave.gradient.x()));
    };
    return result;
}
}
