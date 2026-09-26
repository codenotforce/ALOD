#pragma once

#include "mesh/types.h"

#include <array>
#include <functional>
#include <vector>
#include <complex>
#include <memory>
#include <string>

namespace lod2d::helmholtz {

struct QuadraturePolicy {
    int base_triangle_order = 12;
    int gaussian_triangle_order = 20;
    int singular_triangle_order = 24;
    int max_recursive_subdivisions = 8;
};

enum class QuadratureClass {
    Regular,
    LocalizedGaussian,
    ReentrantSingular
};

struct QuadratureContext {
    QuadratureClass integrand_class = QuadratureClass::Regular;
    Point2 feature_point = Point2::Zero();
    double feature_scale = 0.0;
};

// Source-specific moments with shared quadrature mass matrices. The caller must
// pair them with the same frozen source; geometry and quadrature are checked.
struct SourceMomentData {
    std::string mesh_identity;
    QuadraturePolicy policy;
    QuadratureContext context;
    std::vector<double> squared;
    std::vector<std::array<std::complex<double>,3>> load;
    std::shared_ptr<std::vector<std::array<double,9>>> mass;
};

struct PhysicalTriangleQuadraturePoint {
    Point2 point = Point2::Zero();
    std::array<double, 3> barycentric{};
    double weight = 0.0;
};

void validate_quadrature_policy(const QuadraturePolicy &policy);

// Reuse caller-owned storage; no mesh-dependent global cache or callback lifetime.
void triangle_quadrature_points_into(
    const TriMesh &mesh, int element, const QuadraturePolicy &policy,
    const QuadratureContext &context,
    std::vector<PhysicalTriangleQuadraturePoint> &result);

std::vector<PhysicalTriangleQuadraturePoint> triangle_quadrature_points(
    const TriMesh &mesh,
    int element,
    const QuadraturePolicy &policy = {},
    const QuadratureContext &context = {});

double integrate_scalar_function(
    const TriMesh &mesh,
    const std::function<double(const Point2 &)> &function,
    const QuadraturePolicy &policy = {},
    const QuadratureContext &context = {});

bool operator==(const QuadraturePolicy &lhs, const QuadraturePolicy &rhs);

} // namespace lod2d::helmholtz
