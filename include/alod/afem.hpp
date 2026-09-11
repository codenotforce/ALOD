#pragma once
#include "helmholtz/operators.h"
#include <array>
namespace lod2d::helmholtz::adaptive {
namespace diagnostics {
struct ResidualEdgeContribution {
    Edge nodes{-1, -1};
    int left_element = -1;
    int right_element = -1;
    int left_parent = -1;
    int right_parent = -1;
    bool robin_boundary = false;
    bool neumann_boundary = false;
    double length = 0.0;
    double residual_l2_squared = 0.0;
    std::array<Complex, 2> residual_nodal{Complex(0.0), Complex(0.0)};
};
struct HelmholtzResidualContributions {
    std::vector<int> fine_element_parent;
    std::vector<double> body_l2_squared;
    std::vector<std::array<Complex, 3>> body_residual_nodal;
    std::vector<ResidualEdgeContribution> edges;
    ComplexVector reconstructed_residual;
    double algebraic_relative_difference = 0.0;
};
struct HelmholtzP1ResidualEstimate {
    std::vector<double> body_squared;
    std::vector<double> interior_jump_squared;
    std::vector<double> neumann_boundary_squared;
    std::vector<double> robin_boundary_squared;
    std::vector<double> element_squared;
    double eta = 0.0;
    double algebraic_relative_difference = 0.0;
};
HelmholtzP1ResidualEstimate estimate_conforming_p1_residual(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexVector &solution,
    const ComplexVector &load,
    const ComplexFunction &source,
    const QuadraturePolicy &quadrature = {},
    const QuadratureContext &quadrature_context = {});
}
std::vector<int> mark_doerfler(const std::vector<double>&, double, const std::vector<char>& = {});
}
