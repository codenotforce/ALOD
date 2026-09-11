#include "alod/problems.hpp"
#include <stdexcept>

namespace alod {
Problem make_problem(std::string_view id) {
    using namespace lod2d::helmholtz::benchmarks;
    if (id == "E1") {
        return make_shifted_r1_paper_case(16.0, lod2d::Point2(0.75, 0.5));
    }
    if (id == "E2") {
        return make_parameterized_boundary_gaussian_s_paper_case(
            16.0, 1.0, 0.5, 0.0, 80.0, lod2d::Point2(-0.5, 0.5));
    }
    throw std::invalid_argument("problem must be E1 or E2");
}

lod2d::helmholtz::QuadraturePolicy paper_quadrature(std::string_view id) {
    if (id != "E1" && id != "E2") throw std::invalid_argument("unknown problem");
    return {12, 16, 24, id == "E1" ? 8 : 6};
}
}
