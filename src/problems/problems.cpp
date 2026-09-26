#include "alod/problems.hpp"
#include <stdexcept>

namespace alod {
Problem make_problem(std::string_view id, double wavenumber) {
    using namespace lod2d::helmholtz::benchmarks;
    if (id == "E1") {
        return make_shifted_r1_paper_case(wavenumber, lod2d::Point2(0.75, 0.5));
    }
    if (id == "E2") {
        return make_parameterized_boundary_gaussian_s_paper_case(
            wavenumber, 1.0, 0.5, 0.0, 80.0, lod2d::Point2(-0.5, 0.5));
    }
    throw std::invalid_argument("problem must be E1 or E2");
}

double mesh_diameter(const lod2d::TriMesh& mesh) {
    double diameter=0;
    for(const auto& t:mesh.elems)for(int j=0;j<3;++j)
        diameter=std::max(diameter,(mesh.nodes[t[j]]-mesh.nodes[t[(j+1)%3]]).norm());
    return diameter;
}

lod2d::helmholtz::QuadraturePolicy paper_quadrature(std::string_view id) {
    if (id != "E1" && id != "E2") throw std::invalid_argument("unknown problem");
    return {12, 16, 24, id == "E1" ? 8 : 6};
}
}
