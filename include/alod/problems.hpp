#pragma once

#include "helmholtz/benchmarks/paper_cases.h"
#include <string_view>

namespace alod {
using Problem = lod2d::helmholtz::benchmarks::PaperCaseData;

// Nominal factory. P2 explicit member tables use the parameterized factories;
// adaptive family scheduling remains a later migration phase.
Problem make_problem(std::string_view id, double wavenumber = 16.0);
double mesh_diameter(const lod2d::TriMesh& mesh);
lod2d::helmholtz::QuadraturePolicy paper_quadrature(std::string_view id);
}
