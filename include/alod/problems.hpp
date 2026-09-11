#pragma once

#include "helmholtz/benchmarks/paper_cases.h"
#include <string_view>

namespace alod {
using Problem = lod2d::helmholtz::benchmarks::PaperCaseData;

// Nominal problems only. Family scheduling is a later migration phase.
Problem make_problem(std::string_view id);
lod2d::helmholtz::QuadraturePolicy paper_quadrature(std::string_view id);
}
