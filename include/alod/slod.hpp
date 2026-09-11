#pragma once
#include "alod/problems.hpp"
#include "helmholtz/operators.h"
namespace alod {
struct SlodSolution {
    lod2d::TriMesh fine;
    lod2d::helmholtz::ComplexVector values;
    double residual=0, patch_residual=0, constraint_residual=0;
};
SlodSolution solve_slod(const lod2d::TriMesh&, const Problem&, const lod2d::helmholtz::QuadraturePolicy&, int ell=3, int reference_gap=4);
}
