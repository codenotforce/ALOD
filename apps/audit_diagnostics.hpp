#pragma once
#include "alod/checkpoint.hpp"
#include "alod/lod.hpp"

// Offline transformations of an owned audit snapshot. These operations never
// update the accepted checkpoint or the adaptive controller's live state.
namespace audit_diagnostics {
void deepen_reference(std::unique_ptr<alod::LodSpace>& space, alod::Checkpoint& snapshot,
                      int levels, int ell, double wavenumber,
                      alod::InterpolationPolicy policy, const alod::LodLimits& limits);

// Energy Frobenius norm of the corrected-basis difference, not a single-RHS
// error or a certified operator-norm upper bound.
double two_level_corrector(const alod::LodSpace& space, int ell);

// Whole triangles selected by their centroid; this is not curved-cell clipping.
lod2d::TriMesh centroid_region(const lod2d::TriMesh& mesh, double radius);
}
