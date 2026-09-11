#pragma once
#include "helmholtz/patch_system.h"
namespace lod2d::helmholtz {
struct HelmholtzPatchSolverConfig {
    int symbolic_cache_slots = 8;
    bool reuse_identical_factorization = false;
};
struct HelmholtzPatchSolveDiagnostics {
    double primal_residual=0, adjoint_residual=0, constraint_residual=0;
    double schur_residual=0, schur_rcond=1;
    bool symbolic_reused=false, factorization_reused=false;
};
struct HelmholtzPatchSolveResult {
    ComplexMatrix corrector, multipliers;
    HelmholtzPatchSolveDiagnostics diagnostics;
};
HelmholtzPatchSolveResult solve_helmholtz_patch(const HelmholtzPatchSystem&, const HelmholtzPatchSolverConfig& = {});
}
