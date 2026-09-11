# Manuscript-to-implementation map

Manuscript SHA-256: `f714fb0ba66c2cc276437002aa2306d16271fd91fc88205ff8c3ca07ccdcbc37`. Exact labels and the current seven figures/eight tables are in `provenance/paper_inventory.json`.

| Definition | Archived implementation | Migration validation |
|---|---|---|
| Mixed-boundary Helmholtz, positive energy norm | `src/helmholtz/operators.cpp`, `boundary.cpp` | P1 FEM oracle, PDE/boundary residuals and integration convergence |
| E1 localized packet, E2 weighted corner plus Gaussian wave | `src/helmholtz/benchmarks/paper_cases.cpp` | Preserve production parameters; retain only E1/E2 factories |
| NVB and nested P1 injection | `src/mesh/refine.cpp` | Stable identities, ancestry, conformity and injection tests |
| Element L2 projection followed by vertex averaging | `src/lod/quasi_interp.cpp` | P0 counterexample below; coarse P1 reproduction |
| Two-sided corrected spaces and PG solve | `patch_system.cpp`, `patch_solver.cpp`, `corrector.cpp`, `model.cpp` | P1 uniform SLOD baseline; general adaptive machinery belongs to P2 |
| Kernel Riesz estimator and Theta | `adaptive/kernel_residual.cpp` and production workflow | P2; do not substitute AFEM strong residual |
| Candidate strong residual, reference-only alternation, lazy ell | E1 runner and `workflow_campaign_20260908.inc` | P0 old 27-state replay and exact marking recipe comparison |
| Regional AS, POD/AOT, inheritance and accepted publication | E2 runner and included workflow/audit helpers | P0 two-state accepted replay; full migration P3/P4 |
| Four E1 nominal/family controls | `family_control_common.h` and control entry points | Frozen 14208-row audit; new control migration deferred |

## Averaging discrepancy

The manuscript defines area-weighted averaging of coarse discontinuous L2 projections. Both archived implementations average equally over incident coarse elements. They agree when incident element areas are equal. Both reproduce coarse conforming P1 functions, so that identity alone cannot detect the difference.

`tests/quasi_average_probe.cpp` uses two triangles of areas 1/2 and 1, refines twice, and compares the operators on every continuous fine nodal basis function. The Frobenius difference is 0.15023130314433289; replacing the second area by 1/2 reduces it to 3.63e-16. The probe was first linked against the archived E2 library. Evidence is in `provenance/p0/quasi-average-counterexample.json`.

P1 preserves the archived numerical definition and only uses SLOD on uniform meshes. Before claiming manuscript equivalence on graded meshes in P2, explicitly choose and validate the operator definition. Do not silently change the archived operator or regenerate an oracle to conceal this discrepancy.
