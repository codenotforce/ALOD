# Manuscript-to-implementation map

Manuscript SHA-256: `f714fb0ba66c2cc276437002aa2306d16271fd91fc88205ff8c3ca07ccdcbc37`. Exact labels and the current seven figures/eight tables are in `provenance/paper_inventory.json`.

| Definition | Archived implementation | Migration validation |
|---|---|---|
| Mixed-boundary Helmholtz, positive energy norm | `src/helmholtz/operators.cpp`, `boundary.cpp` | P1 FEM oracle, PDE/boundary residuals and integration convergence |
| E1 localized packet, E2 weighted corner plus Gaussian wave | `src/helmholtz/benchmarks/paper_cases.cpp` | Preserve production parameters; retain only E1/E2 factories |
| NVB and nested P1 injection | `src/mesh/refine.cpp` | Stable identities, ancestry, conformity and injection tests |
| Element L2 projection followed by vertex averaging | `src/lod/quasi_interp.cpp` | P0 counterexample below; coarse P1 reproduction |
| Two-sided corrected spaces and PG solve | `patch_system.cpp`, `patch_solver.cpp`, `corrector.cpp`, `model.cpp` | P1 SLOD; P2 fixed nested pairs and batch PG, including graded grids |
| Kernel Riesz estimator and Theta | `adaptive/kernel_residual.cpp`, `additive_kernel_riesz.inc`, `certificates.cpp` | P2 explicit kernel inverses, geometric boundary vertices, per-element conservation, old/new fixtures and two-sided spectra |
| Candidate strong residual, reference-only alternation, lazy ell | E1 runner and `workflow_campaign_20260908.inc` | P3 m_ref controller and six-state archive comparison across promotion; preserves proposed-coarse closure in reference-only sweeps |
| Regional AS, POD/AOT, inheritance and accepted publication | E2 runner and included workflow/audit helpers | P4 mathematical gates, same-mesh 2-to-3-to-4 and strict base-target/Theta comparison; corrected-bootstrap enriched proximity is labelled separately |
| Four E1 nominal/family controls | `family_control_common.h` and control entry points | P3 explicit presets and four short 48-member contracts; full new campaigns deferred to P8 |

## Averaging discrepancy

The manuscript defines area-weighted averaging of coarse discontinuous L2 projections. Both archived implementations average equally over incident coarse elements. They agree when incident element areas are equal. Both reproduce coarse conforming P1 functions, so that identity alone cannot detect the difference.

`tests/quasi_average_probe.cpp` uses two triangles of areas 1/2 and 1, refines twice, and compares the operators on every continuous fine nodal basis function. The Frobenius difference is 0.15023130314433289; replacing the second area by 1/2 reduces it to 3.63e-16. The probe was first linked against the archived E2 library. Evidence is in `provenance/p0/quasi-average-counterexample.json`.

P1 preserves the archived numerical definition and only uses SLOD on uniform meshes. P2 explicitly defaults to manuscript area weighting and retains an arithmetic compatibility policy. Both assembly paths support the policy. The graded probe checks the manuscript operator against independently assembled projection and averaging factors. Historical fixtures remain unchanged.

## Vertex-patch discrepancy

The manuscript defines N⁰(z)={z}, N(z) as the vertex star and N²(z) as one further shared-vertex expansion. P2's manuscript Riesz policy uses this geometry, including every geometric Dirichlet vertex as a center.

The archived implementation first measures interpolation support through the complete coarse-parent supports of fine nodal hats, then uses `patch_layers = interpolation_support_layers + 1` expansions from the star. On the E2 regression grids the measured support distance is one, giving N³(z). This differs from the manuscript even on equal-area meshes. For the H2/h5, ell=1 nominal E2 state, both modes have exactly the same LOD solution and strong residual, while their kernel estimators and Theta differ. Numerical evidence is in `provenance/p2_validation.json`.

P2 preserves the complete archived support-expansion algorithm as an explicit compatibility policy. Old/new equality tests select both archived policies; manuscript tests use independent local kernel bases and the stated N² geometry. Historical curves have not been reinterpreted as manuscript-mode trajectories.

## Localization mapping

P4 adds an explicit E2 lazy policy with E1's initial/terminal/2^ell cadence and E2's unchanged 0.1 ratio threshold. The every-state policy remains available. The new E2 cold start also follows the manuscript's deflated normalized phi recurrence; the archived nominal raw-kernel bootstrap is a documented discrepancy, not silently carried into manuscript runs. See [P3/P4](p3_p4.md).

For real interpolation and complex-symmetric A, the adjoint defect `A* B_test` is the conjugate of the primal defect `A B_trial`. Real energy constraints imply equal spectra, checked on uniform and graded E1/E2 states. Taking their maximum is redundant in this setting; adding their energies is incorrect. P2 uses one explicitly labelled side and the uncorrected coarse energy denominator. It retains the actual archived small-space iteration and matrix-free four-vector Ritz core, with strict state identities and failure on nonconvergence. Ritz values and residuals are implementation-study evidence, not rigorous upper bounds. Details are in `p2.md`.
