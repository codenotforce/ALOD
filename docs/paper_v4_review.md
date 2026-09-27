# Version 4 implementation review

Reviewed against `version4/arxiv.tex` and the supplied September 14–20 result
packages. Source hashes and the externally reported E3 endpoints are recorded
in [the evidence inventory](provenance/v4_review.json). These are supplied
results, not production runs made by this checkout.

## Implemented changes

| Item | Current implementation |
| --- | --- |
| Enriched tests | `J_ell phi = phi - T_ell I_H phi`, using the free coarse coefficients of the actual trial basis. Both growing and POD dictionaries use the same map. |
| Coupled solve | Full projected Petrov–Galerkin blocks and their Schur complement. Kernel tests do not satisfy the old adjoint triangular identities; no lower-left block is discarded. |
| Inheritance | Raw regional representatives remain kernel constrained and support restricted. After repair, projection and POD, their equality to the lifted tests is checked by regression tests. |
| Global adjoint solves | Zero enrichment AOT factorizations or solves in kernel mode. The archived adjoint implementation remains explicitly selectable. Reference FEM audits still require their own reference solve. |
| Localization | `q = max(norm_k(U_nom), 1e-12) * Theta / eta_nom`, threshold 0.3 in current main presets, alongside the manual absolute gate tau(k)=3.2/k. Same-mesh rechecks follow every promotion. Lazy timing remains initial, terminal, extra states, and every `2^ell` accepted states. |
| Ritz convergence | The normalized policy uses a residual relative to the eigenvalue, rather than `max(1, lambda)`. Its dense path computes the largest Hermitian eigenpair directly. Archived stopping remains available. |
| Degenerate local spaces | Exactly zero local kernels return exactly zero Riesz solutions. Poorly conditioned unconstrained patch blocks can fall back to the same constrained saddle equations, with residual checks and iterative refinement. |
| E1/E2 presets | Both use two reference sweeps, fractions 0.3 and 0.2, H6/h10, ell 2–4, normalized balance and kernel tests. E1 retains extra checks 35/44; E2 uses nestedness without the E1 local gap reserve. |
| E3 presets | `k = 8,16,32,64,128`; coarse levels 4,6,8,10,12 and reference gaps 3,4,5,6,7, giving `k H0 = 2 sqrt(2)` and `k^3 h0^2 = 8`. Each run owns its operators and correctors. |
| Uniform LOD baseline | `reference_gap` is configurable in 1–8, with legacy default 4 and fixed ell 3. This permits the E3 reference scaling. |
| Restart/audit isolation | Test policy and normalized balance settings participate in restart validation. Independent audits select the frozen policy; a conflicting test policy is rejected for enriched checkpoints. |

`raw_base_block` and `raw_dictionary_block` are retained for archive diagnostics.
They are not zero/identity acceptance gates in kernel mode. The relevant fields
are `kernel_lift_residual`, `pg_residual` and the dictionary invariants.

## Configuration and compatibility

Use `configs/adaptive/e1_main.json`, `e2_main.json`, and `e3_k*.json` for the
current paper policies. Old configuration files without the new fields are
interpreted by the Python runner as `enrichment_tests="adjoint"`,
`ell_ratio_mode="raw"`, `ell_threshold=0`, preserving historical replay.
The library and direct adaptive CLI default to kernel tests; the direct CLI
retains raw balance unless `--ell-ratio-mode=solution_scaled` is supplied.
Never resume a trajectory with changed mathematics; start a new output directory.

```sh
python3 tools/run_adaptive.py --config configs/adaptive/e2_main.json --output results/e2-v4
build/alod_geometry --paper-v4
python3 tools/make_paper_v4_campaign.py --output results/v4-campaign
```

The last command only prepares a generic seven-job sequential specification;
it is not the current ordered 32+32 production launcher. Use
`tools/make_production_campaign.py` for the current dependency barriers, audit
budgets and tau(k)=3.2/k. Current E3 supports exact-error target stopping at 1%
with exact-only audits, subject to its horizon/resource limits. A horizon stop
without a target crossing is not target completion. Current launchers impose no
CPU affinity or single-thread-BLAS requirement; earlier pinned measurements
remain historical evidence.

## Implementation status and remaining validation (2026-09-28)

| Original review item | Current disposition |
|---|---|
| Generic unconditional pre-training defect loop | Current production uses the user-selected lazy absolute and relative gates. An unconditional generic-paper controller is not a separate selectable mode; do not claim the two schedules are identical. |
| Exact-error target and exact-only audits | Implemented, including terminal/resume handling; used by E3. |
| Recovery and large-reference support | Periodic checkpoints, shared geometry, recovery, large-file loading and the 64-bit UMFPACK reference path are implemented. Sustained resource/fault testing remains separate. The original 1 GiB loader restriction is not the current contract. |
| Scientific diagnostic menu | Implemented: see [scientific workflows](scientific_workflow.md) for supported controls and their precise definitions. |
| Portable paper delivery | Five-figure/five-table general pipeline and explicit historical adapters are implemented. The later 25-cycle plot/28-cycle fit export has its own scoped entry point. |
| Main production | E1/E2 completed 28 cycles, 57 accepted states each; E3 k=8,16,32,64,128 completed with audits on frozen deployments. |
| Remaining experiments | New-binary family/nominal controls, full lazy/every and inherit/reset comparisons, sustained failure tests and full campaigns after subsequent optimizations remain separate. |

These implementation and production statements supersede the original pending
list without changing its supplied numerical evidence. The historical k=128
LOD control endpoint at error 0.01153597 was resource-stopped; it is distinct
from the new E3 ALOD target result.



The old 33-state E2/AOT campaigns remain historical evidence. They are not the
current kernel-lifted E2 paper target. Old triangular AOT performance candidates
must not be used to simplify the kernel-lifted PG equations.

## Bounded validation

[Validation evidence](provenance/v4_validation.json) records 33 broad regression
cases and 10 final targeted checks, all passing (36 distinct CTest cases), plus
20 Python unit tests. The final checks cover the repaired lifting invariant,
checkpoint resume/audits, E3 mesh diameters, eigenvalue-relative Ritz versus
dense eigenpairs, and the supplied E1/E2 initial states. Their comparison uses
relative tolerance 2e-7 and absolute tolerance 1e-10; no full-trajectory equality
is claimed. The two-load kernel profile has maximum relative discrepancy
2.154e-15 and zero enrichment AOT solve columns. Repository integrity checks
preserve the 14,208 E1 and 1,650 E2 historical audit rows and reject local
absolute paths. All new documentation is English; the original plan remains
an intentionally Chinese historical document with a versioned update.

## 2026-09-27 implementation update

Exact-target termination, exact-only audits, the scientific diagnostic menu and the version 4 result-delivery pipeline are now implemented; see [scientific workflows](scientific_workflow.md). These supersede the corresponding pending implementation statements above. The update itself provided bounded validation. Later main production completion and remaining new-binary comparisons are distinguished in the status table above. Delivery explicitly records the supplied E1 family-control interpolation/first-crossing discrepancy.
