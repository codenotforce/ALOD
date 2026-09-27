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
| Localization | `q = max(norm_k(U_nom), 1e-12) * Theta / eta_nom`, threshold 0.3 in all current paper presets. Same-mesh rechecks follow every promotion. Lazy timing remains initial, terminal, extra states, and every `2^ell` accepted states. |
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

The last command only prepares seven adaptive jobs. Its conservative resource
allowances require server inspection before launch. Follow the plan's PowerShell
SSH, fresh server directory, tmux, distinct physical cores and BLAS-one rules.
E3 presets use the supplied accepted-state horizons; this is not an automatic
exact-error target stop. Do not infer target completion from the horizon alone.

## Remaining implementation and validation

Execution update: [lazy absolute checks and asynchronous recovery](async_execution.md) now implement the requested manual lazy absolute gate, concurrent auditing, automatic checkpoint boundaries/recovery, and large-file reader support. The numbered list below records the original review; these completed portions supersede its first and third items.

1. **Generic paper Algorithm 1 versus the reported production controller.**
   The paper has an unconditional pre-training loop for `Theta > tau_Theta`.
   The supplied production code uses scheduled nominal balance checks, with
   E1's absolute acceptance branch after the ratio decision. This review follows
   those result-producing policies. The generic pre-training defect loop is
   not implemented as a separate selectable controller; the two must not be
   described as identical.
2. **Adaptive exact-error target termination and exact-only audits.** The supplied
   later E3 workflow can stop after reaching 1% without assembling a full
   reference FEM audit. This checkout still uses bounded accepted-state horizons
   and its independent full audits. A resumable target-terminal check and an
   exact-only audit mode remain to be integrated.
3. **Large-run recovery and checkpoint throughput.** The later source archive has
   UMFPACK allocation/order/index-width recovery and larger checkpoint capacity.
   These have not been imported wholesale. The current checkpoint loader's
   default 1 GiB limit, large reference factors, late-run disk use and scheduler
   fault recovery need server validation before full E3/family campaigns.
4. **Remaining diagnostics.** Automated E2 same-final-mesh ell/rank/radius
   controls, reference deepening, physical-region diagnostics and optional
   two-level corrector norms are not all available as dedicated audited commands.
   Rank-zero, fresh and pure-RHS audits exist; E1 ell and quadrature overrides
   exist. Kernel-mode performance measurements must replace old AOT claims.
5. **Portable paper results and figures.** P7–P9 still need a versioned importer,
   complete audit/coverage checks and generators for the current manuscript's
   figure/table inventory. The small initial-state fixtures added here are not
   a substitute for that pipeline or for complete trajectory comparisons.
6. **Deferred server experiments.** Full E1/E2 51-state curves, family/nominal
   controls, E2 lazy/every and inherit/reset comparisons, E3 through k=128,
   common-budget endpoint controls and late-state resource measurements remain
   unvalidated for this checkout. The supplied k=128 LOD endpoint is explicitly
   resource-stopped at error 0.01153597; it did not reach 1%.

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

Exact-target termination, exact-only audits, the scientific diagnostic menu and the version 4 result-delivery pipeline are now implemented; see [scientific workflows](scientific_workflow.md). These supersede the corresponding pending implementation statements above. Full new-binary production experiments and late-state resource validation remain deferred. Delivery explicitly records the supplied E1 family-control interpolation/first-crossing discrepancy.
