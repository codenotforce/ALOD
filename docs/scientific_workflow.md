# Scientific workflow

This guide covers diagnostics and paper data. [Configuration](configuration.md)
defines controls; [runtime](runtime.md) defines the retained inputs they need.
Frozen scientific spaces and new diagnostic spaces must be labelled separately.

## Exact-error termination

Adaptive configurations accept these additional fields:

```json
{
  "exact_target": 0.01,
  "exact_scope": "nominal",
  "audit_mode": "exact"
}
```

`exact_target` is a positive relative energy-error target, or `-1` to disable it. `exact_scope` is `nominal` or `training_max`; the latter requires every frozen training member to satisfy the target. Held-out members do not enter stopping or marking. A near-zero exact norm is an error, not successful termination. The check runs on each accepted candidate using the manufactured exact solution; it does not solve a reference FEM system or wait for asynchronous audits.

The accepted transaction records `exact_target_ratio` and `target_reached`. A target-terminal checkpoint remains terminal on resume, including when the outer horizon is increased. Changing the target or its member scope is a new mathematical experiment. `termination_reason`, `target_reached` and `complete_horizon` distinguish a reached target from an exhausted horizon. Existing node, cycle and state limits remain active: an experiment can finish its configured horizon without reaching its target. The final scheduled localization policy is unchanged; an unexpected error crossing does not introduce an extra localization promotion.

E3 presets now request the one-percent nominal target and exact-only asynchronous audits. Their existing horizons and resource ceilings remain bounds, not guarantees of reaching one percent.

## Exact-only audits

```sh
python tools/run_audit.py --checkpoint results/e3/checkpoints/latest --output results/e3-exact --audit-mode exact
```

This mode publishes the exact error, exact norm, relative error and solution energy. Reference error, reference gap and their normalized values are `null`, with `reference_status=not_computed` for ALOD. It neither allocates nor factors an auxiliary reference FEM system. Frozen members already present in an unmodified checkpoint reuse their accepted values and PG evidence. Other members are solved in the frozen accepted space. AFEM needs its ordinary FEM solve only for members whose accepted values are absent; this is its solution, not an auxiliary ALOD reference audit.

Full audits remain the default. They retain all e/f/g fields. Exact-only audits do not satisfy a full-audit publication claim. Audit cache identities distinguish full and exact requests, so resume cannot silently substitute one for the other. Timing output reports reference factor counts and whether accepted values were reused.

## Scientific diagnostic menu

```sh
python tools/scientific_diagnostics.py --menu
python tools/scientific_diagnostics.py --checkpoint results/e2/checkpoints/latest --output results/e2-diagnostics --diagnostics all --threads 32
python tools/scientific_diagnostics.py --checkpoint results/e2/checkpoints/latest --output results/e2-controls --diagnostics ell,rank,radius --ells 2,3 --ranks 1,2,4 --radii 0.3,0.6,0.9
```

The menu covers full/exact audits, rank-zero, independent per-RHS training, pure singular/smooth loads, same-mesh ell/rank/radius controls, reference deepening, physical-region error, two-level corrector norms, quadrature-order controls, E1 iterative residual refinement, localization defects and mesh export. Unsupported combinations are explicitly skipped in an `all` request; numerical failures are retained and cause a nonzero exit. Each subrun records its options, source checkpoint hash and numerical output checksum. Source checkpoints are immutable.

Same-mesh E2 controls retrain the complete frozen training family with the requested parameters, rather than reinterpreting old enrichment coefficients in another LOD space. Rank-zero explicitly removes enrichment. Reference deepening injects the frozen trial and dictionary into a uniformly refined reference mesh and resolves that fixed space; retraining is a separate option. Node limits are checked before refinement allocation. Two-level output is the energy Frobenius norm of the corrected-basis difference between ell and ell+1, not a certified localization bound. Localization diagnostics use the eigenvalue-relative Ritz residual check and the configured tolerance.

The physical-region diagnostic uses the union of triangles whose centroids lie in the origin-centred disk of the requested radius. It reports the selected element count and region exact norm. It is an explicitly defined discrete region; it is not curved-cell integration over the exact circular boundary.



## Error metrics

For one accepted state and RHS case,
\[
e=\|u_\mu-U_\mu\|_k,\quad f=\|u_\mu-u_{h,\mu}\|_k,\quad
g=\|u_{h,\mu}-U_\mu\|_k.
\]
`E`, `F`/`E_ref`, `G` normalize these by exact_norm; near-zero denominators
have explicit null/status handling. Only e<=f+g follows, not equality.
`PG_residual` is algebraic, not a solution error. AFEM reference-gap fields are
not applicable; exact-only reference fields are not computed, not zero.
Audits solve additional RHS cases in the frozen space without retraining;
`--fresh` deliberately trains from empty at the same mesh/ell, and E2 `--pure`
adds the two pure-component cases. These diagnostics never feed adaptation.

## One-command version 4 delivery

The plot renderer now preserves the tuned parameters of
`generate_kernel_lift_figures.py` (E1/E2) and `plot_results.py` (E3).
E3 writes `kappa_error_dof.pdf` and `.png`. Its ALOD curve now uses the same
initial-state plus completed-H-step convention as E1/E2. The one-percent cutoff
is evaluated on all accepted states before filtering; AFEM/LOD keep every point
in that prefix. If the target or resource stop occurs during an incomplete cycle,
that intermediate endpoint is omitted from the ALOD curve but retained in the
scientific tables and raw data. Displayed points are joined on logarithmic axes.
Rendering, table arithmetic, data import and packaging are separate modules;
see [architecture](architecture.md).

Install the optional plotting dependencies from `requirements-paper.txt`. The versioned inventory in `configs/paper_v4_inventory.json` is tied to the reviewed manuscript hash and contains five figures and five principal tables.

For native results, prepare a schema-1 input file mapping scientific roles to run directories, relative to that file:

```json
{
  "schema": 1,
  "runs": {
    "e1-ALOD": "runs/e1_main",
    "e1-AFEM": "runs/e1_afem",
    "e2-ALOD": "runs/e2_main",
    "k8-ALOD": "runs/e3_k8"
  }
}
```

```sh
python tools/paper_delivery.py --spec results/paper-input.json --executable build/alod_run --output results/paper-delivery
```

The minimal example intentionally omits roles: a complete delivery also requires E1/E2 AFEM, UFEM and LOD; ALOD/AFEM/LOD at all five E3 wavenumbers; and four E1 nominal/family controls. Missing evidence yields a report and nonzero exit. `--allow-partial` produces an explicitly partial package. Adaptive meshes can be exported directly from the final checkpoint without an additional scientific solve; baselines export their terminal meshes during execution.

The supplied historical export layout has a dedicated versioned importer:

```sh
python tools/paper_delivery.py --archive-e1 results/archive-q03 --archive-e2 results/archive-kernel-lift --archive-e3 results/archive-k8-128 --output results/paper-delivery
```

The package includes PDF/PNG figures, CSV/LaTeX tables, compressed canonical data, a portable generator, declared manuscript settings, observed run settings, checksums and a completion report. It can be regenerated without the original experiment directories:

```sh
python tools/paper_delivery.py --canonical results/paper-delivery/data/runs.json.gz --output results/paper-regenerated
```

Imports validate state continuity, unique state/member audits, frozen member coverage, finite positive errors and mesh connectivity. The original historical file hashes are retained, while host paths are excluded by a whitelist. Imported results remain labelled historical; they are not speed or accuracy evidence for the new executable. Resource-stopped E3 endpoints are explicitly distinguished from one-percent crossings. Tail windows and displayed E2 endpoints follow the reviewed manuscript/export conventions and are recorded in the report.

The supplied manuscript has a reproducible inconsistency: its E1 family-control prose describes the first state with online dimension at least 2000, but the published values are log-budget interpolation at 2000. Delivery therefore includes both `e1-family-control` (published interpolated values) and `e1-family-first-crossing` (actual states), with an explicit scientific note. It does not silently edit the manuscript or claim the two quantities are equal.



## Scoped figures and convergence fits

`tools/paper_snapshot.py` implements the September 28 export request. It merges
current ALOD runs with explicitly historical controls, plots initial state plus
E1/E2 cycle endpoints 1–25, and requires ALOD meshes from cycle 25. E3 retains
its completed-cycle/first-1%-prefix selection. The separate convergence report
fits nominal relative energy error against online DOFs on cycles 19–28 of the
full 28-cycle E1/E2 runs. It does not fit the truncated figure curves.

The general `paper_delivery.py` terminal tables retain their historical fitting
windows (including a five-cycle E1 ALOD tail). Use the scoped exporter for the
ten-cycle fit request; these two CLIs are not interchangeable. Both preserve
source identities and historical/current origin labels. The scoped package has
`full-runs.json.gz`, `plotted-runs.json.gz`, fit points and `reproduce.py`.


## Fine-grid stability diagnostic

With free-DOF Helmholtz matrix A and energy matrix E, the optional diagnostic
estimates gamma_h(k)^2=lambda_min(A* E^-1 A,E), using inverse actions and checked
Ritz residuals. This is fine-grid stability, not online LOD stability or a
certified extreme bound. It never edits production thresholds.

```sh
cmake --build build --target alod_stability
python3 tools/run_stability_calibration.py --executable build/alod_stability --output results/gamma-calibration
```

The resumable sweep uses k=8,16,32,64,128, reference levels 7,10,13,16,19,
independent starts and mesh-sensitivity offsets. Use a fresh output after
changing executable/definition. `--problem E2` defines a separate-domain study;
E1 estimates cannot be silently reused for E2. Gamma-based production calibration
was cancelled; the deployed threshold is the prescribed 3.2/k rule.
