# Scientific workflows and current-paper delivery

The pre-existing reuse work was committed and pushed as `1c8410f`. This follow-up adds exact-target termination, exact-only audits, offline scientific controls, the version 4 result-delivery pipeline, cross-refinement integration reuse and background checkpoint geometry preparation.

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

## Incremental integration

Source assembly retains the current family's element load moments, squared-source integrals and quadrature mass matrices. After refinement, an unchanged ordered physical triangle reuses those values. New children are integrated normally. Frozen family identity, quadrature policy and localized-feature context are checked. Global node numbers do not authorize reuse. Reductions into the new load vector retain their original element order. The existing source-moment budget accounts conservatively for old and new buffers during the transition.

Exact-error integration can retain quadrature points and analytic values/gradients in a bounded 64 MiB cache per run target checker or resident audit worker. Geometry, family/member selection, wavenumber and quadrature settings participate in the identity. Discrete values and gradients are always recomputed, including on unchanged cells: they change when the numerical solution changes. A sampled working-set estimate bypasses retention when the mesh/family would exceed the budget, avoiding LRU churn on large grids. Concurrent audit workers own separate caches. This retention budget is additional to shared-snapshot and numerical-workspace budgets.

`ALOD_INTEGRATION_COLD=1` disables both incremental paths for numerical comparisons. Counters report reused source elements and exact-jet hits/retained bytes. Restart rebuilds these optional caches from the durable scientific state; cached integrals are not checkpoint requirements.

## Lower checkpoint blocking

The adaptive process starts a single background task to fingerprint, serialize, verify and durably publish the current immutable mesh object while the main thread prepares the state. Load integration proceeds before the before-training checkpoint joins the task. An explicit pause-before-training request still publishes its restart point immediately. Later phase checkpoints for that same live mesh reuse the prepared object without rescanning it. Geometry views are joined before refinement or destruction, and background failures propagate to the numerical controller.

Checkpoint payload commit, checksum verification, `latest` publication and journal publication retain their original order. A crash may leave an unused mesh object but cannot publish an uncommitted accepted state. Dense payload serialization remains synchronous; this optimization overlaps geometry work, not every checkpoint byte. `ALOD_CHECKPOINT_SYNC=1` restores synchronous geometry publication. Timers distinguish `checkpoint_geometry_background`, `checkpoint_geometry_wait`, `checkpoint_restart` and `checkpoint_accepted`; overlapping times must not be added as wall time.

## One-command version 4 delivery

The plot renderer now preserves the tuned parameters of
`generate_kernel_lift_figures.py` (E1/E2) and `plot_results.py` (E3).
E3 uses all accepted states through the first one-percent crossing and writes
`kappa_error_dof.pdf` and `.png`. E1/E2 retain their completed-H-step convention.
Rendering, table arithmetic, data import and packaging are separate modules;
see [the code structure review](code_structure_review.md).

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

## Validation scope

Bounded server tests exercise target crossing before the horizon, terminal resume, training-family targets, exact/full error agreement, zero auxiliary reference factors, cold/incremental parity, background/synchronous checkpoint parity and every applicable scientific diagnostic. A separate local-refinement test checks exact source/mass equality and error-jet reuse with changed numerical values. The historical delivery is also exercised against the supplied E1/E2/E3 exports. Full new-binary production campaigns and high-k late-state resource measurements remain separate work.

### Bounded performance evidence

The alternating E2 benchmark uses three accepted states, 16 training/audit members, 32 adaptive threads and 32 threads per asynchronous audit worker. It measures two repetitions after warmup. Cold mode disables incremental integration and background geometry preparation; both full-audit modes retain all earlier optimizations.

| Metric (seconds) | Cold, full audit | Incremental, full audit | Incremental, exact audit |
| --- | ---: | ---: | ---: |
| End-to-end wall time | 13.574 | 12.498 | 12.106 |
| Adaptive subprocess | 11.917 | 10.863 | 11.062 |
| Audit drain after adaptation | 1.404 | 1.381 | 0.793 |
| Load integration | 1.095 | 0.493 | 0.498 |
| Restart + accepted checkpoint blocking | 0.966 | 0.662 | 0.664 |

The optimized full-audit run reuses 49,082 source-element records. Background geometry work totals about 0.256 seconds; only 0.036 seconds remains in the join. The large exact-jet working set is bypassed, keeping target-check cost at 0.159 seconds rather than retaining an ineffective cache. Exact auditing reuses all three accepted solution blocks and performs zero auxiliary reference factorizations. Full auditing performs three. All three modes produce equivalent accepted marks/ranks and exact errors. These small-run medians establish bounded evidence, not universal production speedups. Reproduce with `tools/benchmark_scientific.py`; machine-independent evidence is retained in `docs/provenance/scientific_workflow_validation.json`.

Final validation passed all 52 server CTest cases (184.75 seconds) and all 29 local Python tests. The local-refinement test reduced source evaluations from 66,048 to 4,096, retained 62 exact-jet hits, and matched cold source moments and exact integrals bit for bit.
