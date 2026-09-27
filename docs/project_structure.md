# ALOD project structure and file catalogue

Audited on 2026-09-28 against the local working tree, including uncommitted and
new non-ignored files. All paths below are repository-relative. This guide maps
current code ownership and entry points; it does not turn historical evidence
into measurements of the latest executable. See [the audit findings](project_audit_20260928.md)
for resolved documentation findings, historical recipes and remaining code cleanup candidates.

## Start here

| Task | Entry point | Supporting guide |
|---|---|---|
| Build and run regression tests | `CMakeLists.txt`, `CMakePresets.json` | `README.md` build section |
| Launch current ordered E1/E2/E3 jobs | `tools/make_production_campaign.py`, then `tools/run_campaign.py` | `docs/production_campaign_20260927.md` |
| Run/resume one adaptive job | `tools/run_adaptive.py` -> `alod_run adaptive` | `docs/async_execution.md` |
| Run a baseline | `tools/run_baseline.py` -> `alod_run` | `docs/baselines.md` |
| Audit a saved state | `tools/run_audit.py` -> shared native audit implementation | `docs/scientific_workflow.md` |
| Inspect/pack checkpoints | `tools/checkpoint_io.py`, `alod_run` checkpoint subcommands | `docs/five_optimizations.md` |
| Measure final-state Theta | `benchmark_theta` optional target | `docs/theta_fused_operator.md` |
| Generate general paper package | `tools/paper_delivery.py` | `docs/scientific_workflow.md` |
| Reproduce requested 25-cycle plots and 28-cycle tail fits | `tools/paper_snapshot.py` | Output package `README.md` and `convergence.md` |
| Check repository integrity | `tools/check_repository.py` | `README.md` provenance section |

## Runtime and mathematical ownership

```text
configs + frozen RHS tables
  -> Python validation / run lease / campaign scheduling
  -> alod_run
       -> adaptive controller -> LOD hierarchy/correctors/reduced solves
                              -> regional training/kernel lifting/POD
                              -> Riesz indicators and lazy Theta checks
                              -> reference/coarse refinement
                              -> immutable audit snapshot + periodic checkpoint
       -> bounded shared audit or persistent disk-backed audit worker
            -> reference solve / batched loads / exact errors / optional diagnostics
  -> accepted events + audit samples + timing and resource logs
  -> paper data adapters -> selection -> figures / tables -> portable package
```

The numerical library is `alod_core`. `include/alod` provides the newer public
interfaces; lower-level `include/helmholtz`, `include/lod`, `include/mesh` and
matching source files are active extracted kernels. They are not a second
abandoned project. `apps` owns CLI/process boundaries; `tools` owns supervision,
portable result handling and rendering. Mathematical kernels do not belong in
Python plotting modules, and server job scheduling does not belong in them either.

`alod_audit` and `alod_run` share numerical audit source files. Current production
can serve concurrent audits through `alod_run`; the standalone frontend remains
useful for explicit offline requests. Theta's default uses cached Ritz images,
parallel global sparse products and deterministic Riesz gathers. The fused
patch-local action is an optional measured alternative. Legacy numerical policies
remain explicit compatibility paths, not current production defaults.

## Current policies versus historical records

Current main E1/E2 presets use 28 cycles, kernel-lifted enrichment tests, lazy
localization, a 0.3 solution-scaled relative threshold, tau(16)=0.2, and Ritz
relative tolerance 1e-9. The production generator uses tau(k)=3.2/k, 32 adaptive
threads and one 32-thread auditor, and the explicit dependency graph. E3 presets
have separate horizons and 1% exact-error stopping. These values are not implied
for every smoke, fixed-state, compatibility or family-control preset.

Periodic restart checkpoints and per-accepted-state audit snapshots are separate
contracts. Five-cycle restart retention does not mean only every fifth cycle has
an audit input. Snapshots and checkpoint geometry can be hard-linked; never sum
logical sizes as physical storage without inode deduplication. Removing snapshots
can prevent all-state re-audits even when formal restart checkpoints remain valid.

Recorded production completion and new optimization acceptance are separate:
completed E1/E2/E3 runs used frozen deployments; later Theta/async optimizations
have their own fixed-state and bounded regression evidence. An old document's
"not implemented" or "pending" paragraph must be read with its date and later
supersession notices. Frozen scientific fixtures must not be updated merely to
make a changed implementation pass.

## Top-level ownership

| Location | Ownership / retention |
|---|---|
| `apps/` | Executable frontends, adaptive/audit orchestration and transport |
| `include/`, `src/` | Numerical APIs, data ownership and implementations |
| `tools/` | Validation, runners, campaign tools, profiling and paper delivery |
| `configs/` | Explicit runnable presets and current paper inventory |
| `data/rhs/` | Frozen scientific sample definitions |
| `tests/` | C++/Python regressions, manual benchmarks and immutable fixtures |
| `docs/` | Current guides, dated design records and provenance evidence |
| `.github/` | Continuous integration workflow |
| `build/`, `build-v4/`, other ignored build trees | Generated binaries and CMake state; not authoritative source |
| `results/` | Ignored experiment evidence and self-contained result packages; potentially unique data |
| `_local/` | Ignored local scratch, transfer scripts and inspection outputs; not a maintained API |
| `.git/`, `__pycache__/` | Version-control internals and generated interpreter caches |

Ignored generated trees are documented by purpose rather than listing every
binary, transient log or downloaded copy. The maintained-file inventory below
is exhaustive for the audited working tree. A copied script inside a result
bundle belongs to that bundle's provenance, not to the live tool implementation.

## File catalogue

Status labels are not deletion instructions. "Historical / scoped" means the page preserves dated evidence and explicitly
distinguishes it from current operational policy.
"Manual historical tool" means intentional oracle/replay tooling outside normal
production. No whole numerical source file was proven abandoned by this audit.


Total maintained files: **342**, including this guide and the audit report.

### root (9 files)

| File | Role | Purpose |
|---|---|---|
| [.editorconfig](../.editorconfig) | Repository infrastructure | Shared whitespace, charset and line-ending editor conventions. |
| [.gitattributes](../.gitattributes) | Repository infrastructure | Git text normalization and repository line-ending policy. |
| [.gitignore](../.gitignore) | Repository infrastructure | Exclude generated builds, results, caches and machine-specific settings. |
| [ALOD_SUBPROJECT_AGENT_PLAN_20260909.md](../ALOD_SUBPROJECT_AGENT_PLAN_20260909.md) | Repository infrastructure | Chinese migration plan and chronological acceptance/change log; preserve history and distinguish later overrides. |
| [CMakeLists.txt](../CMakeLists.txt) | Repository infrastructure | Numerical library, executable targets, optional benchmarks, dependencies and CTest registration. |
| [CMakePresets.json](../CMakePresets.json) | Repository infrastructure | Release configure/build/test presets and compiler-job settings; not a numerical runtime thread cap. |
| [NOTICE.md](../NOTICE.md) | Repository infrastructure | Source ownership, inherited-code provenance and dependency licensing notices. |
| [README.md](../README.md) | Repository infrastructure | Project landing page, build/run examples, current production scope and documentation links. |
| [requirements-paper.txt](../requirements-paper.txt) | Repository infrastructure | Python dependencies for scientific figure/table generation. |

### .github (1 files)

| File | Role | Purpose |
|---|---|---|
| [.github/workflows/ci.yml](../.github/workflows/ci.yml) | Repository infrastructure | Ubuntu CI dependency installation, Release build and CTest execution. |

### apps (19 files)

| File | Role | Purpose |
|---|---|---|
| [apps/adaptive_run.cpp](../apps/adaptive_run.cpp) | Active / compatibility | Adaptive CLI, state loop, solution/training/check overlap, snapshot publication and shared-audit lifetime. |
| [apps/adaptive_run.hpp](../apps/adaptive_run.hpp) | Active / compatibility | Declaration of the adaptive CLI entry point used by alod_run. |
| [apps/alod_audit.cpp](../apps/alod_audit.cpp) | Active / compatibility | Small standalone audit executable frontend; delegates to the shared audit implementation. |
| [apps/alod_fixed.cpp](../apps/alod_fixed.cpp) | Active / compatibility | Fixed nested-mesh LOD/estimator/localization probe and structured numerical output. |
| [apps/alod_geometry.cpp](../apps/alod_geometry.cpp) | Active / compatibility | Initial mesh-resolution diagnostic; historical four-k mode and --paper-v4 five-k mode. |
| [apps/alod_profile.cpp](../apps/alod_profile.cpp) | Active / compatibility | Phase/allocation profiling entry for fixed-space, training and audit workloads. |
| [apps/alod_run.cpp](../apps/alod_run.cpp) | Active / compatibility | Unified command dispatch plus AFEM/UFEM/SLOD baseline trajectory execution. |
| [apps/alod_stability.cpp](../apps/alod_stability.cpp) | Active / compatibility | Optional global fine-grid inf-sup diagnostic and dense/sparse validation modes. |
| [apps/audit_diagnostics.cpp](../apps/audit_diagnostics.cpp) | Active / compatibility | Optional audit controls, fresh training, rank-zero and related diagnostic calculations. |
| [apps/audit_diagnostics.hpp](../apps/audit_diagnostics.hpp) | Active / compatibility | Declarations and data contracts for audit diagnostics. |
| [apps/audit_run.cpp](../apps/audit_run.cpp) | Active / compatibility | Immutable-checkpoint audit computation, batched RHS/error integration and scientific outputs. |
| [apps/audit_run.hpp](../apps/audit_run.hpp) | Active / compatibility | Audit options, reusable execution context and shared numerical audit API. |
| [apps/audit_worker.cpp](../apps/audit_worker.cpp) | Active / compatibility | Persistent audit request loop hosted by the unified executable. |
| [apps/fem_smoke.cpp](../apps/fem_smoke.cpp) | Active / compatibility | Small nominal FEM oracle probe used by migration regression tests. |
| [apps/fixed_support.hpp](../apps/fixed_support.hpp) | Active / compatibility | Shared fixed-state CLI validation, explicit member loading and structured-output helpers. |
| [apps/mesh_export.hpp](../apps/mesh_export.hpp) | Active / compatibility | JSON mesh export helpers for checkpoint/figure workflows. |
| [apps/shared_audit.cpp](../apps/shared_audit.cpp) | Active / compatibility | Bounded shared-memory audit transport, snapshot ownership, requests and fallback coordination. |
| [apps/shared_audit.hpp](../apps/shared_audit.hpp) | Active / compatibility | Shared audit endpoint/snapshot interfaces and ownership contracts. |
| [apps/supervised.hpp](../apps/supervised.hpp) | Active / compatibility | Linux parent-death handling to prevent orphan numerical workers after supervisor failure. |

### include (32 files)

| File | Role | Purpose |
|---|---|---|
| [include/alod/adaptive.hpp](../include/alod/adaptive.hpp) | Active / compatibility | Adaptive configuration, cursor/events, cycle and localization scheduling interfaces. |
| [include/alod/afem.hpp](../include/alod/afem.hpp) | Active / compatibility | Conforming FEM/AFEM solve, residual-estimator and marking interfaces. |
| [include/alod/batch.hpp](../include/alod/batch.hpp) | Active / compatibility | Shared batch quadrature, load/error integration and reusable audit/reference solve interfaces. |
| [include/alod/bordered_solve.hpp](../include/alod/bordered_solve.hpp) | Active / compatibility | Incremental bordered/Schur coupled-system solves with checked fallback. |
| [include/alod/checkpoint.hpp](../include/alod/checkpoint.hpp) | Active / compatibility | Checkpoint state/schema, shared geometry, atomic persistence and restore interfaces. |
| [include/alod/estimator.hpp](../include/alod/estimator.hpp) | Active / compatibility | Additive constrained-kernel Riesz applications, regional masks and family marking contracts. |
| [include/alod/execution.hpp](../include/alod/execution.hpp) | Active / compatibility | Scoped thread budgets, budget propagation, background tasks and bounded persistent task lane. |
| [include/alod/kernel_defect.hpp](../include/alod/kernel_defect.hpp) | Active / compatibility | Borrowed-defect operator API for parallel global and patch-local fused Theta actions. |
| [include/alod/legacy_case_id.hpp](../include/alod/legacy_case_id.hpp) | Active / compatibility | Stable identifier mapping for the retained legacy benchmark cases. |
| [include/alod/local_factor_cache.hpp](../include/alod/local_factor_cache.hpp) | Active / compatibility | Bounded exact-identity cache of immutable local numerical factors. |
| [include/alod/localization.hpp](../include/alod/localization.hpp) | Active / compatibility | Theta eigenproblem configuration, convergence records and warm-start/result interfaces. |
| [include/alod/lod.hpp](../include/alod/lod.hpp) | Active / compatibility | LOD hierarchy/space ownership, assembly, reduced solves, caches and resource limits. |
| [include/alod/mesh_state.hpp](../include/alod/mesh_state.hpp) | Active / compatibility | Adaptive mesh identity, lineage, marks and geometry state structures. |
| [include/alod/meshes.hpp](../include/alod/meshes.hpp) | Active / compatibility | Initial E1/E2 geometry factories and mesh utilities. |
| [include/alod/patch_cache.hpp](../include/alod/patch_cache.hpp) | Active / compatibility | Bounded local-corrector reuse with dependency identities and retained shared results. |
| [include/alod/problems.hpp](../include/alod/problems.hpp) | Active / compatibility | Current E1/E2 manufactured-problem factory and explicit family parameters. |
| [include/alod/regional.hpp](../include/alod/regional.hpp) | Active / compatibility | Regional kernel enrichment, training, POD, kernel lifting and archived adjoint-test APIs. |
| [include/alod/slod.hpp](../include/alod/slod.hpp) | Active / compatibility | Uniform-grid two-sided SLOD baseline interface. |
| [include/alod/timing.hpp](../include/alod/timing.hpp) | Active / compatibility | Phase timing, counters and thread-local timing destinations. |
| [include/helmholtz/benchmarks/paper_cases.h](../include/helmholtz/benchmarks/paper_cases.h) | Active / compatibility | Retained manufactured benchmark-case declarations used by current E1/E2 factories. |
| [include/helmholtz/boundary.h](../include/helmholtz/boundary.h) | Active / compatibility | Boundary classification and mixed-boundary data structures/functions. |
| [include/helmholtz/corrector.h](../include/helmholtz/corrector.h) | Active / compatibility | Localized corrector equation and corrected-basis assembly interfaces. |
| [include/helmholtz/manufactured.h](../include/helmholtz/manufactured.h) | Active / compatibility | Manufactured Helmholtz data/evaluation interfaces. |
| [include/helmholtz/operators.h](../include/helmholtz/operators.h) | Active / compatibility | Fine-grid Helmholtz/energy assembly, load integration and operator structures. |
| [include/helmholtz/patch_solver.h](../include/helmholtz/patch_solver.h) | Active / compatibility | Constrained patch linear-solver interfaces and factorization policies. |
| [include/helmholtz/patch_system.h](../include/helmholtz/patch_system.h) | Active / compatibility | Local patch assembly, constraints and reduced patch-system structures. |
| [include/helmholtz/quadrature.h](../include/helmholtz/quadrature.h) | Active / compatibility | Triangle/edge quadrature and prepared integration context interfaces. |
| [include/helmholtz/types.h](../include/helmholtz/types.h) | Active / compatibility | Shared complex scalar, matrix, function and Helmholtz parameter types. |
| [include/lod/patches.h](../include/lod/patches.h) | Active / compatibility | Geometric patch/support construction declarations. |
| [include/lod/quasi_interp.h](../include/lod/quasi_interp.h) | Active / compatibility | Coarse prolongation and area/arithmetic quasi-interpolation policies. |
| [include/mesh/refine.h](../include/mesh/refine.h) | Active / compatibility | NVB refinement, nested embeddings and ancestry output declarations. |
| [include/mesh/types.h](../include/mesh/types.h) | Active / compatibility | Triangular mesh, vertices, edges and boundary representation. |

### src (31 files)

| File | Role | Purpose |
|---|---|---|
| [src/adaptive/refine.cpp](../src/adaptive/refine.cpp) | Active / compatibility | Accepted adaptive mesh transactions, coarse/reference refinement and state transfer. |
| [src/adaptive/schedule.cpp](../src/adaptive/schedule.cpp) | Active / compatibility | Cycle/reference-sweep scheduling, lazy checks and promotion decisions. |
| [src/enrichment/aot.cpp](../src/enrichment/aot.cpp) | Active / compatibility | Cached global adjoint-test solves retained for compatibility/diagnostics; not the current kernel-lift default. |
| [src/enrichment/regional.cpp](../src/enrichment/regional.cpp) | Active / compatibility | Regional greedy family training, reduced coordinates, kernel lifting, POD and dictionary inheritance. |
| [src/estimator/kernel_defect.cpp](../src/estimator/kernel_defect.cpp) | Active / compatibility | Parallel D X / D* Y products and optional patch-local fused Gram actions. |
| [src/estimator/marking.cpp](../src/estimator/marking.cpp) | Active / compatibility | Explicit training-family normalization and mean/worst-member marking. |
| [src/estimator/riesz.cpp](../src/estimator/riesz.cpp) | Active / compatibility | Patch preparation/factor reuse and full/selected/lightweight kernel Riesz applications, including ordered gathers. |
| [src/estimator/riesz_internal.hpp](../src/estimator/riesz_internal.hpp) | Active / compatibility | Shared private Riesz storage and exact constrained solve used by estimator and defect paths. |
| [src/fem/afem.cpp](../src/fem/afem.cpp) | Active / compatibility | Conforming Helmholtz baseline solves, strong residuals and adaptive marking implementation. |
| [src/fem/batch.cpp](../src/fem/batch.cpp) | Active / compatibility | Batched reference solves, accepted-solution reuse, quadrature/load/error integration and incremental data reuse. |
| [src/fem/boundary.cpp](../src/fem/boundary.cpp) | Active / compatibility | Mixed Dirichlet/Neumann/Robin boundary treatment. |
| [src/fem/manufactured.cpp](../src/fem/manufactured.cpp) | Active / compatibility | Manufactured solution, derivatives and forcing support. |
| [src/fem/operators.cpp](../src/fem/operators.cpp) | Active / compatibility | Fine-grid sparse Helmholtz and energy matrices, loads and integration operations. |
| [src/fem/quadrature.cpp](../src/fem/quadrature.cpp) | Active / compatibility | Triangle/edge rules and reusable quadrature data construction. |
| [src/helmholtz/corrector.cpp](../src/helmholtz/corrector.cpp) | Active / compatibility | Localized corrector computation and corrected trial/test basis construction. |
| [src/helmholtz/patch_solver.cpp](../src/helmholtz/patch_solver.cpp) | Active / compatibility | Checked constrained patch solves, direct Schur/saddle and supported reference paths. |
| [src/helmholtz/patch_system.cpp](../src/helmholtz/patch_system.cpp) | Active / compatibility | Local matrix/constraint assembly and patch algebra representation. |
| [src/io/checkpoint.cpp](../src/io/checkpoint.cpp) | Active / compatibility | Versioned binary checkpoint/geometry formats, checksums, atomic writes, hard-link retention and readers. |
| [src/localization/dense.cpp](../src/localization/dense.cpp) | Active / compatibility | Dense generalized spectral reference/whitening routines for bounded cross-checks. |
| [src/localization/theta.cpp](../src/localization/theta.cpp) | Active / compatibility | Defect construction, projected Gram action, block Ritz iterations, image reuse and fresh residual acceptance. |
| [src/lod/fingerprint.hpp](../src/lod/fingerprint.hpp) | Active / compatibility | Internal stable numerical/geometry fingerprint construction. |
| [src/lod/lod.cpp](../src/lod/lod.cpp) | Active / compatibility | LOD hierarchy ownership, patch reuse, corrected basis/reduced operator assembly and solves. |
| [src/lod/patches.cpp](../src/lod/patches.cpp) | Active / compatibility | Patch geometry and support expansion implementation. |
| [src/lod/quasi_interp.cpp](../src/lod/quasi_interp.cpp) | Active / compatibility | Coarse P1 embeddings and manuscript/archived averaging implementations. |
| [src/lod/slod.cpp](../src/lod/slod.cpp) | Active / compatibility | Uniform two-sided SLOD baseline space and solution implementation. |
| [src/mesh/edges.cpp](../src/mesh/edges.cpp) | Active / compatibility | Edge connectivity and boundary edge construction. |
| [src/mesh/initial.cpp](../src/mesh/initial.cpp) | Active / compatibility | Retained E1 square and E2 L-shaped initial meshes. |
| [src/mesh/nested.cpp](../src/mesh/nested.cpp) | Active / compatibility | Nested mesh mappings, prolongations and fine/coarse relationship construction. |
| [src/mesh/refine.cpp](../src/mesh/refine.cpp) | Active / compatibility | NVB refinement, conformity closure and parent/child transfer data. |
| [src/problems/legacy_cases.cpp](../src/problems/legacy_cases.cpp) | Active / compatibility | Extracted E1/E2 exact-solution and forcing formulas; active numerical code despite the filename. |
| [src/problems/problems.cpp](../src/problems/problems.cpp) | Active / compatibility | Public current-problem/family factory wrapping the retained formulas. |

### tools (45 files)

| File | Role | Purpose |
|---|---|---|
| [tools/adaptive_config.py](../tools/adaptive_config.py) | Shared module | Compatibility fields shared by new runs and checkpoint validation. |
| [tools/async_audit.py](../tools/async_audit.py) | Shared module | Durable checkpoint-backed audit queue; workers never mutate solver state. |
| [tools/benchmark_async.py](../tools/benchmark_async.py) | Manual profiling | Compare bounded audit drain scheduling with numerical equivalence gates. |
| [tools/benchmark_deep_reuse.py](../tools/benchmark_deep_reuse.py) | Manual profiling | Alternating bounded E1/E2 comparisons against the prior optimized executable. |
| [tools/benchmark_integrations.py](../tools/benchmark_integrations.py) | Manual profiling | Alternating before/after E2 workflow benchmark with semantic output checks. |
| [tools/benchmark_reuse.py](../tools/benchmark_reuse.py) | Manual profiling | Alternating workflow benchmark with numerical and marking equivalence gates. |
| [tools/benchmark_scientific.py](../tools/benchmark_scientific.py) | Manual profiling | Bounded alternating cold/incremental and full/exact workflow measurements. |
| [tools/benchmark_spaces.py](../tools/benchmark_spaces.py) | Manual profiling | Alternate identical phase probes linked against before/after core libraries (Linux). |
| [tools/build_baseline_oracles.py](../tools/build_baseline_oracles.py) | Manual historical tool | Explicitly build a probe against a freshly built external legacy source tree. |
| [tools/build_legacy_oracle.py](../tools/build_legacy_oracle.py) | Manual historical tool | Build old kernels in an isolated directory and produce FEM oracle fixtures. |
| [tools/build_p2_oracles.py](../tools/build_p2_oracles.py) | Manual historical tool | Developer-only creation of new P2 fixtures from the archived numerical library. |
| [tools/campaign_job.py](../tools/campaign_job.py) | Current CLI/tool | Execute one isolated campaign job and its same-checkpoint diagnostics. |
| [tools/campaign_schedule.py](../tools/campaign_schedule.py) | Shared module | Validate dependency barriers independently of process/resource scheduling. |
| [tools/check_repository.py](../tools/check_repository.py) | Current CLI/tool | Check portable text, frozen source hashes and canonical audit completeness. |
| [tools/checkpoint_io.py](../tools/checkpoint_io.py) | Shared module | Validated checkpoint metadata, atomic files and accepted-journal recovery. |
| [tools/compare_adaptive.py](../tools/compare_adaptive.py) | Comparison tooling | Compare a bounded new trajectory with immutable archived P3/P4 evidence. |
| [tools/compare_baselines.py](../tools/compare_baselines.py) | Comparison tooling | Compare a complete small baseline trajectory with an old-kernel fixture. |
| [tools/compare_csv.py](../tools/compare_csv.py) | Comparison tooling | Compare selected numeric and identity columns of two archived CSV files. |
| [tools/compare_oracle.py](../tools/compare_oracle.py) | Comparison tooling | Compare a fresh nominal FEM solve with the frozen legacy-kernel output. |
| [tools/compare_p2.py](../tools/compare_p2.py) | Comparison tooling | Compare a P2 fixed state with frozen archived-kernel evidence. |
| [tools/execution.py](../tools/execution.py) | Shared module | Unpinned runtime configuration and an OS-released exclusive run lease. |
| [tools/export_results.py](../tools/export_results.py) | Current CLI/tool | Portable canonical tables with bounded CSV fields and hashed marking sets. |
| [tools/make_kappa_campaign.py](../tools/make_kappa_campaign.py) | Historical recipe — review A2 | Create the predeclared E1 kappa controls and revised E1/E2 production jobs. |
| [tools/make_paper_v4_campaign.py](../tools/make_paper_v4_campaign.py) | Alternate preparation | Prepare seven v4 adaptive jobs; this command never launches experiments. |
| [tools/make_production_campaign.py](../tools/make_production_campaign.py) | Current CLI/tool | Prepare the requested E1 -> E2 -> E3 barriers; never launch implicitly. |
| [tools/native_audit.py](../tools/native_audit.py) | Shared module | One bounded, cancellable native worker per active audit thread. |
| [tools/paper/__init__.py](../tools/paper/__init__.py) | Shared module | Version 4 paper delivery: data, selection, figures, tables and packaging. |
| [tools/paper/bundle.py](../tools/paper/bundle.py) | Shared module | Orchestrate validation, rendering and a self-contained delivery package. |
| [tools/paper/data.py](../tools/paper/data.py) | Shared module | Explicit native/historical adapters and portable scientific-data validation. |
| [tools/paper/figures.py](../tools/paper/figures.py) | Shared module | The user's tuned September plotting parameters, without host-specific paths. |
| [tools/paper/selection.py](../tools/paper/selection.py) | Shared module | Keep each manuscript's sampling convention explicit and shared with tables. |
| [tools/paper/tables.py](../tools/paper/tables.py) | Shared module | Scientific table calculations, independent of Matplotlib styling. |
| [tools/paper_delivery.py](../tools/paper_delivery.py) | Current CLI/tool | Generate the version 4 paper package from native, archive or canonical data. |
| [tools/paper_snapshot.py](../tools/paper_snapshot.py) | Specific 25/28-cycle export | Render current ALOD runs with existing controls and a separate full-run tail fit. |
| [tools/profile_fixed.py](../tools/profile_fixed.py) | Manual profiling | Run one fixed-state profile without changing CPU affinity. |
| [tools/replay_legacy.py](../tools/replay_legacy.py) | Manual historical tool | Rebuild an explicitly supplied archived production snapshot and replay P0. |
| [tools/run_adaptive.py](../tools/run_adaptive.py) | Current CLI/tool | Run or resume a checkpointed trajectory, with concurrent checkpoint-backed audits. |
| [tools/run_audit.py](../tools/run_audit.py) | Current CLI/tool | Audit immutable accepted checkpoints without replaying an adaptive prefix. |
| [tools/run_baseline.py](../tools/run_baseline.py) | Current CLI/tool | Validate an explicit P1 baseline preset and write a portable run artifact. |
| [tools/run_campaign.py](../tools/run_campaign.py) | Current CLI/tool | Linux tmux supervisor: restartable jobs and memory guards, no CPU binding. |
| [tools/run_fixed.py](../tools/run_fixed.py) | Current CLI/tool | Run one bounded P2 fixed state using the frozen explicit member table. |
| [tools/run_stability_calibration.py](../tools/run_stability_calibration.py) | Current CLI/tool | Measure fine-grid inf-sup constants and anchor tau(16)=0.2. |
| [tools/runtime_provenance.py](../tools/runtime_provenance.py) | Shared module | Path-free runtime and numerical build provenance for portable result bundles. |
| [tools/sample_run_resources.py](../tools/sample_run_resources.py) | Current CLI/tool | Sample an existing supervisor and its descendants without changing affinity. |
| [tools/scientific_diagnostics.py](../tools/scientific_diagnostics.py) | Current CLI/tool | Explicit offline scientific controls on immutable accepted checkpoints. |

### configs (36 files)

| File | Role | Purpose |
|---|---|---|
| [configs/adaptive/e1_family_afem.json](../configs/adaptive/e1_family_afem.json) | Explicit preset/inventory | E1 fixed-policy family-marked AFEM control; cycles=96, ell mode=fixed. Compatibility/comparison presets are not automatically current production settings. |
| [configs/adaptive/e1_family_alod.json](../configs/adaptive/e1_family_alod.json) | Explicit preset/inventory | E1 fixed-ell family-marked ALOD control; cycles=25, ell mode=fixed. Compatibility/comparison presets are not automatically current production settings. |
| [configs/adaptive/e1_main.json](../configs/adaptive/e1_main.json) | Production preset | Current E1 production ALOD preset; cycles=28, tau=0.2, five-cycle restart interval. |
| [configs/adaptive/e1_nominal_afem.json](../configs/adaptive/e1_nominal_afem.json) | Explicit preset/inventory | E1 nominal-marked AFEM control; cycles=96, ell mode=fixed. Compatibility/comparison presets are not automatically current production settings. |
| [configs/adaptive/e1_nominal_alod.json](../configs/adaptive/e1_nominal_alod.json) | Explicit preset/inventory | E1 fixed-ell nominal-marked ALOD control; cycles=25, ell mode=fixed. Compatibility/comparison presets are not automatically current production settings. |
| [configs/adaptive/e1_server32.json](../configs/adaptive/e1_server32.json) | Production preset | Explicit 32-thread E1 production preset; cycles=28, tau=0.2, five-cycle restart interval. |
| [configs/adaptive/e1_smoke.json](../configs/adaptive/e1_smoke.json) | Explicit preset/inventory | Bounded E1 controller smoke fixture preset; cycles=3, ell mode=lazy. Compatibility/comparison presets are not automatically current production settings. |
| [configs/adaptive/e2_every_inherit.json](../configs/adaptive/e2_every_inherit.json) | Explicit preset/inventory | E2 every-state-check/inherited-dictionary comparison preset; cycles=32, ell mode=every. Compatibility/comparison presets are not automatically current production settings. |
| [configs/adaptive/e2_every_reset.json](../configs/adaptive/e2_every_reset.json) | Explicit preset/inventory | E2 every-state-check/reset-dictionary comparison preset; cycles=32, ell mode=every. Compatibility/comparison presets are not automatically current production settings. |
| [configs/adaptive/e2_forced_smoke.json](../configs/adaptive/e2_forced_smoke.json) | Explicit preset/inventory | Forced same-state ell-promotion smoke preset; cycles=0, ell mode=lazy. Compatibility/comparison presets are not automatically current production settings. |
| [configs/adaptive/e2_lazy_inherit.json](../configs/adaptive/e2_lazy_inherit.json) | Explicit preset/inventory | E2 lazy-check/inherited-dictionary comparison preset; cycles=32, ell mode=lazy. Compatibility/comparison presets are not automatically current production settings. |
| [configs/adaptive/e2_lazy_reset.json](../configs/adaptive/e2_lazy_reset.json) | Explicit preset/inventory | E2 lazy-check/reset-dictionary comparison preset; cycles=32, ell mode=lazy. Compatibility/comparison presets are not automatically current production settings. |
| [configs/adaptive/e2_main.json](../configs/adaptive/e2_main.json) | Production preset | Current E2 kernel-lift production ALOD preset; cycles=28, tau=0.2, five-cycle restart interval. |
| [configs/adaptive/e2_server32.json](../configs/adaptive/e2_server32.json) | Production preset | Explicit 32-thread E2 production preset; cycles=28, tau=0.2, five-cycle restart interval. |
| [configs/adaptive/e2_smoke.json](../configs/adaptive/e2_smoke.json) | Explicit preset/inventory | Bounded E2 regional controller smoke preset; cycles=2, ell mode=lazy. Compatibility/comparison presets are not automatically current production settings. |
| [configs/adaptive/e3_k128.json](../configs/adaptive/e3_k128.json) | Production preset | Current E3 k=128 ALOD preset; cycles=100, tau=0.025, five-cycle restart interval. |
| [configs/adaptive/e3_k16.json](../configs/adaptive/e3_k16.json) | Production preset | Current E3 k=16 ALOD preset; cycles=25, tau=0.2, five-cycle restart interval. |
| [configs/adaptive/e3_k32.json](../configs/adaptive/e3_k32.json) | Production preset | Current E3 k=32 ALOD preset; cycles=22, tau=0.1, five-cycle restart interval. |
| [configs/adaptive/e3_k64.json](../configs/adaptive/e3_k64.json) | Production preset | Current E3 k=64 ALOD preset; cycles=20, tau=0.05, five-cycle restart interval. |
| [configs/adaptive/e3_k8.json](../configs/adaptive/e3_k8.json) | Production preset | Current E3 k=8 ALOD preset; cycles=27, tau=0.4, five-cycle restart interval. |
| [configs/full/e1_afem.json](../configs/full/e1_afem.json) | Explicit preset/inventory | Historical-horizon E1 AFEM baseline preset; explicit input to run_baseline.py. |
| [configs/full/e1_slod.json](../configs/full/e1_slod.json) | Explicit preset/inventory | Historical-horizon E1 SLOD baseline preset; explicit input to run_baseline.py. |
| [configs/full/e1_ufem.json](../configs/full/e1_ufem.json) | Explicit preset/inventory | Historical-horizon E1 UFEM baseline preset; explicit input to run_baseline.py. |
| [configs/full/e2_afem.json](../configs/full/e2_afem.json) | Explicit preset/inventory | Historical-horizon E2 AFEM baseline preset; explicit input to run_baseline.py. |
| [configs/full/e2_slod.json](../configs/full/e2_slod.json) | Explicit preset/inventory | Historical-horizon E2 SLOD baseline preset; explicit input to run_baseline.py. |
| [configs/full/e2_ufem.json](../configs/full/e2_ufem.json) | Explicit preset/inventory | Historical-horizon E2 UFEM baseline preset; explicit input to run_baseline.py. |
| [configs/p2/e1_fixed.json](../configs/p2/e1_fixed.json) | Explicit preset/inventory | Fixed-state e1 fixed LOD/estimator/localization probe configuration, not an adaptive campaign. |
| [configs/p2/e2_fixed.json](../configs/p2/e2_fixed.json) | Explicit preset/inventory | Fixed-state e2 fixed LOD/estimator/localization probe configuration, not an adaptive campaign. |
| [configs/p2/e2_medium.json](../configs/p2/e2_medium.json) | Explicit preset/inventory | Fixed-state e2 medium LOD/estimator/localization probe configuration, not an adaptive campaign. |
| [configs/paper_v4_inventory.json](../configs/paper_v4_inventory.json) | Explicit preset/inventory | Current version-4 five-figure/five-table inventory and manuscript identity. |
| [configs/smoke/e1_afem.json](../configs/smoke/e1_afem.json) | Explicit preset/inventory | Bounded smoke E1 AFEM baseline preset; explicit input to run_baseline.py. |
| [configs/smoke/e1_slod.json](../configs/smoke/e1_slod.json) | Explicit preset/inventory | Bounded smoke E1 SLOD baseline preset; explicit input to run_baseline.py. |
| [configs/smoke/e1_ufem.json](../configs/smoke/e1_ufem.json) | Explicit preset/inventory | Bounded smoke E1 UFEM baseline preset; explicit input to run_baseline.py. |
| [configs/smoke/e2_afem.json](../configs/smoke/e2_afem.json) | Explicit preset/inventory | Bounded smoke E2 AFEM baseline preset; explicit input to run_baseline.py. |
| [configs/smoke/e2_slod.json](../configs/smoke/e2_slod.json) | Explicit preset/inventory | Bounded smoke E2 SLOD baseline preset; explicit input to run_baseline.py. |
| [configs/smoke/e2_ufem.json](../configs/smoke/e2_ufem.json) | Explicit preset/inventory | Bounded smoke E2 UFEM baseline preset; explicit input to run_baseline.py. |

### data (2 files)

| File | Role | Purpose |
|---|---|---|
| [data/rhs/rhs_e1.json](../data/rhs/rhs_e1.json) | Frozen scientific input | E1 frozen formal RHS/sample identities, roles and manufactured-family parameters. |
| [data/rhs/rhs_e2.json](../data/rhs/rhs_e2.json) | Frozen scientific input | E2 frozen formal RHS/sample identities, roles and manufactured-family parameters. |

### tests (71 files)

| File | Role | Purpose |
|---|---|---|
| [tests/benchmark_reference_solvers.cpp](../tests/benchmark_reference_solvers.cpp) | Manual benchmark | Optional offline direct/iterative reference FEM comparison on one saved state; not a production backend. |
| [tests/benchmark_spaces.cpp](../tests/benchmark_spaces.cpp) | Manual benchmark | Explicit before/after phase benchmark for space construction and training. |
| [tests/benchmark_theta.cpp](../tests/benchmark_theta.cpp) | Manual benchmark | Same-checkpoint cold/warm Theta benchmark and global/fused/parallel action probe. |
| [tests/check_adaptive_cli.py](../tests/check_adaptive_cli.py) | Regression | Exercise complete/partial sweeps, fixed controls and same-state E2 promotion. |
| [tests/check_async_execution.py](../tests/check_async_execution.py) | Regression | Compare mathematical trajectories and multi-batch audits across scheduling modes. |
| [tests/check_async_restart.py](../tests/check_async_restart.py) | Regression | Actual asynchronous overlap, absolute lazy checks, and crash recovery. |
| [tests/check_campaign.py](../tests/check_campaign.py) | Regression | An unpinned campaign worker can run and resume without duplicate work. |
| [tests/check_campaign_dependencies.py](../tests/check_campaign_dependencies.py) | Regression | Dependency barriers include successful audit completion through job exit. |
| [tests/check_checkpoint_cli.py](../tests/check_checkpoint_cli.py) | Regression | Checkpoint recovery, immutable audits and numerical trajectory equivalence. |
| [tests/check_checkpoint_interval.py](../tests/check_checkpoint_interval.py) | Regression | Periodic restarts preserve all-state audits and replay superseded tails. |
| [tests/check_cli.py](../tests/check_cli.py) | Regression | Exercise preset execution, artifact status and preflight errors end to end. |
| [tests/check_fixed_cli.py](../tests/check_fixed_cli.py) | Regression | Exercise P2 artifact status, explicit IDs and failure reporting end to end. |
| [tests/check_large_training.py](../tests/check_large_training.py) | Regression | Physical 16-load E2 training where changing local ordering changes POD rank. |
| [tests/check_normalized_ell.py](../tests/check_normalized_ell.py) | Regression | Exercise initial checks, repeated promotions, cap reporting and resume. |
| [tests/check_optimization_workflow.py](../tests/check_optimization_workflow.py) | Regression | Optimized/reference E1, E2 and variable-k training and audit comparisons. |
| [tests/check_scientific_workflow.py](../tests/check_scientific_workflow.py) | Regression | Target stopping, resume, exact-only audits and frozen scientific controls. |
| [tests/check_shared_audit.py](../tests/check_shared_audit.py) | Regression | Shared immutable snapshots, disk spill, persistent fallback and numerical parity. |
| [tests/check_stability.py](../tests/check_stability.py) | Regression | Cross-check the inf-sup diagnostic against dense SVD on both domains. |
| [tests/check_v4_geometry.py](../tests/check_v4_geometry.py) | Regression | Check actual NVB diameters and configurable LOD reference refinement. |
| [tests/check_v4_initial.py](../tests/check_v4_initial.py) | Regression | Compare bounded initial states with the supplied E1/E2 production evidence. |
| [tests/check_wavenumber.py](../tests/check_wavenumber.py) | Regression | Variable-kappa operators, immutable audits, geometry and restart isolation. |
| [tests/check_workflow_reuse.py](../tests/check_workflow_reuse.py) | Regression | Accepted basis restoration, independent rebuild and runtime timing coverage. |
| [tests/fixtures/baselines/E1_AFEM.json](../tests/fixtures/baselines/E1_AFEM.json) | Evidence — retain | Frozen E1_AFEM small old-kernel baseline trajectory and numerical comparison oracle. |
| [tests/fixtures/baselines/E1_SLOD.json](../tests/fixtures/baselines/E1_SLOD.json) | Evidence — retain | Frozen E1_SLOD small old-kernel baseline trajectory and numerical comparison oracle. |
| [tests/fixtures/baselines/E1_UFEM.json](../tests/fixtures/baselines/E1_UFEM.json) | Evidence — retain | Frozen E1_UFEM small old-kernel baseline trajectory and numerical comparison oracle. |
| [tests/fixtures/baselines/E2_AFEM.json](../tests/fixtures/baselines/E2_AFEM.json) | Evidence — retain | Frozen E2_AFEM small old-kernel baseline trajectory and numerical comparison oracle. |
| [tests/fixtures/baselines/E2_SLOD.json](../tests/fixtures/baselines/E2_SLOD.json) | Evidence — retain | Frozen E2_SLOD small old-kernel baseline trajectory and numerical comparison oracle. |
| [tests/fixtures/baselines/E2_UFEM.json](../tests/fixtures/baselines/E2_UFEM.json) | Evidence — retain | Frozen E2_UFEM small old-kernel baseline trajectory and numerical comparison oracle. |
| [tests/fixtures/fem_E1.json](../tests/fixtures/fem_E1.json) | Evidence — retain | Frozen nominal fem_E1 FEM oracle, including solved numerical observables. |
| [tests/fixtures/fem_E2.json](../tests/fixtures/fem_E2.json) | Evidence — retain | Frozen nominal fem_E2 FEM oracle, including solved numerical observables. |
| [tests/fixtures/legacy/e1-family-afem.csv](../tests/fixtures/legacy/e1-family-afem.csv) | Evidence — retain | Immutable historical e1-family-afem accepted audit table; used by provenance/completeness checks, never rewritten to hide a regression. |
| [tests/fixtures/legacy/e1-family-alod.csv](../tests/fixtures/legacy/e1-family-alod.csv) | Evidence — retain | Immutable historical e1-family-alod accepted audit table; used by provenance/completeness checks, never rewritten to hide a regression. |
| [tests/fixtures/legacy/e1-nominal-afem.csv](../tests/fixtures/legacy/e1-nominal-afem.csv) | Evidence — retain | Immutable historical e1-nominal-afem accepted audit table; used by provenance/completeness checks, never rewritten to hide a regression. |
| [tests/fixtures/legacy/e1-nominal-alod.csv](../tests/fixtures/legacy/e1-nominal-alod.csv) | Evidence — retain | Immutable historical e1-nominal-alod accepted audit table; used by provenance/completeness checks, never rewritten to hide a regression. |
| [tests/fixtures/legacy/e2-accepted.csv](../tests/fixtures/legacy/e2-accepted.csv) | Evidence — retain | Immutable historical e2-accepted accepted audit table; used by provenance/completeness checks, never rewritten to hide a regression. |
| [tests/fixtures/p2/E1_graded.json](../tests/fixtures/p2/E1_graded.json) | Evidence — retain | Frozen E1_graded fixed-state old-library LOD/indicator/Theta oracle. |
| [tests/fixtures/p2/E1_ritz.json](../tests/fixtures/p2/E1_ritz.json) | Evidence — retain | Frozen E1_ritz fixed-state old-library LOD/indicator/Theta oracle. |
| [tests/fixtures/p2/E1_uniform.json](../tests/fixtures/p2/E1_uniform.json) | Evidence — retain | Frozen E1_uniform fixed-state old-library LOD/indicator/Theta oracle. |
| [tests/fixtures/p2/E2_graded.json](../tests/fixtures/p2/E2_graded.json) | Evidence — retain | Frozen E2_graded fixed-state old-library LOD/indicator/Theta oracle. |
| [tests/fixtures/p2/E2_ritz.json](../tests/fixtures/p2/E2_ritz.json) | Evidence — retain | Frozen E2_ritz fixed-state old-library LOD/indicator/Theta oracle. |
| [tests/fixtures/p2/E2_uniform.json](../tests/fixtures/p2/E2_uniform.json) | Evidence — retain | Frozen E2_uniform fixed-state old-library LOD/indicator/Theta oracle. |
| [tests/fixtures/p34/e1_prefix.json](../tests/fixtures/p34/e1_prefix.json) | Evidence — retain | Frozen e1_prefix adaptive prefix comparison evidence. |
| [tests/fixtures/p34/e2_initial.json](../tests/fixtures/p34/e2_initial.json) | Evidence — retain | Frozen e2_initial adaptive prefix comparison evidence. |
| [tests/fixtures/v4/initial_states.json](../tests/fixtures/v4/initial_states.json) | Evidence — retain | Supplied version-4 E1/E2 initial-state numerical evidence for bounded validation. |
| [tests/legacy_p2_probe.cpp](../tests/legacy_p2_probe.cpp) | Manual benchmark | Manual external-library fixed-state oracle probe built by tools/build_p2_oracles.py, not ordinary CMake. |
| [tests/quasi_average_probe.cpp](../tests/quasi_average_probe.cpp) | Regression | Independent graded-mesh counterexample for area versus arithmetic vertex averaging. |
| [tests/stress_checkpoint.cpp](../tests/stress_checkpoint.cpp) | Manual benchmark | Optional large checkpoint stress executable, excluded from ordinary builds. |
| [tests/test_adaptive_config.py](../tests/test_adaptive_config.py) | Regression | Unit tests for adaptive configuration validation and supported options. |
| [tests/test_async_execution.cpp](../tests/test_async_execution.cpp) | Regression | Scoped execution-budget, background context and task-lane regression. |
| [tests/test_baseline_tools.py](../tests/test_baseline_tools.py) | Regression | Unit tests for baseline validation, command construction and comparisons. |
| [tests/test_bordered_solve.cpp](../tests/test_bordered_solve.cpp) | Regression | Bordered-system updates, multi-RHS equivalence and unsafe-update fallback checks. |
| [tests/test_cross_state.cpp](../tests/test_cross_state.cpp) | Regression | Cross-state local reuse and identity invalidation regression. |
| [tests/test_execution.py](../tests/test_execution.py) | Regression | Unit tests for unpinned environments and execution/run-lease helpers. |
| [tests/test_fixed_config.py](../tests/test_fixed_config.py) | Regression | Unit tests for fixed-state parameters and explicit sample membership. |
| [tests/test_incremental_integration.cpp](../tests/test_incremental_integration.cpp) | Regression | Cross-refinement integration reuse versus fresh integration. |
| [tests/test_integration.cpp](../tests/test_integration.cpp) | Regression | Shared/fused integration equivalence and numerical checks. |
| [tests/test_kernel_defect.cpp](../tests/test_kernel_defect.cpp) | Regression | Fused/global actions, Hermitian behavior and ordered serial/parallel regional gather equivalence. |
| [tests/test_kernel_lift.cpp](../tests/test_kernel_lift.cpp) | Regression | Kernel-lifted enrichment test-space mathematical identities and coupled behavior. |
| [tests/test_mathematics.cpp](../tests/test_mathematics.cpp) | Regression | Independent FEM, boundary, quadrature, mesh and marking mathematical regression. |
| [tests/test_optimization_reuse.cpp](../tests/test_optimization_reuse.cpp) | Regression | Cache and prepared-product reuse versus fresh/reference numerical paths. |
| [tests/test_p2.cpp](../tests/test_p2.cpp) | Regression | Fixed-state LOD, constrained Riesz, interpolation, marking and spectral mathematics. |
| [tests/test_p34.cpp](../tests/test_p34.cpp) | Regression | Adaptive controller and regional-enrichment mathematical regression. |
| [tests/test_paper_delivery.py](../tests/test_paper_delivery.py) | Regression | Unit tests for paper data adapters, selection, calculations and delivery contracts. |
| [tests/test_patch_recovery.cpp](../tests/test_patch_recovery.cpp) | Regression | Patch reuse/recovery and invalidation regression. |
| [tests/test_relative_theta.cpp](../tests/test_relative_theta.cpp) | Regression | Relative Theta residual, dense comparisons, cold/warm starts and independent optimization switches. |
| [tests/test_result_io.py](../tests/test_result_io.py) | Regression | Long legacy fields and canonical marking-set integrity. |
| [tests/test_space_assembly.cpp](../tests/test_space_assembly.cpp) | Regression | Optimized space assembly versus reference assembly. |
| [tests/test_structural_reuse.cpp](../tests/test_structural_reuse.cpp) | Regression | Structural reuse and reduced training checks, including k=8 and k=128 variants. |
| [tests/test_tools.py](../tests/test_tools.py) | Regression | Unit tests for portability checks, frozen comparisons and audit completeness. |
| [tests/test_v4_config.py](../tests/test_v4_config.py) | Regression | Mathematical policy changes must never become implicit checkpoint resumes. |
| [tests/test_zero_kernel.cpp](../tests/test_zero_kernel.cpp) | Regression | Exactly trivial local kernels and zero-defect behavior. |

### docs (96 files)

| File | Role | Purpose |
|---|---|---|
| [docs/async_execution.md](../docs/async_execution.md) | Documentation | Lazy absolute defect checks and asynchronous execution; module/workflow documentation and its stated validation scope. |
| [docs/async_training_optimization.md](../docs/async_training_optimization.md) | Historical / scoped | Asynchronous audits, tighter localization, and bordered training solves; dated evidence with current-policy distinctions and follow-up links. |
| [docs/baselines.md](../docs/baselines.md) | Documentation | P1 baseline contract; module/workflow documentation and its stated validation scope. |
| [docs/code_structure_review.md](../docs/code_structure_review.md) | Historical / scoped | Code structure review and tuned figures; dated evidence with current-policy distinctions and follow-up links. |
| [docs/deep_reuse.md](../docs/deep_reuse.md) | Documentation | Cross-state reuse, coordinate training and shared asynchronous audits; module/workflow documentation and its stated validation scope. |
| [docs/e2_server32.md](../docs/e2_server32.md) | Historical / scoped | E2 server deployment (2026-09-22); dated evidence with current-policy distinctions and follow-up links. |
| [docs/experiments.md](../docs/experiments.md) | Historical / scoped | Experiment contracts; dated evidence with current-policy distinctions and follow-up links. |
| [docs/five_optimizations.md](../docs/five_optimizations.md) | Documentation | Five workflow optimizations; module/workflow documentation and its stated validation scope. |
| [docs/integration_optimization.md](../docs/integration_optimization.md) | Documentation | Integration optimization; module/workflow documentation and its stated validation scope. |
| [docs/kappa_campaign_20260912.md](../docs/kappa_campaign_20260912.md) | Historical / scoped | Revised policies and the E1 wavenumber campaign; dated evidence with current-policy distinctions and follow-up links. |
| [docs/migration_status.md](../docs/migration_status.md) | Historical / scoped | Migration status; dated evidence with current-policy distinctions and follow-up links. |
| [docs/p2.md](../docs/p2.md) | Historical / scoped | Fixed-state LOD and estimators; dated evidence with current-policy distinctions and follow-up links. |
| [docs/p3_p4.md](../docs/p3_p4.md) | Historical / scoped | P3/P4 adaptive controllers and regional enrichment; dated evidence with current-policy distinctions and follow-up links. |
| [docs/p5_p6.md](../docs/p5_p6.md) | Historical / scoped | P5/P6 checkpoints, independent audits and measured optimization; dated evidence with current-policy distinctions and follow-up links. |
| [docs/paper_implementation_map.md](../docs/paper_implementation_map.md) | Historical / scoped | Manuscript-to-implementation map; dated evidence with current-policy distinctions and follow-up links. |
| [docs/paper_v4_review.md](../docs/paper_v4_review.md) | Historical / scoped | Version 4 implementation review; dated evidence with current-policy distinctions and follow-up links. |
| [docs/production_campaign_20260927.md](../docs/production_campaign_20260927.md) | Documentation | Ordered E1/E2/E3 production campaign; module/workflow documentation and its stated validation scope. |
| [docs/project_audit_20260928.md](../docs/project_audit_20260928.md) | Documentation | Original documentation/source audit, completed documentation corrections and remaining historical-tool cleanup. |
| [docs/project_structure.md](../docs/project_structure.md) | Documentation | This architecture guide and exhaustive maintained-file purpose catalogue. |
| [docs/provenance/async_training_validation.json](../docs/provenance/async_training_validation.json) | Evidence — retain | Bounded audit-drain and bordered-training validation. |
| [docs/provenance/async_validation.json](../docs/provenance/async_validation.json) | Evidence — retain | Initial asynchronous execution/recovery validation. |
| [docs/provenance/baseline_contracts.json](../docs/provenance/baseline_contracts.json) | Evidence — retain | Frozen historical baseline settings and comparison contracts. |
| [docs/provenance/configs/e1-control-afem.json](../docs/provenance/configs/e1-control-afem.json) | Evidence — retain | Frozen historical e1-control-afem configuration; evidence, not a current production preset. |
| [docs/provenance/configs/e1-control-alod.json](../docs/provenance/configs/e1-control-alod.json) | Evidence — retain | Frozen historical e1-control-alod configuration; evidence, not a current production preset. |
| [docs/provenance/configs/e1-main.json](../docs/provenance/configs/e1-main.json) | Evidence — retain | Frozen historical e1-main configuration; evidence, not a current production preset. |
| [docs/provenance/configs/e2-main.json](../docs/provenance/configs/e2-main.json) | Evidence — retain | Frozen historical e2-main configuration; evidence, not a current production preset. |
| [docs/provenance/deep_reuse_validation.json](../docs/provenance/deep_reuse_validation.json) | Evidence — retain | Cross-state reuse, reduced training and shared-audit validation. |
| [docs/provenance/dependency_inventory.json](../docs/provenance/dependency_inventory.json) | Evidence — retain | Imported-source dependency closure and ownership inventory. |
| [docs/provenance/e1-all-controls-completed.json](../docs/provenance/e1-all-controls-completed.json) | Evidence — retain | Historical completion evidence for four E1 controls. |
| [docs/provenance/e1-baseline-reuse.json](../docs/provenance/e1-baseline-reuse.json) | Evidence — retain | Historical E1 baseline reuse provenance. |
| [docs/provenance/e1-build-provenance.json](../docs/provenance/e1-build-provenance.json) | Evidence — retain | Historical E1 compiler/link/binary/source identity record. |
| [docs/provenance/e1-family-afem-started.json](../docs/provenance/e1-family-afem-started.json) | Evidence — retain | Historical E1 family AFEM launch record. |
| [docs/provenance/e1-family-alod-started.json](../docs/provenance/e1-family-alod-started.json) | Evidence — retain | Historical E1 family ALOD launch record. |
| [docs/provenance/e1-nominal-afem-started.json](../docs/provenance/e1-nominal-afem-started.json) | Evidence — retain | Historical E1 nominal AFEM launch record. |
| [docs/provenance/e1-nominal-alod-started.json](../docs/provenance/e1-nominal-alod-started.json) | Evidence — retain | Historical E1 nominal ALOD launch record. |
| [docs/provenance/e1-paper-summary.json](../docs/provenance/e1-paper-summary.json) | Evidence — retain | Historical E1 paper result summary. |
| [docs/provenance/e1-portable-data-validation.json](../docs/provenance/e1-portable-data-validation.json) | Evidence — retain | Historical E1 exported-data validation. |
| [docs/provenance/e1-server-environment.json](../docs/provenance/e1-server-environment.json) | Evidence — retain | Historical server/compiler/runtime environment. |
| [docs/provenance/e2-baseline-reuse.json](../docs/provenance/e2-baseline-reuse.json) | Evidence — retain | Historical E2 baseline reuse provenance. |
| [docs/provenance/e2-paper-summary.json](../docs/provenance/e2-paper-summary.json) | Evidence — retain | Historical E2 accepted paper result summary. |
| [docs/provenance/evidence_verification.json](../docs/provenance/evidence_verification.json) | Evidence — retain | Original supplied evidence inventory verification. |
| [docs/provenance/five_optimizations_validation.json](../docs/provenance/five_optimizations_validation.json) | Evidence — retain | Five workflow optimizations: numerical and timing comparisons. |
| [docs/provenance/integration_validation.json](../docs/provenance/integration_validation.json) | Evidence — retain | Load/residual/error integration optimization comparison. |
| [docs/provenance/p0/dependencies.json](../docs/provenance/p0/dependencies.json) | Evidence — retain | P0 dependency versions and reconstruction inputs. |
| [docs/provenance/p0/e1-branch-inventory.json](../docs/provenance/p0/e1-branch-inventory.json) | Evidence — retain | E1: Literal historical source feature/environment branch inventory. |
| [docs/provenance/p0/e1-build-commands.json](../docs/provenance/p0/e1-build-commands.json) | Evidence — retain | E1: Recorded historical compiler and linker commands. |
| [docs/provenance/p0/e1-coarse_mark_recipe.csv](../docs/provenance/p0/e1-coarse_mark_recipe.csv) | Evidence — retain | E1: Archived coarse marking recipe and selected identities. |
| [docs/provenance/p0/e1-lazy_global_ell.csv](../docs/provenance/p0/e1-lazy_global_ell.csv) | Evidence — retain | E1: Archived E1 lazy localization decision trace. |
| [docs/provenance/p0/e1-link-closure.json](../docs/provenance/p0/e1-link-closure.json) | Evidence — retain | E1: Linked object/library dependency closure. |
| [docs/provenance/p0/e1-reference_policy.csv](../docs/provenance/p0/e1-reference_policy.csv) | Evidence — retain | E1: Archived E1 reference-refinement policy trace. |
| [docs/provenance/p0/e1-replay-comparison.json](../docs/provenance/p0/e1-replay-comparison.json) | Evidence — retain | E1: Rebuilt old-source replay versus supplied evidence comparison. |
| [docs/provenance/p0/e1-replay-config.json](../docs/provenance/p0/e1-replay-config.json) | Evidence — retain | E1: Explicit archived replay configuration. |
| [docs/provenance/p0/e1-resource-patch.json](../docs/provenance/p0/e1-resource-patch.json) | Evidence — retain | E1: Resource-only changes applied when replaying old code. |
| [docs/provenance/p0/e1-sc_lod_rhs_family.csv](../docs/provenance/p0/e1-sc_lod_rhs_family.csv) | Evidence — retain | E1: Archived E1 per-state family solution/error evidence. |
| [docs/provenance/p0/e1-source-hashes.json](../docs/provenance/p0/e1-source-hashes.json) | Evidence — retain | E1: Hashes of the exact historical sources used for replay. |
| [docs/provenance/p0/e2-accepted-prefix.csv](../docs/provenance/p0/e2-accepted-prefix.csv) | Evidence — retain | E2: Accepted E2 prefix identity/state records. |
| [docs/provenance/p0/e2-accepted_control_samples.csv](../docs/provenance/p0/e2-accepted_control_samples.csv) | Evidence — retain | E2: Accepted E2 nominal/pure control-sample evidence. |
| [docs/provenance/p0/e2-accepted_workflow_training.csv](../docs/provenance/p0/e2-accepted_workflow_training.csv) | Evidence — retain | E2: Accepted E2 training workflow evidence. |
| [docs/provenance/p0/e2-branch-inventory.json](../docs/provenance/p0/e2-branch-inventory.json) | Evidence — retain | E2: Literal historical source feature/environment branch inventory. |
| [docs/provenance/p0/e2-build-commands.json](../docs/provenance/p0/e2-build-commands.json) | Evidence — retain | E2: Recorded historical compiler and linker commands. |
| [docs/provenance/p0/e2-coarse_mark_recipe.csv](../docs/provenance/p0/e2-coarse_mark_recipe.csv) | Evidence — retain | E2: Archived coarse marking recipe and selected identities. |
| [docs/provenance/p0/e2-ell_events.csv](../docs/provenance/p0/e2-ell_events.csv) | Evidence — retain | E2: Archived E2 localization/promotion events. |
| [docs/provenance/p0/e2-link-closure.json](../docs/provenance/p0/e2-link-closure.json) | Evidence — retain | E2: Linked object/library dependency closure. |
| [docs/provenance/p0/e2-replay-comparison.json](../docs/provenance/p0/e2-replay-comparison.json) | Evidence — retain | E2: Rebuilt old-source replay versus supplied evidence comparison. |
| [docs/provenance/p0/e2-replay-config.json](../docs/provenance/p0/e2-replay-config.json) | Evidence — retain | E2: Explicit archived replay configuration. |
| [docs/provenance/p0/e2-resource-patch.json](../docs/provenance/p0/e2-resource-patch.json) | Evidence — retain | E2: Resource-only changes applied when replaying old code. |
| [docs/provenance/p0/e2-source-hashes.json](../docs/provenance/p0/e2-source-hashes.json) | Evidence — retain | E2: Hashes of the exact historical sources used for replay. |
| [docs/provenance/p0/quasi-average-counterexample.json](../docs/provenance/p0/quasi-average-counterexample.json) | Evidence — retain | Independent graded averaging discrepancy evidence. |
| [docs/provenance/p0/README.md](../docs/provenance/p0/README.md) | Evidence — retain | P0 reconstruction, source provenance and comparison instructions. |
| [docs/provenance/p1_floating_point_sensitivity.json](../docs/provenance/p1_floating_point_sensitivity.json) | Evidence — retain | Architecture/roundoff sensitivity in near-tied adaptive marks. |
| [docs/provenance/p1_validation.json](../docs/provenance/p1_validation.json) | Evidence — retain | P1 core, baseline and portability acceptance record. |
| [docs/provenance/p2_validation.json](../docs/provenance/p2_validation.json) | Evidence — retain | P2 LOD/estimator/localization and fixed-state acceptance record. |
| [docs/provenance/p34_validation.json](../docs/provenance/p34_validation.json) | Evidence — retain | P3/P4 adaptive and regional mathematical/CLI acceptance record. |
| [docs/provenance/p56_profile.json](../docs/provenance/p56_profile.json) | Evidence — retain | P5/P6 fixed-state resource and timing profile. |
| [docs/provenance/p56_validation.json](../docs/provenance/p56_validation.json) | Evidence — retain | Checkpoint/audit and measured P6 optimization acceptance record. |
| [docs/provenance/paper_inventory.json](../docs/provenance/paper_inventory.json) | Evidence — retain | Frozen earlier manuscript figure/table inventory; not the current v4 inventory. |
| [docs/provenance/plan_evidence.json](../docs/provenance/plan_evidence.json) | Evidence — retain | Migration-plan evidence and deliverable mapping. |
| [docs/provenance/plot_style_validation.json](../docs/provenance/plot_style_validation.json) | Evidence — retain | Historical supplied-script pixel/style comparison; predates later E3 sampling revision. |
| [docs/provenance/repeated_work_validation.json](../docs/provenance/repeated_work_validation.json) | Evidence — retain | Repeated-construction removal benchmark and equivalence record. |
| [docs/provenance/scientific_workflow_validation.json](../docs/provenance/scientific_workflow_validation.json) | Evidence — retain | Exact targets, audit modes, diagnostics and incremental workflow validation. |
| [docs/provenance/source_differences.json](../docs/provenance/source_differences.json) | Evidence — retain | Recorded differences between supplied source variants. |
| [docs/provenance/source_manifest.json](../docs/provenance/source_manifest.json) | Evidence — retain | Original imported hashes plus current normalized content hashes; enforced by repository integrity checks. |
| [docs/provenance/space_validation.json](../docs/provenance/space_validation.json) | Evidence — retain | Space/training construction optimization comparison. |
| [docs/provenance/structural_validation.json](../docs/provenance/structural_validation.json) | Evidence — retain | Structural ownership/cache/training optimization comparison. |
| [docs/provenance/v4_review.json](../docs/provenance/v4_review.json) | Evidence — retain | Machine-readable version-4 manuscript/source review. |
| [docs/provenance/v4_validation.json](../docs/provenance/v4_validation.json) | Evidence — retain | Kernel-lift and normalized-localization validation. |
| [docs/provenance/validation.json](../docs/provenance/validation.json) | Evidence — retain | Initial repository/source validation record. |
| [docs/reference_solver_benchmark.md](../docs/reference_solver_benchmark.md) | Documentation | Reference FEM solver comparison; module/workflow documentation and its stated validation scope. |
| [docs/repeated_work_review.md](../docs/repeated_work_review.md) | Historical / scoped | Repeated-work review after the E1 timing run; dated evidence with current-policy distinctions and follow-up links. |
| [docs/scientific_workflow.md](../docs/scientific_workflow.md) | Documentation | Scientific workflows and current-paper delivery; module/workflow documentation and its stated validation scope. |
| [docs/space_optimization.md](../docs/space_optimization.md) | Documentation | LOD, enrichment and training optimization; module/workflow documentation and its stated validation scope. |
| [docs/stability_calibration.md](../docs/stability_calibration.md) | Documentation | Fine-grid stability and localization-threshold calibration; module/workflow documentation and its stated validation scope. |
| [docs/structural_optimization.md](../docs/structural_optimization.md) | Historical / scoped | Structural optimization review and implementation; dated evidence with current-policy distinctions and follow-up links. |
| [docs/theta_fused_operator.md](../docs/theta_fused_operator.md) | Documentation | Parallel and patch-local localization operators; module/workflow documentation and its stated validation scope. |
| [docs/theta_optimization.md](../docs/theta_optimization.md) | Documentation | Localization eigenproblem operator reuse; module/workflow documentation and its stated validation scope. |
| [docs/workflow_reuse.md](../docs/workflow_reuse.md) | Historical / scoped | Workflow reuse and E1 server deployment; dated evidence with current-policy distinctions and follow-up links. |
