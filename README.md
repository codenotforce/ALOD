# ALOD

ALOD is an independent C++20 migration of the E1/E2 Helmholtz experiments. It provides mixed-boundary P1 FEM, NVB meshes, AFEM/UFEM/SLOD baselines, adaptive two-sided LOD, explicit family marking, reference sweeps, regional AS/POD enrichment, dictionary inheritance and lazy localization checks.

For the current module boundaries and a file-by-file inventory, see [project structure](docs/project_structure.md). The [2026-09-28 source/documentation audit](docs/project_audit_20260928.md) records the documentation corrections and remaining historical-tool cleanup. Dated migration measurements remain evidence for their own binaries.

P0 source discovery and old-program replays are recorded in [the P0 evidence](docs/provenance/p0/README.md). P1-P6 implementation and bounded validation are summarized in [migration status](docs/migration_status.md). The [migration plan](ALOD_SUBPROJECT_AGENT_PLAN_20260909.md) preserves the chronological implementation and acceptance history; the latest status below supersedes early pending lists.

Version 4 now uses kernel-lifted enrichment tests and normalized lazy localization in the current paper presets. See [the version 4 review](docs/paper_v4_review.md) for the implemented changes, historical compatibility and remaining work. The September 27 production deployments completed E1/E2 at 28 cycles (57 accepted states each) and E3 at k=8,16,32,64,128, including their audits. These runs used frozen deployed sources. Later Theta and asynchronous changes have separate bounded/checkpoint validation; they have not all been rerun through the full campaign. Supplied historical controls remain labelled historical.

Current execution supports a manual lazy absolute-defect gate, concurrent audits through the same numerical executable, automatic checkpoint recovery, and unpinned runtime-selected threads. See [execution and recovery](docs/async_execution.md). The localization eigencheck uses cached operator images, parallel sparse products and deterministic Riesz gathers; see [operator mathematics and benchmark modes](docs/theta_fused_operator.md).

Current Linux runs also support bounded shared-memory audit snapshots, resident audit workers, exact cross-state local reuse and reduced-coordinate training. See [deep reuse and audit ownership](docs/deep_reuse.md) for defaults, fallback modes and validation.

## Build and test

Linux or WSL with GCC 11 or newer, CMake 3.20 or newer, Eigen 3.3 or newer, OpenMP, SuiteSparse/UMFPACK, BLAS and Python 3.10 or newer is required. The tested environment is Ubuntu 22.04. UMFPACK is a required dependency. Native Windows builds have not been validated.

```sh
sudo apt-get install cmake g++ libeigen3-dev libsuitesparse-dev python3 pkg-config
cmake --preset release
cmake --build --preset release
ctest --preset release
```

The build uses two compiler jobs and portable architecture flags. It does not require an old checkout, old static libraries, a server connection or Python third-party packages. Tests compare frozen old-kernel results, check the manufactured PDE/boundaries/quadrature and NVB invariants, and inspect repository portability.

## Run a baseline

```sh
python3 tools/run_baseline.py --config configs/smoke/e1_afem.json --output results/e1-afem-smoke
python3 tools/run_baseline.py --config configs/smoke/e2_ufem.json --output results/e2-ufem-smoke
python3 tools/run_baseline.py --config configs/smoke/e2_slod.json --output results/e2-slod-smoke
build/alod_run --help
build/alod_fem_smoke E1 6
```

The output directory must be new. `run.json` records the effective configuration, hashes, completion status and stop reason. `states.jsonl` records every solved state, errors, residuals, node counts, mesh fingerprints, element identities and marked sets. `stderr.log` records failures. These runtime files are ignored by Git. Thread settings are user-controlled; production launchers do not force single-thread BLAS or CPU affinity. The optional `emit_solution` field includes the complete nodal complex solution.

Each problem has AFEM, UFEM and SLOD presets under `configs/smoke/` and `configs/full/`. Full presets expose historical horizons and AFEM error targets; they can consume substantial resources and were not run as new production campaigns in P1. The schema, solver choices and baseline limitations are documented in [baselines](docs/baselines.md).

## Numerical scope

Run a P2 fixed state with explicit training-member IDs:

```sh
python3 tools/run_fixed.py --config configs/p2/e1_fixed.json --output results/e1-fixed
python3 tools/run_fixed.py --config configs/p2/e2_fixed.json --output results/e2-fixed
```

See [P2 modules and validation](docs/p2.md) for configuration, API contracts, server profiling and the separate manuscript/archived numerical policies.

E1 is the localized oscillatory packet with Dirichlet, Neumann and Robin boundaries. E2 is the weighted corner singularity plus an alpha=80 Gaussian wave on the L-shaped domain. The E1/E2 main presets use kappa=16; E3 uses the E1 problem at kappa=8 through 128. Formal load tables and imported canonical audits are in `data/rhs/` and `tests/fixtures/legacy/`. Baseline exact-error stopping is separate from strong-residual marking.

SLOD and LOD use two corrected spaces, direct Schur patch solves and an Eigen SparseLU coarse solve. New LOD defaults to manuscript area-weighted interpolation and geometric N² vertex patches. Explicit archive policies preserve the different averaging and patch expansion found in production sources; see [the implementation map](docs/paper_implementation_map.md). P1 SLOD remains restricted to uniform coarse meshes. Checkpoint restart and independent audits are described in [P5/P6](docs/p5_p6.md).

## Run an adaptive experiment

```sh
python3 tools/run_adaptive.py --config configs/adaptive/e1_smoke.json --output results/e1-adaptive
python3 tools/run_adaptive.py --config configs/adaptive/e2_smoke.json --output results/e2-adaptive
```

[P3/P4 documentation](docs/p3_p4.md) describes the state machine, all presets and deferred experiments. Current main presets use normalized `max(norm_k(U_nom), 1e-12) * Theta/eta_H > 0.3`; the old raw E2 threshold 0.1 is retained only for historical configurations. E1 fixed-ell nominal/family ALOD and AFEM controls use the same 48 members. The current 28-cycle E1/E2 main runs are complete on their recorded deployment. New-binary family/nominal controls and full lazy/every or inherit/reset comparison campaigns remain separate validation work. The E2 cold-start recurrence follows the plan and differs from the archive; its regression report distinguishes numerical proximity from equality.

## Resume and audit

```sh
python3 tools/run_adaptive.py --config configs/adaptive/e1_smoke.json --output results/e1-restart --pause-state 1
python3 tools/run_adaptive.py --config configs/adaptive/e1_smoke.json --output results/e1-restart --resume results/e1-restart/checkpoints/latest
python3 tools/run_audit.py --checkpoint results/e1-restart/checkpoints/latest --output results/e1-audit --batch-size 8
```

Current production presets retain restart checkpoints every five completed
cycles and at the final accepted state. Every accepted state still supplies an
immutable audit snapshot containing both meshes, the raw kernel/dictionary,
solutions, counters, warm starts and the frozen sample table. Periodic restart
points reuse those files without another serialization. The compatibility
setting `checkpoint_interval_cycles=0` retains the earlier every-boundary policy.
Audits run concurrently through the same executable without replaying the adaptive prefix.
Optional E2 `--fresh` and `--pure` diagnostics do not change the saved state.
Canonical CSV tables use external hashed marking sets; solver, validation,
audit and paper-data completion have separate status flags. See [P5/P6](docs/p5_p6.md)
for recovery rules, performance measurements and remaining P7/P8 scope.

The [integration optimization notes](docs/integration_optimization.md) describe
shared quadrature storage, joint exact-solution evaluation and reproducible
before/after benchmarks.
The [space and training optimization notes](docs/space_optimization.md) cover
LOD patch assembly, enrichment factors, temporary storage and phase benchmarks.

## Provenance and contribution checks

```sh
python3 tools/check_repository.py
python3 tools/check_repository.py --staged
```

The staged check inspects exact Git blobs. Documentation uses English, except for `ALOD_SUBPROJECT_AGENT_PLAN_20260909.md`, which remains in Chinese at the user's request. All committed files, including the plan, must contain no local absolute paths. Use relative paths or symbolic environment roots in documentation. Original source hashes remain in the provenance manifest even when an imported file is deliberately edited.

Rebuilding old-source references is an explicit developer action, never part of normal builds. See `tools/replay_legacy.py`, `tools/build_legacy_oracle.py` and `tools/build_baseline_oracles.py`. Historical fixtures must not be regenerated to hide a failed comparison. Source ownership and dependency licensing are described in [NOTICE](NOTICE.md).

Remote repository: [codenotforce/ALOD](https://github.com/codenotforce/ALOD).

The [September 12 four-wavenumber protocol](docs/kappa_campaign_20260912.md) is a historical recipe. For new ordered E1/E2/E3 deployments, use `tools/make_production_campaign.py` and [the current production protocol](docs/production_campaign_20260927.md).

See [workflow reuse and E1 timing](docs/workflow_reuse.md) for accepted-basis audit reuse, memory lifetime and the 32+32-thread server configuration.

See [five workflow optimizations](docs/five_optimizations.md) for bounded cross-state reuse, shared checkpoint geometry, and portable checkpoint export.

See [structural optimizations and remaining work](docs/structural_optimization.md) for shared state ownership, POD screening, early patch reuse, and the localization eigensolver.

See [audit scheduling, training updates, and global solve inventory](docs/async_training_optimization.md) for current runtime behavior.

Scientific workflow additions: [exact targets, lightweight audits, diagnostics, current-paper delivery and incremental I/O/integration](docs/scientific_workflow.md).

See [the code structure review and tuned plot conventions](docs/code_structure_review.md) for module boundaries and maintenance notes.

See [fine-grid stability calibration](docs/stability_calibration.md) for the
inf-sup diagnostic and experimental localization-threshold scaling.

See [the ordered production campaign](docs/production_campaign_20260927.md) for
the requested inverse-wavenumber thresholds, 28-cycle E1/E2 runs, 32+32 threads,
E3 dependency barriers and five-cycle restart policy.
