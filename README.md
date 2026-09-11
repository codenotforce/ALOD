# ALOD

ALOD is an independent C++20 migration of the E1/E2 Helmholtz experiments. It currently provides manufactured problems, mixed-boundary P1 FEM, conforming NVB meshes, a strong-residual AFEM baseline, uniform FEM, and uniform-grid SLOD with ell=3 and four additional reference levels.

P0 source discovery and old-program replays are recorded in [the P0 evidence](docs/provenance/p0/README.md). P1 implementation and validation are summarized in [migration status](docs/migration_status.md). The full [migration plan](ALOD_SUBPROJECT_AGENT_PLAN_20260909.md) also covers the remaining adaptive LOD, regional enrichment, checkpoint and production phases.

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

The output directory must be new. `run.json` records the effective configuration, hashes, completion status and stop reason. `states.jsonl` records every solved state, errors, residuals, node counts, mesh fingerprints, element identities and marked sets. `stderr.log` records failures. These runtime files are ignored by Git. BLAS uses one thread; the preset controls OpenMP threads. The optional `emit_solution` field includes the complete nodal complex solution.

Each problem has AFEM, UFEM and SLOD presets under `configs/smoke/` and `configs/full/`. Full presets expose historical horizons and AFEM error targets; they can consume substantial resources and were not run as new production campaigns in P1. The schema, solver choices and baseline limitations are documented in [baselines](docs/baselines.md).

## Numerical scope

E1 is the localized oscillatory packet with Dirichlet, Neumann and Robin boundaries. E2 is the weighted corner singularity plus an alpha=80 Gaussian wave on the L-shaped domain. Both use kappa=16. Formal load tables and imported canonical audits are in `data/rhs/` and `tests/fixtures/legacy/`. Baseline exact-error stopping is separate from strong-residual marking.

SLOD uses two corrected spaces, direct Schur patch solves and an Eigen SparseLU coarse solve. Adaptive quasi-interpolation has a known manuscript/source averaging discrepancy, reproduced by a regression test and documented in [the implementation map](docs/paper_implementation_map.md). P1 SLOD is restricted to uniform coarse meshes. General adaptive LOD/POD/AOT workflows and restart support are later phases.

## Provenance and contribution checks

```sh
python3 tools/check_repository.py
python3 tools/check_repository.py --staged
```

The staged check inspects exact Git blobs. Committed documents use English and contain no local absolute paths. Use relative paths or symbolic environment roots in documentation. Original source hashes remain in the provenance manifest even when an imported file is deliberately edited.

Rebuilding old-source references is an explicit developer action, never part of normal builds. See `tools/replay_legacy.py`, `tools/build_legacy_oracle.py` and `tools/build_baseline_oracles.py`. Historical fixtures must not be regenerated to hide a failed comparison. Source ownership and dependency licensing are described in [NOTICE](NOTICE.md).

Remote repository: [codenotforce/ALOD](https://github.com/codenotforce/ALOD).
