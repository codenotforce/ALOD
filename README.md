# ALOD

ALOD is a C++20 implementation of adaptive two-sided localized orthogonal
decomposition for the E1/E2 Helmholtz problems, with AFEM/UFEM/SLOD baselines,
regional enrichment, family marking, numerical audits and checkpoint recovery.
E3 uses the E1 family at five wavenumbers.

## Guides

| Guide | Purpose |
|---|---|
| [Architecture](docs/architecture.md) | Module ownership, entry points, data flow and dependencies |
| [Algorithm](docs/algorithm.md) | Current mathematics, adaptive state, estimators and paper/archive differences |
| [Configuration](docs/configuration.md) | Configuration fields, presets, operational switches and defaults |
| [Runtime](docs/runtime.md) | Concurrent audits, thread budgets, checkpoint publication, retention and recovery |
| [Experiments](docs/experiments.md) | Current E1/E2/E3 protocol, frozen results and historical comparisons |
| [Scientific workflow](docs/scientific_workflow.md) | Diagnostics, error metrics, figures, tables and portable delivery |
| [Validation](docs/validation.md) | Tests, provenance rules, benchmarks and evidence limits |

These seven guides and this README are the eight maintained user documents.
Optimization mechanisms are summarized where they affect mathematics or use;
the former per-optimization and dated review pages have been consolidated.

## Build and test

Linux or WSL requires GCC 11+, CMake 3.20+, Eigen 3.3+, OpenMP,
SuiteSparse/UMFPACK, BLAS and Python 3.10+. Ubuntu 22.04 was tested;
native Windows numerical builds have not been validated.

```sh
sudo apt-get install cmake g++ libeigen3-dev libsuitesparse-dev python3 pkg-config
cmake --preset release
cmake --build --preset release
ctest --preset release
```

The release preset uses two compiler jobs and portable architecture flags.
UMFPACK is required; CHOLMOD is optional for the reference-solver benchmark.
Core run tools require no Python third-party packages. Paper rendering uses
[requirements-paper.txt](requirements-paper.txt).

## Run

```sh
python3 tools/run_baseline.py --config configs/smoke/e1_afem.json --output results/e1-baseline
python3 tools/run_fixed.py --config configs/p2/e2_fixed.json --output results/e2-fixed
python3 tools/run_adaptive.py --config configs/adaptive/e1_smoke.json --output results/e1-adaptive
python3 tools/run_adaptive.py --config configs/adaptive/e1_smoke.json --output results/e1-adaptive --resume auto
python3 tools/run_audit.py --checkpoint results/e1-adaptive/checkpoints/latest --output results/e1-audit
```

Fresh output directories must be new; resume uses the original adaptive run.
`threads=0` inherits the OpenMP runtime; a positive count is explicit.
There is no hardware-based automatic thread-sizing or CPU-reserve mode.
See the runtime guide before choosing checkpoint cleanup or audit concurrency.

## Evidence and instructions

The September 27 frozen deployments completed E1/E2 at 28 cycles and 57
accepted states each, plus E3 at k=8,16,32,64,128. Later optimizations have
their own bounded/checkpoint validation; they were not all rerun through the
full campaign. Source/configuration identity determines which results apply.

The [P0 provenance guide](docs/provenance/p0/README.md) preserves historical
acceptance context. The original migration plan is operator-local, ignored by
Git and not required by a checkout. Current guides supersede early pending lists.
Machine-readable evidence
is under `docs/provenance/`; frozen numerical fixtures are under `tests/fixtures/`.
Ownership and dependencies are described in [NOTICE](NOTICE.md).

[AGENTS](AGENTS.md) and the [terminology skill](.agents/skills/alod-terminology/SKILL.md)
govern agent work. Its [mathematical notation](.agents/skills/alod-terminology/references/mathematical-notation.md)
and [implementation correspondence](.agents/skills/alod-terminology/references/implementation-correspondence.md)
are explanatory references. Project guides/JSON use English; agent guidance
may be bilingual. Use relative paths or symbolic
environment roots in maintained text.

Remote repository: [codenotforce/ALOD](https://github.com/codenotforce/ALOD).
