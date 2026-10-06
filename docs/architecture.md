# Architecture

This guide describes maintained module ownership, not every historical file.
See [algorithm](algorithm.md) for mathematical contracts and
[runtime](runtime.md) for snapshot ownership and failure boundaries.

## Entry points

| Task | Python entry point | Numerical entry point |
|---|---|---|
| AFEM/UFEM/SLOD baseline | `tools/run_baseline.py` | `alod_run E1|E2 AFEM|UFEM|SLOD` |
| Fixed nested-pair study | `tools/run_fixed.py` | `alod_fixed` |
| Adaptive run/resume | `tools/run_adaptive.py` | `alod_run adaptive E1|E2` |
| Audit accepted input | `tools/run_audit.py` | Shared audit implementation in `alod_run` |
| Ordered campaign | `tools/make_production_campaign.py`, `tools/run_campaign.py` | Per-job processes |
| Scientific diagnostics | `tools/scientific_diagnostics.py` | Audit/reconstruction/diagnostic commands |
| Paper delivery | `tools/paper_delivery.py`, `tools/paper_snapshot.py` | Saved metrics and optional mesh export |
| Fine-grid stability | `tools/run_stability_calibration.py` | `alod_stability` |

`alod_audit` is a compatibility wrapper. Launching the C++ adaptive CLI alone
does not run the Python audit scheduler or post-audit state reclamation.

## Module ownership

| Location | Responsibility |
|---|---|
| `apps/adaptive_run.cpp` | Controller, accepted state, within-state overlap and publication |
| `apps/audit_run.*`, `apps/shared_audit.*` | Numerical audit context, shared service and immutable state transfer |
| `apps/adaptive_*.hpp`, `apps/fixed_support.hpp` | CLI support, state/member identities and application helpers |
| `src/mesh/`, `include/alod/mesh.hpp` | NVB refinement, boundary tags, stable identities and ancestry |
| `src/fem/`, `src/helmholtz/`, `src/problems/` | P1 operators, quadrature, strong residuals and RHS families |
| `src/lod/`, `include/alod/lod.hpp` | Interpolation, local correctors, two-sided spaces and reduced PG solves |
| `src/estimator/`, `include/alod/estimator.hpp` | Constrained local Riesz solves, family marking and defect actions |
| `src/localization/` | Dense or matrix-free localization eigenproblem |
| `src/enrichment/`, `include/alod/regional.hpp` | Regional training, kernel lifting, POD and inheritance |
| `src/adaptive/` | Reference-sweep schedule and nested mesh transitions |
| `src/io/`, `include/alod/checkpoint.hpp` | Versioned state/geometry formats, validation and durable writes |
| `include/alod/execution.hpp`, `include/alod/timing.hpp` | Task budgets, completion guards and phase/team observations |
| `tools/checkpoint_io.py`, `tools/async_audit.py`, `tools/state_retention.py` | Journal recovery, audit queue/receipts and state cleanup |
| `tools/runtime_provenance.py`, `tools/check_repository.py` | Runtime/build identity and repository integrity |
| `tools/paper/` | Data selection, importing, plotting, tables and packaging |

## Data flow

Validated configuration and a frozen RHS table define an experiment. The
numerical controller constructs the nested hierarchy and spaces, trains the
permitted RHS cases, verifies an accepted state and durably publishes it.
The supervisor schedules read-only audits and exports metrics; audit outputs
never drive marking or localization decisions. Retention can remove numerical
inputs only after validated durable completion. Resume validates retained
state and its committed mathematical journal.

Local reuse requires exact operator/constraint/geometry/RHS identities, not
global numbering or a hash alone. Optional source moments and analytic jets
reuse unchanged physical cells under the frozen family/quadrature identity;
discrete solution values are reevaluated. Restart reconstructs optional caches
from durable state. Reduced-coordinate training still verifies terminal
results through the original estimator. Cache controls are in configuration.

## Source, configuration and generated data

| Location | Meaning |
|---|---|
| `configs/smoke/`, `configs/full/` | Baseline smoke and historical full horizons |
| `configs/p2/`, `configs/adaptive/` | Fixed-state/adaptive presets; mathematical policy must be read explicitly |
| `data/rhs/` | Formal E1/E2 RHS identities and roles |
| `docs/provenance/`, `tests/fixtures/` | Source/replay/benchmark evidence and immutable numerical controls |
| `tests/test_*.py`, `tests/check_*.py`, C++ tests | Tooling, end-to-end, mathematical and concurrency checks |
| `build*/` | Ignored binaries and build state |
| `results/`, `out/`, `artifacts/` | Ignored run outputs and exports, potentially unique evidence |
| `_local/` | Ignored local scratch and transfer helpers |

Core numerical dependencies are Eigen, OpenMP, SuiteSparse/UMFPACK and BLAS.
Dependencies are discovered, not vendored. Legacy builders/replay tools are
explicit developer actions; normal builds do not need an external old checkout.
