# Workflow reuse and E1 server deployment

This page records the initial 32+32-thread E1 deployment. Current defaults also
include the [five follow-up optimizations](five_optimizations.md), including
v3 checkpoint geometry sharing and reduced-operator reuse.

## Repeated work removed

Accepted checkpoints now include the sparse LOD trial basis. Audits validate and
restore this basis, skipping all local corrector solves. Checkpoint writing
borrows the live basis rather than first copying it into the checkpoint object.
Legacy checkpoints, changed ell, and explicit `--rebuild-lod` audits reconstruct
the basis. The latter provides an independent numerical comparison.

Frozen E2 audits reuse the accepted enrichment dictionary and construct the
coupled evaluator once per state, not once per sample batch. They do not retrain
unless fresh-training diagnostics are explicitly requested. E1 has rank zero,
so there is no enrichment training to reuse. Operators, interpolation, reduced
factors and reference FEM factors are still reconstructed once per independent
audit process/state; mesh changes prevent indiscriminate reuse across states.

Mesh and dictionary fingerprints are computed once per audit state. Per-sample
solution fingerprints remain necessary. Accepted bases increase checkpoint disk
usage in exchange for avoiding repeated local solves.

## Memory lifetime

Adaptive LOD spaces, Riesz caches, loads, solution workspaces and training
scratch are released before mesh transition. Required inherited kernels and
warm-start state remain alive. Each LOD worker releases its last patch-factor
cache after construction. Audits release checkpoint hierarchy matrices, the
serialized basis, warm starts and raw kernels after restoration. The live fine
mesh and energy matrix are referenced rather than copied; temporary reference
residuals are released before error integration. E1 no longer duplicates its
solution solely for zero-enrichment diagnostics.

## Threads

RHS-only partitioning previously limited 16-member training loads to 16 workers
and eight-member audit batches to eight. Small batches now partition elements
as well as members. Deterministic element-order accumulation preserves the
numerical reduction order. Tests observe all 32 workers even with two members.

Explicit thread requests disable OpenMP dynamic teams and remove inherited
OMP_THREAD_LIMIT and affinity overrides. There is no fixed 16-thread production
cap. Task-count bounds remain where fewer independent tasks exist. Sparse
factorization and serial assembly phases do not necessarily occupy 32 cores.
Runtime timing logs record the observed OpenMP team and reject a smaller team.

## Timing and deployment

`configs/adaptive/e1_server32.json` requests the complete E1 horizon: k=16,
51 accepted states, 16 marking members and 48 audited members. The deployment
uses 32 adaptive threads and one asynchronous audit worker with 32 threads.
The numerical solver continues while checkpoint-backed audits are pending.

Run from the deployment directory:

```sh
python3 tools/run_adaptive.py --config configs/adaptive/e1_server32.json --output results/e1-main --executable build/alod_run --audit-workers 1 --audit-threads 32 --timeout 604800
```

The server deployment directory is `e1-32-20260922` under the existing project
root; its tmux session is `alod-e1-32-20260922`. The launch wrapper and source
snapshot are retained with the deployment. No existing experiment is replaced.

Timing evidence is kept separately from mathematical journals:

- `results/e1-main/timings.jsonl`: adaptive phases, including hierarchy/operators,
  local correctors, basis assembly, reduced assembly/factorization, loads,
  training/solve, Riesz preparation/estimation, Theta, strong residuals,
  checkpoint I/O and mesh transitions.
- `results/e1-main/audits/*/timings.jsonl`: checkpoint loading, basis restoration,
  operator/factor preparation, reference and frozen solves, loads and errors.
- `results/e1-main/audit_queue.json`: audit scheduling and completion status.
- `results/e1-resources.jsonl`: process-tree CPU ticks, active threads, resident
  memory, swap and allowed CPU lists, sampled every five seconds.
- `results/e1-supervisor.log`: supervisor output; `results/e1-exit-code` appears
  when the supervised run exits.

Phase timers are nested wall-clock intervals; do not sum overlapping parent
and child intervals. Resource sampling observes the processes without binding
cores or changing their thread settings. A running experiment is not a completed
paper result. Detailed timing analysis is deferred until the user requests it.

## Validation

The Release server build passes all 41 CTest tests. Restored and independently
rebuilt E1/E2 audits match exactly on the numerical sample records. A separate
three-state E1 rehearsal with all 48 members completed with 32 adaptive threads
and 32 audit threads; every audit reported basis reuse and zero patch rebuilds.
