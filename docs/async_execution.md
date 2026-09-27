# Lazy absolute defect checks and asynchronous execution

## Localization policy

`ell_absolute_threshold` is the manually supplied value of tau(k) for the
configured wavenumber. Nonnegative finite values enable the gate; `-1` disables
it for archived configurations. Current production presets explicitly store
`3.2/k`, as requested: 0.4, 0.2, 0.1, 0.05 and 0.025 for k=8,16,32,64,128.
The runtime still consumes a manual scalar; it does not infer a scaling law.

At each scheduled lazy check, promote if either `Theta > tau(k)` or the existing
balance ratio exceeds its threshold. Equality does not promote. Both gates use
the same freshly computed Theta. After a promotion, rebuild and recheck on the
same mesh until both gates pass or `maximum_ell` is reached. The cap is reported
explicitly and never presented as satisfying the criteria. `ell_mode="fixed"`
continues to forbid promotions. Unscheduled states do not compute Theta.

This is the requested lazy variant, not an unconditional pre-training check on
every mesh. Changing tau(k) changes the mathematics and is rejected on resume.
An E3 configuration can set a different manual value for each k.

## One workflow, concurrent audits

```sh
python3 tools/run_adaptive.py --config configs/adaptive/e2_main.json --output results/e2-async --audit-workers 1 --audit-threads 0
python3 tools/run_adaptive.py --config configs/adaptive/e2_main.json --output results/e2-async --resume auto --audit-workers 1 --audit-threads 0
```

Set `audit=true` in the configuration. The supervisor discovers atomically
published accepted checkpoints while the adaptive solver continues. Audit
workers invoke the **same** `alod_run` executable with its `audit` subcommand.
The old `alod_audit` executable is only a compatibility wrapper; the integrated
workflow does not require it. The user does not launch a second audit program.

Workers are isolated processes, not threads sharing mutable meshes or solver
factors. A bounded worker pool controls simultaneous audit memory. Pending work
stays backed by immutable checkpoints on disk; it does not queue dense snapshots
in memory and does not block adaptive steps. The supervisor waits for outstanding
audits only after the numerical horizon ends. `run.json` then distinguishes
`auditing`, `complete`, `audit_failed` and numerical failure. The durable
`audit_queue.json` records pending/running/completed/failed jobs and timestamps.

Ordinary frozen E2 audits no longer construct regional Riesz patches/factors or
compute unused training seeds. Fresh retraining diagnostics still construct
them. Current accepted v3 checkpoints store the LOD trial basis and reduced
operator; ordinary audits restore both. They still assemble fine operators and
factor the reduced and reference FEM systems once per state. Numerical factors
are not shared between processes. Legacy v1/v2 inputs fall back to constructing
whichever cached objects are absent. See [five optimizations](five_optimizations.md)
for shared geometry storage and standalone checkpoint export.
See [workflow review](workflow_reuse.md) for lifetime, thread and timing details.

## Recovery contract

- With `checkpoint_interval_cycles=0`, training boundaries, ell promotions and
  accepted states retain the compatibility checkpoint policy. Current production
  presets use 5: only completed cycles divisible by five and the final accepted
  state become automatic restart points. Exact-target and state-limit endings
  count as final states. Explicit diagnostic pauses also retain a restart point.
  Partial sparse factorizations and individual greedy iterations are not serialized.
- Periodic-mode per-state audit transactions live in `audit_snapshots`; the
  automatic recovery search uses only `checkpoints`. Scheduled points hard-link
  their transaction and geometry into `checkpoints`, with a copy fallback on
  filesystems without hard links. A restart preserves superseded later snapshots
  in recovery folders before recomputing the interval. Audit coverage is not
  reduced to checkpoint states. Retained audit inputs still require disk I/O
  and storage; reducing restart frequency does not eliminate that cost.
- `--resume auto` searches newest transactions first, verifies the complete
  checksum, recovers a missing/stale pointer and repairs an uncommitted journal
  tail. A corrupt newest snapshot is retained in `checkpoints/recovery` while
  recovery uses the newest valid boundary. A mismatched committed journal fails
  closed instead of silently changing history.
- Mathematical settings and the frozen member table must match. Operational
  thread counts, audit concurrency and permitted resource/horizon increases do
  not redefine the mathematical checkpoint. Extending a previously terminal
  horizon retains the terminal check already performed; it is not claimed to be
  identical to a run that originally declared a longer horizon.
- Completed audits are reused only after checksum and sample-coverage checks.
  Interrupted, failed or invalid outputs are preserved; retries use new output
  directories. Resuming a completed numerical trajectory can therefore finish
  audits and exports without replaying the solver.
- An OS-backed exclusive lease prevents two supervisors writing one run. Linux
  numerical children terminate when their supervisor dies, preventing orphan
  writes after a hard kill. Graceful interruption cancels outstanding workers.
- The old 1 GiB file cap and two-million-row dense reader restriction are gone.
  Allocation checks still enforce file-backed lengths, finite values, valid
  indices, dimensions and checksums. Readers support v1, v2 and v3. Current shared snapshots require their
  referenced geometry objects; pack-checkpoint exports a self-contained file.
  Older binaries cannot read v3 checkpoints.
  Metadata inspection verifies the file without reconstructing its dense arrays.
- A startup interruption before the first checkpoint restarts at state zero
  only when the frozen configuration/member table match. In periodic mode,
  any preliminary journal and audit snapshots are preserved as recovery evidence
  before replay; compatibility mode requires an empty journal.
  Checkpointing does not make disk exhaustion or loss of all checkpoint files
  recoverable. Resource failures preserve the last completed transaction.

The campaign supervisor supports `--resume` for unfinished adaptive jobs and
retains memory/disk guards. It refuses duplicate live workers. Standalone
UFEM/SLOD baseline controllers do not have adaptive checkpoint continuation.

## Threads and affinity

`threads=0` inherits the OpenMP runtime; positive counts are explicit user
choices. The application no longer imposes a 64-thread ceiling. The current
paper presets use zero. `--audit-threads` can independently select audit threads;
omitting it uses the run's setting. The application no longer forces BLAS to one
thread, applies `taskset`, or calls `sched_setaffinity`. Launchers remove inherited
OpenMP affinity directives but cannot override scheduler/cgroup CPU restrictions.
Explicit user BLAS/OpenMP thread environment variables remain effective.

Regression tests and archived oracle builders still request fixed thread counts
for repeatability; those are test settings, not production limits. Parallelism
also depends on available work: an eight-RHS audit batch has at most eight outer
RHS tasks, and sparse factorization backends do not necessarily scale with the
OpenMP thread count. Removing a cap does not parallelize a serial algorithm.

## Bottlenecks and measurements

Use accepted-state phase timers for LOD/Riesz preparation, training, estimator,
Theta, load assembly and strong residuals. Audit completion records separate
preparation, reference factor/solve, coupled solves and error integration.
`solver_wall_seconds` and `audit_drain_seconds` describe the concurrent workflow;
overlapping inclusive times must not be added as if they were serial phases.

Important remaining costs are rebuilding correctors after mesh changes, sparse
reference factorizations, quadrature, and checkpoint transfer/checksum work.
Checkpoint construction still copies state arrays before serialization. The
latest-state serializer is synchronous; only audit computation is asynchronous.
Use measured phase and memory data before choosing further optimizations or
increasing audit concurrency. Large paper trajectories remain separate server
experiments; bounded tests do not establish their speedup or scalability.

## Server validation (2026-09-22)

[Evidence](provenance/async_validation.json) records 37 broad regression cases,
five final targeted cases and two campaign/tooling cases, all passing (38
unique CTest cases), including 23 Python unit tests. The asynchronous test holds
an audit worker until the final solver checkpoint appears, then checks equality
with the unaudited mathematical journal. It also covers startup recovery, an
actual supervisor kill, corruption quarantine, torn journals and audit retry.
The CLI accepts 128 threads; the server exposes 128 CPUs with no applied binding.

The opt-in large-file test successfully wrote, checksummed and reloaded a
1,552,584,656-byte checkpoint with 1,050,625 nodes and 48 solution columns. Peak
RSS was 1,911,996 KiB. Write plus verification took 12.43 seconds, metadata
inspection 7.32 seconds and full loading 8.55 seconds. This tests serialization
capacity, not a million-node Helmholtz solve.

```sh
cmake --build build --target stress_checkpoint
mkdir -p results
build/stress_checkpoint results/checkpoint-large
```

A three-state E2 H6/h10 probe with 16 training/audit members, eight explicitly
selected threads and one audit worker completed the solver in 6.24 seconds and
the full workflow in 10.09 seconds, with 3.73 seconds of final audit draining.
Two audits started before the numerical solver finished. No BLAS thread
environment variable or affinity was applied. Another server experiment was
present, so these are bounded observations, not an isolated speedup benchmark.

In this probe, summed strong-residual and load times were 2.44 and 2.22 seconds,
about 78% of accepted-state wall time (5.97 seconds). Theta took 0.27 seconds and
LOD/Riesz preparation 0.49 seconds. Audit error integration took 4.39 seconds
out of 7.35 seconds summed audit wall time. Quadrature dominates at this size;
large-grid factorization, memory and checkpoint scaling require separate data.

The subsequent [integration optimization](integration_optimization.md) reduced
audit error integration by 44% and complete workflow time by 18% on the same
three-state probe. Load and strong-residual times changed little at that size;
see the alternating-run measurements before extrapolating to larger jobs.


## Overlap inside an adaptive state

The current implementation starts load assembly while constructing the LOD space
and preparing kernel Riesz factors. Once both are ready, a scheduled localization
check can run concurrently with solution/training and the coarse estimator:

- Theta uses only the current immutable LOD basis, operators and kernel factors.
  Concurrent Riesz applications have private workspaces and atomic usage counters.
  Factors and patch geometry are not duplicated.
- The foreground still waits for both Theta and eta before deciding whether to
  promote ell. It joins the old task before replacing the space, publishing an
  accepted state or refining a mesh. The greedy training order is unchanged.
- While both computations are active, they divide the configured OpenMP budget.
  When the foreground finishes, subsequent background patch applications reclaim
  the full budget. A 32-thread adaptive run does not launch two 32-worker patch
  teams. Audit workers retain their separate configured budget.

`ALOD_ASYNC_DISABLE=1` selects sequential execution for controlled comparisons.
It is an operational switch, not a mathematical checkpoint parameter. Timing logs
include `load_wait`, `theta_wait`, `theta_denominator`, `theta_defect`,
`theta_cancellation_projection`, `theta_ritz`, `riesz_patch_apply` and
`riesz_scatter`. Threaded writes are serialized as complete JSONL records.
`theta_seconds` measures the full Theta task, rather than only the foreground's
wait. Concurrent phase durations must not be added to infer wall time.

## Bounded audit batch pipeline

Full audits overlap reference/frozen solves with quadrature for adjacent batches.
A persistent background lane performs load and error integration in order,
reusing its OpenMP team. The foreground alone accesses mutable reference/reduced
solver contexts. Numerical factors are constructed once per audited state.

The lane retains at most one running and one queued task. At most one next load
and one previous batch await consumption; outputs preserve sample order. Error
integration caches have a single writer. Ownership survives exceptions, and the
lane drains before borrowed mesh/operator objects are destroyed. Worker exceptions
propagate through futures and prevent an audit-complete record.

The foreground reserves one worker and the quadrature lane uses the remaining
workers (31 for a 32-thread audit). Load and integration deliberately share this
lane: competing quadrature teams caused a regression in the initial timing tests.
Integral work is partitioned using the actual task budget, avoiding unequal
chunk counts after a team is resized. The final drain returns the full team
to the lane. This controls ALOD's OpenMP
teams; separately threaded BLAS libraries still require an appropriate external
thread policy. Validation uses single-thread BLAS, as in existing regression tests.

The pipeline requires at least three threads and more than one batch. Fresh
training diagnostics and exact audits that directly reuse all accepted values
stay sequential. Additional retained dense buffers are checked against
`ALOD_ASYNC_BUFFER_BYTES`, default 1073741824 bytes. Above this budget the audit
uses its sequential path; reducing batch size or explicitly increasing this
budget can enable overlap on larger meshes. This is a retained-buffer budget,
not a cap on total process RSS or on existing sparse factors.

Regression coverage compares E1/E2 serial and concurrent trajectories, promotions,
resume, full/exact multi-batch audits, memory-budget fallback, concurrent Riesz
applications, worker reuse, exception propagation and queue draining. New binaries
are validated independently; already running production processes keep their
original executable and scheduling behavior.
