# Lazy absolute defect checks and asynchronous execution

## Localization policy

`ell_absolute_threshold` is the manually supplied value of tau(k) for the
configured wavenumber. Nonnegative finite values enable the gate; `-1` disables
it for archived configurations. Current paper presets explicitly set 0.2.
There is no inferred formula for tau(k).

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

- Training boundaries, every ell promotion and accepted states are automatically
  checkpointed, not just explicit test pauses. Work interrupted inside a solve
  restarts at the latest completed boundary; partial sparse factorizations and
  individual greedy iterations are not serialized.
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
  only when the frozen configuration/member table match and the numerical
  journal is empty. Checkpointing does not make disk exhaustion or loss of all checkpoint files
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
