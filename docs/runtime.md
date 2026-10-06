# Runtime, audits and recovery

See [configuration](configuration.md) for settings and [algorithm](algorithm.md)
for adaptive-state semantics. This guide owns runtime behavior and storage.

## Run artifacts and completion

An adaptive output contains frozen `members.txt`/configuration, `run.json`,
mathematical `solver.jsonl`, derived `events.jsonl`, timing/error logs,
checkpoint/audit inputs, audit outputs and `audit_queue.json`.
`solver_completed`, `validation_passed`, `audit_complete` and `paper_complete`
are separate flags. An accepted state or horizon stop is not automatically
target attainment, audit completion or full paper delivery.

Standalone baselines stream `states.jsonl`, `stderr.log`, `run.json` and mesh
exports; fixed probes produce `state.json`. Fresh outputs must be new.
Baseline thread metadata records requested environment, not measured teams.

## Concurrent work and ownership

Load assembly overlaps LOD/Riesz preparation. Scheduled Theta overlaps
solution/training and estimation. Completion guards return a task's budget to
its surviving partner after either completion or exception. Later parallel
regions adopt it; a running OpenMP team cannot resize. Both futures are joined
before borrowed state changes or is destroyed. Mathematical decisions wait
for the required Theta and estimator.

Audits consume every accepted audited state while adaptation continues. Linux
ALOD runs normally use bounded resident threads in the numerical process with
immutable shared hierarchy/state ownership and audit-owned mutable solver
contexts. At most two unclaimed snapshots are retained within the configured
allowance; active snapshots also count. Evicted/oversized inputs use disk.
Durable input integrity is checked even when memory supplies numerical arrays.
Adaptation does not wait for an audit slot.

The fallback uses persistent separate workers, with optional one-shot mode.
Ordinary audit exceptions do not change adaptive decisions. A native crash or
process OOM in an integrated auditor affects the shared numerical process;
disk-only execution supplies additional process isolation. Worker-owned reference
FEM factors are not reused across changed meshes.

`--audit-workers` controls active concurrency. Drain concurrency starts after
numerical completion; the integrated process can remain alive to serve jobs.
Default drain expansion reuses a known explicit released thread budget;
unknown/inherited counts do not trigger expansion. Ordered production explicitly
keeps one drain worker. Concurrent state audits can finish out of order;
RHS solves sharing one mutable factor are serialized. Audit batches overlap
solves with bounded quadrature tasks and preserve output order.

## Threads and measurements

Positive `threads` is an adaptive OpenMP budget; auditors have separate budgets.
There is no application 64-thread ceiling, hardware sizing or CPU reserve mode.
Launchers clear OpenMP affinity directives and, for explicit thread counts,
an inherited OMP_THREAD_LIMIT. They cannot override scheduler/cgroup CPU limits.
BLAS variables are preserved. BLAS does not obey a universal process-wide cap
merely because an OpenMP budget was selected.

Element-partitioned batched quadrature can use more workers than RHS cases;
some per-case paths, long regions and serial factorization still limit use.
`team` timing events contain actual `workers` and `requested_workers`.
`budget_return` contains an allocation, not an observed team.
Startup team probes and CPU percentage do not prove continuous useful work.
Nested/overlapping phase times must not be added as disjoint wall time.

## Serialized state and publication

Audit snapshots and checkpoints share serialization but have different roles:
accepted states are auditable; before-training/after-ell boundaries can be
resumable without being audit inputs. They contain meshes, boundary/stable
identities, ancestry, embeddings, raw kernel/enrichment basis, computed values
and IDs, warm blocks, ell/cursor/marks, frozen RHS/configuration and journal
consistency metadata. V3 also stores a trial basis and base reduced operator.
They do not persist every fine-grid matrix or LU factor. New processes rebuild
required operators/factors; the shared memory path can share extra live objects.

Readers support v1/v2/v3 and validate sizes, indices, finite values, ordering,
dimensions and checksums. The old fixed 1 GiB cap is gone. V3 may reference
shared geometry objects; older binaries cannot read v3. Export a self-contained
copy when transferring one state:

```sh
build/alod_run pack-checkpoint "$CHECKPOINT" results/portable-state.bin
```

Background geometry preparation joins before refinement/destruction. Dense
payload serialization and durable commit remain synchronous. Temporary files
are verified/synchronized before atomic publication; latest and committed
journal rows follow the defined order. Partial factorizations/greedy iterations
are not serialized. Linux durability is implemented; native Windows durability
has not been validated.

Positive `checkpoint_interval_cycles` publishes scheduled completed-cycle restart
points plus final accepted state; every accepted audited state still publishes
an input. Scheduled points hard-link inputs/geometry into `checkpoints`, with
copy fallback. Interval zero preserves boundary publication. Diagnostic pauses
publish their requested restart boundary. Publication frequency is separate
from final retention.

| Retention | Successful-audit input lifetime |
|---|---|
| `all` (default) | Keep every published input |
| `delete_after_audit` | Remove successful inputs; final retained only with `keep_final=true` |
| `periodic` | Keep positive cycle boundaries divisible by the interval and always keep the final accepted state |

Both deletion modes require audit=true; periodic requires a positive interval.
Final includes early exact target and partial-cycle cap; state zero is retained
by periodic only if terminal. `keep_final` defaults false and only controls
delete_after_audit. The Python runner performs cleanup across both hard links.

## Validated cleanup

The state file is published and verified before audit execution. On success,
the runner validates the audit results, synchronizes the result files, and
atomically commits a completion receipt in `audit_queue.json`. Only then may it
remove the corresponding input. The receipt binds the checkpoint hash, output
hash, RHS-case coverage and accepted journal row. It is separate from the large
numerical payload.

Reclamation also waits for the accepted row to appear completely in
`solver.jsonl` and synchronizes its journal prefix before removal. Pending and
failed audit inputs remain available for retries.
The current `checkpoints/latest` target is temporarily protected while the
producer runs and while the supervisor uses it to recover the committed journal.
At audit drain completion, it is repointed to a surviving restart point or
removed before its old target is deleted. Thus strict deletion can finish with
no state files and no `latest`. Explicit pause checkpoints remain protected so
the requested pause can be resumed.

Training/ell restart boundaries do not themselves enter the audit queue. In the
reclamation modes they are removed only after the accepted state of that same
step has passed audit; an interrupted step's restart boundaries remain intact.
Shared mesh objects are collected after the producer has stopped, when no
background geometry publication can be in flight. Objects referenced by
retained, pending, failed or quarantined state files are preserved. Unreadable
retained metadata defers geometry collection. Recovery evidence is not deleted.

## Recovery and limitations

Completed audits can be recovered from their validated receipts even after
input deletion. Missing or corrupt results for a reclaimed state fail closed;
the runner cannot silently invent an audit or mark it complete. Rollback/replay
invalidates receipts whose accepted journal rows differ. Retained inputs permit
ordinary retries of incomplete or invalid audits.

In strict `delete_after_audit`, no restart point may remain. `--resume auto` then
replays from the beginning, preserving superseded journal/state evidence, rather
than reconstructing numerical state from audit metrics. With `keep_final=true`,
the final checkpoint remains resumable. Periodic retention permits replay from
the most recent retained cycle boundary. The original committed journal and RHS
table are still required for trajectory resume.

These modes reduce long-term storage, not every per-state disk write. Pending
inputs accumulate if audits run more slowly than adaptation; this change does
not introduce an admission budget. `retention.removed_logical_bytes` counts
unlinked paths and can double-count hard-linked data; it is not a physical disk
reclamation measurement.

## Resume contract

```sh
python3 tools/run_adaptive.py --config configs/adaptive/e1_smoke.json --output results/pause --pause-state 1
python3 tools/run_adaptive.py --config configs/adaptive/e1_smoke.json --output results/pause --resume auto
```

Resume needs the original run directory, committed journal and frozen RHS table.
Automatic recovery verifies newest transactions, repairs stale pointers and
uncommitted tails, and quarantines corrupt/newer superseded inputs before replay.
If interruption precedes the first scheduled boundary, replay starts at zero
under matching configuration/table. A changed committed journal fails closed.
OS-backed exclusive leases prevent duplicate supervisors; Linux numerical
children terminate with their supervisor. Graceful interruption cancels workers.

Mathematical settings must match. Operational threads, audits, retention and
permitted resource/horizon increases can change. Extending a terminal horizon
keeps the terminal check already performed and is not claimed equal to a run
that originally declared that longer horizon. Pause at an unchanged horizon
supports uninterrupted-run comparison. Resource failure preserves completed
transactions; loss of all inputs cannot be repaired from error metrics.
