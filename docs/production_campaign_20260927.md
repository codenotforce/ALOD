# Ordered E1/E2/E3 production campaign

## Recorded completion and source scope

The September 27 deployments completed E1 and E2 at 28 cycles (57 accepted
states each), followed by E3 k=8,16,32,64,128 and their audits. E3 recorded
54,50,41,36,30 accepted states respectively and reached the nominal 1% target.
E1's successful result belongs to the original campaign; successful E2/E3
results belong to the replacement campaign. The earlier failed E2 attempt is
not the completed E2 result.

These are measurements of frozen deployment sources. Later asynchronous and
Theta optimizations have separate bounded or same-checkpoint validation, not a
full rerun of every production trajectory. Current source availability does not
change which executable produced an existing result. The protocol below is the
current preparation/retention policy, not an assertion that all binaries in old
server directories are synchronized with the working tree.

## Production policy

The numerical gamma calibration was cancelled at the user's request. Production
uses the chosen rule tau(k)=0.2*(16/k)=3.2/k, including k=8. This is a prescribed
policy, not a claim that the cancelled calibration established that law.

| k | Absolute defect threshold |
|---:|---:|
| 8 | 0.4 |
| 16 | 0.2 |
| 32 | 0.1 |
| 64 | 0.05 |
| 128 | 0.025 |

The existing solution-scaled relative threshold 0.3 and lazy check schedule
remain active. Either absolute or relative failure can request promotion.

The revised dependency order is E1, E2, E3 k=8/16 together, then k=32,
k=64 and k=128 sequentially. Use `--start-at-e2` to prepare a fresh replacement
campaign while retaining the completed E1 and the original E2 outputs.
A dependency completes only after the adaptive solver, all audits and exports
succeed. A failed prerequisite blocks its dependents. Resume retries unfinished
jobs and retains successful jobs. Memory/disk guards can defer admission or
stop a worker; they never silently reduce thread counts or mark a failed job
complete. The concurrent k=8/16 jobs request up to 128 compute threads in total.

The reference FEM auditor uses UMFPACK's 64-bit index interface and METIS
ordering for large systems. This avoids the 32-bit factor workspace limit and
reduces fill; the existing relative residual acceptance gate remains unchanged.

Every job uses 32 adaptive threads and one 32-thread asynchronous auditor. Drain
also uses one auditor; the default extra drain worker is explicitly disabled.
CPU affinity and inherited OpenMP thread limits are not imposed by the launcher.

E1 and E2 run 28 cycles with two reference sweeps per cycle: 57 accepted states,
including state zero. Exact-target stopping remains disabled for these two runs.
E3 retains nominal exact-error stopping at 1%, exact-only audits, and its existing
bounded horizons:

| k | Cycles | State cap | Maximum accepted states |
|---:|---:|---:|---:|
| 8 | 27 | 55 | 55 |
| 16 | 25 | 51 | 51 |
| 32 | 22 | 45 | 45 |
| 64 | 20 | 38 | 38 |
| 128 | 100 | 30 | 30 |

Reaching a horizon/resource limit without reaching 1% is not target completion.

## Restart and audit storage

`checkpoint_interval_cycles=5` retains automatic restart points only at cycle
5,10,15,... and at normal termination, including an early exact target or a
partial-cycle state cap. E1/E2 therefore retain states 10,20,30,40,50,56. No
initial checkpoint is added unless the initial state itself is terminal.
Explicit diagnostic pause requests are an exception and save their requested
state/phase. A crash may require replay from the preceding five-cycle boundary.

All accepted states still have audit inputs in `audit_snapshots`. They are
separate from the automatic restart index and remain available for offline
diagnostics. Scheduled restart points link the immutable snapshot and geometry
into `checkpoints`; no second serialization is needed. On replay, later audit
inputs are quarantined before publication, preventing stale-audit reuse and
immutable-name collisions. Before the first scheduled checkpoint, replay starts
at state zero and preserves the superseded journal and audit inputs.

Reducing restart frequency suppresses intermediate training/ell restart writes
and non-audited state writes. It does not remove the per-state disk snapshots
needed by the asynchronous audit fallback. Snapshot and restart timing labels
are separate. Logical byte counts for these directories may share physical
storage through hard links.

## Prepare, launch and resume

```sh
python3 tools/make_production_campaign.py --output results/production-campaign
python3 tools/run_campaign.py --campaign results/production-campaign --build build
python3 tools/run_campaign.py --campaign results/production-campaign --build build --resume
```

Use tmux for the long-running supervisor. The generated campaign holds its own
frozen configurations; edits to repository presets do not mutate running jobs.
Per-job `run.json`, `solver.jsonl`, `timings.jsonl` and audit timing files retain
the mathematical settings, accepted states, wall time and stage timings.
The campaign also records process-group memory and CPU usage.
