# Cross-state reuse, coordinate training and shared asynchronous audits

This update supersedes the pending items in the earlier structural and asynchronous optimization notes. The comparison baseline is commit `566d98b`. Mathematical policies, greedy block order, localization tolerance, integration rules and audit sample coverage are unchanged.

## Local correctors and Riesz factors

LOD cache identities now contain the actual local unknown geometry, the restricted Helmholtz operator, restricted interpolation constraints and the target RHS contributions. Global node/element numbers and irrelevant entries outside the local restriction are excluded. Coordinates are relative to the target triangle, so translated patches are also reusable when their discrete operators, constraints and target contributions match exactly. Interpolation dependency rows are ordered once per hierarchy, rather than repeatedly for each patch. Exact byte comparison establishes equality; a hash alone never authorizes reuse. Local solve ordering remains unchanged for historical compatibility. Cache identities also retain the local unknown order, constraint discovery order and target-element accumulation order. This prevents a parallel insertion race from selecting a mathematically equivalent solve with a different floating-point rounding path. Cached corrector rows use geometric order and are explicitly mapped back to the current local numbering. A cache hit skips local matrix assembly, constraint QR, factorization and solve.

The implementation does not require a global refinement-lineage graph to recognize identical local problems. It still constructs the current patch geometry and inspects its local dependencies. Changed constraints, operator entries, target contributions or unknown sets invalidate reuse. This is exact reuse of an unchanged discrete problem, not approximate reuse across genuinely changed patches.

A run-owned Riesz cache retains immutable local LDLT/Schur factors under exact restricted energy/constraint keys. Production Riesz DOF and constraint ordering is preserved for both manuscript and archived policies; changing that ordering was found to perturb a larger E2 greedy trajectory. A changed Riesz ordering can conservatively miss the cache. Large saddle factors are not retained in this cache. Independent dense-RHS solves share read-only factors and own their workspaces; no global solve lock serializes identical patches. Active contexts keep shared ownership after eviction. Local restrictions and constraint reduction still run before a factor-cache lookup.

`ALOD_PATCH_CACHE_BYTES` and `ALOD_RIESZ_CACHE_BYTES` each default to 64 MiB; zero disables the respective cache. Corrector retention is bounded over two generations. Riesz retention uses a conservative dense-fill allowance. These are retention budgets, not whole-process RSS limits; live factors and temporary assembly buffers also consume memory. Counters distinguish corrector and Riesz hits, misses and retained bytes.

## Training in reduced coordinates

Growing training dictionaries retain a residual Gram matrix and Riesz images. Intermediate greedy evaluations compute member indicators in reduced coordinates and reconstruct only the selected member's seed. They no longer materialize the entire family of fine-grid solutions, residuals and seeds at every block. Near ties, stopping boundaries, cancellation and final publication still use the direct fine-grid evaluation. Memory-budget exhaustion selects the direct path.

POD screening now solves transformed small Schur systems and uses residual and correction-energy Gram matrices. Rejected candidates need neither a fine-grid candidate dictionary nor a full family solution. A potentially accepted candidate is materialized and checked with the original fine-grid PG, error-distance and Riesz gates. Greedy selection is still updated after every enrichment block; blocks are not batched across decisions.

Timing counters include `training_coordinate_evaluations`, `training_fine_evaluations` and `pod_screen_factorizations`. Existing Schur counters describe the growing/direct path; POD screening has its own small-factor count. The training coordinate workspace remains limited by `ALOD_TRAIN_RIESZ_BYTES`, with a more conservative five-buffer allowance. POD workspace remains controlled by `ALOD_POD_GRAM_BYTES` and the configured dense allocation ceiling.

## Shared snapshots and resident auditors

On Linux, the Python runner starts a private local audit service inside the adaptive numerical process. Accepted snapshots can transfer their immutable hierarchy directly to an audit thread. Fine meshes, interpolation, embeddings, energy and Helmholtz operators share ownership, so an audit does not rebuild them on this path. Trial/reduced basis data and accepted values have snapshot-owned lifetimes. Auditing restores the accepted space; it does not retrain its dictionary or recompute correctors.

The scheduler still limits active and drain concurrency. A typical configuration has one 32-thread adaptive team and one 32-thread audit team, then two 32-thread auditors after adaptation completes. The adaptive completion marker is published before audit drain. Solver wall time excludes drain even though the numerical process stays alive to serve pending audit work. Per-thread timing destinations keep audit and adaptive logs separate. The runtime does not impose CPU affinity or a 16-thread ceiling.

`ALOD_SHARED_SNAPSHOT_BYTES` defaults to 256 MiB. At most two unclaimed snapshots are retained, and active shared snapshots also count against the allowance. Old unclaimed snapshots are evicted when necessary. Oversized or evicted states use their durable checkpoints; adaptation never waits for a snapshot slot or an audit. Immutable disk integrity is still checked before consuming an in-memory snapshot. This avoids decoding and operator reconstruction, not every checkpoint read. The allowance estimates retained payloads and hierarchies; audit workspaces and reference factors have separate costs.

The fallback is a persistent native worker per scheduled audit slot. It decodes a checkpoint once for metadata inspection and transfers that owned object to the numerical audit. A 64 MiB `ALOD_AUDIT_HIERARCHY_BYTES` allowance can retain an identical hierarchy between standalone requests. Changed meshes are never treated as matching factors. Queueing, cancellation, completed-audit validation and recovery remain checkpoint-backed.

Operational controls:

- `ALOD_AUDIT_DISK_ONLY=1`: use resident separate processes without sharing adaptive memory.
- `ALOD_AUDIT_ONESHOT=1`: use the previous process-per-snapshot execution path.
- `ALOD_SHARED_SNAPSHOT_BYTES=0`: exercise disk fallback while retaining the integrated audit service.
- `ALOD_REFERENCE_EXECUTION=1`: retain the numerical reference comparison path and separate audit workers.

Ordinary audit exceptions remain isolated from adaptive decisions. A native crash or process-level OOM in an integrated auditor affects the shared process; durable checkpoints support recovery. Use the separate-process mode when that additional isolation is required. Full reference FEM factors remain audit-owned and are not shared across changed meshes.

## Validation and remaining scope

The server regression suite includes actual local refinement and global node renumbering, concurrent shared-factor solves, reduced/direct training at k=8/32/128, audit on/off equivalence, shared-memory/disk/standalone audit equivalence, forced snapshot spill, worker cancellation and checkpoint recovery. A local refinement probe reused 90 correctors and 52 Riesz factors; subsequent fine-node renumbering and a translated-mesh probe each reused 130 correctors. The five-state audit probe consumed five shared snapshots and matched standalone audit observables.

All 50 server CTest cases passed, including an isolated integrity rerun after excluding server-only diagnostic logs; 26 local Python unit tests also passed. The final cached-versus-cold structural probes at k=8, 32 and 128 passed the original strict bitwise-equality assertions. The physical 16-load E2 comparison also passed three consecutive reruns after the cache ordering fix. Every one of 48 workflow/audit thread probes observed 32 threads.

The bounded comparison used three states, 16 training/audit members, k=16, fixed ell=2, 32 adaptive threads and 32 threads per auditor. One warm-up pair was excluded; two measured pairs alternated execution order. E1 used level 10 and radius 0.6; E2 used level 8 and diagnostic radius 10, rather than the production region.

| Median seconds | Previous binary | Current binary | Change |
| --- | ---: | ---: | ---: |
| E1 workflow wall | 22.422 | 22.748 | +1.5% |
| E1 solver wall | 20.693 | 21.256 | +2.7% |
| E1 audit drain | 1.452 | 1.213 | -16.4% |
| E2 workflow wall | 14.459 | 13.676 | -5.4% |
| E2 solver wall | 12.664 | 11.916 | -5.9% |
| E2 training phase | 4.226 | 3.084 | -27.0% |

All compared decisions and audit observables matched within the recorded checks. Each current run consumed three shared snapshots. E2 performed 34 coordinate evaluations and nine fine-grid evaluations, with 18 small POD screening factorizations per measured run. Audit hierarchy-operator preparation summed across the three states fell from approximately 0.219 to 0.002 seconds for E1 and from 0.160 to 0.002 seconds for E2. Audit stage sums overlap adaptation and must not be added to workflow wall time.

Cross-state reuse is now demonstrated, but it is not a universal speedup. On these small patches, key construction and Riesz preparation overhead exceed the saved local solve work; E1 remains dominated by approximately 15.9 seconds of localization checks. Aggregate hit counters include within-state equivalent patches; the separate refinement test establishes actual cross-state reuse. Cache budgets remain configurable, including zero to disable retention. Full-scale tests are needed before claiming a net corrector-construction speedup.

Performance measurements and binary identities are recorded in `provenance/deep_reuse_validation.json`. Bounded timings do not establish full E1/E2/E3 production speedups. Full current-tolerance trajectories, high-k late-state memory pressure and sustained failure/resource tests remain deferred. Remaining optimization candidates include incremental global hierarchy assembly, bypassing unchanged Riesz constraint reduction, cross-refinement quadrature updates, asynchronous durable checkpoint publication and memory-aware audit admission. These are distinct from the exact local reuse and shared snapshot paths implemented here.
