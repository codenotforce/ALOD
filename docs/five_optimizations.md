# Five workflow optimizations

This page records the September 22 implementation and measurements. The later
[structural optimization](structural_optimization.md) moves patch lookup before
assembly, shares hierarchy ownership, and adds POD and analytic-batch reuse.
The initial measurements below remain historical baseline evidence.

Source moments, reduced operators, incremental energy products and shared
geometry are enabled by default in newly built adaptive ALOD E1, E2 and
variable-wavenumber E3 workflows. Cross-state patch caching is implemented
but opt-in because the measured low-hit-rate workloads regressed. Existing binaries and running processes
are unchanged. Mathematical thresholds, quadrature orders, greedy rank limits,
compression gates and lazy ell scheduling are unchanged.

## 1. Cross-state local corrector reuse

A solver-owned cache retains local corrector solutions and their diagnostics,
not LU factors. Full byte identities include local coordinates, wavenumber,
sparse operator entries, constraint entries and RHS entries. Hash-table hash
collisions cannot cause reuse because the complete keys are compared. Changed
operators, constraints or RHS invalidate reuse automatically. Assembly and rank
revelation still precede lookup; a hit skips factorization and solution.

The default budget is zero (disabled). Set `ALOD_PATCH_CACHE_BYTES=268435456`
for a 256 MiB accounted storage budget. Oversized entries are skipped.
The previous generation is read-only during the next construction while new
results enter a separate generation, each with half the total budget. This
prevents early misses from evicting old entries before their turn for lookup.
Keys use sparse constraint entries rather than storing their many zeros. Reused
results are promoted into the current generation. The cache is released when
the adaptive process ends and is not serialized. Cache misses after restart do
not change the numerical solution. Changing every patch can legitimately yield
zero hits; reuse is not an unconditional speedup.

`ALOD_PATCH_CACHE_BYTES` supplies the budget; zero disables it. Timings record
cumulative `patch_cache_hits`, `patch_cache_misses` and accounted
`patch_cache_bytes` per accepted state. Accounted storage is not process RSS.

## 2. Load and strong-residual source moments

Load assembly can also retain element integrals of |f| squared, f times each
P1 basis function, and the quadrature mass matrix. Identical quadrature contexts
share mass matrices. The subsequent strong residual evaluates the expanded
squared norm of f + k^2 n u and its nodal moments without evaluating f again.
The original quadrature rule and points define these moments.

The expanded norm is combined in extended precision. If it is nonfinite or
smaller than 1e-6 times the sum of absolute component terms, that element falls
back to direct residual quadrature. Mesh and quadrature identities are checked.
The frozen source/member association belongs to the caller. Numerical results
are compared with tolerances, not promised to be bitwise equal after changing
floating-point evaluation order.

The default moment budget is 768 MiB. `ALOD_SOURCE_MOMENT_BYTES=0` disables it;
a batch whose estimated buffers exceed the budget uses the original path.
Moments are released immediately after strong residual estimation. AFEM does
not allocate these optional moments. The `source_moment_members` counter shows
whether an adaptive ALOD batch used them.

## 3. Accepted reduced operators in audits

Accepted checkpoints carry both the LOD trial basis and its reduced operator.
Audits restore both and factor the small reduced matrix. They validate dimensions,
finite entries, space identity and three deterministic operator actions against
the restored basis and fine operator. These probes supplement the complete
checkpoint checksum; they are not a full independent matrix reconstruction.
`lod_reduced_restore` replaces `lod_reduced_assembly` in the normal audit trace.
Forced rebuilds and legacy checkpoints still assemble the matrix. The regional
E2 evaluator also uses the LodSpace reduced matrix instead of independently
assembling the same base block.

## 4. Incremental enrichment energy products

E2/E3 orthogonalization caches E*phi columns, computes E*v once at the start and
updates it as each dictionary column is subtracted. It refreshes E*v after each
of the two passes. A final orthogonality gate triggers the original direct
product loop when drift is too large. Existing norm, kernel-lifting, training
stop and compression gates remain in force. The cache is scoped to the training
model and released with it. E1 has no enrichment training.

## 5. Shared immutable checkpoint geometry

Format v3 stores mesh states, lineage and P_node/P_elem/P_dg embeddings in an
immutable content-addressed object under `checkpoints/meshes`. Training,
ell-decision and accepted snapshots on unchanged geometry reference the same
object. The object and each referring checkpoint have complete checksums.
Objects are verified and synced before the referring checkpoint is atomically
published; `latest` advances last. Inspection validates shared dependencies as
well as the checkpoint itself. Missing or corrupt dependencies invalidate a
snapshot so automatic recovery can try an older state. If live reconstructed
geometry regenerates a damaged object, the damaged file is quarantined before
publishing a verified replacement. No automatic geometry deletion is performed.

Copy the entire checkpoint directory to retain resumability. To transfer a
single independent file, first export it with:

```sh
build/alod_run pack-checkpoint results/run/checkpoints/state-000050-phase-0-ell-3.bin results/portable
```

The command prints the exported filename. That file embeds its geometry and
can be copied alone for auditing or recovery. Readers retain v1/v2 support;
older binaries cannot read v3. Checkpoint size accounting includes shared
objects. Sharing reduces repeated writes and storage but still verifies and
hashes geometry; I/O savings are not a guarantee of lower checkpoint wall time.

## Reference execution and validation

`ALOD_REFERENCE_EXECUTION=1` disables these five optimizations for controlled
comparisons; it is a diagnostic environment setting, not a mathematical preset.
Thread counts and CPU affinity behavior are unchanged.

The expanded server test suite covers cache hits and invalidation, cache bounds,
source-moment equivalence and cancellation fallback, reduced-operator rejection,
shared geometry deduplication, corruption and recovery, and portable single-file
audits. E1 and E2 at k=8,32,128 compare training decisions and marking sets exactly,
and numerical observables with relative tolerance 1e-8 and absolute tolerance
1e-11. This is bounded validation, not a new complete paper campaign.

## Initial measurements with patch caching enabled

Two measured alternating pairs followed a discarded warm-up pair. Both cases
use three accepted states, 16 members, 32 adaptive threads and one 32-thread
audit worker. E1 uses level 11 and gap 4. E2 uses level 5, gap 3 and diagnostic
radius 10 to exercise nonzero enrichment; it is not the production radius.
The following measurements used a 256 MiB patch cache, before the decision to
make that cache opt-in:

| Metric | E1 before | E1 optimized | E2 before | E2 optimized |
| --- | ---: | ---: | ---: | ---: |
| Workflow wall seconds | 37.6484 | 36.9899 | 1.8759 | 1.8344 |
| Strong residual seconds | 2.2829 | 0.5105 | 0.1577 | 0.0106 |
| Local corrector phase seconds | 3.8831 | 4.6626 | 0.1438 | 0.1536 |
| Training/LOD solve seconds | 0.8746 | 0.8821 | 0.5412 | 0.5205 |
| Audit reduced assembly/restore seconds | 1.1464 | 0.1095 | 0.0183 | 0.0028 |
| Checkpoint bytes | 335374628 | 293457198 | 9083899 | 7828400 |

E1 had about 400 hits in 16436 lookups; E2 had 108 in 976. Identity construction
and copying outweighed the saved solves on these runs, so patch caching is not
part of the default fast path. Shared geometry reduced storage by 12.5%/13.8%,
but checkpoint wall time increased due to hashing/verification and reduced-matrix
serialization. The new format is retained for its storage benefit and recovery
contract, not presented as a measured checkpoint-time speedup.

See [machine-readable validation](provenance/five_optimizations_validation.json)
for individual measurements and binary hashes. No full-horizon speedup is
claimed from these bounded runs.

## Final default validation (patch cache disabled)

All 43 CTest tests pass for the final build. Two measured alternating pairs,
plus a discarded warm-up pair, use the default no-patch-cache path. E1 uses
level 8/gap 4; E2 uses level 5/gap 3/radius 10. Both use three states and 16
members with 32+32 threads. These differ from the larger cached E1 case above.

| Metric | E1 before | E1 final default | E2 before | E2 final default |
| --- | ---: | ---: | ---: | ---: |
| Workflow wall seconds | 2.3192 | 2.3233 | 1.8814 | 1.7410 |
| Adaptive wall seconds | 1.7814 | 1.6309 | 1.5272 | 1.3746 |
| Strong residual seconds | 0.2696 | 0.0533 | 0.1613 | 0.0108 |
| Training/LOD solve seconds | 0.0462 | 0.0476 | 0.5464 | 0.5074 |
| Checkpoint bytes | 28817823 | 23354998 | 9083860 | 7828497 |

E1 adaptive time falls by 8.4%, while total wall time is effectively unchanged
at this scale. E2 adaptive time falls by 10.0% and total wall time by 7.5% in the
diagnostic training case. The numerical comparison gates pass; marking and
training decisions agree. These measurements do not establish full-horizon or
large E3 speedups. Original full E1 artifacts remain unchanged.
