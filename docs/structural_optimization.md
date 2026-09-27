# Structural optimization review and implementation

Current follow-up: [cross-state reuse, coordinate training and shared asynchronous audits](deep_reuse.md). This update supersedes earlier pending statements about local factor reuse, intermediate fine-grid training vectors and resident/shared audit execution. Full production validation remains separate.

The pre-change baseline is commit `8b4e03f` on `codex/p0-p1`. The isolated server validation directory is `deep-20260926` under the project root. Previously deployed experiments and their binaries are unchanged.

## Implemented changes

1. **One reduced factor per LOD state.** Regional evaluators call `LodSpace::solve_reduced` instead of factoring the same reduced operator again. Access is sequential: a factor is not concurrently solved from multiple workers. The diagnostic reference path retains its independent factor.
2. **Incremental training and POD residual coordinates.** Growing kernel-lift dictionaries append only new lower-left block rows and coarse responses; POD transforms the retained small blocks instead of multiplying through the fine operator again. Within a fixed state, training retains the Riesz image of the base residual and adds images only for new correction directions. Later samples obtain their seeds and estimator actions by linear combination. Cancellation, nearly tied worst members, and near-threshold stopping decisions use direct evaluations. Terminal results are verified with the original local-energy estimator. `ALOD_TRAIN_RIESZ_BYTES` bounds additional estimated dense workspace (default 256 MiB; zero disables it), and batches with at most two members use the direct path. POD can reuse the training images. **POD screening:** For more than two candidate dimensions, construct a residual basis from the base residual and `A*(Phi - B*response)`. Apply the fixed selected-region Riesz operator to that basis in bounded blocks, then form its Gram matrix. Candidate estimator norms use small quadratic forms. Cancellation or proximity to an acceptance threshold falls back to the original evaluation. Every candidate accepted by screening is verified through the original Riesz path. Fine-grid solution reconstruction and PG residual verification remain in the greedy loop. `ALOD_POD_GRAM_BYTES` bounds the additional estimated dense workspace (default 256 MiB; zero disables it); the original path remains available when the budget is exceeded.
3. **Dependency versions before patch assembly.** Exact node/operator/interpolation and element/embedding records are compared once per hierarchy construction. Stable interleaved record slots distinguish nodes from elements as the mesh grows. Patch identities contain dependency version tokens. A hit skips matrix assembly, constraint QR, factorization and solve; only patch geometry is reconstructed to map the cached rows. Cached results are immutable shared objects, so lookup/promotion no longer copies numerical matrices under the cache lock. Reindexing or changes outside a patch can conservatively invalidate it. At this historical stage, reuse was conservative and did not persist Riesz factors across changed meshes; the cache default was zero. Current [deep reuse](deep_reuse.md) provides exact local identities across refinement/reindexing and separate corrector/Riesz retention budgets of 64 MiB each, with zero disabling retention. Its budget bounds retained records/results, not transient dependency-building scratch or process RSS. The measured E1 refinement workload still produced no hits with the conservative dependency identities, so this is not presented as a production speedup.
4. **Shared immutable hierarchy.** Meshes, interpolation, embeddings, fine operators and energy matrices have a common owner. Ell-only LOD reconstruction and the Riesz context share that owner. Destroying the old ell space cannot leave the estimator with dangling data. Riesz preparation uses the existing energy matrix and skips factors for exactly zero-dimensional constrained spaces. The legacy `factorizations()` diagnostic still counts distinct local groups.
5. **Borrowed checkpoint geometry.** Synchronous publication borrows the current immutable geometry instead of copying it into each checkpoint. Geometry identity and mesh fingerprints are computed once per view. A new view is created after each state transition. Existing geometry files are still fully verified, and readers still check complete checksums. Dense state payloads and final durable publication remain synchronous; subsequent [scientific workflow changes](scientific_workflow.md) overlap geometry preparation/publication with useful work.
6. **State-scoped audit geometry and shared analytic evaluation.** Audit batches reuse element gradients without storing all quadrature points. E2 family callbacks have identifiable callable types, allowing the corner singularity and boundary weight to be evaluated once per quadrature point across a compatible batch. Member-specific waves remain separate. Both load/source-moment assembly and exact error integration use this path. Arbitrary replacement callbacks and distinct quadrature contexts retain the general path. Integration rules and within-member accumulation order are unchanged.
7. **Run-scoped asynchronous provenance.** Audit jobs share one executable/source provenance snapshot rather than rescanning the source tree and querying build information for every state. Executable file identity is checked before and after each audit. Checkpoint validation, per-state process isolation, sample coverage, and crash recovery remain independent.

`ALOD_REFERENCE_EXECUTION` retains the earlier numerical fallback paths where needed for comparisons. It does not undo immutable ownership or all housekeeping changes; baseline-binary comparisons are required to measure the entire implementation.

## What is not replaced

The reference FEM solver remains a direct UMFPACK solve with one factor per audited state and batched right-hand sides. The completed production E2 evidence shows 153.67 seconds of reference factorization, compared with 1162.39 seconds of audit load integration and 1259.46 seconds of error integration. A new iterative Helmholtz backend is therefore not introduced without a separate conditioning and high-wavenumber study.

Since this stage, exact local reuse across refinement, reduced-coordinate training, persistent/shared audits and incremental quadrature have been implemented; see [deep reuse](deep_reuse.md) and [scientific workflows](scientific_workflow.md). Further candidates include incremental global hierarchy assembly, eliminating remaining fallback/verification fine-grid materializations, full asynchronous dense checkpoint publication and memory-aware audit admission. Global reference FEM factors are not reused across changed meshes.

## Theta algorithm

Let `B_H` denote the uncorrected coarse basis, `E` the fine energy matrix, `Y_ell` the localized test basis, and `R` the additive local kernel Riesz map. The adaptive check uses

```text
M_H = B_H^* E B_H
D   = A^* Y_ell
G x = D^* R(D x)
Theta_ell = sqrt(lambda_max(G, M_H))
```

For dimensions at most the configured dense threshold (normally 64), the code assembles the Gram matrix in blocks of 16 and solves the dense generalized Hermitian eigenproblem. Larger problems use a matrix-free, energy-orthonormal block Ritz iteration with four vectors, inverse-energy residual preconditioning, and previous search directions. Its expanded search space contains current Ritz vectors, preconditioned residuals and history directions. The coarse energy preconditioner is a sparse LDLT factorization. Warm blocks are transported between states through coarse nodal injection and reorthogonalized.

Current normalized configurations use the relative residual `sqrt(r^* M_H^{-1} r) / max(1e-30, abs(lambda))`, with tolerance `1e-9` and at most 750 iterations. The archived residual normalization uses `max(1, abs(lambda))`. Dense fallback is subject to its configured dimension/resource limit and is disabled by default for larger problems. Failure to converge is reported, not silently accepted. A converged Ritz value is an approximation, not a certified upper bound. Lazy scheduling controls when this computation runs; the absolute threshold and normalized balance check use the resulting value.

## Asynchronous work still worth considering

Production E2 completed 51 accepted states and 51 audits with exit code zero. Its wall time was 3553.91 seconds, adaptive subprocess time 2870.63 seconds, and final audit drain 669.47 seconds. Phase sums are nested and must not be added as disjoint costs. These numbers describe the earlier deployed binary, not this optimization build.

The next scheduling opportunity is a bounded increase in audit concurrency after adaptation finishes, when its 32-thread allocation is free. Persistent workers and immutable in-memory snapshots could additionally remove repeated checkpoint decoding and fine-operator assembly, but would need bounded queues, spill-to-disk, cancellation and lifetime rules. Increasing concurrent audit jobs throughout adaptation can instead create CPU and memory-bandwidth contention. Theta-dependent state transitions must still wait for their localization decision; independent audit work must not drive adaptation.


## Validation scope

The final server Release build passed all 46 CTest checks, including crash recovery, shared checkpoint corruption handling, arbitrary callback fallback, cache-hit resource ceilings, immutable hierarchy lifetime, and direct-versus-optimized training at k=8,32,128. In the nontrivial k=32 fixture, Riesz RHS columns fell from 77 to 33 with working rank 10 and six POD trials. The k=8 and k=128 fixtures fell from 44 to 24 columns. Numerical results, stopping decisions and compressed ranks passed comparison gates.

Bounded workflow benchmarks alternate the baseline and optimized binaries, discard one warm-up pair, and report medians of two measured pairs. Both use three accepted states, 16 members, 32 adaptive threads and one asynchronous 32-thread audit worker. E2 uses level 8/gap 4/radius 10 to exercise training; this is a diagnostic radius, not the production radius 0.6. E1 uses level 10/gap 4. Both versions run under the same updated Python supervisor, so these timings do not isolate the provenance-cache improvement.

For E2, total wall time changed from 27.7031 to 25.6323 seconds (7.5% lower). Audit load integration changed from 2.4607 to 1.4794 seconds and error integration from 2.5576 to 1.6043 seconds. Training changed from 10.8381 to 10.7488 seconds: this small difference is within shared-server variability and is not claimed as a measured training speedup. Fewer Riesz RHS operations do not imply a proportional end-to-end speedup. No full-horizon E1/E2/E3 speedup is claimed.


E1 wall time was effectively unchanged (11.0214 to 10.9972 seconds). With a 256 MiB patch cache enabled in both binaries, the local-corrector phase changed from 0.9669 to 0.4873 seconds, but the new cache still had zero hits in this workload and remained slower than the new no-cache phase (0.4315 seconds). The redesign reduces cache overhead; it does not establish useful cross-state reuse for this refinement sequence. The cache was disabled by default for that measured version; current defaults and later reuse evidence are documented in [deep reuse](deep_reuse.md).

Individual measurements, configurations and binary hashes are retained in [validation evidence](provenance/structural_validation.json). The new binary was tested in an isolated server directory; the completed production E2 deployment was not replaced or rerun.

See [asynchronous audit and training updates](async_training_optimization.md) for the tighter tolerance, bordered solves, audit drain scheduling, and global-solve inventory.

See [asynchronous audits and training updates](async_training_optimization.md) for current scheduling, bordered solves, and the global-solve inventory.
