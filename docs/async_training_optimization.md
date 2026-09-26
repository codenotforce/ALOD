# Asynchronous audits, tighter localization, and bordered training solves

## Audit scheduling

The active solver never waits for an audit. Checkpoints remain immutable and durable. The scheduler submits only as many jobs as can run; remaining states stay on disk rather than entering an unbounded executor backlog. Completed sample records are retained for the final export.

`--audit-workers` controls concurrency while the adaptive process is running. After successful solver exit, `--audit-drain-workers` controls concurrency while the remaining audits finish. By default, explicit solver thread counts are reused, adding at most the active number of audit workers without exceeding the released thread budget. Inherited/unknown thread counts do not trigger an automatic increase. A 32-thread adaptive solver and one 32-thread auditor therefore become two 32-thread auditors after solver exit. Set `--audit-drain-workers=1` to retain one auditor and its memory envelope. Each additional auditor owns its own mesh, factors and integration workspace; CPU budget reuse does not make that memory free.

Cancellation stops active subprocesses and leaves pending snapshots available for recovery. A completed audit is reused only after its checkpoint and sample hashes and member coverage match. Run-scoped executable provenance remains shared. The native numerical worker is still a process per snapshot: persistent native workers, shared in-memory snapshots, and factor reuse across different audit states are not implemented. Simply keeping a process alive would not permit reuse across changed meshes or operators. The existing saved LOD basis and reduced matrix are restored instead of rebuilding correctors.

## Coupled training family

Each block adds at most two modes, then solves the whole training family and updates its regional indicators before selecting the next worst member. This greedy order is unchanged; blocks are not batched across selections.

The state-owned coarse reduced matrix is factorized once. All loads share that factorization and are solved as matrix right-hand sides. Eliminating coarse coordinates produces the enrichment Schur matrix. Its growing test blocks now reuse old entries and calculate only new rows/columns. A new border reuses previous factors, solves for the response to the new columns, and factorizes only the new Schur complement. No explicit inverse is formed.

A changed/POD basis rebuilds the small pivoted LU. An ill-conditioned new border or a failed small-system residual check falls back to full LU; the coupled Petrov-Galerkin residual gate remains active. `ALOD_REFERENCE_EXECUTION` retains the independent full-factor path for numerical comparison. Timing logs record `training_schur_border`, `training_schur_factor`, and actual `training_schur_updates`/`training_schur_rebuilds` counters. The rank cap is 24, so these dense factorizations were already small: reduced operation counts should not be interpreted as a large end-to-end training speedup.

An iterative warm start is not added on top of these direct solves: direct back substitution does not use an initial guess. Switching the coupled solve to Krylov iteration would need a separate convergence and cost comparison.

## Localization algorithm and tolerance

All current adaptive/fixed presets and C++/Python defaults now use `ritz_tolerance = 1e-9`. Historical fixture configurations remain frozen. Resuming an old experiment preserves its original mathematical configuration; changing tolerance requires a new run.

For coarse dimension above the configured dense threshold (64 by default), Theta uses a custom matrix-free, preconditioned four-vector block Ritz method of the LOBPCG type. Its search space contains current Ritz vectors, preconditioned residuals and previous search directions, with double energy orthogonalization. It is not a call to an external LOBPCG library. Smaller spaces use the dense path; normalized configurations use the largest dense Hermitian eigenpair.

The generalized problem is `G x = lambda M_H x`, where `G = D* R D`, `D = A* B_test`, and `M_H = B_H* E B_H` uses the uncorrected coarse basis. Theta is `sqrt(lambda_max)`. Applying R uses local constrained energy solves. The preconditioner factors the coarse energy matrix, not a global fine-grid matrix. Normalized configurations stop on `sqrt(r* M_H^-1 r) / max(1e-30, abs(lambda)) <= 1e-9`. Archived residual normalization remains distinct.

At tight tolerances, the constraint-range component is removed from D using `D <- (I - I_H* P_H*) D`. Since `R I_H* = 0`, this preserves `D* R D` in exact arithmetic and reduces cancellation when the defect is small. Small residual directions are normalized before orthogonalization rather than discarded by an absolute cutoff. Dense residuals are also checked against the requested tolerance. No tolerance relaxation or unrestricted dense fallback is used.

The inspected legacy source header has a localization default of `1e-7`; the supplied E2 figure dataset contains Ritz residuals around `1e-5` to `1e-4`. These inspected files do not establish a legacy production setting of `1e-9`. The new value is applied explicitly here regardless of those older settings.

## Global fine-grid solve inventory

| Path | Global fine-grid linear solve? | Purpose |
| --- | --- | --- |
| Current kernel-lift ALOD adaptivity, E1/E2/E3 | No | Correctors and additive Riesz solves are local; coupled and energy solves are coarse/reduced. |
| Theta localization check | No | Local Riesz action, coarse energy preconditioning, small Ritz eigenproblems. |
| Strong residual, load integration, error integration itself | No | Element/boundary quadrature and sparse products. |
| Standard audit reference solution | Yes | One UMFPACK factorization per snapshot, reused across batched RHS, for reference error floor and reference-to-ALOD gap. Exact-solution error alone does not require this solve. |
| Optional audit reference iterative refinement | Yes | Additional back substitutions with the same reference factor. |
| Archived adjoint enrichment tests | Yes | Solve `A* y = E phi`; cached full fine-grid adjoint factor, shared across RHS. This is outside reference-error computation and is absent from kernel lifting. |
| Optional fresh training audit with kernel lifting | No additional global fine solve | Fresh local training/coarse solves; the standard reference solve still supplies audit comparison metrics. |
| Adaptive AFEM or fixed AFEM/UFEM baseline | Yes, on the baseline's own conforming mesh | This is the baseline numerical solution, not an auxiliary reference solve for ALOD. |

Relevant implementation locations are `apps/adaptive_run.cpp`, `apps/audit_run.cpp`, `src/enrichment/aot.cpp`, `src/enrichment/regional.cpp`, `src/localization/theta.cpp`, and `src/fem/batch.cpp`.

## Validation

The isolated server Release build passed all 47 CTest tests in 42.24 seconds. Tests include a complex nonsymmetric 24-by-24 growing matrix with 16 RHS, unsafe-border fallback, full-factor versus optimized training at k=8/32/128, cold/warm matrix-free eigenpairs versus dense results at 1e-9, equivalence of the dual projection, repeated ell promotion, restart, cancellation and active/drain concurrency bounds. A one-iteration nonconvergence test uses a genuinely nontrivial spectrum; a tiny exactly zero defect can correctly converge immediately.

The bounded E2 benchmark uses five accepted states, 48 members per state, production radius 0.6, 32 adaptive threads, and one active 32-thread auditor. One warm-up pair is discarded, followed by two alternating measured pairs. The only comparison variable is drain concurrency (one versus two auditors), using the same optimized binary and 1e-9 tolerance. Median wall time was 6.6638 versus 5.7015 seconds (14.4% lower); drain wait was 3.3524 versus 2.3494 seconds (29.9% lower). Solver time was 3.1454 versus 3.1972 seconds. The measured peak audit concurrency was one before solver exit and at most two afterwards. Every run recorded eight bordered Schur updates and fifteen full small factorizations, including POD basis changes. Those counts do not isolate the speed benefit of bordered solves.

All marking and training decisions matched; floating observables passed rtol=1e-8 and atol=1e-10. Maximum recorded Ritz residual was below 1e-9. Individual measurements and the executable hash are retained in [validation evidence](provenance/async_training_validation.json). Reproduce with `python tools/benchmark_async.py --executable build/alod_run --output results/async-drain --repeats 2`. These short-run measurements are not a full-horizon speedup claim; full E1/E2/E3 runs at the new tolerance remain future work.
