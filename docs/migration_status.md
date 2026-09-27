# Migration status

Current follow-up: [cross-state reuse, coordinate training and shared asynchronous audits](deep_reuse.md). This update supersedes earlier pending statements about local factor reuse, intermediate fine-grid training vectors and resident/shared audit execution. Full production validation remains separate.

Version 4 update (2026-09-22): see [the current implementation review](paper_v4_review.md). Numerical policies and performance statements below describe their original migration stage unless explicitly updated.

Updated: 2026-09-12. **P0-P6 implementation and bounded acceptance are complete.** Full production experiments remain deferred to P8; P7-P9 remain open. Historical full campaigns remain imported evidence; this repository has not generated new 51-state E1 or 33-state E2 production trajectories.

## P0: source discovery and old replay

- Verified the 58 originally cited evidence items and froze the current manuscript inventory (seven figures, eight tables), configurations, source differences, sample roles and canonical audits.
- Recovered the actual archived E1 production source and verified runner, workflow, executable and static-library hashes against production checksums. E2 uses its accepted production snapshot, including the selected-value kernel API absent from the main legacy worktree.
- Rebuilt both old programs from source without reusing old libraries. E1 completed the historically validated 27-state H2/h4 smoke; its 170 sampled family rows and 120 coarse-mark entries agree with the archive. E2 completed two H6/h10 states; all 100 accepted audit rows agree with the archived prefix. Tolerance is 1e-12 + 1e-8*abs(expected); identities and marks are exact.
- Recorded resource-only replay patches, effective configurations, build-hash/environment branch inventories, source hashes, compiler/link commands, actual linked archive members and dependency-version evidence. Historical package locking was not archived; observed host versions are labelled by date.
- Reproduced the manuscript/source averaging discrepancy on a graded mesh. This is an explicit P2 scientific decision. P0 establishes provenance and the counterexample; it does not claim adaptive manuscript equivalence.

See `provenance/p0/README.md` and the accompanying machine-readable evidence. The replay tools require explicitly supplied legacy sources and never run server launchers.

## P1: independent minimal core and baselines

- Independent CMake build with mandatory UMFPACK, Eigen and OpenMP; portable compiler flags and a self-contained CI smoke workflow.
- E1 and alpha=80 E2 manufactured factories only; historical R2 and alternate singular/cutoff profiles removed. Retained numerical formulas agree with the old fixed-grid FEM oracles.
- Mixed-boundary FEM, quadrature, NVB, nested nodal/DG prolongations and stable mesh/element identities. Ancestry records parent IDs and NVB generation changes at refinement-transaction boundaries.
- `alod_run` provides AFEM, UFEM and fixed uniform-grid SLOD. AFEM uses the old conforming strong residual and exact-index Doerfler ordering; UFEM refines all elements. SLOD uses ell=3, H+4 reference levels, direct Schur patch equations and two-sided PG assembly. Exact-error stopping is confined to the baseline controller.
- Strict smoke/full JSON presets, a streaming run wrapper, portable configuration/executable hashes, explicit failure status and state records. Full presets expose historical horizons; their complete campaigns have not been rerun.
- Six old-kernel baseline fixtures, each with three states and complete complex solutions. With matching portable compiler settings, all 18 new states have maximum absolute difference 0, including AFEM indicators/marks and mesh identities. Original SLOD fixtures call the complete old model; the new implementation uses the extracted minimal core.
- Manufactured PDE identities, independent finite-difference derivatives, boundary residuals, quadrature convergence, mesh conformity/area/boundary preservation, nodal injection, lineage, residual reconstruction and Doerfler scaling/tie/error cases are tested. CLI artifact status, invalid inputs, exact stopping and thread reproducibility are covered.

The final clean-checkout CTest and staged portability results are recorded in `provenance/p1_validation.json`. The clean build contains no legacy checkout, old library or preexisting build cache. No new remote CI pass or server performance result is inferred from local tests.

## P2: fixed-state LOD, estimators and localization

- General fixed nested LOD spaces, including locally refined coarse meshes, with direct-Schur correctors, conjugate two-sided PG and a reused coarse factor for multiple RHS.
- Batched constrained energy Riesz contexts, all geometric coarse vertices, per-element mass conservation, shared identical factors with correct multiplicity, and complete-support selected masks.
- Explicit training-member IDs, frozen computed energy normalizers, mean plus worst-member Doerfler marking, nominal `[0]` selection and rejection of audit-only members.
- A separately typed full-reference strong residual evaluated on the embedded LOD solution; AFEM continues to use its own conforming FEM solution.
- Retained small-space and four-vector matrix-free localization iterations, the correct uncorrected coarse-energy denominator, state-bound warm starts, convergence evidence and explicit zero-ratio/nonconvergence handling.
- Explicit manuscript and archive policies for the two discovered numerical differences: area versus arithmetic averaging and N² versus support-expanded vertex patches. Archive equality and manuscript mathematical validation are reported separately.
- Six frozen old-library fixed states with 16 RHS, complete solution/indicator comparisons, exact marking/geometry checks, independent dense-kernel and spectral checks, two-thread reproducibility and CLI failure tests.
- One isolated server window for H6/h10 E2, ell=2, 16 RHS and 16 physical cores, followed by a peak-memory instrumentation confirmation. All workers participated; the matrix-free Ritz iteration converged. This is a bounded fixed-state profile, not a production campaign or speedup claim.

See `p2.md` and `provenance/p2_validation.json` for the acceptance evidence, limits and reproducible commands.

## Reproducibility limits

P3/P4 add the unified reference-sweep controller, E1 fixed/lazy policies and four explicit 48-member controls, E2 regional budget/POD/inheritance, cached AOT factors and every-state/lazy ell policies. E2 lazy reuses E1's timing only and retains the 0.1 threshold. Short tests cover full/partial cycles, same-state promotion and mathematical invariants. The E1 prefix agrees across an ell promotion. E2 target/Theta equality and corrected-training numerical proximity are separate checks. See [P3/P4](p3_p4.md) and `provenance/p34_validation.json` for scope, evidence and pending campaigns.

Native vectorization changes the ordering of nearly equal indicators on an E1 symmetric coarse mesh. `provenance/p1_floating_point_sensitivity.json` records the first observed difference. Matching compiler flags give exact old/new agreement; arbitrary cross-architecture adaptive trajectory identity is not claimed. The original comparison rule remains unchanged, with no new tie tolerance.

The archived quasi-interpolator averages incident elements arithmetically, while the manuscript uses area weighting. Both reproduce coarse P1 functions and coincide on the uniform meshes allowed by P1 SLOD. The graded counterexample has operator difference 0.15023130314433289. P2 also verified that the archived interpolation-support rule can produce N³ vertex patches where the manuscript specifies N². P2 exposes both policies; future trajectory comparisons must select them explicitly. Old archive equality does not imply manuscript-mode trajectory equality.

## P5: checkpoints and independent audits

Accepted and transitional snapshots persist both meshes/lineage, embeddings,
raw kernel/phi, solutions, state counters, lazy checks and Ritz warm blocks.
Atomic commit precedes journal publication; recovery checks the committed
prefix and never repeats an accepted sweep. Three pause positions match an
uninterrupted E2 trajectory. An isolated accepted state-16 snapshot successfully
runs fresh RHS audit without earlier states. Mathematical overrides and altered
frozen members are rejected; new mathematical runs have distinct experiment IDs.

Independent audits provide 48-member e/f/g and E/E_ref/G, AFEM not-applicable
semantics, optional rank-zero, fresh and pure E2 diagnostics. Batch sizes 1/8
agree; audit switches leave production unchanged. Long fields are imported with
an explicit 64 MiB bound and marking sets are external hashed files. Solver,
validation, audit and paper-data completion are separate; injected audit/export
failures preserve the numerical result. See [P5/P6](p5_p6.md).

## P6: measured performance work

Production computes only training plus nominal RHS. Fused quadrature, prepared
regional PG/AOT products, POD linear combinations, accepted training solutions,
ell-only reference/Riesz reuse and direct NVB warm injection reduce repeated
work. The existing factor-cache identity and residual checks remain active.

An isolated PowerShell/SSH/tmux server profile uses H6/h10 E2, ell=2, 16 RHS,
16 physical cores, BLAS one thread and three repetitions per comparison. It
records raw times, medians/ranges, equivalence, live worker activity, RSS/swap
and named dense-buffer allocations. Raw triangular PG and unified symmetric
forward/adjoint factor variants remain measured candidates. No complete
trajectory speedup or 32-thread gain is inferred. Evidence and numerical limits
are in `provenance/p56_validation.json`.

## Next phases

| Phase | Remaining work |
|---|---|
| P7 | Broader integration, resource-failure supervision, portable result packages, optional diagnostic coverage and late-state allocation/timing attribution |
| P8 | E1 main and all four controls; E2 every/lazy and inherit/reset campaigns; complete audits and reproducible paper figures/tables |
| P9 | Final release documentation and full publication acceptance |

Read this file, the implementation map, source manifest and phase plan before continuing. Do not replace accepted-source headers with main-worktree versions, infer new performance from historical timing, or overwrite regression fixtures to hide a failed comparison.
