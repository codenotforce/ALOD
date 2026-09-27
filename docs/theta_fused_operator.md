# Parallel and patch-local localization operators

## Meaning of the defect matrix

Let `A_h` be the fine-grid Helmholtz matrix and let the columns of `B_test`
represent the LOD test basis in fine-grid coordinates. The adjoint defect is
`D = A_h.adjoint() * B_test`. Thus `D.adjoint() * v` gives the tests of the
fine-grid vector `v` against the discrete Helmholtz operator. Its dimensions
are fine-grid nodes by free coarse coefficients. It is independent of the
current load, coupled solution and enrichment training dictionary.

For accurate small defects, the code uses
`D_tilde = D - I_H.adjoint() * (P_H.adjoint() * D)`.
The additive kernel inverse annihilates the range of `I_H.adjoint()`, so this
subtraction leaves `D.adjoint() * R_AS * D` unchanged in exact arithmetic.
It reduces cancellation in the constrained kernel calculations.

The generalized eigenproblem is `G x = lambda M_H x`, where
`G = D_tilde.adjoint() * R_AS * D_tilde` and
`M_H = B_coarse.transpose() * E_h * B_coarse`.
Theta is the square root of the largest eigenvalue. The relative Ritz residual
requirement remains 1e-9 in current production configurations.

## Two equivalent execution paths

The parallel global path stores a CSR view for `D X` and uses the existing CSC
matrix for `D.adjoint() Y`. Each worker owns complete output rows, with a fixed
nonzero summation order. Dynamic scheduling balances rows of different lengths;
it introduces no concurrent writes to the same output coefficient. This path
retains the ordinary global Riesz application. Its overlapping patch sums use
an inverse incidence map to gather each fine-grid output row in parallel.
Entries retain the original group order, including duplicate multiplicities
and selected-region masks. The optional map is capped at 512 MiB and the
existing resource limit; larger maps fall back to the original serial scatter.
On the final k=128 checkpoint this map used 310,316,808 bytes (296 MiB)
and took 0.387 seconds to construct.
`LodLimits::parallel_riesz_gather` disables this map for reference comparisons.

The fused path uses `R_AS = sum_j J_j R_j J_j.adjoint()` to evaluate
`G X = sum_j D_j.adjoint() R_j (D_j X)`, with `D_j = J_j.adjoint() D_tilde`.
It performs each constrained solve and its coarse projection inside a patch.
Duplicate-patch multiplicities are preserved. Patch results are gathered by
coarse output row in fixed group order, in parallel, without atomics.

Both paths share the same immutable local factors and the same Schur/saddle
solve implementation. No global fine-grid Helmholtz solve or additional local
factorization is introduced. The fused path has no global fine-grid dense RHS
or global Riesz result for ordinary iterations. A single shared CSR matrix is
used; patch-local copies of its nonzeros are not retained. Per-patch coarse
support maps retain unique column indices. CSR storage and projected block
allocations are checked against the resource policy. When the extra CSR view
exceeds the policy, Theta falls back to its original sparse products.

After apparent convergence, the eigensolver recomputes an action using the
global Riesz/scatter path and checks the original generalized residual. A failed
check refreshes the block and discards the search direction. Small dense
cross-checks also assemble their Gram matrix with the global path.

## Lifetime, concurrency and measurement

The operator borrows the immutable CSC defect and Riesz context; both outlive
all actions. Temporaries are rejected by the constructor. Fine-grid and patch
workspaces are call-local. Scoped execution budgets also apply when Theta runs
concurrently with training. The local solve implementation lives in an internal
header shared by the estimator and defect operator rather than being duplicated.

`LocalizationEigenConfig::parallel_defect` and `fused_defect` select the two
paths. Parallel global execution is enabled by default; local fusion is optional.
On the final k=128 checkpoint, a fixed four-column action averaged 8.703 seconds
with the earlier global path, 4.901 seconds with fusion and 2.287 seconds with
parallel global products. The latter was therefore selected as the default.
Fusion differed from the reference by 4.60e-15 relatively; parallel products
matched exactly in that probe. These are action timings, not full-solve speedups. Existing `reuse_operator_actions` and `lightweight_riesz` switches remain
available for reference comparisons. Explicit benchmark modes keep historical
comparisons independent of changed defaults:

```sh
cmake --build build --target benchmark_theta
build/benchmark_theta "$CHECKPOINT" 128 32 optimized cold
build/benchmark_theta "$CHECKPOINT" 128 32 parallel cold
build/benchmark_theta "$CHECKPOINT" 128 32 fused cold
build/benchmark_theta "$CHECKPOINT" 128 32 probe cold
```

`optimized` selects the earlier image-reuse/lightweight implementation.
`probe` compares repeated actions on an identical deterministic four-column
block using global, fused and parallel paths. It does not run the eigensolver.
Cold eigensolver comparisons use identical initial vectors; warm mode uses the
saved checkpoint block, which may already approximate the converged state.
Restoration and Theta computation are timed separately.

Timers distinguish operator preparation, fused patch work, fused coarse gather,
parallel global RHS and dual products. Sparse products used to construct the
defect matrix and cancellation projection still have separate outer timers.
Neither matrix-free iteration nor local solves imply that these global data
operations are free. Runtime and peak-memory measurements determine which path
is appropriate for a given large state.

## Full checkpoint validation

On the final E3 k=128 checkpoint (state 29, ell 3), identical deterministic
cold starts gave 70 iterations and Theta 0.0051645682358023703 with relative
residual 8.9932201220623722e-10. The earlier image-reuse/lightweight version
took 660.786 seconds; parallel sparse products took 218.433 seconds (3.03x).
Checkpoint restoration, measured separately, took approximately 151 seconds.
These measurements precede the additional parallel Riesz gather. They are
single runs, with some shared-server build/probe activity during the earlier
run, rather than an isolated-machine repeated benchmark. The production warm
start is not used as the speedup denominator. No complete experiment rerun is
needed to validate this fixed-state operator change.

### Final gather-enabled result

The final implementation completed the same 70 iterations with exactly the
same reported Theta and residual. Theta wall time was 188.741 seconds, a 3.50x
speedup over the 660.786-second comparison. Operator applications took 137.883
seconds cumulatively. The deterministic gather took 2.596 seconds across 71
applications, versus 33.127 seconds for the preceding serial scatter.

The final run's checkpoint restoration took 205.998 seconds and is excluded
from Theta time. Process peak RSS was 62,604,472 KiB, including restoration.
The extra CSR defect view requires approximately 5.01 GiB for this matrix's
268,558,320 nonzeros. Sparse construction of the defect and cancellation
projection still took 16.837 and 25.520 seconds respectively; these are not
global fine-grid Helmholtz solves.

All eight targeted regressions passed: zero kernels, relative Theta, shared
audit transport, fused/global and serial/parallel gather equivalence, execution
context, asynchronous execution equivalence, normalized ell and P2 mathematics.
The fixed-state records are retained under `results/theta-fused-optimization/`
(ignored experiment artifacts). This validation does not predict whole-campaign
speedup or replace a new campaign when new experimental results are required.
