# Reference FEM solver comparison

The optional `benchmark_reference_solvers` target compares solvers on the fine
mesh in an existing E1 k=16 checkpoint. It reads the checkpoint and assembles
the reference operator and the selected manufactured RHS members directly; it
does not construct a LOD space or modify production outputs.

```sh
cmake --build build --target benchmark_reference_solvers
build/benchmark_reference_solvers smoke 2 1 all
build/benchmark_reference_solvers "$CHECKPOINT" 32 1 ilut
build/benchmark_reference_solvers "$CHECKPOINT" 32 1 energy 1e-14
```

The target is excluded from the default build and is available when CHOLMOD is
installed. This optional diagnostic does not add a production dependency.

## Compared methods

- Direct: complex UMFPACK LU, 64-bit indices and METIS ordering. One factorization
  is shared by all requested RHS members.
- GMRES with ILUT: drop tolerance 1e-3 and fill factor 10.
- GMRES with energy preconditioning: `P = K + k^2 M`, factored with real CHOLMOD
  supernodal Cholesky and METIS. Its inverse is applied to the real and imaginary
  RHS parts separately. This is a one-level energy preconditioner, not AMG or a
  two-level Helmholtz preconditioner.

Both iterative variants use Eigen GMRES, zero initial guesses, restart length
80 and a maximum of 400 iterations. Internal tolerance defaults to 1e-12 and
can be overridden by the final argument. In this implementation the internal
preconditioned residual and the original system residual are different; the
benchmark checks the latter explicitly. See the [Eigen GMRES documentation](https://eigen.tuxfamily.org/dox/unsupported/classEigen_1_1GMRES.html)
for the solver interface and the [ILUT documentation](https://eigen.tuxfamily.org/dox/classEigen_1_1IncompleteLUT.html)
for its two threshold parameters.

Acceptance requires a finite solution, `||Ax-b||_2 / ||b||_2 <= 1e-9`, and an
energy-norm difference from the direct solution at most 1e-7. A preconditioned
stopping flag alone never counts as success. An unconverged run reports its
measured time and residual but cannot be described as a completed solve.

## Timing scope

Common checkpoint reading, matrix reduction and load assembly are reported
separately. The compared total is setup/factorization plus numerical solution.
Residual and energy-difference checks are outside the solution timer. No CPU
affinity is imposed; 32 requested OpenMP threads do not imply that every backend
phase can use 32 cores efficiently.

The first comparison uses nominal member 0. A single-RHS timing must not be
multiplied by 48 to predict the audit family time: both direct and iterative
methods can reuse setup across that family. The benchmark accepts an RHS count
up to 48 for explicit follow-up measurements. Results from these preconditioners
do not establish a universal comparison with multigrid, deflated, or two-level
methods, and must not automatically replace the production solver.

## E1 state 56 measurement, 2026-09-27

On the final E1 reference mesh (667,556 free DOFs), nominal member 0 and 32
requested threads, the paired direct solve took 15.942 seconds (15.400 setup,
0.543 solve). GMRES with energy Cholesky took 69.307 seconds (5.802 setup,
63.505 solve), 233 iterations and achieved a true residual of 5.070e-11; its
relative energy difference from the direct solution was 3.128e-12. The direct
solver was 4.347x faster in this comparison. Its three baseline runs ranged
from 15.942 to 17.255 seconds while production continued concurrently.

ILUT exhausted 400 iterations in 92.317 setup/solve seconds with true residual
3.022e-2, so it was not an acceptable solve. Energy preconditioning at internal
tolerance 1e-12 also missed the true-residual gate (2.462e-9); the accepted run
used internal tolerance 1e-14. These failed precision checks remain in the raw
records. No full 48-RHS speedup or general conclusion about multigrid follows
from this nominal-member test.
