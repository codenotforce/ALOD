# Localization eigenproblem operator reuse

The localization calculation solves the same generalized eigenproblem
`D* R_AS D x = lambda M_H x` and reports `Theta = sqrt(lambda_max)`.
The tolerance, lazy schedule and promotion criteria are unchanged.

The default implementation now carries the operator images of the Ritz block
and previous search directions through orthogonalization and the Rayleigh-Ritz
update. Only new preconditioned residual directions require an operator action.
Before accepting convergence, it recomputes the action on the Ritz block and
checks the original generalized residual. If cached arithmetic fails this check,
it refreshes the block and discards the old search direction before continuing.

`AdditiveKernelRieszContext::apply_action` uses the same local factors, local
solves, patch multiplicities and deterministic scatter as `apply`. It omits the
second global result array, regional indicators, energy products and diagnostic
residuals. The full estimator API is retained for marking, training and diagnosis.
This change does not introduce a reference-grid Helmholtz solve.

`LocalizationEigenConfig::reuse_operator_actions` and `lightweight_riesz` default
to true. Tests can disable them independently to isolate either optimization.
Adaptive localization events report operator application count, solved columns
and operator-action wall time. Those counters exclude coarse energy/defect
assembly, orthogonalization and the small eigenproblems; the outer `theta` timer
still includes the whole localization calculation.

## Validation

The mathematical regression compares the lightweight action with the full
estimator and an independently assembled dense kernel inverse, including zero
kernel and duplicate patch cases. The localization test checks each switch
separately, compares the combined result with a dense generalized eigensolve,
checks cold and warm starts, and rejects an unconverged iteration.

Build the offline Linux benchmark explicitly:

```sh
cmake --build build --target benchmark_theta
build/benchmark_theta "$CHECKPOINT" 16 32 baseline cold
build/benchmark_theta "$CHECKPOINT" 16 32 optimized cold
```

It requires an accepted checkpoint with a stored LOD basis and uses manuscript
area interpolation and N2 Riesz patches. Supply the checkpoint's wavenumber.
`light` and `reuse` select either optimization separately. `warm` uses the saved
checkpoint block; this is the accepted state's block, not necessarily the block
that preceded its original localization check. A cold comparison uses identical
deterministic initial vectors and is independent of that distinction.

Restore/preparation time and Theta time are separate. Peak RSS includes the
whole benchmark process, including restoration; it is not a Theta-only memory
measurement. Use separate processes and alternate baseline/optimized order.
Concurrent production jobs can affect wall time, so reduced operator-column
counts and numerical agreement should accompany any speedup claim.

The current E2/E3 production campaign retains its frozen executable. Offline
checkpoint validation does not establish a new full production wall time.

## Final E1 checkpoint measurement

The 2026-09-27 benchmark used state 56, ell=4, k=16, 32 threads and identical
cold starts. Two independent runs per mode were executed in AB/BA order.

| Measurement | Reference path | Both optimizations |
|---|---:|---:|
| Mean Theta wall time | 132.885 s | 61.464 s |
| Mean operator-action wall time | 113.240 s | 42.723 s |
| Operator columns per solve | 144 | 56 |
| Ritz iterations | 13 | 13 |
| Relative Ritz residual | 1.590e-10 | 1.590e-10 |

The measured mean speedup was 2.162x. The maximum relative Theta difference was
3.331e-16. Checkpoint restoration/preparation averaged about 47 seconds and is
excluded above. Process peak RSS remained approximately 20 GiB in both modes.
The existing production campaign continued running during measurement; these
are paired observations, not isolated-server or full-trajectory performance.


## Parallel products and local fusion

The next optimization adds parallel global sparse block products and an optional
patch-local fused action. See [operator implementation](theta_fused_operator.md)
for the definition of D, memory/lifetime contracts and benchmark modes. Explicit
`optimized` benchmark mode continues to select the earlier image-reuse and
lightweight-Riesz path, allowing independent comparisons with `parallel` and
`fused`. Production defaults select parallel products based on the final k=128
same-block action measurements, while preserving the generalized eigenproblem
and tolerance.
