# Fine-grid stability and localization-threshold calibration

The diagnostic measures the manuscript's fine-grid inf-sup constant, not the
online LOD inf-sup constant or a load-dependent error ratio. With homogeneous
Dirichlet degrees of freedom removed, let

\[
A=K-k^2M-ikB_R,\qquad E=K+k^2M.
\]

Then

\[
\gamma_h(k)^2=\lambda_{\min}(A^*E^{-1}A,E).
\]

The executable `alod_stability` computes the largest eigenvalue of
\(T=A^{-1}EA^{-*}E\), which is self-adjoint in the energy inner product.
Its reciprocal square root is the reported inf-sup estimate. A single sparse
UMFPACK factorization supplies forward and adjoint solves because the assembled
operator is complex symmetric. The diagnostic uses the 64-bit sparse-index
interface and METIS ordering above 256 free degrees of freedom to avoid the
large-workspace limits encountered with the 32-bit AMD path. No inverse or normal matrix is formed. A block
Rayleigh-Ritz iteration uses the energy inner product and checks its relative
residual against `1e-9`. One inverse polishing step is followed by an independent
residual check in the original generalized eigenproblem. Dimensions at most 256
are also checked against the smallest singular value of the Cholesky-whitened
operator. A converged Ritz pair is not a certified bound on the global extreme.

Build and validate:

```sh
cmake --build build --target alod_stability
python3 tests/check_stability.py build/alod_stability
```

Run the resumable E1/E3 calibration:

```sh
python3 tools/run_stability_calibration.py \
  --executable build/alod_stability --output results/gamma-calibration
```

The five wavenumbers are 8, 16, 32, 64 and 128. Initial reference levels are
7, 10, 13, 16 and 19, matching the E3 presets. The campaign repeats these with
two independent random starts, then adds two uniform NVB levels to test mesh
sensitivity. Each job records its mesh size, free dimension, both residuals,
iteration count, assembly time, factorization time and eigensolve time.
`summary.json` is updated atomically after each successful job. The runner
resumes completed result files in the same dedicated output directory.
Use a fresh directory after changing the executable or experiment definition.

For each mesh offset separately, the proposed relative calibration is

\[
\tau(k)=0.2\,\frac{\gamma_{h(k)}(k)}{\gamma_{h(16)}(16)}.
\]

This preserves the sampled ratio of defect tolerance to fine-grid stability.
It does not prove that the anchor 0.2 meets every theorem constant. Other
stability, interpolation and decomposition constants can also depend on the
mesh or geometry. The diagnostic never edits production threshold presets.

The optional log-log power fit is descriptive. Inspect per-wavenumber values,
independent-start agreement and mesh sensitivity before using a fitted law;
mixed boundary conditions can produce resonance-dependent, nonmonotone values.
Uniform-mesh checks do not certify the smallest stability constant along every
adaptive trajectory. E2 has a different domain and boundary partition and
requires a separate calibration (`--problem E2`); E1 values must not be silently
reused for E2. The optional E2 sweep uses the same refinement levels for
comparison, not a claim that its physical mesh diameters match E3.
