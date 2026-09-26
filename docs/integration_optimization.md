# Integration optimization

The load, strong-residual and error integration paths share reusable quadrature
storage. Each OpenMP worker owns its scratch vector; no mesh or source samples
remain cached after the call. Recursive Duffy integration grows storage
geometrically instead of reserving exactly enough space at every leaf. Peak
scratch capacity follows the largest element assigned to that worker, rather
than the total mesh size or the number of right-hand sides.

E1 and E2 provide an optional joint value/gradient callback, `exact_jet`.
Error integration evaluates the analytic envelope, singular corner terms and
complex phase once per point. The existing fused audit integration still
computes solution error, reference error and exact norm together. Ordinary
batched error integration now uses that fused path too. Problems without the
joint callback retain the separate-callback fallback. Callers that replace
`exact` or `exact_gradient` must also replace or clear `exact_jet`.

The E2 Gaussian normalization depends only on member parameters and is computed
once when the problem is created, rather than at every quadrature point. This
also benefits load and strong-residual source evaluation.

Quadrature orders, recursive subdivision decisions, point order, accumulation
order and mathematical definitions are unchanged. There is no fast-math mode,
reduced quadrature order, thread-count cap or CPU binding. The strong residual
continues to evaluate the pointwise squared residual directly, avoiding
cancellation from expanding it into differences of large integral moments.

## Validation and reproducibility

`v4_integration` compares joint and separate analytic evaluations, fused and
independent scalar error passes, one/two-thread loads and the fallback without
`exact_jet`, for E1/E2 at wavenumbers 8 and 128. Existing mathematical tests and
historical fixtures remain unchanged. All 39 CTest cases passed on the server,
including asynchronous restart, checkpoint recovery, lazy localization and
historical numerical regressions. Repository integrity checks also passed.

For a before/after benchmark, retain an old executable under a separate
directory with the basename `alod_run`, then run:

```sh
python3 tools/benchmark_integrations.py \
  --before results/integration-before/alod_run --after build/alod_run \
  --output results/integration-benchmark
```

The benchmark runs one warmup per executable followed by three alternating
measurements. It uses E2 H6/h10, 16 training/audit members, three accepted states,
eight explicitly selected threads, one audit worker and eight audit threads.
It compares all accepted-state fields except timings and checkpoint accounting,
and all audit sample fields exactly, including mesh/solution/dictionary hashes,
marking sets, estimators and errors. `--level 3 --gap 2` exercises more recursive
quadrature on a coarse mesh; `--problem E1` selects the smooth case.

No BLAS environment or CPU affinity is imposed by the benchmark. Phase times
are sums across accepted states or audit jobs; they overlap in an asynchronous
workflow and must not be added to infer total elapsed time. These short runs do
not establish scaling on the largest planned meshes or complete paper campaigns.

## Server measurements

Median seconds from three alternating measured runs on 2026-09-22 (one warmup
per executable excluded):

| E2 H6/h10, three states | Before | After | Time reduction |
| --- | ---: | ---: | ---: |
| Load integration | 2.223 | 2.217 | 0.3% |
| Strong residual | 2.436 | 2.396 | 1.7% |
| Audit error integration | 4.395 | 2.459 | 44.0% |
| Numerical solver elapsed | 6.243 | 6.244 | No improvement |
| Complete asynchronous workflow | 10.095 | 8.263 | 18.1% |

On E2 H3/h5, where recursive quadrature requires larger per-element buffers,
load integration decreased from 0.1768 to 0.1570 seconds (11.2%), strong residual
from 0.2169 to 0.2116 seconds (2.5%), and audit error integration from 0.2950 to
0.1649 seconds (44.1%). Complete elapsed time decreased by 7.7%.

The main gain is error integration. Load and strong-residual gains on the finer
mesh are small and should not be treated as a robust general speedup. Expensive
source evaluations and residual topology/assembly remain; the numerical solver
itself has not measurably accelerated. A separate user experiment was present
on the server, so these are controlled alternating runs on a shared machine,
not isolated hardware measurements. Both mesh cases produced exactly identical
accepted-state semantics and audit sample records across all before/after runs.

Raw measurements and executable hashes are in
[`provenance/integration_validation.json`](provenance/integration_validation.json).
