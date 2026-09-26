# LOD, enrichment and training optimization

The current kernel-lifting implementation avoids constructing the test-energy
factor used only by the archived adjoint-test policy. The trial-energy
projection factor is built lazily when a training vector is actually added.
Frozen audit evaluations and empty training regions do not need that factor.
Both policies retain their original coupled Petrov-Galerkin equations and
residual gates.

Regional evaluators reference immutable trial, test and Helmholtz sparse
matrices owned by `LodSpace`, instead of copying all three. As before, the
`LodSpace` must outlive each evaluator. The product of enriched tests with the
Helmholtz operator is reused for the lower-left block and diagnostics, without
changing multiplication order. Dictionary diagnostics remain cached until the
dictionary changes; changing the load block still refreshes the coarse solve.

LOD patch assembly classifies natural-boundary contact once per fine element.
Target patches reuse those flags, and an explicit all-Dirichlet boundary needs
no edge-count table. Each fine element's three-by-three coarse embedding is
read once for its corrector right-hand sides, replacing repeated sparse-entry
lookups. Element correctors are released before extracting the trial basis.
The full corrected basis and extraction triplets are released before allocating
the test basis and factorizing the reduced operator.

Selected-region Riesz applications used in training assemble their common
result once in the original deterministic order. A contiguous copy supplies
the second public result field, replacing a duplicate scattered accumulation.
Full-estimator assembly is unchanged.

## Verification

The opt-in `benchmark_spaces` target measures LOD construction, Riesz
preparation, training, and frozen enrichment setup/cold/hot solves separately.
The same probe source can be linked against retained before/after core
libraries. `tools/benchmark_spaces.py` alternates executables, discards one
warmup each, records peak process RSS with the Linux time utility, and requires
identical solution/dictionary/kernel hashes, estimators, rank, stop reason and
training/POD counts.

```sh
cmake --build build --target benchmark_spaces
python3 tools/benchmark_spaces.py \
  --before results/spaces-before/benchmark_spaces --after build/benchmark_spaces \
  --members results/completed-run/members.txt --output results/spaces-probe \
  --problem E2 --level 10 --gap 4 --ell 2 --wavenumber 64 --threads 8
```

The probe trains on every member in the supplied table. Use a training-only
table to reproduce the measurements. It uses regional radius 0.9, rank cap 24
and no inherited dictionary to exercise growing and compressed spaces.
Explicit probe thread counts do not impose a production cap or CPU binding.

The patch assembly regression independently reconstructs right-hand sides and
boundary contact on graded E1/E2 meshes, including implicit boundary tags.
Kernel-lifting tests check cached block diagnostics, new right-hand sides with
the same dictionary, changed dictionaries and independent coupled solves.
Existing tests cover the archived adjoint policy, empty regions, inheritance,
mesh repair, checkpoint recovery and asynchronous audits. All 40 CTest cases
passed on the server; source provenance and absolute-path checks also passed.

These changes reduce repeated work and overlapping allocations. Local patch
factorizations, coarse sparse factorizations and regional Riesz solves remain
part of the algorithm; bounded probes do not establish million-node scaling
or completion of the deferred paper campaigns.

## Server measurements (2026-09-22)

Three alternating measured runs per executable, excluding one warmup each,
used 16 E2 training members and eight explicitly requested threads. Both grids
used gap 4 and ell 2. The H8 grid used k=32 and 12,545 reference nodes; H10 used
k=64 and 49,665 reference nodes. Each run reached working rank 16, compressed
to rank 6, with nine training evaluations and six compression trials.

| Phase / process metric | H8 before | H8 after | H10 before | H10 after |
| --- | ---: | ---: | ---: | ---: |
| LOD construction (s) | 0.649 | 0.629 | 2.946 | 2.894 |
| Riesz preparation (s) | 0.132 | 0.132 | 0.559 | 0.560 |
| Training including POD (s) | 2.110 | 1.962 | 10.213 | 9.168 |
| Frozen evaluator initialization (s) | 0.156 | 0.054 | 1.500 | 0.508 |
| Frozen cold evaluation (s) | 0.077 | 0.068 | 0.365 | 0.344 |
| Frozen hot evaluation (s) | 0.043 | 0.037 | 0.182 | 0.159 |
| Complete probe elapsed (s) | 6.233 | 5.987 | 28.088 | 26.094 |
| Peak process RSS (KiB) | 343120 | 316304 | 1377964 | 1254144 |

For H10, training time decreased by 10.2%, evaluator initialization by 66.2%,
and peak RSS by 9.0%. Complete probe time decreased by 7.1%. LOD construction
improved only 1.8%; it should not be described as a major construction speedup.
The initialized evaluator's coarse factor is separate from the enriched Schur
block prepared during its cold evaluation. RSS covers the whole single probe
process, not only the LOD phase or concurrent production audit processes.

No thread affinity or BLAS thread environment was imposed. These measurements
come from a shared server; small differences near 1-3% are not reliable general
scaling claims. The complete probe also includes unchanged load integration.
All reported numerical outputs matched exactly across before/after runs. A
separate three-state asynchronous workflow comparison checked accepted-state
semantics and every audit sample, including refinement marks and hashes.

Executable hashes and unrounded measurements are recorded in
[`provenance/space_validation.json`](provenance/space_validation.json).
