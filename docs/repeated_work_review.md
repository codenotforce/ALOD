# Repeated-work review after the E1 timing run

The subsequent [five-optimization implementation](five_optimizations.md)
implements the follow-up opportunities below; this page retains the earlier
review and its benchmark evidence.

## Changes implemented

### Shared quadrature within a batch

Load and audit-error batches previously rebuilt the same physical quadrature
points for every member. Batches whose quadrature class, feature location and
feature scale are identical now generate points once per element. Error batches
also share element geometry and P1 gradients. Source and exact-solution callbacks
remain member-specific. Distinct feature contexts keep separate quadrature.

Elements are partitioned across the requested OpenMP workers. Tiny meshes keep
the member/element split so that sharing does not reduce available parallelism.
Each member retains its original quadrature-point and element accumulation
order. No quadrature order, recursive refinement policy or analytic source is
changed. Scratch holds element contributions, not a global quadrature cache:
48 bytes per member/element for loads and 24 for error energies. For batches
with at least as many members as threads, sharing can use more scratch than the
previous member-only path; the production 16/32 and 8/32 cases already staged
these contributions before this change.

### Residual geometry scoped to one state

Strong residual estimation previously rebuilt element geometry and the edge
incidence map for each of the 16 training members, then recomputed element
diameters in the estimator. A ResidualMeshContext now prepares them once and
serves all members of an accepted state. It is destroyed before mesh transition.
It borrows an immutable mesh; callers must not modify or destroy that mesh while
using the context. Source integration and solution-dependent jumps still run
per member. The standalone scalar estimator keeps its original API.

### Direct patch operator extraction

The direct LOD corrector path previously restricted stiffness, mass and Robin
matrices separately and recombined them into the Helmholtz operator for every
patch. It now restricts the already assembled global Helmholtz operator directly.
The unused component matrices are omitted on this path. Other callers retain
components by default. Constraint rank detection, Schur solves, saddle fallback,
residual gates and basis assembly remain unchanged.

## Disposition of the original optimization opportunities

The following list distinguishes implementation progress after this dated review
from the measurements below. It is not a new performance claim.

| Original opportunity | Current disposition |
|---|---|
| Cross-state corrector and local factor reuse | Implemented with exact identities and bounded retention; see [deep reuse](deep_reuse.md). |
| Load/source-moment reuse in strong residuals | Implemented with checked fallback; see [five optimizations](five_optimizations.md) and [integration](integration_optimization.md). |
| Accepted reduced-operator reuse in audits | Implemented in checkpoint restoration; see [workflow reuse](workflow_reuse.md). |
| Shared immutable checkpoint geometry | Implemented; background geometry preparation is also available. Dense state payload publication remains synchronous. |
| Repeated energy products during E2 orthogonalization/training | Prepared/reduced-coordinate products and POD screening now cover this work; see [space optimization](space_optimization.md) and [deep reuse](deep_reuse.md). |
| Sparse/iterative alternatives to expensive local rank/Schur work | Still require conditioning, acceptance and performance evidence; not an automatic replacement. |

Checksums and restart/audit transaction contracts remain enforced. E1 has no
regional enrichment training, so its timings cannot establish E2 training speedups.
The original experiment logs and measurements below remain unchanged.

## Measured validation

All 41 server Release CTest tests pass, including mixed quadrature contexts,
32-worker participation, residual-context equivalence and direct/component
patch-operator equivalence. Alternating old/new binaries produced exactly equal
accepted-state and audit-sample records after excluding timing fields.

Benchmarks use three accepted E1 states and 16 members, with 32 adaptive threads
and one asynchronous 32-thread audit worker. One warm-up pair is excluded. The
smaller case uses level 8, gap 4 and three measured pairs; the larger uses level
11, gap 4 and two measured pairs. These are not full 51-state reruns. Medians:

| Metric | Smaller before | Smaller after | Larger before | Larger after |
| --- | ---: | ---: | ---: | ---: |
| Workflow wall seconds | 2.4244 | 2.3260 | 38.3515 | 37.7438 |
| Adaptive wall seconds | 1.8275 | 1.8272 | 34.1856 | 33.5844 |
| Adaptive load seconds | 0.2396 | 0.2476 | 1.6745 | 1.6787 |
| Strong residual seconds | 0.3659 | 0.2790 | 2.9647 | 2.2846 |
| Audit error integration seconds | 0.2947 | 0.2546 | 1.7549 | 1.7091 |

Strong-residual improvement persists at about 23%. Larger-case workflow time
falls by about 1.6%; the remaining factorization work dominates. Load integration
has no measured improvement, and the small error-integration gain in the larger
case should not be generalized beyond these runs. No speedup is claimed for E2
training or the full E1 horizon. Machine-readable samples and binary hashes are
in [validation evidence](provenance/repeated_work_validation.json).
