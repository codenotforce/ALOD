# P1 baseline contract

Version 4 update (2026-09-22): see [the current implementation review](paper_v4_review.md). Numerical policies and performance statements below describe their original migration stage unless explicitly updated.

`alod_run` implements nominal AFEM, UFEM and SLOD for E1/E2. `tools/run_baseline.py` validates JSON presets before execution and streams run artifacts. Neither entry point interprets historical build-hash tokens or family-control environment variables.

| Field | Meaning |
|---|---|
| schema_version | Exactly 1 |
| problem / method | E1 or E2; AFEM, UFEM or SLOD |
| initial_level | Number of uniform NVB refinements from the problem's initial mesh |
| states | Maximum number of solved states, including state 0 |
| theta | Fraction of the total squared strong residual selected by AFEM; defaults to 0.15 |
| target | Exact relative energy error stop; 0 disables it. Evaluated after solving and independently of marking. |
| maximum_nodes | Mesh allocation guard; SLOD also checks a conservative reference-mesh bound |
| threads | OpenMP team size; BLAS remains at one thread |
| emit_solution | Include all nodal complex values in state records |

Unknown fields, missing fields, duplicate CLI options, malformed/nonfinite values and unsupported problems/methods fail explicitly. JSON numbers used as integer counts must be actual integers, not booleans. Output directories must not already exist. A nonzero child exit records a failed run and retains completed state records.

AFEM solves with UMFPACK, computes the conforming P1 strong residual, selects a Doerfler set, then calls the unchanged NVB implementation. Body, interior jump, Neumann and Robin terms remain separate in the estimator. Interior edges contribute half to each neighbor. Element-index tie breaking follows the archived exact floating-point comparison. No exact error, test or shift sample enters marking.

UFEM marks all elements. SLOD also uniformly refines H, rebuilds h at H+4, and uses ell=3. It retains the original patch geometry, constrained direct-Schur equations, correction packing cutoff, corrected trial/test basis assembly and complex-symmetric Petrov-Galerkin convention. The P1 patch loop is sequential; FEM and integration retain their existing OpenMP code. Shared patch scheduling and cross-state caches are later optimizations. Sparse factorization failures and residual violations are explicit failures, without hidden fallback.

Stable element IDs are monotone and never reused. An unchanged element keeps its ID. Each new element records its parent from the preceding refinement transaction and the total NVB generation increment; closure may contain several bisections within that transaction. The in-memory ancestry registry retains retired parents. Node indices remain stable under the archived append/merge ordering. FNV-1a mesh fingerprints encode IEEE-754 coordinates, connectivity and boundary tags in fixed byte order; they are diagnostic identities, not cryptographic hashes. The standalone baseline driver does not serialize checkpoints. P5 restart support is available in the adaptive ALOD/AFEM controller; see [P5/P6](p5_p6.md).

Historical effective configurations and full horizons are in `provenance/baseline_contracts.json`. The migration plan explicitly selects UMFPACK and direct Schur; historical baseline runs sometimes used Eigen SparseLU or direct saddle. Old/new extraction tests use the same explicit backend, while historical timing remains separate. The full E1 targets/horizons are AFEM 0.008 with at most 180 updates, UFEM H6-H20, SLOD H6-H16; E2 uses AFEM 0.003 with at most 160 updates, UFEM H6-H19, SLOD H6-H15. Full campaigns were not rerun in P1.

## Comparison scope and floating-point sensitivity

Six old-kernel fixtures cover three solved states per baseline and problem at H2, including complete solutions, norms, residuals, indicators and exact marked identities. The old SLOD probe calls the complete original model, while ALOD uses its extracted independent core. The probes share deterministic scheduling and output serialization; mathematical and mesh tests independently validate these shared parts.

Old and new libraries must use the same architecture flags for a meaningful extraction comparison. Native CPU vectorization changes rounding: in a symmetric E1 AFEM state, two maximal indicators differ by about 1e-16 and can exchange order. This is recorded in the validation evidence. The migration retains the original strict comparison rather than introducing an undocumented tie tolerance. Identical adaptation across arbitrary architectures is not claimed. Portable old/new comparisons use the same compiler, default architecture and BLAS settings.

The manuscript specifies area-weighted quasi-interpolation, whereas the archived implementation uses arithmetic averaging. P1 restricts SLOD to uniform meshes where the definitions agree. The graded-mesh counterexample remains a required test and a P2 scientific decision.
