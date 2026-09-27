# Code structure review and tuned figures

This review covers the adaptive/audit entry points, checkpoint ownership,
integration caches, Python experiment orchestration and paper delivery. The
numerical library already separates mesh operations, FEM assembly, LOD spaces,
regional training and localization. Those boundaries remain useful; the main
maintenance problems were growing application functions and mixed reporting
responsibilities.

## Changes

| Area | Finding | Change |
| --- | --- | --- |
| Paper delivery | One function combined data selection, plotting, table arithmetic and packaging. | `tools/paper_delivery.py` is a CLI; `paper/data.py`, `selection.py`, `figures.py`, `tables.py` and `bundle.py` have separate responsibilities. |
| Plot conventions | E3 incorrectly reused the completed-cycle sampling used for E1/E2. | E3 keeps all accepted states through the first one-percent crossing, including intermediate reference sweeps. Unreached endpoints remain visible. |
| Plot styles | Generic styles discarded the supplied, tuned parameters. | Separate Matplotlib contexts preserve the two reference scripts' fonts, colors, sizes, line widths, markers, annotations and axis settings. E3 is named `kappa_error_dof.pdf`. |
| Scientific table arithmetic | A large reporting function mixed unrelated metrics. | Terminal fits, family first crossings, family budget interpolation and kappa endpoints each have a named function. |
| Paper validation | The PG residual gate checked only the nominal member. | Every imported member is checked for a finite, admissible residual and a positive finite error. |
| Final meshes | Filesystem traversal order could select an intermediate audit mesh. | Native delivery uses the terminal mesh export or exports the accepted checkpoint. Arbitrary audit meshes are not used as terminal evidence. |
| Audit diagnostics | Reference deepening and two-level calculations expanded the main solve loop. | Their implementation is in `apps/audit_diagnostics.cpp`, with an explicit, small interface. |
| Region diagnostics | The same fine mesh was copied and filtered once per RHS. | Select the region once per audit and reuse it for all member batches. |
| Reference deepening | Reading the previous embeddings made a deep copy. | Borrow them until construction of the new space completes. |
| Configuration compatibility | New runs and resume validation duplicated migration defaults. | `tools/adaptive_config.py` is the shared definition, distinct from current production presets. |

The rendering module ports parameters from `generate_kernel_lift_figures.py` and
`plot_results.py`. Their SHA-256 identities are included in the delivery report.
The external scripts and source archives are read only; host-specific imports
and absolute paths are not copied into the repository. The portable package
includes every Python module needed for regeneration.

## Comments and invariants

Comments explain contracts that are easy to break during optimization:

- A checkpoint geometry view borrows immutable meshes. The background future
  must be joined before reading its prepared metadata, refining or destroying
  the meshes. Destruction order also matters during exception unwinding.
- Exact-integration caches retain analytic jets, never the current discrete
  solution or its gradients.
- Reference deepening injects the frozen space; retraining is a separate choice.
- The two-level diagnostic is an energy Frobenius norm, not a certified operator
  norm or a single-RHS error.
- Region diagnostics select whole cells by their centroids. They do not integrate
  clipped curved cells.
- Compatibility defaults preserve old adjoint-test checkpoints; current presets
  use kernel lifting. Combining those defaults would change old experiments.
- E1/E2 and E3 intentionally have different sampling rules. Tables keep their
  documented interpolation and first-crossing definitions.

## Remaining boundaries

`adaptive_run.cpp` and the audit CLI still coordinate several phases. Their
ordering determines checkpoint transactions, greedy training, lazy checks and
audit lifetimes. Future features should add narrowly scoped helpers beside the
existing diagnostic module rather than more nested policy branches. A generic
plugin or workflow framework is not needed for the current experiment set.

Some imported numerical kernels have long loops and compact formatting. They
are covered by mathematical regressions, and loop/reduction order affects
reproducibility. This review does not claim that every legacy function has been
rewritten. Further extraction should preserve arithmetic order and have a
specific ownership or dependency benefit, rather than merely increasing the
number of files.

Dense checkpoint payload publication remains synchronous. Full refinement
lineage scheduling, memory-aware audit admission and large high-k production
validation remain the separately documented future work.

## Validation

All five PNGs match fresh renders from the supplied original plotting functions
pixel for pixel with the same installed Matplotlib version. The E1/E2 scripts
save PDF before PNG; the E3 script saves PNG before PDF. Keeping this order also
avoids constrained-layout rounding differences. All scientific CSV/LaTeX tables
remain byte-identical to the preceding delivery.

The affected server regressions passed 6/6: scientific controls, asynchronous
restart, shared snapshots, incremental integration, checkpoint audits and source
integrity. Local Python tests cover the distinct E3 sampling, promotion-event
roundtrip, non-nominal residual rejection and terminal mesh selection. Evidence
is retained in `provenance/plot_style_validation.json`.

## Subsequent E3 display revision

At the user's request, E3 ALOD now shows only the initial state and completed
H-step endpoints, after applying the existing one-percent prefix cutoff. AFEM
and LOD sampling is unchanged. Intermediate terminal states remain in tables
and raw data, even when omitted from the curve. The earlier five-figure pixel
comparison documents the previous display policy; E3 intentionally differs now.
