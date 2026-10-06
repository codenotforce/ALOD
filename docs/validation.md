# Validation and provenance

Validation is scoped to named source/configuration/workload identities.
An accepted state, Ritz residual gate, successful bounded test and completed
production campaign establish different facts. See [experiments](experiments.md)
for frozen completion and [algorithm](algorithm.md) for mathematical differences.

## Run checks

```sh
ctest --preset release
python3 -m unittest discover -s tests -p 'test_*.py'
python3 tools/check_repository.py
python3 tools/check_repository.py --staged
```

| Area | Principal coverage |
|---|---|
| Mesh/FEM | Manufactured PDE/boundaries/quadrature, NVB conformity, ancestry and nested injection |
| LOD/Riesz/Theta | PG/local residuals, kernel constraints, mass conservation, explicit local inverses and dense spectra |
| Archive comparison | FEM/baseline/fixed-state fixtures, E1 prefix and E2 scoped comparisons with explicit policy |
| Enrichment | Regional support, deflation, coupled solves, POD gates and inherited/reset behavior |
| Concurrency | Actual team expansion both directions, exceptions, shared/disk parity, cancellation and immutable lifetime |
| Recovery/retention | Pause/resume, process interruption, checksums, corrupt quarantine, torn journals, validated receipts and geometry collection |
| Tools/reporting | Schemas, requested child thread metadata, campaigns, diagnostics and paper selection |

Key executable tests include `async_execution_context`,
`async_execution_equivalence`, `periodic_checkpoint`, `audited_state_retention`,
`asynchronous_restart`, `shared_audit_transport`, `fused_kernel_defect` and
`baseline_cli`. Python tooling discovery includes the metadata and repository
policy regressions. Fixed thread/BLAS settings in comparison tests are
repeatability choices, not production defaults. Full suite durations/results
are dated evidence, not a timeless count of current tests.

October 6 pre-publication verification rebuilt all default targets in local
WSL and exercised all 59 CTests. The full run passed 58 tests; the optimization
workflow test inherited the production preset's five-cycle restart interval
while expecting every state's checkpoint. Its test configuration now explicitly
requests every state, and the complete numerical/corruption-recovery comparison
passed on rerun. Python tooling discovery passed 46 tests. Worktree and exact
staged checks validated 330 text files, all 64 required imported hashes and the
canonical audit rows. Local Markdown targets and `git diff --check` passed.
These local checks did not restart production or rerun the full E1/E2/E3 campaigns.

Publication review found no absolute local paths in maintained Markdown/JSON,
no recognized credential signatures or generated binaries in the candidate
inventory. One local SSH alias was removed from historical access metadata;
numerical results were preserved. A public-file preview without the local plan
passed the same repository checker. Only subsequent commits exclude the plan;
existing Git history still contains earlier versions, personal author metadata
and the old SSH alias. Ordinary commits do not remove those historical records.

## Provenance policy

`provenance/source_manifest.json` separates `origin_sha256` (supplied origin)
from `sha256` (reviewed current content normalized to LF). Keep original identity
when editing imported code. Review intentional changes, retain the preceding
reviewed hash and update current hash plus modification description. Do not
silently refresh all hashes or skip mismatches merely to make a test pass.
Immutable historical numerical fixtures keep their acceptance identity.

The October 6 thread-team instrumentation in `src/lod/lod.cpp`,
`src/estimator/riesz.cpp`, `src/estimator/riesz_internal.hpp` and
`src/fem/batch.cpp` was reviewed: explicit parallel regions expose actual team
sizes while preserving loops, accumulation ordering and exception propagation.
Their current hashes are registered, and origin/prior hashes remain available.
The full checker still enforces each current imported hash and canonical audit
completeness (14,208 E1 and 1,650 E2 archived rows).

The checker validates portable text, guide language, local Markdown file targets,
source hashes, RHS identities and archived audit rows. Markdown targets must
belong to the publishable file inventory, so an ignored local file cannot hide
a broken link in a fresh checkout. Worktree deletions are
excluded from the live text inventory; deleting an imported manifest entry's
required file still fails. `--staged` checks exact index blobs. External URL
availability and Markdown heading anchors are not checked.

Baseline manifests now keep legacy BLAS_threads=null and a whitelist of actual
child thread-environment settings. Unset means runtime default; requested settings
are not measured backend team sizes. Tests cover inherited conflicting BLAS
requests and explicit OpenMP override without changing BLAS environment.

The eight core guides consolidate old phase/optimization/audit narratives.
Machine-readable provenance/fixtures, P0 guide, NOTICE and agent guidance remain.
The original Chinese plan is retained only locally and ignored by Git. Its
historical hashes remain in `external_evidence`, outside the required source
file list, so a public checkout does not need that private planning document.
Prior Markdown worktree contents were saved in an
ignored local ZIP before consolidation; committed historical text remains in Git.
One short `structural_optimization.md` redirect preserves the local plan's
historical hyperlink; it is not an additional maintained optimization guide.

## Benchmark tools

Use separate baseline/candidate builds and fresh outputs; report binary hashes,
configuration, order, repetitions, CPU/thread settings and numerical gates.

```sh
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 tools/benchmark_thread_budgets.py \
  --baseline "$BASELINE_EXE" --candidate "$CANDIDATE_EXE" \
  --output results/thread-comparison --threads 32 --repeats 3 --compute-only --level 10 --gap 4
```

This alternates original async, updated async and sequential internal execution
over two accepted states. E1/E2 use 16 training cases; E3 is nominal E1 at k=8.
The compute-only path excludes checkpoints/exports/audits and compares numerical
observables/decisions; full-vector and audit/recovery checks are separate.
The supervisor path is available without `--compute-only`.

| October 6 server case | Original async median | Updated async median | Internal overlap off |
|---|---:|---:|---:|
| E1 | 4.275 s | 4.225 s | 4.826 s |
| E2 | 17.511 s | 17.009 s | 22.824 s |
| E3 nominal k=8 | 3.824 s | 3.825 s | 3.974 s |

AMD EPYC 9554: 64 physical cores/128 logical CPUs; 32 requested threads,
level 10/gap 4/ell 2, rank cap 8 and three repetitions. Discrete decisions
matched, numeric gates passed and actual 32-worker teams were observed; CPU
activity averaged roughly 1600–1900%. This excludes production-scale audits/I/O
and shows no material E3 improvement. The server evidence lives under
`$SERVER_PROJECT_ROOT/thread-budget-20261006-1503`, including `compute32/report.json`,
`benchmark32-v2/report.json` and `benchmark-binaries.sha256`. Binary identities
precede the final exception-only guard-lifetime adjustment; six relevant final
server CTests subsequently passed. Existing production was not restarted.

Other manual tools are `tools/benchmark_async.py` (drain),
`tools/benchmark_spaces.py` (space phases), `tools/profile_fixed.py` (fixed-state)
and the optional `benchmark_theta` target (same-checkpoint defect execution).
Run each tool's `--help` for arguments. Older raw measurements are preserved in
`provenance/async_validation.json`, `provenance/deep_reuse_validation.json`,
`provenance/async_training_validation.json` and other dated provenance files.
Full dense publication, global hierarchy assembly and memory-aware audit
admission remain possible optimization work; bounded evidence is not a full
campaign speedup claim.

```sh
cmake --build build --target benchmark_theta
build/benchmark_theta "$CHECKPOINT" 128 32 parallel cold
build/benchmark_theta "$CHECKPOINT" 128 32 fused cold
build/benchmark_theta "$CHECKPOINT" 128 32 probe cold
```

`parallel` selects current parallel global actions; `fused` selects patch-local
fusion. `optimized` selects the earlier image-reuse/lightweight path. `probe`
compares actions on the same deterministic four-column block without running
the eigensolver. Cold modes use identical initial vectors; warm uses the saved
block. Restoration/action/eigensolve times have different scopes. Retained
operator images and exact-identity local reuse preserve the mathematical
operator; numeric gates still apply to every optimized comparison.

## Optional reference-solver comparison

```sh
cmake --build build --target benchmark_reference_solvers
build/benchmark_reference_solvers smoke 2 1 all
build/benchmark_reference_solvers "$CHECKPOINT" 32 1 energy 1e-14
```

This CHOLMOD-dependent optional target compares direct UMFPACK, GMRES/ILUT and
GMRES with one-level energy Cholesky preconditioning on checkpoint reference
meshes. It does not change production backends. Acceptance checks true
relative residual <=1e-9 and relative energy difference from direct <=1e-7;
an internal preconditioned stop flag is insufficient. Count setup/factor plus
solve, separate common preparation and checks, and measure multi-RHS reuse.

The September 27 E1 state-56 nominal-case observation had 667,556 free DOFs:
direct setup/solve 15.942 s; accepted energy-GMRES 69.307 s (233 iterations,
true residual 5.070e-11); ILUT failed after 400 iterations with residual 3.022e-2.
Concurrent production and one RHS limit interpretation. These results do not
compare every multigrid method or establish 48-RHS speedups.

## Legacy replay and storage checks

`tools/replay_legacy.py`, `tools/build_legacy_oracle.py`,
`tools/build_baseline_oracles.py` and related builders require supplied external
source and are explicit developer actions, never normal build steps. Do not
combine mismatched old headers/libraries or regenerate fixtures to hide failure.
The [P0 guide](provenance/p0/README.md) describes their source/replay provenance.
The optional `stress_checkpoint` target exercises large-file serialization,
checksumming and loading, not a million-node Helmholtz solve.
