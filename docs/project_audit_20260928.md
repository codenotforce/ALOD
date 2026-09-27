# Project source and documentation audit — 2026-09-28

## Scope and method

This is a local working-tree audit, including tracked files and non-ignored new
files. The baseline inventory contains 340 files before the two documents added
by this audit. It is not an audit of server deployments or ignored experiments.
Build registration, direct header includes, Python imports/entry points, presets,
CTest registration, documentation claims and recorded provenance were inspected.
No numerical implementation, fixture, experiment or server file was removed.
Static references establish a known use; they do not prove every function is
reachable for every input. This is not a new full numerical validation campaign.

## Documentation correction status — 2026-09-28

The findings below preserve the audit-time observations. Documentation fixes
have now been applied in place; their original wording is not a current backlog.

| Finding | Resolution |
|---|---|
| A1 | README, migration status and experiment contracts distinguish completed frozen deployments, historical horizons and later-optimization validation. |
| A2 documentation | Old campaign instructions are explicitly historical and warn about mutable presets; current production links point to `make_production_campaign.py`. The old generator code is unchanged and remains unsuitable for a frozen replay claim. |
| A3 | Obsolete feature/default claims were replaced with current behavior or marked as dated measurements beside the relevant text. Implemented diagnostics/reuse/integration are no longer listed as unimplemented. |
| A4 documentation | P2 and P5/P6 profiler instructions state no affinity or forced single-thread BLAS. Old pinned timings remain historical. Unused profiler code/metadata is unchanged and remains a separate cleanup item. |
| A5 | Historical/current inventories, E3 cycle sampling, generic table windows and the scoped 25-cycle/28-cycle export are explicitly distinguished. |
| A6–A7 | Retention and generated-tree guidance remains valid; no implementation or result was deleted. |

Only documentation and its README navigation were changed. Numerical code,
configuration presets, immutable provenance and fixtures remain untouched.

## Original findings and recommended disposition

### A1. Current-status prose contains obsolete stage conclusions — high priority

`README.md` still describes full current E1/E2 campaigns as pending and links a
four-wavenumber protocol as if it were the current campaign. The dated body of
`docs/migration_status.md` says no new production trajectories have been created;
`docs/experiments.md` ends with the same historical claim and describes the old
51/33-state contracts. These are valid stage records but misleading current
summaries. Completed production now includes 28-cycle E1/E2 (57 accepted states
each) and E3 k=8,16,32,64,128. This does NOT mean those experiments used every
subsequent optimization in the current working tree: the latest Theta changes
were validated separately on a saved final checkpoint.

Recommendation: make the latest production protocol and this structure guide
the navigation entry points; label old stage summaries explicitly as historical.
Preserve old measured counts and numerical evidence instead of rewriting history.

### A2. The old campaign generator mixes old metadata with current presets — high priority

`tools/make_kappa_campaign.py` still describes k=8,16,32,64, old memory limits,
fixed-ell controls, and `e2_threshold=0.1` metadata. It loads today's mutable
`e1_main.json` / `e2_main.json`, whose production relative threshold is 0.3 and
which now carry absolute thresholds and later workflow fields. Thus it is neither
the current ordered campaign nor a frozen reproduction of the old campaign.
`apps/alod_geometry.cpp` deliberately omits k=128 unless `--paper-v4` is passed;
this is an old-mode distinction, not missing k=128 support.

Use `tools/make_production_campaign.py` for the current ordered 32+32-thread
campaign. Keep the older generator only as a historical recipe until it is
explicitly retired or given immutable historical inputs. Do not silently replace
its outputs or delete its evidence. `tools/make_paper_v4_campaign.py` is a valid
seven-job preparation utility, but does not implement the requested dependency
barriers or explicit production audit budgets; it is not the production entry.

### A3. Earlier algorithm/default descriptions are superseded — high priority

- `docs/p3_p4.md`, Localization policies: raw E1/E2 thresholds and one-promotion
  wording describe older controller policies. Current main presets use kernel
  lifting, the solution-scaled relative gate 0.3, and the lazy absolute gate
  tau(k)=3.2/k. Compatibility presets intentionally remain different.
- `docs/async_training_optimization.md` says persistent workers and shared
  in-memory snapshots are not implemented. Its follow-up banner supersedes that
  paragraph; `apps/audit_worker.cpp`, `apps/shared_audit.cpp` and
  `tools/native_audit.py` implement those paths.
- `docs/structural_optimization.md` still contains zero-cache-default and pending
  reuse/worker/integration statements. `docs/deep_reuse.md` supersedes much of
  this, and `docs/scientific_workflow.md` supersedes later incremental-integration
  and checkpoint-preparation items. Current adaptive corrector and Riesz cache
  defaults are 64 MiB each, with explicit disable/fallback options.
- `docs/repeated_work_review.md` lists reduced-operator persistence, shared
  geometry and source-moment reuse as opportunities that later work implemented.
- `docs/paper_implementation_map.md` still lists the optional diagnostic menu as
  pending; its source names are explicitly names in the external archive, not
  missing files to recreate in this repository.

Recommendation: retain these pages as dated design/measurement records. Move
current behavior into a small authoritative guide and put supersession notices
next to obsolete operational claims, rather than relying only on late appendices.

### A4. Profiling documentation still promises CPU pinning — medium priority

`docs/p2.md`, Validation and server profile, says `tools/profile_fixed.py` pins
CPU affinity and fixes BLAS to one thread. The current script explicitly does
not set affinity; `--cpus` is descriptive only and execution uses the shared
runtime environment helper. Its `SYS_CPU` constant is unused, and
`physical_cores` is populated from an always-empty local list. These are small
code/metadata leftovers from the former affinity profiler, not an active pinning
restriction or a 16-thread cap. The old pinned measurement remains historical
measurement evidence, but its instructions no longer describe the current tool.

Recommendation: update the operational paragraph and either remove obsolete
metadata or implement explicitly descriptive topology observation. No change to
numerical execution is needed to resolve this documentation inconsistency.

### A5. Paper inventories and plotting scopes need explicit version labels — medium priority

`docs/provenance/paper_inventory.json` and the early implementation map record a
seven-figure/eight-table manuscript. `configs/paper_v4_inventory.json` defines the
current five-figure/five-table delivery. Both are intentional evidence; they are
not interchangeable current inventories. `docs/code_structure_review.md` first
states all-state E3 plotting, then supersedes it with completed-cycle plotting.
The actual selector is `tools/paper/selection.py` and its use in `figures.py`.

`tools/paper_delivery.py` builds the broader paper package; its table fitting
windows include a five-cycle E1 tail. `tools/paper_snapshot.py` is a deliberately
specific export for the latest request: E1/E2 plots through cycle 25, matching
cycle-25 meshes, and separate ten-point fits on cycles 19–28 of the full 28-cycle
runs. They are not duplicate interchangeable CLIs. The latter combines current
ALOD with explicitly labelled historical controls and is not a new control run.

### A6. No whole compiled source file is proven abandoned

All 29 numerical translation units under `src/` are registered in `alod_core`.
All application translation units are attached to executable targets. Every
maintained C/C++ header has at least one direct include in maintained C/C++ code.
The CMake `foreach` registers several tests indirectly; absence of their literal
filenames from CMake is not evidence of dead code. `tests/legacy_p2_probe.cpp`
is deliberately compiled by `tools/build_p2_oracles.py`, outside normal CMake.

Do not remove these merely because their names look old:

- `src/enrichment/aot.cpp`: old adjoint-test compatibility/diagnostic path;
  kernel lifting is the current production choice, not the only supported one.
- `src/problems/legacy_cases.cpp` and lower-level `include/helmholtz`, `lod`,
  `mesh` headers: active extracted numerical implementations, not stale copies.
- `apps/alod_audit.cpp`: small standalone frontend reusing the same audit
  implementation as `alod_run`; it is not a second independent audit algorithm.
- `tests/fixtures`, `docs/provenance`, oracle builders and replay tools: frozen
  regression/provenance evidence and explicitly manual reconstruction utilities.
- `apps/alod_stability.cpp`: optional diagnostic remains useful even though the
  earlier calibration experiment was cancelled.
- The optional fused Theta path: tested alternative, not the fastest default on
  the measured final k=128 state; retained for comparison and different workloads.

This conclusion is at file and documented-interface level. Removing apparently
unused individual functions requires a more specific API/dead-code review.

### A7. Generated trees and older result bundles are not maintained source

`build/`, `build-v4/`, `_local/`, `results/`, Python caches and Git internals are
excluded from the maintained file catalogue. They contain binaries, temporary
scripts, downloaded records and successive paper/benchmark outputs. Several
paper bundles repeat renderers and canonical data intentionally for portability.
The newest 25-cycle figure deliverable is `results/paper-figures-25cycles-20260928`;
older bundles are historical, not instructions for the current figure cutoff.

Do not equate an ignored tree with safe-to-delete data. Results/checkpoints may
be unique; archival and retention decisions are separate from this audit.
A local build directory also does not establish that its binaries include the
latest uncommitted source changes.

## Maintenance priorities

1. Documentation separation for A1–A5 is complete. Keep current policy changes
   synchronized with operational guides; retain immutable fixtures and provenance.
2. Retire or freeze the inputs of the old campaign generator before using it for
   new experiments; use the explicit ordered production generator meanwhile.
3. Remove the small descriptive profiler leftovers in a focused cleanup.
4. If further modularization is needed, split scheduling/publication concerns
   out of `apps/adaptive_run.cpp`, and selected-region training/coordinate/POD
   concerns out of `src/enrichment/regional.cpp`. `src/estimator/riesz.cpp` can
   separate patch preparation from application. File length alone is not a reason
   to rewrite stable numerical mesh and quadrature kernels.
5. Record a source/version identity with each future benchmark or production
   launch. The working tree includes changes newer than the last Git commit;
   commit identity alone currently does not capture all validated optimizations.

## Validation of this audit

The accompanying structure document has one entry for every maintained file,
including new non-ignored files and both audit documents. Relative file links
and inventory coverage are checked. Repository portability/provenance validation
is run after documentation changes. No full solver experiment is needed for
this read-only source audit and documentation addition.

Completed checks: 342 unique catalogue entries match all maintained working-tree
files exactly; all catalogue/local documentation links resolve; all 29 core
translation units are registered. `tools/check_repository.py` validated 342 text
files, 65 imported hashes, 14,208 E1 and 1,650 E2 archived audit rows with no local
absolute paths. `git diff --check` passed. Numerical C++ code and frozen fixtures
were not changed, and no new numerical performance claim is made by this audit.

The documentation follow-up checked 19 revised pages, 432 resolving local links
and exact coverage of all 342 maintained files. Repository integrity again
validated 65 imported hashes and the unchanged historical audit rows; English
and no-local-absolute-path checks and `git diff --check` passed. The old campaign
generator and profiler code leftovers remain explicitly separate code tasks.
