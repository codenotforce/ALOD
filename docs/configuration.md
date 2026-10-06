# Configuration

Start from a complete checked-in preset. The Python validators reject unknown
fields, malformed/nonfinite values and booleans used as integer counts.
Only stated compatibility fields are filled when absent; default dictionaries
are not a promise that arbitrary partial JSON is accepted.
See [algorithm](algorithm.md) for mathematical meaning and [runtime](runtime.md)
for operational effects.

## Presets and common fixed/adaptive fields

Baseline presets live in `configs/smoke/` and `configs/full/`; fixed probes in
`configs/p2/`; adaptive/main/control/E3 presets in `configs/adaptive/`.
Campaign generators freeze copies. Editing a repository preset does not mutate
an existing campaign or its numerical history.

| Field | Contract |
|---|---|
| `problem`, `wavenumber` | E1/E2 and positive finite k; omitted wavenumber is 16 |
| `level`, `gap`, `ell`, `graded` | Coarse NVB level 0..12, reference gap 1..8, ell 1..4, boolean graded probe |
| `interpolation`, `riesz_patches` | `area`/`arithmetic`, `n2`/`archive`; inspect presets explicitly |
| `member_ids`, `training_ids` | Nonempty distinct existing IDs; training subset must be included and have train role |
| `theta` | Coarse marking fraction in (0,1] |
| `ritz_tolerance`, `ritz_iterations`, `dense_threshold` | Positive tolerance <=1, iteration cap 1..10000, dense cutoff 0..512 |
| `threads` | Integer 0..2147483647; zero inherits runtime, positive requests a budget |
| `maximum_nodes` | Fixed probe 1..200000; adaptive 1..4000000 |

The fixed-state defaults in `tools/run_fixed.py` are level=2, gap=3, ell=1,
area/n2, threads=0, theta=0.15, tolerance=1e-9 and 16 training cases.
Formal E1/E2 tables have IDs 0..47: 16 train, 24 test, eight shift; ID 0 is
nominal. Adaptive configurations must include ID 0. Diagnostic E2 pure cases
48/49 are separate from the formal table.

## Adaptive additions

| Field | Contract |
|---|---|
| `method` | ALOD/AFEM; AFEM requires `m_ref=1` |
| `cycles`, `m_ref`, `reference_theta` | Cycles 0..1000, sweeps 1..16, exactly one positive fraction <=1 per sweep |
| `state_limit`, `minimum_gap` | Accepted-state cap 0..10001 (zero uses horizon); reference reserve 0..8 |
| `ell_mode`, `maximum_ell`, `extra_checks` | lazy/every/fixed; cap between initial ell and 4; distinct state IDs 0..10000 |
| `ell_ratio_mode`, `ell_threshold` | raw/solution_scaled; nonnegative finite relative threshold |
| `ell_absolute_threshold` | Nonnegative finite manual scalar; -1 disables the absolute gate |
| `force_promotions` | Explicit E2 ALOD diagnostic states; forbidden in fixed mode; production uses empty list |
| `enrichment_tests` | kernel_lift/adjoint; changing this changes mathematics |
| `radius`, `rank_cap`, `inherit` | Nonnegative finite region radius, rank cap 1..24, boolean dictionary transfer |
| `exact_target`, `exact_scope` | Positive exact relative-error target or -1 disabled; nominal/training_max |
| `audit`, `audit_mode`, `emit_solution` | Boolean audit, full/exact metrics, boolean full solution emission |
| `maximum_patch_entries`, `maximum_dense_entries` | Positive allocation ceilings, at most 64000000000 and 1000000000 respectively |
| `checkpoint_interval_cycles` | Integer 0..1000; zero uses compatibility boundary publication; positive schedules cycle/final restart points |
| `checkpoint_retention`, `keep_final` | all/delete_after_audit/periodic; final preservation boolean for delete_after_audit only |

Both reclamation modes require `audit=true`; periodic requires a positive
checkpoint interval and always retains the final accepted state.
See runtime for validation receipts, pauses and failure preservation.

Compatibility fill values are interval=0, retention=all, keep_final=false,
exact_target=-1, exact_scope=nominal, audit_mode=full, absolute threshold=-1,
wavenumber=16, ratio_mode=raw, relative threshold=0 and enrichment_tests=adjoint.
This last default preserves old checkpoints; new production presets explicitly
select kernel_lift, solution_scaled and their actual thresholds. Do not treat
compatibility defaults as today's production policy.

## Baseline schema

`tools/run_baseline.py` requires schema_version=1, problem, method (AFEM/UFEM/SLOD),
initial_level (0..24), states (1..10000 including state zero), theta (0,1],
target (nonnegative exact relative energy error; zero disables), maximum_nodes
(4..2000000000), threads and boolean emit_solution. Only reference_gap (1..8,
default 4) and wavenumber (default 16) may be omitted. SLOD adds a conservative
reference allocation guard. `run.json` now records the child thread environment;
legacy `BLAS_threads` is null because actual BLAS team size is not measured.
Unset environment values mean runtime defaults, not one thread.

## Operational CLI and environment

The adaptive runner supports `--resume auto|CHECKPOINT`, `--pause-state`,
`--pause-phase accepted|training|ell`, `--audit-workers`, `--audit-threads`,
`--audit-drain-workers`, `--checkpoint-retention`, `--keep-final` and
`--no-keep-final`. Use `--help` for complete CLI spelling. Thread selection
remains runtime-inherited or explicit; hardware sizing and CPU reservation
were not added. BLAS settings remain user-controlled.

| Environment setting | Current use / default |
|---|---|
| `ALOD_PATCH_CACHE_BYTES`, `ALOD_RIESZ_CACHE_BYTES` | Each 64 MiB retained local-result budget; zero disables |
| `ALOD_SOURCE_MOMENT_BYTES` | 768 MiB source-moment budget; zero disables |
| `ALOD_TRAIN_RIESZ_BYTES`, `ALOD_POD_GRAM_BYTES` | Each 256 MiB estimated workspace allowance; zero disables optimized path |
| `ALOD_SHARED_SNAPSHOT_BYTES` | 256 MiB shared payload/hierarchy allowance; zero exercises disk fallback |
| `ALOD_AUDIT_HIERARCHY_BYTES` | 64 MiB persistent-worker identical-hierarchy allowance |
| `ALOD_AUDIT_DISK_ONLY=1` | Resident separate processes instead of shared adaptive memory |
| `ALOD_AUDIT_ONESHOT=1` | Previous process-per-input audit mode |
| `ALOD_REFERENCE_EXECUTION=1` | Numerical comparison paths and separate audit workers |
| `ALOD_ASYNC_DISABLE=1` | Sequential within-state execution for comparison |
| `ALOD_INTEGRATION_COLD=1` | Disable incremental integration for comparison |
| `ALOD_CHECKPOINT_SYNC=1` | Synchronous geometry publication for comparison |
| `ALOD_TIMING_FILE` | JSONL timing/team destination; supervisor supplies it |

Budgets above are retention/workspace allowances, not whole-process RSS limits.
Live factors, temporary buffers and each auditor's workspaces add memory.
Process affinity, scheduler/cgroup restrictions and backend algorithms also
affect observed teams and utilization. Test/oracle thread settings are scoped
repeatability choices and do not define production defaults.
