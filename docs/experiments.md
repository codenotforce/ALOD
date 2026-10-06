# Experiments and frozen evidence

This guide owns the current ordered campaign and the distinction between
historical evidence and current source. See [configuration](configuration.md)
for fields and [runtime](runtime.md) for audit/restart/retention behavior.

## Current ordered E1/E2/E3 protocol

```sh
python3 tools/make_production_campaign.py --output results/production-campaign
python3 tools/run_campaign.py --campaign results/production-campaign --build build
python3 tools/run_campaign.py --campaign results/production-campaign --build build --resume
```

Run the long supervisor in tmux. Use `--start-at-e2` on the generator for a fresh
replacement campaign preserving completed E1/earlier output directories.
The dependency order is E1, E2, E3 k=8/16 together, then k=32,64,128 sequentially.
A prerequisite completes only after solver, audits and exports succeed; failure
blocks dependents. Resume retries unfinished jobs and retains successful jobs.
Campaign memory/disk guards can defer/stop admission; they never silently reduce
thread counts or label a failed job complete.

Every job requests 32 adaptive threads and one 32-thread auditor, including
drain. The simultaneous k=8/16 stage requests up to 128 compute threads overall.
No launcher CPU affinity is imposed. E1/E2 use two reference sweeps per cycle,
28 cycles and 57 accepted states including zero; exact-target stopping is off.
Production presets explicitly select manuscript area/n2, kernel_lift,
solution_scaled threshold 0.3, lazy timing and the manual absolute rule 3.2/k.
The proposed gamma calibration was cancelled; this rule is prescribed policy.

| k | Absolute threshold | E3 cycles | State cap / maximum states |
|---:|---:|---:|---:|
| 8 | 0.4 | 27 | 55 |
| 16 | 0.2 | 25 | 51 |
| 32 | 0.1 | 22 | 45 |
| 64 | 0.05 | 20 | 38 |
| 128 | 0.025 | 100 | 30 |

E3 uses nominal exact relative-energy target 1% and exact-only audits.
A horizon/resource stop without that accuracy is not target attainment.
Reference FEM audits use the UMFPACK 64-bit interface and METIS ordering for
large systems, retaining the relative-residual acceptance gate.

The checkpoint interval is five completed cycles plus final state. E1/E2
scheduled restart states are 10,20,30,40,50,56; no initial point unless terminal
or explicitly paused. Frozen deployments used retention=all. Other current
retention modes do not change those already frozen configurations or directories.

## Recorded completion

September 27 frozen deployments completed E1/E2 at 57 states each, followed by
E3 k=8,16,32,64,128 at 54,50,41,36,30 states respectively, with their audits.
E3 reached the nominal target. E1 belongs to the original campaign; successful
E2/E3 belong to the replacement campaign, not the earlier failed E2 attempt.
Later Theta/async/thread changes have separate checkpoint/bounded validation.
These records do not establish a full campaign of the latest working tree.

## RHS cases and controls

Formal E1/E2 tables contain 48 fixed cases: 16 train, 24 test, eight shift;
nominal ID is zero. E1 nominal/family ALOD and AFEM controls use the same table
and explicit marking set: nominal [0], family 0..15. Fixed-ell control policy,
inherit/reset, lazy/every and archive/manuscript modes must be labelled separately.
Audits do not drive adaptation. E2 optional pure cases 48/49 use coefficients
(c,a)=(1,0)/(0,0.5) at the nominal wave location/phase.

| Historical record | Scope |
|---|---|
| Imported E1 production | 51 accepted states; original policy/source identity |
| Imported E1 controls | 51 states per ALOD arm, 97 per AFEM arm; frozen canonical tables |
| Imported E2 production | 33 accepted states and pure-component diagnostics |
| Early September 22 E2 32-thread deployment | Earlier horizon and binary, superseded by final 28-cycle deployment |
| September 12 four-wavenumber recipe | Superseded launch design, not the current ordered campaign |

`tools/make_kappa_campaign.py` is historical tooling that reads mutable current
presets while retaining old metadata; it is not a frozen replay tool. Use the
production generator for new runs. The [P0 guide](provenance/p0/README.md)
preserves acceptance context. The original planning document is operator-local
and is not required to run or validate a public checkout.
Machine-readable contracts and source identities remain under `provenance/`.
Historical curves/fixtures must not be regenerated to hide a failed comparison.

## Results and reproducibility

Keep `run.json`, solver/accepted journal, effective configuration/RHS table,
audit metrics/completion records, hashes and timing logs with a result package.
Rerendering metrics, reauditing a frozen space and resuming adaptation need
different retained data. Baseline target/horizon presets are historical controls,
not automatic replacements for current adaptive experiment definitions.
Architecture/compiler rounding can change near-tied exact marking order;
archive equality is tested with matching policies/toolchains, not asserted
across arbitrary machines. No historical dependency lock guarantees bitwise
toolchain reconstruction.
