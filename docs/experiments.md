# Experiments and frozen evidence

This guide distinguishes the version 7 manuscript experiments, checked-in
campaign recipes, and frozen historical records. See [configuration](configuration.md)
for fields and [runtime](runtime.md) for audit/restart/retention behavior.

## Version 7 manuscript scope

The checked sources are `../LOD_paper/ALOD/version7/arxiv.tex` and `main.tex`;
their numerical sections agree. The paper uses descriptive experiment titles:
**A reentrant corner** (repository/server E2) and **Wavenumber dependence**
(E3, using the E1 problem). E1 remains supported in the code and historical
records, but is not a separate experiment in this manuscript.

| Setting | A reentrant corner | Wavenumber dependence |
| --- | --- | --- |
| Load scope | 16 training loads; 24 independent test loads | One fixed-center load for each k |
| Initial maximum element diameters | \(H_0=\sqrt2/8,h_0=\sqrt2/32\) | \(kH_0=2\sqrt2,k^3h_0^2=8\) |
| Initial / maximum ell | 2 / 4 | 2 / 5 |
| Cycle limit | 28 | 40 |
| Absolute / balance thresholds | 0.2 / 0.3 | \(3.2/k\) / 0.3 |
| Regional enrichment | \(D=\Omega\cap B_{0.6}(0)\); final rank 2 | None |

Both columns use \(\theta_H=0.15\), \(m_{\rm ref}=2\), fine marking
fractions \((0.3,0.2)\) and \(d_{\min}=10^{-12}\). The paper's cycle
description differs from the code's accepted-sweep convention; see
[algorithm](algorithm.md).

The corner parameter is \(\mu=(\mu_c,\mu_\alpha,\mu_x,\mu_y,\mu_\varphi)\),
with reference \((1,1/2,-1/2,1/2,0)\). It varies
\(\mu_c\in[0.25,1.75]\), \(\mu_\alpha\in[0.15,1]\),
\(\mu_\varphi\in[-\pi,\pi]\), and the center in the radius-0.2 disk
about \((-1/2,1/2)\). Training uses the reference and 15 randomized loads;
the 24 test loads do not drive training or marking. The manuscript does not
specify a seed or generator; read those from the frozen run inputs rather than
assuming the checked-in RHS table generated these results.

The independent single-load/family-aware controls each run 25 cycles at fixed
\(\ell=3\). Their final family-aware mesh pair is also used for the fixed-mesh
Base / Retrained shared / Adaptive shared / Per-load fresh controls; it is
distinct from the main 28-cycle endpoint. Retrained shared is rebuilt from the
whole training family for each ell, Adaptive shared freezes the adaptive space,
and Per-load fresh uses each evaluated load to rebuild its own space. The latter
is a diagnostic, not shared-space reuse. The reported test-load maxima at
\(N_{\rm on}=10^3\) use logarithmic interpolation.

Wavenumber dependence fixes \(\mu_0=(3/4,1/2)\), uses separate operators at
\(k=8,16,32,64,128\), and reports \(E=E_{\mu_0}\). ALOD reaches 1% at all
five wavenumbers; the LOD k=128 table entry is its final memory-limited state
at 1.1536%, not a successful target crossing. These statements summarize the
manuscript; this documentation update does not establish a fresh server audit.

## Checked-in ordered E1/E2/E3 recipe

The following generator retains the earlier ordered recipe. Its E3 horizons,
ell caps, RHS tables and dependency barriers are not a version-7 reproduction
specification. Freeze and verify explicit run settings for a new experiment;
do not infer the latest server deployment from this repository recipe.

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

E3 uses reference-load relative energy error target 1% and exact-only audits.
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
E3 reached the reference-load target. E1 belongs to the original campaign; successful
E2/E3 belong to the replacement campaign, not the earlier failed E2 attempt.
Later Theta/async/thread changes have separate checkpoint/bounded validation.
These records do not establish a full campaign of the latest working tree.

## RHS cases and controls

Formal E1/E2 tables contain 48 fixed cases: 16 train, 24 test, eight shift;
reference ID is zero (literal `nominal`). E1 single-load/family-aware ALOD and AFEM controls use the same table
and explicit marking set: single load [0], family 0..15. Fixed-ell control policy,
inherit/reset, lazy/every and archive/manuscript modes must be labelled separately.
Audits do not drive adaptation. E2 optional pure cases 48/49 use coefficients
(c,a)=(1,0)/(0,0.5) at the reference wave location/phase.

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
