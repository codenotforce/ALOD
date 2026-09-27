# Experiment contracts

**Scope:** the contracts in the table below are the original imported experiments. Current production uses E1/E2 28-cycle runs (57 states each) and E3 k=8,16,32,64,128; see [the production protocol](production_campaign_20260927.md). All seven production jobs and audits completed on their frozen deployments. Historical control results and the latest optimized working tree must not be conflated.

The full migration specification is `ALOD_SUBPROJECT_AGENT_PLAN_20260909.md`. This document identifies effective experiment roles; frozen JSON and CSV evidence is under `provenance/` and `tests/fixtures/legacy/`.

| Experiment | Accepted scope | Marking and diagnostics |
|---|---|---|
| E1 production | 51 accepted states, 48 fixed loads, kappa=16 | Nominal adaptation, periodic reference-only steps, lazy global ell policy; retain configured extra checks and endpoint check. |
| E1 controls | Nominal/family ALOD: 51 states each; nominal/family AFEM: 97 states each | Fixed ell=3; nominal sample 0 or mean of 16 training members marks. All 48 members are audited; test and shift never mark. |
| E2 production | 33 accepted states; 48 formal plus pure-corner and pure-wave audits | Regional AS, normalized training selection, POD/AOT, adaptive ell and inheritance. Only accepted artifacts define the canonical trajectory. |

E1 has 16 train, 24 test and 8 shift samples. Imported CSV uses `shift-audit`; the portable parameter table calls the same role `shift`. E2 uses the same formal role counts and two additional diagnostic loads. Sample identities, coordinates, phase and amplitudes are frozen in `data/rhs/`.

Historical flags encoded in build-hash substrings are part of the effective configuration. The P0 branch inventories enumerate literal `wf_has` tests and environment lookups with source locations; replay configurations preserve the original production hash. E1 controls additionally use `E1_FAMILY_CONTROL`, `E1_FAMILY_MARKING` and `E1_CONTROL_OUTPUT`; these control nominal/family marking and artifact destinations. P1 baselines use explicit validated options and do not interpret those hidden controls.

Nominal baselines use AFEM theta=0.15 with exact-error stopping targets 0.008 for E1 and 0.003 for E2, UFEM uniform NVB, and SLOD ell=3 with four additional reference levels. Exact error may stop these baselines but is never a marking indicator. Historical baseline manifests, including their quadrature and horizons, are the source of the full presets. Small smoke runs only establish numerical equivalence.

The raw floor variable F and normalized reference error are distinct quantities. Do not relabel them or mix pre-acceptance diagnostics with accepted errors. Keep production adaptation, repeated-RHS audit, training, and timing roles separate. The original contracts above predate the completed September 27 server campaign. They are retained for interpreting imported evidence, not as current run horizons.
