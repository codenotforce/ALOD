# ALOD Subproject Migration Plan

Original date: 2026-09-09. Paper verification: 2026-09-11. English edition: 2026-09-11. The filename is retained for existing links. Execution status and measured results belong in [docs/migration_status.md](docs/migration_status.md); this document describes the complete migration contract, not completed implementation.

Paper source: `<PAPER_ROOT>/helmholtz_lod_certified_amsart_revised.tex`, SHA-256 `f714fb0ba66c2cc276437002aa2306d16271fd91fc88205ff8c3ca07ccdcbc37`. Recheck it before executing a later phase. If it changes, update the difference inventory before using old figure numbers or tables.

The target is a small, independently buildable ALOD project extracted from `LOD2d_C++`: E1 production, all four E1 nominal/family marking controls, E2 regional additive-Schwarz enrichment with inheritance, and AFEM/UFEM/SLOD baselines. Family-AFEM is an AFEM marking-member configuration, not another algorithm or problem. Do not build a general FEM or Helmholtz platform.

The September 11 paper requires the four E1 marking controls and complete state-by-load audits. The removed E2 shared/fresh figure and five-sample table are optional historical diagnostics. Preserve their evidence without automatically restoring them to the paper.

Use `<ALOD_ROOT>`, `<LEGACY_ROOT>`, `<PAPER_ROOT>` and `<REMOTE_ALOD_ROOT>` as symbolic roots. Committed files must not contain machine-specific absolute paths. Normal builds and runs must not reference legacy headers, libraries, build directories or executables.

## 1. Execution contract

Follow P0–P9 dependencies. Freeze sources and obtain comparable old outputs before extracting numerical modules; preserve mathematical behavior before optimizing. Record each phase's code, actual commands, inputs, hashes, validation results, resource observations and remaining work. Continue from `docs/migration_status.md` across sessions.

- Preserve the original project, manuscript and all original experiment data. Work in the new repository; never overwrite a running server deployment.
- Import necessary source, configuration, explicitly listed regression inputs and provenance. Do not copy whole result trees, executables or static libraries into the release project.
- Replace behavior selected by directory names, `build_hash` substrings or `diagnostic_sc_lod_campaign` with explicit, typed and validated configuration.
- Do not change local operators, AOT, normalization, tie ordering, training/compression thresholds or quadrature to obtain speedups. Changes to the mathematical method require separate comparisons.
- The original planning task did not authorize experiments. Subsequent implementation tasks authorize their stated migration stages; full production campaigns and server occupancy remain separate from local smoke validation.
- On a source/paper mismatch, record the exact difference and a minimal reproducer. Do not silently choose a definition and claim complete reproduction.

Final acceptance requires a clean WSL/Linux checkout to build and run independently; E1/E2, three baselines and four E1 controls; independent diagnostics from saved states; regeneration of the current figures/tables in section 11.3; real resume; resource protection and auditable timings. Assess online dimension, equal-budget error and wall time separately.

## 2. Authoritative sources and version issues

Paths in the following table are relative to `<LEGACY_ROOT>` unless marked as manuscript inputs. Original source evidence is retained under `docs/provenance/`. Historical evidence identifies files; it is not automatically a complete dependency closure.

| Source | Purpose and restrictions |
|---|---|
| Manuscript and actually included `figures/` and `tables/` | Follow LaTeX labels and includes. Unreferenced files are not required deliverables. |
| `results/E1-periodic-lazy12-51-20260909/data/config.json` and `runs/` | Final 51-state E1 configuration, meshes, lazy checks and errors. Some behavior is in the production driver, not the JSON. |
| `results/E1-production51-paper-20260909/` | E1 paper numbers, plotting inputs and baseline audit. Fit complete-cycle endpoints. |
| `results/E1-family-marking-controls-20260910/CONFIGURATION.md`, `configs/`, `data/`, `paper-update/summary.json` | Fixed-ell=3 controls, Figure 3 and matched budgets; 51/51/97/97 states and 48 loads each. Do not replace variable-ell production or the original 87-state AFEM. |
| Same campaign's `reproducibility/production-source/`, build evidence and recovery evidence | Actual production sources. Ordinary AFEM has a complete audit. Recovery executables are archival replay tools, not new defaults. |
| `results/E2-regional-adaptiveell-inherit-controls-20260909/deployment/source/` and `deployment/configs/inheritell.json` | Authoritative E2 production snapshot. It differs from the main working tree. |
| Same campaign's `outputs/inheritell/` | Use `accepted_control_*.csv`; unprefixed control files may describe pre-promotion solutions. |
| `results/E2-AOT-factor-reuse-20260909/` | Verified AOT factor cache, runner and includes. Previously linked snapshot libraries; extract and rebuild all dependencies. |
| `results/E2-production-inherit-paper-20260909/` | Production/baseline/frozen-space evidence and 75 historical independent-load audits. Old timing excludes the later AOT cache. Removed shared/fresh artifacts are optional. |
| `results/E2-paper-revision-20260909/CHANGELOG.md` | Historical revision only; the current manuscript takes precedence. |
| `results/E2-regional-budget-inherit-controls-20260909/analysis/REPORT.md` | Budget/inheritance costs and limitations. Distinguish global and regional budgets. |
| `results/E2-regional-adjacent-fast-20260909/README.md` | Existing parallelism/caches and adjacent-level diagnostics; no expensive global-corrector probe is required. |

The working-tree `include/helmholtz/adaptive/kernel_residual.h` lacks production `AdditiveKernelRieszContext::apply_selected` and `selected_values`. Copying a current runner with that header is not a complete regional-AS migration.

P0 must record each imported file's original relative path, hash, snapshot, copyright/license evidence and modification rationale. Record compiler, CMake, Eigen, SuiteSparse, BLAS and OpenMP versions, and compile/link commands with symbolic roots. `git_commit=WORKTREE` is insufficient. If E1 source is missing locally, retrieve its archived server deployment read-only and verify source and binary provenance.

The earlier `docs/E2_regional_AS_AOT_control_plan_20260909.md` contains obsolete probes/radii; treat it as history.

### 2.1 Paper-to-implementation map

P0 records labels, source functions, configuration fields and regression inputs in `docs/paper_implementation_map.md`.

| Labels | Required correspondence |
|---|---|
| `eq:local-kernel-space`, `eq:as-local-riesz` | Geometric coarse vertices, `W_h(N^2(z))`, positive local Riesz problems, constraints, multiplicity and indicator redistribution |
| `eq:certifier`, `eq:localized-coarse-smallness`, `eq:reference-interval` | Operator-norm Theta and Ritz approximation; distinguish theoretical solvability constants from practical promotion thresholds |
| `eq:regional-as`, `eq:mode-deflation` | Whole-patch regional selection, unweighted additive sum and two-pass energy deflation; distinguish kernel representatives from projected trial modes |
| `eq:aot-matrix`, `eq:coupled-block-system` | Global inverse-adjoint response, retained base tests, block identities and numerical residuals |
| `eq:regional-training-target`, `eq:pod-dictionary`, `eq:distance-gate`, `eq:residual-gate` | Frozen base targets, normalized correction-coordinate POD and actual per-member coupled solves |
| `alg:singular`, `alg:smooth` | Fixed-mesh solve/train/promote and `m_ref` reference sweeps |
| `fig:e1-family-marking`, `tab:e1-family-control` | Four fixed-ell controls, per-sample matched-budget interpolation and reference accuracy |

A maximum of two defect estimates, a different local subspace or a different boundary-vertex set requires an explicit mapping or fixed-grid comparison; matching function names proves nothing.

**P0 finding:** the archived quasi-interpolator uses arithmetic vertex averaging, while the manuscript specifies area weighting. They coincide on equal-area coarse meshes. A graded-mesh counterexample must remain documented; migration equivalence must not silently replace the archived operator. Uniform-grid SLOD is unaffected by this distinction. Resolve the graded-mesh definition before accepting P2 as manuscript-equivalent.

## 3. Retained and excluded scope

| Component | Decision |
|---|---|
| Triangular meshes, nested refinement, closure, boundary tags, injection | Keep used paths; retain stable identities and parent/child relationships. |
| Conforming P1 Helmholtz, positive energy matrix, mixed boundaries, manufactured data and quadrature | Keep, including corner and localized-Gaussian integration. |
| Quasi-interpolation, kernel constraints, two-sided correctors, PG | Keep direct-Schur and its necessary internal dependencies. |
| Local positive Riesz, multiple RHS, global estimator, regional AS | Core; share patch geometry and factorizations. |
| Regional budgets, Krylov directions, POD dictionary inheritance, global AOT | Core E2 modules, not temporary injected diagnostics. |
| Matrix-free Ritz Theta, E1 lazy and E2 statewise promotion | Keep the used numerical paths; not the whole certification framework. |
| Strong-residual reference marking, global Doerfler and `m_ref` reference refinement | Core scheduler; replace the old reference-only/balance/epoch strategy combination. |
| AFEM/UFEM/SLOD | Minimal drivers for both problems. SLOD means the existing standard fixed-ell LOD baseline. |
| E1 nominal/family marking and all-state 48-load audit | Required; use shared drivers and explicit member IDs. |
| Mesh/state snapshots, checkpoint, independent-load/region/floor/rank diagnostics | Retain as needed without changing production states. |
| Balance gate, `C_rel` reference switching, old reference-only controller | Remove from production, including supposedly disabled branches that still affect behavior. |
| Old balance/candidate-dual observations | Optional diagnostics, off by default; never feed the controller. |
| Equilibrated flux, global-kernel inverse enrichment, explicit singular modes, old fixed-rank training | Do not migrate production implementations; retain only proven internal dependencies if needed. |
| hp/DG platforms, pollution studies, patch GMRES/multigrid/Schwarz variants, E0/E3 | Outside scope; derive dependency closure before pruning by filename. |
| MPFR/MPFI interval certification | Outside scope; keep numerical spectral diagnostics and failure detection without claiming rigorous certification. |
| Date-specific deployment patchers, duplicate giant runners, complete historical archives | Replace by structured run/audit/deployment tools. |

Validated baseline data may be imported instead of rerunning solely to plot, but the new project must be capable of regenerating all three baselines. The 75 E2 fresh audits, two-level corrector remainder estimation and multiple independent singular components are not required production work.

## 4. Mathematical semantics that must remain unchanged

### 4.1 Solutions, residuals and indicators

Use `A_h=K_h-kappa^2 M_h-i*kappa M_R` and positive `E_h=K_h+kappa^2 M_h`. Inner products are linear in the first argument: `a(v,w)=w* A_h v`, `b_kappa(v,w)=w* E_h v`, with `*` the conjugate transpose.

The paper defines `I_H=E_H Pi_H^dg`, with element L2 projection and area-weighted vertex averaging, zero at Dirichlet vertices. Verify `I_H v_H=v_H` and the P0 arithmetic-average discrepancy. Estimator patches use **all geometric coarse vertices** and `W_{h,z}={w in ker I_H: supp(w) subset closure(N^2(z))}`. Do not confuse free coarse unknowns or element corrector patches `N^ell(T)` with this set. Implicit constrained representations are allowed.

- `U_base`: independent rank-zero LOD PG solution. `U`: accepted rank-zero E1 or coupled enriched E2 solution. `u_h`: current reference FEM solution for error diagnostics; `u_c`: optional candidate-mesh diagnostic solution.
- Residual `r_mu=f_mu-A_h U_mu`.
- Solve each constrained positive local Riesz problem and compute `eta_H,mu^2=sum_z ||xi_mu,z||_kappa^2`. This is neither a single global-kernel inverse norm nor the energy norm of the additive vector sum.
- Distribute each nodal squared indicator equally among incident coarse elements, preserving the total.
- Normalize a marking decision using frozen `d_mu=max(||U_mu||_kappa,1e-12)`. Mark the mean normalized squared indicator, then supplement the worst normalized member's Doerfler set. Never use the exact-solution norm for production normalization.
- E2 coarse marking uses `eta_H(U_coupled)`, not the base estimator.
- Reference marking uses the embedded accepted solution's strong residual: body `h_K^2 ||f+kappa^2 U||_K^2`; half of `h_e ||jump(partial_n U)||_e^2` to each interior neighbor; `h_e ||partial_n U||_e^2` on Neumann; `h_e ||partial_n U-i*kappa U||_e^2` on Robin. P1 element Laplacians vanish; the jump sums outward derivatives. Do not double-count edges.
- This strong residual is for marking. Do not claim a reliability theorem for the embedded non-Galerkin solution. Mark globally in both experiments; the AS mask does not restrict reference marking.
- Mean-plus-worst enforces those two bulk conditions, not every member's. Preserve stable tie ordering and redistribution. Nominal-only is the single-member case.

### 4.2 Regional AS, training and compression

`D_R=Omega intersect B_R(0)`, production `R=0.6`. Select a patch only when its entire support is contained in the region. An empty mask produces zero enrichment and an explicit status; do not enlarge the radius or fall back to global AS.

`M_D g=sum_selected I_z B_z^{-1} I_z* g`, with `B_z` the positive energy restriction. Reuse the first seed `M_D r` from estimator local solves only for the identical residual and state. For subsequent accepted normalized directions `phi_j`, form `M_D A_h phi_j`, then deflate the base trial, retained dictionary and current block. Apply energy Gram–Schmidt twice, detect dependence, and normalize. If old coordinates use `A_W`, prove equivalence; do not substitute raw `psi_j` for projected `phi_j` or introduce global `B_W^{-1}`.

Verify `Re(r* M_D r)=eta_D^2`. There is no partition-of-unity averaging. Merged identical patches must retain separate full/selected multiplicities. Raw kernel representatives satisfy regional support and interpolation constraints; projected `phi` may have global LOD components. Do not reject it for nonzero values outside D or claim AS cannot affect localization error.

E2 training:

1. Solve the base training family on the current mesh/ell; freeze `t_mu=max(1e-12,0.1*eta_D(U_base,mu))`.
2. Test the inherited dictionary first. If the maximum regional budget ratio is at most 1, do not grow it.
3. Otherwise choose the worst budget member; add at most 2 directions per block, with working-rank cap 24. Record discarded dependent directions.
4. Exit under the old stagnation rule after two consecutive relative worst-ratio reductions below 1%. Distinguish target attained, stagnation, rank cap, empty region and numerical failure.
5. Set `Dtilde_mu=(U_full,mu-U_base,mu)/d_mu` and `S=Phi_full* E_h [Dtilde_mu]`. SVD of S gives `Phi_r=Phi_full U_r`. For a nonempty dictionary, test ranks from 1 upward by recomputing actual coupled solutions. Every training member must satisfy `||U_full-U_r||_kappa<=0.2*t_mu` and `eta_D(U_r)<=max(1.1*t_mu,1.03*eta_D(U_full))`. If no rank qualifies, retain the full dictionary; an empty dictionary has rank 0. POD discarded energy is not PG solution error; do not restore the old 5% snapshot threshold.
6. Output working rank, independent raw direction count and compressed online rank separately. Production rank 1–2 is an observation, not a hardcoded answer; do not infer dimensions from `16*4`.

Fresh single-RHS retraining starts empty and is diagnostic only. Keep train/test/shift/pure identities distinct. Compression bounds the tested training members' **absolute** reference-relative perturbation by `epsilon_dist*t_mu`; it is not a relative guarantee or a guarantee for unseen loads. Training exit, budget success and compression acceptance are independent fields.

### 4.3 Inheritance and AOT

- Ell-only changes retain same-mesh kernel representatives but refresh LOD projection, global AOT, current targets and budget checks; neither clear them unconditionally nor reuse old projected trial/test pairs.
- H/h changes inject old representatives and apply the production repair `M_D,new E_new P psi_old`, then verify constraints/support. Repair changes directions; it does not exactly preserve the previous space.
- Inherit representatives of the compressed dictionary. Saving representatives and coordinates during construction/POD may avoid a later `C T` recovery solve, but requires equivalence checks. A projected trial vector is not a kernel representative.
- AOT solves the **global reference-domain** equation `A_h* z=E_h phi`, then projects/orthogonalizes against the base test space as in production. A local trial seed does not justify a local test solve.
- Cache the AOT factorization across RHS and ell-only changes; different RHS still need solves. Invalidate for changed reference operator or free-node set. Use one entry or an explicit memory budget.

With base trial/test B,C and `Phi* E_h B=0`, `Phi* E_h Phi=I`, exact raw tests satisfy `Zraw* A_h B=0`, `Zraw* A_h Phi=I`. The raw block system is `[C* A_h B, C* A_h Phi; 0, I]`. Projection/orthogonalization is a test-basis change and must preserve `Y^0`. Record block defects, base-test and actual PG residuals. Approximate AOT does not make the identities exact. Invertibility follows from the base block; the enriched quasi-optimality constant uses its own `alpha_r`, not a presumed `alpha_r>=gamma_h`. Trial inclusion does not guarantee monotone PG error.

### 4.4 Theta and ell promotion

`G_ell(v_H)(w)=conj(a(w,T^*_(ell,h) v_H))` and `Theta_ell=sup eta(G_ell(v_H))/||v_H||_kappa`, using the same local positive Riesz operators. The denominator is the uncorrected coarse function's energy norm. It is not a LOD-function norm, load norm or adjacent-level difference. With real interpolation and complex-symmetric coefficients, primal/adjoint norms agree; validate that equivalence before using their maximum, and never sum them.

Retain the numerical paths in `reference_defect_spectrum_matrix_free` / `build_reference_corrector_certificate`. Record Ritz value, residual, convergence and warm-start identity. Neither a Ritz value nor a small iteration residual alone supplies a rigorous upper bound. Production tolerance is `1e-4`, maximum 750 iterations; unconverged estimates cannot be accepted as normal results.

- E1 checks initially, at prescribed endpoints, or after `2^ell` solved states since the last check. Ratio above 1.2 promotes first; otherwise `Theta<0.2` accepts, and other cases hold. Intermediate reference sweeps count as states. Verify promote/re-solve/counting semantics.
- E2 checks `Theta/eta_H(U)>0.1` after training/compression at every state, raises ell on fixed H/h, rebuilds/retrains/rechecks until accepted or ell=4. Only the final accepted branch enters production curves.
- Reset the lazy counter after every check, not only promotion. Fixed-ell controls may observe Theta but cannot promote.
- Zero denominator with positive Theta triggers promotion up to the cap; zero Theta needs none; both zero records `zero_defect` with NA ratio. Uncomputed values are NA, not cached values pretending to be current.
- Record `ell_cap_reached` when the criterion still fails at 4. These are practical heuristics, not certified accuracy guarantees.

`q_ell=Theta_ell/eta_H(U_nominal)` depends on load amplitude: for fixed space/operator, scaling the load by c scales q by `1/|c|`. Absolute training floors can also change the adaptive trajectory. Keep original thresholds; do not add normalization. This rule does not replace `C_a*delta_ell^2 < alpha_hat_(H,infinity)` or directly control the reference floor.

## 5. Production, baseline and control presets

| Setting | E1 main | E2 main |
|---|---|---|
| Problem | Localized smooth R1 on the unit square | Corner plus Gaussian S on the L domain |
| kappa; initial H/h | 16; 6/10 | 16; 6/10 |
| Boundary/profile | D at y=0/1, N at x=0, Robin at x=1 | D on notch edges; Robin outside; `boundary-weight-gaussian-alpha80`, nominal amplitude 0.5 |
| Train/test/shift | 16/24/8; seeds 20260901/20260902 | 16/24/8, seed 20260831; 2 separate pure diagnostics |
| Coarse theta | 0.15, training family | 0.15, coupled training family |
| Reference sweeps/theta | 2; [0.3,0.2] | 1; [0.2] |
| Reference closure | Local reserve gap 2, active patch layers 0 | No old physical-radius matching or extra E1 reserve |
| Enrichment | None; rank 0 | Regional AS R=0.6, budget, POD and inheritance |
| Ell | 2 to at most 4; lazy 1.2, tau 0.2 | 2 to at most 4; each-state ratio 0.1 |
| Horizon | 25 complete coarse cycles, 51 states | 32 cycles, 33 states |
| Historical extra checks | 35,44,50 | Side branches/pre-promotion solves are not accepted states |
| Quadrature triangle/Gaussian/singular/recursion | 12/16/24/8 | 12/16/24/6 |
| Local/reference solver | Direct Schur / UMFPACK | Same, plus verified AOT factor cache |
| Reference/work DoF caps | 500000/600000 | 300000/400000 |
| Historical threads | 16; BLAS=1 | 16; BLAS=1 |

Extract source/exact/boundary code and numerically verify it, including corner-angle branches, Gaussian normalization, Robin sign and singular integration. E1 nominal center is (0.75,0.5), center disk radius 0.05. E2 domain is `(-1,1)^2 minus ([0,1]x[-1,0])`, with angle `0<theta<3*pi/2`. Its fixed singular shape is `(1-x^2)^2(1-y^2)^2 r^(2/3) sin(2 theta/3)`. Family parameters: singular coefficient [0.75,1.25], center within radius 0.05 of (-0.5,0.5), amplitude [0.25,0.75], phase [-pi/4,pi/4]. Follow the source/`eq:e2-family` envelope and phase shift; shift samples use nominal amplitude/phase. Exact shapes supply loads/audits, never explicit enrichment seeds.

Export complete fixed sample tables `rhs_e1.json`/`rhs_e2.json` with IDs, roles and every parameter. Seeds are provenance, not a replacement for numeric tables. Nominal is already in the 16 training members: there are 48, not 49, formal loads. IDs 0–15/16–39/40–47 are train/test/shift; E1 nominal is 0. E2's 2 pure loads remain separate.

Baseline presets retain their effective configuration: nominal AFEM, theta=0.15, exact-error targets 0.008 (E1) and 0.003 (E2); uniformly refined FEM; SLOD ell=3 with four reference levels beyond H. Full E2 UFEM uses H6–H19 and SLOD H6–H15. Obtain remaining horizons/quadrature from baseline manifests, not plot cutoffs. Exact error may stop these specified baselines but may not drive marking or ALOD decisions. Historical complete baseline states are E1 87/15/11 and E2 88/14/10 for AFEM/UFEM/SLOD.

### 5.1 Replay versus production

`paper-replay` retains E1 35/44 checks and historical audit schedules; E2 may reproduce both promotion diagnostic branches while accepting only inherit. `production` performs identical mathematical updates using only accepted branches, with audits from checkpoints. Expose any extra checks explicitly; removing them creates a new trajectory, not a claimed pointwise replay. Observational balance switches must not change solutions, meshes or ell; if they do, find the hidden state effect.

### 5.2 Four controlled E1 experiments

Generate four full validated configurations from shared settings plus explicit overrides; do not copy four drivers. The controlled choice is marking membership, affecting ALOD coarse/reference marking or AFEM strong-residual marking. Resulting meshes may differ naturally.

| Setting | nominal-ALOD | family-ALOD | nominal-AFEM | family-AFEM |
|---|---|---|---|---|
| Marking IDs | [0] | [0..15] | [0] | [0..15] |
| Initial H/h | 6/10 | 6/10 | 6/no reference | 6/no reference |
| Ell/rank | Fixed 3/0 | Fixed 3/0 | N/A | N/A |
| theta_H | 0.15 | 0.15 | 0.15 | 0.15 |
| m_ref/theta_h | 2/[0.3,0.2] | Same | N/A | N/A |
| Reference closure | Local reserve 2, active layers 0 | Same | N/A | N/A |
| Updates/states | 25/51 | 25/51 | 96/97 | 96/97 |
| Audit per state | Same 48 RHS | Same | Same | Same |

Common conditions: kappa=16, quadrature 12/16/24/8, 16 threads, BLAS=1. ALOD uses direct-Schur/UMFPACK; AFEM uses batched UMFPACK. Disable balance control/default observers and candidate FEM solves. Fixed-ell ALOD observes Theta at initial/final/every 8 states without promotion; do not import main-production extra checks 35/44. Audit reference/exact errors never alter marks or fixed horizons.

Control resource starting points: ALOD reference/work caps 600000/700000; AFEM 500000 unknowns and 1200000 elements; 14-day wall limit; historical server RSS caps 125 GiB (ALOD) and 24 GiB (AFEM). Reassess resources on the current host; do not silently substitute the smaller E1-main limits.

Acceptance: identical samples; identical initial ALOD solution/error/reference accuracy; exactly 48 audit rows per state; 1/16 marking members; mean/worst bulk; at least one different marking set; `48*(51+51+97+97)=14208` rows. In nominal arms the other training loads are audited only. Validate/import existing data for plotting; short nominal/family numerical regressions and full-horizon capability are still required. Do not rerun four expensive campaigns merely to generate plots.

### 5.3 Separate identities and purposes

`experiment_id` separates E1 variable-ell production, fixed-ell controls and E2 AS production. `method` and `marking.member_ids` are independent. `stopping` explicitly selects fixed cycles or baseline-only exact targets; `audit` selects states/loads; `plot` never controls solving. Main E1 terminal DoF 2422 must not replace control DoF 2522/2878. Compare family benefit at common budgets, not unequal final endpoints.

## 6. Unified reference-refinement schedule

Use positive integer `reference.sweeps_per_coarse_update=m_ref` and a theta array of exactly that length; a scalar may expand explicitly. Export the effective array. A sweep is neither a mesh level nor a quadrature recursion.

```text
initialize H,h -> solve_accept(initial) -> save
for each coarse cycle:
    form coarse and first-reference marks from current accepted family solutions
    refine both meshes; apply compatibility/nesting/gap closure
    promote work mesh to reference
    solve_accept(cycle, sweep=1) -> save
    for sweep=2..m_ref:
        retain H; recompute strong residual from the preceding NEW accepted solution
        mark/refine reference, close and promote
        solve_accept(cycle, sweep) -> save
```

`solve_accept` completes space updates, training, compression and ell checks. E1 can promote ell at its second sweep without refining H. E2 must not commit pre-promotion solutions. Separate `coarse_cycle`, `reference_sweep`, `state_id`, `ell_event_id`, `mesh_revision`. Initial state is 0; K complete cycles give `1+K*m_ref` accepted states. Ell iterations/diagnostics do not add accepted states. A stop may leave a partial cycle; save its next stage.

Test m_ref=1/2/3, including interrupted cycles. Functional m_ref=3 support is not evidence of better science. `[0.3,0.2]` is neither one theta=0.5 mark nor two markings of an old solution.

## 7. Architecture and extraction boundary

Aim for a few libraries, two main C++ entry points and shared tools:

```text
include/alod/       typed public interfaces
src/mesh/          geometry, nesting, closure, interpolation
src/fem/           P1, boundaries, energy, quadrature
src/lod/           correctors, PG, localization
src/as/            local Riesz, region, training, POD, AOT
src/adapt/         family marking, m_ref and ell state machine
src/problems/      E1/E2 and their pure components only
src/io/            states, checkpoints, provenance, mesh export
src/diagnostics/   optional observers with no control feedback
apps/alod_run.cpp  methods and explicit configurations
apps/alod_audit.cpp independent saved-state diagnostics
configs/           production, baselines, controls, smoke and audit templates
data/rhs/          fixed numeric sample tables
tests/fixtures/    small old numerical/matrix/mesh references
scripts/run_server.py
tools/analyze.py
docs/
```

Extraction starts at mesh refine/edges, quasi_interp and adaptive hierarchy; Helmholtz operators/boundary/quadrature/manufactured/paper_cases; model/corrector/patch_system/patch_solver/corrector_pipeline/patches; production kernel_residual; certificates' used Ritz paths; optimized E2 runner/helper/includes; reference_epoch driver marking/scheduler; adaptive-paper baseline driver; and the E1/E2 paper-analysis scripts. Trace transitive includes before pruning. Do not retain all hp/Schwarz alternatives, the whole epoch controller, giant duplicated runners, dated data roots, or obsolete E2 fresh overlays.

`Problem` supplies loads and an optional exact evaluator. Adaptive controllers do not receive exact errors, `u_h` or `u_c`. Mesh/operator contexts are shared read-only; caches belong to explicit mesh/operator/space versions. Diagnostics consume immutable accepted states. Baseline exact stopping remains outside the marker. The E2 base LOD solve is part of training, not removable audit overhead.

| Change | Reusable | Must update/recheck |
|---|---|---|
| RHS only | Mesh/operator/corrector/Riesz/AOT factorizations | Load, solution, residual, targets; exact/floor cache is RHS-specific |
| Ell only | Reference operator/FEM solves/Riesz geometry/AOT factor | Trial/test, projections, Theta, base solve, budgets, coupled system |
| Reference h | Identity-verified unchanged local geometry | Operator/free nodes, Riesz, AOT, correctors, injection/repair, mesh indicators |
| Coarse H | Only dependencies proved unaffected | Interpolation/kernel/patches/correctors/repair/projection/control values |
| kappa, coefficients or boundaries | Valid pure geometry | All numerical factors; equal dimensions/mesh IDs are insufficient |

Start with exact sparse structure/coefficient/free-node comparisons. A future revision-ID cache may avoid O(nnz) comparisons only after every modifying entry point is proven to advance identity. Bound caches and report their statistics.

Dependencies: C++20, CMake, Eigen, OpenMP, and the SuiteSparse components actually needed. Do not silently equate an Eigen fallback with production UMFPACK. Remove unused TBB/MPFR/MPFI. Python analysis uses an explicit small dependency set such as numpy/matplotlib; no web/UI/service framework is required.

## 8. Results, diagnostics, checkpoints and imports

### 8.1 Always-recorded data

`manifest.json`: effective config, source/new commit, dirty-patch and binary hashes, backend/dependency versions, threads/CPU set, sample hash, run status and reason.

`states.csv`: all counters, mesh/dictionary hashes, ell, N_H/N_reference/N_work/N_online, marks/closure increments, working/online rank, training exit, worst budget ratio, compression status, PG/AS/AOT residual summaries. `N_online=N_H+r_online`. Count free unknowns after boundary elimination; store node/element counts separately.

`samples.csv`: state/load role, nominal/marking flags, parameter and solution identities, exact norm and errors. `marking.csv`: worst member, frozen scales, mean/worst totals and captured masses. `plots.json`: selected states, cutoffs, fit windows, interpolation endpoints/weights and source hashes. E1 controls cannot be `paper_complete` until all 48 loads per state have been validated, even if audits are later batches.

Old `E=e/||u_exact||`, `F=f/||u_exact||`, `G=g/||u_exact||`; manuscript F is now `E_ref`. Preserve source columns and map them explicitly. AFEM reference/gap are NA; legacy N_reference=0 means inapplicable, not zero floor.

`ell_checks.csv`: Theta, denominator, ratio, state identity, reason, trigger/cap and spectral convergence. Unchecked states use NA; displayed cached observations must name their computation state and mesh/ell.

`timing.csv`: unique event and parent IDs and phase. Separate mesh/closure, assembly, correctors, Ritz, base solve, AS prepare/apply, training evaluation, POD, AOT factor/solve/project, coupling, load, integration, reference solve and I/O. Nested times cannot simply be added; shared construction cannot be charged twice.

### 8.2 Diagnostic menu

| Diagnostic | Requirements/default |
|---|---|
| exact/floor/gap | Absolute and exact-normalized e/f/g for the same RHS; production decisions do not use them; paper nominal/all-family observations may be offline. |
| E1 controls | Required all-state held-out/reference errors, members/marks, bulk and matched-budget paired statistics. Terminal-only old AFEM audits cannot manufacture curves. |
| reference/candidate | Optional solution gaps, candidate mesh identity and DoF; not required for strong marking. |
| old balance | Off; preserve old local candidate-Riesz meaning if requested, without C_rel or control feedback. |
| localization | Distinguish coupled/base denominators, adjacent-level solution changes and actual error; required production Theta frequency remains. |
| regional error | Inside/outside base/coupled gap, pure components, support/kernel/AS checks. Handle circle-cut cells by integration, not barycenter classification. |
| training/POD/rank | Worst ratios/member, growth, independent directions, ranks, compression distances/residuals; fixed-rank and disabled compression are separate controls. |
| inheritance | Incoming/repaired/dropped ranks, support/kernel defects and time; separate ell-only and mesh changes. |
| shared/fresh | Optional identical-space single-load retraining; record signed changes and budgets/ranks. The 75 historical audits are not a plotting prerequisite. |
| rank-zero versus AS | Same H/h/ell/RHS, preferably reusing base solves; unequal adaptive meshes do not isolate enrichment effects. |
| scaling/two-level corrector | Optional q amplitude test and specified L>ell difference/remainder; current paper does not claim evaluated two-level certification. |
| reference deepening | Independent saved-state refinement with actual re-solves, bounded rounds and resource/timing records; no feedback into adaptation. |
| ell/radius sensitivity | Small explicit controls, not the full historical script platform. |
| memory/performance | Process-group RSS, MemAvailable, thread affinity/participation, factor/cache counters and scratch estimates; detailed profiling optional. |

Global AS is an economical all-patches mask diagnostic, not permission to restore the global-kernel inverse. Adjacent-level diagnostics distinguish frozen-representative reprojection from retraining to budget. Reuse existing per-member/base quantities.

Diagnostics must not alter dictionary, counters, Ritz warm starts or marking caches. Prefer independent processes from immutable checkpoints; in-process audits must prove unchanged state and future trajectory. Near-zero denominators give NA with a reason. Floor/e and gap/e do not add; `eta_D^2/eta_H^2` is an overlapping patch-indicator ratio, not a physical-region error fraction. Adjacent-level differences are sensitivity, not an isolated localization component or monotonic-error guarantee.

### 8.3 Genuine resume

Avoid old prefix replays such as E1 45-to-51 and audits at H0/H16/H32. Save:

- Mesh nodes/elements/boundaries, stable IDs, parent/child/refinement levels and closure/injection state; a final VTU alone is insufficient.
- Accepted counters, next stage/sweep, ell, last lazy check, extra check schedule and Ritz warm starts.
- Final raw kernel representatives, compression coordinates/rank, accepted solution or deterministically reconstructible data, fixed samples and needed RNG state.
- Config/operator/sample/mesh/dictionary/format hashes, stop reason and committed log position.

LU factors may be rebuilt. First prioritize portable correctness. Optional space snapshots must state storage costs and version checks. Write temporary files, validate, atomically commit the checkpoint, then commit accepted-state logs. Preserve the preceding valid checkpoint after interrupted writes. Resume without duplicate rows or recomputing committed sweeps. Only horizon/resource increases may be controlled overrides; mathematical changes create a new branch identity.

Compare uninterrupted and resumed runs after the first reference sweep, ell promotion and training boundary. Require matching final meshes/ell/rank/errors and direct H16 diagnostics without H0–15 replay.

### 8.4 Completeness, large fields and failure isolation

Old marking `element_ids` exceeded Python's 131072-character CSV default. This postprocessing failure was not OOM. New formats store versioned/checksummed mark sets separately; CSV stores set ID/count/hash. Legacy readers accept at least 64 MiB per field and explicitly reject overflow/truncation.

Separate `solver_completed`, `validation_passed`, `audit_complete`, `paper_complete`. Overall success requires horizon, records and checks. Validation failure must not kill unrelated authorized jobs. Distinguish shared-resource watchdog failures. Resume missing audits from committed checkpoints without replaying adaptive prefixes.

For the four controls, scientific input is `data/<arm>/family_control_samples.csv` plus completion/hash evidence. Family-ALOD states 0–48 are original and 49–50 are verified recovery; family-AFEM 0–95 original and 96 recovery. Recovery `runs/*` prefixes can skip audit rows; never overwrite merged canonical tables with them. Preserve prefix error/mark verification. The ordinary AFEM complete 97-state audit matches the old first 87 DoF states with nominal error relative difference about `5.06e-13`.

Portable result packages contain effective configs, sample tables, completion markers, canonical data, plot inputs and script dependencies. Unpack into another directory and validate/plot there, including entry markers such as `all-controls-completed.json`. Use a separate LaTeX artifact map, not dated exact-string manuscript replacement scripts.

## 9. Performance evidence and priorities

### 9.1 Already implemented in the legacy sources

- Local AS factor reuse, multiple RHS and first-seed reuse with correct weights/residuals.
- OpenMP load/error/regional integration and thread-local quadrature caches; fixed sample output order and safe initialization.
- Same-state reference/exact-norm/floor caches and adjacent-level construction reuse, keyed by operator/RHS.
- AOT factor reuse: 26 calls, 1 factorization, 25 hits; difference about `2.45e-14`, ell/mesh comparisons about `4.21e-14`. Preserve exact free-node/structure/value invalidation and a bounded cache.
- Whole-process affinity before OpenMP initialization and BLAS=1; verify actual workers, not just environment variables.
- E1 four-arm batched FEM/audit solves and parallel integration from archived production sources, with stable sample order.

The AOT cache evidence is small-scale correctness/hit validation: factorization about 0.005 s and hit checking 0.00034 s, not a measured large-production speedup.

### 9.2 Remaining optimization work

Historical E2 inherit (33 states): 58 min 08.5 s wall, family training about 105.14 s (~3%), error audit about 2076.22 s (~59%), load assembly about 346.94 s (~10%), peak RSS 31.14 GiB, average CPU 1076% (~10.76 CPU equivalents). Phase conventions overlap and include side diagnostics/old AOT behavior; do not sum percentages as a complete cost decomposition.

Historical E1 nominal-ALOD all-load audit: 31 min 50.33 s, 1227% CPU, 81.52 GiB peak; nominal-AFEM: 23 min 38.68 s, 1217%, 1.425 GiB. Different work/meshes prevent a direct speedup claim. Recovery family jobs skipped audited prefixes and are not comparable whole-run timings. Profile before attributing the 81.52 GiB peak to a particular cache.

Optimization priority P1 (distinct from migration phase P1): separate production/audit; direct checkpoint access; bounded E1 audit batches with complete 14208-row semantics; stream/externally store long marks; reuse RHS/geometry/integration/reference factors; execute only accepted E2 production branches. Preserve trajectories and required paper data; report production versus audit wall, peak memory, storage and restore cost.

Priority P2: block RHS AOT/PG solves; cache basis products and small matrices; evaluate raw-AOT block solves and POD linear-response reuse only after block/PG equivalence; save raw representatives/coordinates with inherited-budget checks; compute regional residuals during training but the final global eta for coarse marking; bound patch/scratch/Nref-by-RHS storage and preserve stable aggregation. Do not assume shared-factor solve thread safety.

Priority P3: profile Theta apply/matvec/warm starts without relaxing frequency/tolerance; evaluate unified forward/adjoint LU and symbolic reuse only after verifying conjugate-transpose support. All are candidates until measured.

Global AOT still requires at least one large factorization. Local AOT, iterative approximations, different R or altered budgets are mathematical variants, not migration speedups. Compare identical H/h/ell/dictionary/backend/diagnostic intensity/threads in an independent resource window. Report cold/warm behavior and median/range of at least 3 short repeats. Test 16 versus 32 threads explicitly; more cores does not prove speedup.

## 10. Phases and acceptance gates

Proposed commands/file names below are interfaces to implement, not evidence that an executable exists.

| Phase | Dependencies | Deliverables and gate |
|---|---|---|
| P0 Sources | None | Manifest, paper labels/includes, effective production/control configs/samples, snapshot differences, small old oracles and dependency closure. Identify algorithm sources and hidden branches; separate production/controls; reproduce old small states. |
| P1 Minimal core | P0 | Independent CMake, mesh/FEM/boundaries/quadrature, both problems, three minimal baseline drivers and mathematical regression. Clean full build without old libraries; problem/boundary/mesh/baseline smoke passes. |
| P2 LOD/estimators | P1 | Two-sided LOD/PG, batched Riesz, explicit marking members, strong residual, Theta. Fixed-grid old/new agreement, kernel/vertex/mass checks; distinguish AFEM/ALOD estimators. |
| P3 E1 scheduling | P2 | m_ref state machine, lazy/fixed ell, main/control presets and events. Complete/interrupted m_ref=1/2/3 cycles and short promotion/control trajectories. |
| P4 E2 method | P2,P3 | Regional AS, budgets, POD, inheritance, AOT cache and ell rules. Empty/full/partial masks, AS identity, forced 2→3→4, inherit/reset and fixed-state regression. |
| P5 State/audit | P3,P4 | Resume, independent audit, all-state loads, long-field import and optional fresh/rank-zero. Three interruption positions, audit isolation, direct H16, no unrelated-job cancellation. |
| P6 Performance | P4,P5 | Carry existing optimizations, profile/implement candidates. Same-grid equality and correct resources; report measured gains or no gain. |
| P7 Integration | P6 | Small/medium E1/E2, baselines, required diagnostics and deployment smoke. No old dependencies/absolute paths; explicit missing-data/resource/failure parsing. |
| P8 Paper | P7 | Main 51/33 and controls 51/51/97/97 reproduced or differences stated; baselines validated/imported/rerun; current 7 figures/8 tables. Complete 14208 rows and traceable budgets/cutoffs; no historical data presented as new output. |
| P9 Release | P8 | README, algorithm/config/audit/resources/provenance/performance docs and portable package. Independent-directory build/smoke, scoped features and explicit remaining limitations. |

Update status with commands, inputs/hashes, maximum numerical difference, resources and next action after each phase. Do not concurrently edit a shared giant runner.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DALOD_BUILD_TESTS=ON
cmake --build build -j <build_jobs>
ctest --test-dir build --output-on-failure -L smoke
build/alod_run --config configs/e1-paper.json --output results/e1
build/alod_run --config configs/e2-paper.json --output results/e2
build/alod_run --config configs/e1-marking-family-alod.json --output results/e1-family-alod
build/alod_run --config configs/e1-marking-nominal-afem.json --output results/e1-nominal-afem
build/alod_run --resume results/e2/checkpoints/latest --max-coarse-updates 40
build/alod_audit --checkpoint <saved-state> --config configs/audit-frozen-family.json
# Optional historical audit:
build/alod_audit --checkpoint <H16-checkpoint> --config configs/audit-e2-fresh.json
python3 tools/analyze.py --run results/e2 --paper-output artifacts/e2
python3 tools/analyze.py --campaign results/e1-marking-controls --paper-output artifacts/e1-controls
```

Reject unknown fields, wrong theta lengths, train/test misuse and mismatched checkpoint/operators. Missing SuiteSparse must not silently change the production backend. Use explicit resource-limit semantics, not ambiguous zero-as-unlimited values.

## 11. Numerical and paper acceptance

### 11.1 Required meaningful tests

1. Both manufactured sources, boundary residuals and quadrature convergence; FEM residual/errors; nested injection of P1 functions.
2. Interpolation projection and geometric/free vertices; local Riesz energy/kernel; mask AS identities, duplicate-patch multiplicity and global/regional eta.
3. Normalized mean mass, worst-member supplement, stable ties and old strong-residual boundary/jump equivalence.
4. Global AOT equations/block/test/PG residuals and cache behavior under RHS, ell, equal-size coefficient changes, sparsity changes and free-boundary changes; correction-coordinate POD.
5. Training success/stagnation/cap, compression fallback, empty region, repaired kernel/support, accepted versus pre-promotion states.
6. m_ref/lazy counters, forced/no/capped promotion, incomplete cycles, resume/audit isolation.
7. One/multiple-thread agreement and observable actual worker affinity/participation.
8. Four control memberships, fixed ell, 48 RHS per state, 87-point exact-stop versus 97-point fixed horizon; per-sample interpolation before quantiles, no extrapolation and actual reference endpoints; no interpolated main target crossing.
9. Legacy marking fields above 131072 characters; explicit missing/duplicate/partial-checkpoint/parse failure; solver completion distinct from complete audit; no unrelated cancellation; independently unpacked portable package.

Start fixed-grid numerical comparisons at `abs<=1e-12+1e-8*scale`, using meaningful scales. PG/AOT checks follow source standards (old audits used `<1e-7`). Near-zero quantities must not create fake relative failures. The measured AOT `1e-14` differences do not require bitwise equality everywhere.

Adaptive equivalence is statewise: compare mesh/hash, ell, rank, worst member and marks before checking terminal errors. Stop at the first divergence. Distinguish floating tie/reduction effects from mathematical differences; if necessary compare a common frozen state and label the new trajectory. Never tune thresholds to force agreement.

### 11.2 Historical paper oracles

These identify wrong versions; never hardcode them as answers.

| Quantity | E1 | E2 inherit |
|---|---:|---:|
| Accepted states | 51 | 33 |
| Final online DoF | 2422 | 8911=8910+1 |
| Final reference DoF | 316767 | 108169 |
| Final ell / promotion | 3; state 12, 2→3 | 3; H_step 8, 2→3 |
| Nominal relative exact energy error | 0.009205652787925123 | 0.008294536369324071 |
| Relative reference floor | about 0.0039472696 | 0.004130545463867472 |
| Relative reference gap | about 0.0083161506 | 0.007192366436660725 |
| Tail DoF order | 0.644557; 5 complete cycles | 0.625629; last 10 states |

Frozen final E2 test median/p90/max are approximately 0.008450/0.009000/0.009205; max/nominal about 1.110. This supports low-rank reuse for one fixed singular shape plus nearby packets, not multiple independent singular shapes or isolated AS effects on unequal meshes.

Optional historical 75 fresh audits: 24-test median `100*(E_shared/E_fresh-1)` at H0/H16/H32 is about -0.0217%/3.0310%/-1.4051%. These removed artifacts are not current-paper acceptance and do not assert shared superiority for every load.

E1 tail states are 42/44/46/48/50, with 4/6-cycle sensitivity. Repeated-online-DoF reference sweeps remain on curves, not as independent coarse-fit points. E2 reports 5/8/10-state windows and R². Fit `log(error)=a-p log(N_online)`; p is a DoF order, not directly an h order.

| Control | Full states | Final N_on | Final N_ref | Final held-out max E | Figure 3 points / final shown state |
|---|---:|---:|---:|---:|---|
| nominal-ALOD | 51 | 2522 | 266247 | 0.01184514730 | 26 / 50 |
| family-ALOD | 51 | 2878 | 320194 | 0.01016350548 | 25 / 48 |
| nominal-AFEM | 97 | 192960 | N/A | 0.00606796517 | 81 / 80 |
| family-AFEM | 97 | 272164 | N/A | 0.00489073616 | 77 / 76 |

Displayed ranges end at each retained curve's first test maximum <=0.012; this is not a stopping rule. At budgets 500/1000/2000, family-minus-nominal max-error improvement is ALOD +7.51%/-2.85%/+2.50%, AFEM +1.03%/-0.24%/+2.94%; negative means worse. At 2000, 20/24 ALOD test loads improve; reference-error p90 is about 0.620% for both. Do not assert that family marking must improve every comparison.

At 2000, table errors are percentages: ALOD nominal/family max 1.352/1.318, AFEM 5.981/5.805. ALOD reference-DoF interpolation endpoint ranges are 149.1–198.3 and 132.9–179.3 thousand. First retained `E<=0.015` occurs at 1816/29700 online DoF for family-ALOD/AFEM, actual maxima 0.01457/0.01497; ALOD reference DoF is 132888. These are online-compression observations, not offline-time or memory advantages.

### 11.3 Current figure/table deliverables

Follow labels if paper numbering changes.

| Artifact | Data/selection requirements |
|---|---|
| Fig.1 `fig:e1-error-dof`, `e1_error_dof.pdf` | Variable-ell E1 four-method nominal curves and ALOD exact/reference/gap; all 51 states; 5 complete-cycle tail fit, not fixed-ell controls. |
| Fig.2 `fig:e1-final-meshes`, `e1_final_meshes.pdf` | Main ALOD H/h and original 87-state AFEM final meshes; distinguish nodes/free DoF. |
| Fig.3 `fig:e1-family-marking`, `e1_family_marking_control.pdf` | Left nominal/family ALOD, right AFEM; 24-test median/p90/max, common 1.2% display rule; retain full data. |
| Fig.4 `fig:e2-error-dof`, `e2_error_dof.pdf` | Four nominal methods and ALOD three error components; no default removed fresh crosses. |
| Fig.5 `fig:e2-final-meshes`, `e2_final_meshes.pdf` | ALOD H/h and full AFEM final mesh; shade/outline D=Omega intersect B_0.6 in the first two. R is not a LOD patch physical radius. |
| Appendix Figs.6/7 `fig:e1-rhs-centres`, `fig:e2-rhs-centres` | E1 test/shift and final errors; E2 projections of all 48 sample centers/errors. Center coordinates do not display E2 amplitude/phase variation. |
| Tab.1 `tab:adaptive-config` | Effective main-production settings; identify fixed ell=3 controls. |
| Tabs.2/3 `tab:e1-crossings`, `tab:e1-terminal-results` | Discrete first crossings, full terminal states/tails; ordinary AFEM has 87 states. |
| Tab.4 `tab:e1-family-control`, `e1_family_control.tex` | Four arms at budget 2000, nominal/test percentage errors, true reference endpoint ranges and E_ref p90, including ordinary AFEM. |
| Tab.5 `tab:enrichment-config` | Region, target, block/working rank, stagnation and both compression tests. |
| Tabs.6/7 `tab:e2-terminal-results`, `tab:e2-rhs-methods` | Full endpoints/fits and frozen final-space test quantiles/max-to-nominal; no efficiency ranking. |
| Appendix Tab.8 `tab:e1-rhs-methods` | Original four final E1 spaces' errors and growth relative to their own nominal errors; not a controlled marking-benefit experiment. |

Statistics:

- Main equal-accuracy comparison uses first recorded E<=target without interpolation, including intermediate E1 sweeps. Main AFEM is nominal-adaptive and ALOD family-adaptive; unequal endpoints are not strict equal-accuracy timing comparisons.
- Fig.3 ALOD keeps initial/completed-cycle states 0,2,...,50; AFEM keeps every state. Cut each at first max-test error <=0.012 and save full tails separately. Do not fabricate AFEM curves from terminal-only old audits.
- Budgets 500/1000/2000: log-N/log-E interpolation for each matched RHS between adjacent retained states, then linear type-7 quantiles of 24 tests. Never interpolate quantile curves first or extrapolate. Treat reference error similarly and report actual N_ref endpoint ranges. Count improvements by matching sample IDs.
- E2 main display ends each curve at its first point <= the actual ALOD terminal error (historically 0.008294536). Full baseline endpoint tables/AFEM mesh use their full horizons. Document different selections; record a changed new-run cutoff instead of forcing the old number.
- Do not restore `e1_rhs_growth.pdf` as Fig.3, or unreferenced `e2_rhs_reuse.pdf`, `e2_rhs_samples.tex` or fresh crosses to default paper outputs. Keep optional diagnostics separate.
- Optional plots include Theta, ell sensitivity, regional gaps, training/rank, same-state rank-zero/AS, reference deepening and old balance. The 75 fresh solves are not prerequisites.
- Raw CSV is read-only. Derived data records source hashes, exact states/loads, windows, interpolation weights/reference endpoints and missing values. NA stays NA; side states do not enter accepted curves. Use current E_ref headers with explicit old F mapping.
- `tools/analyze.py` emits data, figures, TeX tables and an artifact manifest, without editing the manuscript by default. Installing into the paper is a separate source-hash-protected copy/build step, not old exact-string replacement.

## 12. Server and resource requirements

Run small implementation tests in a local WSL Ubuntu 22.04 workspace. Full production may use the existing SSH host in an isolated `<REMOTE_ALOD_ROOT>` directory and distinct tmux sessions; do not reuse old output folders/sessions.

Discover CPU topology and allowed CPUs. Prefer disjoint physical-core allocations, initially 16 physical cores per job or an explicit smaller count; BLAS=1. Set whole-process affinity before OpenMP workers exist. Begin with the validated `OMP_PROC_BIND=false`, removing incompatible OMP_PLACES, and inspect actual worker sets. Increasing thread count in a single-core-bound shell is insufficient.

- Read topology, MemAvailable, existing jobs and comparable measured peaks before launch. Multiply peak estimates by at least 1.5, account for remaining growth of active jobs and retain host headroom.
- On the historical large server, initial headroom is max(80 GiB,20% physical memory), not an OOM guarantee. Small WSL uses smaller jobs and explicit appropriate limits.
- Monitor the entire process group, MemAvailable and swap. Bound reference/work DoF, patch work, dense batch sizes and wall time separately.
- E1 nominal-ALOD all-state audit measured 81.52 GiB; do not use E2's old 31.14 GiB estimate for it. Production and audit need different models; AFEM has its own factor-fill behavior.
- Estimate/check before expensive allocation/factorization and supplement with a group watchdog. DoF alone does not predict LU fill; reserve factors cannot guarantee every input is safe.
- Do not launch the next batch without resources. Support safe checkpoint-and-exit. Failures must not cause unlimited restart/limit escalation; retain a clear reason.
- `completed` requires horizon, valid complete data and successful termination. Resource/user stops are `paused_resource`/`paused_user`, not convergence.
- Postprocessing failure records `validation_failed` while preserving committed solver output. Pause dependent work only; unrelated jobs require a resource or user reason to cancel.

Begin with one medium resource/performance measurement. Parallel E1/E2 or immutable-checkpoint audit batches follow measured headroom. Avoid one large AOT factor per process/RHS; prefer shared factors and bounded RHS batches. Additional processes require measured justification.

## 13. Release checklist

- [ ] Clean independent-path build, smoke, import and plotting without legacy paths/libraries.
- [ ] Explicit algorithm selection rather than hashes/directories; all three baseline drivers and validated imports.
- [ ] Four E1 arms share marking interfaces and can produce/validate 51/51/97/97 states and 14208 rows.
- [ ] Separate variable-ell main E1, fixed-ell controls, and 87/97-state AFEM conventions.
- [ ] Distinguish eta_H, eta_D, Theta and adjacent-level sensitivity.
- [ ] Correct Theta denominator and documented amplitude sensitivity/Ritz/heuristic limitations.
- [ ] Fresh-solution reference sweeps and correct m_ref=1/2/3 behavior.
- [ ] No balance production control; reference/candidate FEM is diagnostic except required base training.
- [ ] Global AOT and correct factor reuse; complete regional kernel implementation.
- [ ] Kernel-representative inheritance with separate ell-only and mesh-repair paths.
- [ ] Direct resume/audit without prefix replay and without production-state mutation.
- [ ] Mathematical equality, wall, peak memory and diagnostic intensity reported separately.
- [ ] Traceable paper data, selections, cutoffs and unequal endpoint conventions.
- [ ] Current 7 figures/8 tables; four-arm absolute held-out errors; removed E2 artifacts optional.
- [ ] Large-field import, failure isolation, canonical recovery data and independent package validation.
- [ ] Explicit remaining unknowns: large AOT benefit, 32-thread benefit, m_ref>2 science and nonrigorous regional/localization separation.

Final delivery distinguishes implemented capabilities, reproduced results and candidate optimizations. Planned gains are not measured conclusions.
