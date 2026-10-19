# Implementation correspondence and operational terminology

Manuscript notation and experiment scope were reviewed on 2026-10-19 against both `../LOD_paper/ALOD/version7/arxiv.tex` and `main.tex`; their mathematical bodies agree. The source map and paper/code differences below were checked against the current working tree. Existing storage observations and frozen production records retain their dated provenance; no server experiment was rerun or reaudited for this documentation update. Repository paths are relative to its root.

Contents: source map; coordinate representations; version 7 symbol mapping; paper experiments and historical identifiers; code/paper differences; audit metrics; serialization and retention.

## Source map

| Topic | Primary source to inspect |
| --- | --- |
| PDE, manufactured E1/E2 functions, boundary policy | `src/problems/legacy_cases.cpp`, `src/problems/problems.cpp`, `src/fem/operators.cpp`, `src/fem/boundary.cpp` |
| Nested meshes, NVB, embeddings and interpolation | `src/mesh/nested.cpp`, `src/mesh/refine.cpp`, `src/adaptive/refine.cpp`, `src/lod/quasi_interp.cpp` |
| Local correctors, trial/test bases, reduced solve | `src/helmholtz/patch_system.cpp`, `src/helmholtz/patch_solver.cpp`, `src/helmholtz/corrector.cpp`, `src/lod/lod.cpp` |
| Riesz solves, coarse/fine family marking | `src/estimator/riesz.cpp`, `src/estimator/marking.cpp`, `src/fem/afem.cpp` |
| Defect operator and largest eigenvalue | `src/estimator/kernel_defect.cpp`, `src/localization/theta.cpp`, `src/localization/dense.cpp` |
| Regional enrichment, lifting, inheritance and POD | `src/enrichment/regional.cpp`, `include/alod/regional.hpp` |
| Archived adjoint-test compatibility | `src/enrichment/aot.cpp`; inspect selected `enrichment_tests` before using its identities |
| Adaptive state, ell policy and target stopping | `apps/adaptive_run.cpp`, `src/adaptive/schedule.cpp`, `include/alod/adaptive.hpp`, `tools/adaptive_config.py` |
| RHS identities and selection | `apps/fixed_support.hpp`, `data/rhs/rhs_e1.json`, `data/rhs/rhs_e2.json`, `tools/run_fixed.py`, run configuration/member table |
| Disk checkpoint/audit state | `include/alod/checkpoint.hpp`, `src/io/checkpoint.cpp`, `tools/checkpoint_io.py` |
| Audits, integration and scientific diagnostics | `apps/audit_run.cpp`, `apps/audit_diagnostics.cpp`, `src/fem/batch.cpp`, `tools/run_audit.py`, `tools/scientific_diagnostics.py` |
| Asynchronous publication, workers, shared memory | `tools/async_audit.py`, `tools/native_audit.py`, `apps/shared_audit.cpp`, `apps/audit_worker.cpp` |
| Presets, production and delivery | `configs/adaptive/`, `tools/make_production_campaign.py`, `docs/experiments.md`, `docs/scientific_workflow.md`, `tools/paper/`, `tools/paper_snapshot.py` |

Regression tests such as `tests/test_mathematics.cpp`, `tests/test_kernel_lift.cpp` and `tests/test_relative_theta.cpp` provide bounded numerical evidence. They do not establish that historical full experiments used today's executable.

## Coordinate representations

Introduce bold symbols for nodal coordinate representations. For the paper's unit coefficients,
\[
\mathbf A_h=\mathbf K_h-k^2\mathbf M_h-\mathrm i k\mathbf M_{\Gamma_R},\qquad
\mathbf E_h=\mathbf K_h+k^2\mathbf M_h.
\]
These are the Helmholtz system matrix and energy matrix; \(\mathbf K_h\) alone is the stiffness matrix. The assembly API also supports elementwise coefficients and boundary beta; do not silently identify these extensions with unit-coefficient paper experiments.

Let \(B,C\) be base trial/test basis matrices. Then
\[
\mathbf A_{\rm red}=C^\dagger\mathbf A_hB,\quad
\mathbf b_{{\rm red},\mu}=C^\dagger\mathbf b_\mu,\quad
\mathbf U_{H,\ell,\mu}=B\mathbf c_\mu.
\]
\(\dagger\) denotes matrix conjugate transpose; it is distinct from the paper's star on adjoint correctors. Current complex-symmetric assembly gives \(C=\overline B\), without making \(\mathbf A_h\) Hermitian. `LodSpace::trial()`, `test()`, `reduced()` represent \(B,C,\mathbf A_{\rm red}\).

The introduced basis matrix \(B\) is distinct from version 7's coarse energy matrix \(B_H\). For a basis \(\{\phi_i\}\) of uncorrected \(V_H\), the paper uses \((B_H)_{ij}=b(\phi_j,\phi_i)\) and \((R_\ell)_{ij}=\sum_z b(\xi_z(G_\ell(\phi_j)),\xi_z(G_\ell(\phi_i)))\), with \(\Theta_\ell^2=\lambda_{\max}(R_\ell,B_H)\). These coarse basis functions \(\phi_i\) are not the enrichment basis columns of the code's `phi` matrix \(\Phi\).

For `phi` (enrichment basis matrix) \(\Phi\), lifted tests have columns \(Z=\Phi-B(I_H\Phi)_{\rm free}\). With \(X=[B,\Phi]\), \(Y=[C,Z]\), solve
\[
(Y^\dagger\mathbf A_hX)\mathbf c_\mu=Y^\dagger\mathbf b_\mu,\qquad
\mathbf U_\mu=X\mathbf c_\mu.
\]
The implementation uses a Schur complement of the full coupled system. The kernel-mode lower-left block may not be discarded. `raw_base_block`/`raw_dictionary_block` are archived diagnostics, not kernel-mode zero/identity gates.

`raw_kernel` stores supported kernel representatives \(\Psi\) for inheritance. After projection/POD the checked relation is \(\mathscr J_\ell\Phi\approx\Psi\); final trial modes need not be locally supported or in \(\ker I_H\). `values` stores fine-coordinate solution columns ordered by `computed_ids`, not reduced coefficients. `warm_full` stores the coarse-coordinate localization eigensolver warm block, not the ALOD solution or LU factors.

Local Riesz and corrector solves impose kernel constraints through saddle/Schur systems. An exactly zero local kernel contributes exactly zero Riesz solutions. Residual/constraint gates test these finite-dimensional equations, not true solution error.

## Version 7 symbol mapping

Use these symbols in new explanations without renaming code fields or changing frozen experiment values:

| Mathematical object | Version 7 notation | Literal implementation names / historical notation |
| --- | --- | --- |
| Reference parameter | \(\mu_{\rm ref}\) for the corner; \(\mu_0\) for Wavenumber dependence | `nominal`, reference member ID; older \(\mu_{\rm nom}\); verify IDs per run |
| Packet centre components | \(\mu_x,\mu_y\) | `wave_x`, `wave_y`; older \(z_x,z_y\) or \(z_\mu\) |
| E2 singular coefficient | \(\mu_c\) | `singular_coefficient`; older \(c_\mu\) |
| E2 wave amplitude / phase | \(\mu_\alpha,\mu_\varphi\) | `wave_amplitude`, `wave_phase`; older \(\alpha_\mu,\varphi_\mu\) |
| Corner polar radius / angle | \(\varrho,\theta\) | Older \(r,\vartheta\); do not confuse radius with enrichment rank \(r\) |
| Base solution and spaces | \(U_{H,\ell,\mu},X_{H,\ell},Y_{H,\ell}\) | Older \(U_\mu^0,X^0,Y^0\) |
| Compressed ALOD solution | \(U_{\mu,r}:=U_\mu(\mathcal E_r)\) | Stored fine-coordinate solution columns; \(r=\dim\mathcal E_r\) |
| Regional direction block length | \(s\) | Older \(s_{\rm AS}\); inspect the actual block/budget option |
| Primal / adjoint correctors | \(\mathcal C_\infty,\mathcal C_\ell,\mathcal C_\ell^*\) | Older \(\mathcal Q\); source macros `\Qir`/`\Qlr` now expand to \(\mathcal C\) |
| Energy form / resolution scale | \(b\), \(\chi\) | Older \(b_k\), \(\chi_{\rm res}\); the energy norm still has subscript \(k\) |
| Balance indicator | \(q_\ell=d_{\mu_{\rm ref}}\Theta_\ell/\eta_{H,\mu_{\rm ref}}\) | Older \(\widetilde q_\ell\); `ell_ratio_mode=solution_scaled` |
| Energy projection / orthogonal complement | \(P_Z\), \(Z^\perp\), defined with \(b\) | Older \(P_Z^k\), \(Z^{\perp_k}\); distinct from the proof's Ritz projections |
| Enriched localization constant | \(C_{\rm enr,\ell}\) | Older \(C_{\rm loc,\ell}^{\rm enr}\) |
| Stability / reliability margin constants | \(\rho_{\rm st},C_{\rm rb}\) | Older \(\rho_{\rm stab},C_{\rm rob}\); not RHS sampling radius |

A symbol mapping does not change a stored family, seed, training/test membership or configuration. Preserve literal `nominal`, `sample`, `training_ids`, `member_ids` and `SLOD` identifiers. In explanatory prose use reference load, test load and LOD; `SLOD` is this repository's LOD baseline name, not the super-localized LOD method cited in the introduction.

## Version 7 experiments and historical identifiers

- **A reentrant corner** corresponds to repository/server E2: \(\mu=(\mu_c,\mu_\alpha,\mu_x,\mu_y,\mu_\varphi)\), reference \((1,1/2,-1/2,1/2,0)\), and \(D=\Omega\cap B_{0.6}(0)\). The 16 training loads comprise the reference and 15 randomized loads; 24 test loads are excluded from training and marking. Amplitude ranges are \(\mu_c\in[0.25,1.75]\), \(\mu_\alpha\in[0.15,1]\); phase is \([-\pi,\pi]\); the center lies in the radius-0.2 disk about \((-1/2,1/2)\). See `sec:num-lshape` and [experiments](../../../../docs/experiments.md).
- **Wavenumber dependence** corresponds to E3 and uses the E1 manufactured problem at fixed center \(\mu_0=(3/4,1/2)\) for \(k=8,16,32,64,128\). Each run is single-load, with its own operator/correctors and no regional enrichment. It is not cross-wavenumber shared-space reuse. The initial scales satisfy \(kH_0=2\sqrt2\), \(k^3h_0^2=8\).
- E1 remains a smooth-packet code problem and historical experiment, but is absent as a separate experiment in version 7. Do not rename old directories or infer their data were deleted.
- Checked-in formal RHS tables contain 48 IDs (16 `train`, 24 `test`, 8 `shift`), plus optional E2 `pure` controls. This does not establish that those exact tables generated version 7. Its manuscript states randomized loads but gives no numerical seed or generator; obtain them from the particular run's frozen table/configuration. `training_ids`, not the `train` label alone, controls training/marking selection; `member_ids` and audit selection have separate roles.

## Fixed-mesh control terminology

All corner controls in `tab:e2-fixedmesh-control` use the final mesh pair from the **25-cycle, fixed-\(\ell=3\) family-aware control**, not the 28-cycle main run's final mesh pair.

| Paper label | What is fixed or rebuilt |
| --- | --- |
| Base | No enrichment, \(\ell=3,r=0\) |
| Retrained shared | Rebuild from all 16 training loads on the fixed mesh pair separately at \(\ell=2\) and \(3\); both retain rank two |
| Adaptive shared | Freeze the original adaptive rank-two enrichment across all evaluated loads at \(\ell=3\); recompute only PG coefficients |
| Per-load fresh | Build a new enrichment using each evaluated load at \(\ell=3\); rank one in these controls; diagnostic rather than shared-space reuse |

The single-load/family-aware trajectories are independent 25-cycle runs at fixed \(\ell=3\), assessed on the same 24 test loads. The reported maximum errors at \(N_{\rm on}=10^3\) use logarithmic interpolation, not first-crossing values. The fixed-mesh oversampling comparison rebuilds the enrichment as well as changing \(\ell\); it does not keep an old enrichment unchanged while changing only patch depth.

## Known paper/current-code distinctions

These are differences in definitions or scope, not authorization to change either code or frozen experiments. Recheck them when the sources change.

| Topic | Manuscript | Reviewed operational implementation |
| --- | --- | --- |
| Cycle counting | Algorithm 2 solves after one coarse update and each of \(m_{\rm ref}\) fine-only updates | There are \(m_{\rm ref}\) accepted sweeps per cycle: first commits coarse+fine updates, subsequent sweeps are fine-only. Initial state plus 28 cycles at `m_ref=2` gives 57 accepted states. The literal paper loop has one extra solve per cycle; do not identify the counts. |
| Absolute localization threshold | Corner: \(0.2\); Wavenumber dependence: \(3.2/k\) | The generator stores the manual scalar \(3.2/k\); no runtime calibration is inferred. These agree at \(k=16\). |
| Localization timing | Algorithm 1 includes a pre-training absolute-defect loop and post-training balance step | Main controller uses lazy scheduled checks, overlaps \(\Theta\) work with RHS/training, then applies absolute/relative gates and same-mesh rechecks. This differs from an unconditional pre-training loop. |
| Compatibility | Paper uses kernel lifting and scaled balance | Main presets select `kernel_lift`, `solution_scaled`, relative threshold 0.3; missing fields in old Python configs preserve archived adjoint/raw policies. Inspect direct CLI/library defaults separately. |
| Accuracy stopping | Experiments report exact-error crossings | E3 requests reference-load relative energy error at most 0.01 (`exact_scope=nominal`) with resource/horizon bounds. A bound reached without a crossing is not target success. Optional manufactured-error stopping is additional to residual-driven marking. |
| Fit/budget convention | Corner tail windows: 10 ALOD cycles, 10 AFEM states, 8 UFEM states, 4 LOD states; family control uses interpolation at \(10^3\) | The version-4 delivery inventory and September-28 scoped exporter retain older roles/windows and are not version-7 reproduction specifications. |
| Run horizon and ell cap | Corner: 28 cycles, cap 4; Wavenumber dependence: at most 40 cycles, cap 5 | Current parsers/checkpoint loading support ell up to 5, but checked-in E3 presets and the production generator retain older horizons/caps. Verify frozen run configurations instead of substituting defaults. |

Do not silently rewrite abstract paper algorithms to match every optional code control, or describe an archived AOT experiment as kernel-lifted. Current lazy checks are initial, terminal, configured extra states and intervals of about \(2^\ell\) accepted states. `ell_cap_reached` reports an unmet promotion request at the cap. Missing `Theta` means not computed for that state, not zero. Reuse needs the source's checked mesh/operator/basis identities.

## Numerical audit metrics

For full ALOD audits in the unchanged accepted state,
\[
e_{n,\mu}=\|u_\mu-U_{n,\mu}\|_k,\quad
f_{n,\mu}=\|u_\mu-u_{h_n,\mu}\|_k,\quad
g_{n,\mu}=\|u_{h_n,\mu}-U_{n,\mu}\|_k.
\]
`e`/`exact_error`, `f`, `g` are absolute quantities. `E`, `F`/`E_ref`, `G` divide each by `exact_norm` \(=\|u_\mu\|_k\). The uppercase audit field `F` is a normalized error, not the RHS functional \(F_\mu\). Do not assert \(e=f+g\); only the norm inequality \(e\leq f+g\) follows.

The paper's reference-load exact error \(E_{\mu_{\rm ref}}\) (\(E=E_{\mu_0}\) in Wavenumber dependence) evaluates the exact-error metric for a specified RHS case. It is not the literal code field `E_ref`, which aliases normalized fine-grid discretization error. Always disambiguate a reference RHS case from a reference fine-grid Galerkin solution.

`energy` is \(\|U\|_k\), despite its field name, not its square or necessarily physical wave energy. `accepted_solution_distance` compares with the saved accepted column if present. `PG_residual` checks projected algebraic equations, not \(e,f,g\).

Exact-only audits avoid auxiliary reference FEM solves. Reference fields are `null`/`not_computed`, not zero; for AFEM some comparisons are not applicable. Saved accepted columns may be reused; additional RHS cases need frozen-space solves. Diagnostics that deepen meshes, change ell/rank or retrain define different controlled quantities and must be labelled accordingly.

The two-level diagnostic is an energy Frobenius norm of corrected-basis differences at \(\ell\) and \(\ell+1\), not automatically \(\delta_\ell\), \(\Theta_\ell\), or a certified bound. Describe region selection and quadrature controls precisely.

## Saved states, disk payloads and retention

An audit snapshot preserves numerical inputs for an accepted \(S_n\); a checkpoint preserves numerical and algorithm state for restart. They share serialization but differ in publication/retention role. `CheckpointPhase::Accepted` is auditable; `BeforeTraining`/`AfterEll` may be resumable without being published audit states.

The durable record contains meshes/boundary tags/stable identities, embeddings `P_node`/`P_elem`/`P_dg` (possibly in a shared geometry object), `raw_kernel`/`phi`/`values`/`warm_full`, `computed_ids`, ell/cursor/phase/check schedule, refinement marks, configuration/RHS table, mathematics/space identity and journal consistency metadata. Optional persisted LOD trial basis and base reduced operator avoid rebuilding accepted correctors.

It does not serialize every assembled fine-grid stiffness, mass, Helmholtz or energy matrix, or their LU factors. Required matrices/factors are reconstructed or factored in a new process. Linux shared-memory audits can share additional assembled in-memory objects; distinguish that transfer path from disk serialization.

Temporary records are verified before publication; restart pointers and journals follow a defined commit order. Publication does not imply a background audit has completed. Periodic restart checkpoints (`checkpoint_interval_cycles`, currently five in main presets) and per-state audit snapshots have different frequencies; the interval does not limit audit snapshot count.

| Requested operation | Required data |
| --- | --- |
| Replot existing audit results | Complete metrics plus configuration/provenance; canonical delivery data may suffice |
| Rerun an unchanged state's audit | Its auditable numerical inputs and matching geometry/configuration in a snapshot or accepted checkpoint |
| Evaluate additional RHS cases | Frozen spaces/bases, operator/geometry and RHS definitions; saved solution columns alone generally do not suffice |
| Change mesh, ell, rank or training policy diagnostically | Inputs and rebuilding/retraining specified by that diagnostic |
| Regenerate deleted intermediate state | Replay from an adequate restart point with matching code/configuration and validation |
| Continue adaptation | Valid resumable checkpoint including cursor, inheritance and consistency data |

Deleting an audit snapshot leaves existing metrics intact and does not remove the mathematical ability to audit. It removes that file's direct input for its state; a retained accepted checkpoint may still supply it. Replay requires computation and matching inputs. Which intermediate states must remain directly re-auditable is a retention decision.

Before estimating reclaimed storage, account for hard links shared with checkpoints, shared geometry dependencies and outstanding workers. Logical sizes cannot be added to get physical use; removing one hard link may free no blocks. Distinguish historical supplied measurements from fresh observations. Determine dominant payloads by measurement. Retention advice is not authorization to delete experiments or change persistence settings.
