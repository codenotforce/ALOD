# Implementation correspondence and operational terminology

Reviewed on 2026-10-05. Treat current source, frozen run configuration and dated production evidence separately. Repository paths below are relative to its root; manuscript labels refer to `../LOD_paper/ALOD/version4/arxiv.tex`.

Contents: source map; coordinate representations; RHS experiments; paper/code differences; audit metrics; serialization and retention.

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
\mathbf U_\mu^0=B\mathbf c_\mu.
\]
\(\dagger\) denotes matrix conjugate transpose; it is distinct from the paper's star on adjoint correctors. Current complex-symmetric assembly gives \(C=\overline B\), without making \(\mathbf A_h\) Hermitian. `LodSpace::trial()`, `test()`, `reduced()` represent \(B,C,\mathbf A_{\rm red}\).

For `phi` (enrichment basis matrix) \(\Phi\), lifted tests have columns \(Z=\Phi-B(I_H\Phi)_{\rm free}\). With \(X=[B,\Phi]\), \(Y=[C,Z]\), solve
\[
(Y^\dagger\mathbf A_hX)\mathbf c_\mu=Y^\dagger\mathbf b_\mu,\qquad
\mathbf U_\mu=X\mathbf c_\mu.
\]
The implementation uses a Schur complement of the full coupled system. The kernel-mode lower-left block may not be discarded. `raw_base_block`/`raw_dictionary_block` are archived diagnostics, not kernel-mode zero/identity gates.

`raw_kernel` stores supported kernel representatives \(\Psi\) for inheritance. After projection/POD the checked relation is \(\mathscr J_\ell\Phi\approx\Psi\); final trial modes need not be locally supported or in \(\ker I_H\). `values` stores fine-coordinate solution columns ordered by `computed_ids`, not reduced coefficients. `warm_full` stores the coarse-coordinate localization eigensolver warm block, not the ALOD solution or LU factors.

Local Riesz and corrector solves impose kernel constraints through saddle/Schur systems. An exactly zero local kernel contributes exactly zero Riesz solutions. Residual/constraint gates test these finite-dimensional equations, not true solution error.

## E1/E2/E3 RHS terminology

- E1: smooth localized oscillations. \(\mu=z=(z_x,z_y)\) varies the packet centre near \((3/4,1/2)\); phase is \(e^{\mathrm i kx}\). `wave_x`/`wave_y` are centre coordinates, not propagation directions. E1 has no regional enrichment.
- E2: an L-shaped domain with a reentrant-corner term \(r^{2/3}\sin(2\vartheta/3)\), boundary factor \((1-x^2)^2(1-y^2)^2\), and a normalized localized oscillatory packet. Parameters vary singular coefficient \(c_\mu\), centre \(z_\mu\), amplitude \(\alpha_\mu\), and phase \(\varphi_\mu\). Main enrichment uses \(D=\Omega\cap B_{0.6}(0)\). Its mask requires complete patches inside \(D\); a physical-region error diagnostic instead selects triangles by centroid. Distinguish these two region definitions.
- The paper has 16 training and 24 held-out RHS cases at \(k=16\). Current formal JSON tables have 48 IDs: 16 `train`, 24 `test`, 8 `shift`. Two additional E2 `pure` controls may be requested. Do not call all 48 cases training or equate every audit case with a held-out case.
- `training_ids` selects training/marking cases; a `split=train` label alone does not prove selection. `member_ids` and audit selection have separate roles. Current native audit uses ID 0 as nominal; verify the specific run's IDs/selection.
- E3: the E1 nominal problem at \(k=8,16,32,64,128\), currently ID 0 per preset. These are separate operators, not a single fixed-operator RHS family or cross-wavenumber reuse. Initial scaling is \(kH_0=2\sqrt2\), \(k^3h_0^2=8\), hence \(h_0\propto k^{-3/2}\).
- Frozen held-out evaluation uses the retained shared space. Independent per-RHS retraining (`fresh`) is a separate diagnostic.

## Known paper/current-code distinctions

Recheck these observations when the code or paper revision changes.

| Topic | Manuscript | Reviewed operational implementation |
| --- | --- | --- |
| Cycle counting | Algorithm 2 solves after one coarse update and each of \(m_{\rm ref}\) fine-only updates | There are \(m_{\rm ref}\) accepted sweeps per cycle: first commits coarse+fine updates, subsequent sweeps are fine-only. Initial state plus 28 cycles at `m_ref=2` gives 57 accepted states. The literal paper loop has one extra solve per cycle; do not identify the counts. |
| Absolute localization threshold | Settings table gives \(\tau_\Theta=0.2\) | Main presets/production use \(3.2/k\), equal to 0.2 only at \(k=16\). |
| Localization timing | Algorithm 1 includes a pre-training absolute-defect loop and post-training balance step | Main controller uses lazy scheduled checks, overlaps \(\Theta\) work with RHS/training, then applies absolute/relative gates and same-mesh rechecks. This differs from an unconditional pre-training loop. |
| Compatibility | Paper uses kernel lifting and scaled balance | Main presets select `kernel_lift`, `solution_scaled`, relative threshold 0.3; missing fields in old Python configs preserve archived adjoint/raw policies. Inspect direct CLI/library defaults separately. |
| Accuracy stopping | Experiments report exact-error crossings | E3 requests nominal relative exact-energy error at most 0.01 with resource/horizon bounds. A bound reached without a crossing is not target success. Optional manufactured-error stopping is additional to residual-driven marking. |
| Fit/budget convention | Specified endpoints/tails; one E1 prose budget differs from published interpolation | Exporters distinguish interpolation from first crossing, and general paper delivery from scoped 25-cycle plots / cycles 19–28 fits. State the actual convention. |

Do not silently rewrite abstract paper algorithms to match every optional code control, or describe an archived AOT experiment as kernel-lifted. Current lazy checks are initial, terminal, configured extra states and intervals of about \(2^\ell\) accepted states. `ell_cap_reached` reports an unmet promotion request at the cap. Missing `Theta` means not computed for that state, not zero. Reuse needs the source's checked mesh/operator/basis identities.

## Numerical audit metrics

For full ALOD audits in the unchanged accepted state,
\[
e_{n,\mu}=\|u_\mu-U_{n,\mu}\|_k,\quad
f_{n,\mu}=\|u_\mu-u_{h_n,\mu}\|_k,\quad
g_{n,\mu}=\|u_{h_n,\mu}-U_{n,\mu}\|_k.
\]
`e`/`exact_error`, `f`, `g` are absolute quantities. `E`, `F`/`E_ref`, `G` divide each by `exact_norm` \(=\|u_\mu\|_k\). The uppercase audit field `F` is a normalized error, not the RHS functional \(F_\mu\). Do not assert \(e=f+g\); only the norm inequality \(e\leq f+g\) follows.

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
