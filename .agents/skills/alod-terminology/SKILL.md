---
name: alod-terminology
description: Explain this ALOD project using English technical terms, Chinese explanatory prose and the version 4 paper's mathematical notation. Use for algorithms, right-hand-side families, enrichment, adaptivity, experiments, error analysis, audits, checkpoints, and storage.
---

# ALOD terminology and mathematical explanations

Use the version 4 paper's mathematical objects and notation when explaining this project. Answer in the user's language, keeping specialized mathematical, algorithmic and project terms in English even in Chinese responses. Use the canonical English terms below consistently. Explicit user terminology takes precedence. This skill governs explanations, not permission to alter algorithms, experiments, or saved data.

## Sources and scope

Locate the repository root containing `include/alod/checkpoint.hpp` and `apps/fixed_support.hpp`. Source paths below are relative to that root. The reviewed manuscript is `../LOD_paper/ALOD/version4/arxiv.tex`; the definitions here remain usable if that sibling checkout is absent.

For mathematics read [mathematical notation](references/mathematical-notation.md), especially the relevant definitions and assumptions. For code fields, experiments, audits or storage read [implementation correspondence](references/implementation-correspondence.md). Verify mutable details against current source or the particular run's frozen configuration. Distinguish the paper, current implementation and historical experiment evidence whenever they differ.

## How to answer

- Define the mathematical object first, then explain its implementation. Prefer a short defining equation followed by plain Chinese interpreting it. Introduce only symbols needed for the question; mathematical language should clarify, not inflate the answer.
- Keep technical names in English: for example, right-hand-side functional, quasi-interpolation, Riesz representative, Petrov–Galerkin (PG), enrichment space, energy-weighted POD, audit snapshot and checkpoint. Explain their meaning in Chinese without repeatedly adding Chinese translated labels. Ordinary connective prose remains Chinese; do not turn an entire answer into English just because it contains technical terms.
- Prefer the paper's English term where available; expand unfamiliar abbreviations once, then use them consistently. Preserve established capitalization, symbols and literal code identifiers. If the user asks for a translation, provide it explicitly rather than changing the default naming convention.
- Follow the paper's symbols: parameter \(\mu\), RHS functional \(F_\mu\), meshes \(\mathcal T_H,\mathcal T_h\), kernel \(W_h\), correctors \(\mathcal Q_\ell,\mathcal Q_\ell^*\), enrichment \(\mathcal E_r\), and ALOD solution \(U_\mu(\mathcal E_r)\). Reserve \(u_{h,\mu}\) for the fine-grid Galerkin solution; do not use it generically for every numerical solution.
- Distinguish a parameter, function, functional, finite element function, coefficient vector, matrix and subspace. Identify introduced matrix notation as a coordinate representation; preserve literal code identifiers.
- Specify the mesh, RHS parameter, norm and normalization when discussing an error or indicator. State assumptions before invoking a theorem. A residual gate, Ritz convergence or empirical error curve does not alone establish the paper's stability/reliability hypotheses.
- For 样本, 状态, 快照 or 残差, identify the actual object before using the word. Avoid multiplying vague synonyms.
- Describe what changes and what is fixed. Use \((n,\mu)\) for adaptive-state/RHS-case combinations. Keep `state_id`, `sample` and other code/CLI identifiers unchanged.

## Mandatory terminology

| Object | Canonical English term and notation |
| --- | --- |
| `sample`, `member`, `FamilyMember` in an RHS table | right-hand-side case (RHS case), indexed by \(\mu\); a parameterized problem instance, not merely a vector |
| parameter set / selected training cases | parameter set \(\mathcal P\) / training parameter set \(\mathcal P_{\rm tr}\); corresponding loads form the training RHS family |
| held-out / nominal load | held-out RHS case / nominal RHS case \(\mu_{\rm nom}\) |
| source / load functional / assembled vector | source function \(f_\mu\) / right-hand-side functional \(F_\mu\) / assembled RHS vector \(\mathbf b_\mu\) |
| exact / reference FEM / ALOD solution | exact solution \(u_\mu\) / fine-grid Galerkin solution \(u_{h,\mu}\) / ALOD solution \(U_\mu(\mathcal E_r)\) |
| coarse / reference fine / auxiliary mesh | coarse mesh \(\mathcal T_H\) / reference fine mesh \(\mathcal T_h\) / auxiliary work mesh \(\mathcal T_{\rm w}\) |
| quasi-interpolation / fine kernel | quasi-interpolation operator \(I_H\) / fine-scale kernel \(W_h=\ker(I_H|_{V_h})\) |
| primal / adjoint corrector | primal corrector operator \(\mathcal Q_\ell\) / adjoint corrector operator \(\mathcal Q_\ell^*\); their values are corrector functions |
| trial / test space | trial space \(X\) / test space \(Y\); unrelated to training/test RHS sets |
| local Riesz representative | local residual Riesz representative \(\xi_z(R_{\mu,U})\) |
| coarse / regional residual indicator | coarse residual indicator \(\eta_{H,\mu}\) / regional residual indicator \(\eta_{D,\mu}\) |
| fine strong-residual indicator | fine-grid strong-residual indicator \(\zeta_{K,\mu}\) |
| oversampling / localization defect indicator | oversampling level \(\ell\) / localization indicator \(\Theta_\ell\) |
| Dörfler fractions | Dörfler marking parameters \(\theta_H,\theta_h\); lowercase \(\theta\) is not \(\Theta_\ell\) |
| shared enrichment / stored dictionary | shared enrichment space \(\mathcal E_r\) / enrichment basis matrix \(\Phi\); avoid an unexplained dictionary |
| `raw_kernel` / kernel-lifted tests | raw kernel correction matrix \(\Psi\) / kernel-lifted test functions \(\mathscr J_\ell\phi\) |
| POD / rank | energy-weighted proper orthogonal decomposition (POD) / enrichment dimension \(r=\dim\mathcal E_r\) |
| adaptive / accepted state | adaptive state \(S_n\) / accepted state |
| coarse cycle / reference sweep | coarse refinement cycle / reference refinement sweep; verify the controller's counting convention |
| audit / audit snapshot / checkpoint | numerical audit / audit snapshot / checkpoint |
| POD snapshot | POD solution-difference snapshot \(w_\mu\); distinct from an audit snapshot |

The RHS functional is general; a volume source gives \(F_\mu(v)=(f_\mu,v)_{L^2(\Omega)}\). A manufactured exact solution is available for E1/E2, not a prerequisite for the general algorithm. Do not define every RHS case as necessarily containing an analytic solution.

## Essential distinctions

- Separate \(u_\mu\), \(u_{h,\mu}\), coarse coefficient function \(u_{H,\ell,\mu}\), and reconstruction \(U_{H,\ell,\mu}=T_\ell u_{H,\ell,\mu}\). Enriched ALOD solves in \(X(\mathcal E_r)\), not the full \(V_h\).
- \(\Theta_\ell\) depends on the mesh pair, operator, interpolation and \(\ell\), not a chosen RHS. The balance ratio \(\widetilde q_\ell=d_{\mu_{\rm nom}}\Theta_\ell/\eta_{H,\mu_{\rm nom}}\) also depends on the nominal solution. A missing/deferred \(\Theta\) is not zero.
- Regional support and kernel membership apply to raw corrections. Energy projection can give final trial modes \(\phi_j\) global support and nonzero \(I_H\phi_j\). Current enriched tests use \(\mathscr J_\ell=I-T_\ell I_H\); solve the full coupled PG system.
- POD minimizes an average projection error. Accepted compression also needs memberwise distance and regional-residual checks. A rank cap or training stop is not proof that all targets were met.
- Separate local corrector equation residuals, PG algebraic residuals, Riesz indicators, localization defects and solution errors. Keep fine-reference error distinct from exact error:
  \[
  \|u_\mu-U_\mu\|_k\leq\|u_\mu-u_{h,\mu}\|_k+\|u_{h,\mu}-U_\mu\|_k.
  \]
- A frozen audit fixes accepted meshes and trial/test/enrichment spaces while evaluating specified RHS cases. It does not halt the adaptive solver or imply all test-case solutions were saved.
- Disk checkpoints and audit snapshots use common serialization with different roles. Do not claim they save every fine-grid matrix or LU factor. Replotting metrics, rerunning audits and replaying adaptation require different retained data.
- Separate 主求解完成, 数值校验通过, 数值审计完成 and 论文交付完成. An accepted state or horizon stop is not automatically a reached accuracy target or completed audit.

Example: 参数 \(\mu\) 确定 right-hand-side functional \(F_\mu\)。在 adaptive state \(S_n\) 上，程序装配 RHS vector \(\mathbf b_{n,\mu}\)，并在固定的 trial space \(X_n(\mathcal E_{r_n})\) 中求得 ALOD solution \(U_{n,\mu}\)。进入下一步时，\(\mu\) 可以不变，但 mesh、space、RHS vector 和 solution 通常会变化。
