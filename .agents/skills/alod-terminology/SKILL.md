---
name: alod-terminology
description: Explain this ALOD project using English technical terms, Chinese explanatory prose and the version 7 arxiv.tex/main.tex paper's mathematical notation. Use for algorithms, right-hand-side families, enrichment, adaptivity, experiments, error analysis, audits, checkpoints, and storage.
---

# ALOD terminology and mathematical explanations

Use the version 7 paper's mathematical objects and notation when explaining this project. Answer in the user's language, keeping specialized mathematical, algorithmic and project terms in English even in Chinese responses. Use the canonical English terms below consistently. Explicit user terminology takes precedence. This skill governs explanations, not permission to alter algorithms, experiments, or saved data.

## Sources and scope

Locate the repository root containing `include/alod/checkpoint.hpp` and `apps/fixed_support.hpp`. Source paths below are relative to that root. The manuscript sources are `../LOD_paper/ALOD/version7/arxiv.tex` and `../LOD_paper/ALOD/version7/main.tex`. Both were checked: their mathematical and experiment bodies agree; `main.tex` uses review mode and BibTeX, while `arxiv.tex` embeds its bibliography. Neither is an obsolete source. Read TeX macro definitions before interpreting notation: `\Qir` and `\Qlr` now expand to `\mathcal C_\infty` and `\mathcal C_\ell`. If later copies diverge, identify the source and difference rather than silently combining them.

For mathematics read [mathematical notation](references/mathematical-notation.md); for code fields, experiments, audits or storage read [implementation correspondence](references/implementation-correspondence.md). These references remain usable if the sibling manuscript is absent. Verify mutable details against current source or the particular run's frozen configuration.

Version 7 has two experiments, **A reentrant corner** and **Wavenumber dependence**. Repository and server identifiers E1/E2/E3 are retained: E2 corresponds to the corner experiment and E3 to the wavenumber study; E1 remains a code problem and historical experiment, not a third experiment in the current paper. Manuscript settings do not update checked-in presets, deployment files or frozen results. Keep those sources distinct.

## How to answer

- Define the mathematical object first, then explain its implementation. Prefer a short defining equation followed by plain Chinese interpreting it. Introduce only symbols needed for the question; mathematical language should clarify, not inflate the answer.
- Keep technical names in English: for example, right-hand-side functional, quasi-interpolation, Riesz representative, Petrov–Galerkin (PG), enrichment space, energy-weighted POD, audit snapshot and checkpoint. Explain their meaning in Chinese without repeatedly adding Chinese translated labels. Ordinary connective prose remains Chinese; do not turn an entire answer into English just because it contains technical terms.
- Prefer the paper's English term where available; expand unfamiliar abbreviations once, then use them consistently. Preserve established capitalization, symbols and literal code identifiers. If the user asks for a translation, provide it explicitly rather than changing the default naming convention.
- Follow the paper's symbols: parameter \(\mu\), RHS functional \(F_\mu\), meshes \(\mathcal T_H,\mathcal T_h\), kernel \(W_h\), correctors \(\mathcal C_\ell,\mathcal C_\ell^*\), base spaces \(X_{H,\ell},Y_{H,\ell}\), enrichment \(\mathcal E_r\), and compressed ALOD solution \(U_{\mu,r}:=U_\mu(\mathcal E_r)\). Use \(U_{H,\ell,\mu}\) for the base reconstruction and \(\mu_{\rm ref}\) for the corner reference RHS case (\(\mu_0\) is the fixed packet-center parameter in the wavenumber study). Reserve \(u_{h,\mu}\) for the fine-grid Galerkin solution. Prefer these current symbols over the older aliases \(X^0,Y^0,U_\mu^0,\mu_{\rm nom}\).
- Distinguish a parameter, function, functional, finite element function, coefficient vector, matrix and subspace. Identify introduced matrix notation as a coordinate representation; preserve literal code identifiers.
- Specify the mesh, RHS parameter, norm and normalization when discussing an error or indicator. State assumptions before invoking a theorem. A residual gate, Ritz convergence or empirical error curve does not alone establish the paper's stability/reliability hypotheses.
- For 样本, 状态, 快照 or 残差, identify the actual object before using the word. Avoid multiplying vague synonyms.
- Describe what changes and what is fixed. Use \((n,\mu)\) for adaptive-state/RHS-case combinations. Keep `state_id`, `sample` and other code/CLI identifiers unchanged.

## Mandatory terminology

| Object | Canonical English term and notation |
| --- | --- |
| `sample`, `member`, `FamilyMember` in an RHS table | right-hand-side case (RHS case), indexed by \(\mu\); a parameterized problem instance, not merely a vector |
| parameter set / selected training cases | parameter set \(\mathcal P\) / training parameter set \(\mathcal P_{\rm tr}\); corresponding loads form the training RHS family |
| test / reference load | test load / reference load, indexed by \(\mu\) / \(\mu_{\rm ref}\); use test-load / reference-load as modifiers; preserve literal `nominal` fields |
| source / load functional / assembled vector | source function \(f_\mu\) / right-hand-side functional \(F_\mu\) / assembled RHS vector \(\mathbf b_\mu\) |
| exact / reference FEM / ALOD solution | exact solution \(u_\mu\) / fine-grid Galerkin solution \(u_{h,\mu}\) / compressed ALOD solution \(U_{\mu,r}=U_\mu(\mathcal E_r)\) |
| coarse / reference fine / auxiliary mesh | coarse mesh \(\mathcal T_H\) / reference fine mesh \(\mathcal T_h\) / auxiliary work mesh \(\mathcal T_{\rm w}\) |
| quasi-interpolation / fine kernel | quasi-interpolation operator \(I_H\) / fine-scale kernel \(W_h=\ker(I_H|_{V_h})\) |
| energy form / resolution scale | energy inner product \(b\), with \(\|v\|_k^2=b(v,v)\) / resolution parameter \(\chi\) |
| energy projection / orthogonal complement | \(P_Z\) / \(Z^\perp\), both defined using \(b\); the stability proof's Ritz projections are different |
| primal / adjoint corrector | primal corrector operator \(\mathcal C_\ell\) / adjoint corrector operator \(\mathcal C_\ell^*\); their values are corrector functions |
| base / enriched trial and test spaces | base spaces \(X_{H,\ell},Y_{H,\ell}\); enriched spaces \(X_{H,\ell}\oplus\mathcal E_r\), \(Y_{H,\ell}\oplus\mathscr J_\ell\mathcal E_r\); unrelated to training/test RHS sets |
| local Riesz representative | local residual Riesz representative \(\xi_z(R_{\mu,U})\) |
| coarse / regional residual indicator | coarse residual indicator \(\eta_{H,\mu}\) / regional residual indicator \(\eta_{D,\mu}\) |
| fine strong-residual indicator | fine-grid strong-residual indicator \(\zeta_{K,\mu}\) |
| oversampling / localization defect indicator | oversampling level \(\ell\) / localization indicator \(\Theta_\ell\) |
| localization defect / threshold / balance | corrector defect \(\mathcal D_\ell=\mathcal C_\infty-\mathcal C_\ell\), its operator norm \(\delta_\ell\) / threshold \(\tau_\Theta(k)\) / balance indicator \(q_\ell\) |
| Dörfler fractions | Dörfler marking parameters \(\theta_H,\theta_h\); lowercase \(\theta\) is not \(\Theta_\ell\) |
| shared enrichment / stored dictionary | shared enrichment space \(\mathcal E_r\) / enrichment basis matrix \(\Phi\); avoid an unexplained dictionary |
| `raw_kernel` / kernel-lifted tests | raw kernel correction matrix \(\Psi\) / kernel-lifted test functions \(\mathscr J_\ell\phi\) |
| POD / rank | energy-weighted proper orthogonal decomposition (POD) / enrichment dimension \(r=\dim\mathcal E_r\) |
| regional direction block | projected Schwarz block of at most \(s\) directions; \(s\) is not the retained rank \(r\) |
| adaptive / accepted state | adaptive state \(S_n\) / accepted state |
| coarse cycle / reference sweep | coarse refinement cycle / reference refinement sweep; verify the controller's counting convention |
| audit / audit snapshot / checkpoint | numerical audit / audit snapshot / checkpoint |
| POD snapshot | POD solution-difference snapshot \(w_\mu\); distinct from an audit snapshot |

In current prose use single load, reference load and test load as nouns, and single-load, reference-load and test-load as modifiers. Use family-aware for the shared training/adaptation policy. Call the comparison method LOD; preserve the literal code method `SLOD`, which does not mean the distinct super-localized LOD method cited in the paper.

The RHS functional is general; a volume source gives \(F_\mu(v)=(f_\mu,v)_{L^2(\Omega)}\). A manufactured exact solution defines the loads and evaluates errors in the paper's experiments (implemented through E1/E2), but is not a prerequisite for the general algorithm. Do not define every RHS case as necessarily containing an analytic solution.

## Essential distinctions

- Separate \(u_\mu\), \(u_{h,\mu}\), coarse coefficient function \(u_{H,\ell,\mu}\), and reconstruction \(U_{H,\ell,\mu}=T_\ell u_{H,\ell,\mu}\). Enriched ALOD solves in \(X_{H,\ell}\oplus\mathcal E_r\), not the full \(V_h\).
- \(\Theta_\ell\) depends on the mesh pair, operator, interpolation and \(\ell\), not a chosen RHS. Its square is \(\lambda_{\max}(R_\ell,B_H)\), with uncorrected coarse energy matrix \(B_H\). The balance ratio \(q_\ell=d_{\mu_{\rm ref}}\Theta_\ell/\eta_{H,\mu_{\rm ref}}\) also depends on the reference-case solution. Distinguish the measured indicator \(\Theta_\ell\) from its threshold \(\tau_\Theta(k)\); a missing/deferred \(\Theta\) is not zero.
- Regional support and kernel membership apply to raw corrections. Energy projection can give final trial modes global support and nonzero coarse interpolation. Current enriched tests use \(\mathscr J_\ell=I-T_\ell I_H\); for the paper's regional construction and its compressed subspaces, \(\widehat{\mathcal E}=\mathscr J_\ell\mathcal E\subseteq W_{h,D}\) recovers regional support. Solve the full coupled PG system; this support statement is not a claim about arbitrary supplied enrichment or archived tests.
- The prescribed region \(D\) restricts where unprojected enrichment corrections are generated; global coarse/fine marking still covers \(\Omega\). It is not a local-error target region chosen by the residual.
- POD minimizes an average projection error. Accepted compression also needs memberwise distance and regional-residual checks. A rank cap or training stop is not proof that all targets were met.
- Separate local corrector equation residuals, PG algebraic residuals, Riesz indicators, localization defects and solution errors. Keep fine-reference error distinct from exact error:
  \[
  \|u_\mu-U_\mu\|_k\leq\|u_\mu-u_{h,\mu}\|_k+\|u_{h,\mu}-U_\mu\|_k.
  \]
- Adaptive shared freezes the enrichment retained from the adaptive run; Retrained shared rebuilds it from the whole training family on a fixed mesh pair; Per-load fresh constructs a new space using each evaluated load and is a diagnostic rather than shared-space reuse.
- A frozen audit fixes accepted meshes and trial/test/enrichment spaces while evaluating specified RHS cases. It does not halt the adaptive solver or imply all test-case solutions were saved.
- Disk checkpoints and audit snapshots use common serialization with different roles. Do not claim they save every fine-grid matrix or LU factor. Replotting metrics, rerunning audits and replaying adaptation require different retained data.
- Separate 主求解完成, 数值校验通过, 数值审计完成 and 论文交付完成. An accepted state or horizon stop is not automatically a reached accuracy target or completed audit.

Example: 参数 \(\mu\) 确定 right-hand-side functional \(F_\mu\)。在 adaptive state \(S_n\) 上，程序装配 RHS vector \(\mathbf b_{n,\mu}\)，并在固定的 trial space \(X_{H_n,\ell_n}\oplus\mathcal E_{r_n}\) 中求得 ALOD solution \(U_{n,\mu}\)。这里的 \(n\) 是为说明 adaptive state 而补充的下标。进入下一步时，\(\mu\) 可以不变，但 mesh、space、RHS vector 和 solution 通常会变化。
