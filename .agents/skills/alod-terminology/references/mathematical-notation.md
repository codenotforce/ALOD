# Version 7 mathematical notation

Reviewed on 2026-10-19 against both complete version 7 sources, `../LOD_paper/ALOD/version7/main.tex` (SHA-256 `ea0fc85d9739001c81e065b6c7f8f65cb9973c1a7e3de1ce396f4c952ca6cae7`) and `../LOD_paper/ALOD/version7/arxiv.tex` (SHA-256 `1febb1ed15177b1833efb4df60b743ee4e2fdbf173bc90ef48daa11e5c4eefbf`). Their mathematical and experiment bodies agree; document class and bibliography packaging differ. This reference records definitions, not a guarantee about future manuscript or code revisions. Use equation labels rather than unstable line numbers. Experiment scope and code differences are in [implementation correspondence](implementation-correspondence.md).

Contents: variational problem; meshes and LOD; residual and localization; enrichment and POD; stability and best approximation; refinement and errors; adaptive states.

## Variational problem and right-hand-side cases

The paper uses complex-valued Helmholtz solutions with homogeneous Dirichlet, Neumann and impedance boundary conditions:
\[
V=\{v\in H^1(\Omega;\mathbb C):v|_{\Gamma_D}=0\},\qquad
a(u_\mu,v)=F_\mu(v)\quad(v\in V).
\]
Inner products are linear in the first argument and conjugate-linear in the second. The sesquilinear form and energy inner product are
\[
a(v,w)=(\nabla v,\nabla w)_\Omega-k^2(v,w)_\Omega
-\mathrm i k(v,w)_{\Gamma_R},\qquad
b(v,w)=(\nabla v,\nabla w)_\Omega+k^2(v,w)_\Omega,
\quad\|v\|_k^2=b(v,v).
\]
The paper's \(V^\star\) is the continuous anti-dual: \(F_\mu\) is conjugate-linear in its test argument. A volume source gives
\[
F_\mu(v)=(f_\mu,v)_\Omega,\qquad -\Delta u_\mu-k^2u_\mu=f_\mu.
\]
For nodal basis \(\{\lambda_i\}\), \((\mathbf b_\mu)_i=F_\mu(\lambda_i)\). The vector changes with the mesh/basis even if \(\mu\) does not. Do not confuse \(\mathbf b_\mu\) with the energy form \(b\).

A right-hand-side case is the parameterized problem indexed by \(\mu\). A family is \(\{F_\mu:\mu\in\mathcal P\}\); its finite training parameter set is \(\mathcal P_{\rm tr}=\{\mu_1,\ldots,\mu_M\}\). Training constructs shared approximation directions and/or drives family marking; it does not imply random sampling or a neural model. A mathematical parameter need not equal its numeric table ID.

Use \(\mu_{\rm ref}\) for the corner reference parameter and \(\mu_0=(3/4,1/2)\) for the fixed packet center in Wavenumber dependence. These specify loads, not fine-grid reference solutions. The corner parameter is \(\mu=(\mu_c,\mu_\alpha,\mu_x,\mu_y,\mu_\varphi)\). The corner polar radius is \(\varrho\), distinct from enrichment rank \(r\); the polar angle \(\theta\) is distinct from marking fractions \(\theta_H,\theta_h\). These symbol mappings do not change any family's stored parameter values or selection.

At fixed \(k\), coefficients, boundary conditions and mesh pair, the corner family shares the operator and varies the RHS; the code also supports an E1 packet-center family. Across changed \(k\) or meshes this fixed-operator statement no longer applies.

## Meshes, spaces and LOD

Use \(\mathcal T_H\preceq\mathcal T_h\), with conforming piecewise-linear spaces \(V_H\subseteq V_h\subseteq V\). The fine-grid reference solves
\[
a(u_{h,\mu},v_h)=F_\mu(v_h)\qquad(v_h\in V_h).
\]
It is distinct from the exact and ALOD solutions even when all have fine-grid coordinate representations.

The projective quasi-interpolation is
\[
I_H=E_H\Pi_H^{\rm dg},\quad I_H|_{V_H}=I,\quad
W_h=\ker(I_H|_{V_h}),\quad V_h=V_H\oplus W_h.
\]
\(E_H\) uses coarse-element area weights to average DG values at vertices. This direct sum need not be energy-orthogonal. The continuous kernel \(W=\ker I_H\) and discrete \(W_h\) must be distinguished.

The coarse-resolution condition and kernel-coercivity constant are
\[
\chi=C_{\rm app}C_{\rm patch}(1)^{1/2}\|kH\|_{L^\infty(\Omega)}
\leq\chi_0<1,\qquad
c_W=\frac{1-\chi^2}{1+\chi^2}.
\]
Fixed-order patch overlap is bounded by shape regularity, but polynomial growth of \(C_{\rm patch}(m)\) in \(m\) is not automatic on graded meshes and is not required for the a posteriori estimates. See `lem:kernel-coercivity`; verify the mesh-resolution assumption separately from numerical equation residuals.

\(N^m(S)\) denotes an \(m\)-layer coarse patch with adjacency through shared vertices. Correctors use element patches \(N^\ell(T)\); residual Riesz problems use vertex patches \(N^2(z)\). \(\ell\) is neither a mesh level nor a physical radius.

The ideal correctors satisfy
\[
a(\mathcal C_\infty v,w)=a(v,w),\qquad
a(w,\mathcal C_\infty^*v)=a(w,v)\qquad(w\in W_h).
\]
The source macros `\Qir`, `\Qirs`, `\Qlr`, `\Qlrs` expand to the corresponding \(\mathcal C\) symbols. Element correctors \(\mathcal C_{T,\ell},\mathcal C_{T,\ell}^*\) use the restricted form \(a_T\); summing them gives \(\mathcal C_\ell,\mathcal C_\ell^*\). Set
\[
T_\ell=I-\mathcal C_\ell,\quad T_\ell^*=I-\mathcal C_\ell^*,\quad
X_{H,\ell}=T_\ell V_H,\quad Y_{H,\ell}=T_\ell^*V_H.
\]
Ideal correctors act on \(V_h\); localized element sums act on \(V_H\). The paper's star has this variational meaning; it is not automatically the conjugate transpose of a nodal matrix representing \(T_\ell\).

The base PG solution \(U_{H,\ell,\mu}=T_\ell u_{H,\ell,\mu}\) satisfies
\[
a(U_{H,\ell,\mu},T_\ell^*v_H)=F_\mu(T_\ell^*v_H)\quad(v_H\in V_H).
\]
\(u_{H,\ell,\mu}\) is the coarse coefficient function, \(U_{H,\ell,\mu}\) its corrected reconstruction. Keep the explicit version 7 symbols \(X_{H,\ell},Y_{H,\ell},U_{H,\ell,\mu}\) when explaining enrichment rather than introducing the older shorthand \(X^0,Y^0,U_\mu^0\).

## Residual indicators and localization

The residual functional is \(R_{\mu,U}=F_\mu-a(U,\cdot)\). On the local kernel \(W_{h,z}\) supported in the closed patch \(N^2(z)\), define
\[
b(\xi_z(G),w)=G(w)\quad(w\in W_{h,z}),\qquad
\eta(G)^2=\sum_z\|\xi_z(G)\|_k^2.
\]
Hence \(\eta_{H,z,\mu}(U)=\|\xi_z(R_{\mu,U})\|_k\) and \(\eta_{H,\mu}(U)=\eta(R_{\mu,U})\). These involve local energy solves, not Euclidean algebraic residual norms. The additive Schwarz correction is \(\mathscr M_hG=\sum_z\xi_z(G)\), with \(G(\mathscr M_hG)=\eta(G)^2\).

The adjoint coarse-space defect functional and localization indicator are
\[
G_\ell(v_H)(w)=\overline{a(w,T_\ell^*v_H)},\qquad
\Theta_\ell=\sup_{v_H\ne0}\frac{\eta(G_\ell(v_H))}{\|v_H\|_k}.
\]
The corrector defects and their energy operator norm are
\[
\mathcal D_\ell=\mathcal C_\infty-\mathcal C_\ell,\quad
\mathcal D_\ell^*=\mathcal C_\infty^*-\mathcal C_\ell^*,\quad
\delta_\ell=\|\mathcal D_\ell\|_{\mathcal L(V_H,V_h)}
=\|\mathcal D_\ell^*\|_{\mathcal L(V_H,V_h)}.
\]
Call \(\Theta_\ell\) localization indicator and \(\delta_\ell\) localization-error operator norm. They are related, not identical. For an uncorrected coarse basis \(\{\phi_i\}\), version 7 explicitly defines
\[
(B_H)_{ij}=b(\phi_j,\phi_i),\qquad
(R_\ell)_{ij}=\sum_z b\bigl(\xi_z(G_\ell(\phi_j)),\xi_z(G_\ell(\phi_i))\bigr),
\qquad \Theta_\ell^2=\lambda_{\max}(R_\ell,B_H).
\]
\(B_H\) is Hermitian positive definite and \(R_\ell\) is Hermitian positive semidefinite. The denominator is uncorrected coarse energy, not the energy Gram matrix of corrected trial functions. These matrices are distinct from the RHS functional \(F_\mu\), residual functional \(R_{\mu,U}\), and introduced basis matrix \(B\) in the implementation reference.

This load-independent computation belongs to the offline adaptive stage and is shared across RHS cases on a fixed base state. Patch contributions can be computed in parallel; the largest generalized eigenvalue problem is on the coarse space and requires no global fine-grid solve. Local Riesz factorizations and coarse energy factorizations/preconditioners can be reused, with a previous eigenvector as a warm start. Changes to the operator, mesh pair or \(\ell\) invalidate the old value; retained enrichment alone does not. Inspect current lazy-check scheduling before saying the implementation recomputes it at every update.

Under the paper's stable decomposition, synthesis, kernel-coercivity and fine-grid stability hypotheses,
\[
\frac{\Theta_\ell}{C_aC_{\rm syn}}\leq\delta_\ell\leq\frac{C_{\rm sd}}{c_W}\Theta_\ell,
\]
and for a solution with base-test orthogonality,
\[
\frac{\eta_{H,\mu}(U)}{C_aC_{\rm syn}}\leq\|u_{h,\mu}-U\|_k
\leq C_{\rm sd}\left(c_W^{-1}+
\frac{C_IC_FC_{\rm sd}}{\gamma_h(k)c_W}\Theta_\ell\right)\eta_{H,\mu}(U).
\]
See `eq:corrector-certificate` and the global-reliability section. \(\gamma_h(k)\) is the fine-grid inf-sup constant; \(c_W>0\) requires the coarse-resolution condition. An inexact base-test equation adds the paper's \(\epsilon_0\) contribution. Do not claim a computable certified bound when these constants/hypotheses are unverified. A small Ritz eigenpair residual is not automatically a certified upper bound on the largest eigenvalue.

The unified oversampling conditions use the version 7 constants \(\rho_{\rm st}\) and \(C_{\rm rb}\):
\[
\frac{C_a^2}{c_W}\left(\frac{C_{\rm sd}}{c_W}\Theta_\ell\right)^2
\leq\rho_{\rm st}\frac{\gamma_h(k)}{C_FC_I^2},\qquad
\frac{C_{\rm sd}}{c_W}\Theta_\ell\leq C_{\rm rb}\gamma_h(k),
\quad 0<\rho_{\rm st}<1,\ C_{\rm rb}>0.
\]
See `eq:stability-margin`. A single threshold \(\tau_\Theta(k)\) can enforce both conditions if the required constants are controlled. Distinguish a prescribed or calibrated threshold from a measured \(\Theta_\ell\); a chosen \(c/k\) threshold does not itself verify the hypotheses.

The positive solution scale and reference-load balance ratio are
\[
d_\mu=\max\{\|U_\mu\|_k,d_{\min}\},\qquad
q_\ell=d_{\mu_{\rm ref}}\Theta_\ell/\eta_{H,\mu_{\rm ref}}.
\]
Hold chosen scales fixed within a marking/compression decision. Evaluate the base scale and residual at \(U_{H,\ell,\mu}\); the accepted coupled or compressed solution may also be used under `rem:enriched-localization-perturbation`'s assumptions. For the accepted coupled solution, `rem:enriched-localization-perturbation` uses \(C_{\rm enr,\ell}=\lambda_W[C_I+C_aL_\ell(1+C_I\delta_\ell)/m_\ell]\) and assumes it remains bounded. This balance ratio motivates oversampling; it is not a directly evaluated true-error ratio. Zero-denominator behavior is an implementation policy to inspect.

## Regional enrichment, kernel lifting and POD

The region \(D\subseteq\Omega\) is prescribed and fixed during the adaptive run. It restricts unprojected corrections, while coarse and fine marking still cover all of \(\Omega\). Define
\[
\mathcal Z_D=\{z:N^2(z)\subseteq\overline D\},\quad
\mathscr M_{h,D}G=\sum_{z\in\mathcal Z_D}\xi_z(G),\quad
\eta_{D,\mu}(U)^2=\sum_{z\in\mathcal Z_D}\|\xi_z(R_{\mu,U})\|_k^2.
\]
For \(Z\subseteq V_h\), \(Z^\perp=\{v:b(v,z)=0\text{ for all }z\in Z\}\), and \(P_Z\) is the \(b\)-orthogonal projection. The proof's Ritz projections \(P,P^\dagger\) are different, defined with \(a\) on a kernel subspace.

Raw corrections \(\psi=\mathscr M_{h,D}R_{\mu,U}\) lie in \(W_{h,D}\). With \(S=X_{H,\ell}\oplus\mathcal E\), projected Schwarz directions are generated through
\[
g=(I-P_S)\mathscr M_{h,D}R_{\mu,U},\quad
K_S=(I-P_S)\mathscr M_{h,D}\mathscr A_h,\quad
\operatorname{span}\{g,K_Sg,\ldots,K_S^{s-1}g\}.
\]
\(P_S\) is energy-orthogonal projection; \((\mathscr A_hv)(w)=a(v,w)\). This operator is not merely a stiffness matrix. Projection may change support and kernel membership of final modes in \(\mathcal E\subseteq X_{H,\ell}^{\perp}\). The block adds at most \(s\) directions (`eq:enrichment-block`); this block length is distinct from full or retained enrichment rank.

The kernel-lifting operator and enriched spaces are
\[
\mathscr J_\ell=I-T_\ell I_H,\quad\widehat{\mathcal E}=\mathscr J_\ell\mathcal E\subseteq W_h,
\]
\[
X_{H,\ell}\oplus\mathcal E=X_{H,\ell}\oplus\widehat{\mathcal E},\qquad
Y_{H,\ell}\oplus\widehat{\mathcal E},
\]
\[
a(U_\mu(\mathcal E),y)=F_\mu(y)\quad(y\in Y_{H,\ell}\oplus\widehat{\mathcal E}).
\]
\(\mathscr J_\ell\) projects onto \(W_h\) along \(X_{H,\ell}\) and is injective on \(\mathcal E\). For the paper's regional construction, \(X_{H,\ell}+\mathcal E\subseteq X_{H,\ell}+W_{h,D}\), so \(\widehat{\mathcal E}\subseteq W_{h,D}\). This also holds for compressed subspaces; locally supported lifted tests do not require locally supported projected trial modes. Do not apply archived adjoint-test triangular identities to these tests. The proof's block diagonalization uses an additionally transformed basis, not the original assembled coupled blocks.

Targets and greedy selection use
\[
t_\mu=\max\{t_{\min},\tau_D\eta_{D,\mu}(U_{H,\ell,\mu})\},\qquad
\mu_D\in\arg\max_{\mu\in\mathcal P_{\rm tr}}
\frac{\eta_{D,\mu}(U_\mu(\mathcal E))}{t_\mu}.
\]
The intended stopping condition is that this maximum is at most one. Describe any actual budget/dependence/low-gain stop separately.

POD compresses \(\mathcal E_{\rm full}\) using the POD solution-difference snapshots
\[
w_\mu=(U_\mu^{\rm full}-U_{H,\ell,\mu})/d_\mu,\qquad\mathcal E_r\subseteq\mathcal E_{\rm full}.
\]
These functions are distinct from disk audit snapshots and need not lie entirely in \(\mathcal E_{\rm full}\). The average POD optimum contains both the irreducible outside-space component and the singular-value tail (`eq:pod-optimality`). Use \(d_\mu\) evaluated at the full coupled solutions for these snapshots. Define \(U_{\mu,r}:=U_\mu(\mathcal E_r)\), recomputed with tests \(Y_{H,\ell}\oplus\mathscr J_\ell\mathcal E_r\). Retain the smallest tested rank satisfying, for every training case,
\[
\|U_\mu^{\rm full}-U_{\mu,r}\|_k\leq\varepsilon_{\rm dist}t_\mu,\qquad
\eta_{D,\mu}(U_{\mu,r})\leq\max\{c_t t_\mu,c_\rho\eta_{D,\mu}(U_\mu^{\rm full})\}.
\]
\(r\) is the enrichment dimension; \(N_{\rm on}=\dim V_H+r\) counts unconstrained coarse DOFs plus enrichment. It is not the fine-grid node count, offline work or memory use.

## Stability and best approximation

Version 7 places kernel-lifted enrichment stability in `sec:enrichment-stability`. Define
\[
\lambda_W=C_a/c_W,\quad L_\ell=C_F+\lambda_W\delta_\ell,\quad
m_\ell=\widehat\alpha_{H,\infty}-(C_a^2/c_W)\delta_\ell^2.
\]
For finite-dimensional \(\mathcal E\subseteq X_{H,\ell}^{\perp}\), under fine-grid stability, kernel coercivity and \(m_\ell>0\), the coupled inf-sup bound
\[
\alpha_{\mathcal E}\geq\frac{\min\{m_\ell,c_W\}}{1+L_\ell^2}
\]
holds independently of enrichment dimension. It applies to every compressed subspace, so compression does not introduce a separate rank-dependent oversampling condition under these hypotheses.

The sharper quasi-optimality statement is
\[
\|u_{h,\mu}-U_\mu(\mathcal E)\|_k
\leq C_{\rm qo,\ell}\inf_{v\in X_{H,\ell}\oplus\mathcal E}\|u_{h,\mu}-v\|_k,
\qquad
C_{\rm qo,\ell}=\lambda_W(\lambda_W+C_I\delta_\ell)
\left(1+\frac{C_aL_\ell\delta_\ell}{m_\ell}\right).
\]
See `prop:augmented-quasioptimality` and `cor:kernel-enrichment-robust-quasioptimality`. The ALOD–Galerkin gap on the left is the actual solution error relative to the fine-grid Galerkin solution. The infimum is the best-approximation error in the trial space, which is at most that gap; the quantities are not equal in general. This bound does not remove fine-grid discretization error or prove observed convergence rates.

## Marking, fine-grid control and errors

For coarse elements,
\[
m_{T,\mu}=\sum_{z\in\mathcal N_H\cap\overline T}\eta_{H,z,\mu}^2/m_z,\qquad
\overline m_T=M^{-1}\sum_{\mu\in\mathcal P_{\rm tr}}m_{T,\mu}/d_\mu^2.
\]
\(m_z\) counts adjacent coarse elements. Dörfler marking captures fraction \(\theta_H\) of the average and supplements it for the worst case in \(\eta_{H,\mu}/d_\mu\); it does not impose a bulk condition on every case.

In version 7, fine marking uses \(\zeta_{K,\mu}^2/d_\mu^2\) for \(K\in\mathcal T_h\). Components are the volume residual \(h_K^2\|f_\mu+k^2U_\mu\|^2\), interior normal-derivative jumps and Neumann/impedance residuals. An auxiliary work mesh \(\mathcal T_{\rm w}\), when used by code, is an implementation object to identify separately. This guides refinement; the paper does not establish it as an unconditional certified exact-error bound for ALOD. It is distinct from the coarse Riesz indicator even if both code structs have a field `eta`.

Use
\[
\|u_\mu-U_\mu\|_k\leq
\underbrace{\|u_\mu-u_{h,\mu}\|_k}_{\text{fine-grid discretization error}}+
\underbrace{\|u_{h,\mu}-U_\mu\|_k}_{\text{error relative to the fine-grid reference}}.
\]
In the wavenumber study, \(E=E_{\mu_0}=\|u_k-U_k\|_k/\|u_k\|_k\) depends on \(k\), even though \(\mu_0\) is fixed; the table reports \(100E\). The corner study uses \(E_{\mu_{\rm ref}}\).

The paper's relative energy error is \(E_\mu=\|u_\mu-U_\mu\|_k/\|u_\mu\|_k\), with a denominator different from \(d_\mu\). A fit \(E\approx C N_{\rm on}^{-p}\) is observed decay versus online dimension (online compression), not a convergence theorem in total unknowns or a CPU complexity claim. Fine-only refinement may improve error without increasing \(N_{\rm on}\).

## Adaptive state notation

For explaining saved states, introduce this explanatory tuple, explicitly not a named paper definition:
\[
S_n=(\mathcal T_{H,n},\mathcal T_{h,n},I_{H,n},\ell_n,
X_{H_n,\ell_n},Y_{H_n,\ell_n},\mathcal E_{r_n},\{U_{n,\mu}\}_{\mu\in\mathcal P_{n,\rm computed}}).
\]
Some components are determined by others plus fixed operator/configuration data; they need not all be separately serialized. Restart additionally needs cursor, phase, marks, inheritance data, localization warm starts, configuration and journal consistency information.

An accepted state is a completed solver transaction at one adaptive step. A frozen audit preserves its approximation spaces while evaluating an audit parameter set. A new test case can be solved within those spaces without retraining them. This is freezing the accepted state: while adaptation advances to \(S_{n+1}\), that audit still refers to \(S_n\).
