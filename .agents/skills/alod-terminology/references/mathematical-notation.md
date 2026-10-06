# Version 4 mathematical notation

Reviewed on 2026-10-05 against the complete `../LOD_paper/ALOD/version4/arxiv.tex`, SHA-256 `6cd4194e28e0cd58ff0d0ff9e6c588e96377502b69e70fa09717f71a9b9026d0`. This reference records definitions, not a frozen guarantee about future manuscript or code revisions. Use manuscript equation labels to locate arguments rather than unstable line numbers.

Contents: variational problem; meshes and LOD; residual and localization; enrichment and POD; refinement and errors; adaptive states.

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
b_k(v,w)=(\nabla v,\nabla w)_\Omega+k^2(v,w)_\Omega,
\quad\|v\|_k^2=b_k(v,v).
\]
The paper's \(V^\star\) is the continuous anti-dual: \(F_\mu\) is conjugate-linear in its test argument. A volume source gives
\[
F_\mu(v)=(f_\mu,v)_\Omega,\qquad -\Delta u_\mu-k^2u_\mu=f_\mu.
\]
For nodal basis \(\{\lambda_i\}\), \((\mathbf b_\mu)_i=F_\mu(\lambda_i)\). The vector changes with the mesh/basis even if \(\mu\) does not. Do not confuse \(\mathbf b_\mu\) with the energy form \(b_k\).

A right-hand-side case is the parameterized problem indexed by \(\mu\). A family is \(\{F_\mu:\mu\in\mathcal P\}\); its finite training parameter set is \(\mathcal P_{\rm tr}=\{\mu_1,\ldots,\mu_M\}\). Training constructs shared approximation directions and/or drives family marking; it does not imply random sampling or a neural model. A mathematical parameter need not equal its numeric table ID.

At fixed \(k\), coefficients, boundary conditions and mesh pair, E1/E2 share the operator and vary the RHS. Across changed \(k\) or meshes this fixed-operator statement no longer applies.

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

\(N^m(S)\) denotes an \(m\)-layer coarse patch with adjacency through shared vertices. Correctors use element patches \(N^\ell(T)\); residual Riesz problems use vertex patches \(N^2(z)\). \(\ell\) is neither a mesh level nor a physical radius.

The ideal correctors satisfy
\[
a(\mathcal Q_\infty v,w)=a(v,w),\qquad
a(w,\mathcal Q_\infty^*v)=a(w,v)\qquad(w\in W_h).
\]
Localizing their element problems and summing them gives \(\mathcal Q_\ell,\mathcal Q_\ell^*\). Set
\[
T_\ell=I-\mathcal Q_\ell,\quad T_\ell^*=I-\mathcal Q_\ell^*,\quad
X_{H,\ell}=T_\ell V_H,\quad Y_{H,\ell}=T_\ell^*V_H.
\]
The paper's star has this variational meaning; it is not automatically the conjugate transpose of a nodal matrix representing \(T_\ell\).

The base PG solution \(U_{H,\ell,\mu}=T_\ell u_{H,\ell,\mu}\) satisfies
\[
a(U_{H,\ell,\mu},T_\ell^*v_H)=F_\mu(T_\ell^*v_H)\quad(v_H\in V_H).
\]
\(u_{H,\ell,\mu}\) is the coarse coefficient function, \(U_{H,\ell,\mu}\) its corrected reconstruction. For enrichment use \(X^0=X_{H,\ell}\), \(Y^0=Y_{H,\ell}\), \(U_\mu^0=U_{H,\ell,\mu}\).

## Residual indicators and localization

The residual functional is \(R_{\mu,U}=F_\mu-a(U,\cdot)\). On the local kernel \(W_{h,z}\) supported in the closed patch \(N^2(z)\), define
\[
b_k(\xi_z(G),w)=G(w)\quad(w\in W_{h,z}),\qquad
\eta(G)^2=\sum_z\|\xi_z(G)\|_k^2.
\]
Hence \(\eta_{H,z,\mu}(U)=\|\xi_z(R_{\mu,U})\|_k\) and \(\eta_{H,\mu}(U)=\eta(R_{\mu,U})\). These involve local energy solves, not Euclidean algebraic residual norms. The additive Schwarz correction is \(\mathscr M_hG=\sum_z\xi_z(G)\), with \(G(\mathscr M_hG)=\eta(G)^2\).

The adjoint coarse-space defect functional and localization indicator are
\[
G_\ell(v_H)(w)=\overline{a(w,T_\ell^*v_H)},\qquad
\Theta_\ell=\sup_{v_H\ne0}\frac{\eta(G_\ell(v_H))}{\|v_H\|_k}.
\]
Call \(\Theta_\ell\) localization indicator and the paper's \(\delta_\ell\) localization-error operator norm. They are related, not identical. \(\Theta_\ell^2\) is a largest generalized eigenvalue with uncorrected coarse energy in the denominator; it is independent of the RHS and retained enrichment on a fixed base state.

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

The positive solution scale and nominal balance ratio are
\[
d_\mu=\max\{\|U_\mu\|_k,d_{\min}\},\qquad
\widetilde q_\ell=d_{\mu_{\rm nom}}\Theta_\ell/\eta_{H,\mu_{\rm nom}}.
\]
Hold chosen scales fixed within a marking/compression decision. This balance ratio motivates oversampling; it is not a directly evaluated true-error ratio. Zero-denominator behavior is an implementation policy to inspect.

## Regional enrichment, kernel lifting and POD

For prescribed \(D\subseteq\Omega\), define
\[
\mathcal Z_D=\{z:N^2(z)\subseteq\overline D\},\quad
\mathscr M_{h,D}G=\sum_{z\in\mathcal Z_D}\xi_z(G),\quad
\eta_{D,\mu}(U)^2=\sum_{z\in\mathcal Z_D}\|\xi_z(R_{\mu,U})\|_k^2.
\]
Raw corrections \(\psi=\mathscr M_{h,D}R_{\mu,U}\) lie in \(W_{h,D}\). With \(S=X^0\oplus\mathcal E\), projected Schwarz directions are generated through
\[
g=(I-P_S^k)\mathscr M_{h,D}R_{\mu,U},\quad
K_S=(I-P_S^k)\mathscr M_{h,D}\mathscr A_h,\quad
\operatorname{span}\{g,K_Sg,\ldots,K_S^{s_{\rm AS}-1}g\}.
\]
\(P_S^k\) is energy-orthogonal projection; \((\mathscr A_hv)(w)=a(v,w)\). This operator is not merely a stiffness matrix. Projection may change support and kernel membership of final modes in \(\mathcal E\subseteq(X^0)^{\perp_k}\).

The kernel-lifting operator and enriched spaces are
\[
\mathscr J_\ell=I-T_\ell I_H,\quad\widehat{\mathcal E}=\mathscr J_\ell\mathcal E\subseteq W_h,
\]
\[
X(\mathcal E)=X^0\oplus\mathcal E=X^0\oplus\widehat{\mathcal E},\quad
Y(\mathcal E)=Y^0\oplus\widehat{\mathcal E},\quad
a(U_\mu(\mathcal E),y)=F_\mu(y)\quad(y\in Y(\mathcal E)).
\]
\(\mathscr J_\ell\) projects along \(X^0\). Do not apply archived adjoint-test triangular identities to these tests. The proof's block diagonalization uses an additionally transformed basis, not the original assembled coupled blocks.

Targets and greedy selection use
\[
t_\mu=\max\{t_{\min},\tau_D\eta_{D,\mu}(U_\mu^0)\},\qquad
\mu_D\in\arg\max_{\mu\in\mathcal P_{\rm tr}}
\frac{\eta_{D,\mu}(U_\mu(\mathcal E))}{t_\mu}.
\]
The intended stopping condition is that this maximum is at most one. Describe any actual budget/dependence/low-gain stop separately.

POD compresses \(\mathcal E_{\rm full}\) using the POD solution-difference snapshots
\[
w_\mu=(U_\mu^{\rm full}-U_\mu^0)/d_\mu,\qquad\mathcal E_r\subseteq\mathcal E_{\rm full}.
\]
These functions are distinct from disk audit snapshots and need not lie entirely in \(\mathcal E_{\rm full}\). The average POD optimum contains both the irreducible outside-space component and the singular-value tail (`eq:pod-optimality`). Retain the smallest tested rank satisfying, for every training case,
\[
\|U_\mu^{\rm full}-U_\mu(\mathcal E_r)\|_k\leq\varepsilon_{\rm dist}t_\mu,\qquad
\eta_{D,\mu}(U_\mu(\mathcal E_r))\leq\max\{c_t t_\mu,c_\rho\eta_{D,\mu}(U_\mu^{\rm full})\}.
\]
\(r\) is the enrichment dimension; \(N_{\rm on}=\dim V_H+r\) counts unconstrained coarse DOFs plus enrichment. It is not the fine-grid node count, offline work or memory use.

## Marking, fine-grid control and errors

For coarse elements,
\[
m_{T,\mu}=\sum_{z\in\mathcal N_H\cap\overline T}\eta_{H,z,\mu}^2/m_z,\qquad
\overline m_T=M^{-1}\sum_{\mu\in\mathcal P_{\rm tr}}m_{T,\mu}/d_\mu^2.
\]
\(m_z\) counts adjacent coarse elements. Dörfler marking captures fraction \(\theta_H\) of the average and supplements it for the worst case in \(\eta_{H,\mu}/d_\mu\); it does not impose a bulk condition on every case.

Fine marking uses \(\zeta_{K,\mu}^2/d_\mu^2\) on an auxiliary conforming refinement \(\mathcal T_{\rm w}\) of \(\mathcal T_h\). Components are the volume residual \(h_K^2\|f_\mu+k^2U_\mu\|^2\), interior normal-derivative jumps and Neumann/impedance residuals. This guides refinement; the paper does not establish it as an unconditional certified exact-error bound for ALOD. It is distinct from the coarse Riesz indicator even if both code structs have a field `eta`.

Use
\[
\|u_\mu-U_\mu\|_k\leq
\underbrace{\|u_\mu-u_{h,\mu}\|_k}_{\text{fine-grid discretization error}}+
\underbrace{\|u_{h,\mu}-U_\mu\|_k}_{\text{error relative to the fine-grid reference}}.
\]
The paper's relative energy error is \(E_\mu=\|u_\mu-U_\mu\|_k/\|u_\mu\|_k\), with a denominator different from \(d_\mu\). A fit \(E\approx C N_{\rm on}^{-p}\) is observed decay versus online dimension (online compression), not a convergence theorem in total unknowns or a CPU complexity claim. Fine-only refinement may improve error without increasing \(N_{\rm on}\).

## Adaptive state notation

For explaining saved states, introduce this explanatory tuple, explicitly not a named paper definition:
\[
S_n=(\mathcal T_{H,n},\mathcal T_{h,n},I_{H,n},\ell_n,
X_n^0,Y_n^0,\mathcal E_{r_n},\{U_{n,\mu}\}_{\mu\in\mathcal P_{n,\rm computed}}).
\]
Some components are determined by others plus fixed operator/configuration data; they need not all be separately serialized. Restart additionally needs cursor, phase, marks, inheritance data, localization warm starts, configuration and journal consistency information.

An accepted state is a completed solver transaction at one adaptive step. A frozen audit preserves its approximation spaces while evaluating an audit parameter set. A new test case can be solved within those spaces without retraining them. This is freezing the accepted state: while adaptation advances to \(S_{n+1}\), that audit still refers to \(S_n\).
