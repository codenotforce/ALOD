# Algorithm and mathematical contracts

This guide follows the notation in the repository
[mathematical reference](../.agents/skills/alod-terminology/references/mathematical-notation.md)
and describes current implementation policies and archive compatibility.
Read [configuration](configuration.md) for selectable policies
and [experiments](experiments.md) for frozen run scope.

## Problem, spaces and solutions

For RHS parameter \(\mu\), the Helmholtz weak problem is
\[
a(u_\mu,v)=F_\mu(v),\qquad
a(v,w)=(\nabla v,\nabla w)-k^2(v,w)-ik(v,w)_{\Gamma_R},
\]
with homogeneous Dirichlet test constraints and prescribed mixed boundary data.
The energy form is \(b(v,w)=(\nabla v,\nabla w)+k^2(v,w)\), with
\(\|v\|_k^2=b(v,v)=\|\nabla v\|^2+k^2\|v\|^2\).
Version 7 presents A reentrant corner (repository E2) and Wavenumber
dependence (E3, using the E1 manufactured problem). E1 remains a supported
localized packet problem. E3 fixes its center and changes k, with separate
operators rather than one fixed-operator load family.

Keep the manufactured exact solution \(u_\mu\), fine-grid Galerkin solution
\(u_{h,\mu}\), and ALOD reconstruction \(U_\mu(\mathcal E_r)\) distinct.
The fine-grid solution is not the exact solution. RHS cases contain parameterized
problem definitions; `sample` is their stable identifier.

On nested meshes \(\mathcal T_H,\mathcal T_h\), let
\(W_h=\ker(I_H|_{V_h})\), with localized primal/adjoint correctors
\(\mathcal C_\ell,\mathcal C_\ell^*\). The un-enriched trial/test spaces are
\(X_{H,\ell}=T_\ell V_H\), \(Y_{H,\ell}=T_\ell^*V_H\), where
\(T_\ell=I-\mathcal C_\ell\) and \(T_\ell^*=I-\mathcal C_\ell^*\).
The base reconstruction is \(U_{H,\ell,\mu}=T_\ell u_{H,\ell,\mu}\).
The PG system tests the actual reconstructed solution against the test space.
In coordinates its base block is \(B_{\rm test}^* A_h B_{\rm trial}\).
Sparse factors are mutable; calls sharing one solver context are serialized.

## Interpolation, patches and archive differences

| Object | Manuscript/current explicit mode | Archive compatibility |
|---|---|---|
| Quasi-interpolation | Element L2 projection, incident-area vertex averaging, zero Dirichlet rows | Equal incident-element averaging |
| Riesz patches | Geometric \(N^2(z)\), including geometric Dirichlet vertices | Interpolation-support-distance expansion; tested E2 grids can produce \(N^3(z)\) |
| Enrichment tests | Kernel lifting | Earlier global adjoint-test policy |

Area and arithmetic averages agree on equal-area meshes and both reproduce
coarse P1 functions. The graded counterexample has operator difference
0.15023130314433289; reproduction alone does not distinguish the definitions.
The LOD baseline (`SLOD` in code, not the separate super-localized LOD method)
uses uniform coarse meshes and the retained archive interpolation;
new fixed/adaptive manuscript modes select area weighting explicitly.
Archive comparisons select arithmetic interpolation plus archived patches.
The new E2 cold-start deflated recurrence differs from the old reference-load
raw-kernel bootstrap. Historical curves are not relabelled as manuscript runs.

## Estimators and marking

For residual \(R_{\mu,U}=F_\mu-a(U,\cdot)\), the constrained local Riesz
representative \(\xi_z\in W_h(N^2(z))\) satisfies
\[
b(\xi_z,w)=R_{\mu,U}(w),\qquad
\eta_{H,\mu}(U)^2=\sum_z\|\xi_z\|_k^2.
\]
This is a sum of local squared norms, not the energy norm of their overlapping
sum. Per-element allocation conserves the local contributions.
Regional selection uses complete patch support inside the prescribed region;
an empty selection reports `empty_region` and does not expand the region.

Family marking freezes \(d_\mu=\max(\|U_\mu\|_k,10^{-12})\), uses mean
normalized squared mass, and supplements for the current worst normalized
training case. Only `training_ids` drive marking and dictionary growth.
Test/shift RHS cases and exact-solution norms never enter those decisions.
Reference refinement uses the distinct conforming strong residual on the
reference mesh, including body, jump, Neumann and Robin terms.

## Localization

The paper defines
\[
G_\ell(v_H)(w)=\overline{a(w,T_\ell^*v_H)},\qquad
\Theta_\ell^2=\lambda_{\max}(R_\ell,B_H),
\]
where \((B_H)_{ij}=b(\phi_j,\phi_i)\) and
\((R_\ell)_{ij}=\sum_z b(\xi_z(G_\ell(\phi_j)),\xi_z(G_\ell(\phi_i)))\).
The denominator uses the uncorrected coarse energy. In implementation
coordinates, with defect matrix \(D=A_h^*B_{\rm test}\), constrained Riesz
action \(R_{\rm AS}\) and coarse energy matrix \(B_H\), this is
\(D^*R_{\rm AS}Dx=\lambda B_Hx\).
For the real interpolation and complex-symmetric operator used here, the
conjugate primal defect has the same spectrum; the two energies are not added.
The corrector defect is \(\mathcal D_\ell=\mathcal C_\infty-\mathcal C_\ell\),
with energy operator norm \(\delta_\ell\). Under the paper's assumptions,
\(\Theta_\ell/(C_aC_{\rm syn})\leq\delta_\ell\leq C_{\rm sd}\Theta_\ell/c_W\).
Dense checks and matrix-free Ritz residual gates verify numerical contracts;
convergence does not provide a rigorous upper bound on the extreme eigenvalue.

Current production promotes on either
\[
\Theta_\ell>3.2/k,\qquad
q_\ell=d_{\mu_{\rm ref}}\Theta_\ell/\eta_{H,\mu_{\rm ref}}>0.3,\qquad
d_\mu=\max\{\|U_\mu\|_k,10^{-12}\}.
\]
The absolute scalar is manually stored; runtime does not infer a scaling law.
Lazy checks are initial, terminal, explicit extra states, and after about
\(2^\ell\) accepted states since the last check. Reference sweeps count.
Every check resets the counter. Promotion rebuilds/rechecks the same mesh pair
until acceptance or `maximum_ell`; cap exhaustion is reported explicitly.
Missing/deferred Theta is null, not zero. `fixed` forbids promotion.
Historical raw thresholds are preserved only by their explicit configurations.

## Regional enrichment and compression

The region \(D\) is prescribed and fixed; only complete patches inside it
contribute unprojected corrections. Global marking is not restricted to \(D\).
Raw regional Schwarz corrections are kernel representatives. Their energy
projection off \(X_{H,\ell}\oplus\mathcal E\) can give final modes global support.
Let \(P_Z\) be the \(b\)-orthogonal projection and \(Z^\perp\) its energy
orthogonal complement; these replace the older superscript-\(k\) notation.
Current kernel lifting uses
\[
\mathscr J_\ell=I-T_\ell I_H,\quad
X(\mathcal E)=X_{H,\ell}\oplus\mathcal E,\quad
Y(\mathcal E)=Y_{H,\ell}\oplus\mathscr J_\ell\mathcal E.
\]
For this regional construction, \(\mathscr J_\ell\mathcal E\subseteq W_{h,D}\),
including compressed subspaces, even when projected trial modes are global.
Solve the full coupled PG system; archived adjoint-test triangular identities
must not be applied to these tests.

Training fixes each target to `max(1e-12, 0.1*eta_D(base))`, selects the largest
training budget ratio and grows by at most two directions per round. Budget
success, two low-gain rounds, rank exhaustion, dependence and empty support are
distinct stops; a rank cap is not proof that every target was reached.
The second seed follows the deflated normalized current direction.

Energy POD uses normalized full-minus-base solution-difference snapshots,
distinct from disk audit snapshots. Candidate ranks solve the actual PG system;
every training distance must be at most `0.2*target`, and residual at most
`max(1.1*target, 1.03*full_residual)`. If no rank qualifies, retain the full
dictionary. `training_met` and `compression_accepted` are independent.
Ell changes reproject retained raw representatives; mesh changes inject and
repair them in the new kernel. Transfer does not guarantee exact old-span
preservation. `inherit=false` selects reset.

The two-sided residual estimate controls \(\|u_{h,\mu}-U_\mu\|_k\),
provided the base-test equations hold and kernel coercivity and fine-grid
stability are satisfied. Kernel-lifted coupling has rank-independent
stability when \(m_\ell>0\), as specified in the
[mathematical reference](../.agents/skills/alod-terminology/references/mathematical-notation.md).
These estimates do not remove fine-grid discretization error.

## Adaptive state and baseline scope

The following is the code counting convention. Paper Algorithm 2 instead
performs one coarse update and then \(m_{\rm ref}\) fine-only updates, solving
after each; do not infer identical state counts from those two descriptions.

State zero is the initial accepted solve. A coarse cycle performs one joint
coarse/reference update followed by `m_ref-1` reference sweeps, so an unbounded
complete horizon contains `1 + cycles*m_ref` accepted states. State/cycle/sweep,
mesh revision and ell event IDs are separate. Partial-cycle caps retain the
next phase. E1 reference-only closure preserves the proposed coarse-refinement
closure while leaving the accepted coarse mesh fixed.

AFEM uses strong-residual marking and UMFPACK; UFEM marks all elements.
LOD (`SLOD` in code) rebuilds a uniform reference mesh with `reference_gap`
(default 4), ell=3,
direct Schur correctors and reduced SparseLU. Its current patch loop is parallel.
Baseline exact-error stopping is separate from residual marking. Standalone
baseline runs do not implement adaptive checkpoint continuation.
