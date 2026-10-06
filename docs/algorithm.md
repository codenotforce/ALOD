# Algorithm and mathematical contracts

This guide describes current manuscript-mode mathematics and explicit archive
compatibility. Read [configuration](configuration.md) for selectable policies
and [experiments](experiments.md) for frozen run scope.

## Problem, spaces and solutions

For RHS parameter \(\mu\), the Helmholtz weak problem is
\[
a(u_\mu,v)=F_\mu(v),\qquad
a(v,w)=(\nabla v,\nabla w)-k^2(v,w)-ik(v,w)_{\Gamma_R},
\]
with homogeneous Dirichlet test constraints and prescribed mixed boundary data.
The positive energy norm is
\(\|v\|_k^2=\|\nabla v\|^2+k^2\|v\|^2\).
E1 is a localized oscillatory packet; E2 combines a weighted corner singularity
and a Gaussian wave on the L-shaped domain. E3 changes k in the E1 family.

Keep the manufactured exact solution \(u_\mu\), fine-grid Galerkin solution
\(u_{h,\mu}\), and ALOD reconstruction \(U_\mu(\mathcal E_r)\) distinct.
The fine-grid solution is not the exact solution. RHS cases contain parameterized
problem definitions; `sample` is their stable identifier.

On nested meshes \(\mathcal T_H,\mathcal T_h\), let
\(W_h=\ker(I_H|_{V_h})\), with localized primal/adjoint correctors
\(\mathcal Q_\ell,\mathcal Q_\ell^*\). The un-enriched trial/test spaces are
\(X^0=T_\ell V_H\), \(Y^0=T_\ell^*V_H\), where
\(T_\ell=I-\mathcal Q_\ell\) and \(T_\ell^*=I-\mathcal Q_\ell^*\).
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
Baseline SLOD uses uniform coarse meshes and the retained archive interpolation;
new fixed/adaptive manuscript modes select area weighting explicitly.
Archive comparisons select arithmetic interpolation plus archived patches.
The new E2 cold-start deflated recurrence differs from the old nominal
raw-kernel bootstrap. Historical curves are not relabelled as manuscript runs.

## Estimators and marking

For residual \(R_{\mu,U}=F_\mu-a(U,\cdot)\), the constrained local Riesz
representative \(\xi_z\in W_h(N^2(z))\) satisfies
\[
(\xi_z,w)_k=R_{\mu,U}(w),\qquad
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

With defect matrix \(D=A_h^*B_{\rm test}\), additive constrained Riesz action
\(R_{\rm AS}\), and uncorrected coarse energy matrix \(M_H\), solve
\[
D^*R_{\rm AS}D x=\lambda M_Hx,\qquad \Theta_\ell=\sqrt{\lambda_{\max}}.
\]
For real interpolation and the complex-symmetric operator, the conjugate
primal defect has the same spectrum. The two energies are never added.
Dense checks and matrix-free Ritz residual gates verify numerical contracts;
convergence does not provide a rigorous upper bound on the extreme eigenvalue.

Current production promotes on either
\[
\Theta_\ell>3.2/k,\qquad
\max(\|U_{\mu_{\rm nom}}\|_k,10^{-12})\Theta_\ell/\eta_{H,\mu_{\rm nom}}>0.3.
\]
The absolute scalar is manually stored; runtime does not infer a scaling law.
Lazy checks are initial, terminal, explicit extra states, and after about
\(2^\ell\) accepted states since the last check. Reference sweeps count.
Every check resets the counter. Promotion rebuilds/rechecks the same mesh pair
until acceptance or `maximum_ell`; cap exhaustion is reported explicitly.
Missing/deferred Theta is null, not zero. `fixed` forbids promotion.
Historical raw thresholds are preserved only by their explicit configurations.

## Regional enrichment and compression

Raw regional Schwarz corrections are kernel representatives. Their energy
projection off \(X^0\oplus\mathcal E\) can give final modes global support.
Current kernel lifting uses
\[
\mathscr J_\ell=I-T_\ell I_H,\quad
X(\mathcal E)=X^0\oplus\mathcal E,\quad
Y(\mathcal E)=Y^0\oplus\mathscr J_\ell\mathcal E.
\]
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

## Adaptive state and baseline scope

State zero is the initial accepted solve. A coarse cycle performs one joint
coarse/reference update followed by `m_ref-1` reference sweeps, so an unbounded
complete horizon contains `1 + cycles*m_ref` accepted states. State/cycle/sweep,
mesh revision and ell event IDs are separate. Partial-cycle caps retain the
next phase. E1 reference-only closure preserves the proposed coarse-refinement
closure while leaving the accepted coarse mesh fixed.

AFEM uses strong-residual marking and UMFPACK; UFEM marks all elements.
SLOD rebuilds a uniform reference mesh with `reference_gap` (default 4), ell=3,
direct Schur correctors and reduced SparseLU. Its current patch loop is parallel.
Baseline exact-error stopping is separate from residual marking. Standalone
baseline runs do not implement adaptive checkpoint continuation.
