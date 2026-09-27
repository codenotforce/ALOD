# Revised policies and the E1 wavenumber campaign

**Historical recipe — not the current production launcher.** This page records the September 12 four-wavenumber protocol. Use [the ordered production protocol](production_campaign_20260927.md) and `tools/make_production_campaign.py` for new E1/E2/E3 runs. The old `make_kappa_campaign.py` reads mutable current presets while retaining old metadata; it is not a frozen replay tool and must not be used to claim reproduction of this protocol. The policies and limits below describe the original design only.

Version 4 update (2026-09-22): see [the current implementation review](paper_v4_review.md). Numerical policies and performance statements below describe their original migration stage unless explicitly updated.

The original campaign proposed new manuscript-mode trajectories rather than reproductions of archived values. References below to main presets describe their September 12 versions, not their current contents.

## Mathematical policies

- `interpolation=area`: area-weighted coarse-vertex averages.
- `riesz_patches=n2`: geometric two-ring vertex patches.
- E2 starts from an empty dictionary, selects the worst regional budget ratio,
  and generates subsequent directions from the deflated, normalized direction.
  The positive energy Riesz operator is reused from the regional estimator.
- `configs/adaptive/e2_main.json` selects compressed-dictionary inheritance and
  lazy checks. The threshold remains `Theta / eta_H(coupled) > 0.1`.
  Checks occur initially, terminally, and after `2^ell` accepted states since
  the previous check. A due E2 check may promote and recheck on the same mesh
  until the threshold passes or ell reaches four. No balance gate is active.
- The old `every` presets remain explicitly named diagnostic controls.

## Runtime wavenumber

`wavenumber` is a positive finite number, defaulting to 16 for older presets.
It enters the manufactured solution and gradient, cancellation-free source,
Helmholtz and positive energy operators, LOD spaces, reference solves and
independent audits. Checkpoint mathematical identities reject a changed
wavenumber. No factor or dictionary is shared between wavenumber runs.

The historical kappa-16 checkpoint key is retained. Missing `wavenumber` in an
older manifest means 16; this does not undo the intentional policy differences.

## E1 controlled experiment

Four methods, ALOD/AFEM/UFEM/SLOD, at each of kappa 8, 16, 32 and 64: 16 runs.
The smooth E1 nominal wave packet retains its envelope, amplitude, centre
(0.75, 0.5), Gaussian coefficient 80 and boundary conditions. Only its phase
wavenumber and matching manufactured source change. ALOD uses only member 0,
with no enrichment, fixed ell=3, theta_H=0.15, two reference sweeps per coarse
cycle with theta_c=[0.3,0.2], and minimum local reserve two.

Initial meshes use the coarsest uniform NVB level whose *measured maximum
triangle diameter* satisfies H0 <= pi/kappa. The initial reference gap is four
NVB levels. `alod_geometry` measures:

| kappa | coarse level | reference level | kappa H0 | kappa h0 | H0/h0 | initial online DoF |
|---|---|---|---|---|---|---|
| 8 | 4 | 8 | 2.828427125 | 0.707106781 | 4 | 15 |
| 16 | 6 | 10 | 2.828427125 | 0.707106781 | 4 | 63 |
| 32 | 8 | 12 | 2.828427125 | 0.707106781 | 4 | 255 |
| 64 | 10 | 14 | 2.828427125 | 0.707106781 | 4 | 1023 |

This geometric condition alone is not a verification of the theoretical
coercivity/resolution condition. Initial and subsequent PG, corrector and
constraint residuals are retained.

ALOD: 25 coarse cycles, 51 accepted states including the initial state. AFEM:
at most 97 states. Uniform FEM extends at most through level 20; SLOD through
coarse level 15/reference level 19. Baselines stop at the first recorded exact
relative error <=0.01. ALOD retains its predeclared full horizon; exact errors
never control its marking or localization checks. All first-target comparisons
use recorded states without interpolation, at targets 0.02 and 0.01. A run
stopped by a resource ceiling is labelled incomplete or not reached.

AFEM starts from the same coarse mesh and uses strong-residual theta=0.15.
UFEM uniformly refines the same initial coarse mesh. SLOD has fixed ell=3 and
uniform coarse/reference meshes with gap four throughout. SLOD uses independent
parallel corrector solves; its uniform-mesh interpolation agrees with both
vertex-averaging conventions.

## Diagnostics and precision

Every ALOD accepted state is checkpointed. Independent audits integrate the
analytic solution and gradient to obtain e=||u_exact-U||, f=||u_exact-u_h||,
g=||u_h-U|| and the analytic energy norm. The primary normalized outputs are
E=e/norm, F=f/norm, G=g/norm. They are not interpolant errors. The accepted
journal also contains eta_H, H_max, h_max, online/reference DoF, residuals,
marking data, ell checks and phase timings. Together these give eta_H/norm and
eta_H/g without an arbitrary calibration factor. There is no extra candidate
or exact solve in the production marking loop.

For kappa 32 and 64, audit three accepted states on *identical H and h*:
initial, middle, and the first state E<=0.02 (otherwise last). Duplicate state
IDs are removed. Compare ell=3 with ell=4, retaining the energy distance to the
accepted solution. This is a localization sensitivity diagnostic, not an exact
localization-error identity. The same snapshots also receive triangle orders
increased by four and two residual-correction steps for the direct FEM/PG
solves, using their existing factors; the strict audited residual bound is
1e-10. No tighter iterative tolerance is claimed for a direct solver.

The main quadrature is unchanged at all kappas: regular/Gaussian/singular
orders 12/16/24, with the E1 recursive limit eight. Precision diagnostics do
not affect the production trajectory.

## Revised E1/E2 production presets and resources

The queue also contains one revised E1 family production trajectory (51 states)
and one E2 lazy/inheritance trajectory (33 states), retaining the frozen
16-member training family and 48-member independent audit table. Total: 18
main runs. E2 retains radius 0.6, regional training/compression budgets, rank
cap 24, theta_H=0.15, global strong-residual candidate marking theta_c=0.2 and
one reference sweep. E1 retains its lazy threshold 1.2 and reference sweeps
[0.3,0.2]. The exact JSON files are authoritative.

The historical scheduling design assigned 16 disjoint physical cores and one
BLAS thread per job, with at most three jobs. Current launchers do not apply that
affinity policy; use the current production protocol for thread/resource settings.
Admission reserves 80 GiB; a live guard terminates a job if available memory
falls below 64 GiB, its process-group RSS exceeds its reservation (112 GiB for
adaptive jobs; 48 GiB SLOD; 32 GiB FEM), or free disk drops below 35 GiB.
Accepted checkpoints remain available after a guard stop. These guards reduce
OOM risk; an abrupt allocator spike is not a formal memory guarantee.

Adaptive node cap: one million. Baseline node cap: two million. The conservative
hierarchy-product ceiling is 16 billion entries; the actual dense-workspace
ceiling remains 64 million entries. Raising the former avoids rejecting sparse
hierarchies merely because reference nodes times coarse nodes is large. The
entry parser uses a 64-bit value; this does not allocate a matrix of that size.

## Historical preparation sequence (not a current launch recipe)

The sequence below records the original design. Reproduction requires its frozen
source and configurations; executing the old generator against current presets
is not equivalent. For new experiments use `tools/make_production_campaign.py`.

1. Build and run CTest, including `variable_wavenumber`.
2. Generate configuration and measured geometry using
   `tools/make_kappa_campaign.py --output "$CAMPAIGN" --build "$BUILD"`.
3. Launch `tools/run_campaign.py --campaign "$CAMPAIGN" --build "$BUILD"`
   inside tmux. Keep the source archive and its SHA-256 beside the campaign.
4. Inspect `status.json`, `resources.jsonl`, `logs`, and each run's manifest.
   `paper_complete` remains specific to the historical coverage schema; it is
   not the completion test for the single-member kappa controls.
5. After completion, prepare error-versus-online-DoF panels for each kappa,
   E versus online-DoF/kappa^2, ALOD E/F/G curves, first-target tables, eta
   effectivities, and localization/precision diagnostic tables. Report online
   and reference sizes separately. No manuscript edits are part of this run.
