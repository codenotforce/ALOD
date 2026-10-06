# ALOD

ALOD is a C++20 implementation of adaptive two-sided localized orthogonal
decomposition for the E1/E2 Helmholtz problems, with AFEM/UFEM/SLOD baselines,
regional enrichment, family marking, numerical audits and checkpoint recovery.
E3 uses the E1 family at five wavenumbers.

## Guides

| Guide | Purpose |
|---|---|
| [Architecture](docs/architecture.md) | Module ownership, entry points, data flow and dependencies |
| [Algorithm](docs/algorithm.md) | Current mathematics, adaptive state, estimators and paper/archive differences |
| [Configuration](docs/configuration.md) | Configuration fields, presets, operational switches and defaults |
| [Runtime](docs/runtime.md) | Concurrent audits, thread budgets, checkpoint publication, retention and recovery |
| [Experiments](docs/experiments.md) | Current E1/E2/E3 protocol, frozen results and historical comparisons |
| [Scientific workflow](docs/scientific_workflow.md) | Diagnostics, error metrics, figures, tables and portable delivery |
| [Validation](docs/validation.md) | Tests, provenance rules, benchmarks and evidence limits |

These seven guides and this README are the eight maintained user documents.
Optimization mechanisms are summarized where they affect mathematics or use;
the former per-optimization and dated review pages have been consolidated.

## Build and test

Linux or WSL requires GCC 11+, CMake 3.20+, Eigen 3.3+, OpenMP,
SuiteSparse/UMFPACK, BLAS and Python 3.10+. Ubuntu 22.04 was tested;
native Windows numerical builds have not been validated.

```sh
sudo apt-get install cmake g++ libeigen3-dev libsuitesparse-dev python3 pkg-config
cmake --preset release
cmake --build --preset release
ctest --preset release
```

The release preset uses two compiler jobs and portable architecture flags.
UMFPACK is required; CHOLMOD is optional for the reference-solver benchmark.
Core run tools require no Python third-party packages. Paper rendering uses
[requirements-paper.txt](requirements-paper.txt).

## Run

```sh
python3 tools/run_baseline.py --config configs/smoke/e1_afem.json --output results/e1-baseline
python3 tools/run_fixed.py --config configs/p2/e2_fixed.json --output results/e2-fixed
python3 tools/run_adaptive.py --config configs/adaptive/e1_smoke.json --output results/e1-adaptive
python3 tools/run_adaptive.py --config configs/adaptive/e1_smoke.json --output results/e1-adaptive --resume auto
python3 tools/run_audit.py --checkpoint results/e1-adaptive/checkpoints/latest --output results/e1-audit
```

Fresh output directories must be new; resume uses the original adaptive run.
`threads=0` inherits the OpenMP runtime; a positive count is explicit.
There is no hardware-based automatic thread-sizing or CPU-reserve mode.
See the runtime guide before choosing checkpoint cleanup or audit concurrency.

## Instructions

### Experiment definitions

E1 and E2 are manufactured Helmholtz problems,
$-\Delta u_\mu-k^2u_\mu=f_\mu$, with sources and boundary data derived from
the specified exact solution $u_\mu$. The parameter $\mu$ identifies a
right-hand-side case. Keep this exact solution distinct from the fine-grid
Galerkin solution and the ALOD reconstruction.

**E1: smooth localized oscillations.** On $\Omega=(0,1)^2$, use

$$
u_\mu(x,y)=x^2(1-x)^2\sin(\pi y)
\exp\!\left(-80\bigl[(x-z_{x,\mu})^2+(y-z_{y,\mu})^2\bigr]\right)
e^{\mathrm i kx}.
$$

The packet centre $z_\mu$ varies near the nominal centre $(3/4,1/2)$.
The top/bottom edges are Dirichlet, the left edge is Neumann and the right
edge is impedance. The main family preset uses $k=16$ without regional
enrichment. Use [e1_main.json](configs/adaptive/e1_main.json) and the
[E1 RHS table](data/rhs/rhs_e1.json); `wave_x` and `wave_y` specify the centre.

**E2: corner singularity with localized oscillations.** On the L-shaped domain
$\Omega=(-1,1)^2\setminus([0,1]\times[-1,0])$, let
$w(x,y)=(1-x^2)^2(1-y^2)^2$ and define

$$
\begin{aligned}
u_\mu(x,y)&=c_\mu w(x,y)r^{2/3}\sin(2\vartheta/3)\\
&\quad+\alpha_\mu p_{z_\mu}(x,y)
e^{\mathrm i[k(x-z_{x,\mu})+\varphi_\mu]},\\
p_z(x,y)&=\frac{xy\,w(x,y)}{z_xz_y\,w(z_x,z_y)}
e^{-80[(x-z_x)^2+(y-z_y)^2]}.
\end{aligned}
$$

Here $(r,\vartheta)$ are polar coordinates about the reentrant corner,
with $0\leq\vartheta\leq3\pi/2$. RHS cases vary the singular coefficient
$c_\mu$, packet centre $z_\mu$, amplitude $\alpha_\mu$ and phase
$\varphi_\mu$. The nominal case has $(c,\alpha,\varphi)=(1,0.5,0)$
and $z=(-0.5,0.5)$. The reentrant edges are Dirichlet and the outer edges
are impedance. The main preset uses $k=16$ and regional enrichment in
$D=\Omega\cap B_{0.6}(0)$. Use
[e2_main.json](configs/adaptive/e2_main.json) and the
[E2 RHS table](data/rhs/rhs_e2.json).

**E3: wavenumber scaling of the E1 nominal problem.** Solve the E1 nominal
RHS case (ID 0, centre $(3/4,1/2)$) separately at
$k\in\{8,16,32,64,128\}$, using the `e3_k*.json` presets in
[configs/adaptive](configs/adaptive/). Each wavenumber defines its own
Helmholtz operator. These presets target nominal relative exact energy error
$\|u_0-U_0\|_k/\|u_0\|_k\leq0.01$, with $U_0$ the nominal ALOD
reconstruction and
$\|v\|_k^2=\|\nabla v\|_{L^2(\Omega)}^2+k^2\|v\|_{L^2(\Omega)}^2$.
Initial meshes and stopping limits are set per preset.

For E1/E2 family runs, use `training_ids` to select RHS cases that drive marking
and enrichment; use `member_ids` to select cases included in the run and audit.
Consult [configuration](docs/configuration.md) for overrides and
[experiments](docs/experiments.md) for the campaign launcher.

### Project conventions

Use the maintained guides for current configuration and workflow instructions.
Consult [NOTICE](NOTICE.md) for source ownership and dependencies.

[AGENTS](AGENTS.md) and the [terminology skill](.agents/skills/alod-terminology/SKILL.md)
govern agent work. Its [mathematical notation](.agents/skills/alod-terminology/references/mathematical-notation.md)
and [implementation correspondence](.agents/skills/alod-terminology/references/implementation-correspondence.md)
are explanatory references. 