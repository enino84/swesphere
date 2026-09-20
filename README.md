# swesphere

[![tests](https://github.com/enino84/swesphere/actions/workflows/tests.yml/badge.svg)](https://github.com/enino84/swesphere/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/swesphere.svg)](https://pypi.org/project/swesphere/)
[![Python](https://img.shields.io/pypi/pyversions/swesphere.svg)](https://pypi.org/project/swesphere/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/enino84/swesphere/blob/main/LICENSE)

Shallow-water equations on the sphere with three documented flow regimes: a small geophysical
flow for numerical experiments, verified against the standard test cases and reproducible from a
container. This package is the model only; assimilation filters live in `pyteda`.

![The two_jets preset on the sphere](https://raw.githubusercontent.com/enino84/swesphere/main/docs/img/hero_two_jets.png)

*A state of the `two_jets` record: zonal wind, meridional wind, and the departure of `h` from its
zonal mean drawn as relief (exaggerated). Polar caps under the sponge are grey. Rendered with
`python figures/make_globe3d.py two_jets`.*

## Install

```
pip install swesphere
```

Documentation and API: https://aml-cs.github.io/software/swesphere/

## Quick start

```python
from swesphere import presets, climatology
model, x0 = presets.two_jets()
S = climatology.build_climatology(model, x0, spinup_days=40, n_snapshots=200, every_days=2)
scales, per_field = climatology.anomaly_scales(model, S)
```

State vector `[u, v, h]` (row-major, north to south), `model.var_blocks`,
`model.interior_mask()` (outside the polar sponge), `model.propagate(x, [t0, t1])`.

## The model

Equations: vector-invariant momentum and flux-form continuity,

    du/dt = eta v - (1 / a cos(lat)) dB/dlon        eta = zeta + f,   B = g h + (u^2 + v^2) / 2
    dv/dt = -eta u - (1 / a) dB/dlat
    dh/dt = -div(h V)

Discretization: regular latitude-longitude (Driscoll-Healy) grid, `Nlat = 2 (LMAX + 1)`,
`Nlon = 4 (LMAX + 1)` (66 x 132 at the default `LMAX = 32`); exact zonal derivative by FFT,
second-order centered meridional difference; RK4 in time (`dt = 120 s`). Dissipation: an
exponential filter in spherical harmonics (`pyshtools`) applied every `filter_every = 3` steps,
and a multiplicative sponge on the wind over the `n_pole_rows = 10` rows next to each pole
(poleward of about 63 deg), which is how the model avoids the polar singularity of the grid.
The rows under the sponge are not meant to be used (`model.interior_mask()`).

![Verification](https://raw.githubusercontent.com/enino84/swesphere/main/docs/img/verification.png)

*Left: normalized l2 error of `h` in Williamson test case 2, with and without the polar sponge.
Right: relative vorticity of the Galewsky instability at day 6, LMAX 64 (compare with the
published reference solution).*

Two consequences that the experiments quantify:

- The scheme is second order (EXP-00: Williamson test case 2 without the sponge, l2 error of h
  2.2e-4 at LMAX 32 and 5.5e-5 at LMAX 64 after 5 days). With the sponge the steady state of
  that test is disturbed from the polar caps and the error in the interior is about 1e-2.
- The dissipation acts per step, not per second, so the eddy amplitude of a regime depends on
  `dt * filter_every`. The statistics quoted for the presets hold for the defaults (EXP-01).

## Presets

| preset | what it is | regime |
|---|---|---|
| `waves` | solid-body rotation (Williamson TC2) + wave packets (wavenumber 4 north, 3 south), no forcing | quasi-linear Rossby waves on a zonal flow, **slowly decaying** (no energy source; the filter is the only sink) |
| `one_jet` | Galewsky jet at 45 N + relaxation of the zonal mean (5 d) | stationary, fairly regular eddies in the north |
| `two_jets` | Galewsky jets at +-45 deg + the same relaxation | stationary irregular turbulence in both hemispheres; errors advected by the flow |

![The three regimes](https://raw.githubusercontent.com/enino84/swesphere/main/docs/img/regimes.png)

The measured statistics of each regime (eddy std of h, wind, anomaly scales, lag correlations,
first-half against second-half amplitude) are written by EXP-02 to `results/EXP-02/regimes.csv`.
Measured with the defaults: eddy std of `h` 8.8 m (`waves`, decaying), 124 m (`one_jet`) and
156 m (`two_jets`); `|U|max` about 42, 76 and 76 m/s. Note the scales in the figure above: the
meridional wind and the depth anomaly of `waves` are an order of magnitude weaker than those of
the forced presets.

![Nine days of the two_jets regime](https://raw.githubusercontent.com/enino84/swesphere/main/docs/img/chaos_vorticity.png)

*Relative vorticity of the `two_jets` regime over nine days, starting from a spun-up state. The jet
meanders, sheds eddies and reforms: the flow is chaotic but statistically stationary, which is what
makes it usable as a reference.*

## The article's experiments (Docker)

```
docker compose build test
docker compose run --rm test              # unit tests
QUICK=1 docker compose run --rm all       # smoke
docker compose run --rm all               # everything, one core
```

or detached and in parallel (`verification`, `stability` and `regimes` are independent;
`footprint` needs `regimes`; `figures` needs all):

```
docker compose run -d verification; docker compose run -d stability; docker compose run -d regimes
docker compose run -d footprint           # when regimes has finished
docker compose run --rm figures           # when everything has finished
```

| experiment | what it measures | output |
|---|---|---|
| `exp00_verification.py` | Williamson TC2 error norms (LMAX 32, 64; with and without sponge), Galewsky instability: invariants and vorticity | `results/EXP-00/` |
| `exp01_stability.py` | the stable range of `dt` for the shipped RK4 (pushed until it blows up) and the eddy amplitude against the filter interval | `results/EXP-01/` |
| `exp02_regimes.py` | 200-state record per preset, per-snapshot series, statistics, stationarity, decorrelation | `results/EXP-02/` |
| `exp03_footprint.py` | balanced against unbalanced perturbations: advection against gravity waves | `results/EXP-03/` |
| `figures/make_figures.py` | figures and tables of the article, one file per panel | `results/figures/` |

Besides those, `figures/make_globe3d.py <preset> [state.npy]` renders the illustrative
three-dimensional views above (it spins the preset up itself if no state is given).

Progress of a detached run: `docker ps`, `docker logs -f <name>`; every script logs the preset
and the step, and rewrites its CSV after every case.

## Using your own integrator

The package publishes the right-hand side and ships RK4 as the tested stepper. Any other
scheme can be plugged in and runs with the same filter, sponge and forcing (its stable range
of `dt` is its own: rerun EXP-01):

```python
from swesphere import presets
from swesphere.integrators import INTEGRATORS

def my_step(rhs, u, v, h, dt):          # rhs(u, v, h) -> (du, dv, dh)
    du, dv, dh = rhs(u, v, h)
    return u + dt * du, v + dt * dv, h + dt * dh

INTEGRATORS["mine"] = my_step
model, x0 = presets.two_jets(scheme="mine")
```

`swesphere.dynamics.rhs(u, v, h, grid)` is the physics alone; `model.rhs(u, v, h)` adds the
preset's forcing. `swesphere.diagnostics` gives mass, energy, potential enstrophy and the
Williamson error norms.

Changes between versions are listed in `CHANGELOG.md` (0.2.0 fixes the continuity equation of 0.1.0).

## Citing, contributing and license

If you use swesphere in your work, please cite it: the `CITATION.cff` file gives the reference, and
GitHub shows it under "Cite this repository". Bug reports and pull requests are welcome; see
`CONTRIBUTING.md`. Questions: enino@uninorte.edu.co.

swesphere is released under the MIT license (`LICENSE`). Its dependencies are NumPy and SciPy
(BSD-3-Clause) and pyshtools (BSD-3-Clause), all compatible with it.

