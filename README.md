# swesphere

Shallow-water equations on the sphere (spectral filter via `pyshtools`, RK4)
with three stationary regimes, meant as a testbed for data assimilation (the model only; the filters live in pyteda):

- `presets.waves()`, Rossby-Haurwitz, no forcing: zonal jets with weak waves.
- `presets.two_jets()`, Galewsky jets at ±45° with zonal relaxation: sustained
  barotropic turbulence, errors advected by the flow.

```python
from swesphere import presets, climatology
model, x0 = presets.two_jets()
S = climatology.build_climatology(model, x0, spinup_days=40, n_snapshots=200, every_days=2)
scales, per_field = climatology.anomaly_scales(model, S)
```

State vector `[u, v, h]` (row-major, north to south), `model.var_blocks`,
`model.interior_mask()` (outside the polar sponge), `model.propagate(x, [t0, t1])`.

Time stepping is RK4 (`dt = 120 s` by default); EXP-01 documents the range of
`dt` and filter cadence for which each preset stays stationary. See `forcing.py`
and `presets.py` for the measurements behind each choice.

## Presets

| preset | what it is | stationary statistics (rows outside the sponge) |
|---|---|---|
| `waves` | Rossby-Haurwitz IC, no forcing | eddy std of h ~20 m, `|U|` ~36 m/s; nearly linear waves |
| `one_jet` | Galewsky jet at 45°N + zonal relaxation (5 d) | eddy std ~110 m, `|U|max` ~70 m/s, regular (wavenumber ~6) |
| `two_jets` | Galewsky jets at ±45° + zonal relaxation | eddy std ~150-160 m, irregular turbulence in both hemispheres |

## The article's experiments (Docker)

```
docker compose build
QUICK=1 docker compose run --rm all   # smoke, ~5 min
docker compose run --rm all           # ~2.5 h on one core: stability, regimes + climatologies, footprints, figures
```

`experiments/exp01_stability.py` (stable envelope: dt and filter cadence, 100 days per preset), `exp02_regimes.py`
(200-state climatology per preset with statistics and decorrelation),
`exp03_footprint.py` (balanced vs unbalanced perturbations: advection vs gravity
waves), `figures/make_figures.py`. Data-assimilation examples on this model are in
`pyteda` (`pip install pyteda[swe]`). Results in `results/`.

## Using your own integrator

The package publishes the right-hand side and ships RK4 as the tested stepper.
Any other scheme can be plugged in and runs with the same filter, sponge and
forcing:

```python
from swesphere import presets
from swesphere.integrators import INTEGRATORS

def my_step(rhs, u, v, h, dt):          # rhs(u, v, h) -> (du, dv, dh)
    du, dv, dh = rhs(u, v, h)
    return u + dt * du, v + dt * dv, h + dt * dh

INTEGRATORS["mine"] = my_step
model, x0 = presets.two_jets(scheme="mine")
```

`swesphere.dynamics.rhs(u, v, h, grid)` is the physics alone; `model.rhs(u, v, h)`
adds the preset's forcing. Multistep schemes that need state from earlier steps
can keep it in the stepper's closure.
