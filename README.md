# swesphere

Shallow-water equations on the sphere (spectral filter via `pyshtools`, RK4)
with two stationary regimes for data-assimilation testbeds:

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

The original explicit-Euler stepper is kept as `scheme="euler"` for reference;
it blows up in 25 days at `dt = 120 s`. See `integrators.py`, `forcing.py` and
`presets.py` for the measurements behind each choice.

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
docker compose run --rm all           # ~4 h on one core: stability, regimes + climatologies, footprints, DA baseline, figures
```

`experiments/exp01_stability.py` (Euler vs RK4 energy series), `exp02_regimes.py`
(200-state climatology per preset with statistics and decorrelation),
`exp03_footprint.py` (balanced vs unbalanced perturbations: advection vs gravity
waves), `exp04_da_baseline.py` (first analysis, h observed at 10%: LETKF vs
modified-Cholesky EnKF, 10 seeds), `figures/make_figures.py`. Results in `results/`.
