# Changelog

## 0.4.0 (2026-10-13)

**`propagate` now gives the same answer however a forecast is split, and does
what the documentation says.** All four changes below are fixes to behaviour
that was wrong in 0.3.0 and in every earlier version.

### Fixed

- **The filter schedule restarted on every call.** The exponential filter
  fired every `filter_every` steps *of the current call*, counted from 1. At
  the defaults (`dt = 120 s`, `filter_every = 3`) a call shorter than 360 s
  therefore never filtered at all, so the same forecast interval gave a
  different answer depending on how it was chopped up: one 360 s call filtered
  once, three 120 s calls filtered never. Measured on `two_jets` at LMAX 8,
  one call and six calls over the same 12 minutes differed by 274 m in `h`;
  with `filter_every = 1`, where both paths filter identically, they agreed to
  2e-12, which is what identified the schedule as the cause. At the default
  LMAX 32 the split path did not merely differ, it blew up within the first
  model day.

  This mattered most for the use the package is built for. An external driver,
  such as a data assimilation loop, calls `propagate` once per forecast window,
  and the length of that window is the experiment's choice: a 6-hour cycling
  experiment and a 10-minute one were running different models.

  The filter now fires on **absolute model time**, every `filter_every` steps
  counted from `t = 0`, so any partition of an interval gives a bit-identical
  result. The schedule is derived from `T`, not stored on the model, so
  propagating the members of an ensemble through one model object stays safe.

- **A zero-length interval advanced the state.** `n_steps` was
  `max(1, round(span/dt))`, so `propagate(x, [0, 0])` took a step; it now
  returns the state unchanged.

- **An interval that was not a whole number of steps landed elsewhere in
  silence.** A 60 s request with `dt = 120 s` returned the state at 120 s. It
  now raises `ValueError` naming the two nearest valid intervals. Integrating
  backwards raises as well.

- **`just_final_state=False` raised `NotImplementedError`** while the API
  documentation (and Table 3 of the submitted manuscript) said it returned
  every step. It now returns the trajectory, shape `(n_steps + 1, dim)`, first
  row `x0`.

### Unchanged

The experiments of the repository and the records they produce are **bit
identical** to 0.3.0: they drive the model in single long calls starting at
`t = 0`, which is exactly the case the old schedule got right. This was
verified by running both versions side by side over 10 model days
(`max |old - new| = 0`) and by rebuilding a climatology record with both
(`max |old - new| = 0`). No figure, table or number in the manuscript changes.

`climatology.build_climatology` now passes the running time to each leg
instead of restarting at zero. For the shipped presets each leg is a whole
number of filter periods, so the records are identical; the change is so that
a user's own `every_days` cannot fall out of step.

- **`swesphere.__version__` was a hand-written string** and still said
  `"0.3.0"` while `pyproject.toml` said otherwise, so the package could report
  a version it was not. It is now read from the installed metadata and cannot
  drift.

### Added

- Six tests covering the above. Five of them fail on 0.3.0.


## 0.3.0 (2026-09-20)

- EXP-03 redesigned. The first version placed the bumps at 30 N, the equator and 30 S (none on the jets of
  `two_jets`, and geostrophic balance is not defined at the equator) and tracked the centroid of |dh| over the
  whole sphere without periodic wrapping, which the radiated gravity waves dominate after a few hours. Now: sites
  on the jet cores and in the quiet subtropics, centroid of dh^2 above half maximum with a circular mean in
  longitude, hourly parcel trajectories, the steering wind averaged over the footprint, and the compact fraction.
- Figures: one file per panel (PNG + PDF) so the manuscript composes them with LaTeX subfigures.
- EXP-01 panels are drawn wide (they are stacked in the manuscript), the filter intervals are labelled explicitly
  instead of by log minor ticks, and labels that would collide are dropped.
- Repository housekeeping: `CITATION.cff`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `.zenodo.json`, a CI workflow
  that runs the tests on Python 3.10 and 3.12, and status badges in the README.
- Packaging for PyPI: complete metadata and project URLs in `pyproject.toml`, absolute image links in the README
  (PyPI does not resolve relative paths), and a GitHub Actions workflow that publishes on release through PyPI
  trusted publishing.
- README: illustrative figures in `docs/img/` (a rendered state, the three regimes, the verification, and nine
  days of vorticity in the `two_jets` regime).
- New `figures/make_globe3d.py`: an illustrative three-dimensional rendering of u, v and h on the sphere
  (shaded globe, depth drawn as exaggerated relief, polar sponge greyed out). Not part of the experiment
  pipeline; `python figures/make_globe3d.py <preset> [state.npy]`.
- Figures: the filter intervals are labelled explicitly on the EXP-01 panels (the log minor ticks were unreadable at
  the printed size), and the displacement panel of EXP-03 shows only the balanced case: once the unbalanced bump has
  radiated it is not a localized object and the displacement of its centroid is not meaningful. Its amplitude, which
  is the informative quantity, stays in the peak panel.

## 0.2.0 (2026-09-19)

**The continuity equation was wrong in 0.1.0 and is fixed.** `dynamics.rhs` computed
`dh/dt = -h*div(V)`, which omits the advection of `h`; the shallow-water equation is
`dh/dt = -div(h V) = -h*div(V) - V.grad(h)`. It is now discretized in flux form. Measured
with the 0.1.0 code: the missing term was 26 % (`waves`) to 87 % (`two_jets`) of the term
that was there (rms), and the global mass of the unforced `waves` preset drifted 0.41 % in
10 days against 0.015 % now. The model went unnoticed because the term vanishes for exactly
geostrophic flow and because the zonal relaxation of the jets restores the mass.
**Every result produced with 0.1.0 (climatologies, stability table, footprints, and any
assimilation experiment on them) must be regenerated.**

Also:
- `initial_condition("rossby", perturbed=True, seed=...)` raised `TypeError`; fixed. New kind `"tc2"`.
- New `swesphere.diagnostics`: mass, total energy, potential enstrophy, Williamson error norms.
- New EXP-00 (verification): Williamson test case 2 and the Galewsky instability.
- EXP-01: dt pushed until RK4 blows up; stationarity judged on the eddy std of h averaged over
  windows (the old criterion on the std of h was true by construction); groups of equal filter
  interval to separate the time step from the dissipation cadence.
- EXP-02: per-snapshot series and first-half / second-half eddy amplitude (the unforced `waves`
  preset decays; it is not stationary).
- Tests for mass conservation, the continuity equation, the TC2 steady state, the perturbed IC.
- `LICENSE` (MIT).
