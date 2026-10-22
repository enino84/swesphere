# Changelog

## 0.6.0 (2026-10-22)

**The polar sponge is off by default.** `n_pole_rows` changes from 10 to 0.
This changes the answer the model gives with default settings, so read the
migration note below before upgrading a running experiment.

### Why

EXP-10 of the model description paper ran 45 configurations for 120 model days
each: three presets, three truncations, five sponge widths from ten rows down
to none. The measurements say the sponge is not worth what it costs.

- **It is not needed for stability.** All 45 runs completed the full 120 days,
  including every sponge-free run at every preset and every resolution.
- **It costs nothing to omit.** Cost per model day with and without it agrees
  to within the scatter of the measurement; it multiplies two arrays per step.
- **It destroys the convergence of the scheme.** On Williamson test case 2 at
  a 90 degree rotation angle, with the sponge covering a fixed latitude band,
  the observed order of the normalized l2 error of the depth is **-0.23**: the
  error grows under grid refinement, from 8.6e-2 at LMAX 16 to 1.3e-1 at
  LMAX 96. With no sponge the same case gives 1.11, and the grid-aligned case
  gives 2.04.
- **At the default resolution it removed more than half of the flow.** The
  standard deviation of the depth anomaly outside 62.7 degrees of latitude,
  after 120 days at LMAX 32: `two_jets` 341 m without the sponge against
  146 m with ten rows; `waves` 34 m against 12 m. Ten rows reach 62.7 degrees
  at LMAX 32, and the damping does not stay inside the band.

What the sponge prevented was grid-scale structure in the polar rows, which
without it reach 0.6 to 1.0 times the midlatitude amplitude in the forced
presets against 0.004 to 0.019 with ten rows. Two things make that acceptable.
It does not grow with resolution: `two_jets` reads 0.59, 0.60 and 0.71 at
LMAX 32, 64 and 96. And those rows were never intended to be used;
`interior_mask()` has always excluded them.

### Migration

**If you want the old behaviour, pass it explicitly:**

    SWEModel(n_pole_rows=10)                 # as before 0.6.0
    presets.two_jets(n_pole_rows=10)

**Check your code for `field[model.n_pole_rows:-model.n_pole_rows]`.** That
idiom was correct while the sponge was always on. With `n_pole_rows = 0` it
returns an **empty array**, because `-0` is `0` and `field[0:0]` selects
nothing, and it fails quietly: the mean and the standard deviation of an empty
array are `nan`, so the first symptom is usually a blank plot rather than an
error. This bit two tests inside this package.

Use the new helper instead, which is correct at every width:

    model.interior(field)        # the rows outside the sponge

`interior_mask()` is unchanged and already returned all-True at zero width.

### Added

- `SWEModel.interior(field)`, the safe form of the row slice above.

### Changed

- `SWEModel(n_pole_rows=...)` defaults to 0 (was 10).
- `grid.make_sponge(n_pole_rows=...)` defaults to 0 (was 10), so the two
  agree.

### Note for anyone reproducing published numbers

Results generated with swesphere 0.5.0 or earlier used a ten-row sponge
unless they said otherwise. To reproduce them, pass `n_pole_rows=10`. The
model description paper reports both configurations and says which is which
in every table.

## 0.5.0 (2026-10-14)

**Standard test cases, and the one thing the model was missing to run them.**
Nothing here changes an existing result: every default is what it was in
0.4.0, and `williamson_tc2(grid)` returns the same three arrays bit for bit.

### Added

- **`SphereGrid.with_coriolis(f)` and `SWEModel.set_coriolis(f)`.** The
  Coriolis field was fixed at the Earth's `f = 2 Omega sin(phi)`, which is
  right for every configuration the package shipped and wrong for the
  standard test case that matters most for verification. `with_coriolis`
  returns a copy of the grid carrying another field; `set_coriolis` does it
  on the model and accepts `None` to go back to the Earth's.

- **Williamson test case 2 at a rotation angle:
  `williamson_tc2(grid, alpha_deg=...)`, `tc2_coriolis`,
  `tc2_sin_rotated_latitude`, and the preset `presets.williamson_tc2`.**

  The grid-aligned case, `alpha_deg = 0`, is a weak test and it was the only
  one available. Its solution has `v == 0` at every point, so the meridional
  momentum equation is only ever evaluated on zero, every meridional
  derivative is taken along a grid line, and the polar caps carry no flow.
  That is the configuration in which this discretization looks its best. At
  `alpha_deg = 90` the same solid-body rotation runs pole to pole, straight
  through the converging meridians and the polar sponge, and the solution is
  still exact and steady, so the whole departure from the initial state is
  the model's own error. Measured over six hours at LMAX 16 with no sponge,
  the depth moves 0.08 percent at `alpha = 0` and 1.0 percent at
  `alpha = 90`; add a four-row polar sponge and the rotated case moves 7.8
  percent instead.

  The rotated test is why `with_coriolis` had to exist. It does not tilt the
  flow and leave the planet alone: it tilts the whole problem, planetary
  rotation included, so it prescribes its own `f = 2 Omega sin(phi')` with
  `phi'` the latitude in the tilted frame. Run the tilted flow against the
  grid's own `f` and it is not a steady solution of anything, because the
  balance

      dB/dlam = a cos(phi) (f + zeta) v,   dB/dphi = -a (f + zeta) u

  then fails its own integrability condition by
  `-2 Omega U0 a sin(alpha) sin(lam) cos^2(phi)`: no height field whatsoever
  holds that flow steady. The drift you would measure would be the drift of a
  problem that was never in balance. Verified symbolically; with the tilted
  `f` both residuals are identically zero. `presets.williamson_tc2` installs
  the right field for you, which is the only way to use this that cannot go
  wrong.

- **Williamson test case 6, the Rossby-Haurwitz wave:
  `rossby_haurwitz(grid, R=4, ...)`, `rh_angular_phase_speed`, the preset
  `presets.rossby_haurwitz`, and `initial_condition("rossby_haurwitz")`
  (also spelled `"tc6"`).** Non-zonal and time dependent, with an analytic
  eastward angular phase speed, `2.46e-6 rad/s` for the standard `R = 4`, a
  pattern that circles the globe in about 29.5 days. Test case 2 asks whether
  the model leaves a structure alone; this one asks whether it carries one at
  the right speed. The `cos^(-2)` term of the published `A(phi)` is folded
  into the power, so nothing divides by `cos(phi)` at the poles.

### Notes

- These are verification cases, not regimes: no spin-up, no forcing, and a
  known answer. The characterization presets `waves`, `one_jet` and
  `two_jets` are unchanged.

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
