"""EXP-01 The stable envelope of the shipped RK4 and the dependence of the climate on the dissipation cadence.

Every preset is run 100 days for a grid of (dt, filter_every). Two different questions are answered:

* stability (a property of RK4 against the fastest gravity waves): does the run reach day 100 with finite,
  physical fields? dt is pushed until it does not (the edge differs between presets: the gravity-wave speed
  of ``waves``, mean depth ~2.8 km, is about half that of the jets, mean depth 10 km).
* climate (a property of the dissipation): the spectral filter acts every ``filter_every`` steps, i.e. every
  ``dt * filter_every`` seconds, and the polar sponge every step, so the eddy amplitude can change with dt even
  where the run is stable. The grid contains groups with the same filter interval (360 s: (60,6), (120,3),
  (360,1); 720 s: (120,6), (240,3), (720,1)) to separate the two effects.

The eddy std of h (rows outside the sponge) is sampled every 5 days; ``eddy_mean`` and ``eddy_sd`` are over days
40-100 and ``drift`` compares days 70-100 with days 40-70 (``stationary`` when |drift| < 0.15).
"""
import time, warnings, numpy as np, pandas as pd
from swesphere import presets
from common import DAY, QUICK, out, log, zonal_stats
warnings.filterwarnings("ignore", category=RuntimeWarning)
days_max, every = (20, 5) if QUICK else (100, 5); t_a, t_b = (5, 10) if QUICK else (40, 70)
grid = [(120.0, 3), (240.0, 3), (2400.0, 3)] if QUICK else [(60.0, 3), (60.0, 6), (120.0, 1), (120.0, 3), (120.0, 6), (240.0, 3), (360.0, 1), (480.0, 3),
                                                             (720.0, 1), (960.0, 3), (1440.0, 3), (1800.0, 3), (2100.0, 3), (2400.0, 3), (3000.0, 3), (3600.0, 3), (4800.0, 3)]
rows, series = [], []
for name in ("waves", "one_jet", "two_jets"):
    for dt, fe in grid:
        model, x = presets.make(name, dt=dt, filter_every=fe); t = 0.0; t0 = time.time(); ok = True; mine = []
        while t < days_max:
            try:
                x = model.propagate(x, [0.0, every * DAY]); st = zonal_stats(model, x)
                if not st["U_max"] < 400.0: raise RuntimeError("unphysical wind")
            except RuntimeError:
                ok = False; break
            t += every; mine.append(dict(preset=name, dt=dt, filter_every=fe, day=t, **st))
        series += mine; e = pd.DataFrame(mine)
        a = e[(e.day > t_a) & (e.day <= t_b)].eddy_std if ok else pd.Series(dtype=float); b = e[e.day > t_b].eddy_std if ok else pd.Series(dtype=float)
        drift = float(b.mean() / a.mean() - 1) if ok else np.nan
        rows.append(dict(preset=name, dt=dt, filter_every=fe, filter_interval_s=dt * fe, reached_day=t, stable=ok, stationary=bool(ok and abs(drift) < 0.15), drift=drift,
                         eddy_mean=(float(pd.concat([a, b]).mean()) if ok else np.nan), eddy_sd=(float(pd.concat([a, b]).std()) if ok else np.nan),
                         h_std_end=(mine[-1]["h_std"] if ok else np.nan), U_max_mean=(float(e[e.day > t_a].U_max.mean()) if ok else np.nan),
                         seconds_per_day=(time.time() - t0) / max(t, every)))
        r = rows[-1]; log(f"{name} dt={dt:g} filter_every={fe}: " + (f"eddy {r['eddy_mean']:.1f} +- {r['eddy_sd']:.1f} m, drift {r['drift']:+.2f}, {r['seconds_per_day']:.2f} s/day" if ok else f"BLEW UP after day {t:g}"))
        pd.DataFrame(rows).to_csv(out("EXP-01", "stability.csv"), index=False); pd.DataFrame(series).to_csv(out("EXP-01", "stability_series.csv"), index=False)
