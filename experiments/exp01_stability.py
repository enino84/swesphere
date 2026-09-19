"""EXP-01 The stable envelope of the model: RK4 with dt in {60, 120, 240, 480} s and the spectral
filter applied every {1, 3, 6} steps, 100 days per preset. Reports whether each setting stays
stationary (the std of h in the last 40 days within 10% of the days 40-60 value) and the cost."""
import time, numpy as np, pandas as pd
from swesphere import presets
from common import DAY, QUICK, out, log, zonal_stats
days_max = 20 if QUICK else 100; every = 10 if QUICK else 20
grid = [(60.0, 3), (120.0, 1), (120.0, 3), (120.0, 6), (240.0, 3), (480.0, 3)] if not QUICK else [(120.0, 3), (240.0, 3)]
rows = []
for name in ("waves", "one_jet", "two_jets"):
    for dt, fe in grid:
        model, x = presets.make(name, dt=dt, filter_every=fe); t = 0.0; series = []; t0 = time.time(); ok = True
        while t < days_max:
            try:
                x = model.propagate(x, [0.0, every * DAY]); t += every; series.append((t, zonal_stats(model, x)))
            except RuntimeError:
                ok = False; break
        h_std = [s["h_std"] for _, s in series]
        stationary = ok and len(h_std) >= 3 and abs(h_std[-1] - h_std[len(h_std)//2]) / max(h_std[len(h_std)//2], 1e-9) < 0.10
        rows.append(dict(preset=name, dt=dt, filter_every=fe, reached_day=t, stable=ok, stationary=stationary,
                         h_std_end=(h_std[-1] if h_std else np.nan), eddy_std_end=(series[-1][1]["eddy_std"] if series else np.nan),
                         U_max_end=(series[-1][1]["U_max"] if series else np.nan), seconds_per_day=(time.time() - t0) / max(t, 1)))
        log(f"{name} dt={dt:g} filter_every={fe}: {'stationary' if stationary else ('stable' if ok else 'BLEW UP at day ' + str(t))}, {rows[-1]['seconds_per_day']:.1f} s/day")
        pd.DataFrame(rows).to_csv(out("EXP-01", "stability.csv"), index=False)
