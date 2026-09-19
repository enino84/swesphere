"""EXP-01 Euler vs RK4: energy series until blow-up or 200 days, for each preset."""
import numpy as np, pandas as pd
from swesphere import presets, SWEModel
from common import DAY, QUICK, out, log, zonal_stats
days_max = 20 if QUICK else 200; every = 5 if QUICK else 5
rows = []
for name in ("waves", "one_jet", "two_jets"):
    for scheme, dt in (("euler", 120.0), ("euler", 60.0), ("rk4", 120.0)):
        model, x = presets.make(name, dt=dt, scheme=scheme)
        t = 0.0
        while t < days_max:
            try:
                x = model.propagate(x, [0.0, every * DAY]); t += every
                rows.append(dict(preset=name, scheme=scheme, dt=dt, day=t, **zonal_stats(model, x), blew_up=False))
            except RuntimeError:
                rows.append(dict(preset=name, scheme=scheme, dt=dt, day=t + every, blew_up=True)); log(f"{name} {scheme} dt={dt}: blew up at ~day {t + every}"); break
        log(f"{name} {scheme} dt={dt}: reached day {t}")
        pd.DataFrame(rows).to_csv(out("EXP-01", "stability.csv"), index=False)
