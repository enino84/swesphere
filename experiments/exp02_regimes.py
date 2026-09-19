"""EXP-02 The three regimes: spin-up, climatology, statistics, decorrelation."""
import numpy as np, pandas as pd
from swesphere import presets, climatology
from common import DAY, QUICK, out, log, zonal_stats
n_snap = 12 if QUICK else 200
rows = []
for name in ("waves", "one_jet", "two_jets"):
    model, x0 = presets.make(name)
    S = climatology.build_climatology(model, x0, spinup_days=(10 if QUICK else presets.SPINUP_DAYS[name]), n_snapshots=n_snap, every_days=presets.SNAPSHOT_DAYS[name], verbose=True)
    np.save(out("EXP-02", f"clim_{name}.npy"), S)
    scales, sd = climatology.anomaly_scales(model, S)
    st = zonal_stats(model, S[-1].astype(float))
    rows.append(dict(preset=name, n_snapshots=len(S), every_days=presets.SNAPSHOT_DAYS[name], anom_u=sd[0], anom_v=sd[1], anom_h=sd[2],
                     corr_lag1=climatology.lag_correlation(model, S.astype(float), 1), corr_lag2=climatology.lag_correlation(model, S.astype(float), 2), **st))
    log(f"{name}: {rows[-1]}")
    pd.DataFrame(rows).to_csv(out("EXP-02", "regimes.csv"), index=False)
