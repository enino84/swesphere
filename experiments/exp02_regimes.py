"""EXP-02 The three regimes: spin-up, 200-state record, statistics, stationarity and decorrelation.

Besides the record itself (``clim_<preset>.npy``) the eddy std of h, the wind maximum and the invariants of every
snapshot go to ``regime_series.csv``; ``regimes.csv`` summarizes them, with ``eddy_first_half`` and
``eddy_second_half`` to show whether the regime is stationary (forced jets) or decaying (unforced waves).
"""
import numpy as np, pandas as pd
from swesphere import presets, climatology, diagnostics
from common import DAY, QUICK, out, log, zonal_stats
n_snap = 12 if QUICK else 200
rows, series = [], []
for name in ("waves", "one_jet", "two_jets"):
    model, x0 = presets.make(name); spin = 10 if QUICK else presets.SPINUP_DAYS[name]; every = presets.SNAPSHOT_DAYS[name]
    S = climatology.build_climatology(model, x0, spinup_days=spin, n_snapshots=n_snap, every_days=every, verbose=True, label=name)
    np.save(out("EXP-02", f"clim_{name}.npy"), S)
    for k, s in enumerate(S):
        series.append(dict(preset=name, snapshot=k, day=spin + k * every, **zonal_stats(model, s.astype(float)), **diagnostics.invariants(model, s.astype(float))))
    e = pd.DataFrame([r for r in series if r["preset"] == name]); half = len(e) // 2
    scales, sd = climatology.anomaly_scales(model, S)
    rows.append(dict(preset=name, n_snapshots=len(S), every_days=every, anom_u=sd[0], anom_v=sd[1], anom_h=sd[2],
                     corr_lag1=climatology.lag_correlation(model, S.astype(float), 1), corr_lag2=climatology.lag_correlation(model, S.astype(float), 2),
                     eddy_mean=e.eddy_std.mean(), eddy_sd=e.eddy_std.std(), eddy_first_half=e.eddy_std[:half].mean(), eddy_second_half=e.eddy_std[half:].mean(),
                     h_std=e.h_std.mean(), U_mean=e.U_mean.mean(), U_max=e.U_max.mean(), mass_drift=e["mass"].iloc[-1] / e["mass"].iloc[0] - 1))
    log(f"{name}: {rows[-1]}")
    pd.DataFrame(rows).to_csv(out("EXP-02", "regimes.csv"), index=False); pd.DataFrame(series).to_csv(out("EXP-02", "regime_series.csv"), index=False)
