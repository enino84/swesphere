# -*- coding: utf-8 -*-
"""Climatology builder: spin-up, snapshots, decorrelation check, normalization."""
from __future__ import annotations

import numpy as np

DAY = 86400.0


def build_climatology(model, x0, spinup_days, n_snapshots, every_days, verbose=False):
    x = model.propagate(x0, [0.0, spinup_days * DAY])
    snaps = [x.astype(np.float32)]
    for k in range(n_snapshots - 1):
        x = model.propagate(x, [0.0, every_days * DAY]); snaps.append(x.astype(np.float32))
        if verbose and (k + 1) % 10 == 0:
            print(f"  {k + 1} snapshots", flush=True)
    return np.array(snaps)


def anomaly_scales(model, snapshots):
    """Per-field standard deviation of the climatological anomalies, repeated over the state."""
    A = snapshots - snapshots.mean(axis=0); fs = model.field_size
    sd = np.array([A[:, k * fs:(k + 1) * fs].std() for k in range(3)])
    return np.repeat(sd, fs), sd


def lag_correlation(model, snapshots, lag=1, field="h"):
    A = snapshots - snapshots.mean(axis=0); blk = model.var_blocks[field]; a = A[:, blk]
    return float(np.mean([np.corrcoef(a[i], a[i + lag])[0, 1] for i in range(len(a) - lag)]))
