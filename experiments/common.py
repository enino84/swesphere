import os, sys, time, json, numpy as np
DAY = 86400.0
RESULTS = os.environ.get("RESULTS_DIR", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results"))
QUICK = os.environ.get("QUICK", "0") == "1"


def out(*parts):
    p = os.path.join(RESULTS, *parts); os.makedirs(os.path.dirname(p), exist_ok=True); return p


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def zonal_stats(model, x):
    u, v, h = model.unpack(x); k = model.n_pole_rows; sl = slice(k, model.Nlat - k)
    hz = h[sl] - h[sl].mean(axis=1, keepdims=True); U = np.hypot(u, v)[sl]
    return dict(h_std=float(h[sl].std()), eddy_std=float(hz.std()), U_mean=float(U.mean()), U_max=float(U.max()))
