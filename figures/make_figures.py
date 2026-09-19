"""Figures of the article from results/: stability series, regime snapshots, footprints, DA baseline."""
import os, sys, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "experiments"))
from common import out, RESULTS
from swesphere import presets
FIG = os.path.join(RESULTS, "figures"); os.makedirs(FIG, exist_ok=True)
d = pd.read_csv(out("EXP-01", "stability.csv"))
t = d.pivot_table(index=["dt", "filter_every"], columns="preset", values="stationary", aggfunc="first")
t.to_csv(os.path.join(FIG, "stability_table.csv")); print("stable envelope (True = stationary over 100 days):"); print(t)
for name in ("waves", "one_jet", "two_jets"):
    S = np.load(out("EXP-02", f"clim_{name}.npy")); model, _ = presets.make(name); u, v, h = model.unpack(S[-1].astype(float)); G = model.grid
    lats = 90 - np.arange(G.Nlat)*180/G.Nlat; lons = np.arange(G.Nlon)*360/G.Nlon; LON, LAT = np.meshgrid(lons, lats); lon_r, lat_r = np.radians(LON-180), np.radians(LAT)
    fig = plt.figure(figsize=(15, 4.5))
    for pos, (F, t, cmap, sym) in enumerate([(np.hypot(u, v), "wind speed (m/s)", "magma", False), (np.gradient(v, axis=1) - np.gradient(u, axis=0), "relative vorticity (grid units)", "RdBu_r", True), (h - h.mean(axis=1, keepdims=True), "h minus zonal mean (m)", "RdBu_r", True)], 1):
        ax = fig.add_subplot(1, 3, pos, projection="mollweide"); vmax = np.nanpercentile(np.abs(F), 99.5) if sym else None
        im = ax.pcolormesh(lon_r, lat_r, F, cmap=cmap, shading="auto", vmin=(-vmax if sym else None), vmax=(vmax if sym else None)); ax.set_title(t, fontsize=10); ax.grid(alpha=.3); ax.set_xticklabels([]); ax.set_yticklabels([]); plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    fig.suptitle(f"preset '{name}': last climatology snapshot"); fig.tight_layout(); fig.savefig(os.path.join(FIG, f"regime_{name}.png"), dpi=130)
f = pd.read_csv(out("EXP-03", "footprint.csv")); g = f.groupby(["preset", "kind", "T_h"])[["shift_cells", "advection_cells", "spread_cells", "peak"]].mean().round(2); g.to_csv(os.path.join(FIG, "footprint_table.csv")); print(g)
print("figures ->", FIG)
