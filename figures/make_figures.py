"""Figures and tables of the article from results/: verification, stable envelope, regimes, footprints."""
import os, sys, glob, re, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "experiments"))
from common import out, RESULTS
from swesphere import presets
FIG = os.path.join(RESULTS, "figures"); os.makedirs(FIG, exist_ok=True)
PRESETS = ("waves", "one_jet", "two_jets")

# ---- EXP-00 verification
if os.path.exists(out("EXP-00", "tc2_norms.csv")):
    d = pd.read_csv(out("EXP-00", "tc2_norms.csv")); g = pd.read_csv(out("EXP-00", "galewsky_invariants.csv"))
    t = d[d.day == d.day.max()].pivot_table(index=["LMAX", "sponge"], values=["h_l1", "h_l2", "h_linf", "drift_mass", "drift_energy"]); t.to_csv(os.path.join(FIG, "tc2_table.csv")); print("Williamson TC2, last day:"); print(t)
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    for (L, sp), e in d.groupby(["LMAX", "sponge"]):
        ax[0].semilogy(e.day, e.h_l2, "o-" if sp else "s--", label=f"LMAX {L}, {'default sponge' if sp else 'no sponge'}")
    ax[0].set_xlabel("day"); ax[0].set_ylabel("normalized $l_2$ error of h"); ax[0].set_title("Williamson test case 2"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
    for L, e in g.groupby("LMAX"):
        for n, ls in (("mass", "-"), ("energy", "--"), ("enstrophy", ":")):
            ax[1].semilogy(e.day, np.abs(e[f"drift_{n}"]) + 1e-12, ls, label=f"{n}, LMAX {L}")
    ax[1].set_xlabel("day"); ax[1].set_ylabel("|relative drift|"); ax[1].set_title("Galewsky jet, unforced: invariants"); ax[1].legend(fontsize=7, ncol=2); ax[1].grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "verification.png"), dpi=130); plt.close(fig)
    files = sorted(glob.glob(out("EXP-00", "galewsky_vorticity_L*_day*.npy")), key=lambda f: [float(s) for s in re.findall(r"L(\d+)_day([\d.]+)\.npy", f)[0]])
    fig, axes = plt.subplots(len(files), 1, figsize=(10, 2.4 * len(files)), squeeze=False)
    for ax_, f in zip(axes[:, 0], files):
        z = np.load(f); N = z.shape[0]; L, day = re.findall(r"L(\d+)_day([\d.]+)\.npy", f)[0]
        ax_.contourf(np.arange(2 * N) * 180 / N, 90 - np.arange(N) * 180 / N, z, levels=np.arange(-1.1e-4, 1.6e-4, 2e-5), cmap="RdBu_r", extend="both"); ax_.set_ylim(10, 80); ax_.set_title(f"relative vorticity, LMAX {L}, day {day}", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "galewsky_vorticity.png"), dpi=130); plt.close(fig)

# ---- EXP-01 stable envelope and climate against the filter interval
d = pd.read_csv(out("EXP-01", "stability.csv"))
for col in ("stable", "stationary", "eddy_mean"):
    t = d.pivot_table(index=["dt", "filter_every", "filter_interval_s"], columns="preset", values=col, aggfunc="first")[[p for p in PRESETS if p in set(d.preset)]]
    t.to_csv(os.path.join(FIG, f"stability_{col}.csv")); print(f"EXP-01 {col}:"); print(t.round(1) if col == "eddy_mean" else t)
fig, ax = plt.subplots(1, 3, figsize=(14, 3.8))
for a, name in zip(ax, PRESETS):
    e = d[(d.preset == name) & d.stable]; sc = a.scatter(e.filter_interval_s, e.eddy_mean, c=np.log2(e.dt), cmap="viridis", s=45, zorder=3); a.errorbar(e.filter_interval_s, e.eddy_mean, yerr=e.eddy_sd, fmt="none", ecolor="0.6", zorder=2)
    a.set_xscale("log"); a.set_xlabel("filter interval dt x filter_every (s)"); a.set_ylabel("eddy std of h, days 40-100 (m)"); a.set_title(name); a.grid(alpha=.3)
    cb = plt.colorbar(sc, ax=a); ticks = sorted(set(np.log2(e.dt))); cb.set_ticks(ticks[::2]); cb.set_ticklabels([f"{2 ** t:g}" for t in ticks[::2]]); cb.set_label("dt (s)")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "stability_climate.png"), dpi=130); plt.close(fig)

# ---- EXP-02 regimes: series and last snapshot
if os.path.exists(out("EXP-02", "regime_series.csv")):
    s = pd.read_csv(out("EXP-02", "regime_series.csv")); fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
    for name, e in s.groupby("preset"):
        ax[0].plot(e.day, e.eddy_std, label=name); ax[1].plot(e.day, e.U_max, label=name)
    ax[0].set_ylabel("eddy std of h (m)"); ax[1].set_ylabel("|U| max (m/s)")
    for a in ax: a.set_xlabel("day"); a.legend(); a.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "regime_series.png"), dpi=130); plt.close(fig)
for name in PRESETS:
    S = np.load(out("EXP-02", f"clim_{name}.npy")); model, _ = presets.make(name); u, v, h = model.unpack(S[-1].astype(float)); G = model.grid
    lats = 90 - np.arange(G.Nlat)*180/G.Nlat; lons = np.arange(G.Nlon)*360/G.Nlon; LON, LAT = np.meshgrid(lons, lats); lon_r, lat_r = np.radians(LON-180), np.radians(LAT)
    from swesphere.dynamics import relative_vorticity
    fig = plt.figure(figsize=(15, 4.5))
    for pos, (F, t, cmap, sym) in enumerate([(np.hypot(u, v), "wind speed (m/s)", "magma", False), (relative_vorticity(u, v, G), "relative vorticity (1/s)", "RdBu_r", True), (h - h.mean(axis=1, keepdims=True), "h minus zonal mean (m)", "RdBu_r", True)], 1):
        ax = fig.add_subplot(1, 3, pos, projection="mollweide"); vmax = np.nanpercentile(np.abs(F[model.n_pole_rows:-model.n_pole_rows]), 99.5) if sym else None
        im = ax.pcolormesh(lon_r, lat_r, F, cmap=cmap, shading="auto", vmin=(-vmax if sym else None), vmax=(vmax if sym else None)); ax.set_title(t, fontsize=10); ax.grid(alpha=.3); ax.set_xticklabels([]); ax.set_yticklabels([]); plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    fig.suptitle(f"preset '{name}': last snapshot of the record"); fig.tight_layout(); fig.savefig(os.path.join(FIG, f"regime_{name}.png"), dpi=130); plt.close(fig)

# ---- EXP-03 footprints
f = pd.read_csv(out("EXP-03", "footprint.csv")); g = f.groupby(["preset", "kind", "T_h"])[["shift_cells", "advection_cells", "spread_cells", "peak"]].mean().round(2); g.to_csv(os.path.join(FIG, "footprint_table.csv")); print(g)
print("figures ->", FIG)
