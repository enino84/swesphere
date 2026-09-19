"""Figures and tables of the article from results/.

Every panel is written as its own file (PNG at 300 dpi and PDF) so that the manuscript composes multi-panel
figures with LaTeX subfigures and subcaptions; no panel carries a title, only axis and colour-bar labels.

  EXP-00  tc2_error, galewsky_invariants, galewsky_L<L>_day<d>
  EXP-01  stability_<preset>                         (+ tables stability_*.csv)
  EXP-02  regime_series_eddy, regime_series_umax, regime_<preset>_{wind,vorticity,eddy}
  EXP-03  footprint_<preset>                          (+ table footprint_table.csv)
"""
import os, sys, glob, re, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "experiments"))
from common import out, RESULTS
from swesphere import presets
from swesphere.dynamics import relative_vorticity

FIG = os.path.join(RESULTS, "figures"); os.makedirs(FIG, exist_ok=True)
PRESETS = ("waves", "one_jet", "two_jets")
COLORS = {"waves": "#0072B2", "one_jet": "#E69F00", "two_jets": "#009E73"}          # Okabe-Ito, colour-blind safe
plt.rcParams.update({"font.size": 10, "axes.labelsize": 10, "legend.fontsize": 8, "xtick.labelsize": 9, "ytick.labelsize": 9,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.03})


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300)
    plt.close(fig)


def globe(F, name, label, cmap, symmetric, mask_rows=0):
    """One Mollweide panel with a horizontal colour bar; rows under the polar sponge are hatched out in grey."""
    N, M = F.shape; lat = np.radians(90 - np.arange(N) * 180 / N); lon = np.radians(np.arange(M) * 360 / M - 180)
    core = F[mask_rows:N - mask_rows] if mask_rows else F
    vmax = float(np.nanpercentile(np.abs(core), 99.5)) if symmetric else float(np.nanpercentile(core, 99.8))
    G = np.roll(F, M // 2, axis=1).astype(float)
    if mask_rows:
        G[:mask_rows] = np.nan; G[N - mask_rows:] = np.nan
    fig = plt.figure(figsize=(3.4, 2.3)); ax = fig.add_subplot(111, projection="mollweide"); ax.set_facecolor("0.85")
    im = ax.pcolormesh(lon, lat, G, cmap=cmap, shading="auto", vmin=(-vmax if symmetric else 0.0), vmax=vmax, rasterized=True)
    ax.grid(alpha=.35, lw=.4); ax.set_xticklabels([]); ax.set_yticklabels([])
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.06, pad=0.04, aspect=30); cb.set_label(label); cb.ax.tick_params(labelsize=8)
    save(fig, name)


# ---- EXP-00 verification
if os.path.exists(out("EXP-00", "tc2_norms.csv")):
    d = pd.read_csv(out("EXP-00", "tc2_norms.csv")); g = pd.read_csv(out("EXP-00", "galewsky_invariants.csv"))
    t = d[d.day == d.day.max()].pivot_table(index=["LMAX", "sponge"], values=["h_l1", "h_l2", "h_linf", "drift_mass", "drift_energy"])
    t.to_csv(os.path.join(FIG, "tc2_table.csv")); print("Williamson TC2, last day:"); print(t)
    fig, ax = plt.subplots(figsize=(3.6, 3.1))
    for (L, sp), e in d.groupby(["LMAX", "sponge"]):
        ax.semilogy(e.day, e.h_l2, marker=("o" if sp else "s"), ls=("-" if sp else "--"), color=("#D55E00" if L == 32 else "#0072B2"), ms=4,
                    label=f"L = {L}, {'sponge' if sp else 'no sponge'}")
    ax.set_xlabel("day"); ax.set_ylabel(r"normalized $l_2$ error of $h$"); ax.legend(frameon=False, ncol=2, fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 1.0)); ax.grid(alpha=.3, which="both", lw=.4); save(fig, "tc2_error")
    fig, ax = plt.subplots(figsize=(3.6, 3.1))
    for L, e in g.groupby("LMAX"):
        for n, c in (("mass", "#0072B2"), ("energy", "#E69F00"), ("enstrophy", "#009E73")):
            ax.semilogy(e.day, np.abs(e[f"drift_{n}"]) + 1e-12, ls=("-" if L == 32 else "--"), color=c, label=f"{n}, L = {L}")
    ax.set_xlabel("day"); ax.set_ylabel("|relative drift|"); ax.legend(frameon=False, ncol=2, fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 1.0)); ax.grid(alpha=.3, which="both", lw=.4); save(fig, "galewsky_invariants")
    for f in sorted(glob.glob(out("EXP-00", "galewsky_vorticity_L*_day*.npy"))):
        L, day = re.findall(r"L(\d+)_day([\d.]+)\.npy", f)[0]; z = np.load(f); N = z.shape[0]
        fig, ax = plt.subplots(figsize=(3.6, 1.7))
        cf = ax.contourf(np.arange(2 * N) * 180 / N, 90 - np.arange(N) * 180 / N, z * 1e5, levels=np.arange(-11, 16, 2), cmap="RdBu_r", extend="both")
        ax.set_ylim(10, 80); ax.set_xlabel("longitude (deg)"); ax.set_ylabel("latitude (deg)"); ax.set_xticks([0, 90, 180, 270, 360])
        cb = fig.colorbar(cf, ax=ax, pad=0.02); cb.set_label(r"$\zeta$ ($10^{-5}$ s$^{-1}$)", fontsize=8); cb.ax.tick_params(labelsize=7)
        save(fig, f"galewsky_L{L}_day{float(day):g}")

# ---- EXP-01 stable envelope and climate against the filter interval
if os.path.exists(out("EXP-01", "stability.csv")):
    d = pd.read_csv(out("EXP-01", "stability.csv"))
    for col in ("stable", "stationary", "eddy_mean"):
        t = d.pivot_table(index=["dt", "filter_every", "filter_interval_s"], columns="preset", values=col, aggfunc="first")[[p for p in PRESETS if p in set(d.preset)]]
        t.to_csv(os.path.join(FIG, f"stability_{col}.csv")); print(f"EXP-01 {col}:"); print(t.round(1) if col == "eddy_mean" else t)
    for name in PRESETS:
        e = d[(d.preset == name) & d.stable]
        if e.empty: continue
        fig, ax = plt.subplots(figsize=(3.4, 2.7)); ax.errorbar(e.filter_interval_s, e.eddy_mean, yerr=e.eddy_sd, fmt="none", ecolor="0.65", lw=.8, zorder=2)
        sc = ax.scatter(e.filter_interval_s, e.eddy_mean, c=np.log2(e.dt), cmap="viridis", s=28, zorder=3, edgecolor="k", linewidth=.3)
        ax.set_xscale("log"); ax.set_xlabel(r"filter interval $\Delta t\,n_\mathrm{f}$ (s)"); ax.set_ylabel(r"eddy std of $h$ (m)"); ax.grid(alpha=.3, which="both", lw=.4)
        ticks = sorted(set(np.log2(e.dt))); ticks = ticks[::max(1, len(ticks) // 5)]; cb = fig.colorbar(sc, ax=ax, pad=0.02); cb.set_ticks(ticks); cb.set_ticklabels([f"{2 ** k:g}" for k in ticks]); cb.set_label(r"$\Delta t$ (s)", fontsize=8); cb.ax.tick_params(labelsize=7)
        save(fig, f"stability_{name}")

# ---- EXP-02 regimes: series and last state of each record
if os.path.exists(out("EXP-02", "regime_series.csv")):
    s = pd.read_csv(out("EXP-02", "regime_series.csv"))
    for col, lab, fname in (("eddy_std", r"eddy std of $h$ (m)", "regime_series_eddy"), ("U_max", r"$|U|_{\max}$ (m s$^{-1}$)", "regime_series_umax")):
        fig, ax = plt.subplots(figsize=(3.6, 2.7))
        for name in PRESETS:
            e = s[s.preset == name]
            if len(e): ax.plot(e.day, e[col], color=COLORS[name], lw=1.1, label=name)
        if col == "eddy_std": ax.set_yscale("log")
        ax.set_xlabel("day"); ax.set_ylabel(lab); ax.legend(frameon=False); ax.grid(alpha=.3, which="both", lw=.4); save(fig, fname)
for name in PRESETS:
    path = out("EXP-02", f"clim_{name}.npy")
    if not os.path.exists(path): continue
    S = np.load(path, mmap_mode="r"); model, _ = presets.make(name); u, v, h = model.unpack(np.asarray(S[-1], dtype=float)); k = model.n_pole_rows
    globe(np.hypot(u, v), f"regime_{name}_wind", r"wind speed (m s$^{-1}$)", "viridis", False, k)
    globe(relative_vorticity(u, v, model.grid) * 1e5, f"regime_{name}_vorticity", r"relative vorticity ($10^{-5}$ s$^{-1}$)", "RdBu_r", True, k)
    globe(h - h.mean(axis=1, keepdims=True), f"regime_{name}_eddy", r"$h$ minus zonal mean (m)", "RdBu_r", True, k)

# ---- EXP-03 footprints
if os.path.exists(out("EXP-03", "footprint.csv")):
    f = pd.read_csv(out("EXP-03", "footprint.csv"))
    g = f.groupby(["preset", "where", "kind", "T_h"])[["peak", "compact", "dx_cells", "steer_dx", "adv_dx"]].mean().round(2)
    g.to_csv(os.path.join(FIG, "footprint_table.csv")); print(g)
    for name, e in f.groupby("preset"):
        jet = e[e["where"] == "jet"]; m = jet.groupby(["kind", "T_h"])[["dx_cells", "adv_dx", "steer_dx", "peak"]].mean().reset_index(); sdv = jet.groupby(["kind", "T_h"])["dx_cells"].std().reset_index()
        fig, ax = plt.subplots(figsize=(3.6, 2.7)); b = m[m.kind == "balanced"]; ax.plot(b.T_h, b.adv_dx, "k:", lw=1.2, label="parcel at the jet core"); ax.plot(b.T_h, b.steer_dx, "k--", lw=1.0, label="steering wind (footprint mean)")
        for kind, c, mk, lab in (("balanced", "#0072B2", "o", "balanced"), ("h_only", "#D55E00", "s", "unbalanced ($h$ only)")):
            q = m[m.kind == kind]; ax.errorbar(q.T_h, q.dx_cells, yerr=sdv[sdv.kind == kind].dx_cells.fillna(0), marker=mk, ms=4, lw=1.1, capsize=2, color=c, label=lab)
        ax.set_xlabel("lead time (h)"); ax.set_ylabel("eastward displacement (cells)"); ax.set_xticks(sorted(m.T_h.unique())); ax.legend(frameon=False); ax.grid(alpha=.3, lw=.4); save(fig, f"footprint_{name}")
        fig, ax = plt.subplots(figsize=(3.6, 2.7)); m2 = e.groupby(["kind", "T_h"])["peak"].mean().reset_index()
        for kind, c, mk, lab in (("balanced", "#0072B2", "o", "balanced"), ("h_only", "#D55E00", "s", "unbalanced ($h$ only)")):
            q = m2[m2.kind == kind]; ax.plot([0] + list(q.T_h), [1.0] + list(q.peak), marker=mk, ms=4, lw=1.1, color=c, label=lab)
        ax.set_ylim(0, None); ax.set_xlabel("lead time (h)"); ax.set_ylabel("peak of $|\\delta h|$ / initial"); ax.set_xticks([0] + sorted(m2.T_h.unique())); ax.legend(frameon=False); ax.grid(alpha=.3, lw=.4); save(fig, f"footprint_peak_{name}")
print("figures ->", FIG)
