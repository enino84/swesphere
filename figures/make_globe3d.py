"""A three-dimensional view of a model state: u, v and h rendered on the sphere.

Not part of the paper's experiment pipeline; this is an illustrative rendering. Usage:

    python figures/make_globe3d.py [preset] [state.npy]

With no state file the script spins the preset up for 40 days first. The height field is drawn as a
radial displacement of the sphere (exaggerated) coloured by the field itself; u and v are drawn on the
undeformed sphere. Rows under the polar sponge are left grey.
"""
import os, sys, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import pyshtools as pysh

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "experiments"))
from swesphere import presets

DAY = 86400.0
name = sys.argv[1] if len(sys.argv) > 1 else "two_jets"
model, x0 = presets.make(name)
if len(sys.argv) > 2:
    x = np.load(sys.argv[2]).astype(float)
else:
    x = model.propagate(x0, [0.0, 40 * DAY])
u, v, h = model.unpack(x); G = model.grid; k = model.n_pole_rows
UP = 4                                                           # render mesh refinement


def refine(F):
    """Synthesize the band-limited field on a finer grid: the picture is smooth without inventing structure."""
    c = pysh.expand.SHExpandDH(F, norm=1, sampling=2, lmax_calc=G.LMAX)
    L = UP * (G.LMAX + 1) - 1
    cc = np.zeros((2, L + 1, L + 1)); cc[:, :G.LMAX + 1, :G.LMAX + 1] = c
    return pysh.expand.MakeGridDH(cc, norm=1, sampling=2, lmax=L)


NLAT = 2 * UP * (G.LMAX + 1); NLON = 2 * NLAT
KK = UP * k                                                      # the sponge rows on the fine mesh
lat = np.radians(90 - np.arange(NLAT) * 180 / NLAT)
lon = np.radians(np.arange(NLON + 1) * 360 / NLON)               # wrap so the sphere closes
LON, LAT = np.meshgrid(lon, lat)


def panel(ax, F, cmap, label, relief=0.0, symmetric=True, title=""):
    """Draw one field on the sphere, with Lambert shading from a fixed light so the globe reads as solid."""
    Ff = refine(np.asarray(F, dtype=float))
    Fw = np.concatenate([Ff, Ff[:, :1]], axis=1)
    core = Fw[KK:NLAT - KK]
    vmax = float(np.nanpercentile(np.abs(core), 99.0)) if symmetric else float(np.nanpercentile(core, 99.5))
    vmin = -vmax if symmetric else float(np.nanpercentile(core, 0.5))
    norm = colors.Normalize(vmin, vmax)

    scaled = np.clip((Fw - core.mean()) / max(vmax - vmin, 1e-9), -1.5, 1.5)
    scaled[:KK] = scaled[NLAT - KK:] = 0.0                       # no relief under the polar sponge
    r = 1.0 + relief * scaled
    X = r * np.cos(LAT) * np.cos(LON); Y = r * np.cos(LAT) * np.sin(LON); Z = r * np.sin(LAT)

    C = matplotlib.colormaps[cmap](norm(Fw))
    C[:KK] = C[NLAT - KK:] = (0.78, 0.78, 0.80, 1.0)             # polar sponge
    light = np.array([-0.55, -0.62, 0.56]); light /= np.linalg.norm(light)
    el, az = np.radians(26.0), np.radians(-58.0)                 # the camera, to match view_init
    eye = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    nrm = np.stack([np.cos(LAT) * np.cos(LON), np.cos(LAT) * np.sin(LON), np.sin(LAT)], axis=-1)
    lam = np.clip((nrm * light).sum(-1), 0.0, 1.0)
    half = light + eye; half /= np.linalg.norm(half)
    spec = np.clip((nrm * half).sum(-1), 0.0, 1.0) ** 28         # a soft highlight
    rim = np.clip(1.0 - (nrm * eye).sum(-1), 0.0, 1.0) ** 3      # limb brightening
    C[..., :3] *= (0.50 + 0.52 * lam)[..., None]                 # ambient + diffuse
    C[..., :3] = np.clip(C[..., :3] + 0.15 * spec[..., None] + 0.07 * rim[..., None], 0.0, 1.0)

    ax.plot_surface(X, Y, Z, facecolors=C, rstride=1, cstride=1, linewidth=0, antialiased=True, shade=False)
    ax.set_box_aspect((1, 1, 1)); ax.set_axis_off(); ax.view_init(elev=26, azim=-58)
    lim = 0.74
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    m = cm.ScalarMappable(norm=norm, cmap=cmap); m.set_array([])
    cb = plt.colorbar(m, ax=ax, orientation="horizontal", fraction=0.04, pad=0.0, aspect=24, shrink=0.82)
    cb.set_label(label, fontsize=15); cb.ax.tick_params(labelsize=13)
    if title and not PAPER:
        ax.set_title(title, fontsize=16, pad=-4)


PAPER = os.environ.get("PAPER", "0") == "1"
DARK = os.environ.get("DARK", "0") == "1"          # transparent background, light labels (for a dark page)
plt.rcParams.update({"font.size": 13})
if DARK:
    plt.rcParams.update({"text.color": "#EAF2FA", "axes.labelcolor": "#EAF2FA", "xtick.color": "#EAF2FA", "ytick.color": "#EAF2FA"})
fig = plt.figure(figsize=(12, 4.6))
for i, (F, cmap, lab, rel, sym, ttl) in enumerate([
        (u, "RdBu_r", r"$u$ (m s$^{-1}$)", 0.0, True, "zonal wind"),
        (v, "PuOr_r", r"$v$ (m s$^{-1}$)", 0.0, True, "meridional wind"),
        (h - h.mean(axis=1, keepdims=True), "Spectral_r", r"$h$ minus zonal mean (m)", 0.09, True,
         "depth anomaly (relief exaggerated)")], 1):
    panel(fig.add_subplot(1, 3, i, projection="3d"), F, cmap, lab, rel, sym, ttl)
fig.subplots_adjust(left=0.0, right=1.0, top=0.97, bottom=0.06, wspace=0.0)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"globe3d_{name}{'_dark' if DARK else ''}.png")
fig.savefig(out, dpi=220, bbox_inches="tight", transparent=DARK)

if not DARK:
    fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight"); print("->", out)
