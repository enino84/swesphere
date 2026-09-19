"""EXP-03 Error dynamics: how a small perturbation of h evolves, balanced against unbalanced.

A Gaussian bump in h (5 % of the climatological std of h, 4 cells wide) is added to the last state of the
record, alone ("h_only") or with the geostrophic wind that balances it ("balanced"), at four sites per preset:
two longitudes on the northern jet core, one on the southern jet core (rows of maximum zonal-mean wind; 45 N and
45 S in the waves preset, whose flow is broad) and one in the quiet subtropics (20 N). The equator is avoided: geostrophic balance does not hold there.

The difference with the unperturbed run is followed hour by hour up to 24 h. Reported at 3, 6, 12, 24 h:
  peak            max |dh| relative to the initial bump
  dx, dy, shift   displacement of the compact part of the perturbation: centroid of dh^2 where |dh| > half its
                  maximum, with a circular mean in longitude (cells; dx eastward, dy northward)
  adv_dx, adv_dy  displacement of a parcel released at the site and advected by the unperturbed wind
  steer_dx, steer_dy  the same with the wind averaged over the footprint of the bump (Gaussian weights, 4 cells):
                  the steering flow that a perturbation of that size feels
  compact         fraction of sum(dh^2) within 8 cells of the tracked centre (the rest has radiated away)
"""
import numpy as np, pandas as pd
from swesphere import presets
from common import DAY, QUICK, out, log

TS = (3, 6) if QUICK else (3, 6, 12, 24)


def centre(d, Nlon):
    w = np.where(np.abs(d) >= 0.5 * np.abs(d).max(), d ** 2, 0.0); w = w / w.sum(); yy, xx = np.mgrid[0:d.shape[0], 0:Nlon]
    ang = np.angle((w * np.exp(2j * np.pi * xx / Nlon)).sum()); return float((w * yy).sum()), float((ang % (2 * np.pi)) * Nlon / (2 * np.pi))


def wrap(dx, Nlon):
    return (dx + Nlon / 2) % Nlon - Nlon / 2


rows = []
for name in ("waves", "two_jets"):
    S = np.load(out("EXP-02", f"clim_{name}.npy"), mmap_mode="r"); model, _ = presets.make(name); G = model.grid; Nlat, Nlon, k = G.Nlat, G.Nlon, model.n_pole_rows
    x0 = np.asarray(S[-1], dtype=float); u0, v0, h0 = model.unpack(x0); fs = model.field_size
    sd = float((np.asarray(S[:, 2 * fs:], dtype=float) - np.asarray(S[:, 2 * fs:], dtype=float).mean(0)).std())
    f = np.asarray(G.f_cor); cosL = np.asarray(G.cos_safe); dlon = 2 * np.pi / Nlon; dlat = np.pi / Nlat; yy, xx = np.mgrid[0:Nlat, 0:Nlon]
    ubar = u0.mean(axis=1); quiet = int(round((90 - 20) * Nlat / 180))
    if name == "waves":                                   # broad zonal flow without a jet core: use 45 N and 45 S
        jetN, jetS = int(round(45 * Nlat / 180)), int(round(135 * Nlat / 180))
    else:
        jetN = k + int(np.argmax(ubar[k:Nlat // 2])); jetS = Nlat // 2 + int(np.argmax(ubar[Nlat // 2:Nlat - k]))
    sites = [("jet", jetN, Nlon // 4), ("jet", jetN, 3 * Nlon // 4), ("jet", jetS, Nlon // 2), ("quiet", quiet, Nlon // 2)]
    for where, r, c in sites:
        lat = 90 - r * 180 / Nlat; bump = 0.05 * sd * np.exp(-((yy - r) ** 2 + wrap(xx - c, Nlon) ** 2) / (2 * 4.0 ** 2))
        du = -G.g * (-np.gradient(bump, axis=0) / (G.a * dlat)) / f[r, c]; dv = G.g * (np.gradient(bump, axis=1) / (G.a * dlon * cosL)) / f[r, c]
        runs = {"balanced": model.pack(u0 + du, v0 + dv, h0 + bump), "h_only": model.pack(u0, v0, h0 + bump)}
        xb = x0.copy(); py, px = float(r), float(c); sy, sx = float(r), float(c)
        for hour in range(1, max(TS) + 1):
            ub, vb, _ = model.unpack(xb); i, j = int(round(py)) % Nlat, int(round(px)) % Nlon                     # parcel moved with the wind at the start of the hour
            px += ub[i, j] * 3600.0 / (G.a * dlon * cosL[i, j]); py -= vb[i, j] * 3600.0 / (G.a * dlat)
            W = np.exp(-((yy - sy) ** 2 + wrap(xx - sx, Nlon) ** 2) / (2 * 4.0 ** 2)); W = W / W.sum()
            sx += float((W * ub / (G.a * dlon * cosL)).sum()) * 3600.0; sy -= float((W * vb).sum()) * 3600.0 / (G.a * dlat)
            xb = model.propagate(xb, [0.0, 3600.0]); runs = {kd: model.propagate(xq, [0.0, 3600.0]) for kd, xq in runs.items()}
            if hour in TS:
                for kd, xq in runs.items():
                    d = model.unpack(xq - xb)[2] / (0.05 * sd); cy, cx = centre(d, Nlon); dx, dy = wrap(cx - c, Nlon), -(cy - r)
                    near = ((yy - cy) ** 2 + wrap(xx - cx, Nlon) ** 2) <= 8.0 ** 2
                    rows.append(dict(preset=name, where=where, row=r, col=c, lat=lat, kind=kd, T_h=hour, peak=float(np.abs(d).max()), dx_cells=float(dx), dy_cells=float(dy),
                                     shift_cells=float(np.hypot(dx, dy)), adv_dx=float(px - c), adv_dy=float(-(py - r)), advection_cells=float(np.hypot(px - c, py - r)), steer_dx=float(sx - c), steer_dy=float(-(sy - r)),
                                     compact=float((d ** 2)[near].sum() / (d ** 2).sum())))
        log(f"{name} {where} site (row {r}, lat {lat:.0f}, col {c}) done")
        pd.DataFrame(rows).to_csv(out("EXP-03", "footprint.csv"), index=False)
