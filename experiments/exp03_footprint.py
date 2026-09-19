"""EXP-03 Error dynamics: balanced vs unbalanced perturbation footprints, and how far the flow carries each point."""
import numpy as np, pandas as pd
from swesphere import presets
from common import DAY, QUICK, out, log
rows = []
for name in ("waves", "two_jets"):
    S = np.load(out("EXP-02", f"clim_{name}.npy")).astype(float); model, _ = presets.make(name); G = model.grid
    x0 = S[-1]; u0, v0, h0 = model.unpack(x0); Nlat, Nlon = G.Nlat, G.Nlon
    f = np.asarray(G.f_cor); cosL = np.asarray(G.cos_safe); dlon = 2*np.pi/Nlon; dlat = np.pi/Nlat; yy, xx = np.mgrid[0:Nlat, 0:Nlon]
    sd = (S[:, 2*Nlat*Nlon:] - S[:, 2*Nlat*Nlon:].mean(0)).std()
    for (r, c) in ((Nlat//3, Nlon//2), (Nlat//2, Nlon//4), (2*Nlat//3, 3*Nlon//4)):
        bump = 0.05*sd*np.exp(-((yy-r)**2 + (xx-c)**2)/(2*4.0**2))
        for kind in ("balanced", "h_only"):
            xp = x0.copy()
            if kind == "balanced":
                fs = np.where(np.abs(f) < 2e-5, 2e-5, f)
                du = -G.g*(-np.gradient(bump, axis=0)/(G.a*dlat))/fs; dv = G.g*(np.gradient(bump, axis=1)/(G.a*dlon*cosL))/fs
                xp = model.pack(u0 + du, v0 + dv, h0 + bump)
            else:
                xp = model.pack(u0, v0, h0 + bump)
            xb, xq, t = x0.copy(), xp.copy(), 0.0
            for T in ((3, 6, 12, 24) if not QUICK else (3, 6)):
                xb = model.propagate(xb, [0.0, (T - t)*3600.0]); xq = model.propagate(xq, [0.0, (T - t)*3600.0]); t = T
                d = model.unpack(xq - xb)[2] / (0.05*sd); ad = np.abs(d); w = ad/ad.sum()
                cy, cx = (w*yy).sum(), (w*xx).sum(); spread = np.sqrt((w*((yy-cy)**2 + (xx-cx)**2)).sum())
                adv = np.hypot(u0, v0)[r, c]*T*3600.0/(G.a*dlon*cosL[r, c])
                rows.append(dict(preset=name, point=f"({r},{c})", kind=kind, T_h=T, shift_cells=float(np.hypot(cx-c, cy-r)), advection_cells=float(adv), spread_cells=float(spread), peak=float(ad.max())))
        log(f"{name} point ({r},{c}) done")
    pd.DataFrame(rows).to_csv(out("EXP-03", "footprint.csv"), index=False)
