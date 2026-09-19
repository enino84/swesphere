"""EXP-00 Verification against standard shallow-water test cases.

A. Williamson et al. (1992) test case 2 (steady zonal geostrophic flow, alpha = 0): normalized
   l1, l2, l-inf errors of h over 5 days, at LMAX 32 and 64, without the polar sponge (accuracy
   and order of the discretization) and with the default sponge (what the polar treatment costs).
   Errors are measured on the rows outside the default sponge in both cases.
B. Galewsky et al. (2004) barotropic instability (northern jet, unforced, 10 days): drift of mass,
   total energy and potential enstrophy, and the relative vorticity at days 6 and 10 for the figure
   (day 6 is the published reference time; at LMAX 32 the instability rolls up later).
"""
import time, numpy as np, pandas as pd
from swesphere import SWEModel, diagnostics as D
from swesphere.dynamics import relative_vorticity
from common import DAY, QUICK, out, log

configs = [(32, 120.0)] if QUICK else [(32, 120.0), (64, 60.0)]
sponge_rows = lambda L: round(10 * (L + 1) / 33)          # the same latitudes (poleward of ~63 deg) at every resolution

# ---- A. Williamson test case 2
rows = []
for L, dt in configs:
    for sponge in (False, True):
        k = sponge_rows(L); m = SWEModel(LMAX=L, dt=dt, n_pole_rows=(k if sponge else 0)); G = m.grid
        x0 = m.initial_condition("tc2", U0=2 * np.pi * G.a / (12 * DAY), H0=2.94e4 / G.g); x = x0.copy(); i0 = D.invariants(m, x0); t0 = time.time()
        for day in range(1, (2 if QUICK else 5) + 1):
            x = m.propagate(x, [0.0, DAY]); m.n_pole_rows = k; e = D.error_norms(m, x, x0, "h", interior=True); m.n_pole_rows = (k if sponge else 0)
            i = D.invariants(m, x)
            rows.append(dict(case="tc2", LMAX=L, dt=dt, sponge=sponge, day=day, **{f"h_{n}": v for n, v in e.items()}, **{f"drift_{n}": (i[n] - i0[n]) / i0[n] for n in i0}))
        log(f"TC2 LMAX={L} sponge={sponge}: day {day} l2(h) = {rows[-1]['h_l2']:.2e}, mass drift {rows[-1]['drift_mass']:.1e} ({time.time() - t0:.0f} s)")
        pd.DataFrame(rows).to_csv(out("EXP-00", "tc2_norms.csv"), index=False)
d = pd.DataFrame(rows); last = d[(d.day == d.day.max()) & (~d.sponge)].set_index("LMAX").h_l2
if len(last) > 1:
    log(f"TC2 without sponge: l2 {last[32]:.2e} -> {last[64]:.2e}, observed order {np.log2(last[32] / last[64]):.2f}")

# ---- B. Galewsky barotropic instability, unforced
rows = []
for L, dt in configs:
    m = SWEModel(LMAX=L, dt=dt, n_pole_rows=sponge_rows(L)); x = m.initial_condition("galewsky", both_hemispheres=False); i0 = D.invariants(m, x)
    for step in range(1, (4 if QUICK else 20) + 1):
        x = m.propagate(x, [0.0, 0.5 * DAY]); i = D.invariants(m, x); u, v, h = m.unpack(x)
        rows.append(dict(case="galewsky", LMAX=L, dt=dt, day=0.5 * step, U_max=float(np.hypot(u, v).max()), **{f"drift_{n}": (i[n] - i0[n]) / i0[n] for n in i0}))
        if 0.5 * step in ((2.0,) if QUICK else (6.0, 10.0)):
            np.save(out("EXP-00", f"galewsky_vorticity_L{L}_day{0.5 * step:g}.npy"), relative_vorticity(u, v, m.grid).astype(np.float32))
    log(f"Galewsky LMAX={L}: day {rows[-1]['day']:g} drifts mass {rows[-1]['drift_mass']:.1e}, energy {rows[-1]['drift_energy']:.1e}, enstrophy {rows[-1]['drift_enstrophy']:.1e}")
    pd.DataFrame(rows).to_csv(out("EXP-00", "galewsky_invariants.csv"), index=False)
