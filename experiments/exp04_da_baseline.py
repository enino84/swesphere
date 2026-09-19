"""EXP-04 Data-assimilation baseline on two_jets: first analysis, h observed at 10%, LETKF vs modified-Cholesky EnKF."""
import numpy as np, pandas as pd, scipy.sparse as sps
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import spsolve
from swesphere import presets
from common import DAY, QUICK, out, log
S = np.load(out("EXP-02", "clim_two_jets.npy")).astype(float); model, _ = presets.make("two_jets"); G = model.grid
Nlat, Nlon = G.Nlat, G.Nlon; fs = model.field_size; n = model.dim; mean = S.mean(0); A = S - mean
sd = np.array([A[:, k*fs:(k+1)*fs].std() for k in range(3)]); scale = np.repeat(sd, fs); Sn = A/scale; M = len(S)
interior = model.interior_mask(); mI = np.tile(interior, 3); rows_ = np.repeat(np.arange(Nlat), Nlon); cols_ = np.tile(np.arange(Nlon), Nlat)
def neighbours(j, r):
    r0, c0 = rows_[j], cols_[j]; rr = np.arange(max(0, r0-r), min(Nlat-1, r0+r)+1); cc = np.arange(c0-r, c0+r+1) % Nlon
    RR, CC = np.meshgrid(rr, cc, indexing="ij"); o = (RR*Nlon + CC).ravel(); return np.sort(o[interior[o]])
def predecessors(i, r):
    f, j = divmod(i, fs); nb = neighbours(j, r); return np.sort(np.concatenate([nb + k*fs for k in range(f)] + [nb[nb < j] + f*fs]))
def build(DX, r, alpha=0.3):
    rows, cols, vals = [], [], []; d = np.empty(n)
    for i in range(n):
        rows.append(i); cols.append(i); vals.append(1.0)
        if not mI[i]: d[i] = 1.0; continue
        idx = predecessors(i, r); y = DX[i]
        if idx.size == 0: v = y.var(); d[i] = 1/v if v > 0 else 0.0; continue
        X = DX[idx].T; Gm = X.T@X; p = idx.size; Gm.flat[::p+1] += alpha*np.trace(Gm)/p; beta = np.linalg.solve(Gm, X.T@y); v = (y - X@beta).var(); d[i] = 1/v if v > 0 else 0.0
        rows.extend([i]*p); cols.extend(idx.tolist()); vals.extend((-beta).tolist())
    L = csr_matrix((vals, (rows, cols)), shape=(n, n)); return (L.T @ diags(d) @ L).tocsr()
def enkf_mc(Xf, Binv, obs, y, R, rng):
    p = obs.size; H = sps.csr_matrix((np.ones(p), (np.arange(p), obs)), shape=(p, n)); Am = (Binv + (H.T@H)/R).tocsc()
    D = (y[:, None] + np.sqrt(R)*rng.standard_normal((p, Xf.shape[1]))) - H@Xf; Z = spsolve(Am, sps.csc_matrix((H.T@D)/R)); return Xf + (Z.toarray() if sps.issparse(Z) else np.asarray(Z))
def letkf(Xf, obs, y, R, r):
    N = Xf.shape[1]; DX = Xf - Xf.mean(1, keepdims=True); mt = Xf.mean(1); yb = Xf[obs]; Yb = yb - yb.mean(1, keepdims=True); innov = y - yb.mean(1)
    where = np.full(n, -1); where[obs] = np.arange(obs.size); out_ = Xf.copy()
    for i in np.flatnonzero(mI):
        f, j = divmod(i, fs); nb = neighbours(j, r); o = where[np.concatenate([nb + k*fs for k in range(3)])]; o = o[o >= 0]
        if o.size == 0: continue
        Yl = Yb[o]; C = Yl.T/R; P = np.linalg.inv((N-1)*np.eye(N) + C@Yl); w = P@(C@innov[o]); ev, V = np.linalg.eigh((N-1)*P); Wsq = (V*np.sqrt(np.maximum(ev, 0)))@V.T
        out_[i] = mt[i] + DX[i]@w + DX[i]@Wsq
    return out_
N = min(20, M - 1); T = 6*3600.0; R = 0.05**2; seeds = range(2 if QUICK else 10); rows = []
for seed in seeds:
    rng = np.random.default_rng(seed); pick = rng.choice(M, size=N+1, replace=False); X0 = Sn[pick[1:]].T; xt0 = Sn[pick[0]]
    Xf = np.stack([(model.propagate(X0[:, e]*scale + mean, [0.0, T]) - mean)/scale for e in range(N)], 1); xt = (model.propagate(xt0*scale + mean, [0.0, T]) - mean)/scale
    DX = Xf - Xf.mean(1, keepdims=True); hs = np.flatnonzero(interior); obs = 2*fs + np.sort(rng.choice(hs, size=int(0.10*hs.size), replace=False)); y = xt[obs] + 0.05*rng.standard_normal(obs.size)
    e = lambda Q: float(np.sqrt(np.mean((Q.mean(1) - xt)[mI]**2))); eh = lambda Q: float(np.sqrt(np.mean((Q.mean(1) - xt)[2*fs:][interior]**2))); bg, bgh = e(Xf), eh(Xf)
    for r in (1, 2, 3): Q = letkf(Xf, obs, y, R, r); rows.append(dict(seed=seed, method=f"LETKF r={r}", total=e(Q)/bg, h=eh(Q)/bgh))
    for r in (1, 2, 3): Q = enkf_mc(Xf, build(DX, r), obs, y, R, np.random.default_rng(1)); rows.append(dict(seed=seed, method=f"EnKF-MC r={r}", total=e(Q)/bg, h=eh(Q)/bgh))
    log(f"seed {seed} done"); pd.DataFrame(rows).to_csv(out("EXP-04", "da_baseline.csv"), index=False)
print(pd.DataFrame(rows).groupby("method")[["total", "h"]].mean().round(3))
