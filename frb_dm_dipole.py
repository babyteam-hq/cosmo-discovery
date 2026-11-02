"""Sweep row 10: a sky dipole in the extragalactic DM of CHIME/FRB Catalog 2 (dm_exc_ne2001, NE2001 subtracted), beyond
the Galactic/halo model. CHIME sees Dec > -11 only, so the fit is on the northern cap; a Galactic-latitude term
(csc|b|, the halo/disk residual) is always in the model, and the dipole is fit in Galactic and equatorial frames with a
look-elsewhere correction for the 3 free components. Robust statistic: median-based (DM distributions are heavy-tailed),
via quantile (median) regression bootstrapped. usage: python frb_dm_dipole.py data/chimefrbcat2.csv
"""
import sys, numpy as np, pandas as pd
from scipy.stats import chi2 as chi2d
from scipy.optimize import minimize

def main(csv):
    d = pd.read_csv(csv, low_memory=False)
    d = d.drop_duplicates('tns_name')   # one row per source (repeaters' first burst)
    d = d[np.isfinite(d.dm_exc_ne2001) & (np.abs(d.gb) > 10)]
    dm = d.dm_exc_ne2001.values.astype(float); l = np.radians(d.gl.values.astype(float)); b = np.radians(d.gb.values.astype(float))
    ra = np.radians(d.ra.values.astype(float)); dec = np.radians(d.dec.values.astype(float))
    print(f"sources: {len(d)} (|b|>10), DM_exc median {np.median(dm):.0f}, 16-84% {np.percentile(dm,16):.0f}..{np.percentile(dm,84):.0f}")
    y = np.log(np.clip(dm, 20, None))   # log-DM: the scatter is multiplicative
    def X_of(frame):
        if frame == 'gal': n = np.stack([np.cos(b)*np.cos(l), np.cos(b)*np.sin(l), np.sin(b)], 1)
        else: n = np.stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)], 1)
        return np.column_stack([np.ones(len(y)), 1/np.abs(np.sin(b)), n])
    def qreg(X, tau=0.5):
        f = lambda c: np.sum(np.abs(y - X@c)*np.where(y - X@c > 0, tau, 1-tau))
        c0 = np.linalg.lstsq(X, y, rcond=None)[0]; return minimize(f, c0, method='Powell').x
    rng = np.random.default_rng(1)
    for frame in ('gal', 'eq'):
        X = X_of(frame); c = qreg(X); amp = np.linalg.norm(c[2:])
        boots = np.array([qreg(X[i := rng.integers(0, len(y), len(y))]) if False else None for _ in range(0)])
        # bootstrap the dipole vector
        bs = []
        for _ in range(60):
            i = rng.integers(0, len(y), len(y)); Xi, yi = X[i], y[i]
            f = lambda cc: np.sum(np.abs(yi - Xi@cc)*np.where(yi - Xi@cc > 0, 0.5, 0.5))
            bs.append(minimize(f, c, method='Powell').x[2:])
        bs = np.array(bs); covd = np.cov(bs.T); chi = c[2:] @ np.linalg.solve(covd, c[2:])
        print(f"  {frame:3s} frame: dipole (log-DM) {np.round(c[2:],3)} |D|={amp:.3f} (x100% in DM), csc|b| coef {c[1]:.3f}; chi2(3 dof, bootstrap cov)={chi:.1f} p={chi2d.sf(chi,3):.3g}")
    # mean-based GLS with log-DM for comparison, and hemisphere medians
    X = X_of('gal'); c, res, *_ = np.linalg.lstsq(X, y, rcond=None); s2 = np.var(y - X@c); cv = s2*np.linalg.inv(X.T@X)
    chi = c[2:] @ np.linalg.solve(cv[2:, 2:], c[2:]); print(f"  OLS log-DM galactic dipole {np.round(c[2:],3)} chi2={chi:.1f} p={chi2d.sf(chi,3):.3g} (2 frames: Bonferroni x2)")
    for name, sel in [("b>0", b > 0), ("b<0", b < 0), ("l<180", l < np.pi), ("l>180", l > np.pi)]:
        print(f"  median DM_exc {name}: {np.median(dm[sel]):.0f} (n={sel.sum()})")
main(sys.argv[1])
