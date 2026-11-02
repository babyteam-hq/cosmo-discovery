"""Sweep row 23: do Pantheon+ Hubble residuals depend on the local matter density at the SN (an environment step), using
the CF4++ reconstructed density contrast d_mean (Hollinger & Courtois 2025 grid, 1000 Mpc box, supergalactic) at the SN
position? Residual y = m_b_corr - mu(zHD) (so the 2M++ velocity correction is in), full covariance GLS; tests: slope on
delta, a step at delta=0, and on the CF4++ radial velocity residual (what the SN's zHD misses if CF4++ were right).
usage: python env_density.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/new_CF4pp_mean_std_grids.npz
"""
import sys, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.ndimage import map_coordinates
from scipy.stats import chi2 as chi2d
from astropy.coordinates import SkyCoord
import astropy.units as u
C = 299792.458; OM = 0.334; N = 128; L = 1000.0; DELTA = L/N
def dc_h(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3 + 1-OM)
    return C/100*np.array([quad(e, 0, zi)[0] for zi in np.atleast_1d(z)])
def mu(z): return 5*np.log10((1+z)*dc_h(z)/0.7) + 25

def main(dat, covf, grid):
    g = np.load(grid); D = g['d_mean']; Ds = g['d_std']
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f: n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    for zmax in (0.1, 0.06):
        m = (d.IS_CALIBRATOR == 0) & (d.zHD > 0.01) & (d.zHD < zmax); idx = np.where(m)[0]
        zh = d.zHD.values[idx]; r = dc_h(d.zCMB.values[idx])
        sc = SkyCoord(ra=d.RA.values[idx]*u.deg, dec=d.DEC.values[idx]*u.deg, distance=r*u.Mpc, frame='fk5', equinox='J2000').supergalactic
        xyz = np.stack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value]); ci = (xyz + L/2)/DELTA - 0.5
        dl = map_coordinates(D, ci, order=1, mode='nearest'); dls = map_coordinates(Ds, ci, order=1, mode='nearest')
        y = d.m_b_corr.values[idx] - mu(zh)
        Lc = np.linalg.cholesky(cov[np.ix_(idx, idx)]); w = lambda X: np.linalg.solve(Lc, X)
        def fit(cols):
            X = np.column_stack([np.ones(len(y))]+cols); Xw, yw = w(X), w(y); b, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
            e = np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))); rr = yw-Xw@b; return b, e, rr@rr
        _, _, c0 = fit([])
        print(f"zHD<{zmax}: n={len(idx)} delta at SNe: median {np.median(dl):.2f}, 10-90% {np.percentile(dl,10):.2f}..{np.percentile(dl,90):.2f}, grid std median {np.median(dls):.2f}")
        for name, cols, dof in [("slope on delta", [dl], 1), ("step at delta>0", [(dl > 0).astype(float)], 1),
                                ("slope on delta, S/N>1 voxels only (others 0)", [np.where(np.abs(dl) > dls, dl, 0.0)], 1),
                                ("slope on log(1+delta)", [np.log(np.clip(1+dl, 0.05, None))], 1)]:
            b, e, c = fit(cols); dc = c0-c
            print(f"  {name:46s} coef={b[1]:7.4f}±{e[1]:6.4f}  dchi2={dc:5.2f} p={chi2d.sf(dc, dof):.3g}")
main(*sys.argv[1:4])
