"""Which reconstructed velocity field explains Pantheon+ Hubble residuals better: 2M++ (Pantheon+'s own VPEC) or CF4++?

Residual y = m_b_corr - mu(zCMB). A flow model with radial velocity v gives zHD via (1+zHD) = (1+zCMB)/(1+v/c), so a
correction c_model = mu(zHD_model) - mu(zCMB). Fit y = M + beta*c_model with the full STAT+SYS covariance, for 2M++,
for CF4++, and jointly; and the dipole left over after each. CF4++ grid: Hollinger & Courtois 2025 (new_CF4pp_mean_std_grids.npz),
128^3 voxels, L=1000 Mpc/h, supergalactic Cartesian, v in km/s.
usage: python cf4_vs_2mpp.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/new_CF4pp_mean_std_grids.npz
"""
import sys, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.ndimage import map_coordinates
from scipy.stats import chi2 as chi2d
from astropy.coordinates import SkyCoord
import astropy.units as u
C = 299792.458; OM = 0.334; N = 128; L = 1000.0; DELTA = L/N

def dc_h(z):  # comoving distance, Mpc/h
    e = lambda x: 1/np.sqrt(OM*(1+x)**3 + 1-OM)
    return C/100*np.array([quad(e, 0, zi)[0] for zi in np.atleast_1d(z)])
def mu(z): return 5*np.log10((1+z)*dc_h(z)/0.7) + 25
def zhd(zcmb, v): return (1+zcmb)/(1+v/C) - 1

def main(dat, covf, grid, zmax=0.1):
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f:
        n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    m = (d.IS_CALIBRATOR == 0) & (d.zCMB > 0.01) & (d.zCMB < zmax)
    idx = np.where(m)[0]
    zc = d.zCMB.values[idx]; zh2 = d.zHD.values[idx]; vp2 = d.VPEC.values[idx]
    r = dc_h(zc)
    sc = SkyCoord(ra=d.RA.values[idx]*u.deg, dec=d.DEC.values[idx]*u.deg, distance=r*u.Mpc, frame='fk5', equinox='J2000').supergalactic
    xyz = np.stack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value])  # (3, n) Mpc/h
    rhat = xyz/np.linalg.norm(xyz, axis=0)
    g = np.load(grid)
    print("grid keys:", {k: g[k].shape for k in g.files})
    ci = (xyz + L/2)/DELTA - 0.5   # continuous voxel index, voxel centres at (i+0.5)*DELTA - L/2
    interp = lambda a: map_coordinates(a, ci, order=1, mode='nearest')
    v = np.stack([interp(g['v_mean'][k]) for k in range(3)])
    vr_dot = np.sum(v*rhat, 0); vr_grid = interp(g['vr_mean']); vr_std = interp(g['vr_std'])
    print(f"N={len(idx)} z<{zmax}. sanity: corr(v.rhat, vr_mean)={np.corrcoef(vr_dot, vr_grid)[0,1]:.3f}; "
          f"corr(CF4 vr, 2M++ VPEC)={np.corrcoef(vr_dot, vp2)[0,1]:.3f}; rms CF4 vr={np.sqrt(np.mean(vr_dot**2)):.0f}, "
          f"rms 2M++={np.sqrt(np.mean(vp2**2)):.0f}, mean CF4 vr_std={vr_std.mean():.0f} km/s")
    vcf = vr_dot
    y = d.m_b_corr.values[idx] - mu(zc)
    c2 = mu(zh2) - mu(zc); ccf = mu(zhd(zc, vcf)) - mu(zc)
    Lc = np.linalg.cholesky(cov[np.ix_(idx, idx)]); w = lambda X: np.linalg.solve(Lc, X)
    ra, dec = np.radians(d.RA.values[idx]), np.radians(d.DEC.values[idx])
    nh = np.stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)], -1)
    def fit(cols, yy=y):
        X = np.column_stack([np.ones(len(yy))] + cols); Xw, yw = w(X), w(yy)
        b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); rr = yw - Xw@b
        return b, np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))), rr@rr
    _, _, c0 = fit([])
    print(f"chi2 no correction        = {c0:.1f}")
    for name, cc in [('2M++', c2), ('CF4++', ccf)]:
        _, _, cf = fit([], y - cc)
        b, e, cb = fit([cc])
        _, _, cd = fit([cc] + [nh[:, k] for k in range(3)])
        print(f"{name:6s} beta=1: chi2={cf:.1f} (dchi2={c0-cf:+.1f}) | free beta={b[1]:.2f}±{e[1]:.2f}: chi2={cb:.1f} | +dipole: dchi2={cb-cd:.1f} (3 dof, p={chi2d.sf(cb-cd,3):.3f})")
    b, e, cj = fit([c2, ccf])
    print(f"joint: beta_2M++={b[1]:.2f}±{e[1]:.2f}, beta_CF4++={b[2]:.2f}±{e[2]:.2f}, chi2={cj:.1f}")
    # which one does the data prefer, at beta=1: direct chi2 difference, and in shells
    print("shells (beta=1 chi2 contributions are not separable with a full cov; so per-shell free betas instead):")
    shells = [(0.01, 0.02), (0.02, 0.03), (0.03, 0.045), (0.045, 0.065), (0.065, zmax)]
    for name, cc in [('2M++', c2), ('CF4++', ccf)]:
        cols = [np.where((zc >= lo) & (zc < hi), cc, 0.0) for lo, hi in shells]
        b, e, _ = fit(cols)
        print(f"  {name:6s} " + "  ".join(f"z[{lo:.3f},{hi:.3f}) {bi:5.2f}±{ei:4.2f}" for (lo, hi), bi, ei in zip(shells, b[1:], e[1:])))

main(*sys.argv[1:4])
