"""Convention scan: for each axis permutation and distance scale of the CF4++ grid, the correlation of its radial velocity
at Pantheon+ positions with Pantheon+'s 2M++ VPEC, and the beta the Hubble residuals want. The right convention should
maximise the correlation; if none gets beta near 1, the field's amplitude really is off.
usage: python cf4_variants.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/new_CF4pp_mean_std_grids.npz
"""
import sys, itertools, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.ndimage import map_coordinates
from astropy.coordinates import SkyCoord
import astropy.units as u
C = 299792.458; OM = 0.334; N = 128; L = 1000.0; DELTA = L/N
def dc_h(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3 + 1-OM)
    return C/100*np.array([quad(e, 0, zi)[0] for zi in np.atleast_1d(z)])
def mu(z): return 5*np.log10((1+z)*dc_h(z)/0.7) + 25
def zhd(zcmb, v): return (1+zcmb)/(1+v/C) - 1

def main(dat, covf, grid, zmax=0.1):
    g = np.load(grid); V = g['v_mean']
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f:
        n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    m = (d.IS_CALIBRATOR == 0) & (d.zCMB > 0.01) & (d.zCMB < zmax); idx = np.where(m)[0]
    zc = d.zCMB.values[idx]; vp2 = d.VPEC.values[idx]; r = dc_h(zc)
    sc = SkyCoord(ra=d.RA.values[idx]*u.deg, dec=d.DEC.values[idx]*u.deg, distance=r*u.Mpc, frame='fk5', equinox='J2000').supergalactic
    xyz0 = np.stack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value])
    y = d.m_b_corr.values[idx] - mu(zc)
    Lc = np.linalg.cholesky(cov[np.ix_(idx, idx)]); w = lambda X: np.linalg.solve(Lc, X)
    def beta(v):
        cc = mu(zhd(zc, v)) - mu(zc)
        X = np.column_stack([np.ones(len(y)), cc]); Xw, yw = w(X), w(y)
        b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); e = np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw)))
        return b[1], e[1]
    print("perm(pos axes -> grid axes)  scale   corr(vr,VPEC)   rms vr   beta")
    rows = []
    for perm in itertools.permutations(range(3)):
        for scale in [0.7, 1.0, 1/0.7]:
            xyz = xyz0[list(perm)]*scale; rhat = xyz/np.linalg.norm(xyz, axis=0)
            ci = (xyz + L/2)/DELTA - 0.5
            v = np.stack([map_coordinates(V[k], ci, order=1, mode='nearest') for k in range(3)])
            vr = np.sum(v*rhat, 0)
            b, e = beta(vr); c = np.corrcoef(vr, vp2)[0, 1]
            rows.append((c, perm, scale, np.sqrt(np.mean(vr**2)), b, e))
    for c, perm, scale, rms, b, e in sorted(rows, reverse=True):
        print(f"  {perm}  {scale:5.3f}   {c:+.3f}   {rms:5.0f}   {b:5.2f}±{e:4.2f}")
    # the well-constrained subset under the default convention
    xyz = xyz0; rhat = xyz/np.linalg.norm(xyz, axis=0); ci = (xyz + L/2)/DELTA - 0.5
    vr = np.sum(np.stack([map_coordinates(V[k], ci, order=1, mode='nearest') for k in range(3)])*rhat, 0)
    vs = map_coordinates(g['vr_std'], ci, order=1, mode='nearest')
    print("\nby CF4++ vr_std tercile (default convention): n, mean std, rms vr, corr with VPEC, beta")
    for lo, hi in [(0, 33), (33, 67), (67, 100)]:
        s = (vs >= np.percentile(vs, lo)) & (vs <= np.percentile(vs, hi))
        vv = np.where(s, vr, 0.0); b, e = beta(vv)
        print(f"  std[{lo},{hi}]%: n={s.sum()} std={vs[s].mean():.0f} rms={np.sqrt(np.mean(vr[s]**2)):.0f} corr={np.corrcoef(vr[s], vp2[s])[0,1]:+.3f} beta={b:.2f}±{e:.2f}")

main(*sys.argv[1:4])
