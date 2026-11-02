"""Robustness of the CF4++ result: beta for the full field, its dipole and its non-dipole rest, under jackknife by survey,
redshift cuts, and the same split applied to 2M++ (to show the test has power to accept non-dipole structure).
usage: python cf4_robust.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/new_CF4pp_mean_std_grids.npz
"""
import sys, numpy as np, pandas as pd
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

def main(dat, covf, grid):
    g = np.load(grid); V = g['v_mean']
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f:
        n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    base = (d.IS_CALIBRATOR == 0) & (d.zCMB > 0.01) & (d.zCMB < 0.1)
    def run(mask, label):
        idx = np.where(mask)[0]
        zc = d.zCMB.values[idx]; vp2 = d.VPEC.values[idx]; r = dc_h(zc)
        sc = SkyCoord(ra=d.RA.values[idx]*u.deg, dec=d.DEC.values[idx]*u.deg, distance=r*u.Mpc, frame='fk5', equinox='J2000').supergalactic
        xyz = np.stack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value]); rhat = xyz/np.linalg.norm(xyz, axis=0)
        ci = (xyz + L/2)/DELTA - 0.5
        vcf = np.sum(np.stack([map_coordinates(V[k], ci, order=1, mode='nearest') for k in range(3)])*rhat, 0)
        y = d.m_b_corr.values[idx] - mu(zc)
        Lc = np.linalg.cholesky(cov[np.ix_(idx, idx)]); w = lambda X: np.linalg.solve(Lc, X)
        def beta(v):
            cc = mu(zhd(zc, v)) - mu(zc)
            X = np.column_stack([np.ones(len(y)), cc]); Xw, yw = w(X), w(y)
            b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); e = np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw)))
            return f"{b[1]:5.2f}±{e[1]:4.2f}"
        A = np.column_stack([r, rhat.T])
        def split(v):
            coef, *_ = np.linalg.lstsq(A, v, rcond=None); return A@np.r_[0, coef[1:]], v - A@coef
        dcf, rcf = split(vcf); d2, r2 = split(vp2)
        print(f"{label:34s} n={len(idx):4d} | CF4++ full {beta(vcf)} dipole {beta(dcf)} rest {beta(rcf)} (rms rest {np.sqrt(np.mean(rcf**2)):3.0f}) | 2M++ full {beta(vp2)} dipole {beta(d2)} rest {beta(r2)} (rms rest {np.sqrt(np.mean(r2**2)):3.0f})")
    run(base, "all 0.01<z<0.1")
    run(base & (d.zCMB < 0.05), "z<0.05")
    run(base & (d.zCMB >= 0.05), "z>=0.05")
    run(base & (d.zCMB < 0.03), "z<0.03")
    run(base & (d.zCMB >= 0.03) & (d.zCMB < 0.06), "0.03<=z<0.06")
    print("jackknife: leave one survey out (IDSURVEY)")
    for sid, cnt in d[base].IDSURVEY.value_counts().items():
        if cnt >= 30:
            run(base & (d.IDSURVEY != sid), f"  without survey {sid} ({cnt} SNe)")

main(*sys.argv[1:4])
