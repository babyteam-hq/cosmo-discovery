"""CF4 follow-up 1: the cf4_robust.py test (beta for the CF4++ field, its dipole and its non-dipole rest, and 2M++ for
comparison) on the Pantheon+ SNe whose hosts are NOT Cosmicflows-4 distance inputs, vs those that are. CF4 inputs are
the table-2 galaxies with a SN Ia distance modulus (DMsnIa; Tully+2023, VizieR J/ApJ/944/94); a Pantheon+ SN is "in
CF4" if its host position (HOST_RA/DEC, else the SN position) lies within `sep` arcsec of such a galaxy.
usage: python cf4_notin.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/new_CF4pp_mean_std_grids.npz data/cf4_table2.csv [sep=30]
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

def main(dat, covf, grid, cf4csv, sep=30.0):
    sep = float(sep)
    g = np.load(grid); V = g['v_mean']
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f:
        n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    t = pd.read_csv(cf4csv); t = t[t.DMsnIa.notna()]
    cf = SkyCoord(t.RAJ2000.values*u.deg, t.DEJ2000.values*u.deg)
    hra = np.where(d.HOST_RA.values > -1, d.HOST_RA.values, d.RA.values); hdec = np.where(d.HOST_RA.values > -1, d.HOST_DEC.values, d.DEC.values)
    _, s1, _ = SkyCoord(hra*u.deg, hdec*u.deg).match_to_catalog_sky(cf)
    _, s2, _ = SkyCoord(d.RA.values*u.deg, d.DEC.values*u.deg).match_to_catalog_sky(cf)
    incf4 = (s1.arcsec < sep) | (s2.arcsec < sep)
    base = (d.IS_CALIBRATOR == 0) & (d.zCMB > 0.01) & (d.zCMB < 0.1)
    print(f"Pantheon+ SNe 0.01<z<0.1: {base.sum()}; hosts within {sep:.0f}\" of a CF4 SN-Ia-distance galaxy: {(base & incf4).sum()}, not: {(base & ~incf4).sum()}"
          f" (CF4 galaxies with DMsnIa: {len(t)}; matched by any Pantheon+ entry: {len(set(s1[incf4]) | set(s2[incf4]))})")
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
    run(base & incf4, "in CF4 (host matched)")
    run(base & ~incf4, "NOT in CF4")
    run(base & ~incf4 & (d.zCMB < 0.05), "NOT in CF4, z<0.05")
    run(base & ~incf4 & (d.zCMB >= 0.05), "NOT in CF4, z>=0.05")
    print("jackknife by survey on the NOT-in-CF4 set (>=20 SNe):")
    for sid, cnt in d[base & ~incf4].IDSURVEY.value_counts().items():
        if cnt >= 20: run(base & ~incf4 & (d.IDSURVEY != sid), f"  without survey {sid} ({cnt})")

main(*sys.argv[1:6])
