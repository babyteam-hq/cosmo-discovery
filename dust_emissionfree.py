"""Does the GNILC-SFD residual slope (results.md row 3, finding.md) mean SFD is wrong where GNILC differs? Test with an
emission-free reddening map: Edenhofer+2023 (Gaia/2MASS stellar extinction, 3D to 1.25 kpc, integrated here to its far
edge, Zenodo; dustmaps 'edenhofer2023'), and Lenz+2017 (HI-based E(B-V), if fetched). Each map is first put on SFD's
scale by a least-squares fit MAP ~ k x SFD over the sample (so the difference has no SFD-scaling term), then the
Hubble residual is regressed on (k MAP - SFD) and on (k MAP - SFD, GNILC - SFD) together, for Pantheon+ (full
covariance) and ZTF DR2 (diagonal errors, alpha/beta free, as ztf_tests.py). Prediction if GNILC is right and SFD wrong:
a positive slope on (k MAP - SFD) of the GNILC size, and (k MAP - SFD) takes over from (GNILC - SFD) in the joint fit.
usage: python dust_emissionfree.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/ztf/ztfsniadr2_lite/tables/snia_data.csv data/dustmaps
"""
import sys, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.stats import chi2 as chi2d
from astropy.coordinates import SkyCoord
import astropy.units as u
C = 299792.458; OM = 0.334; AL, BE = 0.15, 3.1

def dl(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3+1-OM); return (1+z)*C/70*np.array([quad(e, 0, zi)[0] for zi in z])

def maps(c, data_dir):
    from dustmaps.config import config; config['data_dir'] = data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.planck import PlanckGNILCQuery
    out = {'SFD': SFDQuery()(c), 'GNILC': PlanckGNILCQuery()(c)}
    try:
        from dustmaps.edenhofer2023 import Edenhofer2023Query
        q = Edenhofer2023Query(integrated=True)
        far = SkyCoord(c.ra, c.dec, distance=1.24*u.kpc, frame='icrs'); out['Edenhofer'] = q(far)
    except Exception as e: print("no Edenhofer map:", repr(e)[:150])
    try:
        from dustmaps.lenz2017 import Lenz2017Query; out['Lenz'] = Lenz2017Query()(c)
    except Exception as e: print("no Lenz map:", repr(e)[:150])
    return out

def report(y, w, base, mp, bdeg, label):
    """w: whitening function on column stacks; base: list of base columns."""
    s, g = mp['SFD'], mp['GNILC']
    def gls(cols, sel):
        X = np.column_stack([b[sel] for b in base] + [cc[sel] for cc in cols]); Xw, yw = w(X, sel), w(y[sel], sel)
        b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); e = np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))); r = yw-Xw@b; return b, e, r@r
    k0 = len(base)
    for name in [n for n in ('Edenhofer', 'Lenz') if n in mp]:
        m = mp[name]; ok = np.isfinite(m) & (m > -0.5)
        for bl, sel in [("all", ok), ("|b|>20", ok & (bdeg > 20)), ("|b|>30", ok & (bdeg > 30))]:
            if sel.sum() < 50: continue
            k = np.sum(s[sel]*m[sel])/np.sum(m[sel]**2); x = k*m - s
            _, _, c0 = gls([], sel); b1, e1, c1 = gls([x], sel); b2, e2, c2 = gls([x, g - s], sel); b3, e3, c3 = gls([g - s], sel)
            print(f"  {label} {name:9s} {bl:7s} n={sel.sum():4d} k={k:.3f} std(kMAP-SFD)={np.std(x[sel]):.4f} | slope kMAP-SFD {b1[k0]:5.2f}±{e1[k0]:4.2f} (p={chi2d.sf(c0-c1,1):.2g})"
                  f" | GNILC-SFD alone {b3[k0]:5.2f}±{e3[k0]:4.2f} | joint: kMAP-SFD {b2[k0]:5.2f}±{e2[k0]:4.2f}, GNILC-SFD {b2[k0+1]:5.2f}±{e2[k0+1]:4.2f}"
                  f" | corr(kMAP-SFD, GNILC-SFD)={np.corrcoef(x[sel], (g-s)[sel])[0,1]:.2f}")

def main(dat, covf, ztfcsv, data_dir):
    # Pantheon+
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f: n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    msk = (d.IS_CALIBRATOR == 0) & (d.zHD > 0.01); idx = np.where(msk)[0]
    z = d.zHD.values[idx]; y = d.m_b_corr.values[idx] - (5*np.log10(dl(z))+25)
    c = SkyCoord(d.RA.values[idx]*u.deg, d.DEC.values[idx]*u.deg, frame='icrs'); mp = maps(c, data_dir)
    Lfull = cov[np.ix_(idx, idx)]
    def w(X, sel):
        Lc = np.linalg.cholesky(Lfull[np.ix_(sel, sel)]); return np.linalg.solve(Lc, X)
    print("Pantheon+ (full covariance):")
    report(y, w, [np.ones(len(y))], mp, np.abs(c.galactic.b.deg), "P+ ")
    # ZTF
    t = pd.read_csv(ztfcsv)
    sel = (t.sn_type == 'snia-cosmo') & (t.lccoverage_flag == 1) & (t.fitquality_flag == 1) & (t.source == 'z_gal') \
        & (t.x1.abs() < 3) & (t.c.abs() < 0.3) & (t.redshift > 0.015) & (t.redshift < 0.1) & t.x0.notna()
    t = t[sel].reset_index(drop=True)
    ct = SkyCoord(t.ra.values*u.deg, t.dec.values*u.deg, frame='icrs'); apex = SkyCoord(l=264.021*u.deg, b=48.253*u.deg, frame='galactic')
    zc = (1+t.redshift.values)/(1 - 369.82*np.cos(ct.galactic.separation(apex).rad)/C) - 1
    mb = -2.5*np.log10(t.x0.values); x1 = t.x1.values; cc = t.c.values
    var = (2.5/np.log(10)*t.x0_err.values/t.x0.values)**2 + AL**2*t.x1_err.values**2 + BE**2*t.c_err.values**2 \
        + 2*(-AL)*(-2.5/np.log(10)/t.x0.values)*t.cov_x0_x1.values + 2*BE*(-2.5/np.log(10)/t.x0.values)*t.cov_x0_c.values \
        + 2*(-AL)*BE*t.cov_x1_c.values + (5/np.log(10)*300/(C*zc))**2 + (5/np.log(10)*t.redshift_err.values/zc)**2 + 0.078**2
    yz = mb - (5*np.log10(dl(zc))+25); wz = 1/np.sqrt(var)
    def w2(X, sel): return X*wz[sel][:, None] if X.ndim == 2 else X*wz[sel]
    mpz = maps(ct, data_dir)
    print("ZTF DR2 (diagonal errors, alpha and beta free):")
    report(yz, w2, [np.ones(len(yz)), -x1, cc], mpz, np.abs(ct.galactic.b.deg), "ZTF")

main(*sys.argv[1:5])
