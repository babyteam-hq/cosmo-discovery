"""Sweep rows 1-4: Pantheon+ Hubble residuals vs Milky-Way dust-map quantities, and residual-scatter anisotropy.
Full STAT+SYS covariance GLS, calibrators out, zHD>0.01, flat LCDM Om=0.334 (as frames.py). Tests (each a 1-dof slope
on top of an intercept, unless noted):
  MWEBV            Pantheon+'s own E(B-V) (0.86 x SFD): a residual slope means the extinction scaling is off
  SFD-CSFD         extragalactic (CIB/LSS) contamination of SFD at the SN position (Chiang 2023)
  GNILC-SFD        Planck-GNILC minus SFD E(B-V): temperature-dependent SFD bias
  SFD local excess SFD at the SN minus the median in a 20-40 arcmin annulus: unresolved emission at the SN (host)
  scatter          variance of whitened residuals by Galactic hemisphere / |b| / E(B-V) (Levene + Breusch-Pagan)
usage: python dust_maps.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/dustmaps
"""
import sys, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.stats import chi2 as chi2d, levene, norm
C = 299792.458; OM = 0.334

def dl(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3+1-OM); return (1+z)*C/70*np.array([quad(e, 0, zi)[0] for zi in z])

def maps(ra, dec, data_dir):
    from dustmaps.config import config; config['data_dir'] = data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.csfd import CSFDQuery; from dustmaps.planck import PlanckGNILCQuery
    from astropy.coordinates import SkyCoord; import astropy.units as u
    c = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    sfd = SFDQuery(); s = sfd(c); cs = CSFDQuery()(c); gn = PlanckGNILCQuery()(c)
    # annulus 20-40 arcmin: 24 points on each of two rings, median
    ann = []
    for r_arcmin in (20, 30, 40):
        for th in np.linspace(0, 2*np.pi, 24, endpoint=False):
            off = c.directional_offset_by(th*u.rad, r_arcmin*u.arcmin); ann.append(sfd(off))
    loc = s - np.median(np.array(ann), axis=0)
    b = c.galactic.b.deg
    return s, cs, gn, loc, b

def main(dat, covf, data_dir):
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f: n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    m = (d.IS_CALIBRATOR == 0) & (d.zHD > 0.01); idx = np.where(m)[0]
    z = d.zHD.values[idx]; y = d.m_b_corr.values[idx] - (5*np.log10(dl(z))+25)
    s, cs, gn, loc, b = maps(d.RA.values[idx], d.DEC.values[idx], data_dir)
    L = np.linalg.cholesky(cov[np.ix_(idx, idx)]); w = lambda X: np.linalg.solve(L, X)
    def fit(cols):
        X = np.column_stack([np.ones(len(z))]+cols); Xw, yw = w(X), w(y)
        bb, *_ = np.linalg.lstsq(Xw, yw, rcond=None); r = yw-Xw@bb
        return bb, np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))), r@r, r
    b0, e0, c0, r0 = fit([])
    print(f"N={len(z)} chi2_0={c0:.1f} (dof {len(z)-1})")
    tests = [('MWEBV (0.86 SFD)', d.MWEBV.values[idx]), ('SFD-CSFD', s-cs), ('GNILC-SFD', gn-s),
             ('SFD local excess (20-40\')', loc), ('SFD', s), ('CSFD', cs), ('GNILC', gn)]
    ntests = 0; rows = []
    for name, x in tests:
        x = np.asarray(x, float); ok = np.isfinite(x)
        if not ok.all(): print(f"  {name}: {(~ok).sum()} non-finite, filled with median"); x = np.where(ok, x, np.nanmedian(x))
        bb, ee, cc, _ = fit([x]); dc = c0-cc; p = chi2d.sf(dc, 1); ntests += 1
        rows.append((name, bb[1], ee[1], dc, p))
        print(f"  {name:26s} slope={bb[1]:8.3f} ± {ee[1]:6.3f} mag/E(B-V)  dchi2={dc:6.2f} p={p:.3g}  "
              f"x: median {np.median(x):.4f} std {np.std(x):.4f}")
    # the two map-difference tests together (2 dof)
    bb, ee, cc, _ = fit([s-cs, gn-s]); print(f"  SFD-CSFD + GNILC-SFD jointly: dchi2={c0-cc:.2f} dof=2 p={chi2d.sf(c0-cc,2):.3g}")
    # scatter anisotropy on the whitened residuals of the intercept-only fit
    print("scatter of whitened residuals (var ~1 if the covariance is right):")
    splits = [('Galactic b>0 vs b<0', b > 0), ('|b|>30 vs |b|<30', np.abs(b) > 30), ('Dec>0 vs Dec<0', d.DEC.values[idx] > 0),
              ('MWEBV above vs below median', d.MWEBV.values[idx] > np.median(d.MWEBV.values[idx]))]
    for name, sel in splits:
        a, bq = r0[sel], r0[~sel]; W, p = levene(a, bq, center='median'); ntests += 1
        print(f"  {name:30s} var {a.var():.3f} (n={sel.sum()}) vs {bq.var():.3f} (n={(~sel).sum()})  Levene p={p:.3g}")
    # Breusch-Pagan: r^2 on E(B-V) and |b|
    for name, x in [('MWEBV', d.MWEBV.values[idx]), ('|b|', np.abs(b)), ('SFD-CSFD', s-cs)]:
        X = np.column_stack([np.ones(len(x)), x]); u2 = r0**2
        bb, *_ = np.linalg.lstsq(X, u2, rcond=None); ss_res = ((u2-X@bb)**2).sum(); ss_tot = ((u2-u2.mean())**2).sum()
        lm = len(x)*(1-ss_res/ss_tot); p = chi2d.sf(lm, 1); ntests += 1
        print(f"  Breusch-Pagan r^2 vs {name:10s} LM={lm:6.2f} p={p:.3g}")
    pmin = min(r[4] for r in rows); print(f"trials: {ntests} tests; smallest slope p={pmin:.3g} -> Bonferroni p={min(1, pmin*ntests):.3g}")

main(*sys.argv[1:4])
