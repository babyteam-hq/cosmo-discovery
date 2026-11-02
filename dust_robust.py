"""Robustness of the dust_maps.py result (Pantheon+ residuals vs Planck-GNILC minus SFD E(B-V)): survey fixed effects,
jackknife by survey, sky/redshift/E(B-V) cuts, clipping of the map difference, alternative Planck maps (2013 tau353 and
radiance E(B-V), if data/dustmaps/planck holds the 2013 model map), and the two-map regression (y on SFD and GNILC
together: if GNILC were the truth, b_GNILC ~ +1 and b_SFD ~ -0.86).
usage: python dust_robust.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/dustmaps
"""
import sys, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.stats import chi2 as chi2d
C = 299792.458; OM = 0.334

def dl(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3+1-OM); return (1+z)*C/70*np.array([quad(e, 0, zi)[0] for zi in z])

def main(dat, covf, data_dir):
    from dustmaps.config import config; config['data_dir'] = data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.planck import PlanckGNILCQuery, PlanckQuery
    from astropy.coordinates import SkyCoord; import astropy.units as u
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f: n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    m = (d.IS_CALIBRATOR == 0) & (d.zHD > 0.01); idx0 = np.where(m)[0]
    c = SkyCoord(d.RA.values[idx0]*u.deg, d.DEC.values[idx0]*u.deg, frame='icrs')
    s = SFDQuery()(c); gn = PlanckGNILCQuery()(c); b = c.galactic.b.deg
    alt = {}
    try:
        alt['P13 tau353 E(B-V)'] = PlanckQuery(component='ebv')(c); alt['P13 radiance E(B-V)'] = PlanckQuery(component='radiance')(c)*1.0
    except Exception as e: print("no Planck 2013 map:", repr(e)[:120])
    z0 = d.zHD.values[idx0]; y0 = d.m_b_corr.values[idx0] - (5*np.log10(dl(z0))+25)
    sv = d.IDSURVEY.values[idx0]; mwebv = d.MWEBV.values[idx0]; x0 = gn - s
    def fit(sel, cols, fixed=False, label=""):
        ii = idx0[sel]; y = y0[sel]
        base = [np.ones(sel.sum())]
        if fixed:
            for sid in np.unique(sv[sel])[1:]: base.append((sv[sel] == sid).astype(float))
        L = np.linalg.cholesky(cov[np.ix_(ii, ii)]); w = lambda X: np.linalg.solve(L, X)
        def gls(cs):
            X = np.column_stack(base+cs); Xw, yw = w(X), w(y); bb, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
            r = yw-Xw@bb; return bb, np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))), r@r
        _, _, c0 = gls([]); bb, ee, c1 = gls([cc[sel] for cc in cols]); k = len(base)
        dc = c0-c1; print(f"  {label:44s} n={sel.sum():4d} " + " ".join(f"b{j}={bb[k+j]:6.2f}±{ee[k+j]:4.2f}" for j in range(len(cols)))
              + f"  dchi2={dc:6.2f} p={chi2d.sf(dc, len(cols)):.2g}")
        return bb[k:], ee[k:]
    allsel = np.ones(len(idx0), bool)
    print("baseline and survey fixed effects:")
    fit(allsel, [x0], label="GNILC-SFD")
    fit(allsel, [x0], fixed=True, label="GNILC-SFD + survey fixed effects")
    fit(allsel, [s, gn], label="y on SFD, GNILC (b0=SFD, b1=GNILC)")
    fit(allsel, [s, gn], fixed=True, label="  same + survey fixed effects")
    fit(allsel, [x0, mwebv], label="GNILC-SFD with MWEBV also (b0=diff)")
    fit(allsel, [x0, np.abs(b)/90], label="GNILC-SFD with |b|/90 also (b0=diff)")
    fit(allsel, [x0, z0], label="GNILC-SFD with z also (b0=diff)")
    print("cuts:")
    for name, sel in [("|b|>20", np.abs(b) > 20), ("|b|>30", np.abs(b) > 30), ("|b|<30", np.abs(b) < 30), ("MWEBV<0.1", mwebv < 0.1),
                      ("MWEBV<0.05", mwebv < 0.05), ("|diff|<0.02", np.abs(x0) < 0.02), ("|diff|<0.01", np.abs(x0) < 0.01),
                      ("z<0.1", z0 < 0.1), ("z>=0.1", z0 >= 0.1), ("Dec>0", d.DEC.values[idx0] > 0), ("Dec<0", d.DEC.values[idx0] < 0),
                      ("b>0", b > 0), ("b<0", b < 0)]:
        fit(sel, [x0], label=name)
    print("jackknife by survey (IDSURVEY, >=30 SNe):")
    bs = []
    for sid, cnt in pd.Series(sv).value_counts().items():
        if cnt >= 30:
            bb, ee = fit(sv != sid, [x0], label=f"without survey {sid} ({cnt})"); bs.append(bb[0])
    bs = np.array(bs); print(f"  jackknife spread of the slope over surveys: {bs.min():.2f} to {bs.max():.2f}")
    if alt:
        print("alternative Planck E(B-V) maps minus SFD:")
        for name, mp in alt.items():
            fit(allsel, [mp - s], label=f"{name} - SFD")
            fit(allsel, [mp - s], fixed=True, label=f"{name} - SFD + survey fixed effects")
        fit(allsel, [gn - alt['P13 tau353 E(B-V)']], label="GNILC - P13 tau353 (both Planck)")

main(*sys.argv[1:4])
