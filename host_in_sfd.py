"""Does the SFD dust map contain the SN host galaxy's own far-IR emission? For every Pantheon+ SN (calibrators included),
E(B-V) at the SN position minus the mean on a ring of control points (radius R arcmin, default 30), for SFD, CSFD,
Planck-2013 tau353 and Planck GNILC. Host emission is unresolved at SFD's 6' beam and CIB-like, so the hypothesis
predicts a positive SN-minus-ring excess in SFD/CSFD/tau353 that grows toward low redshift and with host mass, and
none (or less) in GNILC, which separates the CIB. Also the excess vs the Pantheon+ Hubble residual with its expected
coefficient. usage: python host_in_sfd.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/dustmaps [R_arcmin]
"""
import sys, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.stats import chi2 as chi2d, ttest_1samp, wilcoxon
C = 299792.458; OM = 0.334

def dl(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3+1-OM); return (1+z)*C/70*np.array([quad(e, 0, zi)[0] for zi in z])

def main(dat, covf, data_dir, R=30.0):
    R = float(R)
    from dustmaps.config import config; config['data_dir'] = data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.csfd import CSFDQuery; from dustmaps.planck import PlanckGNILCQuery, PlanckQuery
    from astropy.coordinates import SkyCoord; import astropy.units as u
    d = pd.read_csv(dat, sep=r'\s+')
    c = SkyCoord(d.RA.values*u.deg, d.DEC.values*u.deg, frame='icrs')
    qs = {'SFD': SFDQuery(), 'CSFD': CSFDQuery(), 'P13 tau353': PlanckQuery(component='ebv'), 'GNILC': PlanckGNILCQuery()}
    rings = [c.directional_offset_by(th*u.rad, R*u.arcmin) for th in np.linspace(0, 2*np.pi, 16, endpoint=False)]
    exc = {}
    for name, q in qs.items():
        at = q(c); ring = np.mean([q(r) for r in rings], axis=0); exc[name] = at - ring
    z = d.zHD.values; lm = d.HOST_LOGMASS.values
    bins = [(0, 0.01), (0.01, 0.03), (0.03, 0.06), (0.06, 0.1), (0.1, 0.3), (0.3, 3)]
    print(f"E(B-V) at the SN minus the mean of 16 points {R:.0f}' away: mean ± error (x1000), Wilcoxon p, by zHD bin")
    for lo, hi in bins:
        s = (z >= lo) & (z < hi); line = f"  {lo:<5g}<=z<{hi:<4g} n={s.sum():4d}"
        for name in qs:
            x = exc[name][s]; x = x[np.isfinite(x)]
            line += f" | {name} {1000*x.mean():6.2f}±{1000*x.std()/np.sqrt(len(x)):4.2f} p={wilcoxon(x).pvalue:.1e}"
        print(line)
    print("by host log mass (0.01<z<0.1):")
    for lo, hi in [(0, 9.5), (9.5, 10.5), (10.5, 13)]:
        s = (z > 0.01) & (z < 0.1) & (lm >= lo) & (lm < hi); line = f"  {lo:<4g}<=logM<{hi:<4g} n={s.sum():4d}"
        for name in qs:
            x = exc[name][s]; x = x[np.isfinite(x)]
            line += f" | {name} {1000*x.mean():6.2f}±{1000*x.std()/np.sqrt(len(x)):4.2f}"
        print(line)
    # the residual test: GLS of the Hubble residual on the SFD excess, and on (GNILC - SFD) with the excess as a covariate
    with open(covf) as f: n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    m = (d.IS_CALIBRATOR == 0) & (d.zHD > 0.01); idx = np.where(m)[0]
    zz = z[idx]; y = d.m_b_corr.values[idx] - (5*np.log10(dl(zz))+25)
    L = np.linalg.cholesky(cov[np.ix_(idx, idx)]); w = lambda X: np.linalg.solve(L, X)
    def fit(cols, label):
        X = np.column_stack([np.ones(len(y))]+cols); Xw, yw = w(X), w(y); bb, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        ee = np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))); r = yw-Xw@bb
        X0 = np.ones((len(y), 1)); b0, *_ = np.linalg.lstsq(w(X0), yw, rcond=None); r0 = yw - w(X0)@b0
        dc = r0@r0 - r@r
        print(f"  {label:50s} " + " ".join(f"b{j}={bb[1+j]:6.2f}±{ee[1+j]:4.2f}" for j in range(len(cols))) + f" dchi2={dc:5.2f} p={chi2d.sf(dc, len(cols)):.2g}")
    gs = qs['GNILC'](c)[idx] - qs['SFD'](c)[idx]
    print("Hubble residual (0.01<z, calibrators out) GLS, full covariance:")
    for name in qs: fit([np.nan_to_num(exc[name][idx])], f"residual vs {name} excess at the SN")
    fit([gs], "residual vs GNILC-SFD")
    fit([gs, np.nan_to_num(exc['SFD'][idx])], "residual vs GNILC-SFD (b0) and SFD excess (b1)")
    fit([gs, np.nan_to_num(exc['P13 tau353'][idx])], "residual vs GNILC-SFD (b0) and tau353 excess (b1)")
    hm = np.where(np.isfinite(lm[idx]) & (lm[idx] > 0), lm[idx], 10.0)
    fit([gs, (hm > 10).astype(float)], "residual vs GNILC-SFD (b0) and host mass step (b1)")
    fit([gs, hm - 10], "residual vs GNILC-SFD (b0) and host logM-10 (b1)")

main(*sys.argv[1:5])
