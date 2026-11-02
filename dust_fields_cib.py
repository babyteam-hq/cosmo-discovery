"""Two checks the third blind referee asked of the dust finding (finding.md, "Third blind review"), for Pantheon+ (full
STAT+SYS covariance, calibrators out, zHD > 0.01, as dust_rotation_null.py) and ZTF DR2 (ztf_tests.py's sample and diagonal
errors, alpha and beta free, sigma_int for chi2/dof = 1).

1. Per-field residuals. The high-z Pantheon+ SNe sit in a few deep fields, so per-SN counts overstate the independent sky
   samples. Cells = survey (IDSURVEY) x HEALPix nside-32 pixel (1.8 deg) of the SN position.
   (a) the GNILC - SFD template replaced by its cell mean, per-SN full-covariance GLS slope (signal from between-cell variation);
   (b) block average: y, template and covariance (A C A^T, A the cell-averaging matrix), GLS on the cells;
   (c) leave-one-cell-out and leave-one-survey-out on the per-cell fit (min and max slope);
   (d) the per-cell fit for z < 0.1 (the all-sky low-z surveys) and z >= 0.1 (the deep fields).
   ZTF: (b) only (one survey; cells = nside-32 pixels; weighted cell means, cell variance 1/sum w).
2. Which map difference carries it. GNILC differs from the Planck 2013 tau353 model in its CIB removal and temperature model.
   Residuals regressed jointly on (GNILC - tau353, tau353 - SFD) and on (GNILC - radiance, radiance - SFD), with the regressor
   correlation at the SN positions; the slope on GNILC - SFD with the Planck 2013 dust temperature T as a second covariate, and
   on T alone. Planck 2013 components from dustmaps' PlanckQuery (tau x 1.49e4, radiance x 5.4e5, temp).
usage: python dust_fields_cib.py data/pantheon/Pantheon+SH0ES.dat data/pantheon/Pantheon+SH0ES_STAT+SYS.cov
       data/ztf/ztfsniadr2_lite/tables/snia_data.csv data/dustmaps
"""
import sys, numpy as np, pandas as pd, healpy as hp
from scipy.integrate import quad
from astropy.coordinates import SkyCoord
import astropy.units as u
C = 299792.458; NSIDE = 32

def dl(z, om):
    e = lambda x: 1/np.sqrt(om*(1+x)**3 + 1 - om); return (1+z)*C/70*np.array([quad(e, 0, zi)[0] for zi in z])

def cells_of(ra, dec, survey):
    pix = hp.ang2pix(NSIDE, ra, dec, lonlat=True); key = np.array([f'{s}:{p}' for s, p in zip(survey, pix)])
    u_, inv = np.unique(key, return_inverse=True); A = np.zeros((len(u_), len(ra)))
    A[inv, np.arange(len(ra))] = 1; A /= A.sum(1, keepdims=True); return A, inv, u_

def gls(Xcols, y, Ci):
    X = np.column_stack(Xcols); Aι = X.T@Ci@X; b = np.linalg.solve(Aι, X.T@Ci@y); e = np.sqrt(np.diag(np.linalg.inv(Aι))); return b, e

def maps_at(c, data_dir):
    from dustmaps.config import config; config['data_dir'] = data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.planck import PlanckGNILCQuery, PlanckQuery
    out = dict(sfd=SFDQuery()(c), gnilc=PlanckGNILCQuery()(c))
    try:
        out['tau'] = PlanckQuery(component='tau')(c); out['rad'] = PlanckQuery(component='radiance')(c); out['temp'] = np.asarray(PlanckQuery(component='temp')(c), float)
        # R-115/R-116 took dustmaps' 'ebv' component as "tau353 E(B-V)"; it is the file's EBV column, which tracks the radiance estimator
        ebv = PlanckQuery(component='ebv')(c)
        print(f"  Planck 2013 EBV column vs estimators: median EBV/(1.49e4 tau) {np.median(ebv/out['tau']):.2f}, corr {np.corrcoef(ebv, out['tau'])[0,1]:.3f}; median EBV/(5.4e5 R) {np.median(ebv/out['rad']):.2f}, corr {np.corrcoef(ebv, out['rad'])[0,1]:.3f}; "
              f"corr(EBV − SFD, GNILC − SFD) {np.corrcoef(ebv - out['sfd'], out['gnilc'] - out['sfd'])[0,1]:+.2f}, corr(1.49e4 tau − SFD, GNILC − SFD) {np.corrcoef(out['tau'] - out['sfd'], out['gnilc'] - out['sfd'])[0,1]:+.2f}")
    except Exception as e:
        print(f'  Planck 2013 model map not available ({str(e)[:80]}); carrier tests skipped'); out['tau'] = None
    return out

def report_cells(label, y, t, Ci_full, A, cov, z=None, surveys=None, cellkeys=None):
    n = len(y); one = np.ones(n)
    # (a) template replaced by its cell mean, per-SN GLS
    cell_of = np.argmax(A > 0, axis=0); tc = (A@t)[cell_of]  # each SN's cell mean
    b, e = gls([one, tc], y, Ci_full); print(f'  {label} (a) cell-mean template, per-SN GLS: slope {b[1]:+.3f} ± {e[1]:.3f} (z = {b[1]/e[1]:+.2f}), cells {A.shape[0]}')
    # (b) block average
    yc, tc2, Cc = A@y, A@t, A@cov@A.T; Cic = np.linalg.inv(Cc); onec = np.ones(len(yc))
    b, e = gls([onec, tc2], yc, Cic); print(f'  {label} (b) block-averaged GLS on {len(yc)} cells: slope {b[1]:+.3f} ± {e[1]:.3f} (z = {b[1]/e[1]:+.2f})')
    res = dict(a=(None), b=(b[1], e[1], len(yc)))
    if surveys is not None:
        # (c) leave-one-cell-out and leave-one-survey-out on the per-cell fit
        sl = []
        for k in range(len(yc)):
            m = np.arange(len(yc)) != k; Ci_ = np.linalg.inv(Cc[np.ix_(m, m)]); bb, ee = gls([onec[m], tc2[m]], yc[m], Ci_); sl.append(bb[1])
        print(f'  {label} (c) leave-one-cell-out: slope min {min(sl):+.3f}, max {max(sl):+.3f}')
        cs = np.array([k.split(":")[0] for k in cellkeys]); sl = {}
        for s_ in np.unique(cs):
            m = cs != s_
            if m.sum() < 5: continue
            Ci_ = np.linalg.inv(Cc[np.ix_(m, m)]); bb, ee = gls([onec[m], tc2[m]], yc[m], Ci_); sl[s_] = (bb[1], ee[1])
        kmin, kmax = min(sl, key=lambda k: sl[k][0]), max(sl, key=lambda k: sl[k][0])
        print(f'  {label} (c) leave-one-survey-out: min {sl[kmin][0]:+.3f} ± {sl[kmin][1]:.3f} (without survey {kmin}), max {sl[kmax][0]:+.3f} ± {sl[kmax][1]:.3f} (without {kmax})')
    if z is not None:
        zc = A@z
        for name, m in (('z < 0.1', zc < 0.1), ('z >= 0.1', zc >= 0.1)):
            Ci_ = np.linalg.inv(Cc[np.ix_(m, m)]); bb, ee = gls([onec[m], tc2[m]], yc[m], Ci_)
            print(f'  {label} (d) per-cell fit, {name}: {m.sum()} cells, slope {bb[1]:+.3f} ± {ee[1]:.3f} (z = {bb[1]/ee[1]:+.2f})')
    return res

def report_carrier(label, y, Ci, M):
    one = np.ones(len(y)); g, s, t, r, T = M['gnilc'], M['sfd'], M['tau'], M['rad'], M['temp']
    if t is None: return
    for name, x in (('tau353 − SFD', t - s), ('radiance − SFD', r - s), ('GNILC − tau353', g - t), ('GNILC − radiance', g - r), ('GNILC − SFD', g - s)):
        b, e = gls([one, x], y, Ci); print(f'  {label} marginal slope on {name:16s} {b[1]:+.3f} ± {e[1]:.3f} (z {b[1]/e[1]:+.2f}); rms {np.std(x):.4f}; corr with GNILC − SFD {np.corrcoef(x, g - s)[0, 1]:+.2f}')
    for name, x1, x2, n1, n2 in (('GNILC − tau353 | tau353 − SFD', g - t, t - s, 'GNILC−τ', 'τ−SFD'), ('GNILC − radiance | radiance − SFD', g - r, r - s, 'GNILC−R', 'R−SFD')):
        b, e = gls([one, x1, x2], y, Ci); rho = np.corrcoef(x1, x2)[0, 1]
        print(f'  {label} carrier: {n1} {b[1]:+.3f} ± {e[1]:.3f} (z {b[1]/e[1]:+.2f}), {n2} {b[2]:+.3f} ± {e[2]:.3f} (z {b[2]/e[2]:+.2f}); corr(regressors) {rho:+.2f}; rms {np.std(x1):.4f}, {np.std(x2):.4f}')
    Tm = T - np.nanmean(T)
    b, e = gls([one, g - s, Tm], y, Ci); print(f'  {label} GNILC − SFD with T: slope {b[1]:+.3f} ± {e[1]:.3f}, T {1e3*b[2]:+.2f} ± {1e3*e[2]:.2f} mmag/K; corr(GNILC−SFD, T) {np.corrcoef(g - s, Tm)[0, 1]:+.2f}')
    b, e = gls([one, Tm], y, Ci); print(f'  {label} T alone: {1e3*b[1]:+.2f} ± {1e3*e[1]:.2f} mmag/K (z {b[1]/e[1]:+.2f})')

def pantheon(dat, covf, data_dir):
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f: n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    m = (d.IS_CALIBRATOR == 0) & (d.zHD > 0.01); idx = np.where(m)[0]; cov = cov[np.ix_(idx, idx)]
    z = d.zHD.values[idx]; y = d.m_b_corr.values[idx] - (5*np.log10(dl(z, 0.334)) + 25); ra, dec, sid = d.RA.values[idx], d.DEC.values[idx], d.IDSURVEY.values[idx]
    c = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs'); M = maps_at(c, data_dir); t = M['gnilc'] - M['sfd']
    Ci = np.linalg.inv(cov); one = np.ones(len(y)); b, e = gls([one, t], y, Ci)
    print(f'Pantheon+: {len(y)} SNe; per-SN slope on GNILC − SFD {b[1]:+.3f} ± {e[1]:.3f} (z = {b[1]/e[1]:+.2f})')
    A, inv, keys = cells_of(ra, dec, sid); print(f'  cells (survey x nside-{NSIDE} pixel): {len(keys)}; SNe per cell median {np.median(np.bincount(inv))}, max {np.bincount(inv).max()}; {np.sum(np.bincount(inv) == 1)} singleton cells')
    report_cells('Pantheon+', y, t, Ci, A, cov, z=z, surveys=sid, cellkeys=keys)
    report_carrier('Pantheon+', y, Ci, M)

def ztf(csv, data_dir):
    AL, BE = 0.15, 3.1; d = pd.read_csv(csv)
    sel = (d.sn_type == 'snia-cosmo') & (d.lccoverage_flag == 1) & (d.fitquality_flag == 1) & (d.source == 'z_gal') \
        & (d.x1.abs() < 3) & (d.c.abs() < 0.3) & (d.redshift > 0.015) & (d.redshift < 0.1) & d.x0.notna()
    d = d[sel].reset_index(drop=True); c = SkyCoord(d.ra.values*u.deg, d.dec.values*u.deg, frame='icrs'); gal = c.galactic
    apex = SkyCoord(l=264.021*u.deg, b=48.253*u.deg, frame='galactic'); zcmb = (1 + d.redshift.values)/(1 - 369.82*np.cos(gal.separation(apex).rad)/C) - 1
    M = maps_at(c, data_dir); t = M['gnilc'] - M['sfd']
    mb = -2.5*np.log10(d.x0.values); x1 = d.x1.values; col = d.c.values
    var_m = (2.5/np.log(10)*d.x0_err.values/d.x0.values)**2 + AL**2*d.x1_err.values**2 + BE**2*d.c_err.values**2 \
        + 2*(-AL)*(-2.5/np.log(10)/d.x0.values)*d.cov_x0_x1.values + 2*BE*(-2.5/np.log(10)/d.x0.values)*d.cov_x0_c.values + 2*(-AL)*BE*d.cov_x1_c.values
    var_v = (5/np.log(10)*300/(C*zcmb))**2 + (5/np.log(10)*d.redshift_err.values/zcmb)**2
    e_ = lambda x: 1/np.sqrt(0.3*(1+x)**3 + 0.7); dch = C/100*np.array([quad(e_, 0, zi)[0] for zi in zcmb]); y = mb - (5*np.log10((1+zcmb)*dch/0.7) + 25)
    one = np.ones(len(y)); base = [one, -x1, col]
    def fitd(cols, s):
        w = 1/(var_m + var_v + s*s); X = np.column_stack(base + cols); Aι = X.T@(X*w[:, None]); b = np.linalg.solve(Aι, X.T@(w*y)); r = y - X@b; return b, np.sqrt(np.diag(np.linalg.inv(Aι))), np.sum(w*r*r)
    lo, hi = 0.0, 0.5
    for _ in range(40):
        s_ = (lo + hi)/2; _, _, chi = fitd([], s_)
        if chi > len(y) - 3: lo = s_
        else: hi = s_
    sig = (lo + hi)/2; var = var_m + var_v + sig*sig; Ci = np.diag(1/var)
    b, e, _ = fitd([t], sig); print(f'ZTF: {len(y)} SNe, sigma_int {sig:.3f}; per-SN slope on GNILC − SFD {b[3]:+.3f} ± {e[3]:.3f} (z = {b[3]/e[3]:+.2f})')
    # (b) block average with diagonal weights: weighted cell means of y, x1, c, t; cell variance 1/sum w
    w = 1/var; pix = hp.ang2pix(NSIDE, d.ra.values, d.dec.values, lonlat=True); u_, inv = np.unique(pix, return_inverse=True)
    W = np.zeros((len(u_), len(y))); W[inv, np.arange(len(y))] = w; sw = W.sum(1); Aw = W/sw[:, None]
    yc, x1c, cc, tc = Aw@y, Aw@x1, Aw@col, Aw@t; Cic = np.diag(sw); onec = np.ones(len(yc))
    bb, ee = gls([onec, -x1c, cc, tc], yc, Cic); print(f'  ZTF (b) block-averaged fit on {len(yc)} nside-{NSIDE} cells (SNe per cell median {np.median(np.bincount(inv))}): slope {bb[3]:+.3f} ± {ee[3]:.3f} (z = {bb[3]/ee[3]:+.2f})')
    # carrier, with x1 and c in the base
    if M['tau'] is not None:
        g, s, tau, r, T = M['gnilc'], M['sfd'], M['tau'], M['rad'], M['temp']
        for name, x in (('tau353 − SFD', tau - s), ('radiance − SFD', r - s), ('GNILC − tau353', g - tau), ('GNILC − radiance', g - r)):
            b, e, _ = fitd([x], sig); print(f'  ZTF marginal slope on {name:16s} {b[3]:+.3f} ± {e[3]:.3f} (z {b[3]/e[3]:+.2f})')
        for name, xa, xb in (('GNILC−τ | τ−SFD', g - tau, tau - s), ('GNILC−R | R−SFD', g - r, r - s)):
            b, e, _ = fitd([xa, xb], sig); print(f'  ZTF carrier {name}: {b[3]:+.3f} ± {e[3]:.3f} (z {b[3]/e[3]:+.2f}), {b[4]:+.3f} ± {e[4]:.3f} (z {b[4]/e[4]:+.2f}); corr {np.corrcoef(xa, xb)[0, 1]:+.2f}')
        Tm = T - np.nanmean(T); b, e, _ = fitd([g - s, Tm], sig); print(f'  ZTF GNILC − SFD with T: slope {b[3]:+.3f} ± {e[3]:.3f}, T {1e3*b[4]:+.2f} ± {1e3*e[4]:.2f} mmag/K')
        b, e, _ = fitd([Tm], sig); print(f'  ZTF T alone: {1e3*b[3]:+.2f} ± {1e3*e[3]:.2f} mmag/K (z {b[3]/e[3]:+.2f})')

if __name__ == '__main__':
    pantheon(sys.argv[1], sys.argv[2], sys.argv[4]); ztf(sys.argv[3], sys.argv[4])
