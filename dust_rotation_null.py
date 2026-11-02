"""The blind referee's objection to the dust finding (finding.md, "Blind review"): GNILC − SFD is a smooth sky pattern and
Hubble residuals are sky-correlated, so a p-value from per-SN errors overstates the number of independent samples. The right
null keeps the residuals (with their full covariance) and the survey footprint where they are and moves the template: here,
rotations of the dust maps relative to the SNe.

Two nulls, both re-sampling SFD and GNILC at the moved positions and refitting the same full-covariance GLS slope:
  A. azimuthal: a random rotation about the Galactic pole (every SN keeps its latitude, so the template's amplitude
     distribution is matched exactly), with a random north-south reflection;
  B. full SO(3) rotations (the footprint lands anywhere; template statistics differ, so B is the looser comparison).
Statistics per draw: the slope, its GLS z = slope/err and Δχ². p = fraction of draws with |z| (or Δχ²) at least the data's.
Also the same for the slope on MWEBV's residual proxy (0.86 SFD) as a control, and a 1°-smoothed-template version if asked.
usage: python dust_rotation_null.py data/pantheon/Pantheon+SH0ES.dat data/pantheon/Pantheon+SH0ES_STAT+SYS.cov data/dustmaps [--n 1000]
"""
import sys, argparse, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.stats import chi2 as chi2d
from astropy.coordinates import SkyCoord
import astropy.units as u
C = 299792.458; OM = 0.334

def dl(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3 + 1 - OM); return (1+z)*C/70*np.array([quad(e, 0, zi)[0] for zi in z])

def gal_unit(l, b):
    l, b = np.radians(l), np.radians(b); return np.stack([np.cos(b)*np.cos(l), np.cos(b)*np.sin(l), np.sin(b)], -1)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('dat'); ap.add_argument('cov'); ap.add_argument('data_dir'); ap.add_argument('--n', type=int, default=1000)
    ap.add_argument('--seed', type=int, default=118); a = ap.parse_args()
    d = pd.read_csv(a.dat, sep=r'\s+')
    with open(a.cov) as f: n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    m = (d.IS_CALIBRATOR == 0) & (d.zHD > 0.01); idx = np.where(m)[0]
    z = d.zHD.values[idx]; y = d.m_b_corr.values[idx] - (5*np.log10(dl(z)) + 25)
    c = SkyCoord(d.RA.values[idx]*u.deg, d.DEC.values[idx]*u.deg, frame='icrs').galactic; l, b = c.l.deg, c.b.deg
    from dustmaps.config import config; config['data_dir'] = a.data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.planck import PlanckGNILCQuery
    sfd, gn = SFDQuery(), PlanckGNILCQuery()
    def template(l_, b_):
        cc = SkyCoord(l=l_*u.deg, b=b_*u.deg, frame='galactic'); s = sfd(cc); g = gn(cc); return g - s, s
    L = np.linalg.cholesky(cov[np.ix_(idx, idx)]); w = lambda X: np.linalg.solve(L, X); yw = w(y)
    one = np.ones(len(z))
    def fit(x):
        X = np.column_stack([one, x]); Xw = w(X); bb, *_ = np.linalg.lstsq(Xw, yw, rcond=None); r = yw - Xw@bb
        e = np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))); return bb[1], e[1], r@r
    Xw0 = w(one[:, None]); b0, *_ = np.linalg.lstsq(Xw0, yw, rcond=None); chi0 = float(((yw - Xw0@b0)**2).sum())
    t, s = template(l, b); sl, er, ch = fit(t); z0 = sl/er; dchi0 = chi0 - ch
    print(f'N={len(z)}; data: slope on GNILC − SFD {sl:+.3f} ± {er:.3f} (z = {z0:+.2f}, Δχ² = {dchi0:.1f}, Gaussian p = {chi2d.sf(dchi0, 1):.2g})')
    sl_s, er_s, ch_s = fit(s); print(f'      control, slope on SFD itself {sl_s:+.3f} ± {er_s:.3f} (z = {sl_s/er_s:+.2f})')
    # the referee's confounding reading: does the slope survive x1, c and host mass as covariates, and survey fixed effects?
    def fitcols(cols):
        X = np.column_stack([one] + cols); Xw = w(X); bb, *_ = np.linalg.lstsq(Xw, yw, rcond=None); e = np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))); return bb[1], e[1]
    x1, col, hm, sid = (d[k].values[idx] for k in ('x1', 'c', 'HOST_LOGMASS', 'IDSURVEY')); hm = np.where(hm > 0, hm, np.nanmedian(hm[hm > 0]))
    surveys = [np.asarray(sid == s_, float) for s_ in np.unique(sid)[1:]]
    for name, cols in (('+ x1, c', [t, x1, col]), ('+ x1, c, host mass', [t, x1, col, hm]), ('+ survey fixed effects', [t] + surveys), ('+ survey FE, x1, c, host mass', [t, x1, col, hm] + surveys)):
        s1, e1 = fitcols(cols); print(f'      slope on GNILC − SFD {name:32s} {s1:+.3f} ± {e1:.3f}')
    rng = np.random.default_rng(a.seed)
    def run(kind):
        zs, dcs, sls, zs_s = [], [], [], []
        for i in range(a.n):
            if kind == 'A':
                l2 = (l + rng.uniform(0, 360)) % 360; b2 = b*(1 if rng.uniform() < 0.5 else -1)
            else:
                v = gal_unit(l, b); q = rng.normal(size=4); q /= np.linalg.norm(q)      # uniform random rotation (unit quaternion)
                w0, x, yq, zq = q; R = np.array([[1-2*(yq*yq+zq*zq), 2*(x*yq-zq*w0), 2*(x*zq+yq*w0)], [2*(x*yq+zq*w0), 1-2*(x*x+zq*zq), 2*(yq*zq-x*w0)], [2*(x*zq-yq*w0), 2*(yq*zq+x*w0), 1-2*(x*x+yq*yq)]])
                v2 = v@R.T; l2 = np.degrees(np.arctan2(v2[:, 1], v2[:, 0])) % 360; b2 = np.degrees(np.arcsin(np.clip(v2[:, 2], -1, 1)))
            t2, s2 = template(l2, b2); t2 = np.where(np.isfinite(t2), t2, np.nanmedian(t2))
            sl2, er2, ch2 = fit(t2); zs.append(sl2/er2); dcs.append(chi0 - ch2); sls.append(sl2)
            sl3, er3, _ = fit(np.where(np.isfinite(s2), s2, np.nanmedian(s2))); zs_s.append(sl3/er3)
        zs, dcs, sls, zs_s = map(np.array, (zs, dcs, sls, zs_s))
        print(f'null {kind} ({a.n} draws): slope {sls.mean():+.3f} ± {sls.std():.3f} (data {sl:+.3f}); z: std {zs.std():.2f} (1 if per-SN errors were right), '
              f'P(|z| ≥ {abs(z0):.2f}) = {np.mean(np.abs(zs) >= abs(z0)):.4f}, P(z ≤ −|z0| or ≥ |z0| one-sided same sign) = {np.mean(zs >= z0) if z0 > 0 else np.mean(zs <= z0):.4f}; '
              f'Δχ²: P(≥ {dchi0:.1f}) = {np.mean(dcs >= dchi0):.4f}; the data z in units of the null std: {z0/zs.std():+.2f}; control SFD slope z std {zs_s.std():.2f}')
        return zs
    zA = run('A'); zB = run('B')
    # a fair trials-corrected reading: the per-SN z of 3.9 against the null's std
    print(f'summary: per-SN-error z = {z0:+.2f}; against rotation null A the equivalent is {z0/zA.std():+.2f}σ (empirical p above), against B {z0/zB.std():+.2f}σ')

main()
