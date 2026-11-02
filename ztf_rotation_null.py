"""The sky-rotation null for the ZTF DR2 half of the dust finding (the blind referee's remaining objection: ZTF's p = 7e-4
rests on diagonal errors and no spatially correlated null). Same sample and model as `ztf_tests.py` (snia-cosmo, coverage
and fit-quality flags, host spectroscopic z, |x1| < 3, |c| < 0.3, 0.015 < z < 0.1; y = mB − μ(z_cmb) = M − α x1 + β c + slope ×
template; diagonal errors with σ_int for χ²/dof = 1 and 300 km/s), the GNILC − SFD template re-sampled at rotated SN positions:
(A) rotations about the Galactic pole with N-S reflection (latitudes kept), (B) uniform SO(3) rotations. An empirical p from a
rotation null is valid whatever the error model, because the residuals and their correlations stay where they are and only
the template moves; a null z spread above 1 measures how much the diagonal errors understate the sky correlations.
usage: python ztf_rotation_null.py data/ztf/ztfsniadr2_lite/tables/snia_data.csv data/dustmaps [--n 2000]
"""
import sys, argparse, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.stats import chi2 as chi2d
from astropy.coordinates import SkyCoord
import astropy.units as u
C = 299792.458; OM = 0.3; AL, BE = 0.15, 3.1

def dc_h(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3 + 1 - OM); return C/100*np.array([quad(e, 0, zi)[0] for zi in np.atleast_1d(z)])
def mu(z): return 5*np.log10((1+z)*dc_h(z)/0.7) + 25

def gal_unit(l, b):
    l, b = np.radians(l), np.radians(b); return np.stack([np.cos(b)*np.cos(l), np.cos(b)*np.sin(l), np.sin(b)], -1)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('csv'); ap.add_argument('data_dir'); ap.add_argument('--n', type=int, default=2000); ap.add_argument('--seed', type=int, default=118)
    a = ap.parse_args(); d = pd.read_csv(a.csv)
    sel = (d.sn_type == 'snia-cosmo') & (d.lccoverage_flag == 1) & (d.fitquality_flag == 1) & (d.source == 'z_gal') \
        & (d.x1.abs() < 3) & (d.c.abs() < 0.3) & (d.redshift > 0.015) & (d.redshift < 0.1) & d.x0.notna()
    d = d[sel].reset_index(drop=True); print(f'ZTF DR2 sample: {len(d)} SNe')
    c = SkyCoord(d.ra.values*u.deg, d.dec.values*u.deg, frame='icrs'); gal = c.galactic; l, b = gal.l.deg, gal.b.deg
    apex = SkyCoord(l=264.021*u.deg, b=48.253*u.deg, frame='galactic'); cosang = np.cos(gal.separation(apex).rad)
    zcmb = (1 + d.redshift.values)/(1 - 369.82*cosang/C) - 1
    from dustmaps.config import config; config['data_dir'] = a.data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.planck import PlanckGNILCQuery
    sfd, gn = SFDQuery(), PlanckGNILCQuery()
    def template(l_, b_):
        cc = SkyCoord(l=l_*u.deg, b=b_*u.deg, frame='galactic'); t = gn(cc) - sfd(cc); return np.where(np.isfinite(t), t, np.nanmedian(t))
    mb = -2.5*np.log10(d.x0.values); x1 = d.x1.values; col = d.c.values
    var_m = (2.5/np.log(10)*d.x0_err.values/d.x0.values)**2 + AL**2*d.x1_err.values**2 + BE**2*d.c_err.values**2 \
        + 2*(-AL)*(-2.5/np.log(10)/d.x0.values)*d.cov_x0_x1.values + 2*BE*(-2.5/np.log(10)/d.x0.values)*d.cov_x0_c.values + 2*(-AL)*BE*d.cov_x1_c.values
    var_v = (5/np.log(10)*300/(C*zcmb))**2 + (5/np.log(10)*d.redshift_err.values/zcmb)**2
    y = mb - mu(zcmb); base = [np.ones(len(y)), -x1, col]
    def gls(cols, sig_int):
        wgt = 1/np.sqrt(var_m + var_v + sig_int**2); X = np.column_stack(base + cols)*wgt[:, None]; yw = y*wgt
        bb, *_ = np.linalg.lstsq(X, yw, rcond=None); r = yw - X@bb; return bb, np.sqrt(np.diag(np.linalg.inv(X.T@X))), r@r
    lo, hi = 0.0, 0.5
    for _ in range(40):
        s_ = (lo + hi)/2; _, _, chi = gls([], s_)
        if chi > len(y) - 3: lo = s_
        else: hi = s_
    sig_int = (lo + hi)/2; b0, e0, c0 = gls([], sig_int)
    t = template(l, b); bb, ee, c1 = gls([t], sig_int); sl, er = bb[3], ee[3]; z0 = sl/er; dchi0 = c0 - c1
    print(f'sigma_int {sig_int:.3f}, alpha {b0[1]:.3f}, beta {b0[2]:.3f}; data slope on GNILC − SFD {sl:+.3f} ± {er:.3f} (z = {z0:+.2f}, Δχ² = {dchi0:.1f}, Gaussian p = {chi2d.sf(dchi0, 1):.2g})')
    rng = np.random.default_rng(a.seed)
    for kind in ('A', 'B'):
        zs, dcs, sls = [], [], []
        for i in range(a.n):
            if kind == 'A': l2 = (l + rng.uniform(0, 360)) % 360; b2 = b*(1 if rng.uniform() < 0.5 else -1)
            else:
                v = gal_unit(l, b); q = rng.normal(size=4); q /= np.linalg.norm(q); w0, x, yq, zq = q
                R = np.array([[1-2*(yq*yq+zq*zq), 2*(x*yq-zq*w0), 2*(x*zq+yq*w0)], [2*(x*yq+zq*w0), 1-2*(x*x+zq*zq), 2*(yq*zq-x*w0)], [2*(x*zq-yq*w0), 2*(yq*zq+x*w0), 1-2*(x*x+yq*yq)]])
                v2 = v@R.T; l2 = np.degrees(np.arctan2(v2[:, 1], v2[:, 0])) % 360; b2 = np.degrees(np.arcsin(np.clip(v2[:, 2], -1, 1)))
            b2_, e2_, c2 = gls([template(l2, b2)], sig_int); zs.append(b2_[3]/e2_[3]); dcs.append(c0 - c2); sls.append(b2_[3])
        zs, dcs, sls = map(np.array, (zs, dcs, sls))
        print(f'null {kind} ({a.n} draws): slope {sls.mean():+.3f} ± {sls.std():.3f} (data {sl:+.3f}); z spread {zs.std():.2f}; P(|z| ≥ {abs(z0):.2f}) = {np.mean(np.abs(zs) >= abs(z0)):.4f}; '
              f'P(same-sign ≥) = {np.mean(zs >= z0) if z0 > 0 else np.mean(zs <= z0):.4f}; Δχ² P(≥ {dchi0:.1f}) = {np.mean(dcs >= dchi0):.4f}; data z / null spread = {z0/zs.std():+.2f}')

main()
