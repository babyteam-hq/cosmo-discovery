"""ZTF SN Ia DR2 (Rigault+2025, ztfcosmo.in2p3.fr 'lite' release) as an independent low-z sample for two tests:
  (a) the CF4++ velocity-field test of finding.md (beta for the field, its dipole, its non-dipole rest), and
  (b) the GNILC-SFD dust-map slope of results.md row 3.
Caveat the release states: DR2 photometry is not cosmology-grade (DR2.5 will be); here only relative magnitudes across
the sky and across redshift matter, with the sample's own scatter. Sample: sn_type snia-cosmo, lccoverage_flag and
fitquality_flag set, host spectroscopic redshift (source z_gal), |x1|<3, |c|<0.3, 0.015<z<0.1 (or 0.06).
Model: y = -2.5 log10 x0 - mu(z_cmb) = M - alpha x1 + beta c [+ test columns], alpha and beta free; diagonal errors:
SALT covariance (alpha 0.15, beta 3.1), sigma_int fit for chi2/dof=1, 300 km/s peculiar velocity, redshift error.
z_cmb from the heliocentric redshift and the CMB dipole (369.82 km/s toward l=264.021, b=48.253).
usage: python ztf_tests.py data/ztf/ztfsniadr2_lite/tables/snia_data.csv data/new_CF4pp_mean_std_grids.npz data/dustmaps
"""
import sys, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.ndimage import map_coordinates
from scipy.stats import chi2 as chi2d
from astropy.coordinates import SkyCoord
import astropy.units as u
C = 299792.458; OM = 0.334; N = 128; L = 1000.0; DELTA = L/N; AL, BE = 0.15, 3.1
def dc_h(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3 + 1-OM)
    return C/100*np.array([quad(e, 0, zi)[0] for zi in np.atleast_1d(z)])
def mu(z): return 5*np.log10((1+z)*dc_h(z)/0.7) + 25
def zhd(zcmb, v): return (1+zcmb)/(1+v/C) - 1

def main(csv, grid, data_dir):
    g = np.load(grid); V = g['v_mean']
    d = pd.read_csv(csv)
    sel = (d.sn_type == 'snia-cosmo') & (d.lccoverage_flag == 1) & (d.fitquality_flag == 1) & (d.source == 'z_gal') \
        & (d.x1.abs() < 3) & (d.c.abs() < 0.3) & (d.redshift > 0.015) & (d.redshift < 0.1) & d.x0.notna()
    d = d[sel].reset_index(drop=True); print(f"ZTF DR2 sample: {len(d)} SNe")
    c = SkyCoord(d.ra.values*u.deg, d.dec.values*u.deg, frame='icrs'); gal = c.galactic
    apex = SkyCoord(l=264.021*u.deg, b=48.253*u.deg, frame='galactic')
    cosang = np.cos(gal.separation(apex).rad); zcmb = (1+d.redshift.values)/(1 - 369.82*cosang/C) - 1
    from dustmaps.config import config; config['data_dir'] = data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.planck import PlanckGNILCQuery, PlanckQuery
    sfd = SFDQuery()(c); gn = PlanckGNILCQuery()(c); p13 = PlanckQuery(component='ebv')(c)
    print(f"ZTF mwebv / SFD: median {np.median(d.mwebv.values/sfd):.3f} (1 = raw SFD, 0.86 = Schlafly-scaled)")
    mb = -2.5*np.log10(d.x0.values); x1 = d.x1.values; cc = d.c.values
    var_m = (2.5/np.log(10)*d.x0_err.values/d.x0.values)**2 + AL**2*d.x1_err.values**2 + BE**2*d.c_err.values**2 \
        + 2*(-AL)*(-2.5/np.log(10)/d.x0.values)*d.cov_x0_x1.values + 2*BE*(-2.5/np.log(10)/d.x0.values)*d.cov_x0_c.values \
        + 2*(-AL)*BE*d.cov_x1_c.values
    var_v = (5/np.log(10)*300/(C*zcmb))**2 + (5/np.log(10)*d.redshift_err.values/zcmb)**2
    y = mb - mu(zcmb)
    base = [np.ones(len(y)), -x1, cc]
    def gls(cols, sig_int):
        wgt = 1/np.sqrt(var_m + var_v + sig_int**2)
        X = np.column_stack(base + cols)*wgt[:, None]; yw = y*wgt
        b, *_ = np.linalg.lstsq(X, yw, rcond=None); r = yw - X@b
        return b, np.sqrt(np.diag(np.linalg.inv(X.T@X))), r@r
    # sigma_int so that chi2/dof = 1 for the base model
    lo, hi = 0.0, 0.5
    for _ in range(40):
        s = (lo+hi)/2; _, _, chi = gls([], s)
        if chi > len(y)-3: lo = s
        else: hi = s
    sig_int = (lo+hi)/2; b0, e0, c0 = gls([], sig_int)
    print(f"sigma_int={sig_int:.3f}; alpha={b0[1]:.3f}±{e0[1]:.3f} beta={b0[2]:.3f}±{e0[2]:.3f}; chi2={c0:.0f}/{len(y)-3}; rms residual {np.std(y - np.column_stack(base)@b0):.3f}")
    def test(cols, label, dof=None):
        b, e, c1 = gls(cols, sig_int); dc = c0 - c1; k = len(base); dof = dof or len(cols)
        print(f"  {label:50s} " + " ".join(f"b{j}={b[k+j]:6.2f}±{e[k+j]:4.2f}" for j in range(len(cols))) + f"  dchi2={dc:5.2f} p={chi2d.sf(dc, dof):.2g}")
        return b[k:], e[k:]
    print("(b) dust maps at the SN:")
    test([gn - sfd], "GNILC - SFD (Pantheon+: +1.75±0.45)")
    test([p13 - sfd], "Planck13 tau353 - SFD (Pantheon+: +0.17±0.41)")
    test([d.mwebv.values], "ZTF mwebv")
    test([sfd, gn], "y on SFD (b0) and GNILC (b1)")
    # cuts on the GNILC-SFD slope (base refit on each subset)
    bdeg = np.abs(gal.b.deg); x = gn - sfd
    for label, s in [("|b|>20", bdeg > 20), ("|b|>30", bdeg > 30), ("|b|<30", bdeg < 30), ("z<0.06", zcmb < 0.06), ("z>=0.06", zcmb >= 0.06),
                     ("Dec>30", d.dec.values > 30), ("Dec<30", d.dec.values < 30), ("c<0.1", cc < 0.1), ("mwebv<0.05", d.mwebv.values < 0.05),
                     ("|GNILC-SFD|<0.02", np.abs(x) < 0.02), ("b>0", gal.b.deg > 0), ("b<0", gal.b.deg < 0)]:
        wgt = 1/np.sqrt(var_m[s] + var_v[s] + sig_int**2)
        X0 = np.column_stack([np.ones(s.sum()), -x1[s], cc[s]])*wgt[:, None]; yw = y[s]*wgt
        X = np.column_stack([X0, x[s]*wgt]); b, *_ = np.linalg.lstsq(X, yw, rcond=None); e = np.sqrt(np.diag(np.linalg.inv(X.T@X)))
        b_, *_ = np.linalg.lstsq(X0, yw, rcond=None); dc = (yw - X0@b_)@(yw - X0@b_) - (yw - X@b)@(yw - X@b)
        print(f"    cut {label:18s} n={s.sum():4d} GNILC-SFD slope {b[3]:5.2f}±{e[3]:4.2f} dchi2={dc:5.1f}")
    print("(a) CF4++ velocity field at the SN positions:")
    for zmax in (0.1, 0.06):
        m = zcmb < zmax; idx = np.where(m)[0]
        r = dc_h(zcmb[idx]); sc = SkyCoord(ra=d.ra.values[idx]*u.deg, dec=d.dec.values[idx]*u.deg, distance=r*u.Mpc, frame='fk5', equinox='J2000').supergalactic
        xyz = np.stack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value]); rhat = xyz/np.linalg.norm(xyz, axis=0)
        ci = (xyz + L/2)/DELTA - 0.5
        vcf = np.sum(np.stack([map_coordinates(V[k], ci, order=1, mode='nearest') for k in range(3)])*rhat, 0)
        A = np.column_stack([r, rhat.T]); coef, *_ = np.linalg.lstsq(A, vcf, rcond=None)
        dip = A@np.r_[0, coef[1:]]; rest = vcf - A@coef
        ysub = y[idx]
        print(f"  z<{zmax}: n={len(idx)}, CF4++ radial v at SNe rms {np.sqrt(np.mean(vcf**2)):.0f} km/s, dipole rms {np.sqrt(np.mean(dip**2)):.0f}, rest rms {np.sqrt(np.mean(rest**2)):.0f}")
        wgt = 1/np.sqrt(var_m[idx] + var_v[idx] + sig_int**2)
        def sub(cols, label):
            X = np.column_stack([np.ones(len(idx)), -x1[idx], cc[idx]] + cols)*wgt[:, None]; yw = ysub*wgt
            b, *_ = np.linalg.lstsq(X, yw, rcond=None); e = np.sqrt(np.diag(np.linalg.inv(X.T@X))); rr = yw - X@b
            X0 = np.column_stack([np.ones(len(idx)), -x1[idx], cc[idx]])*wgt[:, None]; b_, *_ = np.linalg.lstsq(X0, yw, rcond=None); r0 = yw - X0@b_
            dc = r0@r0 - rr@rr; print(f"    {label:30s} beta={b[3]:5.2f}±{e[3]:4.2f}  dchi2={dc:5.1f}")
        for v, label in [(vcf, "CF4++ full field"), (dip, "CF4++ dipole part"), (rest, "CF4++ non-dipole rest")]:
            sub([mu(zhd(zcmb[idx], v)) - mu(zcmb[idx])], label)

main(*sys.argv[1:4])
