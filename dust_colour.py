"""Two tests of the SFD error with the SN data themselves (sweep-2.md row 40 and finding.md's "to settle it"):
 (1) SALT2 colour c vs the map differences: if the Milky-Way E(B-V) used in the fit (0.86 SFD) is low by dE, the fitted c
     is redder by 1.1 dE (dust_salt_response.py), independently of any distance. Regress c on 0.86 (Zhou - SFD), on
     (GNILC - 0.86 SFD) and on SFD itself, diagonal errors (cERR plus the sample's intrinsic colour scatter), Pantheon+ and
     ZTF; the expected slope on a true E(B-V) error is +1.1. Also x1 (expected +0.2) and mB (expected +4.0).
 (2) Correct the residuals with the DESI-stars map: dE = 0.86 (Zhou - SFD), dmu = SALT2 response x dE, y' = y - dmu; the
     slopes on (Zhou - SFD) and on (GNILC - SFD) before and after. If the DESI map is the truth, the first goes to zero and
     the second drops by 0.65 x cov/var ~ 0.3-0.4; what is left on GNILC - SFD is the part no stellar map explains.
usage: python dust_colour.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/ztf/ztfsniadr2_lite/tables/snia_data.csv data/dustmaps data/zhou/desi_dust_gr_256.fits
"""
import sys, numpy as np, pandas as pd, healpy as hp
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import dust_settle as ds

def wls(y, cols, w, sel):
    X = np.column_stack([np.ones(sel.sum())] + [c[sel] for c in cols])*w[sel][:, None]; yw = y[sel]*w[sel]
    b, *_ = np.linalg.lstsq(X, yw, rcond=None); e = np.sqrt(np.diag(np.linalg.inv(X.T@X))); return b[1:], e[1:]

def main(dat, covf, ztfcsv, data_dir, zhou):
    t = fits.getdata(zhou, 1); okz = t['N_STAR'] > 0; ze = np.full(hp.nside2npix(256), np.nan); zs = np.full(hp.nside2npix(256), np.nan)
    ze[t['HPXPIXEL'][okz]] = t['EBV_GR'][okz]; zs[t['HPXPIXEL'][okz]] = t['EBV_SFD'][okz]
    from dustmaps.config import config; config['data_dir'] = data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.planck import PlanckGNILCQuery
    sq, gq = SFDQuery(), PlanckGNILCQuery()
    # Pantheon+
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f: n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    msk = (d.IS_CALIBRATOR == 0) & (d.zHD > 0.01); idx = np.where(msk)[0]
    z = d.zHD.values[idx]; y = d.m_b_corr.values[idx] - (5*np.log10(ds.dl(z))+25)
    c = SkyCoord(d.RA.values[idx]*u.deg, d.DEC.values[idx]*u.deg, frame='icrs')
    Lfull = cov[np.ix_(idx, idx)]
    def w(X, sel):
        Lc = np.linalg.cholesky(Lfull[np.ix_(sel, sel)]); return np.linalg.solve(Lc, X)
    P = ds.Sample("P+ ", y, w, [np.ones(len(y))], c, z, d.MWEBV.values[idx], 3.04)
    P.c_, P.cerr, P.x1, P.x1err, P.mB, P.mBerr = d.c.values[idx], d.cERR.values[idx], d.x1.values[idx], d.x1ERR.values[idx], d.mB.values[idx], d.mBERR.values[idx]
    # ZTF
    tz = pd.read_csv(ztfcsv)
    sel = (tz.sn_type == 'snia-cosmo') & (tz.lccoverage_flag == 1) & (tz.fitquality_flag == 1) & (tz.source == 'z_gal') \
        & (tz.x1.abs() < 3) & (tz.c.abs() < 0.3) & (tz.redshift > 0.015) & (tz.redshift < 0.1) & tz.x0.notna()
    tz = tz[sel].reset_index(drop=True)
    ct = SkyCoord(tz.ra.values*u.deg, tz.dec.values*u.deg, frame='icrs'); apex = SkyCoord(l=264.021*u.deg, b=48.253*u.deg, frame='galactic')
    zc = (1+tz.redshift.values)/(1 - 369.82*np.cos(ct.galactic.separation(apex).rad)/ds.C) - 1
    mb = -2.5*np.log10(tz.x0.values); x1 = tz.x1.values; cc = tz.c.values; AL, BE = ds.AL, ds.BE
    var = (2.5/np.log(10)*tz.x0_err.values/tz.x0.values)**2 + AL**2*tz.x1_err.values**2 + BE**2*tz.c_err.values**2 \
        + 2*(-AL)*(-2.5/np.log(10)/tz.x0.values)*tz.cov_x0_x1.values + 2*BE*(-2.5/np.log(10)/tz.x0.values)*tz.cov_x0_c.values \
        + 2*(-AL)*BE*tz.cov_x1_c.values + (5/np.log(10)*300/(ds.C*zc))**2 + (5/np.log(10)*tz.redshift_err.values/zc)**2 + 0.078**2
    yz = mb - (5*np.log10(ds.dl(zc))+25); wz = 1/np.sqrt(var)
    def w2(X, sel): return X*wz[sel][:, None] if X.ndim == 2 else X*wz[sel]
    Z = ds.Sample("ZTF", yz, w2, [np.ones(len(yz)), -x1, cc], ct, zc, tz.mwebv.values, 2.47)
    Z.c_, Z.cerr, Z.x1, Z.x1err, Z.mB, Z.mBerr = cc, tz.c_err.values, x1, tz.x1_err.values, mb, 2.5/np.log(10)*tz.x0_err.values/tz.x0.values
    for S in (P, Z):
        s, g = sq(S.c), gq(S.c); pix = hp.ang2pix(256, S.c.ra.deg, S.c.dec.deg, lonlat=True); zh = ze[pix] - zs[pix]; okz_ = np.isfinite(zh)
        print(f"{S.label}: n={len(S.y)}, Zhou coverage {okz_.mean()*100:.0f}%")
        print("  (1) light-curve parameters vs map differences (expected per unit true E(B-V) error: c +1.1, x1 +0.2, mB +4.0; 0.86 (Zhou-SFD) and GNILC-0.86 SFD are both on the applied scale)")
        for pname, val, err, sig_int in [("c", S.c_, S.cerr, 0.08), ("x1", S.x1, S.x1err, 0.9), ("mB", S.mB - (5*np.log10(ds.dl(S.z))+25), S.mBerr, 0.25)]:
            wgt = 1/np.sqrt(err**2 + sig_int**2); ok = np.isfinite(val) & np.isfinite(err) & (err > 0)
            for xl, x, selx in [("0.86 (Zhou - SFD)", 0.86*zh, ok & okz_), ("GNILC - 0.86 SFD", g - 0.86*s, ok), ("SFD (x0.86)", 0.86*s, ok)]:
                cols = [x] if pname != "mB" else [x, -S.x1 if False else S.x1*0]  # mB alone: its residual has x1, c structure, so also fit with x1 and c as covariates
                if pname == "mB":
                    b, e = wls(val, [x, S.x1, S.c_], wgt, selx); print(f"    {pname:2s} on {xl:18s} n={selx.sum():4d}  {b[0]:+6.2f} ± {e[0]:4.2f}  (with x1, c as covariates)")
                else:
                    b, e = wls(val, [x], wgt, selx); print(f"    {pname:2s} on {xl:18s} n={selx.sum():4d}  {b[0]:+6.2f} ± {e[0]:4.2f}")
        print("  (2) residual corrected with the DESI-stars map (dE = 0.86 (Zhou - SFD), SALT2 response):")
        dE = np.where(okz_, 0.86*zh, 0.0); dmu = ds.salt_response(S.z, S.beta)*dE; yc = S.y - dmu
        for xl, x in [("Zhou - SFD", zh), ("GNILC - SFD", g - s), ("GNILC - 0.86 SFD", g - 0.86*s)]:
            b0, e0, d0 = S.slope(S.y, [x], okz_); b1, e1, d1 = S.slope(yc, [x], okz_)
            print(f"    slope on {xl:18s} (Zhou-covered SNe, n={okz_.sum()}): before {b0[0]:+5.2f}±{e0[0]:4.2f} (dchi2 {d0:4.1f}) -> after {b1[0]:+5.2f}±{e1[0]:4.2f} (dchi2 {d1:4.1f})")
        b, e, dd = S.slope(S.y, [dmu], okz_); print(f"    the SNe want {b[0]:.2f}±{e[0]:.2f} times the DESI-map correction (1 = exactly)")

main(*sys.argv[1:6])
