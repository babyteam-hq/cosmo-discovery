"""Follow-ups to dust_colour.py: (a) Pantheon+'s BBC bias correction (biasCor_m_b) and its raw and corrected magnitudes vs
the map differences, to see which ingredient of m_b_corr carries the GNILC slope; (b) SALT2 c vs 1-degree averages of
(Zhou - SFD) and (GNILC - SFD) (the SN-distance slope grew with averaging; does the colour one?); (c) c vs the map
differences by survey group in Pantheon+ (low-z BVRI surveys vs griz surveys).
usage: python dust_colour2.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/ztf/ztfsniadr2_lite/tables/snia_data.csv data/dustmaps data/zhou/desi_dust_gr_256.fits
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
    d = pd.read_csv(dat, sep=r'\s+'); msk = (d.IS_CALIBRATOR == 0) & (d.zHD > 0.01); idx = np.where(msk)[0]; p = d.iloc[idx].reset_index(drop=True)
    c = SkyCoord(p.RA.values*u.deg, p.DEC.values*u.deg, frame='icrs'); s, g = sq(c), gq(c); pix = hp.ang2pix(256, c.ra.deg, c.dec.deg, lonlat=True); zh = ze[pix] - zs[pix]; okz_ = np.isfinite(zh)
    mu = 5*np.log10(ds.dl(p.zHD.values)) + 25
    print("(a) Pantheon+: ingredients of m_b_corr vs the map differences (diagonal weights 1/m_b_corr_err_DIAG)")
    wgt = 1/p.m_b_corr_err_DIAG.values
    raw = p.mB.values + 0.148*p.x1.values - 3.04*p.c.values - mu   # un-bias-corrected standardised residual, alpha=0.148, beta=3.04
    for lbl, val in [("m_b_corr - mu (y)", p.m_b_corr.values - mu), ("mB + 0.148 x1 - 3.04 c - mu", raw), ("biasCor_m_b", p.biasCor_m_b.values), ("mB - mu", p.mB.values - mu), ("c", p.c.values), ("x1", p.x1.values)]:
        for xl, x, selx in [("GNILC - 0.86 SFD", g - 0.86*s, np.ones(len(p), bool)), ("0.86 (Zhou - SFD)", 0.86*zh, okz_)]:
            b, e = wls(val, [x], wgt, selx); print(f"    {lbl:28s} on {xl:18s}: {b[0]:+6.2f} ± {e[0]:4.2f}")
    print("(b) colour vs 1-degree averages:")
    sm = ds.smooth(c, {'SFD': 1, 'GNILC': 1}, data_dir); vec = hp.ang2vec(c.ra.deg, c.dec.deg, lonlat=True)
    zsm = np.array([np.nanmean(ze[hp.query_disc(256, v, np.radians(1.0))] - zs[hp.query_disc(256, v, np.radians(1.0))]) if np.isfinite(ze[hp.query_disc(256, v, np.radians(1.0))]).sum() > 10 else np.nan for v in vec])
    wc = 1/np.sqrt(p.cERR.values**2 + 0.08**2)
    for xl, x in [("0.86 (Zhou - SFD) 1deg", 0.86*zsm), ("GNILC - 0.86 SFD 1deg", sm['GNILC_sm'] - 0.86*sm['SFD_sm'])]:
        selx = np.isfinite(x); b, e = wls(p.c.values, [x], wc, selx); print(f"    c on {xl:24s} n={selx.sum()}: {b[0]:+6.2f} ± {e[0]:4.2f} (expected +1.1 for a true E(B-V) error)")
    print("(c) colour vs the map differences by Pantheon+ survey group:")
    lowz = p.IDSURVEY.isin([5, 50, 51, 56, 57, 58, 61, 62, 63, 64, 65, 66]).values; griz = p.IDSURVEY.isin([1, 4, 10, 15, 150]).values
    for gl, sel in [("low-z BVRI-like (CfA, CSP, LOSS, ...)", lowz), ("griz (SDSS, SNLS, DES, PS1, Foundation)", griz)]:
        for xl, x, selx in [("GNILC - 0.86 SFD", g - 0.86*s, sel), ("0.86 (Zhou - SFD)", 0.86*zh, sel & okz_)]:
            if selx.sum() < 50: continue
            b, e = wls(p.c.values, [x], wc, selx); b2, e2 = wls(p.m_b_corr.values - mu, [x], wgt, selx)
            print(f"    {gl:42s} n={selx.sum():4d}  c on {xl:18s}: {b[0]:+6.2f} ± {e[0]:4.2f} | y on it: {b2[0]:+6.2f} ± {e2[0]:4.2f}")

main(*sys.argv[1:6])
