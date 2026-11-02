"""How much do the dust maps agree about where SFD is wrong? Correlation matrix, at the Pantheon+ SN positions (calibrators
out, z > 0.01, |b| > 20), of the differences from SFD of GNILC, Zhou+2024 (DESI stars), Mudur+2023 (PS1/2MASS/Gaia stars),
Bayestar 2019 (far edge), Lenz 2017 (HI), Edenhofer 2023 (1.25 kpc), each scaled to SFD by one factor k; point values.
Also each difference's rms and, for Zhou, its quoted pixel error (so the intrinsic rms of the SFD error can be read off).
usage: python dust_mapcorr.py data/Pantheon+SH0ES.dat data/dustmaps data/zhou/desi_dust_gr_256.fits data/mudur/recon_fw15_final_mult.fits
"""
import sys, numpy as np, pandas as pd, healpy as hp
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import dust_settle as ds

def main(dat, data_dir, zhou, mudur):
    d = pd.read_csv(dat, sep=r'\s+'); m = (d.IS_CALIBRATOR == 0) & (d.zHD > 0.01)
    c = SkyCoord(d.RA.values[m]*u.deg, d.DEC.values[m]*u.deg, frame='icrs'); b = c.galactic.b.deg
    mp = ds.maps(c, data_dir)
    t = fits.getdata(zhou, 1); ok = t['N_STAR'] > 0; ze = np.full(hp.nside2npix(256), np.nan); zs = np.full(hp.nside2npix(256), np.nan); zerr = np.full(hp.nside2npix(256), np.nan)
    ze[t['HPXPIXEL'][ok]] = t['EBV_GR'][ok]; zs[t['HPXPIXEL'][ok]] = t['EBV_SFD'][ok]; zerr[t['HPXPIXEL'][ok]] = t['EBV_GR_ERR'][ok]
    pix = hp.ang2pix(256, c.ra.deg, c.dec.deg, lonlat=True); mp['Zhou'] = ze[pix]; zhou_sfd = zs[pix]; zhou_err = zerr[pix]
    mu = fits.getdata(mudur, 1)['Recon_Mean'].ravel().astype(float); mu[mu < -1] = np.nan; mp['Mudur'] = mu[hp.ang2pix(2048, c.galactic.l.deg, c.galactic.b.deg, lonlat=True)]
    s = mp['SFD']; diffs = {}
    for name in ('GNILC', 'Zhou', 'Mudur', 'Bayestar19', 'Lenz17', 'Edenhofer'):
        if name not in mp: continue
        v = mp[name]; okv = np.isfinite(v) & (v > -0.5) & (v < 5) & (np.abs(b) > 20)
        if name == 'GNILC': diffs[name] = np.where(okv, v - s, np.nan)
        elif name == 'Zhou': diffs[name] = np.where(okv, v - zhou_sfd, np.nan)
        else:
            k = np.sum(s[okv]*v[okv])/np.sum(v[okv]**2); diffs[name] = np.where(okv, k*v - s, np.nan)
    names = list(diffs); print("rms of each difference from SFD at the SN positions (|b|>20), and the number of SNe:")
    for n in names: x = diffs[n]; print(f"  {n:11s} rms {np.nanstd(x):.4f}  n={np.isfinite(x).sum()}")
    print(f"  Zhou pixel error median {np.nanmedian(zhou_err[np.abs(b) > 20]):.4f} -> intrinsic rms of (Zhou - SFD) ~ {np.sqrt(max(np.nanstd(diffs['Zhou'])**2 - np.nanmedian(zhou_err[np.abs(b) > 20])**2, 0)):.4f}")
    print("correlation matrix (pairwise complete):")
    print("  " + " ".join(f"{n:>11s}" for n in names))
    for a in names:
        row = []
        for bb in names:
            okp = np.isfinite(diffs[a]) & np.isfinite(diffs[bb]); row.append(np.corrcoef(diffs[a][okp], diffs[bb][okp])[0, 1] if okp.sum() > 30 else np.nan)
        print(f"  {a:11s}" + " ".join(f"{r:11.2f}" for r in row))

main(*sys.argv[1:5])
