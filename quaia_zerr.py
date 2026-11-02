"""Quaia redshift errors across the sky (sweep-2.md row 23): Quaia's spectrophotometric redshifts were validated against SDSS
(northern cap only). DESI DR1's QSO clustering catalogues (both caps, 0.8 < z < 2.1 plus the Lya sample, 1.2 M redshifts)
give a truth on both caps: cross-match Quaia to DESI within 1", then the error dz = (z_quaia - z_desi)/(1 + z_desi), its
median and scaled MAD, by cap, by |b|, by SFD E(B-V), by the map differences, and the fraction of catastrophic errors
(|dz| > 0.1) by the same splits. A sky-dependent bias in the Quaia z is a systematic for its tomographic dipole.
usage: python quaia_zerr.py data/quaia/quaia_G20.5.fits data/desi_lss data/dustmaps data/zhou/desi_dust_gr_256.fits
"""
import sys, numpy as np, healpy as hp
from astropy.io import fits
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u

def mad(x): return 1.4826*np.median(np.abs(x - np.median(x)))

def main(cat, lss, data_dir, zhou):
    q = fits.getdata(cat, 1)
    parts = [fits.getdata(f"{lss}/QSO_{cap}_clustering.dat.fits", 1) for cap in ('NGC', 'SGC')]
    ra = np.concatenate([p['RA'] for p in parts]); dec = np.concatenate([p['DEC'] for p in parts]); zd = np.concatenate([p['Z'] for p in parts]); cap = np.concatenate([np.full(len(p), i) for i, p in enumerate(parts)])
    cq = SkyCoord(q['ra']*u.deg, q['dec']*u.deg); cd = SkyCoord(ra*u.deg, dec*u.deg)
    i, sep, _ = match_coordinates_sky(cq, cd); ok = sep.arcsec < 1.0
    zq = q['redshift_quaia'][ok]; zt = zd[i[ok]]; capm = cap[i[ok]]; cm = cq[ok]; b = cm.galactic.b.deg
    dz = (zq - zt)/(1 + zt); cat_ = np.abs(dz) > 0.1
    print(f"Quaia x DESI DR1 QSO: {ok.sum()} matches within 1\" (NGC {np.sum(capm == 0)}, SGC {np.sum(capm == 1)}); overall median dz {np.median(dz):+.4f}, sMAD {mad(dz):.4f}, catastrophic (|dz|>0.1) {cat_.mean()*100:.1f}%")
    from dustmaps.config import config; config['data_dir'] = data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.planck import PlanckGNILCQuery
    s = SFDQuery()(cm); g = PlanckGNILCQuery()(cm)
    t = fits.getdata(zhou, 1); okz = t['N_STAR'] > 0; ze = np.full(hp.nside2npix(256), np.nan); zs = np.full(hp.nside2npix(256), np.nan)
    ze[t['HPXPIXEL'][okz]] = t['EBV_GR'][okz]; zs[t['HPXPIXEL'][okz]] = t['EBV_SFD'][okz]; pix = hp.ang2pix(256, cm.ra.deg, cm.dec.deg, lonlat=True); zh = ze[pix] - zs[pix]
    def line(lbl, sel):
        if sel.sum() < 200: return
        x = dz[sel]; n = sel.sum(); print(f"  {lbl:28s} n={n:7d} median dz {np.median(x):+.4f}±{mad(x)/np.sqrt(n)*1.25:.4f}  sMAD {mad(x):.4f}  catastrophic {cat_[sel].mean()*100:5.1f}%  (median z_desi {np.median(zt[sel]):.2f})")
    print("by cap:"); line("NGC", capm == 0); line("SGC", capm == 1)
    print("by |b|:")
    for lo, hi in [(15, 30), (30, 45), (45, 60), (60, 90)]: line(f"{lo}<=|b|<{hi}", (np.abs(b) >= lo) & (np.abs(b) < hi))
    print("by SFD E(B-V):")
    for lo, hi in [(0, 0.02), (0.02, 0.04), (0.04, 0.08), (0.08, 1)]: line(f"{lo:.2f}<=E<{hi:.2f}", (s >= lo) & (s < hi))
    for name, x in [("GNILC - SFD", g - s), ("Zhou - SFD", zh)]:
        okx = np.isfinite(x); qs = np.nanpercentile(x[okx], [0, 25, 50, 75, 100]); print(f"by {name} quartiles:")
        for k in range(4): line(f"{qs[k]:+.4f}..{qs[k+1]:+.4f}", okx & (x >= qs[k]) & (x <= qs[k+1]))
    print("by DESI redshift:")
    for lo, hi in [(0.8, 1.2), (1.2, 1.6), (1.6, 2.1), (2.1, 3.5)]: line(f"{lo}<=z<{hi}", (zt >= lo) & (zt < hi))
    # the dipole of the bias: median dz in the two hemispheres about the CMB apex and about the Galactic poles
    n_eq = np.stack([np.cos(np.radians(cm.dec.deg))*np.cos(np.radians(cm.ra.deg)), np.cos(np.radians(cm.dec.deg))*np.sin(np.radians(cm.ra.deg)), np.sin(np.radians(cm.dec.deg))], -1)
    cmb = np.array([np.cos(np.radians(-6.944))*np.cos(np.radians(167.942)), np.cos(np.radians(-6.944))*np.sin(np.radians(167.942)), np.sin(np.radians(-6.944))]); pr = n_eq@cmb
    print("hemispheres:"); line("toward the CMB apex", pr > 0); line("away from the CMB apex", pr < 0); line("Galactic north", b > 0); line("Galactic south", b < 0)

main(*sys.argv[1:5])
