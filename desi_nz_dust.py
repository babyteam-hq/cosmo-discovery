"""DESI DR1 LRG clustering catalogues (LSScats v1.5, NGC + SGC, 2.2 M galaxies): does the redshift distribution shift with
Galactic latitude, with SFD E(B-V), and with the dust-map differences (GNILC - SFD, Zhou+2024 DESI stars - SFD) at the
galaxy position (sweep-2.md row 4)? The LRG colour selection uses SFD-dereddened magnitudes, so where the true reddening
differs from SFD the selected colours, and so n(z), shift. Weighted (WEIGHT) mean z and quantiles in bins, and a
regression of z on [1, |b|/90, SFD, GNILC - SFD, Zhou - SFD, south-cap flag], with errors from the per-galaxy scatter.
usage: python desi_nz_dust.py data/desi_lss data/dustmaps data/zhou/desi_dust_gr_256.fits
"""
import sys, numpy as np, healpy as hp
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

def main(lss, data_dir, zhou):
    parts = []
    for cap in ('NGC', 'SGC'):
        t = fits.getdata(f"{lss}/LRG_{cap}_clustering.dat.fits", 1); parts.append((t['RA'], t['DEC'], t['Z'], t['WEIGHT'], t['PHOTSYS'], np.full(len(t), cap == 'SGC')))
    ra, dec, z, w, ps, sgc = [np.concatenate([p[i] for p in parts]) for i in range(6)]
    ok = (z > 0.4) & (z < 1.1); ra, dec, z, w, ps, sgc = ra[ok], dec[ok], z[ok], w[ok], ps[ok], sgc[ok]
    print(f"LRG: {len(z)} galaxies at 0.4<z<1.1; mean z {np.average(z, weights=w):.4f}, std {z.std():.3f}; PHOTSYS: {dict(zip(*np.unique(ps, return_counts=True)))}")
    c = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs'); b = c.galactic.b.deg
    from dustmaps.config import config; config['data_dir'] = data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.planck import PlanckGNILCQuery
    # maps at nside 256 pixel centres (13.7'), to match Zhou's resolution and keep the queries cheap
    pix = hp.ang2pix(256, ra, dec, lonlat=True); up = np.unique(pix); pra, pdec = hp.pix2ang(256, up, lonlat=True); pc = SkyCoord(pra*u.deg, pdec*u.deg, frame='icrs')
    sfd_p, gn_p = SFDQuery()(pc), PlanckGNILCQuery()(pc); lut = dict(zip(up, range(len(up)))); ii = np.array([lut[p] for p in pix])
    s, g = sfd_p[ii], gn_p[ii]
    t = fits.getdata(zhou, 1); okz = t['N_STAR'] > 0; ze = np.full(hp.nside2npix(256), np.nan); zs = np.full(hp.nside2npix(256), np.nan)
    ze[t['HPXPIXEL'][okz]] = t['EBV_GR'][okz]; zs[t['HPXPIXEL'][okz]] = t['EBV_SFD'][okz]; zh = ze[pix] - zs[pix]
    print(f"Zhou coverage of the LRGs {np.isfinite(zh).mean()*100:.1f}%; rms(GNILC-SFD) {np.std(g-s):.4f}, rms(Zhou-SFD) {np.nanstd(zh):.4f}")
    def wmean(sel):
        zz, ww = z[sel], w[sel]; m = np.average(zz, weights=ww); e = np.sqrt(np.average((zz-m)**2, weights=ww)/sel.sum()); return m, e, np.percentile(zz, [25, 50, 75])
    print("weighted mean z by |b| (and the quartiles):")
    for lo, hi in [(15, 25), (25, 35), (35, 50), (50, 65), (65, 90)]:
        sel = (np.abs(b) >= lo) & (np.abs(b) < hi); m, e, q = wmean(sel); print(f"  {lo:2d}<=|b|<{hi:2d} n={sel.sum():7d} mean z {m:.4f}±{e:.4f} quartiles {q[0]:.3f} {q[1]:.3f} {q[2]:.3f} mean SFD {s[sel].mean():.3f}")
    print("by SFD E(B-V):")
    for lo, hi in [(0, 0.02), (0.02, 0.04), (0.04, 0.07), (0.07, 0.12), (0.12, 1)]:
        sel = (s >= lo) & (s < hi); m, e, q = wmean(sel); print(f"  {lo:.2f}<=E<{hi:.2f} n={sel.sum():7d} mean z {m:.4f}±{e:.4f} median {q[1]:.3f} mean |b| {np.abs(b[sel]).mean():.0f}")
    for lbl, x in [("GNILC - SFD", g - s), ("Zhou - SFD", zh)]:
        print(f"by {lbl} (|b|>30):")
        okx = np.isfinite(x) & (np.abs(b) > 30); qs = np.nanpercentile(x[okx], [0, 20, 40, 60, 80, 100])
        for i in range(5):
            sel = okx & (x >= qs[i]) & (x <= qs[i+1]); m, e, q = wmean(sel); print(f"  {qs[i]:+.4f}..{qs[i+1]:+.4f} n={sel.sum():7d} mean z {m:.4f}±{e:.4f} median {q[1]:.3f}")
    okx = np.isfinite(zh)
    X = np.column_stack([np.ones(okx.sum()), np.abs(b[okx])/90, s[okx], (g-s)[okx], zh[okx], sgc[okx].astype(float)]); sw = np.sqrt(w[okx])
    bb, *_ = np.linalg.lstsq(X*sw[:, None], z[okx]*sw, rcond=None); r = z[okx] - X@bb; cov = np.linalg.inv((X*w[okx][:, None]).T@X)*np.average(r**2, weights=w[okx]); e = np.sqrt(np.diag(cov))
    names = ["1", "|b|/90", "SFD", "GNILC-SFD", "Zhou-SFD", "SGC"]
    print("z regressed on [" + ", ".join(names) + f"] (n={okx.sum()}):"); [print(f"  {n:10s} {bb[i]:+.4f} ± {e[i]:.4f}  ({bb[i]/e[i]:+.1f} sigma)") for i, n in enumerate(names)]
    X2 = X[:, [0, 1, 2, 5]]; bb2, *_ = np.linalg.lstsq(X2*sw[:, None], z[okx]*sw, rcond=None); r2 = z[okx] - X2@bb2; e2 = np.sqrt(np.diag(np.linalg.inv((X2*w[okx][:, None]).T@X2)*np.average(r2**2, weights=w[okx])))
    print(f"  without the map differences: |b|/90 {bb2[1]:+.4f}±{e2[1]:.4f}, SFD {bb2[2]:+.4f}±{e2[2]:.4f}, SGC {bb2[3]:+.4f}±{e2[3]:.4f}")
    # what a map error does to the LRG selection: per 0.01 in E(B-V), the dereddened r-z moves by 0.01 x (R_r - R_z) = 0.01 x (2.165-1.211) = 0.0095 and z-W1 by 0.01 x (1.211-0.184)

main(*sys.argv[1:4])
