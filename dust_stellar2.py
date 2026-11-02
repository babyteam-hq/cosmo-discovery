"""Two more emission-free reddening maps at the SN positions (sweep-2.md row 38), in the frame of dust_settle.py part (a):
  Zhou+2024 (arXiv:2409.05140): DESI stellar reddening map, 2.6 M stars, nside 256 (13.7'), E(B-V) from g-r on the SFD scale,
    with SFD sampled at the same stars (EBV_SFD), so the per-pixel difference EBV_GR - EBV_SFD is the map's own; pixels with
    N_STAR > 0 only (no inpainting). DESI footprint (~1/3 of the sky, mostly |b| > 20).
  Mudur+2023 (arXiv:2212.04514): PS1/2MASS/Gaia stellar reddening map built for cosmology, 15' FWHM version, Dec > -30.
For each: the residual slope on (k MAP - SFD) alone and jointly with (GNILC - SFD), point values and 1-degree disk averages,
and the slope expected if GNILC were right and the map a noisy measure of the same truth (slope_GNILC x cov/var).
usage: python dust_stellar2.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/ztf/ztfsniadr2_lite/tables/snia_data.csv data/dustmaps data/zhou/desi_dust_gr_256.fits data/mudur/recon_fw15_final_mult.fits
"""
import sys, os, numpy as np, pandas as pd, healpy as hp
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import dust_settle as ds

def zhou_maps(path):
    t = fits.getdata(path, 1); n = 256; ebv = np.full(hp.nside2npix(n), np.nan); sfd = np.full(hp.nside2npix(n), np.nan); err = np.full(hp.nside2npix(n), np.nan)
    ok = t['N_STAR'] > 0; ebv[t['HPXPIXEL'][ok]] = t['EBV_GR'][ok]; sfd[t['HPXPIXEL'][ok]] = t['EBV_SFD'][ok]; err[t['HPXPIXEL'][ok]] = t['EBV_GR_ERR'][ok]
    return ebv, sfd, err

def mudur_map(path):
    """Mudur+2023 recon map: nside 2048 ring, column Recon_Mean on the SFD E(B-V) scale; masked pixels hold -1.6e30."""
    m = fits.getdata(path, 1)['Recon_Mean'].ravel().astype(float); m[m < -1] = np.nan; return m

def main(dat, covf, ztfcsv, data_dir, zhou, mudur):
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f: n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    msk = (d.IS_CALIBRATOR == 0) & (d.zHD > 0.01); idx = np.where(msk)[0]
    z = d.zHD.values[idx]; y = d.m_b_corr.values[idx] - (5*np.log10(ds.dl(z))+25)
    c = SkyCoord(d.RA.values[idx]*u.deg, d.DEC.values[idx]*u.deg, frame='icrs')
    Lfull = cov[np.ix_(idx, idx)]
    def w(X, sel):
        Lc = np.linalg.cholesky(Lfull[np.ix_(sel, sel)]); return np.linalg.solve(Lc, X)
    P = ds.Sample("P+ ", y, w, [np.ones(len(y))], c, z, d.MWEBV.values[idx], 3.04)
    t = pd.read_csv(ztfcsv)
    sel = (t.sn_type == 'snia-cosmo') & (t.lccoverage_flag == 1) & (t.fitquality_flag == 1) & (t.source == 'z_gal') \
        & (t.x1.abs() < 3) & (t.c.abs() < 0.3) & (t.redshift > 0.015) & (t.redshift < 0.1) & t.x0.notna()
    t = t[sel].reset_index(drop=True)
    ct = SkyCoord(t.ra.values*u.deg, t.dec.values*u.deg, frame='icrs'); apex = SkyCoord(l=264.021*u.deg, b=48.253*u.deg, frame='galactic')
    zc = (1+t.redshift.values)/(1 - 369.82*np.cos(ct.galactic.separation(apex).rad)/ds.C) - 1
    mb = -2.5*np.log10(t.x0.values); x1 = t.x1.values; cc = t.c.values; AL, BE = ds.AL, ds.BE
    var = (2.5/np.log(10)*t.x0_err.values/t.x0.values)**2 + AL**2*t.x1_err.values**2 + BE**2*t.c_err.values**2 \
        + 2*(-AL)*(-2.5/np.log(10)/t.x0.values)*t.cov_x0_x1.values + 2*BE*(-2.5/np.log(10)/t.x0.values)*t.cov_x0_c.values \
        + 2*(-AL)*BE*t.cov_x1_c.values + (5/np.log(10)*300/(ds.C*zc))**2 + (5/np.log(10)*t.redshift_err.values/zc)**2 + 0.078**2
    yz = mb - (5*np.log10(ds.dl(zc))+25); wz = 1/np.sqrt(var)
    def w2(X, sel): return X*wz[sel][:, None] if X.ndim == 2 else X*wz[sel]
    Z = ds.Sample("ZTF", yz, w2, [np.ones(len(yz)), -x1, cc], ct, zc, t.mwebv.values, 2.47)
    from dustmaps.config import config; config['data_dir'] = data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.planck import PlanckGNILCQuery
    sq, gq = SFDQuery(), PlanckGNILCQuery()
    # Zhou: point values (pixel at the SN) and 1-degree disk means, using the map's own SFD column for the point difference
    ze, zs, zerr = zhou_maps(zhou)
    mud = None
    if os.path.exists(mudur):
        try: mud = mudur_map(mudur); print("Mudur map pixels", mud.size, "nside", hp.npix2nside(mud.size), "finite", np.isfinite(mud).mean(), "median", np.nanmedian(mud))
        except Exception as e: print("Mudur map unreadable:", repr(e)[:200])
    for S in (P, Z):
        pix = hp.ang2pix(256, S.c.icrs.ra.deg, S.c.icrs.dec.deg, lonlat=True)
        mp = {'SFD': sq(S.c), 'GNILC': gq(S.c), 'Zhou': ze[pix], 'ZhouSFD': zs[pix]}
        print(f"{S.label}: Zhou coverage {np.isfinite(mp['Zhou']).mean()*100:.1f}%; Zhou pixel error median {np.nanmedian(zerr[pix]):.4f}; corr(Zhou SFD column, SFD at SN) {np.corrcoef(mp['ZhouSFD'][np.isfinite(mp['Zhou'])], mp['SFD'][np.isfinite(mp['Zhou'])])[0,1]:.3f}")
        # 1-degree disks: healpy query_disc on the nside-256 map (NaN-aware), and ring sampling for SFD/GNILC as dust_settle.smooth
        vecs = hp.ang2vec(S.c.icrs.ra.deg, S.c.icrs.dec.deg, lonlat=True)
        zsm = np.array([np.nanmean(ze[hp.query_disc(256, v, np.radians(1.0))]) if np.isfinite(ze[hp.query_disc(256, v, np.radians(1.0))]).sum() > 10 else np.nan for v in vecs])
        zssm = np.array([np.nanmean(zs[hp.query_disc(256, v, np.radians(1.0))]) if np.isfinite(zs[hp.query_disc(256, v, np.radians(1.0))]).sum() > 10 else np.nan for v in vecs])
        sm = ds.smooth(S.c, {'SFD': 1, 'GNILC': 1}, data_dir); mp['SFD_sm'] = sm['SFD_sm']; mp['GNILC_sm'] = sm['GNILC_sm']; mp['Zhou_sm'] = zsm; mp['ZhouSFDsm_sm'] = zssm
        if mud is not None:
            nm = hp.npix2nside(mud.size); gal = S.c.galactic; eq = S.c.icrs
            mg = mud[hp.ang2pix(nm, gal.l.deg, gal.b.deg, lonlat=True)]; me = mud[hp.ang2pix(nm, eq.ra.deg, eq.dec.deg, lonlat=True)]
            okg, oke = np.isfinite(mg), np.isfinite(me)
            cg = np.corrcoef(mg[okg], mp['SFD'][okg])[0, 1] if okg.sum() > 50 else -9; ce = np.corrcoef(me[oke], mp['SFD'][oke])[0, 1] if oke.sum() > 50 else -9
            galactic = cg > ce; print(f"  Mudur coordinate check: corr with SFD if galactic {cg:.3f} (coverage {okg.mean():.2f}), if equatorial {ce:.3f} (coverage {oke.mean():.2f}) -> {'galactic' if galactic else 'equatorial'}")
            mp['Mudur'] = mg if galactic else me
            vv = hp.ang2vec(gal.l.deg, gal.b.deg, lonlat=True) if galactic else hp.ang2vec(eq.ra.deg, eq.dec.deg, lonlat=True)
            mp['Mudur_sm'] = np.array([np.nanmean(mud[hp.query_disc(nm, v, np.radians(1.0))]) for v in vv])
        tests = [("Zhou - its own SFD (point)", mp['Zhou'] - mp['ZhouSFD'], mp['GNILC'] - mp['SFD']), ("Zhou - its own SFD (1 deg)", mp['Zhou_sm'] - mp['ZhouSFDsm_sm'], mp['GNILC_sm'] - mp['SFD_sm'])]
        if mud is not None:
            for lbl, m_, s_, g_ in [("Mudur (point)", mp['Mudur'], mp['SFD'], mp['GNILC']), ("Mudur (1 deg)", mp['Mudur_sm'], mp['SFD_sm'], mp['GNILC_sm'])]:
                ok = np.isfinite(m_); k = np.sum(s_[ok]*m_[ok])/np.sum(m_[ok]**2); tests.append((f"k Mudur - SFD {lbl[6:]} k={k:.3f}", k*m_ - s_, g_ - s_))
        for lbl, xm, xg in tests:
            ok = np.isfinite(xm) & np.isfinite(xg)
            for bl, sel in [("all", ok), ("|b|>30", ok & (np.abs(S.b) > 30))]:
                if sel.sum() < 50: continue
                b1, e1, d1 = S.slope(S.y, [xm], sel); b3, e3, _ = S.slope(S.y, [xg], sel); b2, e2, _ = S.slope(S.y, [xm, xg], sel)
                lam = np.cov(xm[sel], xg[sel])[0, 1]/np.var(xm[sel]); ext = np.median(ds.salt_response(S.z[sel], S.beta))
                print(f"  {S.label} {lbl:34s} {bl:7s} n={sel.sum():4d} rms {np.std(xm[sel]):.4f} | slope {b1[0]:5.2f}±{e1[0]:4.2f} (if GNILC true {b3[0]*lam:5.2f}; pure extinction {ext:4.2f}) | GNILC-SFD alone {b3[0]:5.2f}±{e3[0]:4.2f} | joint {b2[0]:5.2f}±{e2[0]:4.2f}, {b2[1]:5.2f}±{e2[1]:4.2f} | corr {np.corrcoef(xm[sel], xg[sel])[0,1]:.2f}")

main(*sys.argv[1:7])
