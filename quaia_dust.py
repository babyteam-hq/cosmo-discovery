"""Quaia (Storey-Fisher+2024, G < 20.5, 1.3 M quasars) number counts vs the dust-map difference (sweep-2.md row 21): the
Quaia selection function (Gaussian process on nside-64 pixels) takes as its dust feature the CSFD map x 0.86 (Chiang 2023,
SFD with the CIB removed), so the residual counts N / (selection x mean) - 1 should carry no CSFD dependence; if the true
extinction is what GNILC (or a stellar map) says, the residual tracks (GNILC - CSFD) with the count slope -ln(10) s R_G,
where s = d log10 N / dG at the faint limit (measured from the catalogue) and R_G = A_G / E(B-V) ~ 2.7. Fit: weighted least
squares of the residual on [1, dipole (3), template], weights = selection (Poisson), pixels with selection > 0.5 and
|b| > 30 (and > 20), LMC/SMC excluded. Templates: GNILC - CSFD, GNILC - SFD, CSFD - SFD, SFD, Zhou+2024 DESI stellar map
minus its SFD (where covered). The fitted dipole is printed for reference (Mittal+2024 found it consistent with the CMB
one at this cut). Also a 3-bin redshift tomography of the dipole with the same selection function (a caveat: the per-bin
selection functions are not public).
usage: python quaia_dust.py data/quaia/quaia_G20.5.fits data/quaia/selection_function_NSIDE64_G20.5.fits data/dustmaps data/zhou/desi_dust_gr_256.fits data/mudur/recon_fw15_final_mult.fits
"""
import sys, numpy as np, healpy as hp
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
NS = 64; R_G = 2.7

def templates(data_dir, zhou, mudur):
    from dustmaps.config import config; config['data_dir'] = data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.planck import PlanckGNILCQuery
    npix = hp.nside2npix(NS); ra, dec = hp.pix2ang(NS, np.arange(npix), lonlat=True)
    # average each map over 16 sub-pixels (nside 256) of every nside-64 pixel
    sub = np.arange(hp.nside2npix(256)); sra, sdec = hp.pix2ang(256, sub, lonlat=True); parent = hp.ang2pix(NS, sra, sdec, lonlat=True)
    c = SkyCoord(sra*u.deg, sdec*u.deg, frame='icrs')
    def avg(v): return np.bincount(parent, weights=np.nan_to_num(v), minlength=npix)/np.maximum(np.bincount(parent, weights=np.isfinite(v).astype(float), minlength=npix), 1)
    out = {'SFD': avg(SFDQuery()(c)), 'GNILC': avg(PlanckGNILCQuery()(c))}
    try:
        from dustmaps.csfd import CSFDQuery; out['CSFD'] = avg(CSFDQuery()(c))
    except Exception as e: print("no CSFD:", repr(e)[:120])
    try:
        t = fits.getdata(zhou, 1); ok = t['N_STAR'] > 0; ze = np.full(hp.nside2npix(256), np.nan); zs = np.full(hp.nside2npix(256), np.nan)
        ze[t['HPXPIXEL'][ok]] = t['EBV_GR'][ok]; zs[t['HPXPIXEL'][ok]] = t['EBV_SFD'][ok]
        cnt = np.bincount(parent, weights=np.isfinite(ze).astype(float), minlength=npix); out['Zhou'] = np.where(cnt >= 12, avg(ze), np.nan); out['ZhouSFD'] = np.where(cnt >= 12, avg(zs), np.nan)
    except Exception as e: print("no Zhou:", repr(e)[:120])
    try:
        mu = fits.getdata(mudur, 1)['Recon_Mean'].ravel().astype(float); mu[mu < -1] = np.nan; sg = c.galactic
        mv = mu[hp.ang2pix(2048, sg.l.deg, sg.b.deg, lonlat=True)]; cnt = np.bincount(parent, weights=np.isfinite(mv).astype(float), minlength=npix)
        out['Mudur'] = np.where(cnt >= 12, avg(mv), np.nan)
    except Exception as e: print("no Mudur:", repr(e)[:120])
    gal = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs').galactic
    return out, gal.b.deg, gal.l.deg

def main(cat, self_, data_dir, zhou, mudur):
    d = fits.getdata(cat, 1); sel = hp.read_map(self_, verbose=False) if False else fits.getdata(self_, 1)['T'].ravel()
    npix = hp.nside2npix(NS); pix = hp.ang2pix(NS, d['ra'], d['dec'], lonlat=True); N = np.bincount(pix, minlength=npix).astype(float)
    G = d['phot_g_mean_mag']; h, e = np.histogram(G, bins=np.arange(19.0, 20.51, 0.25)); s_cnt = np.polyfit(0.5*(e[1:-1]+e[2:]), np.log10(h[1:]), 1)[0]
    print(f"Quaia G<20.5: {len(d)} sources; faint-end count slope d log10 N / dG = {s_cnt:.2f} (19.25-20.5); expected residual slope per unit E(B-V) error = -ln10 s R_G = {-np.log(10)*s_cnt*R_G:.2f}")
    tm, b, l = templates(data_dir, zhou, mudur)
    vec = np.array(hp.pix2vec(NS, np.arange(npix))).T
    lmc = hp.query_disc(NS, hp.ang2vec(280.47, -32.89, lonlat=True), np.radians(8)); smc = hp.query_disc(NS, hp.ang2vec(302.8, -44.3, lonlat=True), np.radians(5))
    bad = np.zeros(npix, bool); bad[lmc] = True; bad[smc] = True
    for bcut in (30, 20):
        m = (sel > 0.5) & (np.abs(b) > bcut) & ~bad & (N > 0)
        nbar = np.sum(N[m])/np.sum(sel[m]); r = N[m]/(sel[m]*nbar) - 1; w = sel[m]
        def fit(cols, label):
            X = np.column_stack([np.ones(m.sum()), vec[m, 0], vec[m, 1], vec[m, 2]] + cols); sw = np.sqrt(w)
            bb, *_ = np.linalg.lstsq(X*sw[:, None], r*sw, rcond=None); res = r - X@bb; chi = np.sum(res**2*w)/(m.sum()-X.shape[1])
            cov = np.linalg.inv((X*w[:, None]).T@X)*chi; ee = np.sqrt(np.diag(cov))
            D = bb[1:4]; De = ee[1:4]; amp = np.linalg.norm(D); ra_, dec_ = hp.vec2ang(D, lonlat=True)
            # dipole direction in galactic
            gd = SkyCoord(ra_[0]*u.deg, dec_[0]*u.deg, frame='icrs').galactic
            print(f"  |b|>{bcut} n={m.sum():5d} {label:28s} " + " ".join(f"t{j}={bb[4+j]:+6.2f}±{ee[4+j]:4.2f}" for j in range(len(cols))) + f" | dipole {amp:.4f} toward (l,b)=({gd.l.deg:.0f},{gd.b.deg:.0f}) (CMB: 264,48; kinematic ~0.007) | rms resid {np.std(res):.3f}")
            return bb[4:], ee[4:]
        fit([], "dipole only")
        s_, g_ = tm['SFD'][m], tm['GNILC'][m]
        fit([s_], "SFD"); fit([g_ - s_], "GNILC - SFD")
        if 'CSFD' in tm:
            c_ = tm['CSFD'][m]; fit([c_], "CSFD (the GP's own feature)"); fit([g_ - c_], "GNILC - CSFD"); fit([c_ - s_], "CSFD - SFD"); fit([g_ - c_, c_], "GNILC - CSFD, with CSFD")
        if 'Mudur' in tm:
            mu_ = tm['Mudur'][m]; okm = np.isfinite(mu_); k = np.sum(s_[okm]*mu_[okm])/np.sum(mu_[okm]**2)
            fit([np.where(okm, k*mu_ - s_, 0.0)], f"k Mudur - SFD (k={k:.2f}, {okm.sum()} px)")
        if 'Zhou' in tm:
            z_ = tm['Zhou'][m] - tm['ZhouSFD'][m]; ok = np.isfinite(z_)
            mm = m.copy(); mm[np.where(m)[0][~ok]] = False
            nbar2 = np.sum(N[mm])/np.sum(sel[mm]); r2 = N[mm]/(sel[mm]*nbar2) - 1; w2 = sel[mm]
            X = np.column_stack([np.ones(mm.sum()), vec[mm, 0], vec[mm, 1], vec[mm, 2], (tm['Zhou'] - tm['ZhouSFD'])[mm], (tm['GNILC'] - tm['SFD'])[mm]]); sw = np.sqrt(w2)
            bb, *_ = np.linalg.lstsq(X*sw[:, None], r2*sw, rcond=None); res = r2 - X@bb; chi = np.sum(res**2*w2)/(mm.sum()-6); ee = np.sqrt(np.diag(np.linalg.inv((X*w2[:, None]).T@X)*chi))
            print(f"  |b|>{bcut} n={mm.sum():5d} Zhou-SFD and GNILC-SFD jointly (DESI footprint): {bb[4]:+.2f}±{ee[4]:.2f}, {bb[5]:+.2f}±{ee[5]:.2f}")
            X = X[:, :5]; bb, *_ = np.linalg.lstsq(X*sw[:, None], r2*sw, rcond=None); res = r2 - X@bb; chi = np.sum(res**2*w2)/(mm.sum()-5); ee = np.sqrt(np.diag(np.linalg.inv((X*w2[:, None]).T@X)*chi))
            print(f"  |b|>{bcut} n={mm.sum():5d} Zhou-SFD alone: {bb[4]:+.2f}±{ee[4]:.2f}")
    # tomography with the same selection function
    z = d['redshift_quaia']; edges = np.quantile(z, [0, 1/3, 2/3, 1]); m = (sel > 0.5) & (np.abs(b) > 30) & ~bad
    print("dipole by redshift tercile (same selection function; its z-dependence is not modelled, so read amplitudes relative to each other):")
    for i in range(3):
        zz = (z >= edges[i]) & (z < edges[i+1] + (1e-9 if i == 2 else 0)); Nz = np.bincount(pix[zz], minlength=npix).astype(float)
        nbar = np.sum(Nz[m])/np.sum(sel[m]); r = Nz[m]/(sel[m]*nbar) - 1; w = sel[m]
        X = np.column_stack([np.ones(m.sum()), vec[m, 0], vec[m, 1], vec[m, 2]]); sw = np.sqrt(w); bb, *_ = np.linalg.lstsq(X*sw[:, None], r*sw, rcond=None)
        res = r - X@bb; chi = np.sum(res**2*w)/(m.sum()-4); ee = np.sqrt(np.diag(np.linalg.inv((X*w[:, None]).T@X)*chi)); D = bb[1:4]; gd = SkyCoord(*[x[0] for x in hp.vec2ang(D, lonlat=True)]*u.deg, frame='icrs').galactic
        print(f"  {edges[i]:.2f}<z<{edges[i+1]:.2f} n={zz.sum():7d}: dipole {np.linalg.norm(D):.4f} ± {np.mean(ee[1:4]):.4f} toward (l,b)=({gd.l.deg:.0f},{gd.b.deg:.0f})")

main(*sys.argv[1:6])
