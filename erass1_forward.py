"""Forward model for the blind referee: does an isotropic universe, observed through eRASS1's exposure map with a
detection limit in counts, produce the Migkas-cone L_X-T deficit of the observed size?
Model: take the real catalogue's sky positions and exposures (all 12,247 clusters: the sky sampling of the survey), draw
T from the observed T distribution and z from the observed z distribution, draw L from the isotropic fitted L-T relation with
the observed 0.47 dex scatter, convert to flux (L/4 pi D_L^2), to counts = flux * EXP * k (k fixed so the detected
fraction matches the real T-sample size when the detection threshold is the real sample's minimum CTS500), keep the
detected, fit the relation to them, measure the cone deficit about (280,-15). Repeat; compare with the observed -0.039 dex.
Also: the same with a threshold on counts for *temperature measurement* (the T sample needs more counts than detection).
usage: python erass1_forward.py data/erass1/erass1cl_primary_v3.2.fits
"""
import sys, numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
rng = np.random.default_rng(5); OM = 0.3; cosmo = FlatLambdaCDM(H0=70, Om0=OM)
def Ez(z): return np.sqrt(OM*(1+z)**3 + 1-OM)

def main(f):
    t = fits.open(f)[1].data
    allok = np.isfinite(t['EXP']) & (t['EXP'] > 0) & np.isfinite(t['BEST_Z']) & (t['BEST_Z'] > 0.05) & (t['BEST_Z'] < 0.6)
    sc = SkyCoord(ra=t['RA'][allok]*u.deg, dec=t['DEC'][allok]*u.deg).galactic; l, b = sc.l.deg, sc.b.deg
    n_all = np.stack([np.cos(np.radians(b))*np.cos(np.radians(l)), np.cos(np.radians(b))*np.sin(np.radians(l)), np.sin(np.radians(b))], -1)
    exp_all = t['EXP'][allok]; z_all = t['BEST_Z'][allok]
    axis = np.array([np.cos(np.radians(-15))*np.cos(np.radians(280)), np.cos(np.radians(-15))*np.sin(np.radians(280)), np.sin(np.radians(-15))])
    # the real T sample, for the T distribution, the relation and the thresholds
    kt, ktl, kth = t['KT'], t['KT_L'], t['KT_H']; L, Ll, Lh = t['L500'], t['L500_L'], t['L500_H']; z = t['BEST_Z']
    ekt = 0.5*(np.abs(kth - kt) + np.abs(kt - ktl)); eL = 0.5*(np.abs(Lh - L) + np.abs(L - Ll))
    okT = allok & np.isfinite(kt) & (kt > 0) & np.isfinite(L) & (L > 0) & np.isfinite(ekt) & np.isfinite(eL) & (ekt/kt < 0.5)
    cts_T = t['CTS500'][okT]; thr_T = np.nanquantile(cts_T, 0.02)
    x = np.log10(kt[okT]/4.0); y = np.log10(L[okT]/1e44) - np.log10(Ez(z[okT]))
    X = np.column_stack([np.ones_like(x), x]); beta = np.linalg.lstsq(X, y, rcond=None)[0]; scat = np.std(y - X@beta)
    print(f"real T sample: {okT.sum()} clusters, relation logL = {beta[0]:.3f} + {beta[1]:.2f} logT, scatter {scat:.3f} dex; CTS500 2nd percentile {thr_T:.0f}; sky sampling {allok.sum()} clusters")
    # observed cone deficit for reference (same estimator as below)
    def cone_deficit(xs, ys, nn):
        Xs = np.column_stack([np.ones_like(xs), xs]); bb = np.linalg.lstsq(Xs, ys, rcond=None)[0]; r = ys - Xs@bb
        inc = nn@axis >= np.cos(np.radians(75)); return r[inc].mean() - r.mean(), inc.sum()
    n_T = np.stack([np.cos(np.radians(SkyCoord(ra=t['RA'][okT]*u.deg, dec=t['DEC'][okT]*u.deg).galactic.b.deg))*np.cos(np.radians(SkyCoord(ra=t['RA'][okT]*u.deg, dec=t['DEC'][okT]*u.deg).galactic.l.deg)),
                    np.cos(np.radians(SkyCoord(ra=t['RA'][okT]*u.deg, dec=t['DEC'][okT]*u.deg).galactic.b.deg))*np.sin(np.radians(SkyCoord(ra=t['RA'][okT]*u.deg, dec=t['DEC'][okT]*u.deg).galactic.l.deg)),
                    np.sin(np.radians(SkyCoord(ra=t['RA'][okT]*u.deg, dec=t['DEC'][okT]*u.deg).galactic.b.deg))], -1)
    obs, nin = cone_deficit(x, y, n_T); print(f"observed cone deficit (unweighted estimator): {obs:+.4f} dex, {nin} in cone")
    # counts-to-flux conversion: calibrate k so that counts = k * flux * EXP reproduces the real sample's median counts
    DL = cosmo.luminosity_distance(z[okT]).to(u.cm).value; flux = L[okT]/(4*np.pi*DL**2)
    k = np.nanmedian(cts_T/(flux*t['EXP'][okT])); print(f"counts conversion k = {k:.3e} counts per (erg/s/cm2 * s)")
    # forward model
    Tpool = kt[okT]; nmock = 400
    out = []; sizes = []
    for i in range(nmock):
        # a large isotropic population on the survey's sky sampling: each real position gets several draws
        rep = 6
        pos = np.repeat(np.arange(allok.sum()), rep); nn = n_all[pos]; ee = exp_all[pos]; zz = z_all[pos]
        TT = rng.choice(Tpool, len(pos)); xx = np.log10(TT/4.0)
        yy = beta[0] + beta[1]*xx + rng.normal(0, scat, len(pos))          # isotropic relation
        LL = 10**yy*1e44*Ez(zz); DLm = cosmo.luminosity_distance(zz).to(u.cm).value; fl = LL/(4*np.pi*DLm**2)
        cts = fl*ee*k
        det = cts >= thr_T
        # keep a T-measurable sample of the real size: top by counts above the threshold would over-select; use the threshold only
        sel = np.where(det)[0]
        if len(sel) > 3*okT.sum(): sel = rng.choice(sel, 3*okT.sum(), replace=False)
        cd, _ = cone_deficit(xx[sel], yy[sel], nn[sel]); out.append(cd); sizes.append(len(sel))
    out = np.array(out)
    print(f"forward model (isotropic L-T, detection at CTS500 >= {thr_T:.0f} with the real exposure map, {nmock} mocks, ~{int(np.mean(sizes))} detected per mock): cone deficit {out.mean():+.4f} ± {out.std():.4f} dex  (observed {obs:+.4f})")
    # exposure dependence in the mocks, for comparison with the observed Spearman
    from scipy.stats import spearmanr
    rep = 6; pos = np.repeat(np.arange(allok.sum()), rep); nn = n_all[pos]; ee = exp_all[pos]; zz = z_all[pos]
    TT = rng.choice(Tpool, len(pos)); xx = np.log10(TT/4.0); yy = beta[0] + beta[1]*xx + rng.normal(0, scat, len(pos))
    LL = 10**yy*1e44*Ez(zz); DLm = cosmo.luminosity_distance(zz).to(u.cm).value; fl = LL/(4*np.pi*DLm**2); cts = fl*ee*k; sel = cts >= thr_T
    Xs = np.column_stack([np.ones(sel.sum()), xx[sel]]); bb = np.linalg.lstsq(Xs, yy[sel], rcond=None)[0]; rr = yy[sel] - Xs@bb
    print(f"  in one mock: residual vs log EXP Spearman {spearmanr(rr, np.log10(ee[sel]))[0]:+.3f}; residual vs log counts Spearman {spearmanr(rr, np.log10(cts[sel]))[0]:+.3f}  (observed: -0.36 and -0.20)")
    # higher threshold (T needs more counts than detection): scan thresholds
    for q in (0.05, 0.1, 0.25):
        thr = np.nanquantile(cts_T, q); o = []
        for i in range(100):
            pos = np.repeat(np.arange(allok.sum()), 6); nn = n_all[pos]; ee = exp_all[pos]; zz = z_all[pos]
            TT = rng.choice(Tpool, len(pos)); xx = np.log10(TT/4.0); yy = beta[0] + beta[1]*xx + rng.normal(0, scat, len(pos))
            LL = 10**yy*1e44*Ez(zz); DLm = cosmo.luminosity_distance(zz).to(u.cm).value; fl = LL/(4*np.pi*DLm**2); cts = fl*ee*k; sel = cts >= thr
            o.append(cone_deficit(xx[sel], yy[sel], nn[sel])[0])
        print(f"  threshold at the real sample's {int(q*100)}th percentile of counts ({thr:.0f}): cone deficit {np.mean(o):+.4f} ± {np.std(o):.4f} dex")

main(sys.argv[1])
