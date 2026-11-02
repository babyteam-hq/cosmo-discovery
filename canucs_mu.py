"""R-128 secondary: the same regressions (R1-R3 of jwst_mu_test.py) on the CANUCS DR1 cluster fields, with their own Lenstool
models (independent modellers: Gledhill, Rihtarsic, Desprez) and Dense Basis photo-z. mu at each source from the best-fit kappa and
gamma maps (Dds/Ds = 1) scaled to the source redshift; the per-source lens spread from the 100 MCMC sample maps (600 mas) where
released (A370, MACS0416, MACS1149). The catalogue has no fluxes: the UV magnitude is rest_FUV (apparent-frame absolute magnitude,
not de-magnified) and the S/N cut is replaced by an apparent UV magnitude cut m_UV < 29 (labelled difference from the primary).
Usage: python canucs_mu.py [--mcut 29]
"""
import os, sys, json, argparse
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jwst_mu_test import COSMO, D, OUT, estimators, ols, split_normal, dm, rng, MU_MIN, MU_MAX, ZMIN, Z16MIN

ZL = {'a370': 0.375, 'macs0416': 0.396, 'macs1149': 0.543, 'macs0417': 0.443, 'macs1423': 0.545}
C = os.path.join(D, 'canucs')


def f_dds_ds(zl, zs):
    """D_ds / D_s for sources zs behind a lens at zl (maps are for Dds/Ds = 1)."""
    zs = np.atleast_1d(zs)
    dds = COSMO.angular_diameter_distance_z1z2(zl, zs).value; ds = COSMO.angular_diameter_distance(zs).value
    return dds / ds


def mu_from(kappa, gamma, f):
    k = kappa * f; g = gamma * f
    return 1.0 / np.abs((1 - k) ** 2 - g ** 2)


def sample_map(path, ra, dec, ext=0):
    h = fits.open(path)
    w = WCS(h[ext].header); data = h[ext].data
    x, y = w.all_world2pix(ra, dec, 0); xi = np.round(x).astype(int); yi = np.round(y).astype(int)
    ok = (xi >= 0) & (yi >= 0) & (xi < data.shape[1]) & (yi < data.shape[0])
    out = np.full(len(ra), np.nan); out[ok] = data[yi[ok], xi[ok]]
    return out, h


def build(cl, mcut=29.0):
    t = Table.read(os.path.join(C, f'hlsp_canucs_jwst-hst_multi_{cl}-clu_multi_v1_densebasis-cat.fits.gz'))
    z = np.array(t['redshift_50_COLOR03'], float); z16 = np.array(t['redshift_16_COLOR03'], float); z84 = np.array(t['redshift_84_COLOR03'], float)
    Mobs = np.array(t['rest_FUV_COLOR03'], float)
    ok = (np.array(t['USE_PHOT']) == 1) & (np.array(t['FLAG_DBFITS_COLOR03']) == 0) & np.isfinite(z) & np.isfinite(Mobs)
    mab = Mobs + dm(np.where(np.isfinite(z), z, 1)) - 2.5 * np.log10(1 + np.where(np.isfinite(z), z, 1))
    sel = ok & (z > ZMIN) & (z16 > Z16MIN) & (mab < mcut)
    ra = np.array(t['RA'], float)[sel]; dec = np.array(t['DEC'], float)[sel]
    kap, _ = sample_map(os.path.join(C, f'hlsp_canucs_jwst-hst_multi_{cl}-kappa-dlsds1-best-100mas_multi_v1_model.fits.gz'), ra, dec)
    gam, _ = sample_map(os.path.join(C, f'hlsp_canucs_jwst-hst_multi_{cl}-gamma-dlsds1-best-100mas_multi_v1_model.fits.gz'), ra, dec)
    f = f_dds_ds(ZL[cl], z[sel])
    mu = mu_from(kap, gam, f)
    out = Table()
    out['ra'] = ra; out['dec'] = dec; out['z'] = z[sel]; out['z16'] = z16[sel]; out['z84'] = z84[sel]
    out['mab'] = mab[sel]; out['snr'] = np.full(sel.sum(), 99.0)
    out['mu'] = mu
    out['lmass'] = np.array(t['log_stellar_mass_50_COLOR03'], float)[sel] - np.log10(mu)   # the catalogue mass is not de-magnified
    out['lmass16'] = np.array(t['log_stellar_mass_16_COLOR03'], float)[sel] - np.log10(mu); out['lmass84'] = np.array(t['log_stellar_mass_84_COLOR03'], float)[sel] - np.log10(mu)
    out['fr'] = np.array(t['FLUX_RADIUS'], float)[sel]
    # lens spread from the sample maps
    ks = os.path.join(C, f'hlsp_canucs_jwst-hst_multi_{cl}-kappa-dlsds1-samples-600mas_multi_v1_model.fits')   # gunzipped: a 100-extension gz re-reads from the start per extension
    gs = os.path.join(C, f'hlsp_canucs_jwst-hst_multi_{cl}-gamma-dlsds1-samples-600mas_multi_v1_model.fits')
    if os.path.exists(ks) and os.path.exists(gs):
        hk = fits.open(ks); hg = fits.open(gs)
        w = WCS(hk[1].header); x, y = w.all_world2pix(ra, dec, 0); xi = np.round(x).astype(int); yi = np.round(y).astype(int)
        inside = (xi >= 0) & (yi >= 0) & (xi < hk[1].data.shape[1]) & (yi < hk[1].data.shape[0])
        mus = np.full((100, len(ra)), np.nan)
        for i in range(1, 101):
            kk = np.full(len(ra), np.nan); gg = np.full(len(ra), np.nan)
            kk[inside] = hk[i].data[yi[inside], xi[inside]]; gg[inside] = hg[i].data[yi[inside], xi[inside]]
            mus[i - 1] = mu_from(kk, gg, f)
        lm = np.log10(mus)
        out['mu16'] = 10 ** np.nanpercentile(lm, 16, axis=0); out['mu84'] = 10 ** np.nanpercentile(lm, 84, axis=0)
        # centre the spread on the best model (the samples' median can differ from the best fit)
        med = 10 ** np.nanmedian(lm, axis=0)
        out['mu16'] = out['mu16'] * mu / med; out['mu84'] = out['mu84'] * mu / med
        out['mu_med_samples'] = med
    else:
        out['mu16'] = mu / 10 ** 0.1; out['mu84'] = mu * 10 ** 0.1    # no samples released: a nominal 0.1 dex
    out['in_mu'] = np.isfinite(mu) & (mu >= MU_MIN) & (mu <= MU_MAX)
    out['quad'] = (ra > np.median(ra)).astype(int) + 2 * (dec > np.median(dec)).astype(int)
    out['field'] = cl
    out.write(os.path.join(OUT, f'canucs_{cl}.fits'), overwrite=True)
    w = np.array(out['in_mu'])
    print(f'{cl}: {sel.sum()} z>8 (m_UV<{mcut}), {w.sum()} in 1.2<=mu<=30 (median mu {np.nanmedian(mu[w]):.2f}; '
          f'median 68% half-width {np.nanmedian(0.5 * np.log10(out["mu84"][w] / out["mu16"][w])):.3f} dex)')
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--mcut', type=float, default=29.0); a = ap.parse_args()
    from jwst_mu_test import fit
    tabs = []
    for cl in ZL:
        if os.path.exists(os.path.join(C, f'hlsp_canucs_jwst-hst_multi_{cl}-kappa-dlsds1-best-100mas_multi_v1_model.fits.gz')):
            tabs.append(build(cl, a.mcut))
    from astropy.table import vstack
    allt = vstack(tabs, metadata_conflicts='silent')
    allt.write(os.path.join(OUT, 'canucs_all.fits'), overwrite=True)
    for cl in ZL:
        if os.path.exists(os.path.join(OUT, f'canucs_{cl}.fits')):
            try:
                fit(snr_min=0, table=f'canucs_{cl}', nboot=300, ndraw=300)
            except Exception as e:
                print(cl, 'fit failed:', e)
    res = fit(snr_min=0, table='canucs_all', nboot=1000, ndraw=500)
    # jackknife over fields
    t = allt[np.array(allt['in_mu'])]
    fields = np.array(t['field'])
    for k in ['R1', 'R2', 'R3']:
        jk = []
        for fl in np.unique(fields):
            w = fields != fl
            jk.append(estimators(np.array(t['z'])[w], np.array(t['mab'])[w], np.array(t['mu'])[w], np.array(t['lmass'], float)[w], np.array(t['fr'], float)[w])[k])
        jk = np.array(jk, float); jk = jk[np.isfinite(jk)]
        sd = np.sqrt((len(jk) - 1) / len(jk) * np.sum((jk - jk.mean()) ** 2))
        print(f'{k} jackknife over fields: {np.round(jk, 3).tolist()}, sd {sd:.3f}')


if __name__ == '__main__':
    main()
