"""R-156: null tests and contamination bounds (after r156_run.py, which writes the templates to run.json).

Substitute marks X (standardised: residual of a cubic in z per cap, unit variance): ln PSFSIZE_R; a stars' PSF map
(HEALPix nside 128 mean FIBERFLUX_R/FLUX_R of Gaia-confirmed stars); FRACFLUX_R; ln GALDEPTH_R; ln PSFDEPTH_R;
E(B-V); stellar density (nside 128 star counts). For each: the dipole per bin against zero (jackknife covariance, since
the maps are spatially coherent), its amplitude A(X) on the bin's Doppler template, the regression slope beta_X of
kappa_hat on X, and the bound dA_X = beta_X A(X) on A_D.
Splits: NGC vs SGC; DR9 north (RELEASE 9011, BASS/MzLS) vs south; FRACFLUX_R < 0.05. Randomised redshifts: the
`zshuf` pair counts with the real kappa_hat.

usage: python -u r156_nulls.py
"""
import json, os
import numpy as np
import healpy as hp
import glob
import r156_measure as M
import r156_run as RU

D = 'data/r156/'
FIT = RU.FIT


def star_maps():
    parts = [dict(np.load(f)) for f in sorted(glob.glob(D + 'tractor/*.npz'))]
    st = {k: np.concatenate([p[k][~p['ISBGS']] for p in parts if len(p['TARGETID'])]) for k in
          ['RA', 'DEC', 'FIBERFLUX_R', 'FLUX_R', 'GAIA_PHOT_G_MEAN_MAG']}
    ok = (st['FLUX_R'] > 0) & (st['FIBERFLUX_R'] > 0) & (st['GAIA_PHOT_G_MEAN_MAG'] > 15) & (st['GAIA_PHOT_G_MEAN_MAG'] < 20)
    pix = hp.ang2pix(128, st['RA'][ok], st['DEC'][ok], lonlat=True)
    n = np.bincount(pix, minlength=hp.nside2npix(128)).astype(float)
    ff = np.bincount(pix, weights=(st['FIBERFLUX_R'] / st['FLUX_R'])[ok], minlength=hp.nside2npix(128))
    with np.errstate(invalid='ignore'):
        ffm = ff / n
    return ffm, n, ok.sum()


def main():
    run = json.load(open(D + 'run.json'))
    s, sel, k, u, rel = M.prepare(False)
    sh = M.Shells('data')
    kid = sh.kidx
    gj = sel[kid]
    zj = s['Z'][kid]
    reg = M.jk_regions(s['RA'][kid], s['DEC'][kid], np.ones(len(kid), bool))
    order = np.argsort(reg, kind='stable')
    bounds = np.searchsorted(reg[order], np.arange(M.NJK + 1))
    ffm, nst, nstar = star_maps()
    pix = hp.ang2pix(128, s['RA'], s['DEC'], lonlat=True)
    ffg = ffm[pix]
    ffg = np.where(np.isfinite(ffg), ffg, np.nanmedian(ffm))
    ebv = -2.5 * np.log10(np.clip(s['MW_TRANSMISSION_R'], 1e-3, None)) / 2.165
    marks = {'ln PSFSIZE_R': np.log(np.clip(s['PSFSIZE_R'], 0.3, None)), 'stars PSF map': ffg,
             'FRACFLUX_R': np.clip(s['FRACFLUX_R'], 0, 1), 'ln GALDEPTH_R': np.log(np.clip(s['GALDEPTH_R'], 1, None)),
             'ln PSFDEPTH_R': np.log(np.clip(s['PSFDEPTH_R'], 1, None)), 'E(B-V)': ebv,
             'stellar density': nst[pix]}
    out = {'n_stars': int(nstar), 'marks': {}, 'splits': {}}
    for name, x in marks.items():
        X = M.standardise(x, s, sel)
        Xj = np.where(gj, X[kid], 0.0)
        out['marks'][name] = {}
        for bi, (lo, hi) in enumerate(M.BINS):
            R = run['bins'][str(bi)]
            B = gj & (zj > lo) & (zj < hi)
            w = np.where(B, u[kid], 0.0)
            num, den = RU.region_sums(sh, order, bounds, w * Xj, w)
            xi, jk = sh.multipoles(num, den, 1)
            C = M.jk_cov(jk)[FIT, FIT]
            tD = np.array(R['tD'])
            A, Ae, _, chi0 = M.gls(xi[FIT], tD[FIT], C, 13)
            m = sel & (s['Z'] > lo) & (s['Z'] < hi)
            beta = np.polyfit(X[m], k[m], 1)[0]
            from scipy.stats import chi2 as C2
            p = 1 - C2.cdf(chi0, 13)
            out['marks'][name][bi] = {'A_X': A, 'A_X_err': Ae, 'chi2': chi0, 'p': p, 'beta': beta, 'dA': beta * A,
                                      'dA_err': abs(beta) * Ae}
            print('%-16s bin %d: chi2 %.1f/13 p %.3f; A(X) %.2f +- %.2f; beta %.4f; bound dA_D = %.3f +- %.3f' % (
                name, bi, chi0, p, A, Ae, beta, beta * A, abs(beta) * Ae), flush=True)

    kj = np.where(gj, k[kid], 0.0)
    rel9011 = s['RELEASE'] == 9011
    splits = {'NGC': s['CAP'] == 0, 'SGC': s['CAP'] == 1, 'north (BASS/MzLS)': rel9011, 'south (DECaLS)': ~rel9011,
              'FRACFLUX_R<0.05': s['FRACFLUX_R'] < 0.05}
    for name, msk in splits.items():
        mj = msk[kid]
        out['splits'][name] = {}
        for bi, (lo, hi) in enumerate(M.BINS[:2]):
            R = run['bins'][str(bi)]
            B = gj & mj & (zj > lo) & (zj < hi)
            w = np.where(B, u[kid], 0.0)
            num, den = RU.region_sums(sh, order, bounds, w * kj, w)
            xi, jk = sh.multipoles(num, den, 1)
            y = xi - np.array(R['tRel'])
            C = M.jk_cov(jk)[FIT, FIT]  # the mock-validated jackknife covariance (results F4)
            A, Ae, _, _ = M.gls(y[FIT], np.array(R['tD'])[FIT], C, 13)
            out['splits'][name][bi] = {'A_D': A, 'err': Ae, 'N': int(B.sum())}
            print('split %-18s bin %d: N %d A_D = %.2f +- %.2f' % (name, bi, B.sum(), A, Ae), flush=True)

    if os.path.exists(D + 'pairs_zshuf.npz'):
        zs = M.Shells('zshuf')
        zid = zs.kidx
        g2 = sel[zid]
        zz = zs.zk
        r2 = M.jk_regions(s['RA'][zid], s['DEC'][zid], np.ones(len(zid), bool))
        o2 = np.argsort(r2, kind='stable')
        b2 = np.searchsorted(r2[o2], np.arange(M.NJK + 1))
        out['zshuffle'] = {}
        for bi, (lo, hi) in enumerate(M.BINS):
            R = run['bins'][str(bi)]
            B = g2 & (zz > lo) & (zz < hi)
            w = np.where(B, u[zid], 0.0)
            num, den = RU.region_sums(zs, o2, b2, w * np.where(g2, k[zid], 0.0), w)
            xi, jk = zs.multipoles(num, den, 1)
            C = M.jk_cov(jk)[FIT, FIT]
            A, Ae, _, chi0 = M.gls(xi[FIT], np.array(R['tD'])[FIT], C, 13)
            from scipy.stats import chi2 as C2
            out['zshuffle'][bi] = {'A': A, 'err': Ae, 'chi2': chi0, 'p': 1 - C2.cdf(chi0, 13)}
            print('randomised z bin %d: chi2 %.1f/13 p %.3f; A = %.2f +- %.2f' % (bi, chi0, 1 - C2.cdf(chi0, 13), A, Ae),
                  flush=True)
    json.dump(out, open(D + 'nulls.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
