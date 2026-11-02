"""R-148 step 0: the feasibility forecast for DESI DR1 Lya forest x ACT DR6 Compton-y.

  python -u lya_y_forecast.py      # needs data/lya/forests_lya.npz (lya_y_prep.py) and data/act_y/{y,y_cib,wide_mask}.fits

1. Footprint: forests whose sightline sits where the ACT wide GAL070 mask is > 0.9.
2. Noise of the annulus-mean y around a sightline, sigma_y(theta), on the standard and CIB-deprojected maps, measured
   at forest positions rotated by +97 deg in RA (a null placement; no forest-y pairing is looked at).
3. The estimator of the brief: xi_b(theta) = sum_q A_qb ybar_q(theta) / sum_q W_qb. Its expected value under the
   model is sum_q S_qb(theta) / sum W (after the continuum projection) or sum U_qb / sum W (if the continuum were
   known). Its noise, sightlines independent (optimistic at theta > 5' where neighbouring annuli overlap):
   sqrt(sum_q A_qb^2) sigma_y(theta) / sum W.
4. A third, hypothetical case: a continuum known to a per-sightline amplitude error sigma_c (Croft+2005 used
   sigma(tau_eff) = 0.1; a red-side PCA continuum is ~5-10 %), no projection: noise sqrt(sum_q (sigma_c W_qb)^2 + A^2).
"""
import numpy as np, os
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, 'data')
EDGES = np.load(os.path.join(D, 'lya/model.npz'))['theta_edges']
PIX = 0.5  # arcmin
NY, NX = 10320, 43200
CRPIX1, CRPIX2 = 21601.0, 7561.0


def radec_to_pix(ra, dec):
    x = CRPIX1 - 1 - ra / (PIX / 60)          # CDELT1 negative, CRVAL1 = 0
    x = np.mod(np.round(x), NX).astype(int)
    y = np.round(CRPIX2 - 1 + dec / (PIX / 60)).astype(int)
    return y, x


def annulus_means(m, ra, dec, edges=EDGES):
    """Mean of map m in annuli around each position (flat-sky offsets, CAR cos(dec) in RA). Python loop."""
    out = np.full((len(ra), len(edges) - 1), np.nan)
    R = edges[-1]
    nr = int(np.ceil(R / PIX)) + 1
    dyg, dxg = np.mgrid[-nr:nr + 1, -nr:nr + 1]
    jy, jx = radec_to_pix(ra, dec)
    for k in range(len(ra)):
        c = np.cos(np.radians(dec[k]))
        nxr = int(np.ceil(nr / max(c, 0.2)))
        yy = jy[k] + dyg[:, 0][:, None]
        xs = np.arange(-nxr, nxr + 1)
        rr = np.hypot(dyg[:, 0][:, None] * PIX, xs[None, :] * PIX * c)
        if yy.min() < 0 or yy.max() >= NY:
            continue
        sub = m[yy[:, 0].min():yy[:, 0].max() + 1][:, np.mod(jx[k] + xs, NX)]
        ib = np.digitize(rr, edges) - 1
        ok = (ib >= 0) & (ib < len(edges) - 1)
        s = np.bincount(ib[ok], weights=sub[ok], minlength=len(edges) - 1)
        n = np.bincount(ib[ok], minlength=len(edges) - 1)
        out[k] = s / np.maximum(n, 1)
    return out


def main():
    F = dict(np.load(os.path.join(D, 'lya/forests_lya.npz')))
    ra = np.degrees(F['ra']); dec = np.degrees(F['dec'])
    mask = fits.open(os.path.join(D, 'act_y/wide_mask.fits'), memmap=True)[0].data
    jy, jx = radec_to_pix(ra, dec)
    inb = (jy >= 0) & (jy < NY)
    mval = np.zeros(len(ra)); mval[inb] = mask[jy[inb], jx[inb]]
    inact = mval > 0.9
    print('forests', len(ra), 'in ACT wide mask (>0.9):', inact.sum())
    for b in range(3):
        s = inact & (F['W'][:, b] > 0)
        print(f'  bin {b}: forests {s.sum()}, pixels {F["npix"][s, b].sum()}')
    # null placement for the noise
    rng = np.random.default_rng(148)
    rra = np.mod(ra + 97.0, 360)
    jy2, jx2 = radec_to_pix(rra, dec)
    ok2 = (jy2 >= 0) & (jy2 < NY)
    mv2 = np.zeros(len(ra)); mv2[ok2] = mask[jy2[ok2], jx2[ok2]]
    cand = np.where(mv2 > 0.9)[0]
    pick = rng.choice(cand, size=min(6000, len(cand)), replace=False)
    sig = {}
    for name in ('y', 'y_cib'):
        m = fits.open(os.path.join(D, f'act_y/{name}.fits'), memmap=True)[0].data
        am = annulus_means(m, rra[pick], dec[pick])
        good = np.all(np.isfinite(am), 1)
        am = am[good]
        # subtract the 45-60' ring as a local background (as the measurement would), and report both
        sig[name] = am.std(0)
        sig[name + '_bg'] = (am - am[:, -1:]).std(0)
        print(name, 'sigma_ybar(theta) [1e-7]:', np.round(sig[name] * 1e7, 3))
        print(name, ' background-subtracted [1e-7]:', np.round(sig[name + '_bg'] * 1e7, 3))
    th = np.sqrt(EDGES[1:] * EDGES[:-1])
    rows = []
    for name in ('y_cib', 'y'):
        for b in range(3):
            s = inact & (F['W'][:, b] > 0)
            W = F['W'][s, b]; A = F['A'][s, b]; S = F['S'][s, b]; U = F['U'][s, b]
            sw = W.sum()
            noise = np.sqrt((A ** 2).sum()) * sig[name] / sw
            xi_p = S.sum(0) / sw; xi_u = U.sum(0) / sw
            for sc in (0.05, 0.1):
                noise_c = np.sqrt(((sc * W) ** 2 + A ** 2).sum()) * sig[name] / sw
                rows.append((name, b, sc, xi_u / noise_c))
            print(f'{name} bin {b}: noise [1e-10] {np.round(noise * 1e10, 3)}')
            print(f'   model unprojected [1e-10] {np.round(xi_u * 1e10, 3)}  S/N per theta {np.round(xi_u / noise, 2)}')
            print(f'   model projected   [1e-10] {np.round(xi_p * 1e10, 5)}  S/N per theta {np.round(xi_p / noise, 5)}')
            print(f'   S/N quadrature over theta (upper bound): unprojected {np.sqrt(((xi_u / noise) ** 2).sum()):.2f}, '
                  f'projected {np.sqrt(((xi_p / noise) ** 2).sum()):.5f}')
    for name, b, sc, sn in rows:
        print(f'{name} bin {b} known continuum to sigma_c={sc}: S/N per theta {np.round(sn, 3)}, quadrature {np.sqrt((sn ** 2).sum()):.2f}')
    np.savez(os.path.join(D, 'lya/forecast.npz'), theta=th, inact=inact, **{k: v for k, v in sig.items()})


if __name__ == '__main__':
    main()
