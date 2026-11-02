"""R-148 step 0/1: reduce the DESI DR1 Lya (or CIII control) deltas to per-forest, per-z-bin sums.

  python -u lya_y_prep.py lya     # data/lya/delta-lya/*.fits.gz  -> data/lya/forests_lya.npz
  python -u lya_y_prep.py ciii    # data/lya/delta-ciii/*.fits.gz -> data/lya/forests_ciii.npz

Per forest q and z bin b (pixel z = lambda/1215.67 - 1, also for the CIII control, so its pixels sit at the same
observed wavelengths as the forest bins):
  W[q,b]    = sum_{i in b} w_i
  A[q,b]    = sum_{i in b} w_i (P delta)_i         (the data; P = picca's projection: weighted mean and slope in
                                                     log10 lambda removed over the whole forest, picca's cf default)
  Araw[q,b] = sum_{i in b} w_i delta_i              (unprojected, for the record)
  S[q,b,t]  = sum_{i in b} w_i (P xi_t)_i           (the model template at theta bin t, after projection)
  U[q,b,t]  = sum_{i in b} w_i xi_t,i               (the same without projection: what a known continuum would give)
  A2[q,b]   = sum_{i in b} w_i^2 (P delta)_i^2 is not needed; the error comes from rotations.
The template xi_t(z) is the ACT-beam model of lya_y_model.py (data/lya/model.npz).
"""
import numpy as np, glob, os, sys
from astropy.io import fits
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ZB = np.array([2.0, 2.5, 3.0, 3.6])
LYA = 1215.67
M = dict(np.load(os.path.join(HERE, 'data/lya/model.npz')))


def one(f):
    with fits.open(f) as h:
        lam = h['LAMBDA'].data.astype(float)
        meta = h['METADATA'].data
        d = h['DELTA_BLIND'].data.astype(float)
        w = h['WEIGHT'].data.astype(float)
    zpix = lam / LYA - 1
    binix = np.digitize(zpix, ZB) - 1           # -1 or 3 = outside
    L = np.log10(lam)
    nt = len(M['theta'])
    tmpl = np.stack([np.interp(zpix, M['z'], M['xi_act'][:, t]) for t in range(nt)])      # [nt, nlam]
    tmplp = np.stack([np.interp(zpix, M['z'], M['xi_planck'][:, t]) for t in range(nt)])
    n = len(meta)
    out = dict(ra=meta['RA'], dec=meta['DEC'], zq=meta['Z'], tid=meta['TARGETID'], snr=meta['MEANSNR'],
               W=np.zeros((n, 3)), A=np.zeros((n, 3)), Araw=np.zeros((n, 3)), npix=np.zeros((n, 3), int),
               S=np.zeros((n, 3, nt)), U=np.zeros((n, 3, nt)), Sp=np.zeros((n, 3, nt)), Up=np.zeros((n, 3, nt)))
    for q in range(n):
        ok = w[q] > 0
        if ok.sum() < 10:
            continue
        ww = w[q, ok]; Lq = L[ok]; bq = binix[ok]
        sw = ww.sum(); Lb = (ww * Lq).sum() / sw; dL = Lq - Lb; sl = (ww * dL * dL).sum()

        def proj(v):  # v [..., npix]
            return v - (v * ww).sum(-1, keepdims=True) / sw - dL * ((v * ww * dL).sum(-1, keepdims=True) / sl)
        dq = d[q, ok]
        pd = proj(dq)
        T = tmpl[:, ok]; PT = proj(T)
        Tp = tmplp[:, ok]; PTp = proj(Tp)
        for b in range(3):
            s = bq == b
            if not s.any():
                continue
            out['W'][q, b] = ww[s].sum()
            out['A'][q, b] = (ww[s] * pd[s]).sum()
            out['Araw'][q, b] = (ww[s] * dq[s]).sum()
            out['npix'][q, b] = s.sum()
            out['S'][q, b] = (PT[:, s] * ww[s]).sum(1)
            out['U'][q, b] = (T[:, s] * ww[s]).sum(1)
            out['Sp'][q, b] = (PTp[:, s] * ww[s]).sum(1)
            out['Up'][q, b] = (Tp[:, s] * ww[s]).sum(1)
    return out


def main():
    tracer = sys.argv[1]
    files = sorted(glob.glob(os.path.join(HERE, f'data/lya/delta-{tracer}/*.fits.gz')))
    print(len(files), 'files')
    parts = []
    with ProcessPoolExecutor(4) as ex:
        for k, o in enumerate(ex.map(one, files, chunksize=4)):
            parts.append(o)
            if k % 100 == 0:
                print(k, flush=True)
    cat = {key: np.concatenate([p[key] for p in parts]) for key in parts[0]}
    keep = cat['W'].sum(1) > 0
    cat = {k: v[keep] for k, v in cat.items()}
    np.savez(os.path.join(HERE, f'data/lya/forests_{tracer}.npz'), **cat)
    W = cat['W']
    print('forests', keep.sum(), 'pixels per bin', cat['npix'].sum(0), 'forests per bin', (W > 0).sum(0))
    for b in range(3):
        sw = W[:, b].sum()
        s = cat['S'][:, b].sum(0) / sw; u = cat['U'][:, b].sum(0) / sw
        print(f'bin {b}: survival S/U per theta', np.round(s / u, 4))
        print(f'        projected xi [1e-9]', np.round(s * 1e9, 4), ' unprojected', np.round(u * 1e9, 4))


if __name__ == '__main__':
    main()
