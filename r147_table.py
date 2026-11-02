"""R-147: tables from r147_leak.py runs: contamination by type, the predicted log g_obs bias per g_bar bin, and Brouwer's
offset statistic log10 <g_obs,E / g_obs,L> with intervals over seeds x bootstrap (results-R147.md, "Fixed before computing").

Usage: python -I r147_table.py <run dir> <model tag, e.g. lorentz> [--morph <morph npz>] [--maglim 17] [--nboot 100]
"""
import argparse, glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r147_leak import GBAR, DSIG_OWN, RISOL, RFRAC, G_SI, MSUN, MPC

MBINS = [(8.5, 10.3), (10.3, 10.6), (10.6, 10.8), (10.8, 11.0)]
HI = GBAR > 1e-13
PUB = {'colour': (0.19, 0.27), 'quiescent': (0.19, 0.27), 'morph': (0.14, 0.17)}   # (g_bar > 1e-13, all bins)


def splits(d, morph):
    out = {'colour': d['gr'] > 0.75, 'quiescent': d['q'].astype(bool)}
    if morph is not None:
        tid, early, okm = morph
        o = np.searchsorted(tid, d['tid']); o = np.clip(o, 0, tid.size - 1); hit = tid[o] == d['tid']
        e = np.zeros(d['tid'].size, bool); valid = np.zeros_like(e)
        e[hit] = early[o[hit]]; valid[hit] = okm[o[hit]]
        out['morph'] = np.where(valid, e, np.nan)
    return out


def match_weights(m, early, sel, lo=8.5, hi=11.0):
    edges = np.arange(lo, hi + 1e-9, 0.05)
    hE = np.histogram(m[sel & early], edges)[0].astype(float); hL = np.histogram(m[sel & ~early], edges)[0].astype(float)
    k = np.clip(np.digitize(m, edges) - 1, 0, edges.size - 2); ov = np.minimum(hE, hL)
    with np.errstate(divide='ignore', invalid='ignore'):
        u = np.where(early, ov[k] / hE[k], ov[k] / hL[k])
    return np.where(sel & (m >= lo) & (m < hi) & np.isfinite(u), u, 0.0)


def bias(u, w, ds):
    """b_k for a weighted lens set: sum u w ds / (DSIG_OWN sum u w); u can be (nboot, n)."""
    num = u @ (w * ds); den = u @ w
    return num / den / DSIG_OWN


def offset(bE, bL, sel):
    return np.log10(np.mean(((1 + bE) / (1 + bL))[..., sel], axis=-1))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('tag')
    ap.add_argument('--morph'); ap.add_argument('--maglim', type=float, default=0)
    ap.add_argument('--nboot', type=int, default=100)
    a = ap.parse_args()
    morph = None
    if a.morph:
        mm = np.load(a.morph); morph = (mm['tid'], mm['early'], mm['valid'])
    files = sorted(glob.glob(os.path.join(a.run, f'r147_{a.tag}_s[0-9]*.npz')))
    print(f'## {a.tag}: {len(files)} seeds{"; lenses r < %g" % a.maglim if a.maglim else ""}\n')
    rng = np.random.default_rng(12345)
    allb = {}; pooled = {}; contam = {}; perbin = {}; comp = {}; nlens = {}; massb = {}
    for f in files:
        d = dict(np.load(f))
        base = np.ones(d['lens'].size, bool)
        if a.maglim: base &= d['mag'] < a.maglim
        for sname, t in splits(d, morph).items():
            valid = np.isfinite(t.astype(float)); early = t.astype(float) == 1
            for R in RISOL:
                S = base & d[f'iso{R:g}'] & valid
                ds = d[f'ds{R:g}']; w = d['w']
                # check (i): spectroscopically isolated lenses carry no contaminating halo
                truly = S & (d['rpmin'] >= R)
                assert np.all(ds[truly] == 0), 'truly isolated lens with a contaminant'
                key = (sname, R)
                if sname == 'colour':
                    Sa = base & d[f'iso{R:g}']; lo = GBAR < 7.56e-14
                    ball = bias(Sa[None].astype(float), w, ds)[0]
                    allb.setdefault(R, []).append((np.mean(ball[lo]), np.mean(ball[~lo]), Sa.sum()))
                # contamination fractions, per type and Brouwer mass bin
                for mb in MBINS + [(8.5, 11.0)]:
                    inb = (d['logmp'] >= mb[0]) & (d['logmp'] < mb[1])
                    for tn, tt in (('E', early), ('L', ~early)):
                        s = S & inb & tt
                        contam.setdefault(key + (mb, tn), []).append(
                            (s.sum(), [np.mean(d['rpmin'][s] < r) if s.any() else np.nan for r in RFRAC]))
                u = match_weights(d['logmp'], early, S)
                nlens.setdefault(key, []).append((np.sum((u > 0) & early), np.sum((u > 0) & ~early)))
                P = rng.poisson(1.0, (a.nboot, u.size)).astype(float)
                bE = bias(np.vstack([u * early, P * u * early]), w, ds)
                bL = bias(np.vstack([u * ~early, P * u * ~early]), w, ds)
                pooled.setdefault(key, []).append(np.c_[offset(bE, bL, HI), offset(bE, bL, slice(None))])
                perbin.setdefault(key, []).append(np.stack([bE[0], bL[0]]))
                parts = d[f'parts{R:g}']
                comp.setdefault(key, []).append(np.array([[bias((u * tt)[None], w, parts[k])[0] for k in range(3)]
                                                          for tt in (early, ~early)]))
                for mb in MBINS:
                    inb = (d['logmp'] >= mb[0]) & (d['logmp'] < mb[1])
                    ub = match_weights(d['logmp'], early, S & inb, *mb)
                    if (ub * early).sum() == 0 or (ub * ~early).sum() == 0: continue
                    be = bias((ub * early)[None], w, ds)[0]; bl = bias((ub * ~early)[None], w, ds)[0]
                    massb.setdefault(key + (mb,), []).append((offset(be, bl, HI), offset(be, bl, slice(None)),
                                                              np.sum((ub > 0) & early), np.sum((ub > 0) & ~early)))
    for R, v in allb.items():
        v = np.array(v).mean(0)
        print(f'All photo-z-isolated lenses, R_isol = {R:g} Mpc (N ~ {v[2]:.0f}, no type split): mean fractional ESD excess '
              f'over the shared own profile {v[0]:.3f} at g_bar < 7.56e-14 (Brouwer: R > 0.3 Mpc/h70 at the mean mass; '
              f'his MICE test ~0.30) and {v[1]:.3f} at larger g_bar.')
    print()
    for key in pooled:
        sname, R = key
        x = np.vstack(pooled[key]); pub = PUB[sname]
        n = np.array(nlens[key]).mean(axis=0)
        print(f'### split {sname}, R_isol = {R:g} Mpc (matched N_E, N_L per seed ~ {n[0]:.0f}, {n[1]:.0f})\n')
        print('Contamination (fraction of photo-z-isolated lenses with a real > 10 % neighbour within R, 500 km/s), '
              'mean over seeds (seed spread in brackets):\n')
        print('| log M*_ph | type | N | ' + ' | '.join(f'{r:g} Mpc' for r in RFRAC) + ' |')
        print('|---|---|---|' + '---|' * len(RFRAC))
        for mb in MBINS + [(8.5, 11.0)]:
            for tn in ('E', 'L'):
                v = contam[key + (mb, tn)]; N = np.mean([q[0] for q in v]); fr = np.array([q[1] for q in v])
                print(f'| {mb[0]}-{mb[1]} | {tn} | {N:.0f} | ' + ' | '.join(
                    f'{m:.3f} ({s:.3f})' for m, s in zip(fr.mean(0), fr.std(0))) + ' |')
        pb = np.array(perbin[key]).mean(0); cp = np.array(comp[key]).mean(0)
        print('\nPredicted log g_obs bias per g_bar bin, mass-matched (mean over seeds); R is the median radius of the bin:\n')
        print('| g_bar (m/s²) | R (Mpc) | ETG log(1+b) | LTG log(1+b) | ETG−LTG | ETG parts a/b/c (b) | LTG parts a/b/c (b) |')
        print('|---|---|---|---|---|---|---|')
        Rmed = np.sqrt(G_SI * 10 ** 10.6 * MSUN / GBAR) / MPC
        for k in range(GBAR.size):
            print(f'| {GBAR[k]:.2e} | {Rmed[k]:.3f} | {np.log10(1 + pb[0, k]):.3f} | {np.log10(1 + pb[1, k]):.3f} | '
                  f'{np.log10((1 + pb[0, k]) / (1 + pb[1, k])):.3f} | ' + '/'.join(f'{c:.3f}' for c in cp[0, :, k]) +
                  ' | ' + '/'.join(f'{c:.3f}' for c in cp[1, :, k]) + ' |')
        print('\nOffset log10<g_E/g_L> predicted by leakage (median; 16-84 %; 2.5-97.5 % over seeds x bootstrap) and as a '
              'fraction of the published value:\n')
        print('| range | predicted δ (dex) | 68 % | 95 % | published | fraction (68 %) |')
        print('|---|---|---|---|---|---|')
        for c, lab in ((0, 'g_bar > 1e-13'), (1, 'all 15 bins')):
            q = np.percentile(x[:, c], [50, 16, 84, 2.5, 97.5])
            print(f'| {lab} | {q[0]:.3f} | {q[1]:.3f} to {q[2]:.3f} | {q[3]:.3f} to {q[4]:.3f} | {pub[c]} | '
                  f'{q[0] / pub[c]:.2f} ({q[1] / pub[c]:.2f}-{q[2] / pub[c]:.2f}) |')
        print('\nPer Brouwer mass bin (types matched within the bin; mean over seeds): δ(g_bar > 1e-13), δ(all), N_E, N_L\n')
        for mb in MBINS:
            v = massb.get(key + (mb,))
            if v: v = np.array(v, float).mean(0); print(f'- {mb[0]}-{mb[1]}: {v[0]:.3f}, {v[1]:.3f}, {v[2]:.0f}, {v[3]:.0f}')
        print()


if __name__ == '__main__':
    main()
