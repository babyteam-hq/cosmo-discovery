"""R-149: tables from r149_mock.py runs (results-R147.md, "R-149 fixed before computing").

Per realisation (line of sight x photo-z seed), each virtual redshift's statistics are combined with the fixed z_v weights;
the 16 sky-plane jackknife regions are deleted together across z_v. Interval: sqrt(mean jackknife variance + variance over
realisations).

Usage: python -I r149_table.py <mock dir> [--riso 3] [--prefix ''] [--desi 0.461,0.314,0.327,0.200,0.310,0.118]
"""
import argparse, glob, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r147_leak import GBAR, G_SI, MSUN, MPC

HI = GBAR > 1e-13
SETS = {'BGS': ({0.10: 0.212, 0.20: 0.580, 0.28: 0.208}, 19.5), 'KiDS': ({0.15: 0.25, 0.25: 0.5, 0.35: 0.25}, 20.0)}
EDGES = np.arange(8.5, 11.0 + 1e-9, 0.05)
RS = (0.3, 0.5, 1.0)


def hist(m, w):
    h = np.histogram(m, EDGES, weights=w)[0]; return h / max(h.sum(), 1e-300)


def weights_to(m, sel, target):
    k = np.clip(np.digitize(m, EDGES) - 1, 0, EDGES.size - 2)
    h = np.histogram(m[sel], EDGES)[0].astype(float); h /= max(h.sum(), 1)
    with np.errstate(divide='ignore', invalid='ignore'):
        u = np.where(h[k] > 0, target[k] / h[k], 0.0)
    return np.where(sel & (m >= 8.5) & (m < 11.0), u, 0.0)


def run_stats(d, riso, keep):
    """All statistics for one file with lenses restricted to `keep` (jackknife). Returns a dict of arrays."""
    q = d['q']; lmp = d['lmp']; row = d['esd_row']; has = row >= 0
    P = d[f'isop{riso:g}'] & keep; T = d[f'isot{riso:g}'] & keep
    # contamination (photo-z-isolated, types mass-matched by overlap as Brouwer / R-147)
    hE = np.histogram(lmp[P & q], EDGES)[0].astype(float); hL = np.histogram(lmp[P & ~q], EDGES)[0].astype(float)
    ov = np.minimum(hE, hL); k = np.clip(np.digitize(lmp, EDGES) - 1, 0, EDGES.size - 2)
    with np.errstate(divide='ignore', invalid='ignore'):
        u = np.where(q, ov[k] / hE[k], ov[k] / hL[k])
    u = np.where(P & np.isfinite(u), u, 0)
    fr = []
    for t in (q, ~q):
        w = u * t; fr.append([np.sum(w * (d['rpmin'] < r)) / w.sum() for r in RS] + [np.sum(w * d['hidsat']) / w.sum()])
    fr = np.array(fr)
    # ESD stacks: a common target mass shape (the minimum over the four samples' normalised M*_ph histograms)
    samp = {('P', 'E'): P & q & has, ('P', 'L'): P & ~q & has, ('T', 'E'): T & q & has, ('T', 'L'): T & ~q & has}
    target = np.min([hist(lmp[s], None) for s in samp.values()], axis=0)
    ds = d['ds']; Rk = d['Rall'][:, :15]
    st = {}
    for key, s in samp.items():
        uu = weights_to(lmp, s, target)[s]; r = row[s]
        w = uu[:, None] * Rk[r] ** 2
        st[key] = (np.sum(w * ds[r, :15], 0) / w.sum(0), np.sum(uu[:, None] * ds[r, 15:], 0) / uu.sum())
    A = {t: st[('P', t)][0] / st[('T', t)][0] for t in 'EL'}
    # Brouwer's App. A benchmark: all types, truly isolated reweighted to the photo-z-isolated M*_ph shape, fixed radii
    sp = P & has; stt = T & has; tg = hist(lmp[sp], None)
    up = weights_to(lmp, sp, tg)[sp]; ut = weights_to(lmp, stt, tg)[stt]
    bp = np.sum(up[:, None] * ds[row[sp], 15:], 0) / up.sum(); bt = np.sum(ut[:, None] * ds[row[stt], 15:], 0) / ut.sum()
    return dict(fr=fr, AE=A['E'], AL=A['L'], bench=bp / bt, ratioT=st[('T', 'E')][0] / st[('T', 'L')][0])


def delta(AE, AL, sel):
    return np.log10(np.mean((AE / AL)[..., sel], axis=-1))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('dir'); ap.add_argument('--riso', type=float, default=3.0)
    ap.add_argument('--prefix', default=''); ap.add_argument('--desi', default='0.461,0.314,0.327,0.200,0.310,0.118')
    a = ap.parse_args()
    desi = dict(zip(['E0.5', 'L0.5', 'E0.3', 'L0.3', 'Esat', 'Lsat'], map(float, a.desi.split(','))))
    Rfix = np.logspace(np.log10(0.03), np.log10(3.0), 11); bsel = Rfix >= 0.299
    Rmed = np.sqrt(G_SI * 10 ** 10.6 * (1 + 10 ** (-0.69 * 10.6 + 6.63)) * MSUN / GBAR) / MPC
    res = {}
    for sname, (zw, mlim) in SETS.items():
        reals = {}
        for f in glob.glob(os.path.join(a.dir, f'r149_{a.prefix}z*_r{mlim:g}_ax*_s*.npz')):
            m = re.search(r'z([\d.]+)_r[\d.]+_ax(\d)_s(\d)\.npz$', f)
            if not m or (a.prefix == '' and os.path.basename(f) != f'r149_z{m.group(1)}_r{mlim:g}_ax{m.group(2)}_s{m.group(3)}.npz'):
                continue
            zv = float(m.group(1))
            if zv not in zw: continue
            reals.setdefault((m.group(2), m.group(3)), {})[zv] = f
        out = []
        for key, fz in sorted(reals.items()):
            if set(fz) != set(zw): continue
            dd = {zv: dict(np.load(f)) for zv, f in fz.items()}
            vals = []
            for jk in range(-1, 16):
                comb = None
                for zv, d in dd.items():
                    keep = d['jk'] != jk
                    s = run_stats(d, a.riso, keep)
                    v = np.r_[s['fr'].ravel(), s['AE'], s['AL'], s['bench'], s['ratioT']]
                    comb = zw[zv] * v if comb is None else comb + zw[zv] * v
                vals.append(comb)
            vals = np.array(vals); n = 16
            out.append((vals[0], (n - 1) / n * np.sum((vals[1:] - vals[1:].mean(0)) ** 2, axis=0), vals[1:]))
        if not out:
            print(f'{sname}: no complete realisation'); continue
        res[sname] = out
        print(f'{sname}: {len(out)} realisations ({", ".join("ax%s s%s" % k for k in sorted(reals) if set(reals[k]) == set(zw))})')
    # unpack: fr (2 x 4) | AE 15 | AL 15 | bench 11 | ratioT 15
    def parts(v):
        return v[..., :8].reshape(v.shape[:-1] + (2, 4)), v[..., 8:23], v[..., 23:38], v[..., 38:49], v[..., 49:64]
    tag = f'R_isol = {a.riso:g} Mpc{", " + a.prefix.rstrip("_") if a.prefix else ""}'
    scale = None
    for sname in ('BGS', 'KiDS'):
        if sname not in res: continue
        full = np.array([o[0] for o in res[sname]]); jkv = np.array([o[1] for o in res[sname]])
        jks = np.array([o[2] for o in res[sname]])                         # (nreal, 16, nv)
        mean = full.mean(0); err = np.sqrt(jkv.mean(0) + (full.var(0, ddof=1) if len(full) > 1 else 0))
        fr, AE, AL, bench, ratioT = parts(mean); efr, eAE, eAL, eb, _ = parts(err)
        print(f'\n### Mock, {sname}-like (r < {SETS[sname][1]:g}), {tag}\n')
        print('Contamination of photo-z-isolated lenses (types mass-matched), mock vs DESI (QUIESCENT, mass-matched):\n')
        print('| statistic | mock Q | mock SF | mock Q − SF | DESI Q | DESI SF |')
        print('|---|---|---|---|---|---|')
        dE = [desi['E0.3'], desi['E0.5'], np.nan, desi['Esat']]; dL = [desi['L0.3'], desi['L0.5'], np.nan, desi['Lsat']]
        for i, lab in enumerate(['0.3 Mpc', '0.5 Mpc', '1 Mpc', 'hidden satellite']):
            dq = full[:, i] - full[:, 4 + i]; ed = np.sqrt(np.mean(jkv[:, i] + jkv[:, 4 + i]) + dq.var(ddof=1))
            print(f'| {lab} | {fr[0, i]:.3f} ± {efr[0, i]:.3f} | {fr[1, i]:.3f} ± {efr[1, i]:.3f} | {dq.mean():.3f} ± {ed:.3f} | '
                  f'{dE[i]:.3f} | {dL[i]:.3f} |'.replace('nan', '—'))
        print(f'\nBrouwer App. A benchmark (all lenses, photo-z / truly isolated ESD, fixed radii): '
              + ', '.join(f'{r:.2f} Mpc {b:.2f}' for r, b in zip(Rfix, bench))
              + f'; mean over 0.3-3 Mpc {bench[bsel].mean():.2f} ± {np.sqrt(np.mean(eb[bsel] ** 2)):.2f} (MICE ~1.30)')
        if sname == 'BGS':
            scale = {'0.5 Mpc': (desi['E0.5'] / fr[0, 1], desi['L0.5'] / fr[1, 1]),
                     '0.3 Mpc': (desi['E0.3'] / fr[0, 0], desi['L0.3'] / fr[1, 0]),
                     'hidden satellite': (desi['Esat'] / fr[0, 3], desi['Lsat'] / fr[1, 3])}
            print('\nRescaling factors f_DESI / f_mock (Q, SF): ' + '; '.join(f'{k} {v[0]:.2f}, {v[1]:.2f}'
                                                                             for k, v in scale.items()))
        print('\nPer g_bar bin (R at M* = 10^10.6): A_Q, A_SF, log10(A_Q/A_SF) ± err, truly isolated log10(ESD_Q/ESD_SF):\n')
        print('| g_bar | R (Mpc) | A_Q | A_SF | log A_Q/A_SF | true-iso log ESD_Q/ESD_SF |')
        print('|---|---|---|---|---|---|')
        lr = np.log10(full[:, 8:23] / full[:, 23:38]); lrj = np.log10(jks[:, :, 8:23] / jks[:, :, 23:38])
        elr = np.sqrt(np.mean((15 / 16) * np.sum((lrj - lrj.mean(1, keepdims=True)) ** 2, 1), 0) + lr.var(0, ddof=1))
        for k in range(15):
            print(f'| {GBAR[k]:.2e} | {Rmed[k]:.3f} | {AE[k]:.3f} | {AL[k]:.3f} | {lr.mean(0)[k]:.3f} ± {elr[k]:.3f} | '
                  f'{np.log10(ratioT[k]):.3f} |')
        print('\nPredicted spurious offset δ = log10 mean_k(A_Q/A_SF) (mean; 68 %; 95 %) and fraction of the published value:\n')
        print('| amplitude | range | δ (dex) | 68 % | 95 % | fraction of 0.19 (u−r) / 0.14 (Sérsic) at g_bar > 1e-13; '
              '0.27 / 0.17 all bins |')
        print('|---|---|---|---|---|---|')
        variants = [('mock (raw)', None)] + ([(f'rescaled to DESI at {k}', v) for k, v in scale.items()] if scale else []) \
            + [('diagnostic: mock excess over a common (SF) reference, R-147-style', 'common')]
        for lab, sc in variants:
            for rlab, sel, pub in (('g_bar > 1e-13', HI, (0.19, 0.14)), ('all 15 bins', slice(None), (0.27, 0.17))):
                def dl(v):
                    ae, al = v[..., 8:23], v[..., 23:38]
                    if isinstance(sc, str): ae = 1 + (ae - 1) * v[..., 49:64]
                    elif sc is not None: ae = 1 + (ae - 1) * sc[0]; al = 1 + (al - 1) * sc[1]
                    return delta(ae, al, sel)
                x = dl(full); xj = dl(jks)
                e = np.sqrt(np.mean((15 / 16) * np.sum((xj - xj.mean(1, keepdims=True)) ** 2, 1)) + (x.var(ddof=1) if x.size > 1 else 0))
                m_ = x.mean()
                print(f'| {lab} | {rlab} | {m_:.3f} | {m_ - e:.3f} to {m_ + e:.3f} | {m_ - 2 * e:.3f} to {m_ + 2 * e:.3f} | '
                      f'{m_ / pub[0]:.2f} ({(m_ - e) / pub[0]:.2f} to {(m_ + e) / pub[0]:.2f}) / {m_ / pub[1]:.2f} '
                      f'({(m_ - e) / pub[1]:.2f} to {(m_ + e) / pub[1]:.2f}) |')
        print()


if __name__ == '__main__':
    main()
