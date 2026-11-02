"""R-150 step 0: the feasibility gate for the halo-dust reddening profile at z ~ 0.8-1.6.

Counts (foregrounds, z > 2 quasars, pairs at 20 kpc - 3 Mpc), the per-quasar A_V-estimator noise measured from the data,
the forecast significance for Menard+2010's z ~ 0.36 profile at the foreground redshifts under (a) no evolution and
(b) dust-to-gas scaling with the gas-phase metallicity (dust/gas ~ Z, -0.2 dex from z = 0.36 to 1.2 at fixed mass,
gas columns at fixed proper r_p unchanged), and the DESI targeting loss of reddened quasars (the DR9 random forest,
desitarget's rf_model_dr9_final, rerun on each quasar's photometry reddened by halo dust at z_f).

No colour is paired with a foreground here: the noise comes from the quasars alone and the pair counts from positions.
usage: python -I dust_z1_feas.py  -> prints the tables (logs/dust_feas.log)
"""
import os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dust_z1_common as C

ZMETAL = 10 ** -0.2


def forecast():
    q = C.load_qso(zmin=C.ZQ_MIN)
    print(f'background quasars: z > {C.ZQ_MIN}, main QSO target, r 16.5-23, W1 > 3 sigma: {len(q["z"])} '
          f'(N {q["north"].sum()}, S {(~q["north"]).sum()}); median z {np.median(q["z"]):.2f}, median r {np.median(q["mag"]["r"]):.2f}')
    res, qr = C.colour_residuals(q)
    ok = np.all(np.isfinite(res), axis=1)
    print(f'colour residuals: {ok.sum()} quasars; robust sigma (1.4826 MAD) g-r, r-z, z-W1 = '
          + ', '.join(f'{1.4826 * np.median(np.abs(res[ok, i])):.3f}' for i in range(3))
          + '; std = ' + ', '.join(f'{res[ok, i].std():.3f}' for i in range(3)))
    out = {}
    for tracer in ('ELG', 'LRG'):
        f = np.load(os.path.join(C.D, f'fg_{tracer}.npz'))
        print(f'\n{tracer}: {len(f["z"])} foregrounds, z {f["z"].min():.2f}-{f["z"].max():.2f}')
        for (zlo, zhi) in C.ZBINS[tracer]:
            zf = 0.5 * (zlo + zhi)
            ah, w = C.av_estimator(res, qr, zf)
            good = np.isfinite(ah)
            sig_q = np.zeros(len(ah))
            for k in range(4):
                s = good & (qr == k)
                sig_q[s] = ah[s].std()
            nf = ((f['z'] >= zlo) & (f['z'] < zhi)).sum()
            qi, fi, rp = C.pairs(q['ra'][good], q['dec'][good], q['z'][good], f['ra'], f['dec'], f['z'], zlo, zhi)
            gi = np.flatnonzero(good)
            qi = gi[qi]
            wg = f['w'][fi]
            wq = 1 / sig_q[qi] ** 2
            b = np.digitize(rp, C.RP_EDGES) - 1
            rows = []
            snr2 = {'a': 0., 'b': 0.}
            for j in range(len(C.RP_EDGES) - 1):
                s = b == j
                if not s.any():
                    continue
                # per-quasar summed pair weight in this bin -> variance of the pair-weighted mean (quasars independent)
                W = np.bincount(qi[s], weights=wg[s] * wq[s], minlength=len(ah))
                tot = W.sum()
                var = np.sum(W ** 2 * sig_q ** 2) / tot ** 2
                rpm = np.average(rp[s], weights=wg[s] * wq[s])
                model = np.average(C.menard_av_rest(rp[s]), weights=wg[s] * wq[s])   # proper kpc
                snr2['a'] += model ** 2 / var
                snr2['b'] += (ZMETAL * model) ** 2 / var
                rows.append((C.RP_EDGES[j], C.RP_EDGES[j + 1], s.sum(), len(np.unique(qi[s])), rpm, np.sqrt(var), model))
            print(f'  {tracer} {zlo}-{zhi}: {nf} foregrounds; per-quasar sigma(A_V,rest) = {np.median(sig_q[good]):.3f} mag '
                  f'(quartiles {", ".join(f"{sig_q[good & (qr == k)][0]:.3f}" for k in range(4))}); '
                  f'pairs {len(rp)}, quasars with a pair {len(np.unique(qi))}')
            print('     r_p bin (kpc)   pairs   quasars  <r_p>   sigma(A_V)   model (a) A_V,rest')
            for r in rows:
                print(f'     {r[0]:6.0f}-{r[1]:6.0f} {r[2]:8d} {r[3]:8d} {r[4]:7.0f}   {r[5]:.2e}   {r[6]:.2e}')
            print(f'     forecast S/N (bins in quadrature, independent quasars): (a) no evolution {np.sqrt(snr2["a"]):.1f}, '
                  f'(b) dust/gas ~ Z {np.sqrt(snr2["b"]):.1f}')
            out[(tracer, zlo)] = (np.sqrt(snr2['a']), np.sqrt(snr2['b']))
    tot_a = np.sqrt(sum(v[0] ** 2 for v in out.values()))
    tot_b = np.sqrt(sum(v[1] ** 2 for v in out.values()))
    eg = [v for k, v in out.items() if k[0] == 'ELG']
    print(f'\nforecast, all six bins in quadrature: (a) {tot_a:.1f} sigma, (b) {tot_b:.1f} sigma; '
          f'ELG alone (a) {np.sqrt(sum(v[0] ** 2 for v in eg)):.1f}, (b) {np.sqrt(sum(v[1] ** 2 for v in eg)):.1f}')
    return q, res, qr


def selection_loss(q, res, qr, nsub=40000, seed=1):
    """Fraction of z > 2 quasars the DR9 QSO random forest loses when reddened by halo dust at zf, and the attenuation of
    the measured mean A_V among the survivors (response = measured / true)."""
    sys.path.insert(0, os.path.join(C.D, 'ext', 'desitarget', 'py'))
    from desitarget.myRF import myRF   # numpy-only module; cuts.py pulls desimodel, so its two helpers are copied below

    def shift_photo_north(g, r, z):   # desitarget.cuts.shift_photo_north (DR9 north -> south colour shift), verbatim maths
        gs = (g * 10 ** (-0.4 * 0.004) * (g / r) ** complex(-0.059)).real
        rs = (r * 10 ** (0.4 * 0.003) * (r / z) ** complex(-0.024)).real
        zs = (z * 10 ** (0.4 * 0.013) * (r / z) ** complex(+0.015)).real
        return gs, rs, zs

    def _getColors(n, nf, g, r, z, w1, w2):   # desitarget.cuts._getColors
        lim = 1e-4
        m = [np.where(f.clip(lim) > lim, 22.5 - 2.5 * np.log10(f.clip(lim)), 0.) for f in (g, r, z, w1, w2)]
        G, R, Z, W1, W2 = m
        photok = (G > 0) & (R > 0) & (Z > 0) & (W1 > 0) & (W2 > 0)
        cols = np.c_[G - R, R - Z, G - Z, G - W1, R - W1, Z - W1, G - W2, R - W2, Z - W2, W1 - W2, R]
        return cols, R, photok
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(q['z']), size=min(nsub, len(q['z'])), replace=False)
    rf_file = os.path.join(C.D, 'ext', 'desitarget', 'py', 'desitarget', 'data', 'rf_model_dr9_final.npz')

    def select(g, r, z, w1, w2, north, ra, dec, gn, rn, zn):
        sel = np.zeros(len(g), bool)
        prob = np.zeros(len(g))
        for nflag in (False, True):
            s = np.flatnonzero(north == nflag)
            gg, rr, zz = (shift_photo_north(g[s], r[s], z[s]) if nflag else (g[s], r[s], z[s]))
            cols, rmag, photok = _getColors(len(s), 11, gg, rr, zz, w1[s], w2[s])
            rf = myRF(cols, os.path.dirname(rf_file), numberOfTrees=500, version=2)
            rf.loadForest(rf_file)
            p = rf.predict_proba()
            if nflag:
                pcut = 0.88 - 0.04 * np.tanh(rmag - 20.5)
            else:
                is_des = (gn[s] > 4) & (rn[s] > 4) & (zn[s] > 4) & ((ra[s] >= 320) | (ra[s] <= 100)) & (dec[s] <= 10)
                pcut = np.where(is_des, 0.7 - 0.05 * np.tanh(rmag - 20.5), 0.84 - 0.04 * np.tanh(rmag - 20.5))
            W1 = 22.5 - 2.5 * np.log10(np.clip(w1[s], 1e-4, None))
            W2 = 22.5 - 2.5 * np.log10(np.clip(w2[s], 1e-4, None))
            sel[s] = (p >= pcut) & (W1 <= 22.3) & (W2 <= 22.3) & (rmag < 23) & (rmag > 16.5) & photok
            prob[s] = p
        return sel, prob

    fl = {b: q[f'flux_{b.lower()}'][idx] / q[f'mw_transmission_{b.lower()}'][idx] for b in ('G', 'R', 'Z', 'W1', 'W2')}
    args = (q['north'][idx], q['ra'][idx], q['dec'][idx], q['nobs_g'][idx], q['nobs_r'][idx], q['nobs_z'][idx])
    t0 = time.time()
    sel0, p0 = select(fl['G'], fl['R'], fl['Z'], fl['W1'], fl['W2'], *args)
    print(f'\nselection: {len(idx)} z > 2 quasars rerun through the DR9 RF ({time.time() - t0:.0f} s): '
          f'{sel0.mean():.3f} pass as observed (the rest were targeted at a different photometry version or in SV)')
    base = np.flatnonzero(sel0)
    for zf in (0.95, 1.2, 1.45, 0.7):
        k = C.band_k(zf)
        ah, _w = C.av_estimator(res, qr, zf)
        ah = ah[idx]
        for av in (0.01, 0.03, 0.1):
            fac = {b: 10 ** (-0.4 * av * k[b]) for b in ('g', 'r', 'z', 'W1', 'W2')}
            sel1, p1 = select(fl['G'] * fac['g'], fl['R'] * fac['r'], fl['Z'] * fac['z'], fl['W1'] * fac['W1'],
                              fl['W2'] * fac['W2'], *args)
            lost = 1 - sel1[base].mean()
            surv = base[sel1[base]]
            fin = np.isfinite(ah)
            # measured shift among survivors: (true A_V added to every estimate) - selection removes some reddened ones
            meas = np.nanmean(ah[surv] + av) - np.nanmean(ah[base])
            print(f'  z_f = {zf}: A_V,rest = {av}: lost {lost:.4f} of selected quasars ({lost / av:.3f} per mag); '
                  f'measured/true mean A_V among survivors = {meas / av:.3f}')


if __name__ == '__main__':
    q, res, qr = forecast()
    if '--no-rf' not in sys.argv:
        selection_loss(q, res, qr)
