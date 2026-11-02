"""R-150: the claim table from the stacks (dust_z1_stack.py): per bin A_V in the 8 r_p bins, the power-law fits, the ratio to
Menard+2010 at 150 kpc, the combined ratios, Mg II with/without, the colour-by-colour and MW-curve versions, the
achromatic term, and the nulls. usage: python dust_z1_table.py  (prints; logs/dust_table.log)
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dust_z1_common as C
from dust_z1_stack import jk_cov, fit_powerlaw, HARTLAP, NJK, REF_RP

BINS = [('ELG', z) for z in C.ZBINS['ELG']] + [('LRG', z) for z in C.ZBINS['LRG']]
MEN_ERR = 0.1 / 1.4


def load(tr, z):
    fn = os.path.join(C.D, f'stack_{tr}_{z[0]}_{z[1]}.npz')
    return dict(np.load(fn)) if os.path.exists(fn) else None


def chi2_zero(a, jk):
    cov = jk_cov(jk)
    g = np.isfinite(a)
    ci = np.linalg.inv(cov[np.ix_(g, g)]) * HARTLAP
    from scipy.stats import chi2
    x = float(a[g] @ ci @ a[g])
    return x, int(g.sum()), float(chi2.sf(x, g.sum()))


def rho_jk(rp, jk, ref, gamma=0.84):
    """the fixed-slope amplitude ratio for every jackknife sample (for cross-bin covariances)."""
    t = (rp / REF_RP) ** -gamma
    return np.array([np.sum(t * j) / np.sum(t * t) for j in jk]) / ref    # unweighted LS per sample (shape only)


def combine(vals, sig, jks):
    """GLS mean of per-bin ratios with the jackknife cross-covariance of their (unweighted) estimates, rescaled so the
    diagonal matches each bin's own fitted error."""
    J = np.array(jks)                    # (nbin, NJK)
    d = J - J.mean(1, keepdims=True)
    cov = (NJK - 1) / NJK * d @ d.T
    sd = np.sqrt(np.diag(cov))
    corr = cov / np.outer(sd, sd)
    cov = corr * np.outer(sig, sig)
    ci = np.linalg.inv(cov)
    one = np.ones(len(vals))
    v = 1 / (one @ ci @ one)
    return v * (one @ ci @ np.asarray(vals)), np.sqrt(v), corr


def main():
    ref = {law: C.menard_amplitude(law, REF_RP) for law in ('SMC', 'MW')}
    print(f'Menard+2010 reference, rest-frame A_V at {REF_RP:.0f} kpc: SMC {ref["SMC"]:.2e}, MW {ref["MW"]:.2e} (gamma 0.84)')
    print(f'r_p bin edges (kpc): {", ".join(f"{e:.0f}" for e in C.RP_EDGES)}')
    summary = {}
    for est in ('SMC', 'noMgII', 'MW', 'g-r', 'r-z', 'z-W1', 'slope', 'lowEBV', 'highEBV'):
        law = 'MW' if est == 'MW' else 'SMC'
        print(f'\n=== estimator {est} ===')
        vals, sigs, jks = {}, {}, {}
        for tr, z in BINS:
            s = load(tr, z)
            if s is None or f'A_{est}' not in s:
                continue
            a, jk, rp = s[f'A_{est}'], s[f'jk_{est}'], s['rp']
            cov = jk_cov(jk)
            err = np.sqrt(np.diag(cov))
            f1 = fit_powerlaw(rp, a, cov, gamma=0.84)
            rho, srho = f1['A'] / ref[law], f1['sA'] / ref[law]
            srho_tot = np.hypot(srho, rho * MEN_ERR)
            x0, n0, p0 = chi2_zero(a, jk)
            line = f'{tr} {z[0]}-{z[1]}: rho(150 kpc, gamma 0.84) = {rho:.2f} +- {srho_tot:.2f} ' \
                   f'(A_150 = {f1["A"]:.2e} +- {f1["sA"]:.2e}, chi2 {f1["chi2"]:.1f}/{f1["dof"]}); ' \
                   f'detection chi2 vs 0 = {x0:.1f}/{n0} (p {p0:.1e})'
            if est in ('SMC', 'noMgII', 'MW'):
                f2 = fit_powerlaw(rp, a, cov)
                line += f'\n      free slope: A_150 = {f2["A"]:.2e} [{f2["A_lo"]:.2e}, {f2["A_hi"]:.2e}], gamma = {f2["gamma"]:.2f} ' \
                        f'[{f2["g_lo"]:.2f}, {f2["g_hi"]:.2f}], chi2 {f2["chi2"]:.1f}/{f2["dof"]}'
                pl = f1['A'] * (np.array([50., 150., 500., 1500.]) / REF_RP) ** -0.84
                spl = f1['sA'] * (np.array([50., 150., 500., 1500.]) / REF_RP) ** -0.84
                line += '\n      fit A_V at 50/150/500/1500 kpc: ' + ', '.join(f'{p:.1e} +- {e:.1e}' for p, e in zip(pl, spl))
                line += '\n      bins (rp, A_V +- err, x1e-3): ' + ' '.join(f'{r:.0f}:{v * 1e3:.2f}+-{e * 1e3:.2f}' for r, v, e in zip(rp, a, err))
                if est == 'SMC':
                    summary[(tr, z)] = dict(rho=rho, srho=srho_tot, A=f1['A'], sA=f1['sA'], free=f2, rp=rp, a=a, err=err, pl=pl, spl=spl)
            print(line)
            vals[(tr, z)], sigs[(tr, z)], jks[(tr, z)] = rho, srho, rho_jk(rp, jk, ref[law])
        for group in ('ELG', 'LRG', 'all'):
            keys = [k for k in vals if group == 'all' or k[0] == group]
            if len(keys) < 2:
                continue
            m, sm, corr = combine([vals[k] for k in keys], [sigs[k] for k in keys], [jks[k] for k in keys])
            print(f'  combined {group}: rho = {m:.2f} +- {sm:.2f} (stat) +- {m * MEN_ERR:.2f} (Menard amplitude) '
                  f'-> +- {np.hypot(sm, m * MEN_ERR):.2f}; max |corr| between bins {np.max(np.abs(corr - np.eye(len(keys)))):.2f}')
    # magnification: achromatic term from the per-band shifts, per r_p bin
    print('\n=== per-band brightness shifts: dm_b = -m_mu + A_V k_b(z_f) per r_p bin (SMC); fit over the 8 bins with a '
          '-0.84 profile ===')
    for tr, z in BINS:
        s = load(tr, z)
        if s is None:
            continue
        zf = 0.5 * (z[0] + z[1])
        k = C.band_k(zf)
        kb = np.array([k['g'], k['r'], k['z'], k['W1']])
        bands = ('g', 'r', 'z', 'W1')
        rp = s['rp']
        t = (rp / REF_RP) ** -0.84
        # stack bands x bins; model: (-mu0 + A0 k_b) t(rp); jackknife covariance of the 32-vector
        y = np.concatenate([s[f'dm_{b}'] for b in bands])
        J = np.concatenate([s[f'jkdm_{b}'] for b in bands], axis=1)
        cov = jk_cov(J)
        ci = np.linalg.pinv(cov) * (NJK - len(y) - 2) / (NJK - 1)
        X = np.c_[np.concatenate([-t for _ in bands]), np.concatenate([kb_i * t for kb_i in kb])]
        F = X.T @ ci @ X
        p = np.linalg.solve(F, X.T @ ci @ y)
        e = np.sqrt(np.diag(np.linalg.inv(F)))
        print(f'{tr} {z[0]}-{z[1]}: achromatic m_mu(150 kpc) = {p[0]:.2e} +- {e[0]:.2e} mag; A_V(150 kpc) from brightness '
              f'= {p[1]:.2e} +- {e[1]:.2e}')
    # nulls
    fn = os.path.join(C.D, 'stack_nulls.npz')
    if os.path.exists(fn):
        n = np.load(fn)
        print('\n=== nulls (A_V,rest estimator; chi2 against zero, jackknife covariance; fixed-slope amplitude in units of '
              'the Menard SMC reference) ===')
        names = sorted({k.split('__')[0] for k in n.files})
        for nm in names:
            rp, a, jk = n[f'{nm}__0'], n[f'{nm}__1'], n[f'{nm}__2']
            x0, n0, p0 = chi2_zero(a, jk)
            cov = jk_cov(jk)
            f1 = fit_powerlaw(rp, a, cov, gamma=0.84)
            print(f'{nm}: chi2 {x0:.1f}/{n0} (p {p0:.2f}); amplitude/ref = {f1["A"] / ref["SMC"]:.2f} +- {f1["sA"] / ref["SMC"]:.2f}')
    return summary


if __name__ == '__main__' and 'leak' not in sys.argv:
    main()


def leak_fit():
    """D2: per-quartile fixed-shape amplitudes alpha_k = a_dust + s_k * d, s_k the within-quartile colour-magnitude slope,
    d the induced magnification magnitude shift at 150 kpc; GLS with the 4x4 jackknife covariance."""
    ref = C.menard_amplitude('SMC', REF_RP)
    print('\n=== D2: dust vs magnification leak through the colour-luminosity relation (per r-quartile amplitudes) ===')
    rows = {}
    for tr, z in BINS:
        s = load(tr, z)
        if s is None or 'A_q0' not in s:
            continue
        rp = s['rp']
        t = (rp / REF_RP) ** -0.84
        al, jks = [], []
        for k in range(4):
            a, jk = s[f'A_q{k}'], s[f'jk_q{k}']
            g = np.all(np.isfinite(jk), axis=0) & np.isfinite(a)   # a bin empty in some jackknife region is dropped
            cov = jk_cov(jk[:, g])
            ci = np.linalg.inv(cov) * HARTLAP
            wv = ci @ t[g] / (t[g] @ ci @ t[g])     # the fixed-shape GLS amplitude is linear in the data: alpha = wv . a
            al.append(wv @ a[g]); jks.append(jk[:, g] @ wv)
        al = np.array(al); J = np.array(jks)
        d = J - J.mean(1, keepdims=True)
        cov = (NJK - 1) / NJK * d @ d.T
        ci = np.linalg.inv(cov) * (NJK - 4 - 2) / (NJK - 1)
        X = np.c_[np.ones(4), s['cm_slope']]
        F = X.T @ ci @ X
        p = np.linalg.solve(F, X.T @ ci @ al)
        e = np.sqrt(np.diag(np.linalg.inv(F)))
        r = al - X @ p
        rows[(tr, z)] = (p[0] / ref, e[0] / ref, p[1], e[1])
        print(f'{tr} {z[0]}-{z[1]}: quartile amplitudes (x1e-3) {np.round(al * 1e3, 2)} +- {np.round(np.sqrt(np.diag(cov)) * 1e3, 2)}; '
              f'slopes {np.round(s["cm_slope"], 3)}\n      a_dust = {p[0]:.2e} +- {e[0]:.2e} (rho {p[0] / ref:.2f} +- {e[0] / ref:.2f}); '
              f'd (induced shift at 150 kpc) = {p[1]:.3f} +- {e[1]:.3f} mag; chi2 {float(r @ ci @ r):.1f}/2')
    for group in ('ELG', 'LRG'):
        ks = [k for k in rows if k[0] == group]
        if ks:
            v = np.array([rows[k][0] for k in ks]); sg = np.array([rows[k][1] for k in ks])
            w = 1 / sg ** 2
            print(f'  combined {group} (bins independent): rho_dust = {np.sum(w * v) / w.sum():.2f} +- {1 / np.sqrt(w.sum()):.2f}')
    return rows


if __name__ == '__main__' and 'leak' in sys.argv:
    leak_fit()
