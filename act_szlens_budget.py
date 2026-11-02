"""The systematic budget of R-154 from `act_szlens_fit.py sys` (data/act_szlens/fits/systematics<tag>.json): per bin the shift
of ln R and of zeta for each variant against the central fit, grouped (two-sided pairs as half the up/down difference, single
variants as the full shift), added in quadrature. usage: python -I act_szlens_budget.py [_final]"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
tag = sys.argv[1] if len(sys.argv) > 1 else '_final'
S = json.load(open(os.path.join(HERE, 'data', 'act_szlens', 'fits', f'systematics{tag}.json')))
c = S['central']; lnR0 = np.log(c['R']); z0 = c['trend']['zeta']
groups = {'tSZ profile (UPP vs Planck-13 truth)': ['l_tsz_p13'], 'tSZ mass scaling alpha': ['l_tsz_alpha'],
          'hardening profile (leakage)': ['l_hard_hard15', 'l_hard_hard30'], 'CIB SED and shape': ['l_cib_sed'], 'CIB flux error': ['l_cib_flux'],
          'response (3\' hardening)': ['r_hard30'], 'response error': ['r_plus1sig'],
          'selection: scatter 0.1 / 0.3': ['B_sint0.1', 'B_sint0.3'], 'selection: sigma8': ['B_s8_0.83'], 'selection: optical confirmation': ['B_optM50_1.5'],
          'contamination 0-4 % (6 % at z > 1)': ['f_max'], 'miscentring 0.5\' -> 0': ['mis0'], 'concentration x0.7 / x1.3': ['c0.7', 'c1.3'],
          'two-halo: inner 30\' only': ['inner30'], 'two-halo: amplitude fixed at 1': ['2h_fixed1'], 'mass dependence eta = -0.18': ['eta-0.18']}
tot_R = np.zeros(4); tot_z = 0.0
print(f'central: R {np.round(c["R"], 3)}, zeta {z0:.2f} [{c["trend"]["lo68"]:.2f}, {c["trend"]["hi68"]:.2f}] 95 % [{c["trend"]["lo95"]:.2f}, {c["trend"]["hi95"]:.2f}]')
print('| systematic | dlnR bin 1 | bin 2 | bin 3 | bin 4 | d zeta |')
print('|---|---|---|---|---|---|')
for g, keys in groups.items():
    ks = [k for k in keys if k in S]
    if not ks:
        continue
    dR = [np.log(S[k]['R']) - lnR0 for k in ks]
    dz = [S[k]['trend']['zeta'] - z0 for k in ks if 'trend' in S[k]]
    if len(ks) == 2:
        sR = np.abs(dR[0] - dR[1]) / 2; sz = abs(dz[0] - dz[1]) / 2 if len(dz) == 2 else 0
        sR = np.maximum(sR, np.maximum(np.abs(dR[0]), np.abs(dR[1])) * 0)  # half-difference for brackets
    else:
        sR = np.abs(dR[0]); sz = abs(dz[0]) if dz else 0
    tot_R += sR ** 2; tot_z += sz ** 2
    print(f'| {g} | ' + ' | '.join(f'{x:.3f}' for x in sR) + f' | {sz:.2f} |')
print(f'| **total (quadrature)** | ' + ' | '.join(f'{x:.3f}' for x in np.sqrt(tot_R)) + f' | {np.sqrt(tot_z):.2f} |')
sz = np.sqrt(tot_z)
s68 = (c['trend']['hi68'] - c['trend']['lo68']) / 2
s95lo = z0 - c['trend']['lo95']; s95hi = c['trend']['hi95'] - z0
lo95 = z0 - np.hypot(s95lo, 1.96 * sz); hi95 = z0 + np.hypot(s95hi, 1.96 * sz)
lo68 = z0 - np.hypot(z0 - c['trend']['lo68'], sz); hi68 = z0 + np.hypot(c['trend']['hi68'] - z0, sz)
print(f'zeta with systematics: {z0:.2f} [{lo68:.2f}, {hi68:.2f}] 68 %, [{lo95:.2f}, {hi95:.2f}] 95 %')
sig18 = (1.8 - z0) / np.hypot(c['trend']['hi68'] - z0, sz)
print(f'zeta = 1.8 is {sig18:.1f} sigma above (upper 68 % half-width with systematics); statistical-only Delta chi2 at 1.8: {c["trend"]["dchi2_zeta18"]:.2f}')
for b in range(4):
    print(f'bin {b + 1}: R {c["R"][b]:.3f}, systematic on ln R {np.sqrt(tot_R[b]):.3f}')
