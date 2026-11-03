"""R-168 (cosmo-cf4-rev): the claim table. Reads the null fits (data/cf4rev/null_<variant>.npz) and the data fit
(data/cf4rev/data.json) and applies the plan's decision rules R1-R3 (plan-cf4-rev.md).
usage: python cf4rev_table.py
"""
import os, glob, json, numpy as np
from scipy.stats import norm
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'cf4rev')
R = json.load(open(os.path.join(D, 'data.json')))
NULL = {os.path.basename(f)[5:-4]: dict(np.load(f)) for f in glob.glob(os.path.join(D, 'null_*.npz'))}
B = NULL['base']; SAMPLES = ['ZTF z<0.06', 'Pantheon+', 'ZTF', 'DES']; R1_SAMPLES = ('ZTF z<0.06', 'Pantheon+')
TWOMPP_VARS = [v for v in ('b10', 'b14', 'lpt2') if v in NULL]
CF4_VARS = [v for v in NULL if v not in ('base', 'b10', 'b14')]

def stat(null, s, f, k, obs):
    v = null.get(f'{s}|{f}|{k}')
    if v is None or len(v) < 10: return None
    m, sd = v.mean(), v.std(ddof=1)
    return dict(m=m, sd=sd, z=(obs - m)/sd, p_lo=norm.cdf((obs - m)/sd), rank=(v <= obs).mean(), n=len(v))

print('## R1: the 2M++ check (observed 2M++ beta against its own 2M++-like baseline null)')
r1_fail = []
for s in SAMPLES:
    for k in ('br', 'bd'):
        obs = R[s]['2mpp'][k]; st = stat(B, s, '2mpp', k, obs)
        p2 = 2*min(st['p_lo'], 1 - st['p_lo'])
        flag = p2 < 0.01 and s in R1_SAMPLES
        if flag: r1_fail.append((s, k))
        others = ', '.join(f"{v}: {NULL[v][f'{s}|2mpp|{k}'].mean():.2f}±{NULL[v][f'{s}|2mpp|{k}'].std(ddof=1):.2f} (z {(obs - NULL[v][f'{s}|2mpp|{k}'].mean())/NULL[v][f'{s}|2mpp|{k}'].std(ddof=1):+.1f})" for v in TWOMPP_VARS)
        print(f"| {s} | 2M++ | {k} | {obs:.2f} ± {R[s]['2mpp']['e' + k[1]]:.2f} | {st['m']:.2f} ± {st['sd']:.2f} (n={st['n']}) | z {st['z']:+.1f}, two-sided p {p2:.4f}, rank {st['rank']:.3f} | {others} | {'FAIL' if flag else 'pass' if s in R1_SAMPLES else '(not in R1)'} |")
print('R1:', 'FAIL ' + str(r1_fail) if r1_fail else 'pass')

print('\n## The table (baseline null; not interpreted if R1 failed): observed beta, calibrated null, ratio')
for s in SAMPLES:
    for f, nf in (('2mpp', '2mpp'), ('cf4', 'cf4'), ('wfreal', 'cf4')):
        for k in ('br', 'bd', 'corr_rest'):
            obs = R[s][f][k]; st = stat(B, s, nf, k, obs)
            if st is None: continue
            err = f" ± {R[s][f]['e' + k[1]]:.2f}" if k != 'corr_rest' else ''
            print(f"| {s} | {f} | {k} | {obs:.2f}{err} | {st['m']:.2f} ± {st['sd']:.2f} | {obs/st['m']:.2f} ± {st['sd']/st['m']:.2f} | z {st['z']:+.1f} (rank {st['rank']:.3f}) |")
print('\n## variants (CF4 rows only where run)')
for s in SAMPLES:
    for k in ('br',):
        obs = R[s]['cf4'][k]; st = stat(B, s, 'cf4', k, obs)
        vz = {v: stat(NULL[v], s, 'cf4', k, obs) for v in CF4_VARS}
        vz = {v: x for v, x in vz.items() if x is not None}
        if not vz: continue
        worst = max(vz, key=lambda v: vz[v]['z']); ratios = [obs/x['m'] for x in vz.values()]
        print(f"| {s} | CF4++ | {k} | variants {', '.join(f'{v} {x[chr(109)]:.2f}±{x[chr(115)+chr(100)]:.2f}' for v, x in vz.items())} | worst: {worst} z {vz[worst]['z']:+.1f} |")
print('\nper-variant CF4-like null of beta_rest (mean ± sd) and the data z:')
for s in ('ZTF z<0.06', 'Pantheon+'):
    obs = R[s]['cf4']['br']
    print(f"  {s}: " + '; '.join(f"{v} {NULL[v][f'{s}|cf4|br'].mean():.2f}±{NULL[v][f'{s}|cf4|br'].std(ddof=1):.2f} (z {(obs - NULL[v][f'{s}|cf4|br'].mean())/NULL[v][f'{s}|cf4|br'].std(ddof=1):+.1f})"
                         for v in ['base'] + sorted(CF4_VARS) if f'{s}|cf4|br' in NULL[v]))
print('\n## R2 (ZTF z<0.06 beta_rest, every variant beyond 2 null-sd below; applied only if R1 passed):', 'R1 FAILED: not applied' if r1_fail else '')
s = 'ZTF z<0.06'; obs = R[s]['cf4']['br']
zs = {v: stat(NULL[v], s, 'cf4', 'br', obs)['z'] for v in ['base'] + CF4_VARS if f'{s}|cf4|br' in NULL[v]}
print('  z per variant:', {v: round(z, 1) for v, z in zs.items()})
print('  R2 claim holds' if all(z < -2 for z in zs.values()) else ('  consistent (baseline within 2 sd)' if abs(zs['base']) < 2 else '  inconclusive'))
print('\n## R3: paired CF4++ - 2M++ difference')
for s in SAMPLES:
    for k, kk in (('br', 'rest'), ('bd', 'dip')):
        d = R[s]['cf4'][k] - R[s]['2mpp'][k]; bt = R[s]['boot']; j4 = R[s]['joint4']
        nd = B.get(f'{s}|diff|{k}')
        lo, hi = np.percentile(nd, [2.5, 97.5])
        bci = (d - 1.96*bt[f'sd_d_{kk}'], d + 1.96*bt[f'sd_d_{kk}'])
        sig = (d < lo or d > hi) and not (bci[0] <= nd.mean() <= bci[1])
        print(f"| {s} | Δβ_{kk} | {d:+.2f} | bootstrap ± {bt[f'sd_d_{kk}']:.2f} | joint-fit Δ {j4[f'd_{kk}']:+.2f} ± {j4[f'e_d_{kk}']:.2f} | null {nd.mean():+.2f} ± {nd.std(ddof=1):.2f} (95% {lo:+.2f} to {hi:+.2f}) | {'significant' if sig else 'not significant'} |")
print('\n## bulk flow')
for s in SAMPLES:
    bt = R[s]['boot']; na = B.get(f'{s}|angle_cf4|deg'); nt = B.get(f'{s}|angle_true|deg')
    print(f"| {s} | SNe B {bt['B']:.0f} ± {bt['sd_B']:.0f} km/s toward (l, b) = ({bt['l']:.0f}, {bt['b']:.0f}) ± {bt['dir_err68']:.0f}° (68%) | CF4++ at the SNe {bt['cf4_B']:.0f} toward ({bt['cf4_l']:.0f}, {bt['cf4_b']:.0f}) | angle {bt['angle_to_cf4']:.0f}° (bootstrap 16-84%: {bt['angle_p16_84'][0]:.0f}-{bt['angle_p16_84'][1]:.0f}) | null angle SN-fit vs CF4-like field: median {np.median(na):.0f}°, 84% {np.percentile(na, 84):.0f}°, 95% {np.percentile(na, 95):.0f}° ; vs truth {np.median(nt):.0f}° | rank {(na <= bt['angle_to_cf4']).mean():.2f} |")
print('\n## data checks'); [print(' ', s, json.dumps(R[s]['checks'])) for s in SAMPLES]
