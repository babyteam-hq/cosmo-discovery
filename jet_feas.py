"""cosmo-deep-28 (R-160) Step 0: feasibility forecast for the QJAR re-test. Uses positions, dates and redshifts only;
no position angle is measured or read here.

Independent candidates per frozen region (cap radius and z range of 2106.06546 Table 1):
  A: NEW sources (every Astrogeo image observed after 2020-06-10) with a spectroscopic z in range   (primary)
  B: A plus NEW sources with an OCARS photometric / Quaia / Makarov z in range                      (extended)
  C: B plus OLD sources whose OCARS z comes from a post-2020 compilation (origin d DESI, e, f, g)    (most generous)

Claim strength: the S-test D of a region behaves like the mean resultant length R of doubled angles, with
<D> = 0.886/sqrt(n), sigma_D = 0.463/sqrt(n) under random PAs, so the published S.L. Z gives R = (0.886 + 0.463 Z)/sqrt(n).
n (the claim's region count) is not published; it is bracketed by 0.40-0.55 x the OLD sources now in the region with a
NED/SDSS z (the claim's 48 % good-jet rate, 7,290/15,235, with some of today's NED z postdating 2020).
On n' new good PAs the V-test against the published axis has E[V] = sqrt(2 n') R c, c = 0.94 for ~10 deg PA errors.
Combined: Stouffer sum V_k / 2. The good-jet rate on the new, fainter sources is bracketed by 0.30-0.48.

Usage: python -I jet_feas.py
"""
import os, sys, csv
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jet_common import DATA, QJARS


def ocars_any():
    oc = {}
    for l in open(os.path.join(DATA, 'ocars.txt'), encoding='latin-1'):
        if l.startswith('#') or len(l) < 130:
            continue
        zs, fl, org = l[122:128].strip(), l[128], l[121]
        if zs and fl in 'sp':
            try:
                oc[l[9:19].strip()] = (float(zs), org, fl)
            except ValueError:
                pass
    return oc


def main():
    R = list(csv.DictReader(open(os.path.join(DATA, 'sources.csv'))))
    oc = ocars_any()
    rows = []
    for q in QJARS:
        cap = [r for r in R if float(r[q[0] + '_sep']) <= q[3]]
        inz = lambda z: q[4] <= z <= q[5]
        A = [r for r in cap if r['new'] == '1' and r[q[0]] == '1']
        B = [r for r in cap if r['new'] == '1' and (r[q[0]] == '1' or (r['name'] in oc and inz(oc[r['name']][0])))]
        C = B + [r for r in cap if r['new'] == '0' and r['name'] in oc and oc[r['name']][1] in 'defg'
                 and oc[r['name']][2] == 's' and inz(oc[r['name']][0])]
        old = [r for r in cap if r['new'] == '0' and r[q[0]] == '1' and r['zsrc'] in ('ocarsN', 'sdss', 'ocars*')]
        rows.append((q, len(A), len(B), len(C), len(old)))
        print(f'{q[0]}: cap {len(cap)}; candidates A {len(A)}, B {len(B)}, C {len(C)}; old with NED/SDSS z {len(old)}')
    print()
    print('forecast combined V (Stouffer) on the new sample, published strength / half strength:')
    for lab, idx in (('A', 1), ('B', 2), ('C', 3)):
        for g in (0.30, 0.48):
            for nf in (0.40, 0.55):
                tot = {1.0: 0, 0.5: 0}; ff = []
                for q, *cnt in rows:
                    n = nf * cnt[3]; npr = g * cnt[idx - 1]
                    Rk = (0.886 + 0.463 * q[6]) / np.sqrt(n)
                    ff.append(npr / n)
                    for s in tot:
                        tot[s] += np.sqrt(2 * npr) * Rk * s * 0.94
                print(f'  sample {lab}, good-jet rate {g:.2f}, claim n = {nf:.2f} x old: f_k = '
                      + ', '.join(f'{x:.2f}' for x in ff)
                      + f'; Z(full) = {tot[1.0] / 2:.2f}, Z(half) = {tot[0.5] / 2:.2f}')


if __name__ == '__main__':
    main()
