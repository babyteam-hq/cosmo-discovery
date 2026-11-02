"""cosmo-deep-28 (R-160): build the source table from the Astrogeo image listings, RFC rfc_2026b and OCARS.

Per source: position, first/last image date, number of map files and bands, redshift (spectroscopic only, see below),
membership of each frozen QJAR, and whether the source is NEW (every Astrogeo image observed after the claim's database
version of 2020-06-10, so it cannot be in the claim sample).

Redshift rule (fixed before any PA is measured): OCARS column 123-128 with quality flag 's', excluding origins 'b'
(Quaia, Gaia spectrophotometric) and 'c' (Makarov); else OCARS's SDSS column (139-144). Photometric, unreliable and limit
values are dropped.

Usage: python -I jet_sample.py   (reads data/jet/{rfc,listings,ocars.txt}, writes data/jet/sources.csv)
"""
import os, re, sys, csv
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jet_common import DATA, QJARS, CLAIM_CUTOFF, angsep


def rfc_positions():
    pos = {}
    for l in open(os.path.join(DATA, 'rfc', 'rfc_2026b_cat.txt')):
        if not l.startswith('RFC '):
            continue
        name = l[4:14].strip()
        rh, rm, rs = int(l[26:28]), int(l[29:31]), float(l[32:41])
        sg = -1 if l[42] == '-' else 1
        dd, dm, ds = int(l[43:45]), int(l[46:48]), float(l[49:57])
        pos[name] = (15 * (rh + rm / 60 + rs / 3600), sg * (dd + dm / 60 + ds / 3600))
    return pos


def ocars():
    out = {}
    for l in open(os.path.join(DATA, 'ocars.txt'), encoding='latin-1'):
        if l.startswith('#') or len(l) < 130:
            continue
        name = l[9:19].strip()
        rh, rm, rs = l[20:37].split(); dd, dm, ds = l[38:55].split()
        ra = 15 * (int(rh) + int(rm) / 60 + float(rs) / 3600)
        sg = -1 if dd.startswith('-') else 1
        dec = sg * (abs(int(dd)) + int(dm) / 60 + float(ds) / 3600)
        org, zs, fl = l[121], l[122:128].strip(), l[128]
        zsd = l[138:144].strip() if len(l) > 144 else ''
        z, zsrc = np.nan, ''
        if zs and fl == 's' and org not in 'bc':
            z, zsrc = float(zs), 'ocars' + (org.strip() or 'N')
        elif zsd:
            try:
                z, zsrc = float(zsd), 'sdss'
            except ValueError:
                pass
        typ = l[151:153].strip() if len(l) > 153 else ''
        out[name] = (ra, dec, z, zsrc, typ)
    return out


FN = re.compile(r'^(J[0-9A-Z+\-]+)_([SCXUKQLW])_(\d{4})_(\d{2})_(\d{2})_(\w+?)_map\.fits$')


def main():
    pos = rfc_positions(); oc = ocars()
    rows = []
    for f in sorted(os.listdir(os.path.join(DATA, 'listings'))):
        if not f.endswith('.txt'):
            continue
        s = f[:-4]
        maps = [m for m in (FN.match(x.strip()) for x in open(os.path.join(DATA, 'listings', f))) if m]
        if not maps:
            continue
        dates = sorted(f'{m.group(3)}_{m.group(4)}_{m.group(5)}' for m in maps)
        bands = ''.join(sorted(set(m.group(2) for m in maps)))
        if s in pos:
            ra, dec = pos[s]
        elif s in oc:
            ra, dec = oc[s][:2]
        else:
            continue
        z, zsrc, typ = (oc[s][2], oc[s][3], oc[s][4]) if s in oc else (np.nan, '', '')
        new = dates[0] > CLAIM_CUTOFF
        r = dict(name=s, ra=ra, dec=dec, first=dates[0], last=dates[-1], nmap=len(maps), bands=bands,
                 new=int(new), z=z, zsrc=zsrc, type=typ)
        for q in QJARS:
            sep = angsep(np.array([ra]), np.array([dec]), q[1], q[2])[0]
            r[q[0] + '_sep'] = round(sep, 3)
            r[q[0]] = int(sep <= q[3] and np.isfinite(z) and q[4] <= z <= q[5])
        rows.append(r)
    with open(os.path.join(DATA, 'sources.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    n = len(rows); nn = sum(r['new'] for r in rows)
    print(f'sources with map files: {n}; new (first image after {CLAIM_CUTOFF}): {nn}')
    print(f'with spec z: {sum(np.isfinite(r["z"]) for r in rows)}; new with spec z: '
          f'{sum(np.isfinite(r["z"]) and r["new"] for r in rows)}')
    for q in QJARS:
        cap = [r for r in rows if r[q[0] + '_sep'] <= q[3]]
        print(f'{q[0]}: in cap {len(cap)} (new {sum(r["new"] for r in cap)}); in cap+z range {sum(r[q[0]] for r in rows)} '
              f'(new {sum(r[q[0]] and r["new"] for r in rows)}); new in cap with any z '
              f'{sum(r["new"] and np.isfinite(r["z"]) for r in cap)}')


if __name__ == '__main__':
    main()
