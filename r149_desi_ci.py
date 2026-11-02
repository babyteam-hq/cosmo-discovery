"""R-149: the DESI photo-z isolation leakage by type with its interval (the headline that stands from R-147).

Reads R-147's per-lens runs (r147_leak.py, 5 seeds) and the group VAC (for RA, Dec), and gives, for photo-z-isolated lenses
mass-matched between types as Brouwer did, the fraction with a real > 10 % neighbour within R and 500 km/s, per type, their
difference and ratio, and the hidden-satellite fractions. Interval: a sky jackknife (~50 regions of equal lens count, RA
slices within Dec bands) combined with the seed scatter.

Usage: python -I r149_desi_ci.py <data dir> <run dir>
"""
import glob, os, sys
import numpy as np, fitsio
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r147_table import match_weights

RS = (0.3, 0.5, 1.0)


def regions(ra, dec, n=50):
    nd = 5; reg = np.zeros(ra.size, int)
    qd = np.quantile(dec, np.linspace(0, 1, nd + 1)); bd = np.clip(np.digitize(dec, qd[1:-1]), 0, nd - 1)
    for b in range(nd):
        s = bd == b; r = np.mod(ra[s] - 100, 360)          # cut the RA wrap away from the footprint
        qr = np.quantile(r, np.linspace(0, 1, n // nd + 1)); reg[s] = b * (n // nd) + np.clip(np.digitize(r, qr[1:-1]), 0, n // nd - 1)
    return reg


def stats(u, early, rpmin, hsat, sel):
    out = []
    for t in (early, ~early):
        ww = u * t * sel
        out.append([np.sum(ww * (rpmin < r)) / ww.sum() for r in RS] + [np.sum(ww * hsat) / ww.sum()])
    E, L = np.array(out)
    return np.r_[E, L, E - L, E / L]


def main():
    data, run = sys.argv[1:3]
    g = fitsio.read(os.path.join(data, 'bgs-groups-dr1-v1.0.fits'), columns=['TARGETID', 'RA', 'DEC'])
    o = np.argsort(g['TARGETID'])
    files = sorted(glob.glob(os.path.join(run, 'r147_lorentz_s[0-9]*.npz')))
    for split in ('colour', 'quiescent'):
        allv = []; jkvar = []
        for f in files:
            d = np.load(f)
            S = d['iso3']; early = (d['gr'] > 0.75) if split == 'colour' else d['q'].astype(bool)
            # hidden satellite: R-147's component (b) at 3 Mpc, a satellite lens whose own group holds a contaminant
            hsat = d['parts3'][1].sum(axis=1) > 0
            k = o[np.searchsorted(g['TARGETID'][o], d['tid'])]
            reg = regions(g['RA'][k].astype(float), g['DEC'][k].astype(float))
            u = match_weights(d['logmp'], early, S)
            full = stats(u, early, d['rpmin'], hsat, S)
            jk = np.array([stats(u, early, d['rpmin'], hsat, S & (reg != r)) for r in range(reg.max() + 1)])
            n = jk.shape[0]
            jkvar.append((n - 1) / n * np.sum((jk - jk.mean(0)) ** 2, axis=0)); allv.append(full)
        allv = np.array(allv); m = allv.mean(0)
        err = np.sqrt(np.mean(jkvar, axis=0) + allv.var(0, ddof=1))
        lab = [f'{r:g} Mpc' for r in RS] + ['hidden satellite']
        print(f'\n### DESI DR1 BGS, {split} split, 3 Mpc photo-z isolation ({len(files)} seeds, sky jackknife + seed scatter)\n')
        print('| statistic | ETG | LTG | ETG − LTG (68 %; 95 %) | ETG / LTG (68 %; 95 %) |')
        print('|---|---|---|---|---|')
        k4 = len(lab)
        for i, l in enumerate(lab):
            dd, ed = m[2 * k4 + i], err[2 * k4 + i]; rr, er = m[3 * k4 + i], err[3 * k4 + i]
            print(f'| {l} | {m[i]:.3f} ± {err[i]:.3f} | {m[k4 + i]:.3f} ± {err[k4 + i]:.3f} | {dd:.3f} ({dd - ed:.3f}-{dd + ed:.3f}; '
                  f'{dd - 2 * ed:.3f}-{dd + 2 * ed:.3f}) | {rr:.2f} ({rr - er:.2f}-{rr + er:.2f}; {rr - 2 * er:.2f}-{rr + 2 * er:.2f}) |')


if __name__ == '__main__':
    main()
