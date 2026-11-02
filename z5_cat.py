"""GOLDRUSH (Harikane et al. 2022, arXiv:2108.01090; Zenodo 5512721) Wide-layer dropout catalogues -> slim .npy (R-146).

usage: python -I z5_cat.py   (reads data/goldrush/goldrush_ver20210826/masked_{gri,riz}_W_W0?.cat, writes data/z5/cat_{g,r}.npy)
Columns kept: ra, dec, field, the aperture-corrected 2" convolvedflux mags (the paper's total magnitudes), the 1.5" aperture
mags and errors (its detection S/N), and the Galactic extinction a_* (the catalogue mags are not dereddened).
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'data', 'goldrush', 'goldrush_ver20210826')
OUT = os.path.join(HERE, 'data', 'z5')
B = 'grizy'
COLS = {'ra': 1, 'dec': 2, **{f'{b}_ap': 45 + i for i, b in enumerate(B)}, **{f'{b}_aperr': 50 + i for i, b in enumerate(B)},
        **{f'a_{b}': 85 + i for i, b in enumerate(B)}, **{f'{b}_tot': 147 + i for i, b in enumerate(B)}}


def main():
    import pandas as pd
    os.makedirs(OUT, exist_ok=True)
    for drop, tag in (('riz', 'r'), ('gri', 'g')):
        parts = []
        for f in ('W02', 'W03', 'W04', 'W05', 'W06', 'W07'):
            fn = os.path.join(SRC, f'masked_{drop}_W_{f}.cat')
            df = pd.read_csv(fn, sep=r'\s+', comment='#', header=None, usecols=list(COLS.values()), engine='c')
            df.columns = [k for k, v in sorted(COLS.items(), key=lambda kv: kv[1])]
            a = np.zeros(len(df), dtype=[(k, 'f8' if k in ('ra', 'dec') else 'f4') for k in df.columns] + [('field', 'i1')])
            for k in df.columns:
                a[k] = pd.to_numeric(df[k], errors='coerce').to_numpy()
            a['field'] = int(f[1:])
            parts.append(a)
            print(tag, f, len(a), flush=True)
        a = np.concatenate(parts)
        np.save(os.path.join(OUT, f'cat_{tag}.npy'), a)
        print(tag, 'total', len(a), flush=True)


if __name__ == '__main__':
    main()
