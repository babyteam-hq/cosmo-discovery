"""R-147: Legacy Survey morphology for the BGS lenses, from the DESI DR1 LSS BGS_BRIGHT full catalogue (matched by TARGETID).

Fixed in results-R147.md: DEV -> bulge-dominated (Sersic n = 4), EXP and REX -> disc-dominated (n = 1), SER by its index if the
catalogue carries one (n > 2 bulge), else dropped; PSF dropped. Writes tid (sorted), early, valid.

Usage: python -I r147_morph.py <LSS full fits> <out npz>
"""
import sys
import numpy as np, fitsio


def main():
    src, out = sys.argv[1], sys.argv[2]
    cols = fitsio.FITS(src)[1].get_colnames()
    want = ['TARGETID', 'MORPHTYPE'] + [c for c in ('SERSIC',) if c in cols]
    print('columns used:', want)
    d = fitsio.read(src, columns=want)
    tid, keep = np.unique(d['TARGETID'], return_index=True); d = d[keep]
    mt = np.char.strip(d['MORPHTYPE'].astype(str))
    print('MORPHTYPE counts:', dict(zip(*np.unique(mt, return_counts=True))))
    early = mt == 'DEV'; valid = np.isin(mt, ['DEV', 'EXP', 'REX'])
    if 'SERSIC' in want:
        ser = mt == 'SER'; early |= ser & (d['SERSIC'] > 2); valid |= ser
    np.savez(out, tid=tid, early=early, valid=valid)
    print(f'{valid.sum()} typed of {tid.size}; bulge fraction {early[valid].mean():.3f}')


if __name__ == '__main__':
    main()
