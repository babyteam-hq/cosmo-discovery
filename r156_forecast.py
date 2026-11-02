"""R-156 step 0: forecast S/N of the Doppler-magnification dipole for DESI DR1 BGS_BRIGHT, from the real n(z).

usage: python r156_forecast.py [sigma_kappa per bin, comma-separated] [eta]
Also reproduces 2511.17858's normalisation (noise per grid cell, 5,900 deg2) for comparison.
"""
import sys
import numpy as np
import r156_theory as T

D = 'data/r156/'
BINS = [(0.15, 0.25), (0.25, 0.35), (0.35, 0.45)]
dcen = np.arange(25, 150, 10.0)  # 20-150 Mpc/h, the pre-registered range


def load_nz():
    rows = []
    area = 0
    for cap in ['NGC', 'SGC']:
        lines = open(D + 'BGS_BRIGHT_%s_nz.txt' % cap).read().splitlines()
        area += float(lines[1].split('is')[1].split('square')[0])
        a = np.loadtxt(D + 'BGS_BRIGHT_%s_nz.txt' % cap)
        rows.append(a)
    return rows, area


def bin_stats(rows, zlo, zhi):
    N = V = n2V = 0
    for a in rows:
        m = (a[:, 0] > zlo) & (a[:, 0] < zhi)
        N += a[m, 4].sum()
        V += a[m, 5].sum()
        n2V += (a[m, 3] ** 2 * a[m, 5]).sum()
    nbar_pair = n2V / (N)  # pair-weighted mean density seen by a galaxy
    zeff = 0
    w = 0
    for a in rows:
        m = (a[:, 0] > zlo) & (a[:, 0] < zhi)
        zeff += (a[m, 0] * a[m, 3] ** 2 * a[m, 5]).sum()
        w += (a[m, 3] ** 2 * a[m, 5]).sum()
    return N, V, nbar_pair, zeff / w


if __name__ == '__main__':
    sig = [float(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [0.45, 0.45, 0.45]
    eta = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    rows, area = load_nz()
    import fitsio
    ACTUAL = {}
    for zlo, zhi in BINS:
        ACTUAL[(zlo, zhi)] = sum(int(((a > zlo) & (a < zhi)).sum()) for a in
                                 [fitsio.read(D + 'BGS_BRIGHT_%s_clustering.dat.fits' % c, columns=['Z'])['Z']
                                  for c in ['NGC', 'SGC']])
    print('effective area (NGC+SGC, header "area"): %.0f deg2' % area)
    tot = 0
    for (zlo, zhi), s in zip(BINS, sig):
        N, V, nb, ze = bin_stats(rows, zlo, zhi)
        # n(z) and Nbin are completeness-weighted; shot noise and size noise come from the observed galaxies
        nobs = ACTUAL[(zlo, zhi)]
        nb *= nobs / N
        N = nobs
        r, H, f, Dz = T.background(ze)
        b = 1.34 / Dz
        sn, sndiag, t, C = T.snr(ze, b, dcen, N, nb, V, s, eta)
        # 2511.17858 normalisation: noise per cell (lp = 4 Mpc/h), i.e. N_k -> V / lp^3
        sn_c, _, _, _ = T.snr(ze, b, dcen, V / 64.0, nb, V, s, eta)
        tot += sn ** 2
        print('z %.2f-%.2f zeff %.3f N %.0f V %.2e (Mpc/h)^3 nbar_pair %.2e nbar*lp^3 %.2f sigma_k %.2f eta %.2f: '
              'S/N %.2f (diag-only %.2f; per-cell normalisation %.2f)' % (zlo, zhi, ze, N, V, nb, nb * 64, s, eta,
                                                                         sn, sndiag, sn_c))
    print('combined S/N %.2f' % np.sqrt(tot))
