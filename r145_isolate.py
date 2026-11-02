"""R-145: isolated DESI DR1 BGS lenses inside the DES Y3 footprint (cosmo-deep-16).

Reads the DESI DR1 BGS group catalogue (bgs-groups-dr1-v1.0.fits: BGS BRIGHT, r < 19.5, 3-pass footprint, fastspecfit
stellar masses, redshifts either spectroscopic or assigned) and the DES Y3 GLIMPSE mass-map mask (footprint only).
Applies the isolation criterion fixed in results-R145.md and writes the lens table.

Units: h70 throughout (H0 = 70, Om = 0.3, flat). The group catalogue's LOGMSTAR is in h=1 units (log Msun h^-2); it is
converted with +2 log10(1/0.7).

Usage: python -I r145_isolate.py <data dir> [--risol 3.0] [--dv 1500] [--fmass 0.1] [--mlim 17.0] [--out name]
"""
import argparse, os, sys
import numpy as np, fitsio, healpy as hp
from scipy.spatial import cKDTree
from scipy.integrate import cumulative_trapezoid

C = 299792.458
H0, OM = 70.0, 0.3
_zg = np.linspace(0, 1.5, 3001)
_dc = cumulative_trapezoid(C / H0 / np.sqrt(OM * (1 + _zg) ** 3 + 1 - OM), _zg, initial=0)  # Mpc, comoving


def dcom(z): return np.interp(z, _zg, _dc)
def dang(z): return dcom(z) / (1 + z)
def dlum(z): return dcom(z) * (1 + z)


def des_footprint(path):
    return fitsio.read(path, ext=1)['T'].ravel() > 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('data'); ap.add_argument('--risol', type=float, default=3.0)
    ap.add_argument('--dv', type=float, default=1500.); ap.add_argument('--fmass', type=float, default=0.1)
    ap.add_argument('--mlim', type=float, default=17.0); ap.add_argument('--zmin', type=float, default=0.05)
    ap.add_argument('--zmax', type=float, default=0.30); ap.add_argument('--mmax', type=float, default=11.0)
    ap.add_argument('--mmin', type=float, default=9.5)
    ap.add_argument('--out', default='r145_lenses.fits'); ap.add_argument('--nofootprint', action='store_true')
    a = ap.parse_args()
    g = fitsio.read(os.path.join(a.data, 'bgs-groups-dr1-v1.0.fits'))
    logm = g['LOGMSTAR'] + 2 * np.log10(1 / 0.7)
    ra, dec, z = g['RA'], g['DEC'], g['Z']
    spec = (g['Z_ASSIGNED_FLAG'] == 0) | (g['Z_ASSIGNED_FLAG'] == -1)
    ok = np.isfinite(logm)
    foot = des_footprint(os.path.join(a.data, 'glimpse_mask.fits'))
    nside = hp.npix2nside(foot.size)
    inDES = foot[hp.ang2pix(nside, ra, dec, lonlat=True)]
    # candidate lenses
    cand = spec & ok & (z > a.zmin) & (z < a.zmax) & (logm > a.mmin) & (logm < a.mmax) & (g['APP_MAG_R'] < a.mlim)
    if not a.nofootprint:
        cand &= inDES
    # neighbour catalogue: everything with a mass, within a padded region around the candidates
    th = np.radians(90 - dec); ph = np.radians(ra)
    xyz = np.c_[np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]
    tree = cKDTree(xyz[ok])
    idx_ok = np.where(ok)[0]
    ci = np.where(cand)[0]
    theta = a.risol / dang(z[ci])  # rad, projected radius at lens distance
    chord = 2 * np.sin(theta / 2)
    iso = np.ones(ci.size, bool); nnb = np.zeros(ci.size, int); nconserv = np.zeros(ci.size, int)
    lists = tree.query_ball_point(xyz[ci], chord)
    for k, (i, L) in enumerate(zip(ci, lists)):
        j = idx_ok[np.asarray(L, int)]
        j = j[j != i]
        if j.size == 0:
            continue
        # spectroscopic neighbours: velocity window, mass as measured
        dv = C * np.abs(z[j] - z[i]) / (1 + z[i])
        js = spec[j]
        big_s = js & (dv < a.dv) & (logm[j] > logm[i] + np.log10(a.fmass))
        # neighbours without a spectrum (assigned z): placed at the lens redshift (conservative), mass rescaled by
        # the luminosity distance squared from their assigned redshift
        mz = logm[j] + 2 * np.log10(dlum(z[i]) / dlum(np.maximum(z[j], 1e-3)))
        big_p = (~js) & (mz > logm[i] + np.log10(a.fmass))
        nnb[k] = big_s.sum(); nconserv[k] = big_p.sum()
        iso[k] = not (big_s.any() or big_p.any())
    out = ci[iso]
    print(f'candidates {ci.size}, isolated {out.size} ({out.size / max(ci.size, 1):.3f}); rejected by spec neighbours '
          f'{(nnb > 0).sum()}, by photometric-only neighbours {((nnb == 0) & (nconserv > 0)).sum()}')
    t = np.zeros(out.size, dtype=[('TARGETID', 'i8'), ('RA', 'f8'), ('DEC', 'f8'), ('Z', 'f8'), ('LOGM', 'f8'),
                                   ('G_R', 'f8'), ('ABS_MAG_R', 'f8'), ('QUIESCENT', 'i2'), ('APP_MAG_R', 'f8'),
                                   ('IS_SAT', 'i2'), ('M_HALO', 'f8')])
    for c, s in [('TARGETID', 'TARGETID'), ('RA', 'RA'), ('DEC', 'DEC'), ('Z', 'Z'), ('G_R', 'G_R'),
                 ('ABS_MAG_R', 'ABS_MAG_R'), ('QUIESCENT', 'QUIESCENT'), ('APP_MAG_R', 'APP_MAG_R'),
                 ('IS_SAT', 'IS_SAT'), ('M_HALO', 'M_HALO')]:
        t[c] = g[s][out]
    t['LOGM'] = logm[out]
    fitsio.write(os.path.join(a.data, a.out), t, clobber=True)
    print('wrote', a.out, 'median z %.3f, median logM %.2f' % (np.median(t['Z']), np.median(t['LOGM'])))


if __name__ == '__main__':
    main()
