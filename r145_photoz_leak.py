"""R-145 by-product: how isolated are photo-z-isolated lenses? Emulates Brouwer et al. 2021's isolation on DESI DR1 BGS.

Uses only galaxies with spectroscopic redshifts in the BGS group catalogue (lenses and neighbours alike, so both views
see the same neighbour set). The "photo-z view" scatters every redshift by 0.02 (1 + z) (Gaussian, the KiDS-bright ANNz2
scatter), rescales stellar masses by D_L(z_phot)^2 / D_L(z_spec)^2, and applies Brouwer's criterion: no neighbour with
more than 10 % of the lens stellar mass within a 3D distance of 3 Mpc (projected separation and line-of-sight distance
from the redshifts, physical, h70). The "spectroscopic view" asks whether the same lens has such a neighbour within 3 Mpc
projected and 500 km/s. Output: the fraction of photo-z-isolated lenses that are spectroscopically non-isolated, split
by type (QUIESCENT flag) in stellar-mass bins, with binomial errors.

Usage: python -I r145_photoz_leak.py <data dir> [--n 60000] [--seed 1] [--risol 3]
"""
import argparse, os
import numpy as np, fitsio
from scipy.spatial import cKDTree
from scipy.integrate import cumulative_trapezoid

C = 299792.458; H0, OM = 70.0, 0.3
_zg = np.linspace(0, 1.5, 3001)
_dc = cumulative_trapezoid(C / H0 / np.sqrt(OM * (1 + _zg) ** 3 + 1 - OM), _zg, initial=0)
def dcom(z): return np.interp(z, _zg, _dc)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('data'); ap.add_argument('--n', type=int, default=60000)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--risol', type=float, default=3.0)
    a = ap.parse_args(); rng = np.random.default_rng(a.seed)
    g = fitsio.read(os.path.join(a.data, 'bgs-groups-dr1-v1.0.fits'),
                    columns=['RA', 'DEC', 'Z', 'LOGMSTAR', 'Z_ASSIGNED_FLAG', 'QUIESCENT', 'APP_MAG_R'])
    g = g[((g['Z_ASSIGNED_FLAG'] == 0) | (g['Z_ASSIGNED_FLAG'] == -1)) & np.isfinite(g['LOGMSTAR'])]
    z = g['Z']; logm = g['LOGMSTAR'] + 2 * np.log10(1 / 0.7)
    zp = np.clip(z + 0.02 * (1 + z) * rng.standard_normal(z.size), 0.002, None)
    dl = dcom(z) * (1 + z); dlp = dcom(zp) * (1 + zp)
    logmp = logm + 2 * np.log10(dlp / dl)
    th = np.radians(90 - g['DEC']); ph = np.radians(g['RA'])
    xyz = np.c_[np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]
    tree = cKDTree(xyz)
    # lenses: Brouwer-like cuts in the photo-z view (0.1 < z_phot < 0.5, M*_phot < 1e11), spectroscopic z < 0.3 so the
    # BGS neighbour catalogue still reaches 10 % of the lens mass at typical masses
    cand = np.where((zp > 0.1) & (zp < 0.5) & (logmp < 11.0) & (logmp > 10.0) & (z < 0.3))[0]
    cand = rng.choice(cand, size=min(a.n, cand.size), replace=False)
    # search radius: 3 Mpc projected at the smaller of the two redshifts
    zmin = np.minimum(z[cand], zp[cand])
    chord = 2 * np.sin(a.risol / (dcom(zmin) / (1 + zmin)) / 2)
    RS = (0.5, 1.0, a.risol)
    iso_p = np.ones(cand.size, bool); iso_s = np.ones((len(RS), cand.size), bool)
    for k, (i, L) in enumerate(zip(cand, tree.query_ball_point(xyz[cand], chord))):
        j = np.asarray(L, int); j = j[j != i]
        if j.size == 0: continue
        ang = 2 * np.arcsin(np.clip(np.linalg.norm(xyz[j] - xyz[i], axis=1) / 2, 0, 1))
        # photo-z view: 3D distance, physical
        rp = ang * dcom(zp[i]) / (1 + zp[i]); los = np.abs(dcom(zp[j]) - dcom(zp[i])) / (1 + zp[i])
        iso_p[k] = not np.any((np.hypot(rp, los) < a.risol) & (logmp[j] > logmp[i] - 1))
        # spectroscopic view: projected radius and velocity window
        rs = ang * dcom(z[i]) / (1 + z[i]); dv = C * np.abs(z[j] - z[i]) / (1 + z[i])
        big = (dv < 500) & (logm[j] > logm[i] - 1)
        for r, R in enumerate(RS): iso_s[r, k] = not np.any(big & (rs < R))
    q = g['QUIESCENT'][cand].astype(bool); m = logmp[cand]
    print(f'lenses tried {cand.size}; photo-z isolated {iso_p.sum()} ({iso_p.mean():.3f}); spectroscopically isolated at '
          f'{RS} Mpc: {iso_s.mean(axis=1).round(4)}')
    print('mass bin    type   N_photo-iso  frac spec-NON-isolated within ' + ', '.join(f'{R} Mpc' for R in RS))
    for lo, hi in [(10.0, 10.5), (10.5, 11.0), (10.0, 11.0)]:
        for name, sel in [('ETG(Q)', q), ('LTG(SF)', ~q)]:
            s = sel & iso_p & (m >= lo) & (m < hi); n = s.sum()
            fs = [1 - iso_s[r][s].mean() if n else np.nan for r in range(len(RS))]
            print(f'{lo:.1f}-{hi:.1f}  {name:7s} {n:8d}     ' + '   '.join(
                f'{f:.3f} +- {np.sqrt(f * (1 - f) / max(n, 1)):.3f}' for f in fs))


if __name__ == '__main__':
    main()
