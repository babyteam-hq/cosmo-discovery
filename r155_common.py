"""R-155 shared pieces: cosmology (flat LCDM, Om 0.27, h 0.7, the Giant Arc papers' choice), comoving coordinates,
sky matching."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.spatial import cKDTree

DATA = "research/cosmo-discovery/data/r155"
OM, H0, C = 0.27, 70.0, 299792.458

_zg = np.linspace(0, 5, 20001)
_E = np.sqrt(OM * (1 + _zg) ** 3 + 1 - OM)
_Dg = np.concatenate([[0], cumulative_trapezoid(C / H0 / _E, _zg)])


def D(z):
    """Comoving distance, Mpc."""
    return np.interp(z, _zg, _Dg)


def zofD(d):
    return np.interp(d, _Dg, _zg)


def dzdD(z):
    """dz per comoving Mpc along the line of sight."""
    return H0 * np.sqrt(OM * (1 + z) ** 3 + 1 - OM) / C


def unit(ra, dec):
    r, d = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)], -1)


def xyz(ra, dec, dist):
    return unit(ra, dec) * np.asarray(dist)[..., None]


def sky_match(ra1, dec1, ra2, dec2, arcsec, many=False):
    """Indices (i1, i2) of pairs within `arcsec`. many=False: nearest in 2 for each 1 within the radius."""
    u1, u2 = unit(ra1, dec1), unit(ra2, dec2)
    chord = 2 * np.sin(np.radians(arcsec / 3600) / 2)
    t = cKDTree(u2)
    if many:
        lists = t.query_ball_point(u1, chord)
        i1 = np.repeat(np.arange(len(u1)), [len(x) for x in lists])
        i2 = np.concatenate([np.asarray(x, int) for x in lists]) if len(i1) else np.zeros(0, int)
        return i1, i2.astype(int)
    d, i2 = t.query(u1, distance_upper_bound=chord)
    ok = np.isfinite(d)
    return np.nonzero(ok)[0], i2[ok]
