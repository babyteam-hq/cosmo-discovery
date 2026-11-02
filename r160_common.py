"""R-160 shared code: polarisation-angle alignment statistics with parallel transport.

Angles are position angles in radians, measured from north through east (IAU). A pair's alignment is
cos 2(psi_i' - psi_j), with psi_i' source i's angle parallel-transported along the great circle to j's position
(Jain, Narain & Sarala 2004; Tiwari & Jain 2013 eqs 2-5).

Statistics:
- pair_C(pairs, ...): the pair-averaged cos 2 dpsi over a pair list (Tiwari & Jain 2013's S'_D with self excluded,
  pair-weighted); its expectation is the pairwise alignment amplitude C (0 for random angles).
- sd_max(neigh, ...): Hutsemekers / Jain S_D, the mean over sites of the resultant length of the transported
  e^{2i psi} of the site and its neighbours (Tiwari & Jain 2013 eqs 2-3).
"""
import numpy as np
from scipy.spatial import cKDTree

H0, OM = 67.7, 0.31  # flat LCDM (Planck 2018-like); distances in Mpc


def unit(ra_deg, dec_deg):
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], -1)


def basis(x):
    """Local north and east unit vectors at unit positions x (N,3)."""
    z = np.array([0.0, 0.0, 1.0])
    east = np.cross(z, x)
    n = np.linalg.norm(east, axis=-1, keepdims=True)
    east = np.where(n > 1e-12, east / np.maximum(n, 1e-12), np.array([0.0, 1.0, 0.0]))
    north = np.cross(x, east)
    return north, east


def transported_angle(xa, xb, psi_a):
    """Angle of a's polarisation after parallel transport along the great circle a->b, in b's local frame."""
    na, ea = basis(xa)
    v = np.cos(psi_a)[:, None] * na + np.sin(psi_a)[:, None] * ea
    axis = np.cross(xa, xb)
    s = np.linalg.norm(axis, axis=-1)
    c = np.clip(np.sum(xa * xb, -1), -1, 1)
    k = axis / np.maximum(s, 1e-15)[:, None]
    vr = v * c[:, None] + np.cross(k, v) * s[:, None] + k * np.sum(k * v, -1)[:, None] * (1 - c)[:, None]
    small = s < 1e-12
    vr[small] = v[small]
    nb, eb = basis(xb)
    return np.arctan2(np.sum(vr * eb, -1), np.sum(vr * nb, -1))


def pair_cos(x, psi, i, j):
    """cos 2(psi_i' - psi_j) per pair."""
    if len(i) == 0:
        return np.zeros(0)
    t = transported_angle(x[i], x[j], psi[i])
    return np.cos(2 * (t - psi[j]))


def pair_C(x, psi, i, j, w=None):
    c = pair_cos(x, psi, i, j)
    if len(c) == 0:
        return np.nan
    return float(np.mean(c)) if w is None else float(np.sum(w * c) / np.sum(w))


class PairGeom:
    """Precomputed transport rotation per pair so that randomisations only need the angles.

    cos 2(psi_i' - psi_j) = cos 2(psi_i + g_ij - psi_j), with g_ij the frame rotation of the transport.
    """

    def __init__(self, x, i, j):
        self.i, self.j = np.asarray(i), np.asarray(j)
        if len(self.i):
            self.g = transported_angle(x[self.i], x[self.j], np.zeros(len(self.i)))
        else:
            self.g = np.zeros(0)

    def C(self, psi):
        if len(self.i) == 0:
            return np.nan
        return float(np.mean(np.cos(2 * (psi[self.i] + self.g - psi[self.j]))))

    def C_many(self, psis):
        """psis: (R, N) array of angle realisations; returns (R,)"""
        d = psis[:, self.i] + self.g[None, :] - psis[:, self.j]
        return np.mean(np.cos(2 * d), axis=1)


def sd_max_from_pairs(n, i, j, g, psi):
    """S_D (maximised) with neighbourhoods given as directed pairs (neighbour i -> site j) plus the site itself.

    Sites with no neighbour are excluded. Returns mean over sites of |(e^{2i psi_j} + sum_i e^{2i(psi_i+g)})/(1+n_j)|.
    """
    z = np.exp(2j * (psi[i] + g))
    s = np.bincount(j, weights=z.real, minlength=n) + 1j * np.bincount(j, weights=z.imag, minlength=n)
    cnt = np.bincount(j, minlength=n)
    m = cnt > 0
    tot = s[m] + np.exp(2j * psi[m])
    return float(np.mean(np.abs(tot) / (cnt[m] + 1)))


def comoving_distance(z):
    from astropy.cosmology import FlatLambdaCDM

    return FlatLambdaCDM(H0=H0, Om0=OM).comoving_distance(z).value


def pairs_within(points, r, min_sep_points=None, x_unit=None, min_ang_deg=None):
    """Unique pairs (i<j) with Euclidean separation <= r. Optionally drop pairs closer than min_ang_deg on the sky."""
    tree = cKDTree(points)
    p = tree.query_pairs(r, output_type="ndarray")
    if len(p) and min_ang_deg is not None and x_unit is not None:
        cosang = np.sum(x_unit[p[:, 0]] * x_unit[p[:, 1]], -1)
        p = p[cosang < np.cos(np.radians(min_ang_deg))]
    return p


def chord(theta_deg):
    return 2 * np.sin(np.radians(theta_deg) / 2)


def wrap_pa(psi):
    """to [0, pi)"""
    return np.mod(psi, np.pi)


def xmatch(ra1, dec1, ra2, dec2, r_arcsec):
    """Nearest match of each source 1 in catalogue 2 within r; returns (idx2 or -1, sep_arcsec)."""
    x1, x2 = unit(ra1, dec1), unit(ra2, dec2)
    tree = cKDTree(x2)
    d, k = tree.query(x1, k=1, distance_upper_bound=chord(r_arcsec / 3600))
    ok = np.isfinite(d)
    sep = np.where(ok, np.degrees(2 * np.arcsin(np.minimum(d, 2) / 2)) * 3600, np.nan)
    return np.where(ok, k, -1), sep
