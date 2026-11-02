"""R-155 engine: comoving path of DESI DR1 independent sightlines inside a template (union of spheres of radius R around
member absorbers, or the convex hull of member spheres), expected Mg II counts from dN/dz measured on the same
sightlines over the whole sky, observed counts, delta = N_obs / N_exp - 1.

Sightline path: z_abs from max(0.30, 1250 (1 + z_q) / 2796 - 1) (redward of Lya: the VAC does not search the forest;
dN/dz at z ~ 0.8 is zero for z_q > 3.4) to (1 + z_q)(1 - 5000 km/s / c) - 1 (proximate zone removed).

Completeness: dN/dz(z) is measured separately in strata of exposure depth (TSNR2_QSO tertiles), quasar g-band flux
(12 quantile bins, finest at the bright end; the slice's Mg II falls at ~5000 A, in g) and z_q (4 bins), on all independent sightlines. The expected
count on a sightline uses its own stratum. Under Poisson counts with mean (1 + delta) e_i, the maximum-likelihood delta is
sum n_i / sum e_i - 1, so this stratified ratio is the TSNR2 (and flux) weighting; a literal TSNR2-weighted ratio
sum w n / sum w e - 1 with w = TSNR2 is reported as a variant.

EW: DESI EW_2796 is rest-frame and runs 1.03x Z&M's REW on 21,789 absorbers both catalogues hold (r155_ewcal in
results-R155.md), so the matched cuts W > 0.3 and > 0.6 A on the Z&M scale are DESI EW_2796 > 0.309 and > 0.618 A.
"""
import numpy as np
from scipy.spatial import ConvexHull, cKDTree

from r155_common import C, DATA, D, unit, zofD

EW_SCALE = 1.03
Z_LO, Z_HI = 0.55, 1.10          # dN/dz table range (templates at R <= 95 Mpc stay inside 0.62-0.98)
ZB = np.arange(Z_LO, Z_HI + 1e-9, 0.01)
DV_PROX = 5000.0
# g-flux strata: completeness rises ~6000x from the faintest to the brightest quintile, so the bins are finest at the top
FLUX_PCT = [30, 50, 60, 70, 77, 83, 88, 92, 95, 97.5, 99]


class Sky:
    def __init__(self, ewcuts=(0.3, 0.6)):
        s = np.load(f"{DATA}/desi_sightlines.npz")
        keep = ~s["sdss"]
        self.ra, self.dec, self.zq, self.tsnr2, self.fg = (s[k][keep] for k in ("ra", "dec", "zq", "tsnr2", "flux_g"))
        idx = -np.ones(len(keep), int)
        idx[keep] = np.arange(keep.sum())
        self.u = unit(self.ra, self.dec)
        self.tree = cKDTree(self.u)
        self.zmin = np.maximum(0.30, 1250 * (1 + self.zq) / 2796 - 1)
        self.zmax = (1 + self.zq) * (1 - DV_PROX / C) - 1
        self.Dmin, self.Dmax = D(self.zmin), D(np.maximum(self.zmax, 0))
        # strata
        t3 = np.digitize(self.tsnr2, np.percentile(self.tsnr2, [33.3, 66.7]))
        f12 = np.digitize(self.fg, np.percentile(self.fg, FLUX_PCT))
        q4 = np.digitize(self.zq, [1.5, 2.2, 2.8])
        nf = len(FLUX_PCT) + 1
        self.stratum = (t3 * nf + f12) * 4 + q4
        self.nstrata = 3 * nf * 4
        # SDSS footprint: share of all DESI quasars within 1 deg that are SDSS quasars (shift templates need > 0.3; the GA has 0.39-0.53)
        allu = unit(s["ra"], s["dec"])
        self.allu_tree, self.all_sdss = cKDTree(allu), s["sdss"]
        a = np.load(f"{DATA}/desi_absorbers_all.npz")
        ok = idx[a["sl"]] >= 0
        sl, z, w = idx[a["sl"][ok]], a["z"][ok], a["w"][ok]
        inpath = (z >= self.zmin[sl]) & (z <= self.zmax[sl])
        self.ewcuts = ewcuts
        self.abs = {}
        for c in ewcuts:
            m = inpath & (w > c * EW_SCALE)
            o = np.lexsort((z[m], sl[m]))
            self.abs[c] = (sl[m][o], z[m][o])
        self._dndz()

    def _dndz(self):
        """dN/dz per stratum and z bin; cumulative F(z) = int dN/dz for exact integration over path intervals."""
        nb = len(ZB) - 1
        path = np.zeros((self.nstrata, nb))
        for b in range(nb):
            lo = np.maximum(self.zmin, ZB[b])
            hi = np.minimum(self.zmax, ZB[b + 1])
            path[:, b] = np.bincount(self.stratum, weights=np.clip(hi - lo, 0, None), minlength=self.nstrata)
        self.path = path
        self.dndz, self.F = {}, {}
        for c in self.ewcuts:
            sl, z = self.abs[c]
            m = (z >= Z_LO) & (z < Z_HI)
            cnt = np.zeros((self.nstrata, nb))
            np.add.at(cnt, (self.stratum[sl[m]], np.digitize(z[m], ZB) - 1), 1)
            # smooth in z within a stratum (5-bin boxcar) to keep the shot noise of the table small
            k = np.ones(5)
            cs = np.array([np.convolve(r, k, "same") for r in cnt])
            ps = np.array([np.convolve(r, k, "same") for r in path])
            d = np.where(ps > 0, cs / np.maximum(ps, 1e-12), 0)
            self.dndz[c] = d
            self.F[c] = np.concatenate([np.zeros((self.nstrata, 1)), np.cumsum(d * np.diff(ZB), 1)], 1)

    def Fz(self, c, st, z):
        """int_{Z_LO}^{z} dN/dz for each sightline stratum st (arrays)."""
        zc = np.clip(z, ZB[0], ZB[-1])
        b = np.clip(np.digitize(zc, ZB) - 1, 0, len(ZB) - 2)
        F = self.F[c]
        return F[st, b] + self.dndz[c][st, b] * (zc - ZB[b])

    def sdss_share(self, uvec, rad_deg=1.0):
        lst = self.allu_tree.query_ball_point(uvec, 2 * np.sin(np.radians(rad_deg) / 2))
        return self.all_sdss[np.asarray(lst, int)].mean() if len(lst) else 0.0

    # ---- geometry -------------------------------------------------------------------------------------------
    def sphere_intervals(self, centres, R):
        """Union over members of each sightline's comoving-distance interval inside the spheres, clipped to the
        sightline's allowed path. Returns (sightline index, D1, D2) arrays of disjoint intervals."""
        Dm = np.linalg.norm(centres, axis=1)
        cu = centres / Dm[:, None]
        ang = np.arcsin(np.clip(R / Dm, 0, 1))
        si, ci = [], []
        for k in range(len(centres)):
            lst = self.tree.query_ball_point(cu[k], 2 * np.sin(ang[k] / 2))
            si.append(np.asarray(lst, int))
            ci.append(np.full(len(lst), k))
        if not si:
            return np.zeros(0, int), np.zeros(0), np.zeros(0)
        si, ci = np.concatenate(si), np.concatenate(ci)
        uc = np.einsum("ij,ij->i", self.u[si], centres[ci])
        disc = R ** 2 - (Dm[ci] ** 2 - uc ** 2)
        ok = disc > 0
        si, ci, uc, h = si[ok], ci[ok], uc[ok], np.sqrt(disc[ok])
        a1, a2 = np.maximum(uc - h, self.Dmin[si]), np.minimum(uc + h, self.Dmax[si])
        ok = a2 > a1
        return merge_intervals(si[ok], a1[ok], a2[ok])

    def hull_intervals(self, centres, r_sph):
        """Convex hull of member spheres (radius r_sph): each sightline's interval inside it."""
        g = fib_sphere(200) * r_sph
        pts = (centres[:, None, :] + g[None]).reshape(-1, 3)
        hull = ConvexHull(pts)
        nrm, off = hull.equations[:, :3], hull.equations[:, 3]
        cen = centres.mean(0)
        Dc = np.linalg.norm(cen)
        rmax = np.max(np.linalg.norm(pts - cen, axis=1))
        si = np.asarray(self.tree.query_ball_point(cen / Dc, 2 * np.sin(np.arcsin(min(1, rmax / Dc)) / 2)), int)
        u = self.u[si]
        un = u @ nrm.T                      # (n, f): ray x = t u, constraint t (u.n) + off <= 0
        with np.errstate(divide="ignore", invalid="ignore"):
            t = -off[None] / un
        lo = np.where(un < 0, t, -np.inf).max(1)
        hi = np.where(un > 0, t, np.inf).min(1)
        bad = ((un == 0) & (off[None] > 0)).any(1)
        a1, a2 = np.maximum(lo, self.Dmin[si]), np.minimum(hi, self.Dmax[si])
        ok = (a2 > a1) & ~bad
        return si[ok], a1[ok], a2[ok], hull.volume

    def close_sightlines(self, X, kpc):
        """Independent sightlines within `kpc` proper kpc of any member's own sightline, at the member's redshift
        (they can cross the member's own absorbing halo: one-halo, not Gpc-scale, signal)."""
        if kpc <= 0:
            return np.zeros(0, int)
        Dm = np.linalg.norm(X, axis=1)
        zm = zofD(Dm)
        th = kpc / 1000 * (1 + zm) / Dm
        out = [np.asarray(self.tree.query_ball_point(X[k] / Dm[k], 2 * np.sin(th[k] / 2)), int) for k in range(len(X))]
        return np.unique(np.concatenate(out)) if out else np.zeros(0, int)

    # ---- counts ---------------------------------------------------------------------------------------------
    def counts(self, iv, c, weight=None):
        """N_obs, N_exp (and their TSNR2-weighted versions) for intervals iv = (si, D1, D2)."""
        si, d1, d2 = iv
        z1, z2 = zofD(d1), zofD(d2)
        st = self.stratum[si]
        e = self.Fz(c, st, z2) - self.Fz(c, st, z1)
        asl, az = self.abs[c]
        lo = np.searchsorted(asl, si, "left")
        hi = np.searchsorted(asl, si, "right")
        n = np.zeros(len(si))
        for k in np.nonzero(hi > lo)[0]:
            zz = az[lo[k]:hi[k]]
            n[k] = np.count_nonzero((zz >= z1[k]) & (zz <= z2[k]))
        w = self.tsnr2[si]
        return dict(nobs=n.sum(), nexp=e.sum(), nobs_w=(w * n).sum(), nexp_w=(w * e).sum(),
                    nsl=len(np.unique(si)), dz=(z2 - z1).sum())


def drop(iv, bad):
    si = iv[0]
    keep = ~np.isin(si, bad)
    return tuple(x[keep] for x in iv)


def intersect(A, B):
    """Per-sightline intersection of two disjoint interval lists (si, D1, D2)."""
    sa, a1, a2 = A[:3]
    sb, b1, b2 = B[:3]
    ob = np.argsort(sb, kind="stable")
    sb, b1, b2 = sb[ob], b1[ob], b2[ob]
    lo, hi = np.searchsorted(sb, sa, "left"), np.searchsorted(sb, sa, "right")
    rs, r1, r2 = [], [], []
    for k in np.nonzero(hi > lo)[0]:
        x1 = np.maximum(a1[k], b1[lo[k]:hi[k]])
        x2 = np.minimum(a2[k], b2[lo[k]:hi[k]])
        ok = x2 > x1
        rs.extend([sa[k]] * ok.sum()); r1.extend(x1[ok]); r2.extend(x2[ok])
    return np.array(rs, int), np.array(r1), np.array(r2)


def clip_D(iv, d1, d2):
    si, a1, a2 = iv[:3]
    x1, x2 = np.maximum(a1, d1), np.minimum(a2, d2)
    ok = x2 > x1
    return si[ok], x1[ok], x2[ok]


def merge_intervals(si, a1, a2):
    if len(si) == 0:
        return si, a1, a2
    o = np.lexsort((a1, si))
    si, a1, a2 = si[o], a1[o], a2[o]
    out_s, out_1, out_2 = [], [], []
    cs, c1, c2 = si[0], a1[0], a2[0]
    for s, x1, x2 in zip(si[1:], a1[1:], a2[1:]):
        if s == cs and x1 <= c2:
            c2 = max(c2, x2)
        else:
            out_s.append(cs); out_1.append(c1); out_2.append(c2)
            cs, c1, c2 = s, x1, x2
    out_s.append(cs); out_1.append(c1); out_2.append(c2)
    return np.array(out_s), np.array(out_1), np.array(out_2)


def fib_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)


def rotation_to(a, b, spin):
    """Rotation taking unit vector a to unit vector b, then spinning by `spin` rad about b."""
    v = np.cross(a, b)
    s, c = np.linalg.norm(v), np.dot(a, b)
    if s < 1e-12:
        R1 = np.eye(3)
    else:
        k = v / s
        K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        R1 = np.eye(3) + s * K + (1 - c) * K @ K
    K = np.array([[0, -b[2], b[1]], [b[2], 0, -b[0]], [-b[1], b[0], 0]])
    R2 = np.eye(3) + np.sin(spin) * K + (1 - np.cos(spin)) * K @ K
    return R2 @ R1


def mst_mean_edge(X):
    from scipy.sparse.csgraph import minimum_spanning_tree
    from scipy.spatial.distance import squareform, pdist
    T = minimum_spanning_tree(squareform(pdist(X)))
    return T.data.mean()
