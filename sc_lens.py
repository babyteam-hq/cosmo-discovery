"""The mass of eRASS1 superclusters outside their member clusters, from ACT DR6 CMB lensing (R-137, cosmo-deep-13).

Statistic (fixed in results-R137.md before any stack): for each supercluster s and aperture A_s, the band-filtered convergence
integrated over the aperture's clean pixels, K_d = sum_p kappa(p) Omega_p, and the same integral of the member model,
K_m = sum over members of (NFW truncated at r200c + linear two-halo), each member at its posterior-averaged eRASS1 mass, filtered
exactly as the map. The excess E = sum_s w_s (K_d - K_m) / sum_s w_s K_m, with w_s = K_m / var(K_d) (the variance from the first
100 simulations, which contain no clusters). The control: the same E on random groupings of isolated eRASS1 clusters, each member
replaced by an isolated cluster matched in M500 and z, and the aperture carried with it (the supercluster's aperture, shifted so
that the member's position falls on the matched cluster), which superposes N isolated clusters at the members' positions: it
carries every cluster's own one- and two-halo term and no supercluster selection. The result is E - E_control.

Apertures: 'hull' = the members' convex hull (comoving, in the tangent plane at the centre) dilated by HULL_MARGIN comoving Mpc;
'circle' = a disc of CIRCLE_R comoving Mpc around the catalogue centre.

usage (python -I; data under data/, gitignored):
  python -I sc_lens.py prep                  # sample, apertures, controls, nulls, models -> data/erass1_sc/apertures.npz
  python -I sc_lens.py data [--map baseline] # the aperture integrals on a data map (baseline, polonly, cibdeproj, tonly, planck)
  python -I sc_lens.py sims [--worker k --nworkers n]   # the same on every streamed simulation (act_dr6_simstream.sh)
  python -I sc_lens.py fit  [--map baseline]  # the table
"""
import argparse, glob, os, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import erass1_cmblens as r121  # noqa: E402

DATA = os.path.join(HERE, 'data')
ACT = os.path.join(DATA, 'act_dr6')
OUT = os.path.join(DATA, 'erass1_sc')
BASE_OUT = OUT
NSIDE = 512
LMIN, LMAX = 20, 2048
HULL_MARGIN = 5.0      # comoving Mpc
CIRCLE_R = 20.0        # comoving Mpc
# aperture-size variants (information rows, defined after the first look): SC_VAR=<hull margin>_<circle radius>, e.g. 10_30;
# their apertures, models and integrals go to data/erass1_sc/var_<SC_VAR>/
VAR = os.environ.get('SC_VAR', '')
if VAR:
    HULL_MARGIN, CIRCLE_R = (float(x) for x in VAR.split('_'))
    OUT = os.path.join(BASE_OUT, f'var_{VAR}')
ZMIN, ZMAX = 0.05, 0.6
CLEAN_MIN = 0.8        # an aperture is kept if >= 80 % of its pixels are clean
NCTRL = 10             # control realisations
NRAND = 40             # random-position null realisations
ROT_DEG = np.arange(30, 331, 15.0)  # rotated nulls: the whole aperture set turned in RA
MATCH_DZ, MATCH_DLNM = 0.02, 0.15
GH_X, GH_W = np.polynomial.hermite.hermgauss(5)
THETA = np.r_[np.linspace(0, 60, 241), np.geomspace(61, 30 * 60, 400)] / 60 * np.pi / 180  # rad


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


# --------------------------------------------------------------------------------------------------------------- model
def band_filter(ell):
    f = np.zeros(len(ell)); sel = (ell >= LMIN) & (ell <= LMAX); f[sel] = 1.0
    tap = (ell > 0.9 * LMAX) & sel
    f[tap] = 0.5 * (1 + np.cos(np.pi * (ell[tap] - 0.9 * LMAX) / (0.1 * LMAX)))
    low = (ell >= LMIN) & (ell < LMIN + 4)
    f[low] *= 0.5 * (1 - np.cos(np.pi * (ell[low] - LMIN + 1) / 5))
    return f


def profile_grid(cscale=1.0):
    """filtered kappa(theta) for the truncated NFW (per log M500, z) and the unit-bias linear two-halo term (per z):
    kappa_1h(L) = M200 u(L/D_A | M) / (Sigma_cr D_A^2), u the Fourier transform of NFW truncated at r200c (Scoccimarro+2001);
    kappa_2h(L) = rho_m0 b P_lin(L/chi, z) / (Sigma_cr D_A^2) (exact in the flat sky, no truncation)."""
    fn = os.path.join(BASE_OUT, f'profgrid_c{cscale:.2f}.npz')
    if os.path.exists(fn):
        g = np.load(fn); return g['lnM'], g['z'], g['p1'], g['p2'], g['bias']
    from scipy.special import sici, j0
    hm = r121.HaloModel(lmax=LMAX, cscale=cscale)
    ell = np.arange(LMAX + 1).astype(float); F = band_filter(ell)
    lnM = np.linspace(np.log(3e12), np.log(4e15), 26)
    zg = np.r_[np.linspace(0.005, 0.1, 9), np.linspace(0.12, 1.0, 23)]
    J = j0(ell[None, :] * THETA[:, None]) * (ell * F / (2 * np.pi))[None, :]   # (ntheta, nell), trapezoid weight 1 per ell
    p1 = np.zeros((len(lnM), len(zg), len(THETA))); p2 = np.zeros((len(zg), len(THETA))); bias = np.zeros((len(lnM), len(zg)))
    for j, z in enumerate(zg):
        Da = hm.co.angular_diameter_distance(z).value; chi = hm.co.comoving_distance(z).value
        scr = hm.sigma_crit(z)
        kk = np.maximum(ell, 0.5) / chi
        k2 = hm.rho_m0 * hm.plin(z, kk) / (Da ** 2 * scr); k2[0] = 0
        p2[j] = J @ k2
        for i, lm in enumerate(lnM):
            M200, c, rs, rho_s = hm.m500_to_nfw(np.exp(lm), z)
            k = np.maximum(ell, 1e-3) / Da
            x = k * rs
            si1, ci1 = sici((1 + c) * x); si0, ci0 = sici(x)
            mfun = np.log(1 + c) - c / (1 + c)
            u = (np.sin(x) * (si1 - si0) - np.sin(c * x) / ((1 + c) * x) + np.cos(x) * (ci1 - ci0)) / mfun
            k1 = M200 * u / (Da ** 2 * scr)
            p1[i, j] = J @ k1
            bias[i, j] = hm.bias(M200, z)
    np.savez(fn, lnM=lnM, z=zg, p1=p1, p2=p2, bias=bias)
    return lnM, zg, p1, p2, bias


def cluster_profiles(M500, sig, z, grid, bscale=1.0):
    """posterior-averaged filtered 1h and 2h profiles per cluster on THETA (5-node Gauss-Hermite over ln M)."""
    lnMg, zg, p1, p2, bias = grid
    n = len(M500)
    P1 = np.zeros((n, len(THETA))); B = np.zeros(n)
    lnM = np.log(M500 * bscale)
    j = np.clip(np.searchsorted(zg, z) - 1, 0, len(zg) - 2); fz = np.clip((z - zg[j]) / (zg[j + 1] - zg[j]), 0, 1)
    for x, w in zip(GH_X, GH_W / np.sqrt(np.pi)):
        lm = lnM + np.sqrt(2) * sig * x
        i = np.clip(np.searchsorted(lnMg, lm) - 1, 0, len(lnMg) - 2); fm = np.clip((lm - lnMg[i]) / (lnMg[i + 1] - lnMg[i]), 0, 1)
        P1 += w * (((1 - fm) * (1 - fz))[:, None] * p1[i, j] + (fm * (1 - fz))[:, None] * p1[i + 1, j]
                   + ((1 - fm) * fz)[:, None] * p1[i, j + 1] + (fm * fz)[:, None] * p1[i + 1, j + 1])
        B += w * ((1 - fm) * (1 - fz) * bias[i, j] + fm * (1 - fz) * bias[i + 1, j] + (1 - fm) * fz * bias[i, j + 1]
                  + fm * fz * bias[i + 1, j + 1])
    P2 = B[:, None] * ((1 - fz)[:, None] * p2[j] + fz[:, None] * p2[j + 1])
    return P1, P2


# ------------------------------------------------------------------------------------------------------------ geometry
def unitvec(ra, dec):
    r, d = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)], -1)


def tangent(ra0, dec0, v):
    """gnomonic coordinates (east, north) in radians of unit vectors v around (ra0, dec0)."""
    c = unitvec(ra0, dec0)
    e = np.array([-np.sin(np.radians(ra0)), np.cos(np.radians(ra0)), 0.0])
    n = np.cross(c, e)
    cz = v @ c
    return (v @ e) / cz, (v @ n) / cz


def from_tangent(ra0, dec0, x, y):
    c = unitvec(ra0, dec0)
    e = np.array([-np.sin(np.radians(ra0)), np.cos(np.radians(ra0)), 0.0]); n = np.cross(c, e)
    v = c[None, :] + np.asarray(x)[..., None] * e[None, :] + np.asarray(y)[..., None] * n[None, :]
    v /= np.linalg.norm(v, axis=-1, keepdims=True)
    dec = np.degrees(np.arcsin(v[..., 2])); ra = np.degrees(np.arctan2(v[..., 1], v[..., 0])) % 360
    return ra, dec


class Hull:
    """points within `margin` of the convex hull of the member positions P (comoving Mpc, tangent plane)."""

    def __init__(self, P, margin):
        from scipy.spatial import ConvexHull
        self.margin = margin
        P = np.asarray(P, float)
        self.eq = None
        try:
            if len(P) >= 3:
                h = ConvexHull(P); self.eq = h.equations; V = P[h.vertices]
                self.segs = np.stack([V, np.roll(V, -1, 0)], 1)
            else:
                raise ValueError
        except Exception:
            d = np.linalg.norm(P[:, None] - P[None], axis=-1); a, b = np.unravel_index(np.argmax(d), d.shape)
            self.segs = np.array([[P[a], P[b]]])
        self.rmax = np.linalg.norm(P, axis=1).max() + margin

    def dist(self, X):
        A = self.segs[:, 0]; B = self.segs[:, 1]; AB = B - A
        L2 = np.maximum((AB ** 2).sum(1), 1e-12)
        t = np.clip(((X[:, None, :] - A[None]) * AB[None]).sum(-1) / L2[None], 0, 1)
        Q = A[None] + t[..., None] * AB[None]
        d = np.linalg.norm(X[:, None, :] - Q, axis=-1).min(1)
        if self.eq is not None:
            inside = (X @ self.eq[:, :2].T + self.eq[:, 2][None]).max(1) <= 0
            d[inside] = 0
        return d

    def contains(self, X):
        return self.dist(X) <= self.margin


def aperture_pixels(ra0, dec0, chi, shape, offset=(0.0, 0.0), rmax=None):
    """pixels (nside NSIDE) of an aperture defined in the comoving tangent plane of a supercluster, placed so that the plane's
    point `offset` (Mpc) falls on (ra0, dec0); chi converts angle to comoving Mpc. shape: ('hull', Hull) or ('circle', R)."""
    import healpy as hp
    kind, obj = shape
    off = np.asarray(offset, float)
    R = obj.rmax if kind == 'hull' else obj
    rad = (R + np.linalg.norm(off)) / chi * 1.05 + np.radians(0.15)
    pix = hp.query_disc(NSIDE, unitvec(ra0, dec0), rad, inclusive=False)
    v = np.stack(hp.pix2vec(NSIDE, pix), -1)
    x, y = tangent(ra0, dec0, v)
    X = np.c_[x, y] * chi + off[None, :]
    ok = obj.contains(X) if kind == 'hull' else (np.hypot(X[:, 0], X[:, 1]) <= obj)
    return pix[ok], v[ok]


# ------------------------------------------------------------------------------------------------------------- catalogue
def load():
    from astropy.io import fits
    p = fits.open(os.path.join(DATA, 'erass1', 'erass1cl_primary_v3.2.fits'))[1].data
    sc = fits.open(os.path.join(DATA, 'erass1', 'erass1sc.fits'))[1].data
    g = np.load(os.path.join(BASE_OUT, 'members_refined.npz'))
    ra = np.asarray(p['RA_XFIT'], float); dec = np.asarray(p['DEC_XFIT'], float)
    bad = ~np.isfinite(ra) | ~np.isfinite(dec); ra[bad] = np.asarray(p['RA'], float)[bad]; dec[bad] = np.asarray(p['DEC'], float)[bad]
    M = np.asarray(p['M500'], float) * 1e13
    sig = 0.5 * (np.log(np.asarray(p['M500_H'], float)) - np.log(np.asarray(p['M500_L'], float)))
    sig = np.where(np.isfinite(sig) & (sig > 0), sig, np.nanmedian(sig))
    cl = dict(ra=ra, dec=dec, z=np.asarray(p['BEST_Z'], float), M=M, sig=sig, pcont=np.asarray(p['PCONT'], float),
              okM=np.isfinite(M) & (M > 0))
    return p, sc, g, cl


def binary_mask(variant='baseline'):
    import healpy as hp
    fn = os.path.join(BASE_OUT, f'mask_{variant}_{NSIDE}.npy')
    if os.path.exists(fn):
        return np.load(fn)
    if variant == 'planck':
        m = planck_mask_eq()
        b = hp.ud_grade(m, NSIDE) >= 0.99
    else:
        m = hp.read_map(os.path.join(ACT, f'mask_act_dr6_lensing_v1_healpix_nside_4096_{variant}.fits'), dtype=np.float32)
        b = hp.ud_grade(m, NSIDE) >= 0.99
    np.save(fn, b)
    return b


# --------------------------------------------------------------------------------------------------------------- prep
def prep(args):
    import healpy as hp
    from scipy import sparse
    os.makedirs(OUT, exist_ok=True)
    p, sc, g, cl = load()
    b = binary_mask('baseline')
    co = r121.cosmo()
    chi_of = lambda z: co.comoving_distance(z).value
    flat, ptr, ok = g['flat'], g['ptr'], g['ok']
    rng = np.random.default_rng(137)
    log(f'catalogue: {len(sc)} systems; members reproduced for {ok.sum()}; {sc["multiplicity"].sum()} catalogued members')

    # ---- the sample: reproduced, 0.05 < z < 0.6, both apertures >= 80 % clean
    rows = []; npix_tot = 0
    count = dict(zrange=0, zrange_ok=0)
    for k in range(len(sc)):
        zk = sc['redshift'][k]
        if not (ZMIN < zk < ZMAX):
            continue
        count['zrange'] += 1
        if not ok[k]:
            continue
        count['zrange_ok'] += 1
        mem = flat[ptr[k]:ptr[k + 1]]
        if not cl['okM'][mem].all():
            continue
        ra0, dec0 = sc['ra'][k], sc['dec'][k]; chi = chi_of(zk)
        x, y = tangent(ra0, dec0, unitvec(cl['ra'][mem], cl['dec'][mem]))
        P = np.c_[x, y] * chi
        H = Hull(P, HULL_MARGIN)
        ap = {}
        good = True
        for kind, obj in (('hull', H), ('circle', CIRCLE_R)):
            pix, _ = aperture_pixels(ra0, dec0, chi, (kind, obj))
            f = b[pix].mean() if len(pix) else 0
            if f < CLEAN_MIN:
                good = False; break
            ap[kind] = pix
        if not good:
            continue
        rows.append(dict(k=k, z=zk, chi=chi, ra=ra0, dec=dec0, mem=mem, P=P, H=H, ap=ap, N=len(mem)))
    log(f'in 0.05 < z < 0.6: {count["zrange"]} catalogued systems, {count["zrange_ok"]} with members reproduced; '
        f'{len(rows)} with both apertures >= 80 % clean in the ACT DR6 mask: '
        f'{sum(r["N"] == 2 for r in rows)} pairs, {sum(r["N"] >= 3 for r in rows)} with >= 3 members')
    # how many catalogued systems in z range lie in the footprint at all (reproduced or not): centre clean
    cen_clean = np.array([b[hp.ang2pix(NSIDE, sc['ra'][k], sc['dec'][k], lonlat=True)] for k in range(len(sc))])
    zin = (sc['redshift'] > ZMIN) & (sc['redshift'] < ZMAX)
    log(f'catalogued systems in z range with a clean centre: {(zin & cen_clean).sum()} ({(zin & cen_clean & (sc["multiplicity"] >= 3)).sum()} '
        f'with >= 3); of them reproduced: {(zin & cen_clean & ok).sum()} ({(zin & cen_clean & ok & (sc["multiplicity"] >= 3)).sum()})')

    # ---- the isolated pool: Liu-like sample members of no rebuilt system, centre clean, at least 1 deg... (no other cut)
    from sc_members import load_sample
    _, sel = load_sample()
    iso = sel.copy(); iso[g['any_member']] = False; iso[g['all_sc_member']] = False; iso &= cl['okM']
    pixc = hp.ang2pix(NSIDE, cl['ra'], cl['dec'], lonlat=True)
    iso &= b[pixc]
    iso_idx = np.where(iso)[0]
    log(f'isolated pool (the supercluster sample selection, no system membership, in the clean footprint): {len(iso_idx)}')

    from sc_members import load_sample as _ls  # noqa
    # ---- every aperture: (set, row, member, realisation) -> pixels, model contributions
    ap_pix, ap_meta, ap_models = [], [], []   # meta: (set_id, row, member_slot, real)
    SETS = ['hull', 'circle', 'ctrl_hull', 'ctrl_circle', 'rand_hull', 'rand_circle', 'rot_hull', 'rot_circle', 'pop_hull', 'pop_circle']
    # the population pool (a check on the control): any cluster of the supercluster sample selection in the clean footprint
    pop_idx = np.where(sel & cl['okM'] & b[pixc])[0]
    # member models are evaluated later for each aperture from a list of (cluster index) contributors
    contrib = []
    for r_i, r in enumerate(rows):
        for kind in ('hull', 'circle'):
            ap_pix.append(r['ap'][kind]); ap_meta.append((SETS.index(kind), r_i, -1, -1)); contrib.append(r['mem'])
    # controls: isolated pool (the control) and population pool (its check)
    nfail = 0
    for pool, pre in ((iso_idx, 'ctrl_'), (pop_idx, 'pop_')):
     for real in range(NCTRL):
        for r_i, r in enumerate(rows):
            for slot, i in enumerate(r['mem']):
                pl = pool if pre == 'ctrl_' else pool[~np.isin(pool, r['mem'])]
                dz = np.abs(cl['z'][pl] - cl['z'][i]); dm = np.abs(np.log(cl['M'][pl] / cl['M'][i]))
                cand = pl[(dz < MATCH_DZ) & (dm < MATCH_DLNM)]
                if len(cand) < 5:
                    cand = pl[np.argsort(dz / MATCH_DZ + dm / MATCH_DLNM)[:20]]
                for attempt in range(20):
                    j = rng.choice(cand)
                    chij = chi_of(cl['z'][j])
                    got = {}
                    for kind, obj in (('hull', r['H']), ('circle', CIRCLE_R)):
                        pix, _ = aperture_pixels(cl['ra'][j], cl['dec'][j], chij, (kind, obj), offset=r['P'][slot])
                        if len(pix) == 0 or b[pix].mean() < CLEAN_MIN:
                            break
                        got[kind] = pix
                    if len(got) == 2:
                        break
                else:
                    nfail += 1; continue
                for kind in ('hull', 'circle'):
                    ap_pix.append(got[kind]); ap_meta.append((SETS.index(pre + kind), r_i, slot, real)); contrib.append(np.array([j]))
    log(f'controls: {NCTRL} realisations, {nfail} member draws failed after 20 tries')
    # random-position nulls: the supercluster's aperture at a random clean position (orientation kept)
    clean_pix = np.where(b)[0]
    for real in range(NRAND):
        for r_i, r in enumerate(rows):
            for attempt in range(50):
                q = rng.choice(clean_pix); ra1, dec1 = hp.pix2ang(NSIDE, q, lonlat=True)
                got = {}
                for kind, obj in (('hull', r['H']), ('circle', CIRCLE_R)):
                    pix, _ = aperture_pixels(ra1, dec1, r['chi'], (kind, obj))
                    if len(pix) == 0 or b[pix].mean() < CLEAN_MIN:
                        break
                    got[kind] = pix
                if len(got) == 2:
                    break
            else:
                continue
            for kind in ('hull', 'circle'):
                ap_pix.append(got[kind]); ap_meta.append((SETS.index('rand_' + kind), r_i, -1, real)); contrib.append(np.array([], int))
    # rotated nulls: every aperture turned together in RA (the members' configuration and the mask re-applied)
    for real, dra in enumerate(ROT_DEG):
        for r_i, r in enumerate(rows):
            got = {}
            for kind, obj in (('hull', r['H']), ('circle', CIRCLE_R)):
                pix, _ = aperture_pixels((r['ra'] + dra) % 360, r['dec'], r['chi'], (kind, obj))
                if len(pix) == 0 or b[pix].mean() < CLEAN_MIN:
                    break
                got[kind] = pix
            if len(got) < 2:
                continue
            for kind in ('hull', 'circle'):
                ap_pix.append(got[kind]); ap_meta.append((SETS.index('rot_' + kind), r_i, -1, real)); contrib.append(np.array([], int))
    log(f'apertures: {len(ap_pix)}, pixels {sum(len(a) for a in ap_pix):,}')

    # ---- models: members' filtered 1h and 2h integrals over each aperture's clean pixels; and every primary cluster
    grid = profile_grid()
    allc = np.where(cl['okM'] & (cl['pcont'] < 0.5) & (cl['z'] > 0.003))[0]
    need = np.unique(np.concatenate([np.concatenate(contrib).astype(int), allc]))
    P1, P2 = cluster_profiles(cl['M'][need], cl['sig'][need], cl['z'][need], grid)
    rowof = {c: i for i, c in enumerate(need)}
    vall = unitvec(cl['ra'][allc], cl['dec'][allc])
    from scipy.spatial import cKDTree
    tree = cKDTree(vall)
    omega = 4 * np.pi / hp.nside2npix(NSIDE)
    K1 = np.zeros(len(ap_pix)); K2 = np.zeros(len(ap_pix)); KA1 = np.zeros(len(ap_pix)); KA2 = np.zeros(len(ap_pix))
    nclean = np.zeros(len(ap_pix), int)
    rows_W, cols_W = [], []
    for a, (pix, cs) in enumerate(zip(ap_pix, contrib)):
        pc = pix[b[pix]]; nclean[a] = len(pc)
        rows_W.append(np.full(len(pc), a, np.int32)); cols_W.append(pc.astype(np.int32))
        v = np.stack(hp.pix2vec(NSIDE, pc), -1)
        for c in cs:
            th = np.arccos(np.clip(v @ unitvec(cl['ra'][c], cl['dec'][c]), -1, 1))
            K1[a] += np.interp(th, THETA, P1[rowof[c]]).sum() * omega
            K2[a] += np.interp(th, THETA, P2[rowof[c]]).sum() * omega
        meta = ap_meta[a]
        if SETS[meta[0]] in ('hull', 'circle', 'ctrl_hull', 'ctrl_circle', 'pop_hull', 'pop_circle'):
            # every primary cluster (PCONT < 0.5, any z) within the aperture's reach: the 'all detected clusters' model
            cen = v.mean(0); cen /= np.linalg.norm(cen)
            reach = np.arccos(np.clip((v @ cen).min(), -1, 1)) + np.radians(1.0)
            near = tree.query_ball_point(cen, 2 * np.sin(reach / 2))
            for t in near:
                c = allc[t]
                th = np.arccos(np.clip(v @ vall[t], -1, 1))
                KA1[a] += np.interp(th, THETA, P1[rowof[c]]).sum() * omega
                KA2[a] += np.interp(th, THETA, P2[rowof[c]]).sum() * omega
        if a % 2000 == 0:
            log(f'  models {a}/{len(ap_pix)}')
    W = sparse.csr_matrix((np.full(sum(len(c) for c in cols_W), omega, np.float32),
                           (np.concatenate(rows_W), np.concatenate(cols_W))), shape=(len(ap_pix), hp.nside2npix(NSIDE)))
    meta = np.array(ap_meta)
    rinfo = dict(k=np.array([r['k'] for r in rows]), z=np.array([r['z'] for r in rows]), N=np.array([r['N'] for r in rows]),
                 Msum=np.array([cl['M'][r['mem']].sum() for r in rows]),
                 npix_hull=np.array([len(r['ap']['hull']) for r in rows]), npix_circle=np.array([len(r['ap']['circle']) for r in rows]))
    sparse.save_npz(os.path.join(OUT, 'W.npz'), W)
    np.savez(os.path.join(OUT, 'apertures.npz'), meta=meta, K1=K1, K2=K2, KA1=KA1, KA2=KA2, nclean=nclean, sets=np.array(SETS),
             **{f'row_{k}': v for k, v in rinfo.items()}, mem_flat=np.concatenate([r['mem'] for r in rows]),
             mem_ptr=np.r_[0, np.cumsum([r['N'] for r in rows])])
    log('wrote apertures.npz and W.npz')
    for s in ('hull', 'circle'):
        sel = meta[:, 0] == SETS.index(s)
        log(f'{s}: median clean pixels {np.median(nclean[sel]):.0f} ({np.median(nclean[sel]) * omega * (180 / np.pi) ** 2:.1f} deg2); '
            f'model 1h/(1h+2h) {K1[sel].sum() / (K1[sel] + K2[sel]).sum():.3f}; all-cluster model / member model '
            f'{(KA1[sel] + KA2[sel]).sum() / (K1[sel] + K2[sel]).sum():.3f}')


# --------------------------------------------------------------------------------------------------------------- maps
def read_alm(fn):
    import healpy as hp
    a = hp.read_alm(fn).astype(np.complex128); a[~np.isfinite(a)] = 0
    return a


def map_from_alm(a):
    import healpy as hp
    lm = hp.Alm.getlmax(len(a))
    if lm > LMAX:
        a = hp.resize_alm(a, lm, lm, LMAX, LMAX); lm = LMAX
    F = band_filter(np.arange(lm + 1).astype(float))
    return hp.alm2map(hp.almxfl(a, F), NSIDE, lmax=lm)


def planck_mask_eq():
    import healpy as hp
    d = os.path.join(DATA, 'pr3lens', 'COM_Lensing_Szdeproj_4096_R3.00')
    m = hp.read_map(os.path.join(d, 'mask.fits.gz'))
    m = hp.ud_grade(m, 1024)
    return hp.Rotator(coord=['G', 'C']).rotate_map_pixel(m)


def data_map(name):
    import healpy as hp
    if name == 'planck':
        d = os.path.join(DATA, 'pr3lens', 'COM_Lensing_Szdeproj_4096_R3.00')
        a = read_alm(os.path.join(d, 'TT', 'dat_klm.fits'))
        hp.Rotator(coord=['G', 'C']).rotate_alm(a, inplace=True)
        return map_from_alm(a)
    return map_from_alm(read_alm(os.path.join(ACT, f'kappa_alm_data_act_dr6_lensing_v1_{name}.fits')))


def data(args):
    from scipy import sparse
    W = sparse.load_npz(os.path.join(OUT, 'W.npz'))
    m = data_map(args.map)
    K = W @ m
    if args.map == 'planck':
        import healpy as hp
        bp = binary_mask('planck')
        Wp = W.multiply(bp[None, :].astype(np.float32)).tocsr()
        K = Wp @ m
        np.save(os.path.join(OUT, 'planck_cleanfrac.npy'), np.asarray(Wp.sum(1)).ravel() / np.maximum(np.asarray(W.sum(1)).ravel(), 1e-30))
    os.makedirs(OUT, exist_ok=True)
    np.save(os.path.join(OUT, f'K_data_{args.map}.npy'), K)
    log(f'data {args.map}: {len(K)} apertures integrated')


def sims(args):
    from scipy import sparse
    outs = [BASE_OUT] + sorted(glob.glob(os.path.join(BASE_OUT, 'var_*')))
    Ws = [sparse.load_npz(os.path.join(o, 'W.npz')) for o in outs]
    for o in outs:
        os.makedirs(os.path.join(o, 'sims'), exist_ok=True)
    simdir = os.path.join(ACT, 'sims')
    idle = 0
    while True:
        files = sorted(glob.glob(os.path.join(simdir, 'kappa_alm_sim_act_dr6_lensing_v1_baseline_*.fits')))
        files = [f for f in files if int(f[-9:-5]) % args.nworkers == args.worker]
        done_any = False
        for f in files:
            n = f[-9:-5]
            if time.time() - os.path.getmtime(f) < 3:
                continue
            try:
                mp = map_from_alm(read_alm(f))
            except Exception as e:  # a partial file; the fetcher writes .part first, so this is a real read error
                log(f'sim {n}: {e}'); continue
            for o, W in zip(outs, Ws):
                np.save(os.path.join(o, 'sims', f'K_{n}.npy'), (W @ mp).astype(np.float64))
            os.remove(f); open(os.path.join(simdir, f'done_baseline_{n}'), 'w').close()
            log(f'sim {n} done'); done_any = True
        if os.path.exists(os.path.join(ACT, 'simstream_baseline.log')) and 'ALLDONE' in open(os.path.join(ACT, 'simstream_baseline.log')).read() \
                and not files:
            break
        idle = 0 if done_any else idle + 1
        if idle > 720:
            break
        time.sleep(5)


# ---------------------------------------------------------------------------------------------------------------- fit
def stats(K, A, rows_sel, w, model='member', apok=None):
    """excess E per aperture kind on one map's aperture integrals K; returns dict."""
    meta = A['meta']; sets = list(A['sets'])
    # 'all': the members' 1h + 2h plus the one-halo term of every other detected cluster (PCONT < 0.5, any z) in the aperture;
    # the other clusters' two-halo terms are not added (summing every cluster's two-halo term counts the same structure many times)
    Km = (A['K1'] + A['K2']) if model == 'member' else (np.maximum(A['KA1'], A['K1']) + A['K2'])
    out = {}
    if apok is not None:
        K = np.where(apok, K, 0.0); Km = np.where(apok, Km, 0.0)
    for kind in ('hull', 'circle'):
        s = meta[:, 0] == sets.index(kind)
        r = meta[s, 1]; use = rows_sel[r]
        ww = w[kind][r][use]
        out[kind] = (ww * (K[s][use] - Km[s][use])).sum() / (ww * Km[s][use]).sum()
        K1 = np.where(apok, A['K1'], 0.0) if apok is not None else A['K1']
        out['h1_' + kind] = (ww * K1[s][use]).sum() / (ww * Km[s][use]).sum()   # the members' one-halo share of the model
        # control: per realisation, sum the shifted apertures of each grouping, same weights
        sc_ = meta[:, 0] == sets.index('ctrl_' + kind)
        Ec = []; H1c = []
        for real in np.unique(meta[sc_, 3]):
            q = sc_ & (meta[:, 3] == real)
            rr = meta[q, 1]; u = rows_sel[rr]
            ww = w[kind][rr][u]
            Ec.append((ww * (K[q][u] - Km[q][u])).sum() / (ww * Km[q][u]).sum())
            H1c.append((ww * K1[q][u]).sum() / (ww * Km[q][u]).sum())
        out['ctrl_' + kind] = np.mean(Ec); out['ctrl_' + kind + '_real'] = np.array(Ec); out['h1_ctrl_' + kind] = np.mean(H1c)
        sp = meta[:, 0] == sets.index('pop_' + kind)
        Ep = []
        for real in np.unique(meta[sp, 3]):
            q = sp & (meta[:, 3] == real)
            rr = meta[q, 1]; u = rows_sel[rr]; ww = w[kind][rr][u]
            Ep.append((ww * (K[q][u] - Km[q][u])).sum() / (ww * Km[q][u]).sum())
        out['pop_' + kind] = np.mean(Ep)
        # nulls: per realisation, normalised by the supercluster's own member model
        for nul in ('rand', 'rot'):
            sn = meta[:, 0] == sets.index(f'{nul}_{kind}')
            vals = []
            Kms = Km[s]; rows_s = meta[s, 1]
            kmap = dict(zip(rows_s, Kms))
            for real in np.unique(meta[sn, 3]):
                q = sn & (meta[:, 3] == real)
                rr = meta[q, 1]; u = rows_sel[rr]
                ww = w[kind][rr][u]; km = np.array([kmap[x] for x in rr[u]])
                if u.sum() < 0.5 * rows_sel.sum():
                    continue
                vals.append((ww * K[q][u]).sum() / (ww * km).sum())
            out[f'{nul}_{kind}'] = np.array(vals)
    return out


def fit(args):
    A = dict(np.load(os.path.join(OUT, 'apertures.npz')))
    meta = A['meta']; sets = list(A['sets'])
    Kd = np.load(os.path.join(OUT, f'K_data_{args.map}.npy'))
    sf = sorted(glob.glob(os.path.join(OUT, 'sims', 'K_*.npy')))
    Ks = np.array([np.load(f) for f in sf])
    log(f'map {args.map}; {len(sf)} simulations')
    N = A['row_N']; z = A['row_z']; nrow = len(N)
    # weights fixed from the first 100 simulations: w = K_m / var(K) per supercluster and aperture (uniform if no sims yet)
    w = {}
    Km = A['K1'] + A['K2']
    wfn = os.path.join(OUT, 'weights_var100.npz')
    if os.path.exists(wfn):
        wv = np.load(wfn)
    else:
        if len(Ks) < 100:
            raise SystemExit('the weights need 100 simulations')
        wv = {k: Ks[:100][:, meta[:, 0] == sets.index(k)].var(0) for k in ('hull', 'circle')}
        np.savez(wfn, sims=np.array([os.path.basename(f) for f in sf[:100]]), **wv)
        log(f'weights frozen from {sf[0]} ... {sf[99]}')
    for kind in ('hull', 'circle'):
        s = meta[:, 0] == sets.index(kind)
        var = wv[kind]
        ww = np.zeros(nrow); ww[meta[s, 1]] = Km[s] / var
        if args.uniform:
            ww[:] = 1.0
        w[kind] = ww
    samples = {'>= 2 members': N >= 2, '>= 3 members': N >= 3, 'pairs': N == 2, '>= 3, z < 0.2': (N >= 3) & (z < 0.2),
               '>= 3, z > 0.2': (N >= 3) & (z >= 0.2), '>= 5 members': N >= 5}
    apok = None
    if args.map == 'planck':
        apok = np.load(os.path.join(OUT, 'planck_cleanfrac.npy')) >= 0.98
        log(f'planck: apertures fully inside the Planck mask: {apok.mean():.3f}')
    base = {}
    if args.map != 'baseline' and os.path.exists(os.path.join(OUT, 'K_data_baseline.npy')):
        Kb = np.load(os.path.join(OUT, 'K_data_baseline.npy'))
    lines = []
    for model in ('member', 'all'):
        for name, rs in samples.items():
            if apok is not None:
                rs = rs.copy()
                for kind in ('hull', 'circle'):
                    sk = meta[:, 0] == sets.index(kind)
                    bad = np.zeros(nrow, bool); bad[meta[sk, 1][~apok[sk]]] = True
                    rs &= ~bad
            d = stats(Kd, A, rs, w, model, apok)
            if args.map != 'baseline':
                db = stats(Kb, A, rs, w, model, apok)
            sim = [stats(k, A, rs, w, model, apok) for k in Ks] if len(Ks) else []
            for kind in ('hull', 'circle'):
                e, c = d[kind], d['ctrl_' + kind]
                if sim:
                    se = np.std([q[kind] for q in sim]); sc = np.std([q['ctrl_' + kind] for q in sim])
                    sdiff = np.std([q[kind] - q['ctrl_' + kind] for q in sim]); mse = np.mean([q[kind] for q in sim])
                    spop = np.std([q[kind] - q['pop_' + kind] for q in sim])
                else:
                    se = sc = sdiff = mse = spop = np.nan
                rnd, rot = d['rand_' + kind], d['rot_' + kind]
                if args.map != 'baseline':
                    # no simulations of this map here: the baseline's simulation errors scaled by the ratio of the two maps'
                    # random-position null spreads on the same apertures
                    rr = rnd.std() / db['rand_' + kind].std()
                    se, sc, sdiff, spop = se * rr, sc * rr, sdiff * rr, spop * rr
                lines.append(f'{model:6s} {name:16s} {kind:6s} n={rs.sum():4d}  E = {e:+.3f} ± {se:.3f}  control = {c:+.3f} ± {sc:.3f} '
                             f'(realisations spread {d["ctrl_" + kind + "_real"].std():.3f})  E - control = {e - c:+.3f} ± {sdiff:.3f} '
                             f'({(e - c) / sdiff:+.1f}σ)  sims mean {mse + 1:+.3f}; random null {rnd.mean():+.3f} ± {rnd.std():.3f} '
                             f'(n {len(rnd)}), rotated {rot.mean():+.3f} ± {rot.std():.3f} (n {len(rot)}); population control '
                             f'{d["pop_" + kind]:+.3f}, E - pop {e - d["pop_" + kind]:+.3f} ± {spop:.3f}')
                h, hc = d['h1_' + kind], d['h1_ctrl_' + kind]
                # the share of the aperture's band-filtered lensing in the member halos (one-halo), h / (1 + E), its 68 % range
                # from E's simulation error
                lines.append(f'       {"":16s} {kind:6s} one-halo share of the model {h:.3f} (control {hc:.3f}); of the measured '
                             f'integral {h / (1 + e):.3f} [{h / (1 + e + se):.3f}, {h / max(1 + e - se, 1e-3):.3f}] (control '
                             f'{hc / (1 + c):.3f} [{hc / (1 + c + sc):.3f}, {hc / max(1 + c - sc, 1e-3):.3f}]); K_d/K_m = {1 + e:.2f} '
                             f'({(1 + e) / se:.1f}σ above zero)')
                log(lines[-2]); log(lines[-1])
    with open(os.path.join(OUT, f'fit_{args.map}{"_uniform" if args.uniform else ""}.txt'), 'w') as f:
        f.write('\n'.join(lines) + '\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd'); ap.add_argument('--map', default='baseline')
    ap.add_argument('--worker', type=int, default=0); ap.add_argument('--nworkers', type=int, default=1)
    ap.add_argument('--uniform', action='store_true')
    a = ap.parse_args()
    {'prep': prep, 'data': data, 'sims': sims, 'fit': fit}[a.cmd](a)


if __name__ == '__main__':
    main()
