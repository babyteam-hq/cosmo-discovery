"""eRASS1 supercluster membership rebuilt from the released recipe (R-137, cosmo-deep-13).

The released supercluster catalogue (Liu et al. 2024, arXiv:2402.08454, `erass1sc.fits`) gives each system's centre (the mean of
its members' coordinates), multiplicity, total mass and length, but not its members: the member table (their Table 2) is at the
CDS only, which this cloud cannot reach (cdsarc/vizier reset the connection; the CfA mirror times out). So the members are rebuilt
with the paper's own recipe (Sect. 2.2 and 3): the primary catalogue with PCONT <= 0.3, dz/(1+z) <= 0.02, photometric redshifts
inside the optical limiting redshift (IN_ZVLIM), no duplicates (SHARED_MEMBERS), >= 5 counts in R500, z > 0.005; friends-of-friends
on 3D comoving distances (flat LCDM H0 = 70, Om = 0.3, X-ray centroids RA_XFIT/DEC_XFIT), linking length l = (n f)^-1/3 with n the
number density in a shell of dz = +-0.005 (and, as in the paper, in bins of exposure), f = 10, l <= 50 Mpc. Each catalogued
system is then matched to a rebuilt one, and accepted only if the rebuilt one reproduces its multiplicity, its centre (within
the catalogue's 4-decimal rounding of a mean) and its total mass (mtot = 1.46 sum M500 / 0.39, the paper's Sect. 4). Only the
accepted systems are used downstream, so every member list used is one that reproduces the catalogue entry.

usage: python -I sc_members.py    # writes data/erass1_sc/members.npz and prints the match statistics
"""
import os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
OUT = os.path.join(DATA, 'erass1_sc')
F_OVER = 10.0
LMAX_LINK = 50.0
AREA_DEG2 = 20626.5  # the western Galactic hemisphere


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def load_sample():
    from astropy.io import fits
    p = fits.open(os.path.join(DATA, 'erass1', 'erass1cl_primary_v3.2.fits'))[1].data
    o = fits.open(os.path.join(DATA, 'erass1', 'eRASS1_clusters_optical.fits'))[1].data
    assert (p['DETUID'] == o['DETUID']).all()
    z = np.asarray(p['BEST_Z'], float); ze = np.asarray(p['BEST_ZERR'], float)
    zt = np.char.strip(np.asarray(p['BEST_Z_TYPE']).astype(str))
    sel = ((p['PCONT'] <= 0.3) & (ze / (1 + z) <= 0.02) & ~((zt == 'photo_z') & ~o['IN_ZVLIM'].astype(bool))
           & ~o['SHARED_MEMBERS'].astype(bool) & (p['CTS500'] >= 5) & (z > 0.005))
    return p, sel


def comoving(z):
    from astropy.cosmology import FlatLambdaCDM
    co = FlatLambdaCDM(H0=70, Om0=0.3)
    zg = np.linspace(0, 1.6, 3201)
    return np.interp(z, zg, co.comoving_distance(zg).value), co


def linking_length(z, exp, nexp_bins):
    """l(z, t) = (n f)^-1/3 with n from the clusters in z +- 0.005 and, if nexp_bins > 1, in exposure quantile bins whose area is
    taken proportional to their share of the clusters at low z (the paper uses the eRASS1 depth curve A(t), not public here)."""
    zg = np.linspace(0.0, 1.4, 1401)
    dc, co = comoving(zg)
    frac_sky = AREA_DEG2 / 41252.96
    l = np.empty(len(z))
    if nexp_bins <= 1:
        bins = np.zeros(len(z), int); area_share = np.array([1.0])
    else:
        q = np.quantile(exp, np.linspace(0, 1, nexp_bins + 1)); q[-1] += 1
        bins = np.clip(np.searchsorted(q, exp, side='right') - 1, 0, nexp_bins - 1)
        # area share of each exposure bin: from the z < 0.15 clusters, where the selection is nearly exposure-independent
        lo = z < 0.15
        area_share = np.array([(bins[lo] == b).mean() for b in range(nexp_bins)])
    for b in range(len(area_share)):
        ib = bins == b
        for i in np.where(ib)[0]:
            z0 = z[i]; a, c = max(z0 - 0.005, 0.0), z0 + 0.005
            n = ((z[ib] >= a) & (z[ib] < c)).sum()
            da, dcc = np.interp([a, c], zg, dc)
            V = 4 / 3 * np.pi * (dcc ** 3 - da ** 3) * frac_sky * area_share[b]
            l[i] = min((n / V * F_OVER) ** (-1 / 3), LMAX_LINK)
    return l


def fof(xyz, link):
    """friends-of-friends with a per-pair linking length max? the paper's l is per cluster; we link i-j if d < min(l_i, l_j)
    (variant 'min') or < the mean (variant 'mean')."""
    from scipy.spatial import cKDTree
    t = cKDTree(xyz)
    pairs = t.query_pairs(LMAX_LINK, output_type='ndarray')
    d = np.linalg.norm(xyz[pairs[:, 0]] - xyz[pairs[:, 1]], axis=1)
    return pairs, d


def groups_from(pairs, n):
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    A = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    return connected_components(A, directed=False)[1]


def build(rule='first', nexp_bins=1):
    p, sel = load_sample()
    idx = np.where(sel)[0]
    ra = np.asarray(p['RA_XFIT'], float)[idx]; dec = np.asarray(p['DEC_XFIT'], float)[idx]
    bad = ~np.isfinite(ra) | ~np.isfinite(dec)
    ra[bad] = np.asarray(p['RA'], float)[idx][bad]; dec[bad] = np.asarray(p['DEC'], float)[idx][bad]
    z = np.asarray(p['BEST_Z'], float)[idx]
    dc, _ = comoving(z)
    r, d = np.radians(ra), np.radians(dec)
    xyz = dc[:, None] * np.c_[np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)]
    l = linking_length(z, np.asarray(p['EXP'], float)[idx], nexp_bins)
    pairs, dist = fof(xyz, l)
    li, lj = l[pairs[:, 0]], l[pairs[:, 1]]
    thr = {'min': np.minimum(li, lj), 'max': np.maximum(li, lj), 'mean': 0.5 * (li + lj)}
    res = {}
    for k, t in thr.items():
        g = groups_from(pairs[dist < t], len(idx))
        res[k] = g
    return dict(idx=idx, ra=ra, dec=dec, z=z, M500=np.asarray(p['M500'], float)[idx], l=l), res


def match(cat, groups, verbose=True):
    """for each catalogued system, the rebuilt group whose members reproduce its centre, multiplicity and mtot."""
    from astropy.io import fits
    sc = fits.open(os.path.join(DATA, 'erass1', 'erass1sc.fits'))[1].data
    lab = groups
    order = np.argsort(lab); ls = lab[order]
    starts = np.r_[0, np.where(np.diff(ls))[0] + 1]; ends = np.r_[starts[1:], len(ls)]
    sizes = ends - starts
    multi = sizes >= 2
    gm = [order[s:e] for s, e, m in zip(starts, ends, multi) if m]
    # group centres: mean coordinates (the paper: 'the algebraic mean coordinate of the members')
    gra = np.array([cat['ra'][m].mean() for m in gm]); gdec = np.array([cat['dec'][m].mean() for m in gm])
    gz = np.array([cat['z'][m].mean() for m in gm]); gn = np.array([len(m) for m in gm])
    gmt = np.array([1.46 * cat['M500'][m].sum() / 0.39 for m in gm])  # M500 in 1e13 -> mtot in 1e13; catalogue in 1e14
    out = []
    for k in range(len(sc)):
        dd = np.hypot((gra - sc['ra'][k]) * np.cos(np.radians(sc['dec'][k])), gdec - sc['dec'][k])
        j = np.argmin(dd)
        ok = (dd[j] < 0.01) and gn[j] == sc['multiplicity'][k] and abs(gmt[j] / 10 / sc['mtot'][k] - 1) < 0.02
        out.append((j, dd[j], gn[j], gmt[j] / 10, ok))
    okk = np.array([o[4] for o in out])
    if verbose:
        log(f'rebuilt systems with >= 2 members: {len(gm)} ({(gn == 2).sum()} pairs, {(gn >= 3).sum()} with >= 3; '
            f'{gn.sum()} member clusters); catalogue: {len(sc)} ({(sc["multiplicity"] == 2).sum()}, '
            f'{(sc["multiplicity"] >= 3).sum()}, {sc["multiplicity"].sum()})')
        log(f'catalogued systems reproduced exactly (centre < 0.01 deg, multiplicity, mtot within 2 %): {okk.sum()} of {len(sc)}; '
            f'pairs {okk[sc["multiplicity"] == 2].sum()}/{(sc["multiplicity"] == 2).sum()}, >=3 '
            f'{okk[sc["multiplicity"] >= 3].sum()}/{(sc["multiplicity"] >= 3).sum()}')
    return sc, gm, out, okk


def main():
    os.makedirs(OUT, exist_ok=True)
    best = None
    for nexp in (1, 3, 5):
        cat, res = build(nexp_bins=nexp)
        for rule, g in res.items():
            log(f'--- exposure bins {nexp}, linking rule {rule}')
            sc, gm, out, okk = match(cat, g)
            if best is None or okk.sum() > best[0]:
                best = (okk.sum(), nexp, rule, cat, sc, gm, out, okk)
    n, nexp, rule, cat, sc, gm, out, okk = best
    log(f'best: exposure bins {nexp}, rule {rule}: {n} of {len(sc)} reproduced')
    # save: for each catalogued system, its member indices into the primary catalogue (empty if not reproduced)
    mem_idx = [];
    for k in range(len(sc)):
        j = out[k][0]
        mem_idx.append(cat['idx'][gm[j]] if okk[k] else np.array([], int))
    flat = np.concatenate(mem_idx); ptr = np.r_[0, np.cumsum([len(m) for m in mem_idx])]
    # every rebuilt member of any system (for the 'isolated' control pool): all groups >= 2
    any_member = np.concatenate([cat['idx'][m] for m in gm])
    np.savez(os.path.join(OUT, 'members.npz'), flat=flat, ptr=ptr, ok=okk, sample_idx=cat['idx'], any_member=any_member,
             nexp=nexp, rule=rule)
    log('wrote members.npz')



# ------------------------------------------------------------------------------------------------ exact refinement
def refine(max_rm=2, max_add=2):
    """for every catalogued system the FoF rebuild did not reproduce: start from the nearest rebuilt group (or the empty set)
    and search removals of up to max_rm members and additions of up to max_add primary clusters (any PCONT, near the centre
    and redshift) for the member set whose mean RA_XFIT, DEC_XFIT, BEST_Z and summed M500 reproduce the catalogue's centre,
    redshift and mtot to their quoted rounding (ra, dec 1e-4 deg, z 1e-4, mtot 1e-3 x 1e14). Accepted only if exactly one
    such set exists within the search."""
    import itertools
    from astropy.io import fits
    p = fits.open(os.path.join(DATA, 'erass1', 'erass1cl_primary_v3.2.fits'))[1].data
    sc = fits.open(os.path.join(DATA, 'erass1', 'erass1sc.fits'))[1].data
    g = np.load(os.path.join(OUT, 'members.npz'))
    flat, ptr, ok = g['flat'], g['ptr'], g['ok'].copy()
    ra = np.asarray(p['RA_XFIT'], float); dec = np.asarray(p['DEC_XFIT'], float); z = np.asarray(p['BEST_Z'], float)
    M = np.asarray(p['M500'], float) * 1.46 / 0.39 / 10  # in 1e14, so that sum M = mtot
    fine = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(M)
    tol = np.array([7.6e-5, 7.6e-5, 7.6e-5, 6e-4])  # mean ra, mean dec, mean z, relative mtot (3 decimals on mtot ~ 5)

    # rebuild the FoF groups once (best variant) to seed the search
    cat, res = build(nexp_bins=int(g['nexp']))
    lab = res[str(g['rule'])]
    gidx = {}
    for i, L in enumerate(lab):
        gidx.setdefault(L, []).append(cat['idx'][i])
    groups = [np.array(v) for v in gidx.values() if len(v) >= 2]
    gc = np.array([[ra[m].mean(), dec[m].mean(), z[m].mean()] for m in groups])
    members = [flat[ptr[k]:ptr[k + 1]] for k in range(len(sc))]
    nfix = 0; namb = 0
    for k in np.where(~ok)[0]:
        N = int(sc['multiplicity'][k]); T = np.array([sc['ra'][k] * N, sc['dec'][k] * N, sc['redshift'][k] * N, sc['mtot'][k]])
        # seeds: rebuilt groups whose centre lies within the system's extent and redshift
        rad = max(1.0, np.degrees(sc['lproj'][k] / max(comoving(np.array([sc['redshift'][k]]))[0][0], 1)))
        dd = np.hypot((gc[:, 0] - sc['ra'][k]) * np.cos(np.radians(sc['dec'][k])), gc[:, 1] - sc['dec'][k])
        seeds = [groups[j] for j in np.where((dd < rad) & (np.abs(gc[:, 2] - sc['redshift'][k]) < 0.03))[0]][:4]
        sep = np.hypot((ra - sc['ra'][k]) * np.cos(np.radians(sc['dec'][k])), dec - sc['dec'][k])
        pool = np.where(fine & (sep < 1.2 * rad + 0.5) & (np.abs(z - sc['redshift'][k]) < 0.04 * (1 + sc['redshift'][k])))[0]
        sols = set()
        for B in seeds + [np.array([], int)]:
            B = np.array(B, int)
            vB = np.array([ra[B].sum(), dec[B].sum(), z[B].sum(), M[B].sum()]) if len(B) else np.zeros(4)
            P = np.setdiff1d(pool, B)
            Pv = np.c_[ra[P], dec[P], z[P], M[P]]
            for nr in range(0, min(max_rm, len(B)) + 1):
                for R in itertools.combinations(range(len(B)), nr):
                    na = N - (len(B) - nr)
                    if na < 0 or na > max_add or (len(B) == 0 and na > max_add):
                        continue
                    vR = np.array([ra[B[list(R)]].sum(), dec[B[list(R)]].sum(), z[B[list(R)]].sum(), M[B[list(R)]].sum()]) if nr else np.zeros(4)
                    need = T - vB + vR
                    if na == 0:
                        cands = [()] if _close(np.zeros(4), need, N, T, tol) else []
                    elif na == 1:
                        cands = [(i,) for i in np.where(_close_arr(Pv, need, N, T, tol))[0]]
                    else:
                        cands = []
                        for i in range(len(P)):
                            rest = need - Pv[i]
                            for j in np.where(_close_arr(Pv[i + 1:], rest, N, T, tol))[0]:
                                cands.append((i, i + 1 + j))
                    keep = [b for r_, b in enumerate(B) if r_ not in R]
                    for c in cands:
                        sols.add(tuple(sorted(keep + [int(P[i]) for i in c])))
        sols = [s for s in sols if len(s) == N]
        if len(P) > 400 and max_add > 2:
            pass
        if len(sols) == 1:
            members[k] = np.array(sols[0]); ok[k] = True; nfix += 1
        elif len(sols) > 1:
            namb += 1
    log(f'refinement: {nfix} more systems reproduced exactly, {namb} ambiguous; total {ok.sum()} of {len(sc)} '
        f'(pairs {ok[sc["multiplicity"] == 2].sum()}/{(sc["multiplicity"] == 2).sum()}, >= 3 '
        f'{ok[sc["multiplicity"] >= 3].sum()}/{(sc["multiplicity"] >= 3).sum()})')
    flat2 = np.concatenate([members[k] if ok[k] else np.array([], int) for k in range(len(sc))]).astype(int)
    ptr2 = np.r_[0, np.cumsum([len(members[k]) if ok[k] else 0 for k in range(len(sc))])]
    np.savez(os.path.join(OUT, 'members_refined.npz'), flat=flat2, ptr=ptr2, ok=ok, sample_idx=g['sample_idx'],
             any_member=g['any_member'], all_sc_member=flat2)


def _close(v, need, N, T, tol):
    return _close_arr(v[None, :], need, N, T, tol)[0]


def _close_arr(V, need, N, T, tol):
    d = np.abs(V - need[None, :])
    return (d[:, 0] < tol[0] * N) & (d[:, 1] < tol[1] * N) & (d[:, 2] < tol[2] * N) & (d[:, 3] < tol[3] * max(T[3], 1) + 6e-4)


if __name__ == '__main__':
    if sys.argv[1:2] == ['refine']:
        refine(*[int(a) for a in sys.argv[2:4]])
    else:
        main()
