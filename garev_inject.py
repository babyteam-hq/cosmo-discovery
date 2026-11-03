"""cosmo-ga-rev (R-171): selection shrinkage of the Giant Arc by injection into the real SDSS catalogue
(plan and the method check's seven points: `ga_rev_plan.md`).

A real overdensity of amplitude A_true is injected at random places in the SDSS footprint in one of three shapes: H, uniform
in the GA's own CHMS shape; S, uniform in its 30.4 Mpc member spheres only (a clumpy arc); L, uniform in the hull of 60.8 Mpc
member spheres (larger, more diffuse). Single linkage at 95 Mpc is re-run on the real Z&M slice plus the injected absorbers.
The recovered group's member overdensity (delta_mem, the estimator also run on the GA rebuild) decides whether it is
GA-like. For GA-like recoveries two DESI quantities are computed on real DESI DR1 data:
  f_rec: the share of the DESI expectation in the recovered hull that lies in the injected region, with both rotated back
         onto the GA's own sightlines (a real structure of amplitude A_true adds A_true * f_rec to the hull's delta);
  bg:    DESI N_obs and N_exp in the recovered hull at the injection site itself (DESI does not see the injected
         absorbers, so this is the chance part of a GA-like hull under that hypothesis, measured, not assumed).
Completeness models for the SDSS placement: "u" two classes (BOSS searched, DR7 only), each with one slice dN/dz; "g"
eight g-magnitude classes from DR16Q PSF magnitudes (the variant).

Run from the repo root after r155_template.py and r155_sightlines.py:
  python -u research/cosmo-discovery/garev_inject.py calib
  nohup python -u research/cosmo-discovery/garev_inject.py run <nsite> <nproc> > .../logs/garev_run.log &
"""
import sys
import time

import numpy as np
from astropy.io import fits
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import ConvexHull, cKDTree
from scipy.spatial.distance import pdist

sys.path.insert(0, "research/cosmo-discovery")
from r155_common import C, D, DATA, sky_match, unit, xyz, zofD  # noqa: E402
from r155_core import Sky, clip_D, fib_sphere, intersect, mst_mean_edge, rotation_to  # noqa: E402

ZC, DZ, LINK = 0.802, 0.060, 95.0
Z1, Z2 = ZC - DZ, ZC + DZ
D1, D2 = D(Z1), D(Z2)
A_GRID = [0.0, 0.25, 0.5, 0.75, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 4.0]
GA_SEP = 742.0
# runs: (profile, completeness model)
RUNS = [("H", "u"), ("S", "u"), ("L", "u"), ("H", "g")]
# loose pre-filter for computing the DESI quantities; the GA-like window itself is applied in garev_infer.py
PRE_SIZE, PRE_SEP, PRE_DMEM = (25, 80), (0.4, 2.2), 0.6


class SDSSRays:
    def __init__(self, sl):
        q7 = fits.getdata(f"{DATA}/Searched_Quasars_SDSS_DR7_107.fits", 1)
        q12 = fits.getdata(f"{DATA}/Searched_Quasars_BOSS_DR12_107.fits", 1)
        r7, d7, z7 = (np.asarray(q7[c], float) for c in ("RA", "DEC", "ZQSO"))
        r12, d12, z12 = (np.asarray(q12[c], float) for c in ("RA", "DEC", "ZQSO"))
        i7, _ = sky_match(r7, d7, r12, d12, 1.5)
        only7 = np.ones(len(r7), bool)
        only7[i7] = False
        self.ra = np.concatenate([r12, r7[only7]])
        self.dec = np.concatenate([d12, d7[only7]])
        self.zq = np.concatenate([z12, z7[only7]])
        boss = np.concatenate([np.ones(len(r12), int), np.zeros(only7.sum(), int)])
        self.zlo = np.maximum(Z1, 1216 * (1 + self.zq) / 2796 - 1)
        self.zhi = np.minimum(Z2, (1 + self.zq) * (1 - 5000 / C) - 1)
        self.Dmin, self.Dmax = D(self.zlo), D(np.maximum(self.zhi, self.zlo))
        self.u = unit(self.ra, self.dec)
        self.tree = cKDTree(self.u)
        ia, js = sky_match(sl["ra"], sl["dec"], self.ra, self.dec, 1.5)
        on = np.full(len(sl["ra"]), -1)
        on[ia] = js
        # g magnitude from DR16Q (a completeness proxy: the spectrum's S/N at the slice's Mg II, ~5000 A)
        q = fits.open(f"{DATA}/DR16Q_v4.fits", memmap=True)[1].data
        i1, i2 = sky_match(self.ra, self.dec, np.asarray(q["RA"], float), np.asarray(q["DEC"], float), 1.5)
        g = np.full(len(self.ra), np.nan)
        g[i1] = np.asarray(q["PSFMAG"][:, 1], float)[i2]
        g[(g < 10) | (g > 30)] = np.nan
        ok = np.isfinite(g)
        gb = np.full(len(g), 7)
        gb[ok] = np.digitize(g[ok], np.nanpercentile(g, [100 / 7 * k for k in range(1, 7)]))
        self.cls = dict(u=boss, g=gb)
        path = np.clip(self.zhi - self.zlo, 0, None)
        self.dndz = {}
        for m, cl in self.cls.items():
            k = cl.max() + 1
            na = np.bincount(cl[on[on >= 0]], minlength=k)
            pa = np.bincount(cl, weights=path, minlength=k)
            self.dndz[m] = na / pa
        inpath = (sl["z"][on >= 0] >= self.zlo[on[on >= 0]] - 1e-6) & (sl["z"][on >= 0] <= self.zhi[on[on >= 0]] + 1e-6)
        self.report = dict(nsl=len(self.ra), nboss=int(boss.sum()), matched=round(float((on >= 0).mean()), 4),
                           inpath=round(float(inpath.mean()), 4), g_matched=round(float(ok.mean()), 3),
                           dndz_u=np.round(self.dndz["u"], 4).tolist(), dndz_g=np.round(self.dndz["g"], 4).tolist(),
                           nexp_all=round(float((path * self.dndz["u"][boss]).sum()), 1), nobs_all=len(sl["ra"]))

    hull_intervals = Sky.hull_intervals
    sphere_intervals = Sky.sphere_intervals

    def expect(self, iv, model="u"):
        si, a1, a2 = iv[:3]
        return (zofD(a2) - zofD(a1)) * self.dndz[model][self.cls[model][si]]


def hull_eq(centres, r):
    pts = (centres[:, None, :] + fib_sphere(200)[None] * r).reshape(-1, 3)
    return ConvexHull(pts).equations


def in_hull(X, eq):
    return (X @ eq[:, :3].T + eq[:, 3]).max(1) <= 0


def link(X):
    pairs = cKDTree(X).query_pairs(LINK, output_type="ndarray")
    n = len(X)
    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    return connected_components(g, directed=False)[1]


def delta_mem(rays, X_all, members, model="u"):
    """All slice absorbers in the group's CHMS (radius half its mean MST edge) over the SDSS expectation in
    CHMS n slice, minus 1."""
    r = 0.5 * mst_mean_edge(members)
    e = rays.expect(clip_D(rays.hull_intervals(members, r)[:3], D1, D2), model).sum()
    n = in_hull(X_all, hull_eq(members, r)).sum()
    return n / e - 1, n, e, r


def region_iv(obj, G, P, r_ga):
    if P == "H":
        return clip_D(obj.hull_intervals(G, r_ga)[:3], D1, D2)
    if P == "L":
        return clip_D(obj.hull_intervals(G, 2 * r_ga)[:3], D1, D2)
    return clip_D(obj.sphere_intervals(G, r_ga), D1, D2)


def region_mask(X, G, P, r_ga):
    if P == "S":
        return cKDTree(G).query(X, distance_upper_bound=r_ga)[0] < np.inf
    return in_hull(X, hull_eq(G, r_ga if P == "H" else 2 * r_ga))


def load():
    sl = dict(np.load(f"{DATA}/zm_slice.npz"))
    ga = dict(np.load(f"{DATA}/ga_members.npz"))
    return sl, ga, xyz(sl["ra"], sl["dec"], D(sl["z"])), xyz(ga["ra"], ga["dec"], D(ga["z"]))


def desi_delta(sky, M, r, c=0.3):
    k = sky.counts(sky.hull_intervals(M, r)[:3], c)
    return k["nobs"], k["nexp"]


def calib(sky=None):
    sl, ga, X, G = load()
    rays = SDSSRays(sl)
    print("SDSS rays:", rays.report)
    out = {}
    for m in ("u", "g"):
        d, n, e, r = delta_mem(rays, X, G, m)
        out[m] = d
        print(f"GA rebuild ({m}): delta_mem = {d:.3f} (N {n} in CHMS n slice, E_SDSS {e:.2f}, r_sph {r:.1f} Mpc); "
              f"published 1.3 +- 0.3")
    lab = link(X)
    sizes = np.bincount(lab)
    rows = []
    for gid in np.nonzero(sizes >= 20)[0]:
        M = X[lab == gid]
        sep = pdist(M).max()
        dm, n, e, r = delta_mem(rays, X, M)
        row = [gid, sizes[gid], sep, dm, np.nan, np.nan]
        if sky is not None and 25 <= sizes[gid] <= 80:
            row[4:] = desi_delta(sky, M, r)
        rows.append(row)
    rows = np.array(rows)
    like = (rows[:, 1] >= 30) & (rows[:, 1] <= 70) & (rows[:, 2] >= 0.5 * GA_SEP) & (rows[:, 2] <= 2 * GA_SEP)
    print(f"slice groups >= 20 members: {len(rows)}; 30-70 members and GA-like extent: {like.sum()}, delta_mem median "
          f"{np.median(rows[like, 3]):.2f} (range {rows[like, 3].min():.2f} to {rows[like, 3].max():.2f}); within 0.3 of "
          f"the GA's: {(like & (np.abs(rows[:, 3] - out['u']) <= 0.3)).sum()}")
    if sky is not None:
        nob, nex = desi_delta(sky, G, 0.5 * mst_mean_edge(G))
        print(f"GA on DESI (whole CHMS, W > 0.3): N_obs {nob:.0f}, N_exp {nex:.2f}")
    return rays, sl, X, G, out, lab, rows


# ---------------------------------------------------------------------------------------------------------------------
_G = {}


def one(args):
    k, A, P, model, seed = args
    g = _G
    rng = np.random.default_rng(seed)
    R = g["sites"][k]
    Gr = g["G"] @ R.T
    rays = g["rays"]
    iv = region_iv(rays, Gr, P, g["r_ga"])
    e = rays.expect(iv, model)
    E = e.sum()
    nin = rng.poisson(A * E) if A > 0 else 0
    if nin:
        pick = rng.choice(len(e), nin, p=e / E)
        si, a1, a2 = iv[0][pick], iv[1][pick], iv[2][pick]
        zz = zofD(a1) + rng.random(nin) * (zofD(a2) - zofD(a1))
        Xall = np.vstack([g["X"], rays.u[si] * D(zz)[:, None]])
    else:
        Xall = g["X"]
    inj = np.r_[np.zeros(len(g["X"]), bool), np.ones(nin, bool)]
    lab = link(Xall)
    sizes = np.bincount(lab)
    if nin:
        score = np.bincount(lab[inj], minlength=len(sizes)).astype(float)
    else:  # chance: the group with the most members inside the region
        score = np.bincount(lab[region_mask(Xall, Gr, P, g["r_ga"])], minlength=len(sizes)).astype(float)
    score[sizes < 20] = -1
    gid = int(np.argmax(score))
    row = dict(k=k, A=A, P=P, model=model, E_sdss=E, nin=nin, found=bool(score[gid] > 0), size=0, sep=0.0,
               dmem=np.nan, ninj=0, merged=0, f_rec=np.nan, bg_obs=np.nan, bg_exp=np.nan, bg6_obs=np.nan,
               bg6_exp=np.nan, f_rec6=np.nan)
    if not row["found"]:
        return row
    inG = lab == gid
    M = Xall[inG]
    # merges: members drawn from one pre-existing real group of >= 20 members
    old = g["lab0"][inG[:len(g["X"])]]
    old = old[g["size0"][old] >= 20]
    row.update(size=int(len(M)), sep=float(pdist(M).max()), ninj=int((inj & inG).sum()),
               merged=int(np.bincount(old).max()) if len(old) else 0)
    if not (PRE_SIZE[0] <= len(M) <= PRE_SIZE[1] and PRE_SEP[0] * GA_SEP <= row["sep"] <= PRE_SEP[1] * GA_SEP):
        return row
    dm, n, ee, rr = delta_mem(rays, Xall, M, model)
    row["dmem"] = float(dm)
    if abs(dm - g["dGA"][model]) > PRE_DMEM:
        return row
    sky = g["sky"]
    # chance part at the site (real DESI absorbers in the recovered hull; the injected ones are not in DESI)
    for c, ko, ke in ((0.3, "bg_obs", "bg_exp"), (0.6, "bg6_obs", "bg6_exp")):
        row[ko], row[ke] = desi_delta(sky, M, rr, c)
    # rotate back onto the GA's sightlines: the injected region's share of the recovered hull's DESI expectation
    Mb = M @ R          # R^-1 = R^T
    rec = sky.hull_intervals(Mb, rr)[:3]
    reg = region_iv(sky, g["G"], P, g["r_ga"])
    ov = intersect(rec, reg)
    for c, kf in ((0.3, "f_rec"), (0.6, "f_rec6")):
        e_rec = sky.counts(rec, c)["nexp"]
        e_ov = sky.counts(ov, c)["nexp"] if len(ov[0]) else 0.0
        row[kf] = float(e_ov / e_rec) if e_rec > 0 else np.nan
    return row


def make_sites(rays, G, r_ga, E_ga, X, lab0, size0, n, rng):
    """Rigid moves of the GA's CHMS to random SDSS sightlines' directions: > 25 deg from the GA, SDSS expectation
    0.75-1.33x the GA's, and the moved hull holding < 10 members of any real slice group of >= 40 members."""
    cu = G.mean(0) / np.linalg.norm(G.mean(0))
    big = size0[lab0] >= 40
    sites, tries = [], 0
    while len(sites) < n and tries < 300 * n:
        tries += 1
        tgt = rays.u[rng.integers(len(rays.ra))]
        if np.degrees(np.arccos(np.clip(tgt @ cu, -1, 1))) < 25:
            continue
        R = rotation_to(cu, tgt, rng.uniform(0, 2 * np.pi))
        Gr = G @ R.T
        e = rays.expect(clip_D(rays.hull_intervals(Gr, r_ga)[:3], D1, D2)).sum()
        if not 0.75 * E_ga <= e <= 1.33 * E_ga:
            continue
        inside = in_hull(X[big], hull_eq(Gr, 2 * r_ga))
        if inside.any() and np.bincount(lab0[big][inside]).max() >= 10:
            continue
        sites.append(R)
    return sites, tries


def run(nsite=200, nproc=3):
    t0 = time.time()
    sky = Sky()
    print(f"DESI sky loaded; {time.time() - t0:.0f} s", flush=True)
    rays, sl, X, G, dGA, lab0, grp = calib(sky)
    np.save(f"{DATA}/garev_groups.npy", grp)
    size0 = np.bincount(lab0)
    r_ga = 0.5 * mst_mean_edge(G)
    E_ga = rays.expect(clip_D(rays.hull_intervals(G, r_ga)[:3], D1, D2)).sum()
    sites, tries = make_sites(rays, G, r_ga, E_ga, X, lab0, size0, nsite, np.random.default_rng(171))
    cen = np.array([(G.mean(0) @ R.T) / np.linalg.norm(G.mean(0)) for R in sites])
    print(f"sites: {len(sites)} accepted of {tries} tries (E_SDSS 0.75-1.33x the GA's {E_ga:.2f}); "
          f"{time.time() - t0:.0f} s", flush=True)
    np.save(f"{DATA}/garev_site_centres.npy", cen)
    _G.update(rays=rays, X=X, G=G, r_ga=r_ga, sites=sites, dGA=dGA, sky=sky, lab0=lab0, size0=size0)
    jobs = [(k, A, P, m) for P, m in RUNS for A in A_GRID for k in range(len(sites))]
    jobs = [j + (1000003 * i + 17,) for i, j in enumerate(jobs)]
    rows = []
    import multiprocessing as mp
    with mp.get_context("fork").Pool(nproc) as pool:
        for i, r in enumerate(pool.imap_unordered(one, jobs, chunksize=8)):
            rows.append(r)
            if (i + 1) % 500 == 0:
                print(f"{i + 1}/{len(jobs)} done, {time.time() - t0:.0f} s", flush=True)
                save(rows, dGA, E_ga, r_ga)
    save(rows, dGA, E_ga, r_ga)
    print(f"wrote {DATA}/garev_runs.npz, {len(rows)} rows, {time.time() - t0:.0f} s")


def save(rows, dGA, E_ga, r_ga):
    np.savez(f"{DATA}/garev_runs.npz", dGA_u=dGA["u"], dGA_g=dGA["g"], E_ga=E_ga, r_ga=r_ga,
             **{k: np.array([r[k] for r in rows]) for k in rows[0]})


if __name__ == "__main__":
    if sys.argv[1] == "calib":
        calib()
    elif sys.argv[1] == "time":
        sky = Sky()
        rays, sl, X, G, dGA, lab0, grp = calib(sky)
        size0 = np.bincount(lab0)
        r_ga = 0.5 * mst_mean_edge(G)
        E_ga = rays.expect(clip_D(rays.hull_intervals(G, r_ga)[:3], D1, D2)).sum()
        sites, tries = make_sites(rays, G, r_ga, E_ga, X, lab0, size0, 5, np.random.default_rng(1))
        _G.update(rays=rays, X=X, G=G, r_ga=r_ga, sites=sites, dGA=dGA, sky=sky, lab0=lab0, size0=size0)
        for A in (0.0, 1.3, 3.0):
            t = time.time()
            print(one((0, A, "H", "u", 5)), f"{time.time() - t:.1f} s")
    else:
        run(*(int(a) for a in sys.argv[2:]))
