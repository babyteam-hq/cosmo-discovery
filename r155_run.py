"""R-155 driver.

  python -I research/cosmo-discovery/r155_run.py feasibility   # step 0: expected counts and null scatter; the template's
                                                               # observed count is neither computed nor printed
  python -I research/cosmo-discovery/r155_run.py measure       # the fixed measurement (after the method check)

Templates: GA (r155_templates.json, frozen by r155_template.py). Statistic per template, R in {47.5, 70, 95} Mpc and
the CHMS hull, EW cut in {0.3, 0.6} A: delta = N_obs / N_exp - 1 (r155_core). Nulls:
  shift  - the template rotated rigidly on the sky (random centre drawn from independent DESI sightlines, random spin;
           redshifts kept), kept only inside the SDSS footprint (every 5th member position has > 30 % SDSS quasars
           among DESI quasars within 1 deg, so the SDSS exclusion thins the sightlines as it does on the GA), dropped if any sphere centre lies within 2R of a GA sphere centre (or the shifted hull
           centre within 2 x the hull radius of the GA's); kept if its N_exp is within x0.75-1.33 of the template's
           (so its scatter is the template's own sigma at the same depth).
  clust  - spheres on random sets of 46 Z&M slice absorbers (GA members excluded) drawn from those with independent
           DESI coverage within R: the ordinary-clustering baseline.
  groups - every other single-linkage group (same recipe, >= 20 members) with N_exp >= 3: what the algorithm's
           groups in general show on new sightlines.
Writes data/r155/<mode>.json.
"""
import json
import sys

import numpy as np

sys.path.insert(0, "research/cosmo-discovery")
from r155_common import DATA, D, xyz  # noqa: E402
from r155_core import Sky, mst_mean_edge, rotation_to  # noqa: E402

RS = (47.5, 70.0, 95.0)
CUTS = (0.3, 0.6)
rng = np.random.default_rng(155)


def members(name="GA"):
    m = np.array(json.load(open("research/cosmo-discovery/r155_templates.json"))[name]["members"])
    return xyz(m[:, 0], m[:, 1], D(m[:, 2]))


def delta(c):
    return c["nobs"] / c["nexp"] - 1 if c["nexp"] > 0 else np.nan


def run_template(sky, X, with_obs):
    out = {}
    for R in RS:
        iv = sky.sphere_intervals(X, R)
        for c in CUTS:
            k = sky.counts(iv, c)
            if not with_obs:
                k = {x: v for x, v in k.items() if not x.startswith("nobs")}
            out[f"R{R}_W{c}"] = k
    r_sph = 0.5 * mst_mean_edge(X)
    si, d1, d2, vol = sky.hull_intervals(X, r_sph)
    for c in CUTS:
        k = sky.counts((si, d1, d2), c)
        if not with_obs:
            k = {x: v for x, v in k.items() if not x.startswith("nobs")}
        k.update(r_sph=r_sph, volume=vol)
        out[f"hull_W{c}"] = k
    return out


def shift_null(sky, X, ref, n_try):
    cen = X.mean(0)
    cu = cen / np.linalg.norm(cen)
    r_hull = np.max(np.linalg.norm(X - cen, axis=1))
    res = {key: [] for key in ref}
    done = 0
    while done < n_try:
        p = rng.integers(0, len(sky.u))
        Rm = rotation_to(cu, sky.u[p], rng.uniform(0, 2 * np.pi))
        Y = X @ Rm.T
        # inside the SDSS footprint, where the SDSS-sightline exclusion bites as it does on the GA
        if min(sky.sdss_share(y / np.linalg.norm(y)) for y in Y[::5]) < 0.3:
            continue
        done += 1
        for R in RS:
            dmin = np.min(np.linalg.norm(Y[:, None] - X[None], axis=2))
            if dmin < 2 * R:
                continue
            iv = sky.sphere_intervals(Y, R)
            for c in CUTS:
                key = f"R{R}_W{c}"
                k = sky.counts(iv, c)
                if 0.75 < k["nexp"] / ref[key]["nexp"] < 1.33:
                    res[key].append([k["nobs"], k["nexp"], k["nobs_w"], k["nexp_w"]])
        if np.linalg.norm(Y.mean(0) - cen) > 2 * r_hull:
            si, d1, d2, _ = sky.hull_intervals(Y, ref["hull_W0.3"]["r_sph"])
            for c in CUTS:
                key = f"hull_W{c}"
                k = sky.counts((si, d1, d2), c)
                if 0.75 < k["nexp"] / ref[key]["nexp"] < 1.33:
                    res[key].append([k["nobs"], k["nexp"], k["nobs_w"], k["nexp_w"]])
    return res


def clust_null(sky, nsets, nmem=46):
    sl = np.load(f"{DATA}/zm_slice.npz")
    ga = np.array(json.load(open("research/cosmo-discovery/r155_templates.json"))["GA"]["members"])
    isga = np.zeros(len(sl["ra"]), bool)
    for a, b, z, _ in ga:
        isga |= (np.abs(sl["ra"] - a) < 1e-4) & (np.abs(sl["dec"] - b) < 1e-4) & (np.abs(sl["z"] - z) < 1e-4)
    P = xyz(sl["ra"], sl["dec"], D(sl["z"]))[~isga]
    res = {}
    for R in RS:
        # per-absorber coverage: sightlines through its sphere
        cov = np.array([len(sky.sphere_intervals(P[i:i + 1], R)[0]) for i in range(len(P))])
        pool = np.nonzero(cov > 0)[0]
        for c in CUTS:
            res[f"R{R}_W{c}"] = []
        for _ in range(nsets):
            pick = rng.choice(pool, nmem, replace=False)
            iv = sky.sphere_intervals(P[pick], R)
            for c in CUTS:
                k = sky.counts(iv, c)
                res[f"R{R}_W{c}"].append([k["nobs"], k["nexp"], k["nobs_w"], k["nexp_w"]])
        res[f"pool_R{R}"] = int(len(pool))
    return res, int(isga.sum())


def groups_null(sky):
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree
    sl = np.load(f"{DATA}/zm_slice.npz")
    P = xyz(sl["ra"], sl["dec"], D(sl["z"]))
    pr = cKDTree(P).query_pairs(95.0, output_type="ndarray")
    g = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(len(P), len(P)))
    _, lab = connected_components(g, directed=False)
    sizes = np.bincount(lab)
    ga = members()
    out = []
    for gi in np.nonzero(sizes >= 20)[0]:
        X = P[lab == gi]
        if np.min(np.linalg.norm(X[:, None] - ga[None], axis=2)) < 1.0:
            continue
        row = {"size": int(sizes[gi])}
        for R in RS:
            iv = sky.sphere_intervals(X, R)
            for c in CUTS:
                k = sky.counts(iv, c)
                row[f"R{R}_W{c}"] = [k["nobs"], k["nexp"]]
        out.append(row)
    return out


def summarize(arr):
    a = np.array(arr)
    if len(a) == 0:
        return {}
    d = a[:, 0] / a[:, 1] - 1
    dw = a[:, 2] / a[:, 3] - 1
    return dict(n=len(a), mean=float(d.mean()), sd=float(d.std(ddof=1)), p16=float(np.percentile(d, 16)),
                p84=float(np.percentile(d, 84)), p2=float(np.percentile(d, 2.5)), p97=float(np.percentile(d, 97.5)),
                ratio_of_sums=float(a[:, 0].sum() / a[:, 1].sum() - 1), median_nexp=float(np.median(a[:, 1])),
                mean_w=float(dw.mean()), sd_w=float(dw.std(ddof=1)))


def main(mode):
    sky = Sky(CUTS)
    print(f"independent sightlines {len(sky.zq)}; all-sky dN/dz at z 0.742-0.862 on them:")
    for c in CUTS:
        sl, z = sky.abs[c]
        m = (z > 0.742) & (z < 0.862)
        p = np.clip(np.minimum(sky.zmax, 0.862) - np.maximum(sky.zmin, 0.742), 0, None).sum()
        print(f"  W > {c}: {m.sum()} absorbers / dz {p:.0f} = {m.sum() / p:.4f}")
    X = members()
    with_obs = mode == "measure"
    ref = run_template(sky, X, with_obs)
    for key, k in ref.items():
        extra = f" N_obs {k['nobs']:.0f} delta {delta(k):+.3f} (TSNR2-weighted {k['nobs_w'] / k['nexp_w'] - 1:+.3f})" \
            if with_obs else ""
        print(f"GA {key}: sightlines {k['nsl']}, path dz {k['dz']:.1f}, N_exp {k['nexp']:.2f}{extra}")
    nshift = 600 if mode == "feasibility" else 3000
    nclust = 200 if mode == "feasibility" else 1000
    sh = shift_null(sky, X, ref, nshift)
    cl, nga = clust_null(sky, nclust)
    out = dict(mode=mode, template=ref, shift={k: summarize(v) for k, v in sh.items()},
               clust={k: (summarize(v) if isinstance(v, list) else v) for k, v in cl.items()}, ga_matched=nga)
    if mode == "measure":
        out["groups"] = groups_null(sky)
        out["shift_raw"] = sh
        out["clust_raw"] = {k: v for k, v in cl.items() if isinstance(v, list)}
    for key in ref:
        s, c = out["shift"].get(key, {}), out["clust"].get(key, {})
        line = f"{key}: shift n={s.get('n')} mean {s.get('mean', np.nan):+.3f} sd {s.get('sd', np.nan):.3f}"
        if c:
            line += f" | clust n={c['n']} mean {c['mean']:+.3f} (ratio of sums {c['ratio_of_sums']:+.3f}) sd {c['sd']:.3f}" \
                    f" median N_exp {c['median_nexp']:.1f}"
        print(line)
    json.dump(out, open(f"{DATA}/{mode}.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main(sys.argv[1])
