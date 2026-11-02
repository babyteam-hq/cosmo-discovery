"""R-155 measurement, revised after the method check (results-R155.md F4).

  python -I research/cosmo-discovery/r155_measure.py predict   # real-structure prediction per cell, sigma, primary cell;
                                                               # no observed count in any GA cell is computed
  python -I research/cosmo-discovery/r155_measure.py measure   # the measurement

Cells: R in {47.5, 70, 95} Mpc spheres and the CHMS hull, x W2796 > 0.3, 0.6 A (Z&M scale). Close-pair exclusion: DESI
sightlines within KPC proper kpc of a member's own sightline (at the member's z) are dropped, in the GA and in every null
(primary 200 kpc; 0 and 500 as variants).

Calibration: every template (GA, each comparison group) is referred to its own shift null (the same shape rotated to
random places in the SDSS footprint, N_exp within 0.75-1.33x its own); delta_cal = delta - mean(shift delta). This removes
the completeness residual of the SDSS-thinned sightlines (feasibility: -0.05 to -0.07) the same way for all templates.

Chance baseline: the other single-linkage groups of the same slice and recipe, matched to the GA in size (23-92 members)
and extent (largest member separation 0.5-2x the GA's). Baseline = mean of their delta_cal; claim per cell
Delta = delta_cal(GA) - baseline, interval from the GA's shift-null percentiles with the baseline's standard error added in
quadrature. The random-46 clustering null of the proposal is kept as a secondary row.

Real-structure prediction ("injection" computed in expectation): an overdensity A uniform in the CHMS cut to the slice
z 0.742-0.862 (profile H; the published A = 1.3) raises a cell's expected count by A x f_H, f_H = e(cell inside CHMS and
slice) / e(cell). Profile S puts the same excess (A V_H / V_S) only in member spheres of radius r_sph cut to the slice.
The fitted amplitude is A = Delta / f_H, the overdensity of a real structure in the CHMS: 0 for a chance pattern, 1.3 as
published, in between for a real structure weaker than its selected members imply.
"""
import json
import sys

import numpy as np
from scipy.spatial import ConvexHull, cKDTree

sys.path.insert(0, "research/cosmo-discovery")
from r155_common import DATA, D, xyz  # noqa: E402
from r155_core import Sky, clip_D, drop, fib_sphere, intersect, mst_mean_edge, rotation_to  # noqa: E402

RS = (47.5, 70.0, 95.0)
CUTS = (0.3, 0.6)
KPC = 200.0
SLICE = (D(0.742), D(0.862))
A_PUB = 1.3
rng = np.random.default_rng(1550)


def members(name="GA"):
    m = np.array(json.load(open("research/cosmo-discovery/r155_templates.json"))[name]["members"])
    return xyz(m[:, 0], m[:, 1], D(m[:, 2]))


def cell_ivs(sky, X, kpc, r_sph):
    bad = sky.close_sightlines(X, kpc)
    out = {f"R{R}": drop(sky.sphere_intervals(X, R), bad) for R in RS}
    si, a1, a2, _ = sky.hull_intervals(X, r_sph)
    out["hull"] = drop((si, a1, a2), bad)
    return out


def deltas(sky, ivs, with_obs=True):
    res = {}
    for cell, iv in ivs.items():
        for c in CUTS:
            k = sky.counts(iv, c)
            res[f"{cell}_W{c}"] = k if with_obs else {x: v for x, v in k.items() if not x.startswith("nobs")}
    return res


def shift_null(sky, X, ref, n_keep, kpc, max_try=None):
    """Shift deltas per cell for template X; ref holds the template's own N_exp per cell key."""
    r_sph = 0.5 * mst_mean_edge(X)
    cen = X.mean(0)
    cu = cen / np.linalg.norm(cen)
    r_h = np.max(np.linalg.norm(X - cen, axis=1))
    res = {k: [] for k in ref}
    tries = 0
    max_try = max_try or 20 * n_keep
    while min(len(v) for v in res.values()) < n_keep and tries < max_try:
        tries += 1
        p = rng.integers(0, len(sky.u))
        Y = X @ rotation_to(cu, sky.u[p], rng.uniform(0, 2 * np.pi)).T
        if min(sky.sdss_share(y / np.linalg.norm(y)) for y in Y[::max(1, len(Y) // 10)]) < 0.3:
            continue
        dmin = np.min(np.linalg.norm(Y[:, None] - X[None], axis=2))
        bad = sky.close_sightlines(Y, kpc)
        for cell in [f"R{R}" for R in RS] + ["hull"]:
            if cell == "hull":
                if np.linalg.norm(Y.mean(0) - cen) < 2 * r_h:
                    continue
                si, a1, a2, _ = sky.hull_intervals(Y, r_sph)
                iv = drop((si, a1, a2), bad)
            else:
                R = float(cell[1:])
                if dmin < 2 * R:
                    continue
                iv = drop(sky.sphere_intervals(Y, R), bad)
            for c in CUTS:
                key = f"{cell}_W{c}"
                if len(res[key]) >= 4 * n_keep:
                    continue
                k = sky.counts(iv, c)
                if ref[key] > 0 and 0.75 < k["nexp"] / ref[key] < 1.33:
                    res[key].append([k["nobs"], k["nexp"], k["nobs_w"], k["nexp_w"]])
    return res


def summ(a):
    a = np.array(a)
    if len(a) < 2:
        return dict(n=len(a))
    d = a[:, 0] / a[:, 1] - 1
    return dict(n=len(a), mean=float(d.mean()), sd=float(d.std(ddof=1)), p16=float(np.percentile(d, 16)),
                p84=float(np.percentile(d, 84)), p2=float(np.percentile(d, 2.5)), p97=float(np.percentile(d, 97.5)),
                mean_w=float((a[:, 2] / a[:, 3] - 1).mean()))


def regions(sky, X, r_sph):
    """CHMS-and-slice and spheres(r_sph)-and-slice: sightline intervals and volumes (Monte Carlo)."""
    si, a1, a2, _ = sky.hull_intervals(X, r_sph)
    H = clip_D((si, a1, a2), *SLICE)
    S = clip_D(sky.sphere_intervals(X, r_sph), *SLICE)
    pts = (X[:, None, :] + fib_sphere(200)[None] * r_sph).reshape(-1, 3)
    hull = ConvexHull(pts)
    lo, hi = pts.min(0), pts.max(0)
    nh = ns = 0
    n = 0
    for _ in range(40):  # 40 x 20k points, chunked (the hull has thousands of facets)
        P = rng.uniform(lo, hi, (20_000, 3))
        r = np.linalg.norm(P, axis=1)
        sl = (r > SLICE[0]) & (r < SLICE[1])
        nh += np.count_nonzero(np.all(P @ hull.equations[:, :3].T + hull.equations[:, 3] <= 0, axis=1) & sl)
        ns += np.count_nonzero((cKDTree(X).query(P, distance_upper_bound=r_sph)[0] < np.inf) & sl)
        n += len(P)
    box = np.prod(hi - lo)
    return H, S, nh / n * box, ns / n * box


def matched_groups():
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial.distance import pdist
    sl = np.load(f"{DATA}/zm_slice.npz")
    P = xyz(sl["ra"], sl["dec"], D(sl["z"]))
    pr = cKDTree(P).query_pairs(95.0, output_type="ndarray")
    g = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(len(P), len(P)))
    _, lab = connected_components(g, directed=False)
    sizes = np.bincount(lab)
    ga = members()
    ext_ga = pdist(ga).max()
    out = []
    for gi in np.nonzero((sizes >= 23) & (sizes <= 92))[0]:
        X = P[lab == gi]
        if np.min(np.linalg.norm(X[:, None] - ga[None], axis=2)) < 1.0:
            continue
        e = pdist(X).max()
        if 0.5 * ext_ga <= e <= 2 * ext_ga:
            out.append((int(gi), X, float(e)))
    return out, float(ext_ga)


def predict():
    sky = Sky(CUTS)
    X = members()
    r_sph = 0.5 * mst_mean_edge(X)
    ivs = cell_ivs(sky, X, KPC, r_sph)
    H, S, VH, VS = regions(sky, X, r_sph)
    print(f"GA: r_sph {r_sph:.1f} Mpc; V(CHMS & slice) {VH:.3e} Mpc^3, V(spheres r_sph & slice) {VS:.3e}; "
          f"profile-S amplitude {A_PUB * VH / VS:.2f}")
    ref = deltas(sky, ivs, with_obs=False)
    out = {"r_sph": r_sph, "VH": VH, "VS": VS, "cells": {}}
    for cell, iv in ivs.items():
        for c in CUTS:
            key = f"{cell}_W{c}"
            e = ref[key]["nexp"]
            fH = sky.counts(intersect(iv, H), c)["nexp"] / e
            fS = sky.counts(intersect(iv, S), c)["nexp"] / e
            out["cells"][key] = dict(nexp=e, nsl=ref[key]["nsl"], fH=fH, fS=fS, predH=A_PUB * fH,
                                     predS=A_PUB * VH / VS * fS)
    sh = shift_null(sky, X, {k: v["nexp"] for k, v in out["cells"].items()}, 300, KPC)
    for key, v in out["cells"].items():
        s = summ(sh[key])
        v["shift"] = s
        v["sepH"] = v["predH"] / s["sd"]
        v["sepS"] = v["predS"] / s["sd"]
        print(f"{key}: sightlines {v['nsl']}, N_exp {v['nexp']:.2f}, f_H {v['fH']:.3f} -> pred {v['predH']:.2f}; "
              f"f_S {v['fS']:.3f} -> pred {v['predS']:.2f}; shift n {s['n']} mean {s['mean']:+.3f} sd {s['sd']:.3f}; "
              f"separation H {v['sepH']:.1f} sigma, S {v['sepS']:.1f} sigma")
    groups, ext = matched_groups()
    print(f"matched comparison groups (23-92 members, extent {0.5 * ext:.0f}-{2 * ext:.0f} Mpc): {len(groups)}, "
          f"sizes {[len(g[1]) for g in groups]}")
    best = max((k for k in out["cells"] if k.endswith("W0.3")), key=lambda k: out["cells"][k]["sepH"])
    out["primary"] = best
    print(f"primary cell (largest profile-H separation, W > 0.3): {best}")
    json.dump(out, open(f"{DATA}/predict.json", "w"), indent=1, default=float)


def measure():
    sky = Sky(CUTS)
    pred = json.load(open(f"{DATA}/predict.json"))
    X = members()
    r_sph = pred["r_sph"]
    out = {"primary": pred["primary"], "ga": {}, "ga_shift": {}, "groups": [], "baseline": {}, "claim": {}}
    for kpc in (0.0, KPC, 500.0):
        res = deltas(sky, cell_ivs(sky, X, kpc, r_sph))
        out["ga"][str(kpc)] = res
    ga = out["ga"][str(KPC)]
    sh = shift_null(sky, X, {k: v["nexp"] for k, v in ga.items()}, 1000, KPC, max_try=6000)
    out["ga_shift"] = {k: summ(v) for k, v in sh.items()}
    out["ga_shift_raw"] = sh
    groups, _ = matched_groups()
    for gi, Xg, ext in groups:
        rg = 0.5 * mst_mean_edge(Xg)
        r = deltas(sky, cell_ivs(sky, Xg, KPC, rg))
        if r["R95.0_W0.3"]["nexp"] < 3:
            continue
        s = shift_null(sky, Xg, {k: v["nexp"] for k, v in r.items()}, 60, KPC, max_try=1500)
        row = dict(group=gi, size=len(Xg), extent=ext, cells={})
        for k, v in r.items():
            ss = summ(s[k])
            d = v["nobs"] / v["nexp"] - 1 if v["nexp"] > 0 else np.nan
            row["cells"][k] = dict(nobs=v["nobs"], nexp=v["nexp"], delta=d, off=ss.get("mean", np.nan), nshift=ss["n"])
        out["groups"].append(row)
        print(f"group {gi} ({len(Xg)} members): R95 W0.3 N_obs {r['R95.0_W0.3']['nobs']:.0f} N_exp "
              f"{r['R95.0_W0.3']['nexp']:.1f}", flush=True)
    for k, v in ga.items():
        dg = v["nobs"] / v["nexp"] - 1
        s = out["ga_shift"][k]
        cal = np.array([g["cells"][k]["delta"] - g["cells"][k]["off"] for g in out["groups"]
                        if g["cells"][k]["nshift"] >= 20 and g["cells"][k]["nexp"] >= 2])
        base, se = float(cal.mean()), float(cal.std(ddof=1) / np.sqrt(len(cal)))
        out["baseline"][k] = dict(n=len(cal), mean=base, se=se, sd=float(cal.std(ddof=1)))
        Dl = dg - s["mean"] - base
        lo68 = np.hypot(s["mean"] - s["p16"], se)
        hi68 = np.hypot(s["p84"] - s["mean"], se)
        lo95 = np.hypot(s["mean"] - s["p2"], se)
        hi95 = np.hypot(s["p97"] - s["mean"], se)
        fH = pred["cells"][k]["fH"]
        out["claim"][k] = dict(delta=dg, delta_w=v["nobs_w"] / v["nexp_w"] - 1, off=s["mean"], base=base,
                               Delta=Dl, lo68=Dl - lo68, hi68=Dl + hi68, lo95=Dl - lo95, hi95=Dl + hi95,
                               z0=Dl / ((lo68 + hi68) / 2), A=Dl / fH, A_lo68=(Dl - lo68) / fH,
                               A_hi68=(Dl + hi68) / fH, A_lo95=(Dl - lo95) / fH, A_hi95=(Dl + hi95) / fH,
                               z_pub=(pred["cells"][k]["predH"] - Dl) / ((lo68 + hi68) / 2),
                               rank=float(np.mean(cal < dg - s["mean"])))
        c = out["claim"][k]
        print(f"{k}: N_obs {v['nobs']:.0f} N_exp {v['nexp']:.2f} delta {dg:+.3f} (TSNR2-w {c['delta_w']:+.3f}); "
              f"shift off {s['mean']:+.3f} sd {s['sd']:.3f}; baseline {base:+.3f}+-{se:.3f} (n {len(cal)}); "
              f"Delta {Dl:+.3f} [{c['lo68']:+.2f},{c['hi68']:+.2f}] 95% [{c['lo95']:+.2f},{c['hi95']:+.2f}]; "
              f"A {c['A']:+.2f} [{c['A_lo95']:+.2f},{c['A_hi95']:+.2f}]; vs 0 {c['z0']:.1f} sig; "
              f"vs published {c['z_pub']:.1f} sig; GA rank among groups {c['rank']:.2f}")
    for kpc in ("0.0", "500.0"):
        r = out["ga"][kpc]
        print(f"kpc {kpc}: " + "; ".join(f"{k} {v['nobs']:.0f}/{v['nexp']:.2f}" for k, v in r.items()))
    json.dump(out, open(f"{DATA}/measure.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    {"predict": predict, "measure": measure}[sys.argv[1]]()
