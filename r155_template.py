"""R-155 step 1 (method step 1): rebuild the Giant Arc from the Zhu & Menard DR7 + DR12 trimmed Mg II catalogues
with the discovery recipe (Lopez, Clowes & Williger 2022, arXiv:2201.06875): absorbers in z = 0.802 +- 0.060, single
linkage at 95 Mpc in comoving coordinates, flat Om = 0.27, h = 0.7. Freezes the members (data/r155/ga_members.npz,
plus a copy in r155_templates.json, which is committed) before any DESI file is opened.

DR7/DR12 deduplication: an absorber in both catalogues (sightlines within 1.5", |dz| < 0.005) is kept once (DR12).
The GA is the linkage group containing the absorber nearest (216.0, 39.0).

Run from the repo root: python -I research/cosmo-discovery/r155_template.py
"""
import json
import sys

import numpy as np
from astropy.io import fits
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

sys.path.insert(0, "research/cosmo-discovery")
from r155_common import D, DATA, xyz, sky_match  # noqa: E402

ZC, DZ, LINK = 0.802, 0.060, 95.0


def load_zm():
    rows = []
    for tag, f in (("DR7", "Trimmed_SDSS_DR7_107.fits"), ("DR12", "Trimmed_BOSS_DR12_107.fits")):
        d = fits.getdata(f"{DATA}/{f}", 1)
        rows.append(dict(ra=d["RA"].astype(float), dec=d["DEC"].astype(float), z=d["ZABS"].astype(float),
                         zq=d["ZQSO"].astype(float), w=d["REW_MGII_2796"].astype(float),
                         src=np.full(len(d), tag)))
    return rows


def main():
    dr7, dr12 = load_zm()
    # dedupe: DR7 absorbers already in DR12 (same sightline within 1.5", |dz| < 0.005)
    i7, i12 = sky_match(dr7["ra"], dr7["dec"], dr12["ra"], dr12["dec"], 1.5, many=True)
    dup = np.zeros(len(dr7["ra"]), bool)
    close = np.abs(dr7["z"][i7] - dr12["z"][i12]) < 0.005
    dup[i7[close]] = True
    cat = {k: np.concatenate([dr12[k], dr7[k][~dup]]) for k in dr7}
    print(f"Z&M: DR7 {len(dr7['ra'])}, DR12 {len(dr12['ra'])}, DR7 duplicates of DR12 {dup.sum()}, merged {len(cat['ra'])}")

    s = (cat["z"] > ZC - DZ) & (cat["z"] < ZC + DZ)
    sl = {k: v[s] for k, v in cat.items()}
    X = xyz(sl["ra"], sl["dec"], D(sl["z"]))
    pairs = cKDTree(X).query_pairs(LINK, output_type="ndarray")
    n = len(X)
    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    ncomp, lab = connected_components(g, directed=False)
    # GA group: the one containing the absorber nearest (216, 39) among groups with >= 20 members
    sizes = np.bincount(lab)
    c = np.cos(np.radians(sl["dec"])) * np.cos(np.radians(sl["ra"] - 216.0)) * np.cos(np.radians(39.0)) + \
        np.sin(np.radians(sl["dec"])) * np.sin(np.radians(39.0))
    order = np.argsort(-c)
    ga = next(lab[i] for i in order if sizes[lab[i]] >= 20)
    m = lab == ga
    print(f"slice {ZC}+-{DZ}: {n} absorbers, {ncomp} groups; largest sizes {sorted(sizes)[-6:]}")
    print(f"GA group: {m.sum()} members, RA {sl['ra'][m].min():.1f}-{sl['ra'][m].max():.1f}, "
          f"Dec {sl['dec'][m].min():.1f}-{sl['dec'][m].max():.1f}, z {sl['z'][m].min():.3f}-{sl['z'][m].max():.3f}, "
          f"W2796 median {np.median(sl['w'][m]):.2f} (10-90%: {np.percentile(sl['w'][m], 10):.2f}-"
          f"{np.percentile(sl['w'][m], 90):.2f})")
    np.savez(f"{DATA}/ga_members.npz", **{k: v[m] for k, v in sl.items()})
    np.savez(f"{DATA}/zm_slice.npz", **sl)
    np.savez(f"{DATA}/zm_merged.npz", **cat)
    tpl = json.load(open("research/cosmo-discovery/r155_templates.json")) if \
        __import__("os").path.exists("research/cosmo-discovery/r155_templates.json") else {}
    tpl["GA"] = dict(recipe=f"Z&M DR7+DR12 trimmed, z {ZC}+-{DZ}, single linkage {LINK} Mpc, Om 0.27 h 0.7",
                     members=[[round(float(a), 5), round(float(b), 5), round(float(z), 5), round(float(w), 3)]
                              for a, b, z, w in zip(sl["ra"][m], sl["dec"][m], sl["z"][m], sl["w"][m])])
    json.dump(tpl, open("research/cosmo-discovery/r155_templates.json", "w"), indent=0)


if __name__ == "__main__":
    main()
