"""R-160 anchor: reproduce Tiwari & Jain 2013's JVAS/CLASS S_D at n_v = 10 and express it as a pairwise amplitude C.

Set 1 of Tiwari & Jain (arXiv:1201.5180): polarised flux sqrt(Q^2+U^2) > 1 mJy (they quote 4,400 sources).
Neighbours in 2D (angle); S_D maximised with parallel transport, site included; null: shuffled angles (their
Procedure I). Then the pair-averaged C over the same neighbour pairs, and over all pairs within fixed angles.

Usage: python -I r160_jvas.py data/r160/jvas_pol.tsv  -> prints the anchor; writes data/r160/jvas_set.fits
"""
import sys

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u
from scipy.spatial import cKDTree

import r160_common as c


def load_jvas(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path) if l and l[0] not in "#\n-\"R"]
    rows = [r for r in rows if len(r) >= 9 and r[0].strip()]
    ra = [r[0] for r in rows]
    dec = [r[1] for r in rows]
    sc = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

    def f(k):
        return np.array([float(r[k]) if r[k].strip() else np.nan for r in rows])

    t = Table(dict(ra=sc.ra.deg, dec=sc.dec.deg, I=f(2), eI=f(3), Q=f(4), eQ=f(5), U=f(6), eU=f(7), PA=f(8)))
    t["P"] = np.hypot(t["Q"], t["U"])
    t["ePA"] = np.degrees(0.5 * np.hypot(t["Q"] * t["eU"], t["U"] * t["eQ"]) / t["P"] ** 2)
    return t


def knn_pairs(x, nv):
    tree = cKDTree(x)
    d, k = tree.query(x, k=nv + 1)
    j = np.repeat(np.arange(len(x)), nv)  # site
    i = k[:, 1:].ravel()  # neighbour
    return i, j, d[:, 1:]


def main(path):
    t = load_jvas(path)
    print(f"JVAS/CLASS rows {len(t)}")
    s1 = t[(t["P"] > 1) & np.isfinite(t["PA"])]
    print(f"set 1 (P > 1 mJy): {len(s1)} (Tiwari & Jain: 4400)")
    s1.write("data/r160/jvas_set.fits", overwrite=True)
    x = c.unit(s1["ra"], s1["dec"])
    psi = np.radians(np.asarray(s1["PA"], float))
    rng = np.random.default_rng(160)
    n = len(s1)
    for nv in (6, 10, 15):
        i, j, d = knn_pairs(x, nv - 1)  # sum includes the site: n_v members = site + (n_v - 1) neighbours
        g = c.transported_angle(x[i], x[j], np.zeros(len(i)))
        sd = c.sd_max_from_pairs(n, i, j, g, psi)
        null = np.array([c.sd_max_from_pairs(n, i, j, g, rng.permutation(psi)) for _ in range(1000)])
        cp = np.mean(np.cos(2 * (psi[i] + g - psi[j])))
        nullc = np.array([np.mean(np.cos(2 * (p[i] + g - p[j]))) for p in (rng.permutation(psi) for _ in range(1000))])
        th = np.degrees(2 * np.arcsin(d[:, -1] / 2))
        print(f"n_v={nv:2d}: S_D={sd:.4f} null {null.mean():.4f}+-{null.std():.5f} -> {(sd-null.mean())/null.std():+.2f} sig;"
              f"  pairwise C={cp:+.4f} null {nullc.mean():+.4f}+-{nullc.std():.4f} -> {(cp-nullc.mean())/nullc.std():+.2f} sig;"
              f" median angle to n_v-th member {np.median(th):.2f} deg")


if __name__ == "__main__":
    main(sys.argv[1])
