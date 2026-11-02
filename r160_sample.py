"""R-160 samples: de-duplicated SPICE-RACS DR2 polarised components, matched to spectroscopic quasars and to JVAS/CLASS,
with a broadband spectral index.

Steps (all fixed in method-R160.md before any alignment statistic was computed):
1. goodRM rows (the catalogue's own subset: goodI, not leakage_flag, not snr_flag => polarised S/N >= 8).
2. De-duplicate overlapping tiles: rows within 10" are one source; keep the row nearest its tile centre.
3. Quasar match: nearest spectroscopic quasar (r160_qso.py) within 3".
4. Spectral index: least-squares slope of log S vs log nu over RACS-low2 (887.5 MHz, summed by source), RACS-mid
   (1367.5 MHz), NVSS (1400 MHz) and RACS-high (1655.5 MHz), each matched within 10" (nearest); needs >= 2 points
   spanning >= 1.3x in frequency. Flat: alpha >= -0.5 (S ~ nu^alpha), the JVAS/CLASS selection's limit.
5. Compact ("core-dominated" proxy): RACS-mid s_code 'S' (single Gaussian) and deconvolved major axis < 10".
6. JVAS/CLASS match within 5" (positions to mas in JVAS; ASKAP ~1").

Usage: python -I r160_sample.py data/r160 -> data/r160/sample.fits and prints the counts
"""
import glob
import sys

import numpy as np
from astropy.table import Table, vstack
from scipy.spatial import cKDTree

import r160_common as c


def dedup(t, r_arcsec=10.0):
    x = c.unit(t["ra"], t["dec"])
    order = np.argsort(np.asarray(t["separation_tile_centre"]))
    keep = np.ones(len(t), bool)
    tree = cKDTree(x)
    pairs = tree.query_pairs(c.chord(r_arcsec / 3600), output_type="ndarray")
    rank = np.empty(len(t), int)
    rank[order] = np.arange(len(t))
    # drop the worse-ranked member of every close pair (greedy by rank)
    worse = np.where(rank[pairs[:, 0]] > rank[pairs[:, 1]], pairs[:, 0], pairs[:, 1])
    better = np.where(rank[pairs[:, 0]] > rank[pairs[:, 1]], pairs[:, 1], pairs[:, 0])
    srt = np.argsort(rank[better])
    for b, w in zip(better[srt], worse[srt]):
        if keep[b]:
            keep[w] = False
    return t[keep]


def load_tap(d, name):
    parts = [Table.read(f, format="ascii.csv") for f in sorted(glob.glob(f"{d}/tap_{name}_*.csv")) if
             open(f).readline().startswith("ra")]
    parts = [p for p in parts if len(p)]
    return vstack(parts)


def load_nvss(d):
    ra, dec, s = [], [], []
    for l in open(f"{d}/nvss.tsv"):
        if not l[:1].isdigit():
            continue
        a = l.rstrip("\n").split("\t")
        if len(a) < 3 or not a[2].strip():
            continue
        h, m, sec = a[0].split()
        dd, dm, ds = a[1].split()
        sgn = -1 if dd.startswith("-") else 1
        ra.append(15 * (int(h) + int(m) / 60 + float(sec) / 3600))
        dec.append(sgn * (abs(int(dd)) + int(dm) / 60 + float(ds) / 3600))
        s.append(float(a[2]))
    return Table(dict(ra=np.array(ra), dec=np.array(dec), total_flux=np.array(s)))


def main(d):
    t = Table.read(f"{d}/spice_pol.fits")
    t = t[t["goodRM_flag"]]
    n0 = len(t)
    t = dedup(t)
    print(f"goodRM rows {n0} -> de-duplicated {len(t)}")
    q = Table.read(f"{d}/qso_all.fits")
    k, sep = c.xmatch(t["ra"], t["dec"], q["ra"], q["dec"], 3.0)
    t["z"] = np.where(k >= 0, q["z"][np.maximum(k, 0)], np.nan)
    t["z_src"] = np.where(k >= 0, q["src"][np.maximum(k, 0)], "")
    t["qso_sep"] = sep
    # spectral index from source-level catalogues: RACS-low2 (887.5 MHz, components summed by source_id),
    # RACS-mid (1367.5), NVSS (1400), RACS-high (1655.5). SPICE's own total_I_flux is a component flux (not used).
    nu, S = [], []
    low = load_tap(d, "racs_low2_v02")
    sid = np.asarray(low["source_id"]).astype(str)
    u_, inv = np.unique(sid, return_inverse=True)
    fl = np.asarray(low["total_flux"], float)
    tot = np.bincount(inv, weights=fl)
    rax = np.bincount(inv, weights=fl * np.asarray(low["ra"], float)) / tot  # fine away from ra = 0/360 wrap
    first = np.unique(inv, return_index=True)[1]
    ra0 = np.asarray(low["ra"], float)[first]
    rax = np.where(np.abs(rax - ra0) > 1, ra0, rax)
    decx = np.bincount(inv, weights=fl * np.asarray(low["dec"], float)) / tot
    kk, _ = c.xmatch(t["ra"], t["dec"], rax, decx, 10.0)
    S.append(np.where(kk >= 0, tot[np.maximum(kk, 0)], np.nan))
    nu.append(np.full(len(t), 887.5))
    for name, f0 in (("racs_mid_sources_v01", 1367.5), ("racs_high_sources_v01", 1655.5)):
        r = load_tap(d, name)
        kk, _ = c.xmatch(t["ra"], t["dec"], r["ra"], r["dec"], 10.0)
        S.append(np.where(kk >= 0, np.asarray(r["total_flux"], float)[np.maximum(kk, 0)], np.nan))
        nu.append(np.full(len(t), f0))
        if name.startswith("racs_mid"):
            t["mid_dcmaj"] = np.where(kk >= 0, np.asarray(r["dc_maj_axis"], float)[np.maximum(kk, 0)], np.nan)
            t["mid_scode"] = np.where(kk >= 0, np.asarray(r["s_code"]).astype(str)[np.maximum(kk, 0)], "")
    nv = load_nvss(d)
    kk, _ = c.xmatch(t["ra"], t["dec"], nv["ra"], nv["dec"], 10.0)
    S.append(np.where(kk >= 0, nv["total_flux"][np.maximum(kk, 0)], np.nan))
    nu.append(np.full(len(t), 1400.0))
    S, nu = np.array(S), np.array(nu)
    ok = np.isfinite(S) & (S > 0)
    lx, ly = np.where(ok, np.log10(nu), 0), np.where(ok, np.log10(np.where(ok, S, 1)), 0)
    n = ok.sum(0)
    mx = lx.sum(0) / np.maximum(n, 1)
    my = ly.sum(0) / np.maximum(n, 1)
    sxx = (np.where(ok, (lx - mx) ** 2, 0)).sum(0)
    sxy = (np.where(ok, (lx - mx) * (ly - my), 0)).sum(0)
    span = np.where(ok, nu, np.nan)
    span = np.nanmax(span, 0) / np.nanmin(span, 0)
    t["alpha"] = np.where((n >= 2) & (span >= 1.3) & (sxx > 0), sxy / np.maximum(sxx, 1e-12), np.nan)
    t["alpha_npts"] = n
    t["flat"] = t["alpha"] >= -0.5
    t["compact"] = (t["mid_scode"] == "S") & (t["mid_dcmaj"] < 10)
    j = Table.read(f"{d}/jvas_set.fits")  # P > 1 mJy set; the full JVAS list for the same-source check below
    kj, sj = c.xmatch(t["ra"], t["dec"], j["ra"], j["dec"], 5.0)
    t["jvas_idx"] = kj
    t["jvas_PA"] = np.where(kj >= 0, np.asarray(j["PA"], float)[np.maximum(kj, 0)], np.nan)
    t["jvas_ePA"] = np.where(kj >= 0, np.asarray(j["ePA"], float)[np.maximum(kj, 0)], np.nan)
    t["jvas_P"] = np.where(kj >= 0, np.asarray(j["P"], float)[np.maximum(kj, 0)], np.nan)
    t.write(f"{d}/sample.fits", overwrite=True)
    qm = np.isfinite(t["z"])
    print(f"alpha measured {np.isfinite(t['alpha']).sum()}; flat {t['flat'].sum()}; flat & compact {(t['flat'] & t['compact']).sum()}")
    print(f"quasar-matched {qm.sum()} (by source {dict(zip(*np.unique(t['z_src'][qm], return_counts=True)))})")
    print(f"  with alpha {np.isfinite(t['alpha'][qm]).sum()}, flat {(t['flat'] & qm).sum()}, flat & compact {(t['flat'] & t['compact'] & qm).sum()}")
    print(f"JVAS (P>1 mJy) matched {np.sum(kj >= 0)}; of them quasar-matched {np.sum((kj >= 0) & qm)}, flat {np.sum((kj >= 0) & t['flat'])}")


if __name__ == "__main__":
    main(sys.argv[1])
