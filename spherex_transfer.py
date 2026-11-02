"""R-143: the full-chain transfer of R-138's SPHEREx x DESI estimator, measured by injection, and the corrected table.

Stages (outputs under data/spherex/r143/):
  maps   : accumulate spherex_inject.py outputs into NSIDE-2048 maps per channel group (R-138's 18) for every kept array
           (real chain, real chain with the star mask, recovered mocks A-E, the mocks on good pixels, the truths).
  pairs  : per NSIDE-8 patch, R-138's gnomonic 1.72' grids and ring convolutions (spherex_xcorr.py), centres: galaxies,
           randoms and Gaia 12 < G < 16 stars per slice; fields: every map above and the galaxy/random count maps.
  result : per slice and group, w̄ of each field pair; transfer T = w̄[recovered] / w̄[truth]; the real signal on this
           subset (b_r and w̄_m as R-138) corrected by T; nulls; star split; jackknife over patches.
  report : the tables for results-R138.md.
Usage: python -u spherex_transfer.py maps|pairs [nproc]|result|report
"""
import glob
import json
import os
import sys
import time

import numpy as np
import healpy as hp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from spherex_xcorr import (SLICES, NS, NBIN, PIXDEG, NSIDE, GROUPS, TRACERS, bin_edges, gnomonic, load_points,  # noqa
                           region, in_region, theory_wm, wbar_weights, group_edges, CIBER)
from spherex_inject import KEYS, load_stars as load_stars_idx  # noqa

D = os.path.join(HERE, "data")
SP = os.path.join(D, "spherex")
OUT = os.environ.get("R143_OUT", os.path.join(SP, "r143"))
NG = len(GROUPS)
GROUP_OF = np.full(128, -1, int)
for gi, chans in enumerate(GROUPS):
    GROUP_OF[chans] = gi


def stars_all():
    s = load_stars_idx()
    return dict(ra=s["ra"], dec=s["dec"], g=s["g"])


# ---------------------------------------------------------------- maps
def stage_maps():
    os.makedirs(OUT, exist_ok=True)
    files = sorted(f for f in glob.glob(os.path.join(SP, "inj", "*.npz")) if ".tmp." not in f)
    if os.environ.get("R143_MAXIMG"):
        files = files[: int(os.environ["R143_MAXIMG"])]
    print("images", len(files), flush=True)
    upix = np.unique(np.concatenate([np.unique(np.load(f)["hp"]) for f in files]))
    npx = len(upix)
    print("map pixels", npx, flush=True)
    lut = np.full(hp.nside2npix(NSIDE), -1, np.int64)
    lut[upix] = np.arange(npx)
    M = np.lib.format.open_memmap(os.path.join(OUT, "maps.npy"), mode="w+", dtype=np.float32,
                                  shape=(len(KEYS), NG, npx))
    acc = M.reshape(len(KEYS), NG * npx)  # file-backed: the maps do not fit in memory beside the reducers
    dets = []
    fsm = []
    t0 = time.time()
    for k, f in enumerate(files):
        d = np.load(f)
        g = GROUP_OF[d["ch"].astype(int)]
        ok = g >= 0
        flat = g[ok] * npx + lut[d["hp"][ok]]
        uf, inv = np.unique(flat, return_inverse=True)
        for i, key in enumerate(KEYS):
            acc[i, uf] += np.bincount(inv, weights=d[key][ok], minlength=len(uf)).astype(np.float32)
        dets.append(int(d["det"]))
        fsm.append(float(d["frac_sm"]))
        if k % 2000 == 0:
            print(k, f"{time.time()-t0:.0f}s", flush=True)
    M.flush()
    np.savez(os.path.join(OUT, "maps_meta.npz"), upix=upix, nimg=len(files), dets=np.array(dets), frac_sm=np.array(fsm))
    print("maps saved; images per detector", np.bincount(dets), "star-mask block fraction", np.mean(fsm), flush=True)


# ---------------------------------------------------------------- pairs
def _setup():
    meta = np.load(os.path.join(OUT, "maps_meta.npz"))
    upix = meta["upix"]
    stars = stars_all()
    # R-138's map-level Gaia mask (amended): G < 10, r = 2' x 10^(0.2 (10 - G))
    b = stars["g"] < 10
    rad = np.maximum(2.0 / 60 * 10 ** (0.2 * (10 - stars["g"][b])), PIXDEG)
    vec = hp.ang2vec(stars["ra"][b], stars["dec"][b], lonlat=True)
    bad = [hp.query_disc(NSIDE, v, np.radians(rr), nest=True, inclusive=True) for v, rr in zip(vec, rad)]
    keep = ~np.isin(upix, np.unique(np.concatenate(bad)))
    return upix, stars, keep


def pairs_patch(args):
    ip, p8 = args
    out = os.path.join(OUT, "pairs", f"p{p8}.npz")
    if os.path.exists(out):
        return
    t0 = time.time()
    G = _G
    upix, stars, keep, pts, M = G["upix"], G["stars"], G["keep"], G["pts"], G["M"]
    edges, half, ng, nfft, kern_f = G["edges"], G["half"], G["ng"], G["nfft"], G["kern_f"]
    ra0, dec0 = hp.pix2ang(8, p8, lonlat=True)

    def grid_index(ra, dec):
        x, y, cc = gnomonic(ra, dec, ra0, dec0)
        i = np.floor((x + half) / PIXDEG).astype(int)
        j = np.floor((y + half) / PIXDEG).astype(int)
        ok = (cc > 0) & (i >= 0) & (i < ng) & (j >= 0) & (j < ng)
        return i * ng + j, ok

    # only points near this patch can land on its grid (the grid's corner is ~1.42 half from the centre)
    vec0 = hp.ang2vec(ra0, dec0, lonlat=True)
    rad = np.radians(half * 1.45 + 4.5)
    near_ring = hp.query_disc(8, vec0, rad, inclusive=True)
    near_nest = hp.query_disc(8, vec0, rad, inclusive=True, nest=True)
    cand = np.flatnonzero(np.isin(upix >> 16, near_nest) & keep)
    gi, ok = grid_index(G["ra_p"][cand], G["dec_p"][cand])
    sel = cand[ok]
    gsel = gi[ok]
    sub = {(t, key): np.flatnonzero(np.isin(pts[t][key]["p8"], near_ring)) for t in TRACERS for key in ("g", "r")}
    nf = len(KEYS) * NG + 2 * NS
    F = np.zeros((nf, ng * ng), np.float32)
    r = 0
    for i in range(len(KEYS)):
        for g in range(NG):
            F[r] = np.bincount(gsel, weights=M[i, g][sel], minlength=ng * ng); r += 1
    for key in ("g", "r"):
        for s, (t, lo, hi) in enumerate(SLICES):
            cat = pts[t][key]
            m = sub[(t, key)]
            m = m[(cat["z"][m] >= lo) & (cat["z"][m] < hi)]
            ii, okk = grid_index(cat["ra"][m], cat["dec"][m])
            F[r] = np.bincount(ii[okk], weights=cat["w"][m][okk], minlength=ng * ng); r += 1
    FT = F.T.astype(np.float64)
    del F
    sets = []
    sc = (stars["g"] >= 12) & (stars["g"] < 16) & G["star_core"] & (G["star_p8"] == p8)
    si, sok = grid_index(stars["ra"][sc], stars["dec"][sc])
    for s, (t, lo, hi) in enumerate(SLICES):
        for key in ("g", "r"):
            cat = pts[t][key]
            m = sub[(t, key)]
            m = m[(cat["z"][m] >= lo) & (cat["z"][m] < hi) & cat["core"][m] & (cat["p8"][m] == p8)]
            ii, okk = grid_index(cat["ra"][m], cat["dec"][m])
            sets.append((s, "gal" if key == "g" else "ran", ii[okk], cat["w"][m][okk]))
        sets.append((s, "star", si[sok], np.ones(sok.sum())))
    res = np.zeros((len(sets), NBIN, nf), np.float64)
    wsum = np.zeros(len(sets))
    for k, (s, key, ii, ww) in enumerate(sets):
        wsum[k] = ww.sum()
        if len(ii) == 0:
            continue
        P = np.zeros((nfft, nfft), np.float32)
        np.add.at(P, (ii // ng, ii % ng), ww)
        Pf = np.fft.rfft2(P)
        rings = np.empty((NBIN, ng * ng), np.float32)
        for b in range(NBIN):
            rings[b] = np.fft.irfft2(Pf * kern_f[s][b], s=(nfft, nfft))[:ng, :ng].ravel()
        res[k] = rings.astype(np.float64) @ FT
    np.savez(out, res=res, wsum=wsum, keys=np.array([f"{s}:{key}" for (s, key, _, _) in sets]))
    print(f"patch {ip+1} p8={p8} {time.time()-t0:.0f}s", flush=True)


_G = {}


def stage_pairs():
    import multiprocessing as mp
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    upix, stars, keep = _setup()
    print("gaia G<10 mask: masked map pixels", (~keep).sum(), "of", len(upix), flush=True)
    reg32 = region()
    r8 = np.unique(hp.ang2pix(8, *hp.pix2ang(32, reg32, lonlat=True), lonlat=True))
    edges = [bin_edges(0.5 * (lo + hi)) for (_, lo, hi) in SLICES]
    tmax_all = max(e[-1] for e in edges)
    half = 3.9 + tmax_all
    ng = int(np.ceil(2 * half / PIXDEG))
    nfft = 1 << int(np.ceil(np.log2(ng + 2 * tmax_all / PIXDEG + 2)))
    kx = (np.arange(nfft) + nfft // 2) % nfft - nfft // 2
    KX, KY = np.meshgrid(kx, kx, indexing="ij")
    dist = np.hypot(KX, KY) * PIXDEG
    kern_f, kern_theta = [], []
    for e in edges:
        row, th = [], []
        for b in range(NBIN):
            k = ((dist >= e[b]) & (dist < e[b + 1])).astype(np.float32)
            row.append(np.fft.rfft2(k))
            th.append(dist[k > 0])
        kern_f.append(row)
        kern_theta.append(th)
    np.save(os.path.join(OUT, "kern_theta.npy"), np.array(kern_theta, dtype=object), allow_pickle=True)
    ra_p, dec_p = hp.pix2ang(NSIDE, upix, nest=True, lonlat=True)
    _G.update(upix=upix, stars=stars, keep=keep, pts=load_points(), M=np.load(os.path.join(OUT, "maps.npy"), mmap_mode="r"),
              edges=edges, half=half, ng=ng, nfft=nfft, kern_f=kern_f, ra_p=ra_p, dec_p=dec_p,
              star_core=in_region(stars["ra"], stars["dec"], reg32),
              star_p8=hp.ang2pix(8, stars["ra"], stars["dec"], lonlat=True))
    os.makedirs(os.path.join(OUT, "pairs"), exist_ok=True)
    print("patches", len(r8), "grid", ng, "fft", nfft, flush=True)
    with mp.get_context("fork").Pool(nproc) as pool:
        todo = list(enumerate(r8))
        if os.environ.get("R143_PATCHES"):
            todo = todo[: int(os.environ["R143_PATCHES"])]
        pool.map(pairs_patch, todo, chunksize=1)


# ---------------------------------------------------------------- result
def stage_result():
    files = sorted(glob.glob(os.path.join(OUT, "pairs", "p*.npz")))
    keys = list(np.load(files[0])["keys"])
    R = np.stack([np.load(f)["res"] for f in files])
    W = np.stack([np.load(f)["wsum"] for f in files])
    npch = len(files)
    ki = {k: i for i, k in enumerate(keys)}
    col = {}
    r = 0
    for key in KEYS:
        col[key] = np.arange(r, r + NG); r += NG
    gcol = np.arange(r, r + NS); r += NS
    rcol = np.arange(r, r + NS)
    kt = np.load(os.path.join(OUT, "kern_theta.npy"), allow_pickle=True)
    wm_file = os.path.join(OUT, "wm_theory.npy")
    if not os.path.exists(wm_file):
        np.save(wm_file, theory_wm(kt))
    wm = np.load(wm_file)
    pts = load_points()
    n2 = np.zeros(NS)
    for s, (t, lo, hi) in enumerate(SLICES):
        g = pts[t]["g"]
        m = (g["z"] >= lo) & (g["z"] < hi) & g["core"]
        h, e = np.histogram(g["z"][m], bins=np.linspace(lo, hi, 21), weights=g["w"][m])
        nz = h / h.sum() / np.diff(e)
        n2[s] = (nz ** 2 * np.diff(e)).sum()
    stars = stars_all()
    reg32 = region()
    sc = (stars["g"] >= 12) & (stars["g"] < 16) & in_region(stars["ra"], stars["dec"], reg32)
    nstar_sr = sc.sum() / (len(reg32) * hp.nside2pixarea(32))
    del pts
    # field pairs: name -> (numerator keys, denominator key, scale)
    PAIRS = {"real": (["s0"], "n0"), "real_sm": (["s0sm"], "n0sm"),
             "recA": (["dA"], "nA"), "recB": (["dB"], "nB"), "recC": (["dC"], "nA"), "recD": (["dD"], "nA"),
             "recE": (["dE"], "nA"), "recA_sm": (["dAsm"], "nAsm"),
             "maskA": (["mC", "mD"], "nA"), "maskC": (["mC"], "nA"), "maskD": (["mD"], "nA"), "maskE": (["mE"], "nA"),
             "truA": (["tC", "tD"], "nt"), "truC": (["tC"], "nt"), "truD": (["tD"], "nt"), "truE": (["tE"], "nt")}

    def estimate(Rs):
        out = {k: np.zeros((NS, NG)) for k in PAIRS}
        out.update(star=np.zeros((NS, NG)), star_sm=np.zeros((NS, NG)), b=np.zeros(NS), wmb=np.zeros(NS),
                   wgs=np.zeros(NS))
        for s in range(NS):
            wb = wbar_weights(s)
            wb = np.where(np.isfinite(wm[s]), wb, 0.0)
            wb = wb / wb.sum()
            V = wb > 0
            S_ = lambda x: (wb[V] * x[V]).sum(-1) if x.ndim == 1 else (wb[V, None] * x[V]).sum(0)
            g, rr, st = Rs[ki[f"{s}:gal"]], Rs[ki[f"{s}:ran"]], Rs[ki[f"{s}:star"]]
            DD, DR, RD, RR = g[:, gcol[s]], g[:, rcol[s]], rr[:, gcol[s]], rr[:, rcol[s]]
            wgg = DD * RR / (DR * RD) - 1
            wms = np.nan_to_num(wm[s])
            out["b"][s] = np.sqrt(max(S_(wgg) / (S_(wms) * n2[s]), 1e-6))
            out["wmb"][s] = S_(wms)
            SD, SR = st[:, gcol[s]], st[:, rcol[s]]
            out["wgs"][s] = S_(SD * RR / (SR * RD) - 1)

            def mI(x, nums, den):
                return sum(x[:, col[k]] for k in nums) / x[:, col[den]]

            for name, (nums, den) in PAIRS.items():
                out[name][s] = S_(mI(g, nums, den) - mI(rr, nums, den))
            out["star"][s] = S_(mI(st, ["s0"], "n0") - mI(rr, ["s0"], "n0"))
            out["star_sm"][s] = S_(mI(st, ["s0sm"], "n0sm") - mI(rr, ["s0sm"], "n0sm"))
        return out

    tot = R.sum(0)
    full = estimate(tot)
    jk = [estimate(tot - R[i]) for i in range(npch)]
    nst = W[:, ki["0:star"]]
    lo_p = nst <= np.median(nst)
    half_est = {"low": estimate(R[lo_p].sum(0)), "high": estimate(R[~lo_p].sum(0))}
    half_jk = {k: [estimate(R[m].sum(0) - R[i]) for i in np.where(m)[0]] for k, m in (("low", lo_p), ("high", ~lo_p))}
    np.savez(os.path.join(OUT, "result.npz"), **{k: full[k] for k in full},
             **{k + "_jk": np.stack([j[k] for j in jk]) for k in full},
             **{f"{h}_{k}": half_est[h][k] for h in half_est for k in full},
             **{f"{h}_{k}_jk": np.stack([j[k] for j in half_jk[h]]) for h in half_jk for k in full},
             nstar_sr=nstar_sr, npatch=npch, nst=nst)
    print("result saved", flush=True)


# ---------------------------------------------------------------- report
def jkerr(a):
    n = a.shape[0]
    return np.sqrt((n - 1) / n * ((a - a.mean(0)) ** 2).sum(0))


def stage_report():
    from scipy.stats import chi2
    J = np.load(os.path.join(OUT, "result.npz"))
    meta = np.load(os.path.join(OUT, "maps_meta.npz"))
    ge = group_edges()
    lamc = np.array([np.sqrt(a * b) for a, b in ge])
    tonw = 2.99792458e14 / lamc * 1e-11
    zc = np.array([0.5 * (lo + hi) for (_, lo, hi) in SLICES])
    L = []
    L.append(f"Images reduced with injection: {int(meta['nimg'])} (D1/D2/D3: {np.bincount(meta['dets'])[1:].tolist()}); "
             f"patches {int(J['npatch'])}; star-mask block fraction {meta['frac_sm'].mean():.2f}.\n")

    def ratio(a, b, pre=""):
        v = J[pre + a] / J[pre + b]
        vj = J[pre + a + "_jk"] / J[pre + b + "_jk"]
        return v, vj

    def ivm(v, vj, m):
        e = jkerr(vj)
        mm = m & np.isfinite(v) & np.isfinite(e) & (e > 0)
        w = 1 / e[mm] ** 2
        mj = (vj[:, mm] * w).sum(1) / w.sum()
        return (v[mm] * w).sum() / w.sum(), jkerr(mj)

    # transfer summary per group: inverse-variance mean over slices, and per tracer
    T = {}
    for k in "ABCDE":
        tru = "truA" if k == "B" else f"tru{k}"
        v, vj = ratio(f"rec{k}", tru)
        if k == "B":
            v, vj = v / 5, vj / 5
        T[k] = (v, vj)
    T["Asm"] = ratio("recA_sm", "truA")
    Mk = {k: ratio(f"mask{k}", f"tru{k}") for k in "ACDE"}
    L.append("Full-chain transfer T = w̄[recovered mock] / w̄[truth] per channel group, inverse-variance mean over the 32 "
             "slices (jackknife σ); the mask alone (mock on good pixels, no strip median) beside it; B is divided by 5:\n")
    L.append("| group µm | T_A (fiducial) | T_B (×5) | T_C (smooth) | T_D (compact) | T_E (rising z) | mask alone A | mask alone C | mask alone D | T_A, star mask |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    allm = np.ones(NS, bool)
    summ = {}
    for g in range(NG):
        cells = []
        for k in ("A", "B", "C", "D", "E"):
            mu, se = ivm(T[k][0][:, g], T[k][1][:, :, g], allm)
            summ[(k, g)] = (mu, se)
            cells.append(f"{mu:.3f} ± {se:.3f}")
        for k in ("A", "C", "D"):
            mu, se = ivm(Mk[k][0][:, g], Mk[k][1][:, :, g], allm)
            cells.append(f"{mu:.2f} ± {se:.2f}")
        mu, se = ivm(T["Asm"][0][:, g], T["Asm"][1][:, :, g], allm)
        cells.append(f"{mu:.3f} ± {se:.3f}")
        L.append(f"| {lamc[g]:.3f} | " + " | ".join(cells) + " |")
    # per slice transfer (A, C, D), averaged over groups per detector
    L.append("\nT per slice (rows), inverse-variance mean over the 6 groups of each detector (columns D1/D2/D3), for A, C, D:\n")
    L.append("| slice | T_A D1 | T_A D2 | T_A D3 | T_C D1 | T_C D2 | T_C D3 | T_D D1 | T_D D2 | T_D D3 |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for s, (t, lo, hi) in enumerate(SLICES):
        cells = []
        for k in "ACD":
            for det in range(3):
                gs = list(range(det * 6, det * 6 + 6))
                v, vj = T[k][0][s, gs], T[k][1][:, s, gs]
                e = jkerr(vj)
                okk = np.isfinite(v) & (e > 0)
                w = 1 / e[okk] ** 2
                mu = (v[okk] * w).sum() / w.sum()
                se = jkerr((vj[:, okk] * w).sum(1) / w.sum())
                cells.append(f"{mu:.2f}±{se:.2f}")
        L.append(f"| {t[:3]} {lo:.2f}-{hi:.2f} | " + " | ".join(cells) + " |")
    # the real signal on this subset: R-138's estimator, strip-only (R-138's convention is not available here), and
    # corrected by T_A and T_C; b_I d(nuI_nu)/dz = w̄ / (b_r w̄_m T)
    b, wmb = J["b"], J["wmb"]
    bj, wmbj = J["b_jk"], J["wmb_jk"]
    real, realj = J["real"], J["real_jk"]
    # rule fixed at 02:00Z: per-cell T if it is precise enough, else T per slice averaged over each detector's groups
    noisy = float((jkerr(T["C"][1]) > 0.2).mean())
    per_det = noisy > 0.25
    L.append(f"\nFraction of (slice, group) cells with σ(T_C) > 0.2: {noisy:.2f}; the correction uses "
             + ("T per slice averaged over each detector's groups (rule of 02:00Z)." if per_det else "T per cell."))

    def t_used(k):
        v, vj = T[k]
        if not per_det:
            return v, vj
        v2, vj2 = np.empty_like(v), np.empty_like(vj)
        for det in range(3):
            gs = list(range(det * 6, det * 6 + 6))
            e = jkerr(vj[:, :, gs])
            w = np.where(np.isfinite(e) & (e > 0), 1 / e ** 2, 0.0)
            vv = np.nan_to_num(v[:, gs])
            v2[:, gs] = ((vv * w).sum(1) / w.sum(1))[:, None]
            vj2[:, :, gs] = ((np.nan_to_num(vj[:, :, gs]) * w).sum(2) / w.sum(1))[:, :, None]
        return v2, vj2

    TU = {k: t_used(k) for k in ("A", "C", "D")}
    out = {}
    for k in ("A", "C", "D"):
        v = real / (b[:, None] * wmb[:, None] * TU[k][0]) * tonw
        vj = realj / (bj[:, :, None] * wmbj[:, :, None] * TU[k][1]) * tonw
        out[k] = (v, vj)
    raw = real / (b[:, None] * wmb[:, None]) * tonw
    rawj = realj / (bj[:, :, None] * wmbj[:, :, None]) * tonw
    L.append("\nClustering bias of each slice on this subset (b_r, jackknife): " +
             ", ".join(f"{t[:3]} {lo:.2f} {b[s]:.2f}±{jkerr(bj[:, s]):.2f}" for s, (t, lo, _) in enumerate(SLICES)) + ".\n")
    L.append("Corrected b_I d(νI_ν)/dz per group (nW m^-2 sr^-1 per unit z), inverse-variance mean over slices (jackknife σ, "
             "which carries the transfer's own error): uncorrected (w̄ / b_r w̄_m, no transfer), corrected by T_A "
             "(fiducial), by T_C (all smooth: the light outside catalogued sources), by T_D (all compact); z < 0.5 and "
             "z ≥ 0.5 for T_C:\n")
    L.append("| group µm | uncorrected | /T_A | /T_C | /T_D | /T_C z < 0.5 | /T_C z ≥ 0.5 | S/N (T_C) |")
    L.append("|---|---|---|---|---|---|---|---|")
    corr_summary = {}
    for g in range(NG):
        cells = []
        mu, se = ivm(raw[:, g], rawj[:, :, g], allm); cells.append(f"{mu:.2f} ± {se:.2f}")
        for k in "ACD":
            mu, se = ivm(out[k][0][:, g], out[k][1][:, :, g], allm)
            cells.append(f"{mu:.1f} ± {se:.1f}")
            corr_summary[(k, g)] = (mu, se)
        for m in (zc < 0.5, zc >= 0.5):
            mu, se = ivm(out["C"][0][:, g], out["C"][1][:, :, g], m); cells.append(f"{mu:.1f} ± {se:.1f}")
        mu, se = corr_summary[("C", g)]
        cells.append(f"{mu/se:.1f}")
        L.append(f"| {lamc[g]:.3f} | " + " | ".join(cells) + " |")
    # per slice corrected table (T_C), and chi2 against zero per group
    L.append("\nPer slice b_I d(νI_ν)/dz corrected by T_C (the light outside catalogued sources), nW m^-2 sr^-1 per unit z, "
             "value ± jackknife σ:\n")
    L.append("| slice | " + " | ".join(f"{x:.3f}" for x in lamc) + " |")
    L.append("|---|" + "---|" * NG)
    eC = jkerr(out["C"][1])
    for s, (t, lo, hi) in enumerate(SLICES):
        L.append(f"| {t[:3]} {lo:.2f}-{hi:.2f} | " + " | ".join(f"{out['C'][0][s, g]:.0f}±{eC[s, g]:.0f}" for g in range(NG)) + " |")
    # nulls: stars, with and without the fixed star mask; galaxy signal with the star mask
    L.append("\nStellar null (12 < G < 16 stars as centres, R-138's estimator and units with T_C), χ² against zero over 32 "
             "slices (jackknife σ), without and with the fixed-flux star mask; the galaxy signal with the star mask; the "
             "galaxy-star cross-correlation w̄_gs per slice and the contamination it implies, max over slices of "
             "|w̄_gs N_* S| / σ(signal):\n")
    L.append("| group µm | star null χ² (R-138 masks) | star null χ² (star mask) | signal /T_C, star mask | max contamination / σ |")
    L.append("|---|---|---|---|---|")
    conv = lambda key: (J[key] / (b[:, None] * wmb[:, None] * TU["C"][0]) * tonw,
                        J[key + "_jk"] / (bj[:, :, None] * wmbj[:, :, None] * TU["C"][1]) * tonw)
    st, stj = conv("star")
    ss, ssj = conv("star_sm")
    rs, rsj = conv("real_sm")
    est, ess = jkerr(stj), jkerr(ssj)
    th_max = np.array([bin_edges(z)[-1] for z in zc])
    nst_in = float(J["nstar_sr"]) * np.pi * np.radians(th_max) ** 2
    cont = J["wgs"][:, None] * nst_in[:, None] * st
    for g in range(NG):
        okk = est[:, g] > 0
        c1 = float(((st[okk, g] / est[okk, g]) ** 2).sum())
        okk2 = ess[:, g] > 0
        c2 = float(((ss[okk2, g] / ess[okk2, g]) ** 2).sum())
        mu, se = ivm(rs[:, g], rsj[:, :, g], allm)
        cm = np.nanmax(np.abs(cont[:, g]) / eC[:, g])
        L.append(f"| {lamc[g]:.3f} | {c1:.0f}/{okk.sum()} (p={chi2.sf(c1, okk.sum()):.1e}) | {c2:.0f}/{okk2.sum()} "
                 f"(p={chi2.sf(c2, okk2.sum()):.1e}) | {mu:.1f} ± {se:.1f} | {cm:.2f} |")
    wgs, wgsj = J["wgs"], J["wgs_jk"]
    L.append("\nw̄_gs per slice (galaxies around 12 < G < 16 stars, Hamilton estimator, same annulus): " +
             ", ".join(f"{t[:3]} {lo:.2f} {wgs[s]:+.4f}±{jkerr(wgsj[:, s]):.4f}" for s, (t, lo, _) in enumerate(SLICES)) +
             f". Stars per sr {float(J['nstar_sr']):.3g}.")
    # star-density split predicted by the transfer
    L.append("\nStar-density split (patches at the median 12 < G < 16 count): T_A per half and the predicted low/high ratio "
             "T_A(low)/T_A(high); the observed low/high ratio of the uncorrected signal on this subset, and R-138's:\n")
    R138_SPLIT = [(14.2, 6.4), (9.7, 3.8), (12.3, 4.1), (13.6, 5.3), (6.2, 2.8), (31.4, 47.7), (9.3, 6.5), (10.8, 7.2),
                  (8.5, 8.9), (9.4, 7.4), (9.4, 5.7), (6.9, 6.0), (7.1, 5.2), (6.8, 5.4), (6.3, 6.1), (5.2, 6.2),
                  (5.4, 6.5), (6.4, 5.4)]
    L.append("| group µm | T_A low-star | T_A high-star | predicted ratio | observed ratio (this subset) | R-138 ratio |")
    L.append("|---|---|---|---|---|---|")
    for g in range(NG):
        cells = []
        tv = {}
        for h in ("low", "high"):
            v = J[f"{h}_recA"][:, g] / J[f"{h}_truA"][:, g]
            vj = J[f"{h}_recA_jk"][:, :, g] / J[f"{h}_truA_jk"][:, :, g]
            mu, se = ivm(v, vj, allm)
            tv[h] = (mu, se, v, vj)
            cells.append(f"{mu:.3f} ± {se:.3f}")
        pr = tv["low"][0] / tv["high"][0]
        pre = pr * np.hypot(tv["low"][1] / tv["low"][0], tv["high"][1] / tv["high"][0])
        cells.append(f"{pr:.2f} ± {pre:.2f}")
        ob = {}
        for h in ("low", "high"):
            v = J[f"{h}_real"][:, g] / (J[f"{h}_b"] * J[f"{h}_wmb"])
            vj = J[f"{h}_real_jk"][:, :, g] / (J[f"{h}_b_jk"] * J[f"{h}_wmb_jk"])
            ob[h] = ivm(v, vj, allm)
        orat = ob["low"][0] / ob["high"][0]
        ore = abs(orat) * np.hypot(ob["low"][1] / ob["low"][0], ob["high"][1] / ob["high"][0])
        cells.append(f"{orat:.2f} ± {ore:.2f}")
        a, c = R138_SPLIT[g]
        cells.append(f"{a / c:.2f}")
        L.append(f"| {lamc[g]:.3f} | " + " | ".join(cells) + " |")
    # CIBER
    L.append("\nCIBER comparison (arXiv:2608.12116, read off its Fig. a2h_vs_redshift as R-138): this subset, inverse-variance "
             "mean of the slices whose centre lies in the CIBER bin, corrected by T_C and by T_D (the range for the total):\n")
    L.append("| λ | z | CIBER | group µm | /T_C | /T_D |")
    L.append("|---|---|---|---|---|---|")
    for lam_c, rows in CIBER.items():
        gsel = [g for g, (a, bb) in enumerate(ge) if a <= lam_c < bb] + ([7 - 1] if lam_c == 1.1 else [])
        gsel = sorted(set(gsel))
        for g in gsel:
            for (lo, hi, v, e) in rows:
                m = (zc >= lo) & (zc < hi)
                cells = []
                for k in "CD":
                    mu, se = ivm(out[k][0][:, g], out[k][1][:, :, g], m)
                    cells.append(f"{mu:.1f} ± {se:.1f}")
                cib = f"{v} ± {e}" if v is not None else f"< {e} (95%)"
                L.append(f"| {lam_c} | {lo}-{hi} | {cib} | {lamc[g]:.3f} | " + " | ".join(cells) + " |")
    txt = "\n".join(L)
    open(os.path.join(OUT, "report.md"), "w").write(txt)
    json.dump({f"{k}_{g}": v for (k, g), v in summ.items()}, open(os.path.join(OUT, "transfer_summary.json"), "w"))
    print(txt)


if __name__ == "__main__":
    {"maps": stage_maps, "pairs": stage_pairs, "result": stage_result, "report": stage_report}[sys.argv[1]]()
