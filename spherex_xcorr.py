"""R-138: SPHEREx QR2 x DESI DR1 clustering-redshift tomography of the near-infrared background (b_I dI/dz per channel).

Stages (each caches its output under data/spherex/):
  maps    : accumulate the reduced images (spherex_reduce.py outputs) into NSIDE-2048 maps per channel and half-survey
            (A/B, split at the median exposure time), the ZODI-model null map, and the strip-filter transfer maps (each
            slice's galaxy overdensity passed through the identical per-(image, channel) mean subtraction).
  pairs   : per NSIDE-8 patch, gnomonic 1.72' grids; ring-convolve each point set (galaxies, randoms, shuffled galaxies,
            Gaia stars) per slice annulus bin, dot with every field; save the pair sums.
  result  : the estimator of arXiv:1810.00885 with the slice's own clustering bias, jackknife over patches, nulls.

Usage: python -u spherex_xcorr.py maps|pairs|result
Fixed choices are those of results-R138.md, "Fixed before any cross-correlation".
"""
import glob
import json
import os
import sys
import time

import numpy as np
import healpy as hp

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")
SP = os.path.join(D, "spherex")
NSIDE = 2048
NCH = 51
PIXDEG = np.degrees(hp.nside2resol(NSIDE))  # 0.0286 deg
NBIN = 4
SEED = 138

SLICES = []  # (tracer, zlo, zhi)
for t, lo, hi in (("BGS_ANY", 0.05, 0.45), ("LRG", 0.4, 1.1), ("ELG_LOPnotqso", 0.8, 1.6), ("QSO", 0.8, 2.1)):
    n = int(round((hi - lo) / 0.1))
    for i in range(n):
        SLICES.append((t, round(lo + 0.1 * i, 2), round(lo + 0.1 * (i + 1), 2)))
NS = len(SLICES)
TRACERS = ["BGS_ANY", "LRG", "ELG_LOPnotqso", "QSO"]


def chi_of_z(z):
    from astropy.cosmology import Planck18
    return Planck18.comoving_distance(z).value  # Mpc


def annulus(zc):
    chi = chi_of_z(zc)
    tmin = max(np.degrees(1.0 / chi), 3.5 / 60)
    tmax = np.degrees(10.0 / chi)
    return tmin, tmax


def bin_edges(zc):
    tmin, tmax = annulus(zc)
    return np.geomspace(tmin, tmax, NBIN + 1)


def region():
    reg32 = np.load(os.path.join(SP, "region32.npy"))
    return reg32


def in_region(ra, dec, reg32):
    m = np.zeros(hp.nside2npix(32), bool)
    m[reg32] = True
    return m[hp.ang2pix(32, ra, dec, lonlat=True)]


def load_points():
    reg32 = region()
    out = {}
    rng = np.random.default_rng(SEED)
    for t in TRACERS:
        g = np.load(os.path.join(D, "desi_dr1", f"{t}_dat.npz"))
        r = np.load(os.path.join(D, "desi_dr1", f"{t}_ran.npz"))
        # partners anywhere (the patch grids reach beyond the field); centres only inside the field
        out[t] = dict(g=dict(ra=g["ra"], dec=g["dec"], z=g["z"], w=g["w"], core=in_region(g["ra"], g["dec"], reg32),
                             p8=hp.ang2pix(8, g["ra"], g["dec"], lonlat=True)),
                      r=dict(ra=r["ra"], dec=r["dec"], z=r["z"], w=r["w"], core=in_region(r["ra"], r["dec"], reg32),
                             p8=hp.ang2pix(8, r["ra"], r["dec"], lonlat=True)))
        zs = g["z"].copy()
        c = out[t]["g"]["core"]
        zs[c] = rng.permutation(zs[c])
        out[t]["g"]["zshuf"] = zs
    return out


def load_stars():
    fs = sorted(glob.glob(os.path.join(D, "gaia", "g_*.csv")))
    if not fs:
        return None
    arr = []
    for f in fs:
        try:
            a = np.loadtxt(f, delimiter=",", skiprows=1, ndmin=2)
        except Exception:
            continue
        if a.size:
            arr.append(a)
    if not arr:
        return None
    a = np.concatenate(arr)
    return dict(ra=a[:, 0], dec=a[:, 1], g=a[:, 2])


# ---------------------------------------------------------------- maps
def stage_maps():
    files = sorted(glob.glob(os.path.join(SP, "red", "*.npz")))
    print("images", len(files), flush=True)
    ts = []
    pix_all = []
    for f in files:
        d = np.load(f)
        ts.append(float(d["t"]))
        pix_all.append(np.unique(d["hp"]))
    tmed = float(np.median(ts))
    upix = np.unique(np.concatenate(pix_all))
    del pix_all
    print("map pixels", len(upix), "area deg2", len(upix) * hp.nside2pixarea(NSIDE, True), flush=True)
    lut = np.full(hp.nside2npix(NSIDE), -1, np.int32)
    lut[upix] = np.arange(len(upix), dtype=np.int32)
    npx = len(upix)
    S = np.zeros((2, NCH, npx), np.float32)
    N = np.zeros((2, NCH, npx), np.float32)
    Z = np.zeros((NCH, npx), np.float32)
    # transfer templates: absolute galaxy overdensity per slice, on the map pixels
    pts = load_points()
    T = np.zeros((NS, npx), np.float32)
    for s, (t, lo, hi) in enumerate(SLICES):
        g, r = pts[t]["g"], pts[t]["r"]
        mg = (g["z"] >= lo) & (g["z"] < hi)
        mr = (r["z"] >= lo) & (r["z"] < hi)
        alpha = g["w"][mg].sum() / r["w"][mr].sum()
        pg = lut[hp.ang2pix(NSIDE, g["ra"][mg], g["dec"][mg], lonlat=True, nest=True)]
        pr = lut[hp.ang2pix(NSIDE, r["ra"][mr], r["dec"][mr], lonlat=True, nest=True)]
        T[s] += np.bincount(pg[pg >= 0], weights=g["w"][mg][pg >= 0], minlength=npx)
        T[s] -= alpha * np.bincount(pr[pr >= 0], weights=r["w"][mr][pr >= 0], minlength=npx)
    del pts
    TF = np.zeros((NS, npx), np.float32)
    TN = np.zeros(npx, np.float32)
    rng = np.random.default_rng(SEED)
    tsub = set(rng.choice(len(files), size=min(3000, len(files)), replace=False).tolist())
    t0 = time.time()
    for k, f in enumerate(files):
        d = np.load(f)
        idx = lut[d["hp"]]
        ch = d["ch"].astype(np.int64)
        half = 0 if ts[k] < tmed else 1
        flat = ch * npx + idx
        np.add.at(S[half].reshape(-1), flat, d["s"])
        np.add.at(N[half].reshape(-1), flat, d["n"])
        np.add.at(Z.reshape(-1), flat, d["z"])
        if k in tsub:
            n = d["n"]
            Tv = T[:, idx]  # (NS, e)
            for c in np.unique(ch):
                m = ch == c
                mean = (Tv[:, m] * n[m]).sum(1) / n[m].sum()
                np.add.at(TF, (slice(None), idx[m]), ((Tv[:, m] - mean[:, None]) * n[m]))
            np.add.at(TN, idx, n)
        if k % 2000 == 0:
            print(k, f"{time.time()-t0:.0f}s", flush=True)
    np.savez(os.path.join(SP, "maps.npz"), upix=upix, S=S, N=N, Z=Z, T=T, TF=TF, TN=TN, tmed=tmed, nimg=len(files))
    print("maps saved", flush=True)


# ---------------------------------------------------------------- pairs
def gnomonic(ra, dec, ra0, dec0):
    ra, dec, ra0, dec0 = map(np.radians, (ra, dec, ra0, dec0))
    cosc = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    x = np.cos(dec) * np.sin(ra - ra0) / cosc
    y = (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cosc
    return np.degrees(x), np.degrees(y), cosc


def stage_pairs():
    M = np.load(os.path.join(SP, "maps.npz"))
    upix = M["upix"]
    S, N, Z, T, TF, TN = M["S"], M["N"], M["Z"], M["T"], M["TF"], M["TN"]
    stars = load_stars()
    # Gaia bright-star map-level mask
    if stars is not None:
        # amended before any cross-correlation (results-R138.md): G < 10, r = 2' x 10^(0.2 (10 - G))
        b = stars["g"] < 10
        rad = np.maximum(2.0 / 60 * 10 ** (0.2 * (10 - stars["g"][b])), PIXDEG) if b.any() else []
        vec = hp.ang2vec(stars["ra"][b], stars["dec"][b], lonlat=True)
        bad = []
        for v, rr in zip(vec, rad):
            bad.append(hp.query_disc(NSIDE, v, np.radians(rr), nest=True, inclusive=True))
        badpix = np.unique(np.concatenate(bad)) if bad else np.array([], int)
        keep = ~np.isin(upix, badpix)
        print("gaia mask: stars", b.sum(), "masked map pixels", (~keep).sum(), flush=True)
    else:
        keep = np.ones(len(upix), bool)
        print("gaia mask: none (no star file)", flush=True)
    ra_p, dec_p = hp.pix2ang(NSIDE, upix, nest=True, lonlat=True)
    if stars is not None:
        star_core = (stars["g"] >= 12) & (stars["g"] < 16) & in_region(stars["ra"], stars["dec"], region())
        star_p8 = hp.ang2pix(8, stars["ra"], stars["dec"], lonlat=True)
    pts = load_points()
    reg32 = region()
    # patches: NSIDE-8 pixels holding field pixels
    r8 = np.unique(hp.ang2pix(8, *hp.pix2ang(32, reg32, lonlat=True), lonlat=True))
    edges = [bin_edges(0.5 * (lo + hi)) for (_, lo, hi) in SLICES]
    tmax_all = max(e[-1] for e in edges)
    half = 3.9 + tmax_all  # deg; an NSIDE-8 pixel lies within ~3.8 deg of its centre
    ng = int(np.ceil(2 * half / PIXDEG))
    nfft = 1 << int(np.ceil(np.log2(ng + 2 * tmax_all / PIXDEG + 2)))
    print("patches", len(r8), "grid", ng, "fft", nfft, flush=True)
    # field list (rows of the field matrix)
    fnames = ([f"SA{c}" for c in range(NCH)] + [f"SB{c}" for c in range(NCH)] + [f"NA{c}" for c in range(NCH)]
              + [f"NB{c}" for c in range(NCH)] + [f"Z{c}" for c in range(NCH)] + [f"T{s}" for s in range(NS)]
              + [f"TF{s}" for s in range(NS)] + ["TN"] + [f"G{s}" for s in range(NS)] + [f"R{s}" for s in range(NS)])
    # ring kernels per slice and bin, and the pair separations they hold (for the theory average)
    kx = (np.arange(nfft) + nfft // 2) % nfft - nfft // 2
    KX, KY = np.meshgrid(kx, kx, indexing="ij")
    dist = np.hypot(KX, KY) * PIXDEG
    kern_f = []
    kern_theta = []
    for e in edges:
        row = []
        th = []
        for b in range(NBIN):
            k = ((dist >= e[b]) & (dist < e[b + 1])).astype(np.float32)
            row.append(np.fft.rfft2(k))
            th.append(dist[k > 0])
        kern_f.append(row)
        kern_theta.append(th)
    np.save(os.path.join(SP, "kern_theta.npy"), np.array(kern_theta, dtype=object), allow_pickle=True)
    outdir = os.path.join(SP, "pairs")
    os.makedirs(outdir, exist_ok=True)
    if os.environ.get("PAIRS_ONLY"):
        r8 = r8[: int(os.environ["PAIRS_ONLY"])]
    for ip, p8 in enumerate(r8):
        out = os.path.join(outdir, f"p{p8}.npz")
        if os.path.exists(out):
            continue
        t0 = time.time()
        ra0, dec0 = hp.pix2ang(8, p8, lonlat=True)

        def grid_index(ra, dec):
            x, y, cc = gnomonic(ra, dec, ra0, dec0)
            i = np.floor((x + half) / PIXDEG).astype(int)
            j = np.floor((y + half) / PIXDEG).astype(int)
            ok = (cc > 0) & (i >= 0) & (i < ng) & (j >= 0) & (j < ng)
            return i * ng + j, ok

        gi, ok = grid_index(ra_p, dec_p)
        ok &= keep
        sel = np.where(ok)[0]
        gsel = gi[sel]
        F = np.zeros((len(fnames), ng * ng), np.float32)

        def put(row, vals):
            F[row] = np.bincount(gsel, weights=vals[sel], minlength=ng * ng)

        r = 0
        for arr in (S[0], S[1], N[0], N[1], Z):
            for c in range(NCH):
                put(r, arr[c]); r += 1
        for s in range(NS):
            put(r, T[s] * TN); r += 1   # unfiltered template, same weights as the filtered one
        for s in range(NS):
            put(r, TF[s]); r += 1
        put(r, TN); r += 1
        # partner galaxy and random count maps per slice
        # rows G0..G{NS-1} then R0..R{NS-1}, in the order of fnames
        for key in ("g", "r"):
            for s, (t, lo, hi) in enumerate(SLICES):
                cat = pts[t][key]
                m = (cat["z"] >= lo) & (cat["z"] < hi)
                ii, okk = grid_index(cat["ra"][m], cat["dec"][m])
                F[r] = np.bincount(ii[okk], weights=cat["w"][m][okk], minlength=ng * ng)
                r += 1
        assert r == len(fnames)
        FT = F.T.astype(np.float64)
        del F
        # centre sets: galaxies, randoms, shuffled galaxies (core of this patch only), Gaia stars 12<G<16
        sets = []
        for s, (t, lo, hi) in enumerate(SLICES):
            g, rr = pts[t]["g"], pts[t]["r"]
            for key, cat, zkey in (("gal", g, "z"), ("ran", rr, "z"), ("shuf", g, "zshuf")):
                m = (cat[zkey] >= lo) & (cat[zkey] < hi) & cat["core"]
                m &= cat["p8"] == p8
                ii, okk = grid_index(cat["ra"][m], cat["dec"][m])
                sets.append((s, key, ii[okk], cat["w"][m][okk]))
            if stars is not None:
                m = star_core & (star_p8 == p8)
                ii, okk = grid_index(stars["ra"][m], stars["dec"][m])
                sets.append((s, "star", ii[okk], np.ones(okk.sum())))
        res = np.zeros((len(sets), NBIN, len(fnames)), np.float64)
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
                # correlation: sum_x P(x) F(x + d) = (P_reflected * K) . F ; K is symmetric so convolve P with K
                rings[b] = np.fft.irfft2(Pf * kern_f[s][b], s=(nfft, nfft))[:ng, :ng].ravel()
            res[k] = rings.astype(np.float64) @ FT
        np.savez(out, res=res, wsum=wsum, keys=np.array([f"{s}:{key}" for (s, key, _, _) in sets]),
                 fnames=np.array(fnames))
        print(f"patch {ip+1}/{len(r8)} p8={p8} sets={len(sets)} {time.time()-t0:.0f}s", flush=True)


# ---------------------------------------------------------------- theory
def theory_wm(kern_theta):
    """w_m(theta, z) of 1810.00885 eq. w_m (per unit z), averaged over each kernel's pixel separations."""
    import camb
    from scipy.special import j0
    from astropy.cosmology import Planck18
    pars = camb.set_params(H0=67.4, ombh2=0.0224, omch2=0.120, ns=0.965, As=2.1e-9, mnu=0.06)
    zs = sorted(set([0.5 * (lo + hi) for (_, lo, hi) in SLICES]))
    PK = camb.get_matter_power_interpolator(pars, zmin=0, zmax=2.5, nz_step=50, kmax=60, nonlinear=True,
                                            hubble_units=False, k_hunit=False)
    k = np.geomspace(1e-4, 50, 6000)
    out = np.zeros((NS, NBIN))
    for s, (_, lo, hi) in enumerate(SLICES):
        zc = 0.5 * (lo + hi)
        chi = Planck18.comoving_distance(zc).value
        hc = (Planck18.H(zc).value / 299792.458)  # 1/Mpc
        P = PK.P(zc, k)
        for b in range(NBIN):
            th = np.radians(np.asarray(kern_theta[s][b]))
            if th.size == 0:  # bin narrower than the pixel lattice: dropped, weights renormalised
                out[s, b] = np.nan
                continue
            tg = np.geomspace(th.min() * 0.999, th.max() * 1.001, 40)
            wg = np.array([np.trapezoid(k * P * j0(k * t * chi), k) for t in tg]) / (2 * np.pi) * hc
            out[s, b] = np.interp(th, tg, wg).mean()
    return out


def wbar_weights(s):
    e = bin_edges(0.5 * (SLICES[s][1] + SLICES[s][2]))
    w = (e[1:] ** 0.2 - e[:-1] ** 0.2)
    return w / w.sum()


def ch_lambda():
    from spherex_reduce import EDGES
    lam = []
    for det in (1, 2, 3):
        lo, hi = EDGES[det]
        e = lo * (hi / lo) ** (np.arange(18) / 17)
        lam += list(np.sqrt(e[:-1] * e[1:]))
    return np.array(lam)


GROUPS = []
for det in range(3):
    for a, b in ((0, 3), (3, 6), (6, 9), (9, 12), (12, 15), (15, 17)):
        GROUPS.append(list(range(det * 17 + a, det * 17 + b)))


def stage_result():
    sys.path.insert(0, HERE)
    files = sorted(glob.glob(os.path.join(SP, "pairs", "p*.npz")))
    d0 = np.load(files[0])
    fn = list(d0["fnames"])
    keys = list(d0["keys"])
    R = np.stack([np.load(f)["res"] for f in files])  # (patch, set, bin, field)
    W = np.stack([np.load(f)["wsum"] for f in files])
    npch = len(files)
    ki = {k: i for i, k in enumerate(keys)}
    fi = {f: i for i, f in enumerate(fn)}
    kt = np.load(os.path.join(SP, "kern_theta.npy"), allow_pickle=True)
    wm_file = os.path.join(SP, "wm_theory.npy")
    if os.path.exists(wm_file):
        wm = np.load(wm_file)
    else:
        wm = theory_wm(kt)
        np.save(wm_file, wm)
    lam = ch_lambda()
    nu = 2.99792458e14 / lam
    # n(z) integral per slice from the core galaxies
    pts = load_points()
    n2 = np.zeros(NS)
    ngal = np.zeros(NS)
    for s, (t, lo, hi) in enumerate(SLICES):
        g = pts[t]["g"]
        m = (g["z"] >= lo) & (g["z"] < hi) & g["core"]
        h, e = np.histogram(g["z"][m], bins=np.linspace(lo, hi, 21), weights=g["w"][m])
        nz = h / h.sum() / np.diff(e)
        n2[s] = (nz ** 2 * np.diff(e)).sum()
        ngal[s] = m.sum()
    del pts
    cols = {}
    for pre in ("SA", "SB", "NA", "NB", "Z"):
        cols[pre] = [fi[f"{pre}{c}"] for c in range(NCH)]

    def estimate(Rs):
        """Rs: (set, bin, field) summed over the patches kept. Returns dict of per-slice arrays."""
        out = dict(bdidz=np.zeros((NS, NCH)), bdidz_grp=np.zeros((NS, len(GROUPS))), b=np.zeros(NS),
                   tr=np.zeros(NS), shuf=np.zeros((NS, len(GROUPS))), star=np.zeros((NS, len(GROUPS))),
                   zodi=np.zeros((NS, len(GROUPS))), half=np.zeros((NS, len(GROUPS))))
        TR_ELG = {}
        for s in range(NS):
            wb = wbar_weights(s)
            wb = np.where(np.isfinite(wm[s]), wb, 0.0)
            wb = wb / wb.sum()
            wms = np.nan_to_num(wm[s])
            g, r, sh = Rs[ki[f"{s}:gal"]], Rs[ki[f"{s}:ran"]], Rs[ki[f"{s}:shuf"]]
            st = Rs[ki[f"{s}:star"]] if f"{s}:star" in ki else None
            # clustering bias (Hamilton estimator)
            DD, DR = g[:, fi[f"G{s}"]], g[:, fi[f"R{s}"]]
            RD, RR = r[:, fi[f"G{s}"]], r[:, fi[f"R{s}"]]
            wgg = DD * RR / (DR * RD) - 1
            V = wb > 0
            S_ = lambda x: (wb[V] * x[V]).sum()
            bb2 = S_(wgg) / (S_(wms) * n2[s])
            out["b"][s] = np.sqrt(max(bb2, 1e-6))
            wmb = S_(wms)
            # transfer of the strip filter
            TNi = fi["TN"]
            wT = g[:, fi[f"T{s}"]] / g[:, TNi] - r[:, fi[f"T{s}"]] / r[:, TNi]
            wTF = g[:, fi[f"TF{s}"]] / g[:, TNi] - r[:, fi[f"TF{s}"]] / r[:, TNi]
            tr = S_(wTF) / S_(wT)
            out["tr"][s] = tr
            if SLICES[s][0] == "QSO":  # amendment: the sparse QSO template is noise; use the same-z ELG slice's
                zl = min(SLICES[s][1], 1.5)
                se = [k for k, (t2, lo2, _) in enumerate(SLICES) if t2 == "ELG_LOPnotqso" and abs(lo2 - zl) < 1e-6][0]
                tr = TR_ELG[se]
                out["tr"][s] = tr
            TR_ELG[s] = tr
            conv = 1.0 / (out["b"][s] * wmb * tr)

            def mean_I(x, chans, pre=("SA", "SB"), den=("NA", "NB")):
                num = sum(x[:, cols[p]][:, chans].sum(1) for p in pre)
                dd = sum(x[:, cols[p]][:, chans].sum(1) for p in den)
                return num / dd

            for c in range(NCH):
                w = mean_I(g, [c]) - mean_I(r, [c])
                out["bdidz"][s, c] = S_(w) * conv
            for gi_, chans in enumerate(GROUPS):
                wbar = lambda x: S_(mean_I(x, chans) - mean_I(r, chans)) * conv
                out["bdidz_grp"][s, gi_] = wbar(g)
                out["shuf"][s, gi_] = wbar(sh)
                if st is not None:
                    out["star"][s, gi_] = wbar(st)
                out["zodi"][s, gi_] = S_(mean_I(g, chans, ("Z",)) - mean_I(r, chans, ("Z",))) * conv
                wa = mean_I(g, chans, ("SA",), ("NA",)) - mean_I(r, chans, ("SA",), ("NA",))
                wbB = mean_I(g, chans, ("SB",), ("NB",)) - mean_I(r, chans, ("SB",), ("NB",))
                out["half"][s, gi_] = S_(wa - wbB) / 2 * conv
        return out

    tot = R.sum(0)
    full = estimate(tot)
    jk = [estimate(tot - R[i]) for i in range(npch)]
    # added after the interim look (labelled so in results-R138.md): split the patches at the median Gaia 12<G<16 star
    # count; a signal carried by star-associated light would scale with stellar density
    if "0:star" in ki:
        nst = W[:, ki["0:star"]]
        lo_p = nst <= np.median(nst)
        split = {"lowstar": estimate(R[lo_p].sum(0)), "highstar": estimate(R[~lo_p].sum(0))}
        split_jk = {k: [estimate(R[m].sum(0) - R[i]) for i in np.where(m)[0]]
                    for k, m in (("lowstar", lo_p), ("highstar", ~lo_p))}
        np.savez(os.path.join(SP, "result_split.npz"), nstar=nst,
                 **{f"{k}": split[k]["bdidz_grp"] for k in split},
                 **{f"{k}_jk": np.stack([j["bdidz_grp"] for j in split_jk[k]]) for k in split})

    def jkerr(key):
        a = np.stack([j[key] for j in jk])
        return np.sqrt((npch - 1) / npch * ((a - a.mean(0)) ** 2).sum(0)), a

    res = dict(slices=SLICES, ngal=ngal.tolist(), npatch=npch, lam=lam.tolist(),
               groups=[[float(lam[g[0]]), float(lam[g[-1]])] for g in GROUPS])
    for key in full:
        e, a = jkerr(key)
        res[key] = full[key].tolist()
        res[key + "_err"] = e.tolist()
    np.savez(os.path.join(SP, "result_jk.npz"), **{k + "_jk": np.stack([j[k] for j in jk]) for k in full},
             **{k: full[k] for k in full})
    json.dump(res, open(os.path.join(SP, "result.json"), "w"))
    print("result saved", flush=True)


# ---------------------------------------------------------------- report
LINES = {"Halpha": 0.6563, "OIII": 0.5007, "Paalpha": 1.8751}
# arXiv:2608.12116 Fig. a2h_vs_redshift (CIBER x DESI-LS), b_I dI/dz in nW m^-2 sr^-1, read off the figure:
# (zlo, zhi, value, err) or (zlo, zhi, None, 95% upper limit)
CIBER = {1.1: [(0.0, 0.2, 100, 30), (0.2, 0.4, 90, 18), (0.4, 0.6, 23, 9), (0.6, 0.8, None, 5.5), (0.8, 1.0, None, 4.0)],
         1.8: [(0.0, 0.2, 140, 40), (0.2, 0.4, 63, 15), (0.4, 0.6, 35, 8), (0.6, 0.8, None, 14), (0.8, 1.0, None, 7.5)]}


def group_edges():
    from spherex_reduce import EDGES
    out = []
    for det in (1, 2, 3):
        lo, hi = EDGES[det]
        e = lo * (hi / lo) ** (np.arange(18) / 17)
        for a, b in ((0, 3), (3, 6), (6, 9), (9, 12), (12, 15), (15, 17)):
            out.append((e[a], e[b]))
    return out


def line_frac(lo_z, hi_z, lam_rest, lam_lo, lam_hi):
    zz = np.linspace(lo_z, hi_z, 201)
    lo_ = lam_rest * (1 + zz)
    return float(((lo_ >= lam_lo) & (lo_ < lam_hi)).mean())


def stage_report():
    sys.path.insert(0, HERE)
    r = json.load(open(os.path.join(SP, "result.json")))
    J = np.load(os.path.join(SP, "result_jk.npz"))
    npch = r["npatch"]
    ge = group_edges()
    lamc = np.array([np.sqrt(a * b) for a, b in ge])
    nu = 2.99792458e14 / lamc
    tonw = nu * 1e-11  # MJy/sr -> nW m^-2 sr^-1
    val = np.array(r["bdidz_grp"]) * tonw
    jk = J["bdidz_grp_jk"] * tonw

    def err_of(a):
        return np.sqrt((npch - 1) / npch * ((a - a.mean(0)) ** 2).sum(0))

    err = err_of(jk)
    zc = np.array([0.5 * (lo + hi) for (_, lo, hi) in SLICES])
    L = []
    L.append(f"Images reduced: see maps.npz; patches (jackknife regions): {npch}.\n")
    L.append("Per slice: tracer, z range, N galaxies in the field, clustering bias b_r (jackknife), strip-filter transfer.\n")
    L.append("| slice | tracer | z | N_gal | b_r | transfer |")
    L.append("|---|---|---|---|---|---|")
    for s, (t, lo, hi) in enumerate(SLICES):
        L.append(f"| {s} | {t} | {lo:.2f}-{hi:.2f} | {int(r['ngal'][s])} | {r['b'][s]:.2f} ± {r['b_err'][s]:.2f} | "
                 f"{r['tr'][s]:.2f} ± {r['tr_err'][s]:.2f} |")
    L.append("\nb_I d(νI_ν)/dz per channel group (columns, central µm) and slice (rows), nW m^-2 sr^-1 per unit z, "
             "value ± jackknife σ:\n")
    L.append("| slice | " + " | ".join(f"{x:.3f}" for x in lamc) + " |")
    L.append("|---|" + "---|" * len(lamc))
    for s, (t, lo, hi) in enumerate(SLICES):
        L.append(f"| {t[:3]} {lo:.2f}-{hi:.2f} | " + " | ".join(f"{val[s,g]:.0f}±{err[s,g]:.0f}" for g in range(len(lamc))) + " |")
    # nulls
    L.append("\nNulls: χ² of each null over the 32 slices per group against zero with jackknife σ "
             "(dof 32; p from χ²):\n")
    from scipy.stats import chi2
    L.append("| group µm | signal χ² | shuffled χ² | stars χ² | zodi χ² | half-diff χ² | half-diff σ / jackknife σ |")
    L.append("|---|---|---|---|---|---|---|")
    stats = {}
    for g in range(len(lamc)):
        row = []
        for key in ("bdidz_grp", "shuf", "star", "zodi", "half"):
            a = J[key + "_jk"][:, :, g]
            v = J[key][:, g]
            e = err_of(a)
            ok = np.isfinite(v) & np.isfinite(e) & (e > 0)
            c2 = float(((v[ok] / e[ok]) ** 2).sum())
            row.append((c2, int(ok.sum())))
            stats[(g, key)] = (c2, int(ok.sum()))
        # half-difference scatter vs jackknife error of the signal
        h = J["half"][:, g]
        e_sig = err_of(J["bdidz_grp_jk"][:, :, g])
        ok = np.isfinite(h) & (e_sig > 0)
        ratio = float(np.sqrt(np.mean((h[ok] / e_sig[ok]) ** 2)))
        L.append(f"| {lamc[g]:.3f} | " + " | ".join(f"{c:.0f}/{n} (p={chi2.sf(c, n):.1e})" for c, n in row) + f" | {ratio:.2f} |")
    # shuffled null, as amended: per tracer the shuffled slices must be flat in z (chi2 about the tracer's
    # inverse-variance mean), and the line model fitted to them must return zero line amplitudes
    L.append("\nShuffled-redshift null (amended): χ² of the shuffled b dI/dz about each tracer's own mean "
             "(dof = slices - 4), per group:\n")
    L.append("| group µm | shuffled flatness χ²/dof (p) |")
    L.append("|---|---|")
    tr_of = np.array([t for (t, _, _) in SLICES])
    shv = J["shuf"] * tonw
    shj = J["shuf_jk"] * tonw
    she = err_of(shj)
    for g in range(len(lamc)):
        c2 = 0.0
        dof = 0
        for t in TRACERS:
            m = (tr_of == t) & np.isfinite(shv[:, g]) & (she[:, g] > 0)
            if m.sum() < 2:
                continue
            w = 1 / she[m, g] ** 2
            mu = (shv[m, g] * w).sum() / w.sum()
            c2 += float((((shv[m, g] - mu) / she[m, g]) ** 2).sum())
            dof += int(m.sum()) - 1
        L.append(f"| {lamc[g]:.3f} | {c2:.0f}/{dof} (p={chi2.sf(c2, dof):.1e}) |")
    # line model per group
    L.append("\nLine model per group: continuum quadratic in z + line amplitudes × the fraction of each slice whose line "
             "λ_obs falls in the group; amplitudes in nW m^-2 sr^-1 per unit z, σ from jackknife refits:\n")
    L.append("| group µm | " + " | ".join(f"{k} (A ± σ, S/N)" for k in LINES) + " | continuum χ²/dof |")
    L.append("|---|---|---|---|---|")
    lines_out = {}
    for g, (a, b) in enumerate(ge):
        F = np.array([[line_frac(lo, hi, lr, a, b) for (_, lo, hi) in SLICES] for lr in LINES.values()]).T
        use = [k for k in range(F.shape[1]) if F[:, k].max() > 0]
        X = np.column_stack([np.ones(NS), zc, zc ** 2] + [F[:, k] for k in use])
        y = val[:, g]
        e = err[:, g]
        ok = np.isfinite(y) & (e > 0)

        def fit(yy):
            w = 1 / e[ok]
            beta, *_ = np.linalg.lstsq(X[ok] * w[:, None], yy[ok] * w, rcond=None)
            return beta

        beta = fit(y)
        bj = np.stack([fit(jk[i, :, g]) for i in range(npch)])
        sb = err_of(bj)
        ok_s = ok & np.isfinite(shv[:, g])
        # shuffled galaxies carry their tracer's projected signal: fit with one offset per tracer
        Xs = np.column_stack([X] + [(tr_of == t).astype(float) for t in TRACERS[1:]])

        def fit_s(yy):
            w = 1 / e[ok_s]
            beta, *_ = np.linalg.lstsq(Xs[ok_s] * w[:, None], yy[ok_s] * w, rcond=None)
            return beta

        bsh = fit_s(shv[:, g])
        bshj = np.stack([fit_s(shj[i, :, g]) for i in range(npch)])
        ssh = err_of(bshj)
        for k, name in enumerate(LINES):
            if k in use:
                j = 3 + use.index(k)
                lines_out[(g, name + "_shuf")] = (bsh[j], ssh[j])
        resid = (y[ok] - X[ok] @ beta) / e[ok]
        cells = []
        for k, name in enumerate(LINES):
            if k in use:
                j = 3 + use.index(k)
                cells.append(f"{beta[j]:.0f} ± {sb[j]:.0f} ({beta[j]/sb[j]:.1f})")
                lines_out[(g, name)] = (beta[j], sb[j])
            else:
                cells.append("—")
        L.append(f"| {lamc[g]:.3f} | " + " | ".join(cells) + f" | {(resid**2).sum():.0f}/{ok.sum()-X.shape[1]} |")
    # joint line significance
    for name in list(LINES) + [n + "_shuf" for n in LINES]:
        v = [lines_out[k] for k in lines_out if k[1] == name]
        if v:
            c2 = sum((x / s) ** 2 for x, s in v)
            mean = sum(x / s ** 2 for x, s in v) / sum(1 / s ** 2 for x, s in v)
            L.append(f"\n{name}: {len(v)} groups, χ² against zero {c2:.1f} (p={chi2.sf(c2, len(v)):.2e}); "
                     f"inverse-variance mean amplitude {mean:.0f} ± {1/np.sqrt(sum(1/s**2 for x, s in v)):.0f}.")
    # headline: inverse-variance mean over slices per group, and over slices in broad z ranges
    L.append("\nSummary per group: inverse-variance mean of b_I d(νI_ν)/dz over slices (σ from the jackknife of the "
             "mean), all slices, z < 0.5 and z > 0.5; and the stellar-density split (patches below/above the median "
             "Gaia 12<G<16 count; added after the interim look):\n")
    L.append("| group µm | all z | z < 0.5 | z ≥ 0.5 | low-star half | high-star half |")
    L.append("|---|---|---|---|---|---|")
    sp = os.path.join(SP, "result_split.npz")
    SPL = np.load(sp) if os.path.exists(sp) else None

    def ivm(v, jj, m):
        e_ = err_of(jj)
        mm = m & np.isfinite(v) & (e_ > 0)
        w = 1 / e_[mm] ** 2
        return (v[mm] * w).sum() / w.sum(), (jj[:, mm] * w).sum(1) / w.sum()

    def ivm_split(v, jj, m):
        n_ = jj.shape[0]
        e_ = np.sqrt((n_ - 1) / n_ * ((jj - jj.mean(0)) ** 2).sum(0))
        mm = m & np.isfinite(v) & (e_ > 0)
        w = 1 / e_[mm] ** 2
        mj = (jj[:, mm] * w).sum(1) / w.sum()
        return (v[mm] * w).sum() / w.sum(), np.sqrt((n_ - 1) / n_ * ((mj - mj.mean()) ** 2).sum())

    summary = {}
    for g in range(len(lamc)):
        cells = []
        for m in (np.ones(NS, bool), zc < 0.5, zc >= 0.5):
            mu, mj = ivm(val[:, g], jk[:, :, g], m)
            cells.append(f"{mu:.1f} ± {err_of(mj):.1f}")
            summary.setdefault(g, []).append((mu, err_of(mj)))
        if SPL is not None:
            for k in ("lowstar", "highstar"):
                mu, se = ivm_split(SPL[k][:, g] * tonw[g], SPL[k + "_jk"][:, :, g] * tonw[g], np.ones(NS, bool))
                cells.append(f"{mu:.1f} ± {se:.1f}")
        L.append(f"| {lamc[g]:.3f} | " + " | ".join(cells) + " |")
    json.dump({f"{lamc[g]:.3f}": summary[g] for g in summary}, open(os.path.join(SP, "summary.json"), "w"))
    # CIBER comparison
    L.append("\nCIBER comparison (arXiv:2608.12116, CIBER × DESI-LS, read off its Fig. a2h_vs_redshift); ours: inverse-"
             "variance mean of the slices whose centre lies in the CIBER bin, σ from jackknife of that mean:\n")
    L.append("| λ | z | CIBER | SPHEREx group (µm) | this work |")
    L.append("|---|---|---|---|---|")
    for lam_ciber, rows in CIBER.items():
        gsel = [g for g, (a, b) in enumerate(ge) if a <= lam_ciber < b]
        for g in gsel:
            for (lo, hi, v, e) in rows:
                m = (zc >= lo) & (zc < hi) & np.isfinite(val[:, g]) & (err[:, g] > 0)
                if not m.any():
                    continue
                wts = 1 / err[m, g] ** 2
                mean = (val[m, g] * wts).sum() / wts.sum()
                mj = (jk[:, m, g] * wts).sum(1) / wts.sum()
                cib = f"{v} ± {e}" if v is not None else f"< {e} (95%)"
                L.append(f"| {lam_ciber} | {lo}-{hi} | {cib} | {lamc[g]:.3f} | {mean:.0f} ± {err_of(mj):.0f} |")
    txt = "\n".join(L)
    open(os.path.join(SP, "report.md"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    {"maps": stage_maps, "pairs": stage_pairs, "result": stage_result, "report": stage_report}[sys.argv[1]]()
