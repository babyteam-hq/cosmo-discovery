"""R-143: the source-mask transfer by injection. Reduce one SPHEREx QR2 L2 image through R-138's chain (spherex_reduce.py:
FLAGS bad bits and the SOURCE bit, 5-robust-sigma clip per channel strip dilated 2 px, ZODI subtracted, strip median of
good pixels subtracted, 8x8 blocks to NSIDE-2048 HEALPix), once on the real image and again with mock galaxy-correlated
light painted into it, and keep what the chain lets through of the mock.

Mock (fixed in results-R138.md, "R-143: injection test", before any transfer was computed). Every DESI DR1 NGC galaxy of
the four tracers inside its tracer's slice range (BGS_ANY 0.05-0.45, LRG 0.4-1.1, ELG_LOPnotqso 0.8-1.6, QSO 0.8-2.1) is
given a flux F_g = WEIGHT_g x A(z) dz / (K(z) n_ts), n_ts the tracer-slice's weighted surface density in the field and K(z)
the number of tracers covering z, so the mock's mean d(nuI_nu)/dz equals A(z) and every tracer carries an equal share.
Its light is split into two components:
  compact : a Gaussian of sigma = 1 native pixel (6.15") at the galaxy, i.e. the galaxy's own light;
  smooth  : a Gaussian of sigma = 1 comoving Mpc at the galaxy's redshift, an unresolved component correlated with it
            (deposited on the 49" block grid with a margin, convolved by FFT in 9 redshift bins, constant within a block).
The mock is grey in nuI_nu (I_nu = nuI_nu / nu per channel). Variants:
  A : A(z) = 6 nW m^-2 sr^-1 per unit z (the R-138 level), half compact, half smooth (fiducial);
  B : A x 5 (the clip's non-linearity);
  C : the smooth half of A alone;  D : the compact half of A alone;
  E : A(z) = 6 ((1+z)/1.8)^3, half/half (a rising redshift shape).
Per image the chain runs on the real image (its own clip: R-138's output), on real + A (A's own clip) and on real + B (B's
own clip); C, D and E go through A's mask with their own strip medians (the clip at A's level is the realistic one, and
C + D = A on that mask exactly). Kept per (HEALPix pixel, channel), as block sums:
  s0, n0       real residual and good-pixel count (R-138's chain);       s0sm, n0sm  the same with the star mask below
  dA, nA       (real + A) residual minus real residual, both on A's mask, and A's good-pixel count (dAsm, nAsm: star mask)
  dB, nB       the same for B on B's mask
  dC, dD, dE   the same for C, D, E on A's mask
  mC, mD, mE   the mock on A's good pixels, no strip median (the mask alone)
  tC, tD, tE   the mock summed over all 64 pixels of the block (truth: no mask, no clip, no filter), nt = 64 per block
Star mask (R-143's fixed-flux mask, for the stellar null): every block whose centre lies within r(G) + 12" of a Gaia DR3
star with G < 16, r(G) = 1' x 10^(0.2 (12 - G)) (first tried 1.5' and + 25": 49 % of blocks in the test image, too many).

Usage: python -u spherex_inject.py <list.csv> <outdir> <worker_index> <n_workers> [max_images]
"""
import glob
import io
import os
import sys
import time
import warnings

import numpy as np
import healpy as hp
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import binary_dilation

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spherex_reduce import BAD_BITS, SOURCE_BIT, BLOCK, NSIDE, channel_of, fetch, EDGES, NCH  # noqa: E402

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")
PIX_AS = 6.15  # native pixel, arcsec
TRACERS = {"BGS_ANY": (0.05, 0.45), "LRG": (0.4, 1.1), "ELG_LOPnotqso": (0.8, 1.6), "QSO": (0.8, 2.1)}
A0 = 6.0  # nW m^-2 sr^-1 per unit z
ZB = np.array([0.05, 0.15, 0.25, 0.35, 0.5, 0.7, 0.9, 1.2, 1.6, 2.1])
IDX_NSIDE = 16
KEYS = ["s0", "n0", "s0sm", "n0sm", "dA", "nA", "dAsm", "nAsm", "dB", "nB", "dC", "dD", "dE", "mC", "mD", "mE",
        "tC", "tD", "tE", "nt"]


def shape_E(z):
    return ((1 + z) / 1.8) ** 3


def load_galaxies():
    """All DESI galaxies inside their tracer range, with flat (A) and rising (E) fluxes in nW m^-2 sr^-1 x sr."""
    from astropy.cosmology import Planck18
    reg = np.load(os.path.join(D, "spherex", "region32.npy"))
    inreg = np.zeros(hp.nside2npix(32), bool)
    inreg[reg] = True
    omega = len(reg) * hp.nside2pixarea(32)
    ra, dec, z, fa, fe = [], [], [], [], []
    for t, (lo, hi) in TRACERS.items():
        g = np.load(os.path.join(D, "desi_dr1", f"{t}_dat.npz"))
        m = (g["z"] >= lo) & (g["z"] < hi)
        gra, gdec, gz, gw = g["ra"][m], g["dec"][m], g["z"][m].astype(float), g["w"][m].astype(float)
        core = inreg[hp.ang2pix(32, gra, gdec, lonlat=True)]
        sl = np.floor((gz - lo) / 0.1 + 1e-9).astype(int)
        nsl = int(round((hi - lo) / 0.1))
        sl = np.clip(sl, 0, nsl - 1)
        n_ts = np.bincount(sl[core], weights=gw[core], minlength=nsl) / omega  # per sr
        K = np.zeros_like(gz)
        for (lo2, hi2) in TRACERS.values():
            K += (gz >= lo2) & (gz < hi2)
        base = gw * 0.1 / (K * n_ts[sl])
        ra.append(gra); dec.append(gdec); z.append(gz)
        fa.append(base * A0); fe.append(base * A0 * shape_E(gz))
    ra, dec, z, fa, fe = map(np.concatenate, (ra, dec, z, fa, fe))
    chi = Planck18.comoving_distance(z).value
    sig_as = np.degrees(1.0 / chi) * 3600  # 1 cMpc in arcsec
    p = hp.ang2pix(IDX_NSIDE, ra, dec, lonlat=True, nest=True)
    o = np.argsort(p)
    cat = dict(ra=ra[o], dec=dec[o], z=z[o], fa=fa[o], fe=fe[o], sig=sig_as[o])
    cat["start"] = np.searchsorted(p[o], np.arange(hp.nside2npix(IDX_NSIDE) + 1))
    return cat


def load_stars():
    fs = sorted(glob.glob(os.path.join(D, "gaia", "g_*.csv")))
    out = []
    for f in fs:
        ra0, dec0 = [float(x) for x in os.path.basename(f)[2:-4].split("_")]
        a = np.loadtxt(f, delimiter=",", skiprows=1, ndmin=2)
        if a.size == 0:
            continue
        keep = (a[:, 0] >= ra0) & (a[:, 0] < ra0 + 10) & (a[:, 1] >= dec0) & (a[:, 1] < dec0 + 10)
        out.append(a[keep])
    a = np.concatenate(out)
    p = hp.ang2pix(IDX_NSIDE, a[:, 0], a[:, 1], lonlat=True, nest=True)
    o = np.argsort(p)
    return dict(ra=a[o, 0], dec=a[o, 1], g=a[o, 2],
                start=np.searchsorted(p[o], np.arange(hp.nside2npix(IDX_NSIDE) + 1)))


def near(cat, ra0, dec0, rad_deg):
    pix = hp.query_disc(IDX_NSIDE, hp.ang2vec(ra0, dec0, lonlat=True), np.radians(rad_deg), nest=True, inclusive=True)
    idx = np.concatenate([np.arange(cat["start"][p], cat["start"][p + 1]) for p in pix])
    return idx


def strip_index(ch_p, chans):
    flat = ch_p.ravel()
    return {c: np.flatnonzero(flat == c) for c in chans}


def clip_mask(res, good, sidx):
    """R-138's clip: pixels > 5 robust sigma above the strip median of good pixels, dilated by 2 px."""
    hot = np.zeros(res.size, bool)
    r = res.ravel()
    g = good.ravel()
    for c, ii in sidx.items():
        sel = ii[g[ii]]
        if sel.size < 2000:
            continue
        v = r[sel]
        med = np.median(v)
        sig = 1.4826 * np.median(np.abs(v - med))
        hot[sel[v > med + 5 * sig]] = True
    hot = binary_dilation(hot.reshape(res.shape), structure=np.ones((5, 5), bool))
    return good & ~hot


def strip_medians(res, good, sidx):
    """Per strip median of good pixels; strips with < 2000 good pixels are dropped (R-138). Returns (med dict, good)."""
    g = good.ravel().copy()
    r = res.ravel()
    med = {}
    for c, ii in sidx.items():
        sel = ii[g[ii]]
        if sel.size < 2000:
            g[ii] = False
            continue
        med[c] = float(np.median(r[sel]))
    return med, g.reshape(good.shape)


def medians_on(res, good, sidx, chans):
    r = res.ravel()
    g = good.ravel()
    return {c: float(np.median(r[sidx[c][g[sidx[c]]]])) for c in chans}


def paint(cat, wc, n, nb, lam_ch_of_block):
    """Mock nuI_nu components (nW m^-2 sr^-1) on the pixel grid: compact_A, compact_E; on the block grid: smooth_A, smooth_E."""
    ra0, dec0 = wc.all_pix2world([[n / 2, n / 2]], 0)[0]
    margin = 0.6  # deg: the smooth grid's margin (45 blocks), > 3 sigma of 1 cMpc at z ~ 0.1
    idx = near(cat, ra0, dec0, n * PIX_AS / 3600 * 0.75 + margin)
    x, y = wc.all_world2pix(cat["ra"][idx], cat["dec"][idx], 0, quiet=True)
    pix_sr = np.radians(PIX_AS / 3600) ** 2
    blk_sr = pix_sr * BLOCK * BLOCK
    out = {}
    # compact: sigma 1 px Gaussian, 7x7 stamp, normalised
    inn = (x > -4) & (x < n + 3) & (y > -4) & (y < n + 3)
    xi, yi = np.round(x[inn]).astype(int), np.round(y[inn]).astype(int)
    dx, dy = x[inn] - xi, y[inn] - yi
    off = np.arange(-3, 4)
    for key, f in (("cA", cat["fa"][idx][inn]), ("cE", cat["fe"][idx][inn])):
        img = np.zeros(n * n)
        acc_w = []
        acc_i = []
        norm = np.zeros(len(xi))
        for ox in off:
            for oy in off:
                w = np.exp(-0.5 * ((ox - dx) ** 2 + (oy - dy) ** 2))
                norm += w
        for ox in off:
            for oy in off:
                w = np.exp(-0.5 * ((ox - dx) ** 2 + (oy - dy) ** 2)) / norm
                px, py = xi + ox, yi + oy
                ok = (px >= 0) & (px < n) & (py >= 0) & (py < n)
                acc_i.append(py[ok] * n + px[ok])
                acc_w.append((f * w)[ok] / pix_sr)
        img = np.bincount(np.concatenate(acc_i), weights=np.concatenate(acc_w), minlength=n * n)
        out[key] = img.reshape(n, n).astype(np.float32)
    # smooth: block grid with margin, per z bin Gaussian convolution by FFT
    pad = 45
    ng = nb + 2 * pad
    nf = 1 << int(np.ceil(np.log2(ng + 2 * pad)))
    bx = (x + 0.5) / BLOCK - 0.5 + pad
    by = (y + 0.5) / BLOCK - 0.5 + pad
    ok = (bx > -0.5) & (bx < ng - 0.5) & (by > -0.5) & (by < ng - 0.5)
    ib, jb = np.round(bx[ok]).astype(int), np.round(by[ok]).astype(int)
    zz = cat["z"][idx][ok]
    kx = np.fft.fftfreq(nf)
    KX, KY = np.meshgrid(kx, np.fft.rfftfreq(nf), indexing="ij")
    k2 = KX ** 2 + KY ** 2
    for key, f in (("hA", cat["fa"][idx][ok]), ("hE", cat["fe"][idx][ok])):
        tot = np.zeros((nf, nf))
        for b in range(len(ZB) - 1):
            m = (zz >= ZB[b]) & (zz < ZB[b + 1])
            if not m.any():
                continue
            sig_b = np.median(cat["sig"][idx][ok][m]) / (PIX_AS * BLOCK)
            P = np.zeros((nf, nf))
            np.add.at(P, (jb[m], ib[m]), f[m] / blk_sr)
            tot += np.fft.irfft2(np.fft.rfft2(P) * np.exp(-2 * np.pi ** 2 * sig_b ** 2 * k2.T.T), s=(nf, nf))
        out[key] = tot[pad:pad + nb, pad:pad + nb].astype(np.float32)
    return out


def star_blocks(stars, wc, n, nb):
    ra0, dec0 = wc.all_pix2world([[n / 2, n / 2]], 0)[0]
    idx = near(stars, ra0, dec0, n * PIX_AS / 3600 * 0.75 + 0.3)
    g = stars["g"][idx]
    x, y = wc.all_world2pix(stars["ra"][idx], stars["dec"][idx], 0, quiet=True)
    r_blk = (60 * 10 ** (0.2 * (12 - g)) + 12) / (PIX_AS * BLOCK)
    bx = (x + 0.5) / BLOCK - 0.5
    by = (y + 0.5) / BLOCK - 0.5
    m = np.zeros((nb, nb), bool)
    ok = (bx > -r_blk) & (bx < nb + r_blk) & (by > -r_blk) & (by < nb + r_blk)
    J, I = np.mgrid[0:nb, 0:nb]
    for cx, cy, rr in zip(bx[ok], by[ok], r_blk[ok]):
        i0, i1 = max(int(cx - rr), 0), min(int(cx + rr) + 2, nb)
        j0, j1 = max(int(cy - rr), 0), min(int(cy + rr) + 2, nb)
        if i0 >= i1 or j0 >= j1:
            continue
        m[j0:j1, i0:i1] |= (I[j0:j1, i0:i1] - cx) ** 2 + (J[j0:j1, i0:i1] - cy) ** 2 <= rr ** 2
    return m


def reduce_inject(data, cat, stars):
    h = fits.open(io.BytesIO(data), memmap=False)
    det = int(h["IMAGE"].header["DETECTOR"])
    im = h["IMAGE"].data.astype(np.float32)
    fl = h["FLAGS"].data
    zo = h["ZODI"].data.astype(np.float32)
    bad = np.zeros(fl.shape, bool)
    for b in BAD_BITS:
        bad |= ((fl >> b) & 1).astype(bool)
    src = ((fl >> SOURCE_BIT) & 1).astype(bool)
    good0 = ~bad & ~src & np.isfinite(im) & np.isfinite(zo)
    ww = WCS(h["IMAGE"].header, fobj=h, key="W")
    n = im.shape[0]
    nb = n // BLOCK
    yc, xc = np.mgrid[0:nb, 0:nb]
    xc = xc * BLOCK + (BLOCK - 1) / 2.0
    yc = yc * BLOCK + (BLOCK - 1) / 2.0
    lam_b, _ = ww.pixel_to_world_values(xc.ravel(), yc.ravel())
    lam_b = lam_b.reshape(nb, nb)
    ch_b = channel_of(det, lam_b)
    ch_p = np.repeat(np.repeat(ch_b, BLOCK, 0), BLOCK, 1)
    ch_p = np.pad(ch_p, ((0, n - ch_p.shape[0]), (0, n - ch_p.shape[1])), mode="edge")
    chans = [int(c) for c in np.unique(ch_b[ch_b >= 0])]
    sidx = strip_index(ch_p, chans)
    res = np.where(np.isfinite(im) & np.isfinite(zo), im - zo, 0.0).astype(np.float32)
    wc = WCS(h["IMAGE"].header)
    # mock, nuI_nu -> I_nu (MJy/sr) per channel: I_nu = nuI_nu / (nu 1e-11)
    lo, hi = EDGES[det]
    e = lo * (hi / lo) ** (np.arange(NCH + 1) / NCH)
    lamc = np.sqrt(e[:-1] * e[1:])
    conv_b = np.where(ch_b >= 0, 1.0 / (2.99792458e14 / lamc[np.clip(ch_b - (det - 1) * NCH, 0, NCH - 1)] * 1e-11), 0.0)
    conv_p = np.repeat(np.repeat(conv_b, BLOCK, 0), BLOCK, 1)
    conv_p = np.pad(conv_p, ((0, n - conv_p.shape[0]), (0, n - conv_p.shape[1])), mode="edge")
    P = paint(cat, wc, n, nb, None)
    up = lambda b: np.pad(np.repeat(np.repeat(b, BLOCK, 0), BLOCK, 1), ((0, n - nb * BLOCK), (0, n - nb * BLOCK)), mode="edge")
    MC = (0.5 * up(P["hA"]) * conv_p).astype(np.float32)
    MD = (0.5 * P["cA"] * conv_p).astype(np.float32)
    ME = (0.5 * (up(P["hE"]) + P["cE"]) * conv_p).astype(np.float32)
    MA = MC + MD
    MB = 5 * MA
    # real chain (R-138)
    g0 = clip_mask(res, good0, sidx)
    med0, g0 = strip_medians(res, g0, sidx)
    # A: own clip on real + A
    gA = clip_mask(res + MA, good0, sidx)
    medA, gA = strip_medians(res + MA, gA, sidx)
    cA = list(medA)
    med0A = medians_on(res, gA, sidx, cA)
    gB = clip_mask(res + MB, good0, sidx)
    medB, gB = strip_medians(res + MB, gB, sidx)
    med0B = medians_on(res, gB, sidx, list(medB))
    medC = medians_on(res + MC, gA, sidx, cA)
    medD = medians_on(res + MD, gA, sidx, cA)
    medE = medians_on(res + ME, gA, sidx, cA)
    m = nb * BLOCK

    def bsum(a):
        return a[:m, :m].reshape(nb, BLOCK, nb, BLOCK).sum(axis=(1, 3))

    def per_block(md):  # strip value per block
        out = np.zeros((nb, nb))
        for c, v in md.items():
            out[ch_b == c] = v
        return out

    f0, fA, fB = g0.astype(np.float32), gA.astype(np.float32), gB.astype(np.float32)
    n0, nA, nB = bsum(f0), bsum(fA), bsum(fB)
    s0 = bsum(f0 * res) - per_block(med0) * n0
    mA_ = bsum(fA * MA)
    out = dict(s0=s0, n0=n0, nA=nA, nB=nB)
    out["dA"] = mA_ - (per_block(medA) - per_block(med0A)) * nA
    out["dB"] = bsum(fB * MB) - (per_block(medB) - per_block(med0B)) * nB
    for k, M, md in (("C", MC, medC), ("D", MD, medD), ("E", ME, medE)):
        mk = bsum(fA * M)
        out["m" + k] = mk
        out["d" + k] = mk - (per_block(md) - per_block(med0A)) * nA
        out["t" + k] = bsum(M)
    out["nt"] = np.full((nb, nb), float(BLOCK * BLOCK))
    sm = star_blocks(stars, wc, n, nb)
    keep_sm = ~sm
    out["s0sm"], out["n0sm"] = out["s0"] * keep_sm, out["n0"] * keep_sm
    out["dAsm"], out["nAsm"] = out["dA"] * keep_sm, out["nA"] * keep_sm
    ra, dec = wc.all_pix2world(xc.ravel(), yc.ravel(), 0)
    pix = hp.ang2pix(NSIDE, ra, dec, lonlat=True, nest=True)
    # R-138 drops blocks with <= 8 good pixels; each chain by its own count. The truth keeps every block.
    for ks, kn in ((("s0", "s0sm", "n0", "n0sm"), n0), (("dA", "dAsm", "nA", "nAsm", "dC", "dD", "dE", "mC", "mD", "mE"), nA),
                   (("dB", "nB"), nB)):
        for k in ks:
            out[k] = np.where(kn > 8, out[k], 0.0)
    keep = ch_b.ravel() >= 0
    key = pix[keep].astype(np.int64) * 128 + ch_b.ravel()[keep]
    uk, inv = np.unique(key, return_inverse=True)
    res_d = dict(hp=(uk // 128).astype(np.int32), ch=(uk % 128).astype(np.int8), det=det,
                 frac_good=float(g0.mean()), frac_sm=float(sm.mean()))
    for k in KEYS:
        res_d[k] = np.bincount(inv, weights=out[k].ravel()[keep], minlength=len(uk)).astype(np.float32)
    return res_d


def main():
    lst, outdir, wi, nw = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    nmax = int(sys.argv[5]) if len(sys.argv) > 5 else 10 ** 9
    os.makedirs(outdir, exist_ok=True)
    rows = [l.strip().split(",") for l in open(lst) if l.strip() and not l.startswith("obs_id")][:nmax]
    cat = load_galaxies()
    stars = load_stars()
    print(f"w{wi}: galaxies {len(cat['ra'])}, stars {len(stars['ra'])}", flush=True)
    t0 = time.time()
    done = 0
    for i, r in enumerate(rows):
        if i % nw != wi:
            continue
        obs, det, tmid, url = r[0], r[1], r[2], r[3]
        out = os.path.join(outdir, f"{obs}_D{det}.npz")
        if os.path.exists(out):
            continue
        try:
            d = reduce_inject(fetch(url), cat, stars)
        except Exception as e:
            print("FAIL", obs, det, repr(e)[:200], flush=True)
            continue
        np.savez_compressed(out + ".tmp.npz", t=float(tmid), **d)
        os.replace(out + ".tmp.npz", out)
        done += 1
        if done % 10 == 0:
            print(f"w{wi} {done} done, {(time.time()-t0)/done:.2f} s/img", flush=True)
    print(f"w{wi} finished {done}", flush=True)


if __name__ == "__main__":
    main()
