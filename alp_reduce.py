"""R-153: reduce one SPHEREx QR2 L2 image to sky-cell x wavelength sums around a dwarf (or a fake dwarf).

For each image (one detector, one exposure) it
  1. reads IMAGE, FLAGS, VARIANCE, ZODI (MJy/sr) and the per-pixel wavelength and bandwidth (the tabular 'W' WCS;
     bandwidth = lambda / R, the FWHM of the pixel's LVF passband);
  2. masks pixels with a bad flag or the SOURCE flag (bit 21, the pipeline's catalogue-source mask), then pixels more
     than 5 robust sigma above their channel-strip median, dilated by 2 pixels (as spherex_reduce.py, R-138);
  3. takes the residual r = IMAGE - ZODI - the fixed detector pattern (alp_pattern.py, built without this dwarf's own
     images; optional) (no strip median: ON and OFF are differenced later at the same wavelength in
     the same exposure, which removes any per-exposure offset);
  4. bins every good pixel within 2 deg of the centre by (sky cell, ln(lambda) bin of width 0.002) and returns, for each
     occupied (cell, bin): N, sum r, sum VARIANCE, sum ZODI, sum bandwidth, sum lambda. Sky cells are 0.05 x 0.05 deg
     squares in the gnomonic plane about the centre (80 x 80 over +-2 deg), so a fixed sky position keeps its cell in
     every exposure while its wavelength changes; that is what lets the analysis fit each cell's own continuum.

ln(lambda) grid: ln(0.74 um) + 0.002 k, k < 958.
Usage: python -u alp_reduce.py [--pattern <dir> <exclude-group> mean|lin] <list.csv> <outdir> <ra> <dec> <worker_index> <n_workers>
list.csv: obs_id,detector,t_mid,url[,sep] ; output <outdir>/<obs_id>_D<det>.npz (existing ones are skipped).
"""
import io
import os
import sys
import time
import urllib.request
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import binary_dilation

warnings.filterwarnings("ignore")

CELL, NC, RMAX = 0.05, 80, 2.0
L0, DL, NB = np.log(0.74), 0.002, 958
BAD_BITS = [0, 1, 2, 4, 5, 6, 7, 9, 10, 11, 15, 17, 19]
SOURCE_BIT = 21


def geometry(h, ra0, de0, step=8):
    """Per-pixel angular distance (deg) and position-angle quadrant from (ra0, de0), and wavelength, bandwidth (um).

    The WCSs are evaluated on a grid every `step` pixels and interpolated bilinearly (they are smooth: lambda changes
    by ~2e-4 in ln per pixel, against the 0.002 bin), with the sky positions as gnomonic offsets about the centre."""
    n = h["IMAGE"].data.shape[0]
    g = np.arange(0, n + step, step).clip(0, n - 1)
    yc, xc = np.meshgrid(g, g, indexing="ij")
    ww = WCS(h["IMAGE"].header, fobj=h, key="W")
    lam, bw = (a.reshape(xc.shape) for a in ww.pixel_to_world_values(xc.ravel(), yc.ravel()))
    wc = WCS(h["IMAGE"].header)
    ra, de = (np.radians(a).reshape(xc.shape) for a in wc.all_pix2world(xc.ravel(), yc.ravel(), 0))
    r0, d0 = np.radians(ra0), np.radians(de0)
    cosc = np.sin(de) * np.sin(d0) + np.cos(de) * np.cos(d0) * np.cos(ra - r0)
    xi = np.cos(de) * np.sin(ra - r0) / cosc
    eta = (np.cos(d0) * np.sin(de) - np.sin(d0) * np.cos(de) * np.cos(ra - r0)) / cosc
    p = np.arange(n) / step
    i0 = np.minimum(p.astype(int), len(g) - 2)
    f = (p - i0).astype(np.float32)

    def up(a):  # separable bilinear upsampling of the coarse grid (i0 = p // step, so rows repeat `step` times)
        a = a.astype(np.float64)
        b = (a[:, i0] * (1 - f) + a[:, i0 + 1] * f).astype(np.float32)
        lo = np.repeat(b[:-1], step, axis=0)[:n]
        hi = np.repeat(b[1:], step, axis=0)[:n]
        return lo * (1 - f[:, None]) + hi * f[:, None]

    xi, eta, lam, bw = up(xi), up(eta), up(lam), up(bw)
    th = np.degrees(np.arctan(np.hypot(xi, eta)))
    cx = np.floor(np.degrees(xi) / CELL + NC / 2).astype(int)
    cy = np.floor(np.degrees(eta) / CELL + NC / 2).astype(int)
    inside = (cx >= 0) & (cx < NC) & (cy >= 0) & (cy < NC) & (th < RMAX)
    cell = np.where(inside, cy * NC + cx, -1)
    return th, cell, lam, bw


PATTERN = {}


def load_pattern(pdir, exclude, model="mean"):
    """Per detector, the fixed residual pattern (nu I_nu nW) from alp_pattern.py sums over every group but `exclude`.
    model 'mean': the mean residual per pixel; 'lin': a + b Z per pixel (Z the image's median zodi), b regularised
    toward zero where the groups' zodi levels do not separate it."""
    import glob
    out = {}
    for det in range(1, 7):
        acc = None
        for f in glob.glob(os.path.join(pdir, f"*.D{det}.npz")):
            if os.path.basename(f).split(".")[0] == exclude:
                continue
            z = np.load(f)
            if acc is None:
                acc = {k: z[k].astype(np.float64) for k in z.files}
            else:
                for k in z.files:
                    acc[k] += z[k]
        if acc is None:
            continue
        n = np.maximum(acc["n"], 1)
        if model == "lin":
            zm = acc["z"] / n
            vz = acc["zz"] / n - zm ** 2
            cov = acc["sz"] / n - zm * acc["s"] / n
            b = cov / (vz + 0.05 * zm ** 2 + 1e-9)
            p = acc["s"] / n + b * 0.0  # stored as (a at mean Z, b): applied as P(Z) = mean + b (Z - zm)
            out[det] = (np.where(acc["n"] >= 10, p, 0.0), np.where(acc["n"] >= 10, b, 0.0), zm)
        else:
            out[det] = (np.where(acc["n"] >= 10, acc["s"] / n, 0.0), None, None)
    return out


def reduce_file(data, ra0, de0):
    h = fits.open(io.BytesIO(data), memmap=False)
    det = int(h["IMAGE"].header["DETECTOR"])
    im = h["IMAGE"].data.astype(np.float64)
    fl = h["FLAGS"].data
    zo = h["ZODI"].data.astype(np.float64)
    var = h["VARIANCE"].data.astype(np.float64)
    bad = np.zeros(fl.shape, bool)
    for b in BAD_BITS:
        bad |= ((fl >> b) & 1).astype(bool)
    src = ((fl >> SOURCE_BIT) & 1).astype(bool)
    good = ~bad & ~src & np.isfinite(im) & np.isfinite(zo) & np.isfinite(var) & (var > 0)
    th, cell, lam, bw = geometry(h, ra0, de0)
    n = im.shape[0]
    res = im - zo
    pscale = np.nan
    if det in PATTERN:
        # the pattern scaled per image: s = <o P> / <P P> over 8x8 blocks of good pixels, o the image's residual with
        # each ln(lambda) bin's median removed (its sky spectrum), as the pattern itself was built
        p = PATTERN[det][0]
        conv = 1e-20 * 2.99792458e14 / lam * 1e9
        lbp = np.floor((np.log(lam) - L0) / DL).astype(int)
        rn = np.where(good, res * conv, np.nan)
        okp = np.isfinite(rn)
        o = np.full(rn.shape, np.nan)
        kk = lbp[okp]
        vv = rn[okp]
        srt = np.argsort(kk, kind="stable")
        uk_, st_ = np.unique(kk[srt], return_index=True)
        en_ = np.r_[st_[1:], len(srt)]
        med = np.empty(len(vv))
        for a_, b_ in zip(st_, en_):
            med[srt[a_:b_]] = np.median(vv[srt[a_:b_]])
        o[okp] = vv - med
        sg = 1.4826 * np.nanmedian(np.abs(o))
        o[np.abs(o) > 5 * sg] = np.nan
        m8 = (n // 8) * 8
        ob = o[:m8, :m8].reshape(m8 // 8, 8, m8 // 8, 8)
        pb = np.where(np.isfinite(ob), p[:m8, :m8].reshape(m8 // 8, 8, m8 // 8, 8), np.nan)
        ob, pb = np.nanmean(ob, axis=(1, 3)), np.nanmean(pb, axis=(1, 3))
        okb = np.isfinite(ob) & np.isfinite(pb)
        pscale = float((ob[okb] * pb[okb]).sum() / (pb[okb] ** 2).sum())
        res = res - pscale * p / conv
    # hot-pixel clip per coarse strip (17 per detector, as R-138)
    lb = np.floor((np.log(lam) - L0) / DL).astype(int)
    strip = lb // 12
    hot = np.zeros(good.shape, bool)
    for s in np.unique(strip[good]):
        sel = good & (strip == s)
        if sel.sum() < 2000:
            continue
        v = res[sel]
        med = np.median(v)
        sig = 1.4826 * np.median(np.abs(v - med))
        hot |= sel & (res > med + 5 * sig)
    good &= ~binary_dilation(hot, structure=np.ones((5, 5), bool))
    good &= (cell >= 0) & (lb >= 0) & (lb < NB)
    key = cell[good].astype(np.int64) * NB + lb[good]
    uk, inv = np.unique(key, return_inverse=True)

    def agg(a):
        return np.bincount(inv, weights=a, minlength=len(uk))

    r = res[good]
    out = dict(cell=(uk // NB).astype(np.int16), lb=(uk % NB).astype(np.int16), n=agg(np.ones_like(r)).astype(np.int32),
               s=agg(r), v=agg(var[good]).astype(np.float32), z=agg(zo[good]).astype(np.float32),
               bw=agg(bw[good]).astype(np.float32), lam=agg(lam[good]).astype(np.float32))
    out["pscale"] = np.array(pscale)
    return out, det, float(good.mean())


def fetch(url, tries=5):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return r.read()
        except Exception:  # network: retry with backoff
            if k == tries - 1:
                raise
            time.sleep(2 ** (k + 1))


def main():
    args = sys.argv[1:]
    if args and args[0] == "--pattern":
        # --pattern <dir> <exclude-group> <model>
        PATTERN.update(load_pattern(args[1], args[2], args[3]))
        print("pattern detectors", sorted(PATTERN), "excluding", args[2], "model", args[3], flush=True)
        args = args[4:]
    lst, outdir, ra0, de0, wi, nw = args[0], args[1], float(args[2]), float(args[3]), int(args[4]), int(args[5])
    os.makedirs(outdir, exist_ok=True)
    rows = [l.strip().split(",") for l in open(lst) if l.strip() and not l.startswith("obs_id")]
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
            d, dd, fg = reduce_file(fetch(url), ra0, de0)
        except Exception as e:
            print("FAIL", obs, det, repr(e)[:200], flush=True)
            continue
        np.savez_compressed(out + ".tmp.npz", t=float(tmid), det=dd, frac_good=fg, **d)
        os.replace(out + ".tmp.npz", out)
        done += 1
        if done % 10 == 0:
            print(f"w{wi} {done} done, {(time.time()-t0)/done:.2f} s/img", flush=True)
    print(f"w{wi} finished {done}", flush=True)


if __name__ == "__main__":
    main()
