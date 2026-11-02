"""R-138: reduce one SPHEREx QR2 L2 spectral image to HEALPix-binned, zodi-subtracted, masked surface brightness.

For each image (one detector, one exposure) it
  1. reads IMAGE, FLAGS, ZODI (MJy/sr) and the per-pixel wavelength (the image's tabular 'W' WCS);
  2. masks pixels with any bad flag or the SOURCE flag (bit 21, the pipeline's known-source map), then pixels more than
     5 robust sigma above their channel-strip median (unflagged sources), dilated by 2 pixels;
  3. subtracts the ZODI model, then each channel strip's own median of good pixels (the per-exposure offset: zodi-model
     error, dark, diffuse light on scales larger than the strip), so only fluctuations on scales below ~3 deg survive;
  4. averages good pixels in 8x8 blocks (49"), assigns each block its HEALPix pixel (NSIDE 2048, nest) and its channel
     (17 log-uniform bins per detector between the band edges, fixed in results-R138.md);
  5. returns (hp, channel, sum of residual, number of good pixels), and the same for the ZODI model alone (the zodi null).

Usage: python -u spherex_reduce.py <list.csv> <outdir> <worker_index> <n_workers>
list.csv: obs_id,detector,t_mid,url (one row per image). Each output is <outdir>/<obs_id>_D<det>.npz; existing ones skip.
"""
import io
import os
import sys
import time
import urllib.request

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import healpy as hp
from scipy.ndimage import binary_dilation
import warnings

warnings.filterwarnings("ignore")

NSIDE = 2048
BLOCK = 8
NCH = 17
# band edges (um) from the obscore em_min/em_max of QR2 (D1..D6)
EDGES = {1: (0.7455, 1.1167), 2: (1.1031, 1.6476), 3: (1.6317, 2.4237), 4: (2.4070, 3.8255),
         5: (3.8020, 4.4234), 6: (4.4104, 5.0089)}
BAD_BITS = [0, 1, 2, 4, 5, 6, 7, 9, 10, 11, 15, 17, 19]
SOURCE_BIT = 21


def channel_of(det, lam):
    lo, hi = EDGES[det]
    x = np.log(lam / lo) / np.log(hi / lo) * NCH
    c = np.floor(x).astype(int)
    c = np.where((c < 0) | (c >= NCH), -1, c)
    return np.where(c < 0, -1, (det - 1) * NCH + c)


def reduce_file(data):
    h = fits.open(io.BytesIO(data), memmap=False)
    det = int(h["IMAGE"].header["DETECTOR"])
    im = h["IMAGE"].data.astype(np.float32)
    fl = h["FLAGS"].data
    zo = h["ZODI"].data.astype(np.float32)
    bad = np.zeros(fl.shape, bool)
    for b in BAD_BITS:
        bad |= ((fl >> b) & 1).astype(bool)
    src = ((fl >> SOURCE_BIT) & 1).astype(bool)
    good = ~bad & ~src & np.isfinite(im) & np.isfinite(zo)
    # wavelength per pixel (tabular WCS 'W'); evaluate at block centres, and at every pixel row/col coarsely
    ww = WCS(h["IMAGE"].header, fobj=h, key="W")
    n = im.shape[0]
    nb = n // BLOCK
    yc, xc = np.mgrid[0:nb, 0:nb]
    xc = xc * BLOCK + (BLOCK - 1) / 2.0
    yc = yc * BLOCK + (BLOCK - 1) / 2.0
    lam_b, _ = ww.pixel_to_world_values(xc.ravel(), yc.ravel())
    lam_b = lam_b.reshape(nb, nb)
    ch_b = channel_of(det, lam_b)
    # per-pixel channel by nearest block (blocks are small vs the ~120-pixel channel strips)
    ch_p = np.repeat(np.repeat(ch_b, BLOCK, 0), BLOCK, 1)
    ch_p = np.pad(ch_p, ((0, n - ch_p.shape[0]), (0, n - ch_p.shape[1])), mode="edge")
    res = im - zo
    out_res = res.copy()
    zres = zo.copy()
    chans = [c for c in np.unique(ch_b[ch_b >= 0])]
    hot = np.zeros(good.shape, bool)
    for c in chans:
        inc = ch_p == c
        sel = inc & good
        if sel.sum() < 2000:
            good &= ~inc
            continue
        v = res[sel]
        med = np.median(v)
        sig = 1.4826 * np.median(np.abs(v - med))
        hot |= sel & (res > med + 5 * sig)
    good &= ~binary_dilation(hot, structure=np.ones((5, 5), bool))
    for c in chans:
        inc = ch_p == c
        sel = inc & good
        if sel.sum() < 2000:
            good &= ~inc
            continue
        out_res[inc] -= np.median(res[sel])
        zres[inc] -= np.median(zo[sel])
    g = good.astype(np.float32)
    m = nb * BLOCK

    def bsum(a):
        return a[:m, :m].reshape(nb, BLOCK, nb, BLOCK).sum(axis=(1, 3))

    cnt = bsum(g)
    s_res = bsum(np.where(good, out_res, 0.0))
    s_zo = bsum(np.where(good, zres, 0.0))
    wc = WCS(h["IMAGE"].header)
    ra, dec = wc.all_pix2world(xc.ravel(), yc.ravel(), 0)
    pix = hp.ang2pix(NSIDE, ra, dec, lonlat=True, nest=True)
    keep = (cnt.ravel() > 8) & (ch_b.ravel() >= 0)
    # combine blocks that fall in the same (HEALPix pixel, channel)
    key = pix[keep].astype(np.int64) * 128 + ch_b.ravel()[keep]
    uk, inv = np.unique(key, return_inverse=True)
    agg = lambda a: np.bincount(inv, weights=a.ravel()[keep], minlength=len(uk)).astype(np.float32)
    return dict(hp=(uk // 128).astype(np.int32), ch=(uk % 128).astype(np.int8),
                s=agg(s_res), z=agg(s_zo), n=agg(cnt), det=det, frac_good=float(good.mean()))


def fetch(url, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read()
        except Exception as e:  # network: retry with backoff
            if k == tries - 1:
                raise
            time.sleep(2 ** (k + 1))


def main():
    lst, outdir, wi, nw = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
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
            d = reduce_file(fetch(url))
        except Exception as e:
            print("FAIL", obs, det, repr(e)[:200], flush=True)
            continue
        np.savez_compressed(out + ".tmp.npz", t=float(tmid), **d)
        os.replace(out + ".tmp.npz", out)
        done += 1
        if done % 20 == 0:
            print(f"w{wi} {done} done, {(time.time()-t0)/done:.2f} s/img", flush=True)
    print(f"w{wi} finished {done}", flush=True)


if __name__ == "__main__":
    main()
