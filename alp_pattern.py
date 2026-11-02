"""R-153: accumulate SPHEREx's fixed detector-frame residual pattern from QR2 images.

Each image's IMAGE - ZODI (masked as alp_reduce.py, nu I_nu nW m^-2 sr^-1) has its own sky spectrum removed by
subtracting the median of each ln(lambda) bin (0.002 wide, ~10 rows of the LVF), then pixels beyond 5 robust sigma
are dropped. What is left in detector coordinates, averaged over many exposures of different sky, is the fixed pattern
(flat-field and offset residuals): half-stacks of 40 Draco D1 images correlate at 0.98 (results-R153.md, Step 0).
Per detector and per image list (group) it stores, per pixel: N, sum o, sum o*Z, sum Z, sum Z^2, with Z the image's
median ZODI in nu I_nu (so the pattern can be fitted as a + b Z: an offset plus a flat-field error times the zodi).
A dwarf's own images are left out of the pattern applied to it (alp_reduce.py --pattern), so its line cannot leak in.

Usage: python -u alp_pattern.py <list.csv> <group> <worker_index> <n_workers>  -> data/alp/pattern/<group>.<w>.D<det>.npz
"""
import io
import os
import sys
import time
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from alp_reduce import BAD_BITS, SOURCE_BIT, L0, DL, fetch  # noqa: E402

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data", "alp", "pattern")
C_UM = 2.99792458e14


def lam_map(h, step=8):
    n = h["IMAGE"].data.shape[0]
    g = np.arange(0, n + step, step).clip(0, n - 1)
    yc, xc = np.meshgrid(g, g, indexing="ij")
    ww = WCS(h["IMAGE"].header, fobj=h, key="W")
    lam, _ = ww.pixel_to_world_values(xc.ravel(), yc.ravel())
    lam = lam.reshape(xc.shape)
    p = np.arange(n) / step
    i0 = np.minimum(p.astype(int), len(g) - 2)
    f = p - i0
    b = lam[:, i0] * (1 - f) + lam[:, i0 + 1] * f
    return b[i0, :] * (1 - f[:, None]) + b[i0 + 1, :] * f[:, None]


def pattern_image(data):
    """(det, o, Z): the sky-spectrum-removed, clipped residual in nW, and the image's median zodi in nW."""
    h = fits.open(io.BytesIO(data), memmap=False)
    det = int(h["IMAGE"].header["DETECTOR"])
    im = h["IMAGE"].data.astype(np.float64)
    zo = h["ZODI"].data.astype(np.float64)
    fl = h["FLAGS"].data
    bad = np.zeros(fl.shape, bool)
    for b in BAD_BITS + [SOURCE_BIT]:
        bad |= ((fl >> b) & 1).astype(bool)
    lam = lam_map(h)
    conv = 1e-20 * C_UM / lam * 1e9
    r = np.where(bad | ~np.isfinite(im) | ~np.isfinite(zo), np.nan, (im - zo) * conv)
    Z = float(np.nanmedian(zo * conv))
    lb = np.floor((np.log(lam) - L0) / DL).astype(int)
    ok = np.isfinite(r)
    o = np.full(r.shape, np.nan)
    order = np.argsort(lb[ok], kind="stable")
    vals = r[ok][order]
    keys = lb[ok][order]
    uk, st = np.unique(keys, return_index=True)
    en = np.r_[st[1:], len(keys)]
    med = np.empty(len(vals))
    for a, b in zip(st, en):
        med[a:b] = np.median(vals[a:b]) if b - a >= 500 else np.nan
    tmp = np.empty(len(vals))
    tmp[order] = vals - med
    o[ok] = tmp
    s = 1.4826 * np.nanmedian(np.abs(o))
    o[np.abs(o) > 5 * s] = np.nan
    return det, o, Z


def main():
    lst, group, wi, nw = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    os.makedirs(OUT, exist_ok=True)
    rows = [l.strip().split(",") for l in open(lst) if l.strip() and not l.startswith("obs_id")]
    acc = {}
    t0 = time.time()
    done = 0
    for i, r in enumerate(rows):
        if i % nw != wi:
            continue
        try:
            det, o, Z = pattern_image(fetch(r[3]))
        except Exception as e:
            print("FAIL", r[0], r[1], repr(e)[:200], flush=True)
            continue
        if det not in acc:
            acc[det] = {k: np.zeros(o.shape, np.float64) for k in ("n", "s", "sz", "z", "zz")}
        a = acc[det]
        ok = np.isfinite(o)
        a["n"][ok] += 1
        a["s"][ok] += o[ok]
        a["sz"][ok] += o[ok] * Z
        a["z"][ok] += Z
        a["zz"][ok] += Z * Z
        done += 1
        if done % 10 == 0:
            print(f"w{wi} {done} done, {(time.time()-t0)/done:.2f} s/img", flush=True)
        if done % 50 == 0:
            for det_, a_ in acc.items():
                np.savez(os.path.join(OUT, f"{group}.{wi}.D{det_}.npz"), **{k: v.astype(np.float32) for k, v in a_.items()})
    for det_, a_ in acc.items():
        np.savez(os.path.join(OUT, f"{group}.{wi}.D{det_}.npz"), **{k: v.astype(np.float32) for k, v in a_.items()})
    print(f"w{wi} finished {done}", flush=True)


if __name__ == "__main__":
    main()
