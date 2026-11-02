"""R-158 (added after the per-galaxy JADES result): every resolved galaxy in the JADES DR2 GOODS-S F200W mosaic, selected
by our own mirror-invariant cuts, through the same prep and classifiers. Detection: 2048-px blocks (128-px overlap),
Gaussian(1.5 px) smoothing, 3 sigma threshold, connected segments with area > 300 px (0.27 arcsec2) whose flux centroid
lies in the block core; flux-weighted r50 > 7 px (0.21" = 3 x F200W FWHM 0.07") and moment axis ratio > 0.4.
Stamps half-size clip(4 r50, 32, 150), prepped as everywhere else. Output: prep_jadesfield.npz + jadesfield_cat.csv."""
import os, numpy as np, pandas as pd
from scipy import ndimage as ndi
from astropy.io import fits
from astropy.wcs import WCS
from r158_common import DATA
import r158_classify as C
od = os.path.join(DATA, "jades")
hdr = fits.Header.fromtextfile(os.path.join(od, "sci_header.txt")); nx, ny = hdr["NAXIS1"], hdr["NAXIS2"]
a = np.memmap(os.path.join(od, "sci.raw"), dtype=">f4", mode="r", shape=(ny, nx)); w = WCS(hdr)
BS, OV = 2048, 128; cat = []
for y0 in range(0, ny, BS):
    for x0 in range(0, nx, BS):
        ya, yb, xa, xb = max(0, y0 - OV), min(ny, y0 + BS + OV), max(0, x0 - OV), min(nx, x0 + BS + OV)
        blk = np.nan_to_num(np.array(a[ya:yb, xa:xb], dtype=np.float32))
        good = blk != 0
        if good.mean() < 0.2: continue
        v = blk[good]; med = np.median(v); mad = 1.4826 * np.median(np.abs(v - med))
        sm = ndi.gaussian_filter(blk - med, 1.5); smad = 1.4826 * np.median(np.abs(sm[good]))
        lab, n = ndi.label((sm > 3 * smad) & good)
        if n == 0: continue
        idx = np.arange(1, n + 1)
        area = ndi.sum(np.ones_like(blk), lab, idx)
        keep = idx[area > 300]
        objs = ndi.find_objects(lab)
        for l in keep:
            sy, sx = objs[l - 1]
            m = lab[sy, sx] == l; im = np.clip((blk - med)[sy, sx], 0, None) * m
            F = im.sum()
            if F <= 0: continue
            yy, xx = np.mgrid[sy, sx]
            cy, cx = (im * yy).sum() / F, (im * xx).sum() / F
            if not (y0 <= cy + ya < y0 + BS and x0 <= cx + xa < x0 + BS): continue
            r = np.hypot(yy - cy, xx - cx)[m]; f = im[m]; o = np.argsort(r); cf = np.cumsum(f[o])
            r50 = r[o][np.searchsorted(cf, 0.5 * F)]
            qxx = (im * (xx - cx) ** 2).sum() / F; qyy = (im * (yy - cy) ** 2).sum() / F; qxy = (im * (xx - cx) * (yy - cy)).sum() / F
            l1, l2 = np.linalg.eigvalsh([[qxx, qxy], [qxy, qyy]])
            ba = np.sqrt(max(l1, 0) / max(l2, 1e-9))
            cat.append((cx + xa, cy + ya, r50, ba, F, m.sum()))
cat = pd.DataFrame(cat, columns=["x", "y", "r50", "ba", "flux", "area"])
cat["ra"], cat["dec"] = w.all_pix2world(cat.x, cat.y, 0)
sel = cat[(cat.r50 > 7) & (cat.ba > 0.4)].reset_index(drop=True)
print("detected", len(cat), "selected", len(sel), flush=True)
det = np.linalg.det(w.pixel_scale_matrix)
X, keep, r50s, snrs, pxs = [], [], [], [], []
for i, r in sel.iterrows():
    hs = int(np.clip(np.ceil(4 * r.r50), 32, 150)); xi, yi = int(round(r.x)), int(round(r.y))
    if xi - hs < 0 or yi - hs < 0 or xi + hs >= nx or yi + hs >= ny: continue
    s = np.array(a[yi - hs:yi + hs + 1, xi - hs:xi + hs + 1], dtype=np.float32)
    if (s != 0).mean() < 0.8: continue
    if det > 0: s = s[:, ::-1]
    x, info = C.prep(s)
    if x is None: continue
    X.append(x); keep.append(i); r50s.append(info["r50_pix"]); snrs.append(info["snr"]); pxs.append(s.shape[0])
sel = sel.loc[keep].reset_index(drop=True); sel["id"] = np.arange(len(sel))
sel.to_csv(os.path.join(DATA, "jadesfield_cat.csv"), index=False)
np.savez(os.path.join(DATA, "prep_jadesfield.npz"), ids=sel.id.values, X=np.array(X), r50=np.array(r50s), snr=np.array(snrs), px=np.array(pxs))
print("prepped", len(sel))
