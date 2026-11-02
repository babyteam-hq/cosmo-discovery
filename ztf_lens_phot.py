"""Scene-model photometry on ZTF difference-image cutouts (R-127, cosmo-deep-7).

Each difference-image cutout (from ztf_lens_fetch.py) is fitted with a fixed-position multi-PSF scene: one PSF (the
pipeline's difference-image PSF stamp for that epoch) per lensed image at the discovery-paper position, a constant
background, and one global astrometric shift (dx, dy) shared by all images; only the fluxes, the background and the
shift are free. Pixel weights are 1/(sigma_sky^2 + |model|/gain), iterated once; pixels off by more than 5 sigma are
masked and the fit repeated. The difference flux is converted to microJansky with the epoch's MAGZP (the difference
image is in the science image's units). The reference flux of each image comes from the reference PSF-fit catalogue
(nearest source within 1.5"), so total flux = reference + difference.

Usage (from research/cosmo-discovery):
    python ztf_lens_phot.py <target>            # writes data/ztf_lens/<target>/phot.csv
    python ztf_lens_phot.py <target> --limit 50 # first 50 epochs, for a look
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage, optimize

from ztf_lens_targets import TARGETS

DATA = Path(__file__).resolve().parent / "data" / "ztf_lens"
warnings.filterwarnings("ignore")


def place_psf(psf: np.ndarray, shape: tuple[int, int], x: float, y: float) -> np.ndarray:
    """The PSF stamp (centre at its middle pixel) placed at fractional pixel (x, y) of an image of `shape`."""
    k = psf.shape[0]
    c = k // 2
    pad = 3
    p = np.pad(psf, pad)
    ix, iy = int(np.floor(x)), int(np.floor(y))
    fx, fy = x - ix, y - iy
    s = ndimage.shift(p, (fy, fx), order=3, mode="constant", cval=0.0)
    out = np.zeros(shape)
    y0, x0 = iy - c - pad, ix - c - pad
    ys, xs = slice(max(y0, 0), min(y0 + s.shape[0], shape[0])), slice(max(x0, 0), min(x0 + s.shape[1], shape[1]))
    sy, sx = slice(ys.start - y0, ys.stop - y0), slice(xs.start - x0, xs.stop - x0)
    out[ys, xs] = s[sy, sx]
    return out


def solve(img, w, mask, psf, pos, shift):
    """Linear LSQ for fluxes + background at a given global shift. Returns params, cov, model, chi2, n."""
    cols = [place_psf(psf, img.shape, x + shift[0], y + shift[1]) for x, y in pos]
    cols.append(np.ones_like(img))
    A = np.stack([c[mask] for c in cols], axis=1)
    sw = np.sqrt(w[mask])
    Aw = A * sw[:, None]
    bw = img[mask] * sw
    p, *_ = np.linalg.lstsq(Aw, bw, rcond=None)
    model = sum(pi * ci for pi, ci in zip(p, cols))
    r = (img - model)[mask] * sw
    chi2 = float(r @ r)
    try:
        cov = np.linalg.inv(Aw.T @ Aw)
    except np.linalg.LinAlgError:
        cov = np.full((len(p), len(p)), np.nan)
    return p, cov, model, chi2, int(mask.sum())


def fit_epoch(diff_path: Path, psf_path: Path, images: dict, fit_shift: bool = True) -> dict | None:
    with fits.open(diff_path) as h:
        hdu = h[-1]
        img = np.array(hdu.data, dtype=float)
        hd = hdu.header
    with fits.open(psf_path) as h:
        psf = np.array(h[0].data, dtype=float)
    psf /= psf.sum()
    wcs = WCS(hd)
    names = list(images)
    sky = np.array([images[n] for n in names])
    px = wcs.all_world2pix(sky, 0)  # 0-based (x, y)
    ny, nx = img.shape
    inside = (px[:, 0] > 4) & (px[:, 0] < nx - 5) & (px[:, 1] > 4) & (px[:, 1] < ny - 5)
    if inside.sum() == 0:
        return None
    finite = np.isfinite(img)
    if finite.sum() < 0.5 * img.size:
        return None
    # sky noise: clipped std away from the images
    yy, xx = np.mgrid[:ny, :nx]
    near = np.zeros_like(img, bool)
    for x, y in px[inside]:
        near |= (xx - x) ** 2 + (yy - y) ** 2 < 36
    bg = img[finite & ~near]
    if bg.size < 100:
        return None
    med = np.median(bg)
    mad = 1.4826 * np.median(np.abs(bg - med))
    bg = bg[np.abs(bg - med) < 4 * mad]
    sig = bg.std() if bg.size > 20 else mad
    if not (0.5 < sig < 200) or abs(med) > 50:  # flat fill (no reference coverage) or garbage at a quadrant edge
        return None
    gain = float(hd.get("GAIN", 6.2))
    w = np.full_like(img, 1.0 / sig**2)
    mask = finite.copy()
    pos = px[inside]
    pos_names = [n for n, i in zip(names, inside) if i]

    def chi2_at(shift):
        return solve(img, w, mask, psf, pos, shift)[3]

    shift = np.zeros(2)
    p, cov, model, chi2, n = solve(img, w, mask, psf, pos, shift)
    # Poisson term, re-weight, clip, refit
    w = 1.0 / (sig**2 + np.abs(model) / gain)
    p, cov, model, chi2, n = solve(img, w, mask, psf, pos, shift)
    resid = (img - model) * np.sqrt(w)
    mask &= np.abs(resid) < 5
    p, cov, model, chi2, n = solve(img, w, mask, psf, pos, shift)
    snr = np.abs(p[:-1]) / np.sqrt(np.abs(np.diag(cov)[:-1]))
    if fit_shift and np.nanmax(snr) > 8:
        r = optimize.minimize(chi2_at, shift, method="Nelder-Mead",
                              options={"xatol": 0.02, "fatol": 0.1, "maxiter": 60})
        if np.all(np.abs(r.x) < 1.5):
            shift = r.x
            p, cov, model, chi2, n = solve(img, w, mask, psf, pos, shift)
    chi2red = chi2 / max(n - len(p), 1)
    err = np.sqrt(np.abs(np.diag(cov))) * max(1.0, np.sqrt(chi2red))
    zp = float(hd["MAGZP"])
    tojy = 10 ** (-0.4 * (zp - 23.9))  # DN -> microJansky (AB)
    out = {
        "jd": float(hd.get("OBSJD", hd.get("OBSMJD", 0) + 2400000.5)),
        "filter": hd.get("FILTER", hd.get("FILTERID", "")), "seeing": float(hd.get("SEEING", np.nan)),
        "maglim": float(hd.get("MAGLIM", np.nan)), "magzp": zp, "magzprms": float(hd.get("MAGZPRMS", np.nan)),
        "airmass": float(hd.get("AIRMASS", np.nan)), "infobits": int(hd.get("INFOBITS", 0)),
        "dx": shift[0], "dy": shift[1], "chi2red": chi2red, "npix": n, "sky_sig": sig, "bg": p[-1],
        "nx": nx, "ny": ny,
    }
    for nme in names:
        out[f"f_{nme}"] = np.nan
        out[f"e_{nme}"] = np.nan
    for nme, pi, ei in zip(pos_names, p[:-1], err[:-1]):
        out[f"f_{nme}"] = pi * tojy
        out[f"e_{nme}"] = ei * tojy
    return out


def ref_mags(target: str, images: dict) -> dict:
    """Reference PSF-fit magnitude per image per (field, filter, ccd, qid), nearest catalogue source within 1.5"."""
    out = {}
    for cat in sorted((DATA / target / "ref").glob("*_refpsfcat.fits")):
        key = cat.name.replace("_refpsfcat.fits", "")
        try:
            t = fits.getdata(cat, 1)
        except Exception:
            continue
        ra, dec = t["ra"], t["dec"]
        res = {}
        for n, (r0, d0) in images.items():
            d = np.hypot((ra - r0) * np.cos(np.radians(d0)), dec - d0) * 3600
            i = int(np.argmin(d))
            res[n] = (float(t["mag"][i]), float(t["sigmag"][i]), float(d[i])) if d[i] < 1.5 else (np.nan, np.nan, float(d[i]))
        out[key] = res
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("target")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--no-shift", action="store_true")
    a = ap.parse_args()
    t = TARGETS[a.target]
    images = t["images"]
    cut = DATA / a.target / "cut"
    diffs = sorted(cut.glob("*_diff.fits"))
    if a.limit:
        diffs = diffs[: a.limit]
    rows = []
    for i, d in enumerate(diffs):
        psf = d.with_name(d.name.replace("_diff.fits", "_psf.fits"))
        if not psf.exists():
            continue
        try:
            r = fit_epoch(d, psf, images, fit_shift=not a.no_shift)
        except Exception as e:  # a broken cutout
            print(f"{d.name}: {e}", file=sys.stderr)
            continue
        if r is None:
            continue
        stem = d.name.replace("_diff.fits", "")
        ffd, field, filt, ccd, qid = stem.split("_")
        r.update({"stem": stem, "filtercode": filt, "field": int(field), "ccdid": int(ccd[1:]), "qid": int(qid[1:])})
        rows.append(r)
        if (i + 1) % 200 == 0:
            print(f"  {i + 1}/{len(diffs)}", flush=True)
    import pandas as pd

    df = pd.DataFrame(rows).sort_values("jd")
    out = DATA / a.target / ("phot.csv" if not a.limit else "phot_sample.csv")
    df.to_csv(out, index=False)
    print(f"{a.target}: {len(df)} epochs fitted -> {out}")
    refs = ref_mags(a.target, images)
    import json

    (DATA / a.target / "refmags.json").write_text(json.dumps(refs, indent=1))
    for k, v in refs.items():
        print(k, {n: f"{m:.2f}+-{e:.2f} ({d:.2f}\")" for n, (m, e, d) in v.items()})


if __name__ == "__main__":
    main()
