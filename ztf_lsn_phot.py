"""Forced photometry on the ZTF difference cutouts at the lens centre and four control positions (R-130, cosmo-deep-10).

Each 30" difference cutout (ztf_lsn_fetch.py) is fitted at the catalogue lens position with a Gaussian PSF of FWHM = the
header SEEING plus a constant background, by weighted least squares in an 11 x 11 pixel box, weights from the clipped sky
variance of the cutout outside 3" of the centre; pixels off by more than 5 sigma are masked and the fit repeated once. The
same at four control positions 10" N, S, E and W. Flux to microJansky (AB) with MAGZP. For the injections the estimator's
response to a unit point source at offsets 0-3" from the lens centre (same Gaussian, four position angles averaged) is
stored per epoch: the estimator is linear in the image, so an injected image at offset r adds flux x resp(r) to the
measured flux.

Usage (from research/cosmo-discovery):
    python ztf_lsn_phot.py [--only ID,ID] [--procs 8] [--redo]      # writes data/ztf_lsn/<id>/phot.csv for every fetched target
"""
from __future__ import annotations

import argparse
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "ztf_lsn"
warnings.filterwarnings("ignore")
BOX = 5  # half-size of the fit box in pixels (11 x 11)
CONTROL_OFFSET = 10.0  # arcsec
RESP_OFFSETS = (0.0, 0.5, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0)  # arcsec


def gauss(xx, yy, x0, y0, sigma):
    g = np.exp(-0.5 * ((xx - x0) ** 2 + (yy - y0) ** 2) / sigma**2)
    return g / (2 * np.pi * sigma**2)  # unit total flux


def fit_at(img, x, y, sigma, sky_sig):
    """Flux, error, background, chi2red, n and the estimator's flux row for a point source at (x, y)."""
    ix, iy = int(round(x)), int(round(y))
    ny, nx = img.shape
    if ix - BOX < 0 or iy - BOX < 0 or ix + BOX >= nx or iy + BOX >= ny:
        return None
    sub = img[iy - BOX:iy + BOX + 1, ix - BOX:ix + BOX + 1]
    yy, xx = np.mgrid[iy - BOX:iy + BOX + 1, ix - BOX:ix + BOX + 1]
    g = gauss(xx, yy, x, y, sigma)
    A = np.stack([g.ravel(), np.ones(g.size)], axis=1)
    d = sub.ravel()
    mask = np.isfinite(d)
    if mask.sum() < 0.8 * d.size:
        return None
    w = np.full(d.size, 1.0 / sky_sig**2)
    for _ in range(2):
        Aw = A[mask] * np.sqrt(w[mask])[:, None]
        p, *_ = np.linalg.lstsq(Aw, d[mask] * np.sqrt(w[mask]), rcond=None)
        resid = (d - A @ p) * np.sqrt(w)
        newmask = mask & (np.abs(resid) < 5)
        if newmask.sum() == mask.sum():
            break
        mask = newmask
    try:
        cov = np.linalg.inv(Aw.T @ Aw)
    except np.linalg.LinAlgError:
        return None
    if mask.sum() < 0.5 * d.size or not np.isfinite(cov[0, 0]) or cov[0, 0] <= 0:
        return None
    chi2 = float(resid[mask] @ resid[mask])
    n = int(mask.sum())
    chi2red = chi2 / max(n - 2, 1)
    # the estimator's linear row for the flux parameter: a_hat = row . d
    row = np.zeros(d.size)
    row[mask] = (cov @ (A[mask] * w[mask][:, None]).T)[0]
    return p[0], float(np.sqrt(cov[0, 0])), p[1], chi2red, n, row, (xx, yy, mask)


def response(row, grid, x, y, sigma):
    """Estimator response to a unit point source at offsets RESP_OFFSETS (arcsec) from (x, y), 4 angles averaged."""
    xx, yy, _ = grid
    out = []
    for off in RESP_OFFSETS:
        vals = []
        for ang in (0, 90, 180, 270):
            dx = off / PIX * np.cos(np.radians(ang))
            dy = off / PIX * np.sin(np.radians(ang))
            vals.append(row @ gauss(xx, yy, x + dx, y + dy, sigma).ravel())
        out.append(float(np.mean(vals)))
    return out


PIX = 1.0126  # ZTF arcsec per pixel


def one_epoch(args):
    path, ra, dec = args
    try:
        with fits.open(path) as h:
            hdu = h[-1]
            img = np.array(hdu.data, dtype=float)
            hd = hdu.header
    except Exception:
        return None
    if img.ndim != 2 or img.size < 400:
        return None
    seeing = float(hd.get("SEEING", np.nan))
    if not np.isfinite(seeing) or seeing <= 0:
        return None
    wcs = WCS(hd)
    (x0, y0), = wcs.all_world2pix([[ra, dec]], 0)
    ny, nx = img.shape
    yy, xx = np.mgrid[:ny, :nx]
    far = (xx - x0) ** 2 + (yy - y0) ** 2 > (3.0 / PIX) ** 2
    bg = img[far & np.isfinite(img)]
    if bg.size < 100:
        return None
    med = np.median(bg)
    mad = 1.4826 * np.median(np.abs(bg - med))
    if mad <= 0:
        return None
    bgc = bg[np.abs(bg - med) < 4 * mad]
    sky = float(bgc.std()) if bgc.size > 20 else mad
    if not (0.3 < sky < 500) or abs(med) > 100:
        return None
    sigma = seeing / 2.3548 / PIX
    zp = float(hd["MAGZP"])
    tojy = 10 ** (-0.4 * (zp - 23.9))
    out = {
        "jd": float(hd.get("OBSJD", hd.get("OBSMJD", 0) + 2400000.5)), "filter": str(hd.get("FILTER", "")),
        "seeing": seeing, "maglim": float(hd.get("MAGLIM", np.nan)), "magzp": zp,
        "magzprms": float(hd.get("MAGZPRMS", np.nan)), "infobits": int(hd.get("INFOBITS", 0)),
        "field": int(hd.get("FIELDID", 0)), "ccdid": int(hd.get("CCDID", 0)), "qid": int(hd.get("QID", 0)),
        "sky_sig": sky * tojy, "nanfrac": float(1 - np.isfinite(img).mean()), "stem": path.name.replace("_diff.fits", ""),
    }
    r = fit_at(img, x0, y0, sigma, sky)
    if r is None:
        return None
    f, e, b, chi2red, n, row, grid = r
    out.update({"f": f * tojy, "e": e * tojy, "bg": b * tojy, "chi2red": chi2red, "npix": n})
    for off, resp in zip(RESP_OFFSETS, response(row, grid, x0, y0, sigma)):
        out[f"resp_{off:.1f}"] = resp
    cosd = np.cos(np.radians(dec))
    for name, (dra, ddec) in {"N": (0, 1), "S": (0, -1), "E": (1, 0), "W": (-1, 0)}.items():
        (xc, yc), = wcs.all_world2pix([[ra + dra * CONTROL_OFFSET / 3600 / cosd, dec + ddec * CONTROL_OFFSET / 3600]], 0)
        rc = fit_at(img, xc, yc, sigma, sky)
        if rc is None:
            out[f"f_{name}"] = np.nan
            out[f"e_{name}"] = np.nan
        else:
            out[f"f_{name}"] = rc[0] * tojy
            out[f"e_{name}"] = rc[1] * tojy
    return out


def process_target(tid: str, ra: float, dec: float, procs: int, redo: bool) -> str:
    d = DATA / tid
    out = d / "phot.csv"
    cuts = sorted((d / "cut").glob("*_diff.fits")) if (d / "cut").exists() else []
    if not cuts:
        return f"{tid}: no cutouts"
    if out.exists() and not redo:
        n_done = sum(1 for _ in open(out)) - 1
        if n_done >= len(cuts) - 5:
            return f"{tid}: cached ({n_done})"
    args = [(c, ra, dec) for c in cuts]
    if procs > 1:
        with ProcessPoolExecutor(procs) as ex:
            rows = [r for r in ex.map(one_epoch, args, chunksize=32) if r is not None]
    else:
        rows = [r for r in map(one_epoch, args) if r is not None]
    if not rows:
        return f"{tid}: nothing fitted"
    df = pd.DataFrame(rows).sort_values("jd")
    df.to_csv(out, index=False, float_format="%.6g")
    return f"{tid}: {len(df)}/{len(cuts)} epochs"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    ap.add_argument("--targets", default="lsn_targets.csv", help="target list (lsn_targets.csv, or val_targets.csv for the validation SNe)")
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--redo", action="store_true")
    ap.add_argument("--wait-done", action="store_true", help="only targets whose fetch marker exists")
    a = ap.parse_args()
    t = pd.read_csv(HERE / a.targets)
    if a.only:
        t = t[t.id.isin(a.only.split(","))]
    for r in t.itertuples():
        if a.wait_done and not list((DATA / r.id).glob("done_*")) if (DATA / r.id).exists() else a.wait_done:
            continue
        print(process_target(r.id, r.ra, r.dec, a.procs, a.redo), flush=True)


if __name__ == "__main__":
    sys.exit(main())
