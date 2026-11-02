"""PSF calibration of the forced photometry (R-130, cosmo-deep-10).

The difference image's PSF is not the science seeing: the pipeline convolves the reference and the science image to a
common PSF, broader than either, and `ztf_lsn_phot.py` fits a Gaussian of FWHM = SEEING. For every k-th epoch of a target
(default every 10th) this fetches the pipeline's difference-image PSF stamp (`diffimgpsf.fits`, 25 x 25), measures its
FWHM, and computes kappa = the flux the Gaussian estimator returns for a unit source with that true PSF (the estimator's
linear row applied to the stamp), and the same at the 1.2" image offset. Writes data/ztf_lsn/<id>/psfcal.csv; the
analysis interpolates kappa to every epoch by the stamp's relation to SEEING within each (field, ccd, quadrant).

Usage (from research/cosmo-discovery):
    python ztf_lsn_psfcal.py [--targets lsn_targets.csv] [--every 10] [--workers 6] [--only ID] [--all-epochs]
"""
from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy import ndimage

from ztf_lsn_fetch import IBE, get, session
from ztf_lsn_phot import BOX, PIX, gauss

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "ztf_lsn"


def psf_url(stem: str) -> str:
    ffd, field, filt, ccd, qid = stem.split("_")
    return (f"{IBE}/data/ztf/products/sci/{ffd[:4]}/{ffd[4:8]}/{ffd[8:]}/"
            f"ztf_{ffd}_{field}_{filt}_{ccd}_o_{qid}_diffimgpsf.fits")


def estimator_row(sigma: float, n: int = 2 * BOX + 1) -> np.ndarray:
    """The flux row of the (flux, background) weighted LSQ estimator for a Gaussian at the box centre."""
    c = n // 2
    yy, xx = np.mgrid[:n, :n]
    g = gauss(xx, yy, c, c, sigma)
    A = np.stack([g.ravel(), np.ones(g.size)], axis=1)
    cov = np.linalg.inv(A.T @ A)
    return (cov @ A.T)[0]


def kappa(psf: np.ndarray, seeing: float, offset_arcsec: float = 0.0) -> float:
    """Flux the Gaussian(seeing) estimator returns for a unit source with PSF `psf` at `offset` (4 angles averaged)."""
    n = 2 * BOX + 1
    c = n // 2
    row = estimator_row(seeing / 2.3548 / PIX, n)
    p = psf / psf.sum()
    k = p.shape[0]
    vals = []
    for ang in ((0,) if offset_arcsec == 0 else (0, 90, 180, 270)):
        dx = offset_arcsec / PIX * np.cos(np.radians(ang))
        dy = offset_arcsec / PIX * np.sin(np.radians(ang))
        big = np.zeros((n + k, n + k))
        big[(n + k) // 2 - k // 2:(n + k) // 2 - k // 2 + k, (n + k) // 2 - k // 2:(n + k) // 2 - k // 2 + k] = p
        s = ndimage.shift(big, (dy, dx), order=1, mode="constant")
        m0 = (n + k) // 2 - c
        sub = s[m0:m0 + n, m0:m0 + n]
        vals.append(float(row @ sub.ravel()))
    return float(np.mean(vals))


def fwhm_of(psf: np.ndarray) -> float:
    """FWHM in arcsec from the second moments of the stamp (Gaussian-equivalent)."""
    p = np.clip(psf, 0, None)
    p /= p.sum()
    yy, xx = np.mgrid[:p.shape[0], :p.shape[1]]
    cx, cy = (p * xx).sum(), (p * yy).sum()
    var = (p * ((xx - cx) ** 2 + (yy - cy) ** 2)).sum() / 2
    return float(2.3548 * np.sqrt(max(var, 1e-6)) * PIX)


def run_target(tid: str, every: int, workers: int, all_epochs: bool, use_cached: bool = False) -> str:
    p = DATA / tid / "phot.csv"
    if not p.exists():
        return f"{tid}: no phot.csv"
    d = pd.read_csv(p)
    d = d[d["filter"].astype(str).str.contains("_r")]
    stems = list(d.stem.values) if all_epochs else list(d.stem.values[::every])
    pdir = DATA / tid / "psf"
    if use_cached:  # R-141: the stamps ztf_lsn_refetch.py queued (every 10th clean r epoch), no new fetches
        have = {f.name[:-len("_psf.fits")] for f in pdir.glob("*_psf.fits") if f.stat().st_size > 0}
        stems = [s for s in d.stem.values if s in have]
    pdir.mkdir(exist_ok=True)
    jobs = [(psf_url(s), pdir / f"{s}_psf.fits") for s in stems]
    jobs = [j for j in jobs if not j[1].with_suffix(".fits.404").exists()]
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(lambda j: get(*j), jobs))
    rows = []
    see = dict(zip(d.stem, d.seeing))
    for s in stems:
        f = pdir / f"{s}_psf.fits"
        if not f.exists() or f.stat().st_size == 0:
            continue
        try:
            psf = np.array(fits.getdata(f), float)
        except Exception:
            continue
        if psf.ndim != 2 or not np.isfinite(psf).all() or psf.sum() <= 0:
            continue
        rows.append({"stem": s, "seeing": see[s], "fwhm_psf": fwhm_of(psf), "kappa0": kappa(psf, see[s], 0.0),
                     "kappa12": kappa(psf, see[s], 1.2), "kappa20": kappa(psf, see[s], 2.0)})
    out = pd.DataFrame(rows)
    out.to_csv(DATA / tid / "psfcal.csv", index=False, float_format="%.5g")
    if len(out):
        return (f"{tid}: {len(out)} stamps; seeing {out.seeing.median():.2f}, PSF FWHM {out.fwhm_psf.median():.2f}, "
                f"kappa0 median {out.kappa0.median():.3f} [{out.kappa0.quantile(0.1):.3f}, {out.kappa0.quantile(0.9):.3f}], "
                f"kappa12 {out.kappa12.median():.3f}")
    return f"{tid}: no stamps"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="lsn_targets.csv")
    ap.add_argument("--every", type=int, default=10)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--only", default="")
    ap.add_argument("--all-epochs", action="store_true")
    ap.add_argument("--redo", action="store_true")
    ap.add_argument("--use-cached", action="store_true", help="measure the stamps already fetched (ztf_lsn_refetch.py)")
    ap.add_argument("--stop", type=int, default=10**6, help="ranks below this only")
    a = ap.parse_args()
    t = pd.read_csv(HERE / a.targets)
    t = t[t["rank"] < a.stop]
    if a.only:
        t = t[t.id.isin(a.only.split(","))]
    for r in t.itertuples():
        if (DATA / r.id / "psfcal.csv").exists() and not a.redo and not a.only:
            continue
        print(run_target(r.id, a.every, a.workers, a.all_epochs, a.use_cached), flush=True)


if __name__ == "__main__":
    sys.exit(main())
