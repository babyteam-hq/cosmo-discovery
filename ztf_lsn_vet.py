"""Vetting material for a candidate (R-130, cosmo-deep-10): the lens and control light curves with the best window
marked, the mean difference image inside the window next to the mean outside it, the window's centroid offset from the
lens, and the per-season flux (persistent variability check). g band, if fetched (`ztf_lsn_fetch.py --only ID --filters
zg` then `ztf_lsn_phot.py --only ID --filters zg`... the photometry reads every cutout in cut/, so g epochs appear in
phot.csv with filter ZTF_g once fetched and re-run with --redo).

Usage (from research/cosmo-discovery):
    python ztf_lsn_vet.py ID [--t0 JD]        # writes data/ztf_lsn/<ID>/vet.png and prints the numbers
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ztf_lsn_detect import WIN, clean, prepare, stat, windows  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "ztf_lsn"


def stack(tid: str, d: pd.DataFrame, ra: float, dec: float, sel: np.ndarray) -> tuple[np.ndarray, tuple[float, float], tuple[float, float]]:
    """Mean of the selected cutouts (sky-subtracted), the lens pixel position and the flux-weighted centroid within 4"."""
    imgs, pos = [], None
    for stem in d.stem.values[sel]:
        p = DATA / tid / "cut" / f"{stem}_diff.fits"
        if not p.exists():
            continue
        with fits.open(p) as h:
            img = np.array(h[-1].data, float)
            if img.shape != (31, 31):
                continue
            if pos is None:
                (x0, y0), = WCS(h[-1].header).all_world2pix([[ra, dec]], 0)
                pos = (float(x0), float(y0))
            imgs.append(img - np.nanmedian(img))
    if not imgs or pos is None:
        return np.zeros((31, 31)), (15, 15), (np.nan, np.nan)
    m = np.nanmean(imgs, axis=0)
    yy, xx = np.mgrid[:31, :31]
    r2 = (xx - pos[0]) ** 2 + (yy - pos[1]) ** 2
    w = np.where((r2 < 16) & np.isfinite(m) & (m > 0), m, 0)
    if w.sum() > 0:
        cx, cy = (w * xx).sum() / w.sum(), (w * yy).sum() / w.sum()
    else:
        cx, cy = np.nan, np.nan
    return m, pos, (cx, cy)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("id")
    ap.add_argument("--t0", type=float, default=None)
    a = ap.parse_args()
    t = pd.read_csv(HERE / "lsn_targets.csv").set_index("id").loc[a.id]
    d0 = pd.read_csv(DATA / a.id / "phot.csv")
    dr = clean(d0[d0["filter"].astype(str).str.contains("_r")])
    dg = clean(d0[d0["filter"].astype(str).str.contains("_g")])
    r, e, jd, scale = prepare(dr, "")
    w = windows(r, e, jd)
    s = stat(w)
    t0 = a.t0 if a.t0 is not None else s["T_t0"]
    print(f"{a.id} ({t['class']}): T={s['T']:.2f} at JD {s['T_t0']:.0f} (n={s['T_n']}), T2={s['T2']:.2f}, scale {scale:.2f}, "
          f"{jd.size} r epochs")
    inwin = (jd >= t0) & (jd < t0 + WIN)
    print(f"window {t0:.0f}-{t0 + WIN:.0f}: epochs {inwin.sum()}, mean flux {np.average(r[inwin], weights=1 / e[inwin] ** 2) if inwin.any() else np.nan:.2f} uJy")
    for name, (lo, hi) in {"-30..-15": (-30, -15), "-15..0": (-15, 0), "0..15": (0, 15), "15..30": (15, 30),
                           "30..45": (30, 45), "45..60": (45, 60), "60..90": (60, 90)}.items():
        m = (jd >= t0 + lo) & (jd < t0 + hi)
        if m.any():
            print(f"  r {name:>8} d: n={m.sum():2d} mean {np.average(r[m], weights=1 / e[m] ** 2):7.2f} +- {1 / np.sqrt((1 / e[m] ** 2).sum()):5.2f} uJy")
    if len(dg):
        rg, eg, jdg, scaleg = prepare(dg, "")
        for name, (lo, hi) in {"-30..-15": (-30, -15), "-15..0": (-15, 0), "0..15": (0, 15), "15..30": (15, 30),
                               "30..45": (30, 45), "45..60": (45, 60)}.items():
            m = (jdg >= t0 + lo) & (jdg < t0 + hi)
            if m.any():
                print(f"  g {name:>8} d: n={m.sum():2d} mean {np.average(rg[m], weights=1 / eg[m] ** 2):7.2f} +- {1 / np.sqrt((1 / eg[m] ** 2).sum()):5.2f} uJy")
    # controls in the same window
    for pos in ("_N", "_S", "_E", "_W"):
        rc, ec, jdc, _ = prepare(dr, pos)
        m = (jdc >= t0) & (jdc < t0 + WIN)
        if m.any():
            print(f"  control {pos}: window mean {np.average(rc[m], weights=1 / ec[m] ** 2):6.2f} +- {1 / np.sqrt((1 / ec[m] ** 2).sum()):5.2f}")
    # seasons
    yrs = np.floor((jd - 2458000) / 365.25)
    for y in np.unique(yrs):
        m = yrs == y
        print(f"  season {int(y)}: n={m.sum():3d}, mean {np.average(r[m], weights=1 / e[m] ** 2):6.2f} +- {1 / np.sqrt((1 / e[m] ** 2).sum()):4.2f}, "
              f"max window S in season {w.S[(w.t0 >= jd[m].min()) & (w.t0 <= jd[m].max())].max() if ((w.t0 >= jd[m].min()) & (w.t0 <= jd[m].max())).any() else np.nan:.2f}")
    # stamps
    sel_in = (dr.jd.values >= t0) & (dr.jd.values < t0 + WIN)
    sel_out = ~sel_in & (np.abs(dr.jd.values - t0) < 400)
    m_in, pos, cen = stack(a.id, dr, t.ra, t.dec, sel_in)
    m_out, _, _ = stack(a.id, dr, t.ra, t.dec, sel_out)
    off = np.hypot(cen[0] - pos[0], cen[1] - pos[1]) * 1.0126 if np.isfinite(cen[0]) else np.nan
    print(f"  window stack centroid offset from the lens: {off:.2f} arcsec (lens at pixel {pos[0]:.1f},{pos[1]:.1f})")
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    ax[0, 0].errorbar(jd - 2458000, r, e, fmt=".", ms=3, alpha=0.5, label="lens, r")
    ax[0, 0].plot(w.t0 + WIN / 2 - 2458000, w.fmean, "r-", lw=1, label="15-d window mean")
    ax[0, 0].axvspan(t0 - 2458000, t0 + WIN - 2458000, color="orange", alpha=0.3)
    ax[0, 0].set_ylim(np.percentile(r, 1) - 10, max(np.percentile(r, 99.5) + 10, (w.fmean.max() + 5)))
    ax[0, 0].set_xlabel("JD - 2458000"); ax[0, 0].set_ylabel("difference flux (uJy)"); ax[0, 0].legend(); ax[0, 0].set_title(f"{a.id}  T={s['T']:.2f}")
    m = np.abs(jd - t0) < 120
    ax[0, 1].errorbar(jd[m] - t0, r[m], e[m], fmt="o", ms=4, label="lens r")
    if len(dg):
        mg = np.abs(jdg - t0) < 120
        ax[0, 1].errorbar(jdg[mg] - t0, rg[mg], eg[mg], fmt="s", ms=3, color="g", alpha=0.6, label="lens g")
    for pos_ in ("_N", "_S", "_E", "_W"):
        rc, ec, jdc, _ = prepare(dr, pos_)
        mc = np.abs(jdc - t0) < 120
        ax[0, 1].plot(jdc[mc] - t0, rc[mc], ".", color="gray", alpha=0.3)
    ax[0, 1].axvspan(0, WIN, color="orange", alpha=0.3); ax[0, 1].axhline(0, color="k", lw=0.5)
    ax[0, 1].set_xlabel(f"days from JD {t0:.0f}"); ax[0, 1].legend(); ax[0, 1].set_title("around the window (gray: controls)")
    vmax = max(3 * np.nanstd(m_out), 1e-3)
    for axx, img, ttl in ((ax[1, 0], m_in, f"mean difference in window ({sel_in.sum()} epochs)"),
                          (ax[1, 1], m_out, f"mean difference outside ({sel_out.sum()} epochs within 400 d)")):
        axx.imshow(img, origin="lower", cmap="gray", vmin=-vmax, vmax=vmax)
        axx.plot(pos[0], pos[1], "r+", ms=14, mew=1.5)
        axx.add_patch(plt.Circle(pos, 1.2 / 1.0126, fill=False, color="cyan", lw=0.8))
        axx.set_title(ttl)
    if np.isfinite(cen[0]):
        ax[1, 0].plot(cen[0], cen[1], "yx", ms=10)
    fig.tight_layout()
    fig.savefig(DATA / a.id / "vet.png", dpi=110)
    print(f"  -> {DATA / a.id / 'vet.png'}")


if __name__ == "__main__":
    sys.exit(main())
