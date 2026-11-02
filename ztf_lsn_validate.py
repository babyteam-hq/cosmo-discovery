"""Validation of the forced photometry and the statistic on real ZTF supernovae (R-130, cosmo-deep-10).

Eight ALeRCE lc_classifier SNe Ia (probability > 0.8, short-lived, peak r 18.4-20.5; `val_targets.csv`), their
cutouts fetched and measured exactly as the lenses' (`--targets val_targets.csv`), compared epoch by epoch with the
alert-stream PSF magnitudes (ALeRCE /objects/<oid>/detections, kept under data/ztf_lsn/val/<oid>.json), and the
statistic's T and its window against the alert peak.

Usage (from research/cosmo-discovery):  python ztf_lsn_validate.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from ztf_lsn_detect import clean, prepare, stat, windows

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "ztf_lsn"


def main() -> None:
    t = pd.read_csv(HERE / "val_targets.csv")
    rows = []
    for r in t.itertuples():
        p = DATA / r.id / "phot.csv"
        if not p.exists():
            continue
        d = clean(pd.read_csv(p))
        f, e, jd, scale = prepare(d, "")
        # kappa0 per epoch (the Gaussian estimator's recovered fraction of a central source), from the PSF stamps
        kap = np.ones_like(f)
        pc = DATA / r.id / "psfcal.csv"
        if pc.exists():
            c = pd.read_csv(pc)
            if len(c) >= 8:
                kap = np.clip(np.polyval(np.polyfit(c.seeing, c.kappa0, 1), d.seeing.values[np.isfinite(d.f.values)]), 0.3, 1.2)
        det = [x for x in json.load(open(DATA / "val" / f"{r.id}.json")) if x["fid"] == 2 and x.get("magpsf")]
        mjd = np.array([x["mjd"] for x in det])
        mag = np.array([x["magpsf"] for x in det])
        emag = np.array([x.get("sigmapsf") or 0.1 for x in det])
        peak_i = int(np.argmin(mag))
        # epoch match
        diffs = []
        for m, mg, em in zip(mjd, mag, emag):
            k = np.argmin(np.abs(jd - (m + 2400000.5)))
            if abs(jd[k] - (m + 2400000.5)) < 2.0 and f[k] > 0:  # the alert exposures are partnership ones, not public: nearest public epoch within 2 d
                diffs.append((23.9 - 2.5 * np.log10(f[k] / kap[k]) - mg, 1.0857 * e[k] / f[k], em))
        w = windows(f, e, jd)
        s = stat(w)
        diffs = np.array(diffs)
        peak_win = (w.t0.values <= mjd[peak_i] + 2400000.5) & (w.t0.values + 15 > mjd[peak_i] + 2400000.5)
        rows.append({
            "id": r.id, "alert peak r": round(float(mag[peak_i]), 2), "alert r dets": len(det),
            "matched epochs": len(diffs), "mine (kappa-corrected, within 2 d) - alert (median)": round(float(np.median(diffs[:, 0])), 3) if len(diffs) else np.nan,
            "rms": round(float(np.std(diffs[:, 0])), 3) if len(diffs) else np.nan,
            "T": round(s["T"], 1), "T window - peak (d)": round(s["T_t0"] - (mjd[peak_i] + 2400000.5), 1),
            "T2": round(s["T2"], 1) if np.isfinite(s["T2"]) else np.nan,
            "max S at peak": round(float(w.S.values[peak_win].max()), 1) if peak_win.any() else np.nan,
            "clean epochs": int(jd.size), "scale": round(scale, 2),
            "window 5-sigma depth": round(float(23.9 - 2.5 * np.log10(5 * np.median(w.ferr))), 2),
        })
    out = pd.DataFrame(rows)
    pd.set_option("display.width", 250)
    print(out.to_string(index=False))
    out.to_csv(DATA / "validation.csv", index=False)


if __name__ == "__main__":
    main()
