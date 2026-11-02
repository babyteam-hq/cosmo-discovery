"""The detection statistic on the forced light curves (R-130, cosmo-deep-10), as fixed in results-R130.md.

Per target: the epoch cuts (INFOBITS != 0, SEEING > 4", MAGLIM < 19, more than half the cutout undefined), the baseline
(median flux per position and (field, ccd, quadrant)), the error scaling (robust std of (f - baseline)/e, 3-sigma clipped,
floored at 1), then 15-day windows every 5 days: S_w = sum(f/sigma^2)/sqrt(sum(1/sigma^2)) over the window's epochs (at
least two); T = max_w S_w. The same for the four control positions. Writes data/ztf_lsn/detect.csv (one row per target and
position) and, per target, windows.csv (the S_w series for the lens position, for the plots and the vetting).

Usage (from research/cosmo-discovery):
    python ztf_lsn_detect.py                 # every target with a phot.csv
    python ztf_lsn_detect.py --only ID       # one target, prints its top windows
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "ztf_lsn"
WIN = 15.0
STEP = 5.0
POSITIONS = ("", "_N", "_S", "_E", "_W")


KEEP_INFOBITS = False  # robustness variant only (--keep-infobits): the pre-registered cut drops INFOBITS != 0


def clean(d: pd.DataFrame) -> pd.DataFrame:
    is_r = d["filter"].astype(str).str.contains("_r")  # g epochs exist only for candidates whose g band was fetched
    ok = is_r & ((d.infobits == 0) | KEEP_INFOBITS) & (d.seeing <= 4.0) & (d.maglim >= 19.0) & (d.nanfrac < 0.5) & np.isfinite(d.f) & (d.e > 0)
    return d[ok].sort_values("jd").reset_index(drop=True)


def prepare(d: pd.DataFrame, pos: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Baseline-subtracted flux, scaled error and jd for one position; also the error scale factor."""
    f, e = d["f" + pos].values.astype(float), d["e" + pos].values.astype(float)
    ok = np.isfinite(f) & np.isfinite(e) & (e > 0)
    f, e, jd = f[ok], e[ok], d.jd.values[ok]
    key = (d.field.astype(str) + "_" + d.ccdid.astype(str) + "_" + d.qid.astype(str)).values[ok]
    base = np.zeros_like(f)
    for k in np.unique(key):
        m = key == k
        base[m] = np.median(f[m]) if m.sum() >= 5 else np.median(f)
    r = f - base
    z = r / e
    for _ in range(3):
        mad = 1.4826 * np.median(np.abs(z - np.median(z)))
        z = z[np.abs(z - np.median(z)) < 3 * max(mad, 1e-6)] if mad > 0 else z
    scale = max(1.0, float(np.std(z))) if z.size > 10 else 1.0
    return r, e * scale, jd, scale


def windows(r: np.ndarray, e: np.ndarray, jd: np.ndarray) -> pd.DataFrame:
    """S_w per 15-day window every 5 days (windows with >= 2 epochs)."""
    if jd.size == 0:
        return pd.DataFrame(columns=["t0", "n", "S", "fmean", "ferr"])
    starts = np.arange(np.floor(jd.min()), jd.max() + 1e-6, STEP)
    w = 1.0 / e**2
    i0 = np.searchsorted(jd, starts)
    i1 = np.searchsorted(jd, starts + WIN)
    cw = np.concatenate([[0], np.cumsum(w)])
    cwf = np.concatenate([[0], np.cumsum(w * r)])
    n = i1 - i0
    sw = cw[i1] - cw[i0]
    swf = cwf[i1] - cwf[i0]
    ok = (n >= 2) & (sw > 0)
    out = pd.DataFrame({"t0": starts[ok], "n": n[ok], "S": swf[ok] / np.sqrt(sw[ok]), "fmean": swf[ok] / sw[ok],
                        "ferr": 1 / np.sqrt(sw[ok])})
    return out


def stat(w: pd.DataFrame) -> dict:
    if len(w) == 0:
        return {"T": np.nan, "T_t0": np.nan, "T_n": 0, "T2": np.nan, "nwin": 0}
    i = int(np.argmax(w.S.values))
    # T2: the best window whose neighbour (the window 15 d later or earlier) is also positive at S > 3 (the shape cut)
    S = w.S.values
    t0 = w.t0.values
    t2 = -np.inf
    for j in range(len(w)):
        nb = (np.abs(t0 - t0[j] - WIN) < 1e-6) | (np.abs(t0 - t0[j] + WIN) < 1e-6)
        if nb.any() and S[nb].max() > 3.0:
            t2 = max(t2, S[j])
    return {"T": float(S[i]), "T_t0": float(t0[i]), "T_n": int(w.n.values[i]), "T2": float(t2) if np.isfinite(t2) else np.nan,
            "nwin": int(len(w))}


def run_target(tid: str, verbose: bool = False) -> list[dict]:
    p = DATA / tid / "phot.csv"
    if not p.exists():
        return []
    d = pd.read_csv(p)
    d = clean(d)
    rows = []
    for pos in POSITIONS:
        r, e, jd, scale = prepare(d, pos)
        w = windows(r, e, jd)
        s = stat(w)
        s.update({"id": tid, "pos": pos or "C", "nep": int(jd.size), "scale": scale,
                  "jd_first": float(jd.min()) if jd.size else np.nan, "jd_last": float(jd.max()) if jd.size else np.nan,
                  "med_e": float(np.median(e)) if e.size else np.nan,
                  "med_win_err": float(np.median(w.ferr)) if len(w) else np.nan})
        rows.append(s)
        if pos == "":
            if not KEEP_INFOBITS:
                w.to_csv(DATA / tid / "windows.csv", index=False, float_format="%.5g")
            if verbose:
                print(w.sort_values("S", ascending=False).head(8).to_string())
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    ap.add_argument("--targets", default="lsn_targets.csv", help="target list (lsn_targets.csv, or val_targets.csv for the validation SNe)")
    ap.add_argument("--keep-infobits", action="store_true", help="robustness variant: keep INFOBITS != 0 epochs")
    a = ap.parse_args()
    global KEEP_INFOBITS
    KEEP_INFOBITS = a.keep_infobits
    t = pd.read_csv(HERE / a.targets)
    if a.only:
        t = t[t.id.isin(a.only.split(","))]
    rows = []
    for r in t.itertuples():
        rows += run_target(r.id, verbose=bool(a.only))
    out = pd.DataFrame(rows)
    if a.only:
        print(out.to_string())
        return
    out = out.merge(t[["id", "rank", "class", "survey", "z_l", "z_s", "mu", "n_cc", "n_ia", "rmag"]], on="id")
    out.to_csv(DATA / ("detect_keepinfobits.csv" if KEEP_INFOBITS else "detect.csv"), index=False, float_format="%.6g")
    lens = out[out.pos == "C"]
    ctrl = out[out.pos != "C"]
    print(f"{lens.id.nunique()} targets; lens T: median {lens['T'].median():.2f}, 99.5th pct {lens['T'].quantile(0.995):.2f}, "
          f"max {lens['T'].max():.2f}; control T: median {ctrl['T'].median():.2f}, 99.5th pct {ctrl['T'].quantile(0.995):.2f}, "
          f"max {ctrl['T'].max():.2f} (n={len(ctrl)})")
    thr = max(5.0, ctrl["T"].quantile(0.995))
    cand = lens[lens["T"] > thr].sort_values("T", ascending=False)
    print(f"threshold {thr:.2f}: {len(cand)} candidates")
    print(cand[["id", "class", "T", "T2", "T_t0", "T_n", "nep", "scale", "med_win_err"]].to_string())


if __name__ == "__main__":
    sys.exit(main())
