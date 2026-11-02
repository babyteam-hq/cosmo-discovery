"""The table (R-130, cosmo-deep-10): per lens class, N lenses, the injection efficiency, the candidates, the forecast
expectation and the detections/forecast ratio or its 95 % upper limit, as fixed in results-R130.md.

E_class = sum over lenses of (N_Ia,i eps_Ia,i + N_cc,i eps_cc,i) x T_i, with T_i the injection window in years (data span
+ 200 d), eps from the lens's own injections, N from Shu et al. 2018 for class S (class medians for the 18 without a rate)
and the SLACS-like medians (N_Ia 0.004, N_cc 0.030 per year) as the proxy for class L, with the BELLS-like medians
(0.0177, 0.149) as the high bracket. Poisson 95 % intervals on the detected count (Garwood): upper 3.00 for zero.

Usage (from research/cosmo-discovery):  python ztf_lsn_rate.py [--detections ID,ID]   # IDs that survived vetting
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "ztf_lsn"
PROXY = {"SLACS-like": (0.0040, 0.030), "BELLS-like": (0.0177, 0.149)}


def poisson_interval(k: int, cl: float = 0.95) -> tuple[float, float]:
    lo = 0.0 if k == 0 else stats.chi2.ppf((1 - cl) / 2, 2 * k) / 2
    hi = stats.chi2.ppf(1 - (1 - cl) / 2, 2 * k + 2) / 2 if k > 0 else -np.log(1 - cl)
    return lo, hi


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--detections", default="")
    a = ap.parse_args()
    det = set(a.detections.split(",")) - {""}
    t = pd.read_csv(HERE / "lsn_targets.csv")
    d = pd.read_csv(DATA / "detect.csv")
    inj = pd.read_csv(DATA / "inject_summary.csv")
    lens = d[d.pos == "C"].merge(inj.drop(columns=["class"]), on="id", how="inner")
    lens = lens.merge(t[["id"]], on="id")
    ctrl = d[d.pos != "C"]
    thr = max(5.0, ctrl["T"].quantile(0.995))
    lens["years"] = (lens.jd_last - lens.jd_first + 200) / 365.25
    lens["candidate"] = lens["T"] > thr
    lens["detection"] = lens.id.isin(det)
    # rates
    s_rates = lens[(lens["class"] == "S") & lens.n_ia.notna()]
    med_ia, med_cc = s_rates.n_ia.median(), s_rates.n_cc.median()
    lens["r_ia"] = lens.n_ia.where(lens["class"] == "S", np.nan).fillna(med_ia)
    lens["r_cc"] = lens.n_cc.where(lens["class"] == "S", np.nan).fillna(med_cc)
    is_l = lens["class"].isin(["L", "X"])
    lens.loc[is_l, "r_ia"], lens.loc[is_l, "r_cc"] = PROXY["SLACS-like"]
    lens["eff_Ia"] = lens.eff_Ia.fillna(lens.eff)
    lens["eff_CC"] = lens.eff_CC.fillna(lens.eff)
    lens["E"] = (lens.r_ia * lens.eff_Ia + lens.r_cc * lens.eff_CC) * lens.years
    lens["E_hi"] = lens.E
    lens.loc[is_l, "E_hi"] = (PROXY["BELLS-like"][0] * lens.eff_Ia + PROXY["BELLS-like"][1] * lens.eff_CC)[is_l] * lens.years[is_l]
    lens["N_source"] = (lens.r_ia + lens.r_cc) * lens.years  # lensed SNe that occurred, forecast, any brightness
    print(f"candidate threshold {thr:.2f} (controls: n={len(ctrl)}, 99.5th pct {ctrl['T'].quantile(0.995):.2f}, max {ctrl['T'].max():.2f})")
    rows = []
    for cls, g in lens.groupby("class"):
        k = int(g.detection.sum())
        lo, hi = poisson_interval(k)
        E = g.E.sum()
        row = {
            "class": cls, "N lenses": len(g), "epochs/lens (median)": int(g.nep.median()),
            "15-d window depth, 5 sigma (median r mag)": round(float(np.median(23.9 - 2.5 * np.log10(5 * g.med_win_err))), 2),
            "eff (all injected)": round(g.eff.mean(), 4), "eff Ia": round(g.eff_Ia.mean(), 4), "eff CC": round(g.eff_CC.mean(), 4),
            "eff at m_peak < 21": round(g.eff_brighter_21.mean(), 3),
            "candidates (T > thr)": int(g.candidate.sum()), "detections": k,
            "lensed SNe forecast to occur (any mag)": round(g.N_source.sum(), 2),
            "E (forecast x efficiency)": round(E, 4), "E high bracket": round(g.E_hi.sum(), 4),
            "detections/E 95 %": f"[{lo / E:.3g}, {hi / E:.3g}]" if k > 0 else f"< {hi / E:.3g}",
            "detections/E_high 95 %": f"< {hi / g.E_hi.sum():.3g}" if k == 0 else f"[{lo / g.E_hi.sum():.3g}, {hi / g.E_hi.sum():.3g}]",
            # forecast-free (added by R-135, a presentation column on the same injections): lensed SNe whose summed
            # images peak brighter than r = 21, per lens per year, from the exposure sum_i T_i eps_i(m_peak < 21)
            "lens-years x eff(m_peak < 21)": round(float((g.years * g.eff_brighter_21.fillna(0)).sum()), 1),
            "rate m_peak < 21 per lens-year, 95 %": (f"< {hi / max((g.years * g.eff_brighter_21.fillna(0)).sum(), 1e-9):.3g}"
                                                     if k == 0 else "see detections"),
        }
        rows.append(row)
    table = pd.DataFrame(rows)
    pd.set_option("display.width", 250)
    print(table.T.to_string())
    table.to_csv(DATA / "rate_table.csv", index=False)
    # efficiency per peak magnitude per class
    cols = [c for c in inj.columns if c.startswith("eff_m")]
    ncols = [c for c in inj.columns if c.startswith("n_m")]
    for cls, g in lens.groupby("class"):
        eff = []
        for c, nc in zip(cols, ncols):
            n = g[nc].sum()
            w = g[nc] / max(n, 1)
            eff.append(f"{c[5:]}: {float((g[c].fillna(0) * w).sum()):.2f} (n={int(n)})")
        print(f"class {cls} efficiency per peak r mag bin (lower edge): " + "; ".join(eff))
    lens.to_csv(DATA / "lens_table.csv", index=False, float_format="%.5g")
    # the per-lens top rows
    print("\ncandidates:")
    print(lens[lens.candidate].sort_values("T", ascending=False)[["id", "class", "T", "T2", "T_t0", "T_n", "nep", "scale", "eff"]].to_string())


if __name__ == "__main__":
    sys.exit(main())
