"""The limits with the vetted efficiency (R-141, cosmo-deep-10b).

Headline (forecast-free): lensed SNe whose summed images peak brighter than r = m_lim (21, also 20 and 20.5), per lens
per year, 95 %, per class: zero detections over the exposure X = sum_i T_i eps_i(m_peak < m_lim), with eps_i the vetted
efficiency of lens i's injections brighter than m_lim (lenses with fewer than 3 such injections contribute 0) and T_i the
injection window (data span + 200 d). The efficiency's uncertainty enters the limit: the upper limit R95 solves
E_eps[exp(-R95 X)] = 0.05 with X lognormal, its width the quadrature sum of the bootstrap over injections and a systematic
(half the relative change of X between the 1.2" and the 2.0" image offset, `--offset 2.0` runs, i.e. 2.0" taken as
a 2-sigma excursion; the limit with every image at 2.0" is printed beside), so a wider PSF-response or
image-offset uncertainty loosens the limit. Also per 0.5-mag bin of summed peak r (differential limits).
Secondary rows (each conditional on the source population it assumes): E = sum_i (N_Ia eps_Ia,i + N_cc eps_cc,i) T_i with
Shu et al. 2018's per-lens rates (class S), and with the SLACS-like and BELLS-like median rates (class L).

Usage (from research/cosmo-discovery):  python ztf_lsn_rate_vetted.py [--stop RANK] [--detections ID,ID]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

from ztf_lsn_rate import PROXY

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "ztf_lsn"
MLIMS = (20.0, 20.5, 21.0)


def ul_with_uncertainty(X: float, rel: float, k: int = 0, cl: float = 0.95) -> float:
    """One-sided upper limit on R for k events (3.00/X for none, 4.74/X for one) when the exposure X has a lognormal
    relative uncertainty rel (averaged over X)."""
    if X <= 0:
        return np.inf
    z = np.random.default_rng(0).normal(size=20000) if rel > 0 else np.zeros(1)
    xs = X * np.exp(rel * z - rel**2 / 2)
    f = lambda R: np.mean(stats.poisson.cdf(k, R * xs)) - (1 - cl)
    return float(optimize.brentq(f, 1e-9, 1e3 / X))


def exposure(inj: dict[str, pd.DataFrame], years: pd.Series, ids, mlim: float, col: str = "vetted",
             lo: float | None = None) -> float:
    x = 0.0
    for i in ids:
        df = inj[i]
        sel = (df.m_peak < mlim) & ((df.m_peak >= lo) if lo is not None else True)
        if sel.sum() >= 3:
            x += years[i] * df[col][sel].mean()
    return x


def boot_rel(inj, years, ids, mlim, nb: int = 300) -> float:
    rng = np.random.default_rng(141)
    vals = []
    for _ in range(nb):
        x = 0.0
        for i in ids:
            df = inj[i]
            sel = np.flatnonzero((df.m_peak < mlim).values)
            if sel.size >= 3:
                pick = rng.choice(sel, sel.size)
                x += years[i] * df.vetted.values[pick].mean()
        vals.append(x)
    vals = np.array(vals)
    return float(vals.std() / max(vals.mean(), 1e-9))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop", type=int, default=10**6)
    ap.add_argument("--detections", default="")
    ap.add_argument("--detections-chi2", default="", help="survivors under the chi2 colour form (robustness row)")
    a = ap.parse_args()
    det = set(a.detections.split(",")) - {""}
    t = pd.read_csv(HERE / "lsn_targets.csv")
    t = t[t["rank"] < a.stop]
    d = pd.read_csv(DATA / "detect.csv")
    lens = d[d.pos == "C"].drop(columns=["n_ia", "n_cc", "survey", "class"], errors="ignore").merge(
        t[["id", "class", "n_ia", "n_cc", "survey"]].rename(columns={"class": "cls"}), on="id")
    ap_var = {"chi2": "inject_vetted_chi2", "muwide": "inject_vetted_muwide"}
    inj, inj2 = {}, {}
    for i in lens.id:
        p = DATA / "inject_vetted" / f"{i}.csv"
        if p.exists():
            inj[i] = pd.read_csv(p)
        p2 = DATA / "inject_vetted_off2.0" / f"{i}.csv"
        if p2.exists():
            inj2[i] = pd.read_csv(p2)
    lens = lens[lens.id.isin(inj)].copy()
    lens["years"] = (lens.jd_last - lens.jd_first + 200) / 365.25
    years = lens.set_index("id").years
    for c in ("detected", "vetted"):
        lens[f"eff_{c}"] = [inj[i][c].mean() for i in lens.id]
        lens[f"effIa_{c}"] = [inj[i][c][inj[i].type == "Ia"].mean() if (inj[i].type == "Ia").any() else np.nan for i in lens.id]
        lens[f"effCC_{c}"] = [inj[i][c][inj[i].type != "Ia"].mean() for i in lens.id]
    for c in ("effIa_detected", "effIa_vetted"):
        lens[c] = lens[c].fillna(lens[c.replace("Ia", "CC")])
    ctrl = d[(d.pos != "C") & d.id.isin(lens.id)]
    print(f"{len(lens)} lenses with vetted injections; controls n={len(ctrl)}, 99.5th pct {ctrl['T'].quantile(0.995):.2f} "
          f"(threshold used 5.31, R-135's)")
    # stage-by-stage efficiency
    allinj = pd.concat([inj[i].assign(cls=c) for i, c in zip(lens.id, lens.cls)])
    for cls, g in allinj.groupby("cls"):
        b = g[g.m_peak < 21]
        print(f"class {cls}: n_inj {len(g)}; all: det {g.s0.mean():.4f} s1 {g.s1.mean():.4f} s2 {g.s2.mean():.4f} "
              f"s3 {g.s3.mean():.4f} vetted {g.s4.mean():.4f};  m<21 (n={len(b)}): det {b.s0.mean():.3f} s1 {b.s1.mean():.3f} "
              f"s2 {b.s2.mean():.3f} s3 {b.s3.mean():.3f} vetted {b.s4.mean():.3f}")
        bins = np.arange(18.0, 23.01, 0.5)
        cut = pd.cut(g.m_peak, bins, right=False)
        tab = g.groupby(cut, observed=True).agg(n=("s0", "size"), det=("s0", "mean"), vetted=("s4", "mean"))
        print("  per summed peak r: " + "; ".join(f"{iv.left:.1f}: {r.det:.2f}->{r.vetted:.2f} (n={r.n})" for iv, r in tab.iterrows()))
    rows = []
    for cls, g in lens.groupby("cls"):
        ids = list(g.id)
        k = int(g.id.isin(det).sum())
        row = {"class": cls, "N": len(g), "eff det": round(g.eff_detected.mean(), 4), "eff vetted": round(g.eff_vetted.mean(), 4),
               "eff vetted Ia": round(g.effIa_vetted.mean(), 4), "eff vetted CC": round(g.effCC_vetted.mean(), 4),
               "detections": k}
        for m in MLIMS:
            X = exposure(inj, years, ids, m)
            Xd = exposure(inj, years, ids, m, "detected")
            rb = boot_rel(inj, years, ids, m) if m == 21.0 else 0.0
            rs = 0.0
            if inj2 and all(i in inj2 for i in ids):
                X2 = exposure(inj2, years, ids, m)
                rs = 0.5 * abs(X2 - X) / max(X, 1e-9)  # 2.0" taken as a 2-sigma excursion of the image offset
                if m == 21.0:
                    row["X<21 at 2.0in"] = round(X2, 1)
                    row["R95<21 at 2.0in"] = round(ul_with_uncertainty(X2, 0.0, k), 5)
            rel = float(np.hypot(rb, rs))
            row[f"X<{m:g} (lens-yr, vetted)"] = round(X, 1)
            row[f"X<{m:g} unvetted"] = round(Xd, 1)
            row[f"R95<{m:g} no unc"] = round(ul_with_uncertainty(X, 0.0, k), 5)
            if m == 21.0:
                row["rel unc X<21 (boot, offset)"] = f"{rb:.3f}, {rs:.3f}"
                row["R95<21 with unc"] = round(ul_with_uncertainty(X, rel, k), 5)
                row["R95<21 if k=0"] = round(ul_with_uncertainty(X, rel, 0), 5)
        # forecast rows
        g = g.copy()
        if cls == "S":
            s = g[g.n_ia.notna()]
            mi, mc = s.n_ia.median(), s.n_cc.median()
            g["r_ia"], g["r_cc"] = g.n_ia.fillna(mi), g.n_cc.fillna(mc)
            pops = {"Shu 2018 per-lens": (g.r_ia, g.r_cc)}
        else:
            pops = {k_: (pd.Series(v[0], index=g.index), pd.Series(v[1], index=g.index)) for k_, v in PROXY.items()}
        for name, (ri, rc) in pops.items():
            E = float(((ri * g.effIa_vetted + rc * g.effCC_vetted) * g.years).sum())
            Ed = float(((ri * g.effIa_detected + rc * g.effCC_detected) * g.years).sum())
            occ = float(((ri + rc) * g.years).sum())
            hi = ul_with_uncertainty(1.0, 0.0, k)  # one-sided 95 % on the count
            row[f"{name}: occur / E_det / E_vet / ratio95"] = f"{occ:.1f} / {Ed:.3f} / {E:.3f} / < {hi / E:.3g}"
        rows.append(row)
    pd.set_option("display.width", 250, "display.max_colwidth", 80)
    out = pd.DataFrame(rows)
    print(out.T.to_string())
    out.to_csv(DATA / "rate_table_vetted.csv", index=False)
    # robustness rows: the first-frozen chi2 colour test (R-135's J1250 also survives it), the wide-magnification draw
    extra = {"chi2": set(a.detections_chi2.split(",")) - {""}}
    for tag, sub in ap_var.items():
        for cls, g in lens.groupby("cls"):
            have = {i: pd.read_csv(DATA / sub / f"{i}.csv") for i in g.id if (DATA / sub / f"{i}.csv").exists()}
            if len(have) < len(g):
                print(f"variant {tag} class {cls}: {len(have)}/{len(g)} lenses run")
                continue
            X = exposure(have, years, list(g.id), 21.0)
            k = int(g.id.isin(det | extra.get(tag, set())).sum())
            effb = np.mean([have[i].vetted[have[i].m_peak < 21].mean() for i in g.id if (have[i].m_peak < 21).sum() >= 3])
            print(f"variant {tag} class {cls}: X<21 {X:.1f} lens-yr, k={k}, R95<21 {ul_with_uncertainty(X, 0.0, k):.4g}, "
                  f"mean vetted eff(m<21) {effb:.3f}")
    # differential forecast-free limits per bin
    for cls, g in lens.groupby("cls"):
        s = []
        for lo in np.arange(18.5, 21.01, 0.5):
            X = exposure(inj, years, list(g.id), lo + 0.5, lo=lo)
            s.append(f"{lo:.1f}-{lo + 0.5:.1f}: < {(-np.log(0.05) / X) if X > 0 else np.inf:.3g}")
        print(f"class {cls} differential R95 per lens-yr: " + "; ".join(s))
    lens.to_csv(DATA / "lens_table_vetted.csv", index=False, float_format="%.5g")


if __name__ == "__main__":
    sys.exit(main())
