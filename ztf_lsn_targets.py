"""The fixed lens list of R-130 (cosmo-deep-10): known galaxy-scale strong lenses in the ZTF footprint, in the order
their ZTF difference-image cutouts are fetched. Written and committed before any light curve was extracted (2026-11-01).

Classes:
  S  spectroscopically confirmed galaxy-galaxy lenses with a per-lens lensed-SN rate forecast: the 128-lens compilation
     of Shu et al. 2018 (arXiv:1803.07569, Table 1: SLACS, S4TM, BELLS; z_L, z_S, total magnification mu, SFR, N_cc and
     N_Ia per year of the lensed source). Dec > -30. Ordered by N_cc + N_Ia descending, lenses without a rate last.
  L  DESI Legacy Surveys grade-A strong-lens candidates: Huang et al. 2020 (ApJ 894, 78; VizieR J/ApJ/894/78/lens, Q=A),
     Huang et al. 2021 (ApJ 909, 27; J/ApJ/909/27/cand, Q=A), Storfer et al. 2022 (arXiv:2206.02764; the project site's
     CSV, grade A). Dec > -30, not within 3" of a class-S lens, duplicates within 2" merged. Ordered by lens r magnitude
     ascending (the brighter the lens galaxy the larger its Einstein radius and magnification, on average).
  X  the one lensed-SN candidate of Sheu et al. 2023 (arXiv:2301.03578) inside the ZTF footprint, DESI-035.1374+00.4676
     (grade C there), as a single check target.

Download order: X, then S and L interleaved one for one, so a download budget that ends early leaves both classes sampled
by their own priority order. The subset actually fetched is whatever prefix of this order the budget allowed; it is
reported, never chosen after the fact.

Inputs under data/ztf_lsn/inputs/ (gitignored; how to refetch is in results-R130.md). Output: lsn_targets.csv (committed).

Usage (from research/cosmo-discovery):  python ztf_lsn_targets.py
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
INP = HERE / "data" / "ztf_lsn" / "inputs"
DEC_MIN = -30.0


def hms(s: str) -> float:
    h, m, sec = s.split(":")
    return 15.0 * (float(h) + float(m) / 60 + float(sec) / 3600)


def dms(s: str) -> float:
    sign = -1.0 if s.startswith("-") else 1.0
    d, m, sec = s.lstrip("+-").split(":")
    return sign * (float(d) + float(m) / 60 + float(sec) / 3600)


def num(s: str) -> float:
    s = s.strip()
    if s.startswith("...") or s == "":
        return np.nan
    return float(s.split("$\\pm$")[0].strip().rstrip("\\%"))


def shu2018() -> pd.DataFrame:
    rows = []
    for line in (INP / "shu2018_table1.tex").read_text().splitlines():
        cells = [c.strip() for c in line.rstrip("\\ \n").split("&")]
        if len(cells) < 10:
            continue
        name = re.sub(r"\\phantom\{B\}|\\,|\$|\^\{\*\}|\{|\}", "", cells[0]).replace("$", "")
        pos = cells[1].replace("$", "")
        m = re.match(r"(\d\d:\d\d:\d\d\.\d+)([+-]\d\d:\d\d:\d\d\.\d+)", pos)
        ra, dec = hms(m.group(1)), dms(m.group(2))
        rows.append({
            "id": name, "class": "S", "ra": ra, "dec": dec, "z_l": num(cells[2]), "z_s": num(cells[3]),
            "mu": num(cells[4]), "f_fiber": num(cells[5]) / 100, "sfr": num(cells[7]),
            "n_cc": num(cells[8]), "n_ia": num(cells[9]), "survey": name[:5].rstrip("J"), "rmag": np.nan,
        })
    df = pd.DataFrame(rows)
    df["rate"] = df.n_cc.fillna(0) + df.n_ia.fillna(0)
    df.loc[df.n_cc.isna(), "rate"] = -1
    return df.sort_values("rate", ascending=False).drop(columns="rate")


def vizier(path: Path) -> pd.DataFrame:
    lines = [l for l in path.read_text().splitlines() if l and not l.startswith("#")]
    head = lines[0].split("\t")
    body = [l.split("\t") for l in lines[3:]]
    return pd.DataFrame(body, columns=head)


def legacy() -> pd.DataFrame:
    h20 = vizier(INP / "J_ApJ_894_78_lens.tsv")
    h20 = h20[h20.Q.str.strip() == "A"]
    a = pd.DataFrame({"id": h20.Name.str.strip(), "ra": h20._RA.astype(float), "dec": h20._DE.astype(float),
                      "rmag": pd.to_numeric(h20.rmag, errors="coerce"), "z_l": pd.to_numeric(h20.z, errors="coerce"),
                      "survey": "Huang2020"})
    h21 = vizier(INP / "J_ApJ_909_27_cand.tsv")
    h21 = h21[h21.Q.str.strip() == "A"]
    zl = pd.to_numeric(h21.zsp, errors="coerce").fillna(pd.to_numeric(h21.zph, errors="coerce"))
    b = pd.DataFrame({"id": h21.Name.str.strip(), "ra": h21._RA.astype(float), "dec": h21._DE.astype(float),
                      "rmag": pd.to_numeric(h21.rmag, errors="coerce"), "z_l": zl, "survey": "Huang2021"})
    st = pd.read_csv(INP / "storfer2022.csv")
    st = st[st.grade == "A"]
    zph = st["phot-z"].astype(str).str.split("±").str[0]
    zl = pd.to_numeric(st["spec-z"], errors="coerce").fillna(pd.to_numeric(zph, errors="coerce"))
    c = pd.DataFrame({"id": st.name, "ra": st.ra.astype(float), "dec": st.dec.astype(float),
                      "rmag": pd.to_numeric(st.mag_r, errors="coerce"), "z_l": zl, "survey": "Storfer2022"})
    df = pd.concat([a, b, c], ignore_index=True)
    df["class"] = "L"
    df = df[df.dec > DEC_MIN].sort_values("rmag").reset_index(drop=True)
    # merge duplicates within 2"
    keep = []
    for i, r in df.iterrows():
        d = np.hypot((df.loc[keep, "ra"] - r.ra) * np.cos(np.radians(r.dec)), df.loc[keep, "dec"] - r.dec) * 3600
        if len(keep) == 0 or d.min() > 2.0:
            keep.append(i)
    return df.loc[keep].reset_index(drop=True)


def main() -> None:
    s = shu2018()
    s = s[s.dec > DEC_MIN].reset_index(drop=True)
    l = legacy()
    d = np.hypot((l.ra.values[:, None] - s.ra.values[None, :]) * np.cos(np.radians(l.dec.values[:, None])),
                 l.dec.values[:, None] - s.dec.values[None, :]) * 3600
    l = l[d.min(axis=1) > 3.0].reset_index(drop=True)
    x = pd.DataFrame([{"id": "DESI-035.1374+00.4676", "class": "X", "ra": 35.1374, "dec": 0.4676, "z_l": 0.269,
                       "survey": "Sheu2023", "rmag": np.nan}])
    order = [x]
    for i in range(max(len(s), len(l))):
        if i < len(s):
            order.append(s.iloc[[i]])
        if i < len(l):
            order.append(l.iloc[[i]])
    out = pd.concat(order, ignore_index=True)
    out.insert(0, "rank", np.arange(len(out)))
    cols = ["rank", "id", "class", "survey", "ra", "dec", "z_l", "z_s", "mu", "f_fiber", "sfr", "n_cc", "n_ia", "rmag"]
    out = out.reindex(columns=cols)
    out.to_csv(HERE / "lsn_targets.csv", index=False, float_format="%.6g")
    print(f"class S: {len(s)} (with rates {s.n_cc.notna().sum()}; surveys {s.survey.value_counts().to_dict()})")
    print(f"class L: {len(l)} ({l.survey.value_counts().to_dict()}); X: {len(x)}; total {len(out)}")
    print(f"class S rates: N_Ia sum {s.n_ia.sum():.3f}/yr, N_cc sum {s.n_cc.sum():.3f}/yr; medians {s.n_ia.median():.4f}, {s.n_cc.median():.3f}")


if __name__ == "__main__":
    main()
