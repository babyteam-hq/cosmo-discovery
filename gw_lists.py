"""Build the per-run candidate lists for the GW capture-recapture test (R-126, cosmo-deep-6).

One row per candidate, with each search's FAR (1/yr), network SNR and p_astro, from the public catalogue papers'
own tables (arXiv LaTeX sources and PDFs under data/gw/src/) and the 4-OGC search output (data/gw/4ogc/):

  O1-O2 : GWTC-1 (PyCBC, GstLAL, cWB; Tables I-II of arXiv:1811.12907), 4-OGC (arXiv:2112.06878, the HDF release),
          IAS O1 (arXiv:1902.10341, 1902.10331) and O2 (arXiv:1904.07214)
  O3    : GWTC-3's event macros (arXiv:2111.03606: cWB, GstLAL, MBTA, PyCBC, PyCBC-BBH for every O3 candidate with
          p_astro > 0.5 and the O3b marginal list), GWTC-2.1's (arXiv:2108.01045) for the O3a candidates GWTC-3 does
          not carry, 4-OGC, IAS O3a (arXiv:2201.02252), IAS O3b (arXiv:2311.06061), IAS-HM (arXiv:2312.06631)
  O4a   : GWTC-4.0 Tables 1-3 (arXiv:2508.18082, PDF text): cWB-BBH, GstLAL, MBTA, PyCBC

Output: data/gw/lists_<run>.csv with columns name, gps_day, inst, run, mchirp (detector frame where the catalogue
gives one, else source frame x (1+z)), bbh (bool), and per list L in LISTS[run]: L_far, L_snr, L_pastro (NaN when
the list did not report the candidate).

Usage: venv/bin/python gw_lists.py            (writes the CSVs and prints the per-run tallies)
"""
from __future__ import annotations

import math
import os
import re
import sys

import h5py
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "gw")
SRC = os.path.join(DATA, "src")

NAN = float("nan")


# ----------------------------------------------------------------------------------------------------------------
# LaTeX macro files: \DeclareRobustCommand{\NAME}[1]{\IfEqCase{#1}{{key}{value}{key}{value}...}}
# ----------------------------------------------------------------------------------------------------------------
def _balanced(s: str, i: int) -> tuple[str, int]:
    """s[i] == '{'; return (content, index after the matching '}')."""
    assert s[i] == "{", (i, s[i : i + 20])
    depth = 0
    j = i
    while j < len(s):
        c = s[j]
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return s[i + 1 : j], j + 1
        j += 1
    raise ValueError("unbalanced")


def parse_macros(path: str) -> dict[str, dict[str, str]]:
    text = open(path, encoding="utf-8", errors="replace").read()
    out: dict[str, dict[str, str]] = {}
    for m in re.finditer(r"\\(?:DeclareRobustCommand|newcommand)\{\\([A-Za-z]+)\}\[1\]\{\\IfEqCase\{#1\}", text):
        name = m.group(1)
        i = m.end()
        body, _ = _balanced(text, i)
        d: dict[str, str] = {}
        k = 0
        while k < len(body):
            if body[k] != "{":
                k += 1
                continue
            key, k = _balanced(body, k)
            while k < len(body) and body[k] != "{":
                k += 1
            if k >= len(body):
                break
            val, k = _balanced(body, k)
            d[key] = val
        out[name] = d
    return out


def tex_number(v: str | None) -> float:
    """A FAR/SNR/p_astro cell: '--', '\\text{--}', '<1.0 \\times 10^{-5}', '0.47', '>0.99', '$0.12$' -> float (NaN for none).
    An upper limit '<x' returns x; '>x' returns x."""
    if v is None:
        return NAN
    s = v.replace("$", "").replace("\\ensuremath", "").replace("\\text", "").replace("\\textrm", "")
    s = s.replace("\\num", "").replace("{", "").replace("}", "").replace("\\,", "").replace(" ", "").replace("\\%", "")
    s = s.replace("\\leq", "<").replace("\\le", "<").replace("\\geq", ">").replace("\\ge", ">").replace("\\lt", "<")
    s = s.replace("<", "").replace(">", "")
    s = re.sub(r"\\d?dagger|\\ast|\*|\^$|\^\\d?dagger", "", s).rstrip("^")
    s = re.sub(r"\\times10\^(-?\d+)", r"e\1", s)  # 1.0\times10^-5 -> 1.0e-5
    s = re.sub(r"\^\+?[\d.]+", "", s)  # 0.13^{+0.02}_{-0.01} -> 0.13
    s = re.sub(r"_-?[\d.]+", "", s)
    if s in ("", "--", "-", "---", "~", "\\bf", "nan", "Veto", "veto"):
        return NAN
    m = re.match(r"^([0-9.]+)(?:\\times10\^(-?[0-9]+))?$", s)
    if m:
        x = float(m.group(1))
        if m.group(2):
            x *= 10 ** int(m.group(2))
        return x
    m = re.match(r"^([0-9.]+)e(-?[0-9]+)$", s)
    if m:
        return float(m.group(1)) * 10 ** int(m.group(2))
    try:
        return float(s)
    except ValueError:
        return NAN


def pm_median(v: str) -> float:
    """'$52_{-12}^{+17}$' or '43.8+... ' -> 52."""
    s = v.replace("$", "").replace("\\!", "").strip()
    m = re.match(r"^\s*([0-9.]+)", s)
    return float(m.group(1)) if m else NAN


def mchirp_from_m1_q(m1: float, q: float) -> float:
    m2 = q * m1
    return (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2


def mchirp_from_m1_m2(m1: float, m2: float) -> float:
    return (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2


def short_name(full: str) -> str:
    """'GW190403\\_051519' / 'GW190412' / '191118N' -> 'GW190403_051519' (GWYYMMDD_hhmmss when known)."""
    s = full.replace("\\_", "_").replace("$", "").strip()
    s = s.replace("\\bfseries", "").replace("\\bf", "").strip()
    return s


def day_of(name: str) -> str:
    m = re.search(r"(\d{6})", name)
    return m.group(1) if m else ""


def gps_from_name(name: str) -> float | None:
    """GWYYMMDD_hhmmss -> GPS seconds (to the second; names are UTC)."""
    m = re.match(r"^(?:GW)?(\d{2})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2})", name)
    if not m:
        return None
    from datetime import datetime, timezone

    yy, mo, dd, hh, mi, ss = (int(x) for x in m.groups())
    t = datetime(2000 + yy, mo, dd, hh, mi, ss, tzinfo=timezone.utc)
    # GPS = UTC - 315964800 + leap seconds (18 since 2017-01-01, 17 in 2015-2016)
    unix = t.timestamp()
    leap = 18 if unix >= 1483228800 else 17
    return unix - 315964800 + leap


# ----------------------------------------------------------------------------------------------------------------
# LVK O3: GWTC-3 and GWTC-2.1 macros
# ----------------------------------------------------------------------------------------------------------------
LVK_O3_PIPES = {"cwb": "CWBALLSKY", "gstlal": "GSTLALALLSKY", "mbta": "MBTAALLSKY", "pycbc": "PYCBCALLSKY", "pycbcbbh": "PYCBCHIGHMASS"}


def lvk_o3() -> pd.DataFrame:
    g3 = parse_macros(os.path.join(SRC, "2111.03606", "event_macros.tex"))
    g21 = parse_macros(os.path.join(SRC, "2108.01045", "event_macros.tex"))
    g21m = parse_macros(os.path.join(SRC, "2108.01045", "marginal_macros.tex"))
    pe21 = parse_macros(os.path.join(SRC, "2108.01045", "data", "pe", "PRD_version", "Tables_macros_funStuff_PE_vPRD", "PE_macros_PE_vPRD_AllEvents.tex"))
    pe3 = parse_macros(os.path.join(SRC, "2111.03606", "PE_macros.tex"))
    g2 = parse_macros(os.path.join(SRC, "2010.14527", "event_macros.tex"))
    rows = {}

    def add(key, macros, suffix=""):
        full = macros["FULLNAME" + suffix].get(key) if ("FULLNAME" + suffix) in macros else None
        if full is None:
            full = macros["NAME" + suffix].get(key)
        name = short_name(full)
        if not re.match(r"^GW\d{6}_\d{6}$", name):
            # GWTC-3 marginal names like '200121\_031748' (no GW prefix) and short names 'GW190412'
            if re.match(r"^\d{6}_\d{6}$", name):
                name = "GW" + name
            elif re.match(r"^GW\d{6}$", name) and ("FULLNAME" + suffix) in macros:
                name = short_name(macros["FULLNAME" + suffix][key])
        r = rows.setdefault(name, {"name": name})
        r["inst"] = macros["INSTRUMENTS" + suffix].get(key, "")
        for pe in (pe21, pe3):
            mc = tex_number(pe.get("chirpmassdetmed", {}).get(key))
            if not math.isnan(mc):
                r["mchirp_pe"] = mc
        for pipe, pre in LVK_O3_PIPES.items():
            far = macros.get(pre + "FAR" + suffix, {}).get(key)
            if far is None:
                continue
            f = tex_number(far)
            if math.isnan(f):
                continue
            r[pipe + "_far"] = f
            r[pipe + "_snr"] = tex_number(macros.get(pre + "SNR" + suffix, {}).get(key))
            r[pipe + "_pastro"] = tex_number(macros.get(pre + "PASTRO" + suffix, {}).get(key))
            mc = tex_number(macros.get(pre + "CHIRPMASS" + suffix, {}).get(key))
            if not math.isnan(mc):
                r.setdefault("mchirp_pipes", []).append(mc)
            pbbh = tex_number(macros.get(pre + "PBBH" + suffix, {}).get(key))
            if not math.isnan(pbbh):
                r.setdefault("pbbh_pipes", []).append(pbbh)

    for key in g3["INSTRUMENTS"]:
        add(key, g3)
    for key in g21["INSTRUMENTS"]:
        add(key, g21)
    for key in g21m["INSTRUMENTSMAR"]:
        add(key, g21m, "MAR")
    # GWTC-2's cWB results for O3a (GWTC-2.1 dropped the cWB columns; GWTC-3 carries cWB for O3b only)
    for key, full in g2["NAME"].items():
        name = short_name(full)
        if not re.match(r"^GW\d{6}_\d{6}$", name):
            continue
        far = tex_number(g2["CWBALLSKYFAR"].get(key))
        if math.isnan(far):
            continue
        r = rows.setdefault(name, {"name": name, "inst": g2["INSTRUMENTS"].get(key, "")})
        if "cwb_far" not in r:
            r["cwb_far"], r["cwb_snr"], r["cwb_pastro"] = far, tex_number(g2["CWBALLSKYSNR"].get(key)), tex_number(g2["CWBALLSKYPASTRO"].get(key))
            r["cwb_src"] = "GWTC-2"
    df = pd.DataFrame(rows.values())
    df["mchirp"] = df.get("mchirp_pipes", pd.Series([None] * len(df))).apply(
        lambda v: float(np.median(v)) if isinstance(v, list) and v else NAN
    )
    if "mchirp_pe" in df.columns:
        df["mchirp"] = df["mchirp_pe"].astype(float).where(df["mchirp_pe"].notna(), df["mchirp"])
    df["pbbh"] = df.get("pbbh_pipes", pd.Series([None] * len(df))).apply(
        lambda v: float(np.max(v)) if isinstance(v, list) and v else NAN
    )
    df = df.drop(columns=[c for c in ("mchirp_pipes", "pbbh_pipes") if c in df.columns])
    df["run"] = np.where(df["name"].str[2:4].astype(int) * 100 + df["name"].str[4:6].astype(int) <= 1910, "O3a", "O3b")
    return df


# ----------------------------------------------------------------------------------------------------------------
# 4-OGC
# ----------------------------------------------------------------------------------------------------------------
def fourogc(ifar_min: float = 0.02) -> pd.DataFrame:
    f = h5py.File(os.path.join(DATA, "4ogc", "4-OGC_small.hdf"), "r")
    ifar = f["ifar"][:]
    sel = ifar > ifar_min
    d = {k: f[k][:][sel] for k in ("name", "obs", "trig", "ifar", "pastro", "stat", "mass1", "mass2", "time", "H1_snr", "L1_snr", "V1_snr")}
    df = pd.DataFrame(d)
    for c in ("name", "obs", "trig"):
        df[c] = df[c].str.decode("ascii")
    snr2 = np.zeros(len(df))
    for ifo in ("H1", "L1", "V1"):
        s = df[ifo + "_snr"].to_numpy()
        snr2 += np.where(s > 0, s**2, 0.0)
    df["ogc_snr"] = np.sqrt(snr2)
    df["ogc_far"] = 1.0 / df["ifar"]
    df["ogc_pastro"] = df["pastro"]
    df["ogc_mchirp"] = mchirp_from_m1_m2(df["mass1"].to_numpy(), df["mass2"].to_numpy())
    df["gps"] = df["time"]
    t = df["time"].to_numpy()
    run = np.full(len(df), "", dtype=object)
    for lo, hi, lab in ((1126051217, 1137254417, "O1"), (1164556817, 1187733618, "O2"), (1238166018, 1253977218, "O3a"), (1256655618, 1269363618, "O3b")):
        run[(t >= lo) & (t < hi)] = lab
    df["run"] = run
    return df


# ----------------------------------------------------------------------------------------------------------------
# IAS
# ----------------------------------------------------------------------------------------------------------------
def _rows_of_tabular(tex: str, start: int, end: int) -> list[list[str]]:
    """Table rows between two line numbers (1-based, inclusive): split on '&', strip \\rowcolor and comments."""
    lines = tex.splitlines()[start - 1 : end]
    body = []
    for ln in lines:
        ln = re.sub(r"(?<!\\)%.*$", "", ln)
        ln = re.sub(r"\\rowcolor\{[^}]*\}", "", ln)
        body.append(ln)
    text = "\n".join(body)
    rows = []
    for raw in text.split("\\\\"):
        if "&" not in raw:
            continue
        cells = [c.strip() for c in raw.split("&")]
        rows.append(cells)
    return rows


def _clean(cell: str) -> str:
    cell = re.sub(r"\\(tabline|hline|cline\{[^}]*\}|centering|twocolumngrid|egroup|bgroup)", "", cell)
    cell = re.sub(r"\\def\\arraystretch\{[^}]*\}", "", cell)
    lines = [ln.strip() for ln in cell.split("\n") if ln.strip()]
    cell = lines[-1] if lines else ""
    cell = re.sub(r"\\footnote\{.*", "", cell)
    cell = re.sub(r"\$\^\{\\ref\{[^}]*\}\}\$", "", cell)
    cell = re.sub(r"\^\{\\ref\{[^}]*\}\}", "", cell)
    cell = cell.replace("\\textbf{", "").replace("\\texttt{", "")
    return cell.strip().strip("}").strip()


def ias_o3() -> pd.DataFrame:
    """IAS O3a (Olsen+ 2022), O3b (Mehta+ 2023) and IAS-HM (Wadekar+ 2023): one list, the best p_astro / IFAR of the
    three searches per candidate. IFAR for O3a/O3b is per bank (no trials factor, as published); HM's 'overall' IFAR
    includes the trials factor."""
    rows: dict[str, dict] = {}

    def put(name, ifar, pastro, mchirp, src, snr=NAN):
        name = short_name(name)
        if not re.match(r"^GW\d{6}_\d{6}$", name):
            return
        r = rows.setdefault(name, {"name": name, "ias_far": NAN, "ias_pastro": NAN, "ias_mchirp": NAN, "ias_src": "", "ias_snr": NAN})
        if not math.isnan(snr) and math.isnan(r["ias_snr"]):
            r["ias_snr"] = snr
        if not math.isnan(ifar) and (math.isnan(r["ias_far"]) or 1 / ifar < r["ias_far"]):
            r["ias_far"] = 1 / ifar
        if not math.isnan(pastro) and (math.isnan(r["ias_pastro"]) or pastro > r["ias_pastro"]):
            r["ias_pastro"] = pastro
        if not math.isnan(mchirp) and math.isnan(r["ias_mchirp"]):
            r["ias_mchirp"] = mchirp
        r["ias_src"] += src

    # O3a new events (Olsen+): Name, Bank, m1, m2, chi_eff, z, lnL, rhoH2, rhoL2, IFAR(yr), pastro
    tex = open(os.path.join(SRC, "2201.02252", "main.tex"), encoding="utf-8", errors="replace").read()
    for c in _rows_of_tabular(tex, 316, 340):
        if len(c) >= 11 and _clean(c[0]).startswith("GW"):
            m1, m2, z = pm_median(c[2]), pm_median(c[3]), pm_median(c[5])
            put(_clean(c[0]), tex_number(c[9]), tex_number(c[10]), mchirp_from_m1_m2(m1, m2) * (1 + z), "a",
                math.sqrt(tex_number(c[7]) + tex_number(c[8])))
    # O3a known events (extra_tables): Name, Bank, rhoH2, rhoL2, pastro IAS|GWTC-2.1|3-OGC, IFAR IAS|GWTC-2.1|3-OGC
    tex = open(os.path.join(SRC, "2201.02252", "extra_tables.tex"), encoding="utf-8", errors="replace").read()
    for c in _rows_of_tabular(tex, 1, 80):
        if len(c) >= 10 and "GW" in c[0]:
            put(_clean(c[0]), tex_number(c[7]), tex_number(c[4]), NAN, "A", math.sqrt(tex_number(c[2]) + tex_number(c[3])))
    # O3b new events (Mehta+): Event, Bank, m1src, m2src, chieff, z, lnL, rhoH, rhoL, IFAR within bank, IFAR overall, pastro
    tex = open(os.path.join(SRC, "2311.06061", "main.tex"), encoding="utf-8", errors="replace").read()
    for c in _rows_of_tabular(tex, 283, 310):
        if len(c) >= 12 and "GW" in c[0]:
            m1, m2, z = pm_median(c[2]), pm_median(c[3]), pm_median(c[5])
            put(_clean(c[0]), tex_number(c[9]), tex_number(c[11]), mchirp_from_m1_m2(m1, m2) * (1 + z), "b",
                math.hypot(tex_number(c[7]), tex_number(c[8])))
    # O3b known events: Event, Bank, rhoH, rhoL, rhoNet(GWTC-3), pastro IAS|GWTC-3|4-OGC, IFAR IAS|GWTC-3|4-OGC
    for c in _rows_of_tabular(tex, 527, 600):
        if len(c) >= 11 and "GW" in c[0]:
            put(_clean(c[0]), tex_number(c[8]), tex_number(c[5]), NAN, "B", math.hypot(tex_number(c[2]), tex_number(c[3])))
    # IAS-HM new events: Sr, name, bank, m1src, q, chieff, z, rho2net, dlog10 HM, dlog10 HM+prec, rhoH2, rhoL2, IFAR per bank, overall, pastro
    tex = open(os.path.join(SRC, "2312.06631", "Events.tex"), encoding="utf-8", errors="replace").read()
    for c in _rows_of_tabular(tex, 288, 320):
        if len(c) >= 15 and "GW" in c[1]:
            m1, q, z = pm_median(c[3]), pm_median(c[4]), pm_median(c[6])
            put(_clean(c[1]), tex_number(c[13]), tex_number(c[14]), mchirp_from_m1_q(m1, q) * (1 + z), "h",
                math.sqrt(tex_number(c[10]) + tex_number(c[11])))
    # IAS-HM, the (2,2) IAS events it also finds: Event, Bank, m1src, q, chieff, z, rho2net, dlog HM, dlog HM+P, rhoH2, rhoL2, IFAR per bank, overall, pastro
    for c in _rows_of_tabular(tex, 386, 408):
        if len(c) >= 14 and "GW" in c[0]:
            m1, q, z = pm_median(c[2]), pm_median(c[3]), pm_median(c[5])
            put(_clean(c[0]), tex_number(c[12]), tex_number(c[13]), mchirp_from_m1_q(m1, q) * (1 + z), "H",
                math.sqrt(tex_number(c[9]) + tex_number(c[10])))
    # IAS-HM, GWTC-3 events: name, bank, rhoH2, rhoL2, rho2net GWTC-3, pastro IAS|GWTC-3|4-OGC, IFAR IAS|GWTC-3|4-OGC
    for c in _rows_of_tabular(tex, 409, 500):
        if len(c) >= 11 and "GW" in c[0]:
            put(_clean(c[0]), tex_number(c[8]), tex_number(c[5]), NAN, "G", math.sqrt(tex_number(c[2]) + tex_number(c[3])))
    return pd.DataFrame(rows.values())


def name_from_gps(gps: float) -> str:
    from datetime import datetime, timezone

    leap = 18 if gps >= 1167264018 else 17
    t = datetime.fromtimestamp(gps + 315964800 - leap, tz=timezone.utc)
    return "GW" + t.strftime("%y%m%d_%H%M%S")


def ias_o1o2() -> pd.DataFrame:
    rows: dict[str, dict] = {}

    def put(name, ifar, pastro, mchirp, src, snr=NAN):
        r = rows.setdefault(name, {"name": name, "ias_far": NAN, "ias_pastro": NAN, "ias_mchirp": NAN, "ias_src": "", "ias_snr": NAN})
        if not math.isnan(snr) and math.isnan(r["ias_snr"]):
            r["ias_snr"] = snr
        if not math.isnan(ifar) and (math.isnan(r["ias_far"]) or 1 / ifar < r["ias_far"]):
            r["ias_far"] = 1 / ifar
        if not math.isnan(pastro) and (math.isnan(r["ias_pastro"]) or pastro > r["ias_pastro"]):
            r["ias_pastro"] = pastro
        if not math.isnan(mchirp) and math.isnan(r["ias_mchirp"]):
            r["ias_mchirp"] = mchirp
        r["ias_src"] += src

    # O1 (Venumadhav+ 2019b, arXiv:1902.10341, Table at l. 898; 1 O1 = 46 days): IFAR in O1 units, p_astro
    O1_YR = 46 / 365.25
    for name, mc, ifar_o1, pastro in (("GW150914_095045", 28, 20000, 1.0), ("GW151012_095443", 18, 20000, 0.9998),
                                      ("GW151226_033853", 9.74, 20000, 1.0), ("GW151216_094428", 29, 52, 0.71),
                                      ("GW151231_123952", 30, 0.98, 0.15), ("GW151011_192806", 58, 1.1, 0.14)):
        put(name, ifar_o1 * O1_YR, pastro, mc, "O1")
    # O2 (Venumadhav+ 2019). LVC events table (336-360): Name, Bank, GPS, rhoH2, rhoL2, IFAR (O2 units ~ 0.32 yr), (W), pastro
    tex = open(os.path.join(SRC, "1904.07214", "main.tex"), encoding="utf-8", errors="replace").read()
    O2_YR = 118 / 365.25  # the IAS O2 analysis time (118 days of coincident data), their IFAR unit
    for c in _rows_of_tabular(tex, 336, 360):
        if len(c) >= 7 and _clean(c[0]).startswith("GW"):
            put(name_from_gps(tex_number(c[2])), tex_number(c[5]) * O2_YR, tex_number(c[7]) if len(c) > 7 else NAN, NAN, "L",
                math.sqrt(tex_number(c[3]) + tex_number(c[4])))
    # new events (401-430): Name, Bank, Mdet, chieff, z, GPS, rhoH2, rhoL2, IFAR (O2), W, pastro
    for c in _rows_of_tabular(tex, 401, 430):
        if len(c) >= 11 and _clean(c[0]).startswith("GW"):
            put(name_from_gps(tex_number(c[5])), tex_number(c[8]) * O2_YR, tex_number(c[10]), pm_median(c[2]), "N",
                math.sqrt(tex_number(c[6]) + tex_number(c[7])))
    return pd.DataFrame(rows.values())


# ----------------------------------------------------------------------------------------------------------------
# GWTC-1 (O1-O2), from the PDF text
# ----------------------------------------------------------------------------------------------------------------
GWTC1_EVENTS = {
    # name: (UTC day, PyCBC FAR, GstLAL FAR, cWB FAR, PyCBC SNR, GstLAL SNR, cWB SNR, Mchirp source, z) Tables I and III
    "GW150914": (1.53e-5, 1.0e-7, 1.63e-4, 23.6, 24.4, 25.2, 28.6, 0.09),
    "GW151012": (0.17, 7.92e-3, NAN, 9.5, 10.0, NAN, 15.2, 0.21),
    "GW151226": (1.69e-5, 1.0e-7, 0.02, 13.1, 13.1, 11.9, 8.9, 0.09),
    "GW170104": (1.37e-5, 1.0e-7, 2.91e-4, 13.0, 13.0, 13.0, 21.4, 0.20),
    "GW170608": (3.09e-4, 1.0e-7, 1.44e-4, 15.4, 14.9, 14.1, 7.9, 0.07),
    "GW170729": (1.36, 0.18, 0.02, 9.8, 10.8, 10.2, 35.4, 0.49),
    "GW170809": (1.45e-4, 1.0e-7, NAN, 12.2, 12.4, NAN, 24.9, 0.20),
    "GW170814": (1.25e-5, 1.0e-7, 2.08e-4, 16.3, 15.9, 17.2, 24.1, 0.12),
    "GW170817": (1.25e-5, 1.0e-7, NAN, 30.9, 33.0, NAN, 1.186, 0.01),
    "GW170818": (NAN, 4.20e-5, NAN, NAN, 11.3, NAN, 26.5, 0.21),
    "GW170823": (3.29e-5, 1.0e-7, 2.14e-3, 11.1, 11.5, 10.8, 29.2, 0.35),
}
# Table II marginal triggers: day, pipeline, FAR, SNR, detector-frame chirp mass
GWTC1_MARGINAL = [
    ("151008", "pycbc", 10.17, 8.8, 5.12), ("151012.2", "gstlal", 8.56, 9.6, 2.01), ("151116", "pycbc", 4.77, 9.0, 1.24),
    ("161202", "gstlal", 6.00, 10.5, 1.54), ("161217", "gstlal", 10.12, 10.7, 7.86), ("170208", "gstlal", 11.18, 10.0, 7.39),
    ("170219", "gstlal", 6.26, 9.6, 1.53), ("170405", "gstlal", 4.55, 9.3, 1.44), ("170412", "gstlal", 8.22, 9.7, 4.36),
    ("170423", "gstlal", 6.47, 8.9, 1.17), ("170616", "pycbc", 1.94, 9.1, 2.75), ("170630", "gstlal", 10.46, 9.7, 0.90),
    ("170705", "gstlal", 10.97, 9.3, 3.40), ("170720", "gstlal", 10.75, 13.0, 5.96),
]
GWTC1_TIMES = {  # UTC hhmmss of the 11 events (Table I)
    "GW150914": "095045", "GW151012": "095443", "GW151226": "033853", "GW170104": "101158", "GW170608": "020116",
    "GW170729": "185629", "GW170809": "082821", "GW170814": "103043", "GW170817": "124104", "GW170818": "022509",
    "GW170823": "131358",
}


def lvk_o1o2() -> pd.DataFrame:
    rows = []
    for name, (fp, fg, fc, sp, sg, sc, mc, z) in GWTC1_EVENTS.items():
        rows.append({
            "name": name + "_" + GWTC1_TIMES[name], "short": name, "inst": "HL" if name[2:4] == "15" or name[2:6] < "1708" else "HLV",
            "pycbc_far": fp, "gstlal_far": fg, "cwb_far": fc, "pycbc_snr": sp, "gstlal_snr": sg, "cwb_snr": sc,
            "mchirp": mc * (1 + z), "lvk_marginal": False,
        })
    for day, pipe, far, snr, mc in GWTC1_MARGINAL:
        r = {"name": "LVT" + day, "short": "LVT" + day, "inst": "HL", "mchirp": mc, "lvk_marginal": True}
        r[pipe + "_far"] = far
        r[pipe + "_snr"] = snr
        rows.append(r)
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------------------------------------
# GWTC-4.0 (O4a) from the PDF text: Tables 1, 2 (per-pipeline) and 3 (chirp masses)
# ----------------------------------------------------------------------------------------------------------------
_CELL = r"(⋯|<?\d[\d.]*(?:\s*×\s*10[−-]\d+)?|>0\.99|<0\.01)"


def _num4(s: str) -> float:
    s = s.replace("⋯", "").replace("<", "").replace(">", "").replace(" ", "")
    if not s:
        return NAN
    m = re.match(r"^([\d.]+)(?:×10[−-](\d+))?$", s)
    if not m:
        return NAN
    x = float(m.group(1))
    if m.group(2):
        x *= 10 ** (-int(m.group(2)))
    return x


def lvk_o4a() -> pd.DataFrame:
    txt = open(os.path.join(SRC, "2508.18082.txt"), encoding="utf-8", errors="replace").read()
    lines = txt.splitlines()
    rows = {}
    pat = re.compile(r"^\s*(GW\d{6}_\d{6})\*?\s+(H|L|HL|HLV|LV|HV)\s+(.*)$")
    for i, ln in enumerate(lines[1100:1420], start=1100):
        m = pat.match(ln)
        if not m:
            continue
        name, inst, rest = m.group(1), m.group(2), m.group(3)
        toks = re.findall(r"<?\d[\d.]*\s*×\s*10[−-]\d+|⋯|<?>?[\d.]+", rest)
        toks = [t.replace(" ", "") for t in toks]
        # 4 pipelines x (FAR, SNR, pastro); a missing pastro at a line break leaves 11 tokens: pad per pipeline block
        vals = []
        k = 0
        for p in range(4):
            block = toks[k : k + 3]
            if len(block) == 3 and block[0] == "⋯":
                vals.append((NAN, NAN, NAN))
                k += 3
            elif len(block) >= 2 and block[0] != "⋯":
                far, snr = _num4(block[0]), _num4(block[1])
                pas = _num4(block[2]) if len(block) == 3 and not block[2].startswith(("⋯",)) and _num4(block[2]) <= 1.0 and "." in block[2] else NAN
                vals.append((far, snr, pas))
                k += 3 if (len(block) == 3 and not math.isnan(pas)) else 2
            else:
                vals.append((NAN, NAN, NAN))
                k += len(block)
        r = rows.setdefault(name, {"name": name, "inst": inst, "run": "O4a"})
        for (far, snr, pas), pipe in zip(vals, ("cwb", "gstlal", "mbta", "pycbc")):
            if not math.isnan(far):
                r[pipe + "_far"], r[pipe + "_snr"], r[pipe + "_pastro"] = far, snr, pas
    df = pd.DataFrame(rows.values())
    # chirp masses from Table 3 (source frame) and redshift. pdftotext puts each row's '+err' and '-err' on separate
    # lines; a value token is 'v+e' when the plus error stayed inline, else 'v' followed by the minus error 'e'.
    mc = {}
    for ln in lines[1580:1900]:
        m = re.match(r"^\s*(GW\d{6}_\d{6})\s+(.*)$", ln)
        if not m:
            continue
        toks = m.group(2).replace(",", "").split()
        meds = []
        i = 0
        while i < len(toks) and len(meds) < 7:
            t = toks[i]
            if "+" in t:
                meds.append(float(t.split("+")[0]) if t.split("+")[0] else NAN)
                i += 1
            else:
                try:
                    meds.append(float(t))
                except ValueError:
                    meds.append(NAN)
                i += 2
        if len(meds) >= 7:
            mc[m.group(1)] = meds[1] * (1 + meds[6])
    df["mchirp"] = df["name"].map(mc)
    return df


# ----------------------------------------------------------------------------------------------------------------
# Cross-matching and assembly
# ----------------------------------------------------------------------------------------------------------------
def match_by_time(names: pd.Series, gps: np.ndarray, tol: float = 1.0) -> dict[str, int]:
    """Map each LVK-style name (GWYYMMDD_hhmmss) to the index of the nearest 4-OGC trigger within tol seconds."""
    out = {}
    order = np.argsort(gps)
    sg = gps[order]
    for n in names:
        t = gps_from_name(n)
        if t is None:
            continue
        j = np.searchsorted(sg, t)
        best = None
        for k in (j - 1, j):
            if 0 <= k < len(sg) and abs(sg[k] - t) <= tol:
                if best is None or abs(sg[k] - t) < abs(sg[best] - t):
                    best = k
        if best is not None:
            out[n] = int(order[best])
    return out


def assemble(run: str) -> pd.DataFrame:
    ogc_all = fourogc()
    if run == "O3":
        lvk = lvk_o3()
        ias = ias_o3()
        ogc = ogc_all[ogc_all["run"].isin(["O3a", "O3b"])].reset_index(drop=True)
    elif run == "O1O2":
        lvk = lvk_o1o2()
        ias = ias_o1o2()
        ogc = ogc_all[ogc_all["run"].isin(["O1", "O2"])].reset_index(drop=True)
    elif run == "O4a":
        lvk = lvk_o4a()
        return lvk
    else:
        raise ValueError(run)

    # the union of names: LVK, IAS, 4-OGC above the loose reporting threshold (IFAR > 0.5 yr or pastro > 0.5)
    ogc_rep = ogc[(ogc["ifar"] > 0.5) | (ogc["pastro"] > 0.5)].copy()
    ogc_rep["name"] = "GW" + ogc_rep["name"]
    names = set(lvk["name"]) | set(ias["name"]) | set(ogc_rep["name"])
    # 4-OGC name is GW<YYMMDD_hhmmss> from its own time; LVK names can differ by a second: match by GPS time
    ogc_gps = ogc["gps"].to_numpy()
    rows = []
    lvk_i = lvk.set_index("name")
    ias_i = ias.set_index("name") if len(ias) else None
    ogc_match = match_by_time(pd.Series(sorted(names)), ogc_gps, tol=1.5)
    seen_ogc = set()
    for n in sorted(names):
        r = {"name": n, "run": run}
        if n in lvk_i.index:
            for c in lvk_i.columns:
                v = lvk_i.at[n, c]
                if isinstance(v, (float, int, str, bool, np.floating, np.bool_)):
                    r[c] = v
        if ias_i is not None and n in ias_i.index:
            for c in ias_i.columns:
                r[c] = ias_i.at[n, c]
        j = ogc_match.get(n)
        if j is not None:
            seen_ogc.add(j)
            r["ogc_far"], r["ogc_snr"], r["ogc_pastro"], r["ogc_mchirp"] = (
                float(ogc.at[j, "ogc_far"]), float(ogc.at[j, "ogc_snr"]), float(ogc.at[j, "ogc_pastro"]), float(ogc.at[j, "ogc_mchirp"]))
            r["ogc_trig"] = ogc.at[j, "trig"]
        rows.append(r)
    df = pd.DataFrame(rows)
    # IAS-only and 4-OGC-only names that are the same event under a 1-s-different name: merge by GPS time
    df["gps"] = df["name"].map(lambda n: gps_from_name(n) or NAN)
    df = df.sort_values("gps").reset_index(drop=True)
    merged = []
    i = 0
    while i < len(df):
        r = df.iloc[i].to_dict()
        j = i + 1
        while j < len(df) and abs(df.at[j, "gps"] - r["gps"]) <= 2.0:
            other = df.iloc[j].to_dict()
            for k, v in other.items():
                if k in ("name", "gps", "run"):
                    continue
                if (k not in r) or (isinstance(r[k], float) and math.isnan(r[k])) or r[k] in ("", None):
                    r[k] = v
            if not str(r["name"]).startswith("GW") or (str(other["name"]).startswith("GW") and other.get("inst", "") not in ("", None) and r.get("inst", "") in ("", None)):
                r["name"] = other["name"]
            j += 1
        merged.append(r)
        i = j
    df = pd.DataFrame(merged)
    # chirp mass: LVK template median, else 4-OGC template, else IAS
    mc = df.get("mchirp", pd.Series([NAN] * len(df))).astype(float)
    for alt in ("ogc_mchirp", "ias_mchirp"):
        if alt in df.columns:
            mc = mc.where(~mc.isna(), df[alt].astype(float))
    df["mchirp"] = mc
    if "inst" not in df.columns:
        df["inst"] = ""
    df["inst"] = df["inst"].fillna("")
    if "ogc_trig" in df.columns:
        df["inst"] = np.where(df["inst"] == "", df["ogc_trig"].fillna("").str.replace("1", ""), df["inst"])
    if "ias_far" in df.columns:  # the IAS searches use Hanford-Livingston coincident data only
        df["inst"] = np.where((df["inst"] == "") & df["ias_far"].notna(), "HL", df["inst"])
    if run == "O3":
        ym = df["name"].str[2:4].astype(int) * 100 + df["name"].str[4:6].astype(int)
        df["run"] = np.where(ym <= 1910, "O3a", "O3b")
    elif run == "O1O2":
        df["run"] = np.where(df["name"].str.contains(r"15\d{4}"), "O1", "O2")
    return df


def tally(df: pd.DataFrame, lists: list[str], far_thr: float = 1.0) -> str:
    out = []
    for L in lists:
        c = L + "_far"
        if c not in df.columns:
            out.append(f"{L}: none")
            continue
        out.append(f"{L}: {int((df[c] < far_thr).sum())} at FAR<{far_thr}/yr, {int(df[c].notna().sum())} reported")
    return "; ".join(out)


if __name__ == "__main__":
    for run, lists in (("O1O2", ["pycbc", "gstlal", "cwb", "ogc", "ias"]),
                       ("O3", ["cwb", "gstlal", "mbta", "pycbc", "pycbcbbh", "ogc", "ias"]),
                       ("O4a", ["cwb", "gstlal", "mbta", "pycbc"])):
        df = assemble(run)
        path = os.path.join(DATA, f"lists_{run}.csv")
        df.to_csv(path, index=False)
        print(f"{run}: {len(df)} candidates -> {path}")
        print("  " + tally(df, lists))
