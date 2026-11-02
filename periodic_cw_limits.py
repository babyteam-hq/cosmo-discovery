"""NANOGrav 15-year continuous-wave limits against the periodic-quasar binary candidates (R-124, cosmo-deep-4).

The cheap pass of the brief: each candidate's implied strain h0 (from its observed period, redshift and black-hole mass,
mass ratio and mass uncertainty marginalised) against the NANOGrav 15-year 95 % upper limit at its sky position and
gravitational-wave frequency; for the CRTS candidates also against the targeted limits NANOGrav published for them
(Agarwal et al. 2025, arXiv:2508.16534, Tables 6-8).

Usage (from research/cosmo-discovery, with a venv holding numpy, scipy, astropy, healpy):
    python periodic_cw_limits.py            # the whole pass; writes data/periodic/candidates_cw.csv and prints the summary
    python periodic_cw_limits.py --check    # the convention check against Agarwal's own h95 / Mc95 pairs only

Data (gitignored, data/):
    ng15/cw/15yr_cw_analysis-main/data/15yr_cw_3d_limits_v4.npz  (Zenodo 8067506: the 15-year sky- and frequency-dependent
        95 % strain upper limits, 22 frequency bins 1-29.7 nHz x 192 healpix pixels)
    periodic/agarwal_tables_raw.txt   (pdftotext of arXiv:2508.16534 Appendix B, Tables 6-8: the 114 targets)
    periodic/src_1604.01020/PeriodicQSOsPTF.tex   (Charisi et al. 2016, PTF: tables of 50 candidates, 33 surviving)
    periodic/src_1906.08315/ms_arxiv_v2.tex       (Liu et al. 2019, PS1 MDS: 26 candidates)
    periodic/src_2206.11497/ms.tex                (Chen et al. 2024, ZTF: the 11 rows of Table 1 printed in the source;
        the full machine-readable table is behind OUP and CDS, both unreachable from this cloud)

Conventions: f_gw = 2 / P_obs (circular binary, light-curve period = orbital period), also f_gw = 1 / P_obs (the
period-doubling case); h0 = 2 (G M_cz)^(5/3) (pi f_gw)^(2/3) / (c^4 D_L) with M_cz = (1+z) M_c the redshifted chirp mass;
flat LCDM H0 = 70, Om = 0.3 (Xin et al. 2021's choice); M_c = M q^(3/5) / (1+q)^(6/5) with q log-uniform on [0.1, 1];
log M ~ Normal(catalogue value, catalogue error or 0.4 dex). A candidate is "excluded at its claimed parameters" when
95 % of its h0 draws exceed the upper limit at its position and frequency; the expected number excluded is the sum of the
per-candidate exclusion probabilities. Above 29.7 nHz (the map's last bin) the last bin's map is scaled by f / 27.5 nHz
(white-noise regime, h_UL proportional to f), flagged in the table as an extrapolation.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import healpy as hp
import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
G, C, MSUN = const.G.value, const.c.value, const.M_sun.value
RNG = np.random.default_rng(124)
NDRAW = 20000
SIG_DEFAULT = 0.4  # dex, single-epoch virial masses

MINUS = "\u2212"


def num(s: str) -> float:
    return float(s.replace(MINUS, "-").replace("--", "-"))


def sci(s: str) -> float:
    """'7.24 × 106' -> 7.24e6, '5.38 × 10−15' -> 5.38e-15 (pdftotext glues the exponent to the 10)."""
    m = re.match(r"\s*([\d.]+)\s*×\s*10(" + MINUS + r"?\d+)\s*$", s)
    if not m:
        raise ValueError(s)
    return float(m.group(1)) * 10 ** int(m.group(2).replace(MINUS, "-"))


def hms(s: str) -> float:
    m = re.match(r"\s*(\d+)h\s*(\d+)m\s*([\d.]+)s?", s)
    h, mi, se = (float(x) for x in m.groups())
    return 15 * (h + mi / 60 + se / 3600)


def dms(s: str) -> float:
    s = s.replace(MINUS, "-")
    m = re.match(r"\s*([+-]?)(\d+)d\s*(\d+)m\s*([\d.]+)s?", s)
    sign = -1 if m.group(1) == "-" else 1
    d, mi, se = (float(x) for x in m.groups()[1:])
    return sign * (d + mi / 60 + se / 3600)


def sex(ra: str, dec: str) -> tuple[float, float]:
    """'00:04:58.79', '-02:26:29.54' -> degrees."""
    h, m, s = (float(x) for x in ra.split(":"))
    sign = -1 if dec.strip().startswith("-") else 1
    d, dm, ds = (abs(float(x)) for x in dec.split(":"))
    return 15 * (h + m / 60 + s / 3600), sign * (d + dm / 60 + ds / 3600)


def from_jname(name: str) -> tuple[float, float] | None:
    """J hhmmss.ss +ddmmss.s (SDSS style) -> degrees; None when the name is not of that form."""
    m = re.search(r"J(\d{2})(\d{2})(\d{2}(?:\.\d+)?)([+-])(\d{2})(\d{2})(\d{2}(?:\.\d+)?)", name)
    if not m:
        return None
    h, mi, s, sg, d, dm, ds = m.groups()
    ra = 15 * (float(h) + float(mi) / 60 + float(s) / 3600)
    dec = (float(d) + float(dm) / 60 + float(ds) / 3600) * (-1 if sg == "-" else 1)
    return ra, dec


# ----------------------------------------------------------------------------------------------------------- catalogues
def read_agarwal() -> list[dict]:
    rows = []
    pat = re.compile(
        r"^\s*(?P<name>.+?)\s+(?P<ra>\d{2}h\s*\d{2}m\s*[\d.]+s?)\s+(?P<dec>[" + MINUS + r"+-]\d{1,2}d\s*\d{1,2}m\s*[\d.]+s?)"
        r"\s+(?P<z>[\d.]+)\s+(?P<P>\d+)\s+(?P<f>[\d.]+)\s+(?P<rest>.*)$"
    )
    for line in (DATA / "periodic/agarwal_tables_raw.txt").read_text(encoding="utf-8").splitlines():
        m = pat.match(line)
        if not m:
            continue
        rest = m.group("rest").strip()
        parts = re.findall(r"[\d.]+\s*×\s*10" + MINUS + r"?\d+|" + MINUS + r"|–", rest)
        if len(parts) != 3:
            raise ValueError(line)
        mtot = None if parts[0] in (MINUS, "–") else sci(parts[0])
        name = re.sub(r"(?<=\d)[ab]$", "", m.group("name").strip())
        rows.append(dict(cat="CRTS", name=name, ra=hms(m.group("ra")), dec=dms(m.group("dec")), z=num(m.group("z")),
                         P_obs=num(m.group("P")), f_tab=num(m.group("f")) * 1e-9, logM=None if mtot is None else np.log10(mtot),
                         sig=SIG_DEFAULT, Mc95=sci(parts[1]), h95=sci(parts[2])))
    return rows


def read_charisi() -> list[dict]:
    tex = (DATA / "periodic/src_1604.01020/PeriodicQSOsPTF.tex").read_text(encoding="utf-8", errors="replace")
    def table(label):
        i = tex.index(label)
        start = tex.rfind(r"\hline\hline", 0, i)
        return tex[start:i]
    props, periods = {}, {}
    for line in table(r"\label{Table:CandidatesQuasarProperties}").splitlines():
        if "&" not in line:
            continue
        cells = [c.strip() for c in line.rstrip("\\").split("&")]
        raw = cells[0]
        bold = r"\textbf" in raw
        name = re.sub(r"\\textbf\{(.*?)\}", r"\1", raw)
        name = re.sub(r"\$\^\{\*+\}\$", "", name).strip()
        mass = cells[5]
        mm = re.search(r"\)\s*([\d.]+)\$\\pm\$([\d.]+)", mass)
        if mm:
            logM, sig = float(mm.group(1)), float(mm.group(2))
        else:
            logM, sig = float(re.search(r"([\d.]+)", mass).group(1)), SIG_DEFAULT
        pos = from_jname(name)
        if pos is None:
            pos = (num(cells[1]), num(cells[2]))
        props[name] = dict(cat="PTF", name=name, ra=pos[0], dec=pos[1], z=num(cells[3]), logM=logM, sig=sig, survives=bold)
    j = tex.index(r"\label{Table:BinaryProperties}")
    for line in tex[j:tex.index(r"\end{tabular}", j)].splitlines():
        if "&" not in line or "P-value" in line or "[days]" in line:
            continue
        cells = [c.strip() for c in line.rstrip("\\").split("&")]
        name = re.sub(r"\\textbf\{(.*?)\}", r"\1", cells[0]).strip()
        try:
            periods[name] = float(cells[2])
        except ValueError:  # header rows
            continue
    rows = []
    for name, p in props.items():
        if name not in periods:
            print("PTF: no period for", name, file=sys.stderr)
            continue
        rows.append(dict(p, P_obs=periods[name]))
    return rows


def read_liu() -> list[dict]:
    tex = (DATA / "periodic/src_1906.08315/ms_arxiv_v2.tex").read_text(encoding="utf-8", errors="replace")
    def block(label):
        i = tex.index(label)
        return tex[tex.index(r"\startdata", i):tex.index(r"\enddata", i)]
    periods = {}
    for line in block(r"\label{tab:P_xi}").splitlines():
        if "&" not in line or line.strip().startswith("%"):
            continue
        cells = [c.strip() for c in line.rstrip("\\").split("&")]
        periods[cells[0].replace("--", "-").strip()] = float(re.match(r"([\d.]+)", cells[1]).group(1))
    rows = []
    for line in block(r"\label{tab:spec}").splitlines():
        if "&" not in line or line.strip().startswith("%"):
            continue
        cells = [c.strip() for c in line.rstrip("\\").split("&")]
        name = cells[0].replace("--", "-").strip()
        key = name if name in periods else next((k for k in periods if k[:14] == name[:14]), None)
        if key is None:
            print("PS1: no period for", name, file=sys.stderr)
            continue
        m = re.match(r"PSO J([\d.]+)([+-][\d.]+)", key)
        rows.append(dict(cat="PS1", name=key, ra=float(m.group(1)), dec=float(m.group(2)), z=float(cells[6]),
                         logM=float(cells[5]), sig=SIG_DEFAULT, P_obs=periods[key]))
    return rows


def read_ztf() -> list[dict]:
    tex = (DATA / "periodic/src_2206.11497/ms.tex").read_text(encoding="utf-8", errors="replace")
    i = tex.index(r"\label{table1}")
    rows = []
    for line in tex[i:tex.index(r"\end{tabular}", i)].splitlines():
        if "&" not in line or "ID &" in line or "Magnitude" in line:
            continue
        cells = [c.strip() for c in line.rstrip("\\").split("&")]
        if len(cells) < 11 or ":" not in cells[1] or "---" in cells[10]:
            continue
        ra, dec = sex(cells[1], cells[2])
        mm = re.search(r"([\d.]+)\\pm([\d.]+)", cells[10])
        pm = re.search(r"([\d.]+)\\pm", cells[5])
        rows.append(dict(cat="ZTF", name=cells[0], ra=ra, dec=dec, z=float(cells[3]), logM=float(mm.group(1)),
                         sig=float(mm.group(2)), P_obs=float(pm.group(1))))
    return rows


def read_des() -> list[dict]:
    """Chen et al. 2020 (arXiv:2008.12329, Tables 'info' and 'parameters'): the five DES+SDSS candidates, transcribed;
    periods are the multi-band GLS periods (J0252's g-band GLS period from Liao et al. 2021, arXiv:2008.12317)."""
    rows = [("J024613.89-004028.2", 1.736, 8.74, 0.07, 1184), ("J024703.24-010032.0", 2.525, 8.60, 0.02, 1705),
            ("J024944.66-000036.8", 1.295, 9.00, 0.04, 1153), ("J025214.67-002813.7", 1.530, 8.40, 0.10, 1607),
            ("J025406.26+002753.7", 1.765, 8.92, 0.04, 1467)]
    out = []
    for name, z, logM, sig, P in rows:
        ra, dec = from_jname(name)
        out.append(dict(cat="DES", name=name, ra=ra, dec=dec, z=z, logM=logM, sig=max(sig, SIG_DEFAULT), P_obs=P))
    return out


# -------------------------------------------------------------------------------------------------------------- physics
class LimitMap:
    def __init__(self):
        d = np.load(DATA / "ng15/cw/15yr_cw_analysis-main/data/15yr_cw_3d_limits_v4.npz")
        self.F = d["F_edges"]
        self.UL = d["UL_skies"]
        self.sig = d["UL_sky_sigmas"]
        self.nside = hp.npix2nside(self.UL.shape[1])
        self.f_last = np.sqrt(self.F[-2] * self.F[-1])

    def __call__(self, f: float, ra: float, dec: float) -> tuple[float, float, bool]:
        pix = hp.ang2pix(self.nside, np.pi / 2 - np.radians(dec), np.radians(ra))
        if f < self.F[0]:
            raise ValueError("below the map")
        if f >= self.F[-1]:
            return self.UL[-1, pix] * f / self.f_last, self.sig[-1, pix] * f / self.f_last, True
        i = np.searchsorted(self.F, f, side="right") - 1
        return self.UL[i, pix], self.sig[i, pix], False


def h0_of(mcz_msun, f, dl_m):
    gm = G * mcz_msun * MSUN
    return 2 * gm ** (5 / 3) * (np.pi * f) ** (2 / 3) / (C ** 4 * dl_m)


def mc_of_h0(h0, f, dl_m):
    return (h0 * C ** 4 * dl_m / (2 * (np.pi * f) ** (2 / 3))) ** (3 / 5) / G / MSUN


def draws(logM, sig, z):
    m = 10 ** RNG.normal(logM, sig, NDRAW)
    q = 10 ** RNG.uniform(-1, 0, NDRAW)
    mc = m * q ** 0.6 / (1 + q) ** 1.2
    return (1 + z) * mc


def q1(logM, z):
    return (1 + z) * 10 ** logM * 0.25 ** 0.6  # equal masses: Mc = M / 2^(6/5) = 0.435 M


# --------------------------------------------------------------------------------------------------------------- driver
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    lim = LimitMap()
    ag = read_agarwal()
    print(f"Agarwal tables: {len(ag)} targets, {sum(r['logM'] is not None for r in ag)} with an EM mass")

    # convention check: Agarwal's Mc95 should follow from h95, f, z and D_L (flat LCDM; they quote D_L fixed per target)
    ratios = []
    for r in ag:
        dl = COSMO.luminosity_distance(r["z"]).to(u.m).value
        mc_from_h = mc_of_h0(r["h95"], r["f_tab"], dl) / (1 + r["z"])
        ratios.append(mc_from_h / r["Mc95"])
    ratios = np.array(ratios)
    zs = np.array([r["z"] for r in ag])
    print(f"Mc(h95, f, D_L)/Mc95 (should be ~1 up to their cosmology): median {np.median(ratios):.3f}, "
          f"16-84 % {np.percentile(ratios, 16):.3f}-{np.percentile(ratios, 84):.3f}; "
          f"times (1+z): median {np.median(ratios * (1 + zs)):.3f}, 16-84 % {np.percentile(ratios * (1 + zs), 16):.3f}-"
          f"{np.percentile(ratios * (1 + zs), 84):.3f} (so their Mc95 is the redshifted chirp mass if this is ~1)")
    # and the map at their positions vs their all-sky-derived improvement factor (median 2.24 in the paper)
    impr = []
    for r in ag:
        ul, _, ex = lim(r["f_tab"], r["ra"], r["dec"])
        impr.append(ul / r["h95"])
    impr = np.array(impr)
    print(f"map UL / targeted h95 at the 114 positions: median {np.median(impr):.2f}, range {impr.min():.2f}-{impr.max():.2f} "
          f"(paper: median 2.24, 0.82-14.94)")
    if args.check:
        return

    cats = {"CRTS": ag, "PTF": read_charisi(), "PS1": read_liu(), "DES": read_des(), "ZTF": read_ztf()}
    out = []
    for cat, rows in cats.items():
        for r in rows:
            if r["logM"] is None:
                continue
            dl = COSMO.luminosity_distance(r["z"]).to(u.m).value
            mcz = draws(r["logM"], r["sig"], r["z"])
            rec = dict(cat=cat, name=r["name"], ra=r["ra"], dec=r["dec"], z=r["z"], P_obs=r["P_obs"], logM=r["logM"],
                       sig=r["sig"], survives=r.get("survives", True))
            for tag, f in (("2P", 2 / (r["P_obs"] * 86400)), ("1P", 1 / (r["P_obs"] * 86400))):
                h = h0_of(mcz, f, dl)
                ul, ulsig, ex = lim(f, r["ra"], r["dec"])
                rec[f"f_{tag}"] = f
                rec[f"h0med_{tag}"] = np.median(h)
                rec[f"h0q1_{tag}"] = h0_of(q1(r["logM"], r["z"]), f, dl)
                rec[f"h0q1p1s_{tag}"] = h0_of(q1(r["logM"] + r["sig"], r["z"]), f, dl)
                rec[f"UL_{tag}"] = ul
                rec[f"ULextrap_{tag}"] = ex
                rec[f"Pexcl_{tag}"] = float(np.mean(h > ul))
                rec[f"y_{tag}"] = ul / np.median(h)
                rec[f"Mcz_lim_{tag}"] = mc_of_h0(ul, f, dl) / (1 + r["z"])  # chirp-mass limit from the all-sky map
            if cat == "CRTS":
                f = r["f_tab"]
                h = h0_of(mcz, f, dl)
                rec["h95_targeted"] = r["h95"]
                rec["Mc95_targeted"] = r["Mc95"]
                rec["Pexcl_targeted"] = float(np.mean(h > r["h95"]))
                rec["y_targeted"] = r["h95"] / np.median(h)
                rec["Mc95_over_McEM_q1"] = r["Mc95"] / (10 ** r["logM"] * 0.25 ** 0.6)
            out.append(rec)

    keys = list(out[0].keys())
    for rec in out:
        for k in rec:
            if k not in keys:
                keys.append(k)
    with open(DATA / "periodic/candidates_cw.csv", "w") as fh:
        fh.write(",".join(keys) + "\n")
        for rec in out:
            fh.write(",".join("" if rec.get(k) is None else (f"{rec[k]:.4g}" if isinstance(rec[k], float) else str(rec[k])) for k in keys) + "\n")

    # summary
    from scipy.stats import beta
    def cp(k, n):
        lo = 0.0 if k == 0 else beta.ppf(0.025, k, n - k + 1)
        hi = 1.0 if k == n else beta.ppf(0.975, k + 1, n - k)
        return lo, hi

    print()
    for cat in cats:
        sub = [r for r in out if r["cat"] == cat]
        if not sub:
            continue
        for tag in ("2P", "1P"):
            n = len(sub)
            k = sum(r[f"Pexcl_{tag}"] >= 0.95 for r in sub)
            expct = sum(r[f"Pexcl_{tag}"] for r in sub)
            ys = np.array([r[f"y_{tag}"] for r in sub])
            nex = sum(r[f"ULextrap_{tag}"] for r in sub)
            fmax = max(r[f"f_{tag}"] for r in sub) * 1e9
            lo, hi = cp(k, n)
            print(f"{cat:5s} f_gw={tag}: N={n:3d}  excluded(95%)={k} [{100*lo:.1f}-{100*hi:.1f} % of N]  expected={expct:.2f}  "
                  f"UL/h0_med: min {ys.min():.1f}, median {np.median(ys):.0f}, max {ys.max():.0f}; "
                  f"frac with P_excl>0.05: {np.mean([r[f'Pexcl_{tag}']>0.05 for r in sub]):.2f}; "
                  f"{nex} above the map ({fmax:.0f} nHz max)")
        if cat == "CRTS":
            k = sum(r["Pexcl_targeted"] >= 0.95 for r in sub)
            expct = sum(r["Pexcl_targeted"] for r in sub)
            ys = np.array([r["y_targeted"] for r in sub])
            rr = np.array([r["Mc95_over_McEM_q1"] for r in sub])
            lo, hi = cp(k, len(sub))
            print(f"CRTS  targeted (Agarwal h95): N={len(sub)}  excluded={k} [{100*lo:.1f}-{100*hi:.1f} %]  expected={expct:.2f}  "
                  f"h95/h0_med: min {ys.min():.1f}, median {np.median(ys):.0f}; Mc95/Mc_EM(q=1): min {rr.min():.1f}, "
                  f"median {np.median(rr):.0f}; n with Mc95 < Mc_EM(q=1): {np.sum(rr < 1)}; "
                  f"n with P_excl>0.05: {np.sum([r['Pexcl_targeted']>0.05 for r in sub])}")
        if cat == "PTF":
            surv = [r for r in sub if r["survives"]]
            k = sum(r["Pexcl_2P"] >= 0.95 for r in surv)
            print(f"PTF   the 33 that survive the extended light curves: N={len(surv)} excluded={k} expected={sum(r['Pexcl_2P'] for r in surv):.2f} "
                  f"UL/h0 min {min(r['y_2P'] for r in surv):.1f}")

    print("\nclosest calls (smallest UL / median h0, f_gw = 2/P):")
    for r in sorted(out, key=lambda r: r["y_2P"])[:12]:
        print(f"  {r['cat']:4s} {r['name']:28s} z={r['z']:.2f} P={r['P_obs']:.0f} d logM={r['logM']:.2f}+-{r['sig']:.2f} "
              f"f={r['f_2P']*1e9:.1f} nHz h0_med={r['h0med_2P']:.1e} (q=1,+1s: {r['h0q1p1s_2P']:.1e}) UL={r['UL_2P']:.1e}"
              f"{'*' if r['ULextrap_2P'] else ''} y={r['y_2P']:.1f} P_excl={r['Pexcl_2P']:.3f}"
              + (f" | targeted h95={r['h95_targeted']:.1e} y={r['y_targeted']:.1f} P_excl={r['Pexcl_targeted']:.3f}" if r["cat"] == "CRTS" else ""))


if __name__ == "__main__":
    main()
