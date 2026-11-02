"""Injection-recovery of lensed supernovae into the real forced light curves (R-130, cosmo-deep-10), as fixed in
results-R130.md.

Per lens, N (default 200) lensed SNe drawn as the forecasts draw them (type, template, peak M_B, source redshift,
2 or 4 images with lognormal magnifications, relative delays 0-30 d, image positions 1.2" from the centre, explosion time
uniform over the data span +- 100 d). Each image's ZTF-r light curve (sncosmo) times the epoch's measured estimator
response at 1.2" is added to the lens's real baseline-subtracted light curve, and the statistic of ztf_lsn_detect.py is
recomputed. Detected = T2 > threshold (a window above the candidate threshold with a neighbouring window at S > 3).

Writes data/ztf_lsn/inject/<id>.csv (one row per injection) and data/ztf_lsn/inject_summary.csv (one row per lens).

Usage (from research/cosmo-discovery):
    python ztf_lsn_inject.py [--n 200] [--thr 5.0] [--only ID] [--seed 130]
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import sncosmo
from astropy.cosmology import FlatLambdaCDM

from ztf_lsn_detect import clean, prepare, stat, windows

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "ztf_lsn"
warnings.filterwarnings("ignore")
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
# Sheu et al. 2023 Table (from Shu et al.): fraction of CC by type, peak M_B and its population scatter, template
CC_TYPES = [
    ("IIP", 0.5583, -16.75, 0.98, "nugent-sn2p"),
    ("Ic", 0.1700, -17.66, 1.18, "nugent-sn1bc"),
    ("IIb", 0.1243, -16.99, 0.92, "v19-2006t-corr"),
    ("Ib", 0.0900, -17.45, 1.12, "nugent-sn1bc"),
    ("IIL", 0.0334, -17.98, 0.86, "nugent-sn2l"),
    ("IIn", 0.0240, -18.53, 1.36, "nugent-sn2n"),
]
IA_MB, IA_SIG = -19.25, 0.50
IMAGE_OFFSET = "resp_1.2"  # the Gaussian-estimator response; superseded by the measured kappa12 (psf_response) where psfcal.csv exists
CLASS_S_IA_FRACTION = 0.106  # sum N_Ia / sum (N_Ia + N_cc) over Shu's 110 lenses with rates
_SOURCES: dict[str, sncosmo.Source] = {}


def model_for(template: str) -> sncosmo.Model:
    if template not in _SOURCES:
        _SOURCES[template] = sncosmo.get_source(template)
    return sncosmo.Model(source=_SOURCES[template])


MU_TOTAL_FIXED: float | None = None  # robustness variant (--mu-total): every injection's images rescaled to this total


def draw(rng: np.random.Generator, n: int, lens: pd.Series) -> pd.DataFrame:
    rows = []
    cls = lens["class"]
    if cls == "S" and np.isfinite(lens.n_ia) and np.isfinite(lens.n_cc) and (lens.n_ia + lens.n_cc) > 0:
        p_ia = lens.n_ia / (lens.n_ia + lens.n_cc)
    else:
        p_ia = CLASS_S_IA_FRACTION
    cc_p = np.array([c[1] for c in CC_TYPES])
    cc_p /= cc_p.sum()
    for _ in range(n):
        if rng.random() < p_ia:
            typ, mb, template = "Ia", rng.normal(IA_MB, IA_SIG), "salt3"
            x1, c = rng.normal(0, 1), rng.normal(0, 0.1)
        else:
            k = rng.choice(len(CC_TYPES), p=cc_p)
            typ, _, m0, s0, template = CC_TYPES[k]
            mb, x1, c = rng.normal(m0, s0), np.nan, np.nan
        if cls == "S" and np.isfinite(lens.z_s):
            zs = float(lens.z_s)
        else:
            zl = float(lens.z_l) if np.isfinite(lens.z_l) else 0.45
            zs = zl * max(1.0, rng.normal(2.0, 0.5))
        nimg = 2 if rng.random() < 0.7 else 4
        mus = np.exp(rng.normal(1.5, 0.35, nimg))
        if MU_TOTAL_FIXED is not None:
            mus *= MU_TOTAL_FIXED / mus.sum()
        elif cls == "S" and np.isfinite(lens.mu) and lens.mu > 0:
            mus *= lens.mu / mus.sum()
        delays = np.concatenate([[0.0], rng.uniform(0, 30, nimg - 1)])
        rows.append({"type": typ, "template": template, "MB": mb, "x1": x1, "c": c, "z_s": zs, "nimg": nimg,
                     "mu_tot": mus.sum(), "mus": mus, "delays": delays})
    return pd.DataFrame(rows)


def light_curve(row, jd: np.ndarray, t0: float) -> np.ndarray:
    """Summed (over images) ZTF-r flux in microJansky at the epochs, before the estimator response."""
    m = model_for(row.template)
    if row.template == "salt3":
        m.set(z=row.z_s, x1=row.x1, c=row.c)
    else:
        m.set(z=row.z_s)
    m.set_source_peakabsmag(row.MB, "bessellb", "ab", cosmo=COSMO)
    tot = np.zeros_like(jd)
    for mu, dt in zip(row.mus, row.delays):
        m.set(t0=t0 + dt)
        tmin, tmax = m.mintime(), m.maxtime()
        sel = (jd >= tmin) & (jd <= tmax)
        if sel.any():
            tot[sel] += mu * m.bandflux("ztfr", jd[sel], zp=23.9, zpsys="ab")
    return tot


def psf_response(tid: str, d: pd.DataFrame) -> np.ndarray | None:
    """kappa12 per epoch: the measured response of the estimator to a unit source 1.2" off-centre with the true
    difference-image PSF (ztf_lsn_psfcal.py, every 10th epoch), carried to every epoch by a linear fit on SEEING per
    (field, ccd, quadrant) when that set has >= 8 stamps, else the target's fit, clipped to [0.05, 1]."""
    p = DATA / tid / "psfcal.csv"
    if not p.exists():
        return None
    c = pd.read_csv(p)
    if len(c) < 8:
        return None
    c = c.merge(d[["stem", "field", "ccdid", "qid"]], on="stem", how="inner")
    if len(c) < 8:
        return None
    glob_fit = np.polyfit(c.seeing, c.kappa12, 1)
    out = np.polyval(glob_fit, d.seeing.values)
    key = d.field.astype(str) + "_" + d.ccdid.astype(str) + "_" + d.qid.astype(str)
    ckey = c.field.astype(str) + "_" + c.ccdid.astype(str) + "_" + c.qid.astype(str)
    for k in key.unique():
        m = (ckey == k).values
        if m.sum() >= 8:
            fit = np.polyfit(c.seeing.values[m], c.kappa12.values[m], 1)
            out[(key == k).values] = np.polyval(fit, d.seeing.values[(key == k).values])
    return np.clip(out, 0.05, 1.0)


def run_lens(lens: pd.Series, n: int, thr: float, rng: np.random.Generator) -> pd.DataFrame | None:
    p = DATA / lens.id / "phot.csv"
    if not p.exists():
        return None
    d = clean(pd.read_csv(p))
    if len(d) < 20:
        return None
    r, e, jd, scale = prepare(d, "")
    resp = psf_response(lens.id, d)
    if resp is None:
        return None  # no PSF calibration yet: the Gaussian response would overstate the efficiency
    assert resp.size == r.size
    draws = draw(rng, n, lens)
    w0_all = windows(r, e, jd)
    t0s = rng.uniform(jd.min() - 100, jd.max() + 100, n)
    out = []
    for (i, row), t0 in zip(draws.iterrows(), t0s):
        lc = light_curve(row, jd, t0)
        inj = lc * resp
        peak_model = lc.max()
        w = windows(r + inj, e, jd)
        # R-135 amendment: a recovery counts only at the injected SN's time (windows starting from 60 d before the first
        # image's model t0 (the template's peak or zero phase) to 200 observed days after the last image's), so a lens whose own light curve already holds
        # an above-threshold event (or a noise excursion) elsewhere does not count every injection as recovered
        lo_t, hi_t = t0 - 60.0, t0 + float(np.max(row.delays)) + 200.0
        w = w[(w.t0 >= lo_t) & (w.t0 <= hi_t)].reset_index(drop=True)
        s = stat(w)
        # and where the lens's own light curve already passes the cut in that span, the injection counts as missed
        # (conservative: no efficiency is credited during a real event)
        w0 = w0_all[(w0_all.t0 >= lo_t) & (w0_all.t0 <= hi_t)].reset_index(drop=True)
        s0 = stat(w0)
        if np.isfinite(s0["T2"]) and s0["T2"] > thr:
            s["T2"] = np.nan
        out.append({
            "id": lens.id, "type": row.type, "template": row.template, "MB": row.MB, "z_s": row.z_s, "nimg": row.nimg,
            "mu_tot": row.mu_tot, "t0": t0,
            "m_peak": 23.9 - 2.5 * np.log10(peak_model) if peak_model > 0 else np.nan,
            "m_peak_meas": 23.9 - 2.5 * np.log10(inj.max()) if inj.max() > 0 else np.nan,
            "n_ep_30d": int(((jd > t0 - 10) & (jd < t0 + 40 * (1 + row.z_s))).sum()),
            "T": s["T"], "T2": s["T2"], "detected": bool(np.isfinite(s["T2"]) and s["T2"] > thr),
            "flagged": bool(s["T"] > thr),
        })
    return pd.DataFrame(out)


def summarize(df: pd.DataFrame, lens: pd.Series, thr: float, scale: float | None = None) -> dict:
    s = {"id": lens.id, "class": lens["class"], "n_inj": len(df), "eff": df.detected.mean(),
         "eff_flag": df.flagged.mean(), "thr": thr}
    for t in ("Ia", "CC"):
        sel = df.type == "Ia" if t == "Ia" else df.type != "Ia"
        s[f"n_{t}"] = int(sel.sum())
        s[f"eff_{t}"] = df.detected[sel].mean() if sel.any() else np.nan
    for lo in np.arange(18.5, 23.0, 0.5):
        sel = (df.m_peak >= lo) & (df.m_peak < lo + 0.5)
        s[f"eff_m{lo:.1f}"] = df.detected[sel].mean() if sel.sum() >= 3 else np.nan
        s[f"n_m{lo:.1f}"] = int(sel.sum())
    s["eff_brighter_21"] = df.detected[df.m_peak < 21.0].mean() if (df.m_peak < 21.0).sum() >= 3 else np.nan
    return s


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--thr", type=float, default=5.0, help="candidate threshold (from the controls, ztf_lsn_detect.py)")
    ap.add_argument("--only", default="")
    ap.add_argument("--seed", type=int, default=130)
    ap.add_argument("--redo", action="store_true")
    ap.add_argument("--mu-total", type=float, default=None, help="robustness variant: fixed total magnification")
    a = ap.parse_args()
    global MU_TOTAL_FIXED
    MU_TOTAL_FIXED = a.mu_total
    tag = "" if a.mu_total is None else f"_mu{a.mu_total:g}"
    t = pd.read_csv(HERE / "lsn_targets.csv")
    if a.only:
        t = t[t.id.isin(a.only.split(","))]
    (DATA / f"inject{tag}").mkdir(exist_ok=True)
    rows = []
    for _, lens in t.iterrows():
        out = DATA / f"inject{tag}" / f"{lens.id}.csv"
        if out.exists() and not a.redo:
            df = pd.read_csv(out)
            if len(df) and abs(df.get("thr", pd.Series([a.thr])).iloc[0] - a.thr) < 1e-9:
                rows.append(summarize(df, lens, a.thr))
                continue
        rng = np.random.default_rng(a.seed + int(lens["rank"]))
        df = run_lens(lens, a.n, a.thr, rng)
        if df is None:
            continue
        df["thr"] = a.thr
        df.to_csv(out, index=False, float_format="%.6g")
        s = summarize(df, lens, a.thr)
        rows.append(s)
        print(f"{lens.id} ({lens['class']}): eff {s['eff']:.3f} (Ia {s['eff_Ia']:.3f}, CC {s['eff_CC']:.3f}); "
              f"m<21: {s['eff_brighter_21']:.2f}; median m_peak {df.m_peak.median():.2f}", flush=True)
    if not a.only:
        pd.DataFrame(rows).to_csv(DATA / f"inject_summary{tag}.csv", index=False, float_format="%.5g")
    elif rows:
        df = pd.read_csv(DATA / f"inject{tag}" / f"{t.iloc[0].id}.csv")
        print(df.groupby(pd.cut(df.m_peak, np.arange(17, 25, 0.5))).agg(n=("detected", "size"), eff=("detected", "mean"),
                                                                           T=("T", "median")).dropna().to_string())
        print(df.groupby("type").agg(n=("detected", "size"), eff=("detected", "mean"), m=("m_peak", "median")).to_string())


if __name__ == "__main__":
    sys.exit(main())
