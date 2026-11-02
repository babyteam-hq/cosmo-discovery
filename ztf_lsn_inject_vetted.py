"""Injection-recovery through the full vetting (R-141, cosmo-deep-10b).

The injection set is R-135's, unchanged and drawn identically (ztf_lsn_inject.draw with the same seed per lens, the same
explosion times, the same kappa12 response at 1.2"), and every injected light curve now goes through
ztf_lsn_vetting.vet, the function the real light curves go through. What the vetting needs beyond R-135's r-band
injection is made here, from a separate random stream so the injection set is untouched:
  - the image position angles (uniform; images 1.2" from the centre, as R-130 fixed), for the centroid stage;
  - a g-band light curve at the lens's real clean g epochs (IBE search rows: INFOBITS == 0, SEEING <= 4", MAGLIM >= 19)
    in the vetting span: the injected SN's ZTF-g flux times the response (the target's kappa12-on-seeing fit evaluated at
    the g seeing), plus Gaussian noise with the per-epoch error the r band shows at that limiting magnitude
    (e = C x F_5sigma(maglim) / 5 x the lens's r error scale, C the lens's median r ratio). g cutouts are not fetched
    for the 168 lenses (IRSA throttled, 2026-11-02); `--check-g` compares this error model with real g photometry
    where it exists.
Writes data/ztf_lsn/inject_vetted/<id>.csv (one row per injection, with the stage flags) and
data/ztf_lsn/inject_vetted_summary.csv. `--real` instead vets every lens's own light curve and writes
data/ztf_lsn/vet_real.csv.

Usage (from research/cosmo-discovery):
    python ztf_lsn_inject_vetted.py [--thr 5.31] [--only ID] [--procs 4] [--real]
"""
from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from ztf_lsn_detect import WIN, clean, prepare, windows
from ztf_lsn_inject import draw, light_curve, model_for
from ztf_lsn_vetting import bandflux as vbandflux
from ztf_lsn_vetting import vet

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "ztf_lsn"
OUT = DATA / "inject_vetted"
MU_WIDE = False
G_ERR_FACTOR = 1.15  # g's e/F_5sigma over r's at the lens centre, measured on the two lenses with real g (1.18, 1.12)
OFFSET = 1.2  # arcsec; --offset 2.0 is the systematic variant (written to inject_vetted_off2.0)


def kappa_fit(tid: str, d: pd.DataFrame, col: str) -> tuple[np.ndarray, np.ndarray] | None:
    """Per-epoch kappa (col = kappa0 or kappa12) for the r epochs of d (fit on seeing per field/ccd/quadrant with >= 8
    stamps, else the target's fit, as ztf_lsn_inject.psf_response), and the target's global fit coefficients."""
    p = DATA / tid / "psfcal.csv"
    if not p.exists():
        return None
    c = pd.read_csv(p).merge(d[["stem", "field", "ccdid", "qid"]], on="stem", how="inner")
    if len(c) < 8:
        return None
    gfit = np.polyfit(c.seeing, c[col], 1)
    out = np.polyval(gfit, d.seeing.values)
    key = (d.field.astype(str) + "_" + d.ccdid.astype(str) + "_" + d.qid.astype(str)).values
    ckey = (c.field.astype(str) + "_" + c.ccdid.astype(str) + "_" + c.qid.astype(str)).values
    for k in np.unique(key):
        m = ckey == k
        if m.sum() >= 8:
            fit = np.polyfit(c.seeing.values[m], c[col].values[m], 1)
            out[key == k] = np.polyval(fit, d.seeing.values[key == k])
    return np.clip(out, 0.05, 1.0), gfit


def g_epochs(tid: str) -> pd.DataFrame:
    e = pd.read_csv(DATA / tid / "epochs.csv")
    e = e[(e.filtercode == "zg") & (e.imgtypecode == "o") & (e.infobits == 0) & (e.seeing <= 4.0) & (e.maglimit >= 19.0)]
    return e.sort_values("obsjd").reset_index(drop=True)


def err_ratio(d: pd.DataFrame) -> float:
    """Median of the raw r error over F_5sigma(maglim)/5 at the lens centre."""
    f5 = 10 ** (-0.4 * (d.maglim.values - 23.9)) / 5
    return float(np.median(d.e.values / f5))


def ctrl_series(d: pd.DataFrame) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    return [prepare(d, p) for p in ("_N", "_S", "_E", "_W")]


def ctrl_S(ctrls, t_w: float) -> float:
    best = -np.inf
    for rc, ec, jdc, _ in ctrls:
        m = (jdc >= t_w) & (jdc < t_w + WIN)
        if m.sum() >= 2:
            wt = 1 / ec[m] ** 2
            best = max(best, float((wt * rc[m]).sum() / np.sqrt(wt.sum())))
    return best


def load(lens: pd.Series):
    p = DATA / lens.id / "phot.csv"
    if not p.exists():
        return None
    d0 = pd.read_csv(p)
    d = clean(d0)
    if len(d) < 20:
        return None
    k12 = kappa_fit(lens.id, d, "kappa12")
    k0 = kappa_fit(lens.id, d, "kappa0")
    if k12 is None or k0 is None:
        return None
    return d0, d, k12, k0


def run_lens(args) -> pd.DataFrame | None:
    lens, n, thr, seed = args
    L = load(lens)
    if L is None:
        return None
    d0, d, (resp, gfit12), (k0, gfit0) = L
    r, e, jd, scale = prepare(d, "")
    ratio = 1.0
    if OFFSET != 1.2:  # the measured response at 2.0" (kappa20 from the same PSF stamps)
        k20 = kappa_fit(lens.id, d, f"kappa{int(round(OFFSET * 10))}")
        if k20 is None:
            return None
        ratio = float(np.median(k20[0] / resp))
        resp = k20[0]
    ctrls = ctrl_series(d)
    ge = g_epochs(lens.id)
    jdg_all, seeg = ge.obsjd.values, ge.seeing.values
    resp_g = np.clip(np.polyval(gfit12, seeg), 0.05, 1.0) * ratio
    k0g = np.clip(np.polyval(gfit0, seeg), 0.05, 1.0)
    eg_all = G_ERR_FACTOR * err_ratio(d) * 10 ** (-0.4 * (ge.maglimit.values - 23.9)) / 5 * scale
    fwhm = float(pd.read_csv(DATA / lens.id / "psfcal.csv").fwhm_psf.median())
    # R-135's injection set, drawn identically
    rng = np.random.default_rng(seed + int(lens["rank"]))
    draws = draw(rng, n, lens)
    w0_all = windows(r, e, jd)
    t0s = rng.uniform(jd.min() - 100, jd.max() + 100, n)
    rng2 = np.random.default_rng(141_000 + int(lens["rank"]))  # the vetting's own randomness
    if MU_WIDE:  # robustness variant: each injection's total magnification times 10^U(-0.5, 1) (0.3-10 x the model)
        rng3 = np.random.default_rng(142_000 + int(lens["rank"]))
        fac = 10 ** rng3.uniform(-0.5, 1.0, len(draws))
        draws["mus"] = [m * f for m, f in zip(draws.mus, fac)]
        draws["mu_tot"] = draws.mu_tot * fac
    out = []
    for (i, row), t0 in zip(draws.iterrows(), t0s):
        lc = light_curve(row, jd, t0)
        inj = lc * resp
        peak_model = lc.max()
        rr = r + inj
        w = windows(rr, e, jd)
        lo_t, hi_t = t0 - 60.0, t0 + float(np.max(row.delays)) + 200.0
        w = w[(w.t0 >= lo_t) & (w.t0 <= hi_t)].reset_index(drop=True)
        w0 = w0_all[(w0_all.t0 >= lo_t) & (w0_all.t0 <= hi_t)].reset_index(drop=True)
        pa = rng2.uniform(0, 2 * np.pi, row.nimg)
        noise_g = rng2.normal(0, 1, jdg_all.size)
        noise_c = rng2.normal(0, 1, 2)
        # g light curve at the real g epochs (only those near the SN matter to the vetting)
        near = (jdg_all > t0 - 200) & (jdg_all < t0 + 600)
        lcg = np.zeros(jdg_all.size)
        if near.any():
            lcg[near] = light_curve_band(row, jdg_all[near], t0, "ztfg")
        gg = lcg * resp_g + noise_g * eg_all

        def centroid(t_w, S, row=row, t0=t0, pa=pa, rad=OFFSET):
            tt = np.arange(t_w, t_w + WIN, 1.0)
            fl = np.array([image_flux(row, tt, t0, k) for k in range(row.nimg)])
            if fl.sum() <= 0:
                return 0.0
            x = (fl * rad * np.cos(pa)).sum() / fl.sum()
            y = (fl * rad * np.sin(pa)).sum() / fl.sum()
            sig = fwhm / 2.3548 / max(S, 1.0)
            return float(np.hypot(x + sig * noise_c[0], y + sig * noise_c[1]))

        v = vet(w, thr, lambda t_w: ctrl_S(ctrls, t_w), centroid, rr, e, jd, k0, gg, eg_all, jdg_all, k0g, lens,
                w0_restricted=w0)
        out.append({
            "id": lens.id, "type": row.type, "template": row.template, "MB": row.MB, "z_s": row.z_s, "nimg": row.nimg,
            "mu_tot": row.mu_tot, "t0": t0,
            "m_peak": 23.9 - 2.5 * np.log10(peak_model) if peak_model > 0 else np.nan,
            "m_peak_meas": 23.9 - 2.5 * np.log10(inj.max()) if inj.max() > 0 else np.nan,
            "n_g_span": int(((jdg_all >= v["t_w"] - 60) & (jdg_all <= v["t_w"] + 120)).sum()) if np.isfinite(v["t_w"]) else 0,
            **{k: v[k] for k in ("T", "T2", "t_w", "s0", "s1", "s2", "s3", "s4", "hyp", "A", "var_amp")},
            "detected": bool(v["s0"]), "vetted": bool(v["s4"]),
        })
    df = pd.DataFrame(out)
    df["thr"] = thr
    df.to_csv(OUT / f"{lens.id}.csv", index=False, float_format="%.6g")  # noqa: OUT is set per process below
    return df


def light_curve_band(row, jd: np.ndarray, t0: float, band: str) -> np.ndarray:
    m = model_for(row.template)
    if row.template == "salt3":
        m.set(z=row.z_s, x1=row.x1, c=row.c)
    else:
        m.set(z=row.z_s)
    m.set_source_peakabsmag(row.MB, "bessellb", "ab", cosmo=_cosmo())
    tot = np.zeros_like(jd)
    for mu, dt in zip(row.mus, row.delays):
        m.set(t0=t0 + dt)
        sel = (jd >= m.mintime()) & (jd <= m.maxtime())
        if sel.any():
            tot[sel] += mu * vbandflux(m, band, jd[sel])
    return tot


def image_flux(row, tt: np.ndarray, t0: float, k: int) -> float:
    m = model_for(row.template)
    if row.template == "salt3":
        m.set(z=row.z_s, x1=row.x1, c=row.c)
    else:
        m.set(z=row.z_s)
    m.set_source_peakabsmag(row.MB, "bessellb", "ab", cosmo=_cosmo())
    m.set(t0=t0 + row.delays[k])
    sel = (tt >= m.mintime()) & (tt <= m.maxtime())
    return float(row.mus[k] * m.bandflux("ztfr", tt[sel], zp=23.9, zpsys="ab").sum()) if sel.any() else 0.0


def _cosmo():
    from ztf_lsn_inject import COSMO
    return COSMO


def summarize(df: pd.DataFrame, lens: pd.Series) -> dict:
    s = {"id": lens.id, "class": lens["class"], "n_inj": len(df), "eff": df.detected.mean(), "eff_vet": df.vetted.mean()}
    for st in ("s1", "s2", "s3", "s4"):
        s[f"eff_{st}"] = df[st].mean()
    for t in ("Ia", "CC"):
        sel = df.type == "Ia" if t == "Ia" else df.type != "Ia"
        s[f"n_{t}"] = int(sel.sum())
        s[f"eff_{t}"] = df.detected[sel].mean() if sel.any() else np.nan
        s[f"effv_{t}"] = df.vetted[sel].mean() if sel.any() else np.nan
    for lo in np.arange(18.5, 23.0, 0.5):
        sel = (df.m_peak >= lo) & (df.m_peak < lo + 0.5)
        s[f"n_m{lo:.1f}"] = int(sel.sum())
        s[f"eff_m{lo:.1f}"] = df.detected[sel].mean() if sel.sum() >= 3 else np.nan
        s[f"effv_m{lo:.1f}"] = df.vetted[sel].mean() if sel.sum() >= 3 else np.nan
    b = df.m_peak < 21.0
    s["n_brighter_21"] = int(b.sum())
    s["eff_brighter_21"] = df.detected[b].mean() if b.sum() >= 3 else np.nan
    s["effv_brighter_21"] = df.vetted[b].mean() if b.sum() >= 3 else np.nan
    s["k_brighter_21"] = int(df.vetted[b].sum())
    return s


def real_one(args) -> dict | None:
    """The vetting on a lens's own light curve (and g, if its g cutouts were fetched and measured)."""
    lens, thr = args
    L = load(lens)
    if L is None:
        return None
    d0, d, (resp, gfit12), (k0, gfit0) = L
    r, e, jd, scale = prepare(d, "")
    ctrls = ctrl_series(d)
    dg = d0[d0["filter"].astype(str).str.contains("_g")]
    dg = dg[(dg.infobits == 0) & (dg.seeing <= 4.0) & (dg.maglim >= 19.0) & np.isfinite(dg.f) & (dg.e > 0)]
    dg = dg.sort_values("jd").reset_index(drop=True)
    if len(dg) >= 5:
        g, eg, jdg, _ = prepare(dg, "")
        k0g = np.clip(np.polyval(gfit0, dg.seeing.values[np.isfinite(dg.f.values)]), 0.05, 1.0)
    else:
        g = eg = jdg = k0g = None
    w = windows(r, e, jd)

    def centroid(t_w, S):
        from ztf_lsn_vet import stack
        sel = (d.jd.values >= t_w) & (d.jd.values < t_w + WIN)
        _, pos, cen = stack(lens.id, d, lens.ra, lens.dec, sel)
        return float(np.hypot(cen[0] - pos[0], cen[1] - pos[1]) * 1.0126) if np.isfinite(cen[0]) else 0.0

    v = vet(w, thr, lambda t_w: ctrl_S(ctrls, t_w), centroid, r, e, jd, k0, g, eg, jdg, k0g, lens)
    # the variability statistic for every lens (at its best window, or mid-span when nothing passes stage 0)
    from ztf_lsn_vetting import season_dev
    tw = v["t_w"] if np.isfinite(v["t_w"]) else float(np.median(jd))
    amp, sig = season_dev(r, e, jd, tw)
    ctrl_amps = [season_dev(rc, ec, jdc, tw)[0] for rc, ec, jdc, _ in ctrls]
    return {"id": lens.id, "class": lens["class"], **v, "var_amp_any": amp, "var_sig": sig,
            "ctrl_var_amps": ";".join(f"{a:.2f}" for a in ctrl_amps), "n_g": 0 if g is None else int(len(g)),
            "err_ratio": err_ratio(d), "scale": scale}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--thr", type=float, default=5.31, help="R-135's candidate threshold (704 controls, 99.5th pct)")
    ap.add_argument("--only", default="")
    ap.add_argument("--seed", type=int, default=130)
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--redo", action="store_true")
    ap.add_argument("--real", action="store_true")
    ap.add_argument("--stop", type=int, default=10**6, help="ranks below this only")
    ap.add_argument("--offset", type=float, default=1.2, help="image offset (arcsec): 1.2 (fixed) or the 2.0 variant")
    ap.add_argument("--mu-wide", action="store_true", help="robustness: total mu x 10^U(-0.5, 1)")
    ap.add_argument("--colour", default="peak", help="peak (row 8) or chi2 (row 7's first-frozen form, robustness)")
    a = ap.parse_args()
    global OFFSET, OUT, MU_WIDE
    MU_WIDE = a.mu_wide
    import ztf_lsn_vetting
    ztf_lsn_vetting.COLOUR_MODE = a.colour
    OFFSET = a.offset
    if a.offset != 1.2:
        OUT = DATA / f"inject_vetted_off{a.offset:.1f}"
    if a.colour != "peak":
        OUT = Path(str(OUT) + f"_{a.colour}")
    if a.mu_wide:
        OUT = Path(str(OUT) + "_muwide")
    t = pd.read_csv(HERE / "lsn_targets.csv")
    t = t[t["rank"] < a.stop]
    if a.only:
        t = t[t.id.isin(a.only.split(","))]
    lenses = [lens for _, lens in t.iterrows()]
    if a.real:
        with ProcessPoolExecutor(a.procs) as ex:
            rows = [r for r in ex.map(real_one, [(l, a.thr) for l in lenses]) if r is not None]
        out = pd.DataFrame(rows)
        if not a.only:
            out.to_csv(DATA / f"vet_real{'' if a.colour == 'peak' else '_' + a.colour}.csv", index=False, float_format="%.5g")
        pd.set_option("display.width", 250)
        print(out[out.s0 | (out["T"] > a.thr)].to_string())
        return
    OUT.mkdir(exist_ok=True)
    todo = [l for l in lenses if a.redo or not (OUT / f"{l.id}.csv").exists()]
    with ProcessPoolExecutor(a.procs) as ex:
        for lens, df in zip(todo, ex.map(run_lens, [(l, a.n, a.thr, a.seed) for l in todo])):
            if df is not None:
                s = summarize(df, lens)
                print(f"{lens.id} ({lens['class']}): eff {s['eff']:.3f} -> vetted {s['eff_vet']:.3f} "
                      f"(s1 {s['eff_s1']:.3f} s2 {s['eff_s2']:.3f} s3 {s['eff_s3']:.3f}); m<21 {s['eff_brighter_21']:.2f} -> "
                      f"{s['effv_brighter_21']:.2f}", flush=True)
    rows = []
    for lens in lenses:
        p = OUT / f"{lens.id}.csv"
        if p.exists():
            rows.append(summarize(pd.read_csv(p), lens))
    if not a.only and rows:
        pd.DataFrame(rows).to_csv(DATA / f"{OUT.name}_summary.csv", index=False, float_format="%.5g")


if __name__ == "__main__":
    sys.exit(main())
