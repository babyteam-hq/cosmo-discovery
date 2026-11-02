"""The candidate vetting as code (R-141, cosmo-deep-10b): the cuts R-135 applied by hand to the real candidates
(results-R130.md rows 4-5), written as one function that is applied, unchanged, to every real light curve and to every
injected one.

A light curve passes when every stage passes, in order:
  0. detection: T2 > thr (a 15-d window above the candidate threshold with a neighbouring window at S > 3, the shape cut);
     the detection window t_w is the best window that satisfies it. (For injections the windows are restricted to the
     injected SN's time span and no credit is given where the lens's own light curve passes there: R-135's amendments.)
  1. controls quiet: no control position (10" N, S, E, W) has S > 3 in the detection window ("a bad subtraction shows
     everywhere").
  2. centroid: the window's flux centroid within 2.5" of the lens. Real data: the stacked window image (ztf_lsn_vet.stack).
     Injections: the flux-weighted mean of the image positions at the window, plus the centroid's noise
     (sigma = PSF sigma / S per axis), since the injections are made in flux, not in pixels.
  3. no persistent variability: the r season means (365.25-d seasons from JD 2458000, >= 5 epochs, outside
     [t_w - 100, t_w + 250]) must not have any season deviating from the median of the season means by more than
     VAR_AMP microJansky at more than VAR_SIG sigma ("a lensed quasar or AGN shows in every season").
  4. colour and magnification: for some SN hypothesis at the source redshift (seven types: SALT3 Ia with x1 = 0, c = 0 at
     M_B = -19.25 and the six core-collapse types at their mean M_B, as the injections draw them; class S at the
     measured z_s, classes L and X on a z_s grid z_l x 1.0-3.0), fitting the hypothesis's r light curve (kappa0-corrected,
     as the candidates' magnitudes were) to the r epochs in [t_w - 60, t_w + 120] over a grid of explosion times gives a
     required total magnification A within MU_LO-MU_HI x the lens model's (class S: Shu's total mu; L, X: the median
     total of the injection draw), and the g epochs in the same span (if any; none = untested, passes) are consistent
     with A times the hypothesis's g light curve: the weighted mean g flux over the r peak windows [t_w, t_w + 30)
     within 3 sigma of the prediction, A's error included (row 8 amendment; the first-frozen form, chi2 p > 0.001 over
     the whole span, is COLOUR_MODE = "chi2", kept as a robustness row).
  5. (not code) no ZTF alert object within 6" at the time: injections have none by construction; for real candidates
     it was checked with ALeRCE (results-R130.md row 4).

The numbers (VAR_AMP, VAR_SIG, the 3-sigma control cut, the 2.5" centroid, MU_LO-MU_HI) were frozen in results-R130.md
(R-141 section) before any injection was vetted; VAR_AMP and VAR_SIG are set so that the two candidates R-135 rejected
as persistent variables fail and the two lens-galaxy transients pass, and their rejection rate on the control positions
is reported.
"""
from __future__ import annotations

import warnings
from functools import lru_cache

import numpy as np
import pandas as pd
import sncosmo
from astropy.cosmology import FlatLambdaCDM
from scipy import stats

from ztf_lsn_detect import WIN, stat

warnings.filterwarnings("ignore")
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
CTRL_S_MAX = 3.0
CENTROID_MAX = 2.5  # arcsec
VAR_AMP = 1.0  # microJansky
VAR_SIG = 5.0
MU_LO, MU_HI = 0.3, 3.0
COLOUR_P = 1e-3  # the chi2 form (COLOUR_MODE = "chi2", row 7 as first frozen; a robustness row)
COLOUR_Z = 3.0  # the peak form (row 8 amendment): |g_obs - g_pred| / sigma < 3 over the r peak windows [t_w, t_w + 30)
COLOUR_MODE = "peak"
HYPOTHESES = [  # (name, template, mean peak M_B), as the injections draw them (ztf_lsn_inject.py)
    ("Ia", "salt3", -19.25), ("IIP", "nugent-sn2p", -16.75), ("Ic", "nugent-sn1bc", -17.66),
    ("IIb", "v19-2006t-corr", -16.99), ("Ib", "nugent-sn1bc", -17.45), ("IIL", "nugent-sn2l", -17.98),
    ("IIn", "nugent-sn2n", -18.53),
]
ZS_FACTORS = np.arange(1.0, 3.01, 0.25)  # z_s / z_l grid for lenses without a measured source redshift
L_MU_REF = None  # set by mu_ref_unknown()


def season_dev(r: np.ndarray, e: np.ndarray, jd: np.ndarray, t_w: float) -> tuple[float, float]:
    """Largest |season mean - median season mean| (microJansky) and its significance, outside the event span."""
    keep = (jd < t_w - 100) | (jd > t_w + 250)
    r, e, jd = r[keep], e[keep], jd[keep]
    yrs = np.floor((jd - 2458000) / 365.25)
    means, errs = [], []
    for y in np.unique(yrs):
        m = yrs == y
        if m.sum() < 5:
            continue
        w = 1 / e[m] ** 2
        means.append(float((w * r[m]).sum() / w.sum()))
        errs.append(float(1 / np.sqrt(w.sum())))
    if len(means) < 3:
        return 0.0, 0.0
    means, errs = np.array(means), np.array(errs)
    med = np.median(means)
    dev = np.abs(means - med)
    sig = dev / errs
    # the season that deviates most in amplitude among those significant
    ok = sig > VAR_SIG
    amp = float(dev[ok].max()) if ok.any() else 0.0
    return amp, float(sig.max())


def variable(r, e, jd, t_w) -> bool:
    amp, _ = season_dev(r, e, jd, t_w)
    return amp > VAR_AMP


@lru_cache(maxsize=None)
def _trimmed(band: str, lo: float, hi: float) -> sncosmo.Bandpass:
    b = sncosmo.get_bandpass(band)
    m = (b.wave > lo + 1) & (b.wave < hi - 1)
    return sncosmo.Bandpass(b.wave[m], b.trans[m], name=f"{band}_{lo:.0f}")


def bandflux(m: sncosmo.Model, band: str, t: np.ndarray) -> np.ndarray:
    """ZTF bandflux in microJansky; where the model's SED does not cover the whole band (the rest-frame far UV of a
    high-z source in g), the flux through the covered part of the band (0 if none)."""
    try:
        return m.bandflux(band, t, zp=23.9, zpsys="ab")
    except ValueError:
        b = sncosmo.get_bandpass(band)
        lo, hi = max(m.minwave(), b.minwave()), min(m.maxwave(), b.maxwave())
        if hi - lo < 200:
            return np.zeros_like(np.asarray(t, float))
        tb = _trimmed(band, float(lo), float(hi))
        # keep the AB zero point of the full band so the trimmed flux is the covered part of the full band's flux
        full = sncosmo.get_magsystem("ab").zpbandflux(b)
        part = sncosmo.get_magsystem("ab").zpbandflux(tb)
        return m.bandflux(tb, t, zp=23.9, zpsys="ab") * part / full


@lru_cache(maxsize=None)
def template_grid(template: str, z: float, mb: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Unit-magnification ZTF r and g fluxes (microJansky) of one SN image vs observer-frame days from t0."""
    src = sncosmo.get_source(template)
    m = sncosmo.Model(source=src)
    if template == "salt3":
        m.set(z=z, x1=0.0, c=0.0)
    else:
        m.set(z=z)
    m.set_source_peakabsmag(mb, "bessellb", "ab", cosmo=COSMO)
    m.set(t0=0.0)
    t = np.arange(np.ceil(m.mintime()), np.floor(m.maxtime()), 1.0)
    fr = bandflux(m, "ztfr", t)
    fg = bandflux(m, "ztfg", t)
    return t, np.clip(fr, 0, None), np.clip(fg, 0, None)


def mu_ref_unknown(rng_seed: int = 141) -> float:
    """Median total magnification of the injection draw for lenses without a model (2 or 4 images, 0.7/0.3,
    ln mu ~ N(1.5, 0.35) each)."""
    rng = np.random.default_rng(rng_seed)
    tot = [np.exp(rng.normal(1.5, 0.35, 2 if rng.random() < 0.7 else 4)).sum() for _ in range(20000)]
    return float(np.median(tot))


def hypotheses_for(lens: pd.Series) -> tuple[list[tuple[str, str, float, float]], float]:
    global L_MU_REF
    if lens["class"] == "S" and np.isfinite(lens.z_s):
        zs = [float(lens.z_s)]
    else:
        zl = float(lens.z_l) if np.isfinite(lens.z_l) else 0.45
        zs = [round(zl * f, 4) for f in ZS_FACTORS]
    if lens["class"] == "S" and np.isfinite(lens.mu) and lens.mu > 0:
        mu = float(lens.mu)
    else:
        if L_MU_REF is None:
            L_MU_REF = mu_ref_unknown()
        mu = L_MU_REF
    return [(n, t, mb, z) for z in zs for (n, t, mb) in HYPOTHESES], mu


def colour_mag(r, e, jd, k0, g, eg, jdg, k0g, t_w: float, lens: pd.Series) -> dict:
    """Stage 4. r, g: baseline-subtracted fluxes; k0, k0g: per-epoch centred-source responses (kappa0)."""
    hyps, mu_ref = hypotheses_for(lens)
    sel = (jd >= t_w - 60) & (jd <= t_w + 120)
    selg = (jdg >= t_w - 60) & (jdg <= t_w + 120) if jdg is not None else None
    best = {"pass": False, "hyp": "", "A": np.nan, "mu_ref": mu_ref, "n_g": int(selg.sum()) if selg is not None else 0,
            "mag_ok_any": False}
    if sel.sum() < 2:
        return best
    rr, ee, tt, kk = r[sel], e[sel], jd[sel], k0[sel]
    w = 1 / ee**2
    t0s = np.arange(t_w - 90, t_w + 31, 2.0)
    for name, tmpl, mb, z in hyps:
        tg, fr, fg = template_grid(tmpl, z, mb)
        bestchi, bA, bsA, bt0 = np.inf, np.nan, np.nan, np.nan
        for t0 in t0s:
            m = np.interp(tt - t0, tg, fr, left=0, right=0) * kk
            sw = (w * m * m).sum()
            if sw <= 0:
                continue
            A = (w * m * rr).sum() / sw
            chi = (w * (rr - A * m) ** 2).sum()
            if chi < bestchi:
                bestchi, bA, bsA, bt0 = chi, A, 1 / np.sqrt(sw), t0
        if not np.isfinite(bA):
            continue
        mag_ok = MU_LO * mu_ref <= bA <= MU_HI * mu_ref
        best["mag_ok_any"] = best["mag_ok_any"] or mag_ok
        if not mag_ok:
            continue
        col_ok = True
        if COLOUR_MODE == "chi2" and selg is not None and selg.sum() > 0:
            mg = np.interp(jdg[selg] - bt0, tg, fg, left=0, right=0) * k0g[selg]
            var = eg[selg] ** 2 + (bsA * mg) ** 2
            chi_g = float(((g[selg] - bA * mg) ** 2 / var).sum())
            col_ok = stats.chi2.sf(chi_g, selg.sum()) > COLOUR_P
        elif COLOUR_MODE == "peak" and jdg is not None:
            pk = (jdg >= t_w) & (jdg < t_w + 2 * WIN)
            if pk.any():
                mg = np.interp(jdg[pk] - bt0, tg, fg, left=0, right=0) * k0g[pk]
                wg = 1 / (eg[pk] ** 2 + (bsA * mg) ** 2)
                zc = ((wg * (g[pk] - bA * mg)).sum() / wg.sum()) * np.sqrt(wg.sum())
                col_ok = abs(zc) < COLOUR_Z
        if col_ok:
            best.update({"pass": True, "hyp": f"{name}@z{z:g}", "A": bA})
            return best
    return best


def vet(w: pd.DataFrame, thr: float, ctrl_S_in_window, centroid_off: float, r, e, jd, k0, g, eg, jdg, k0g,
        lens: pd.Series, w0_restricted: pd.DataFrame | None = None) -> dict:
    """All stages for one light curve. w: its window table (already restricted in time for injections);
    ctrl_S_in_window(t_w) -> max control S in that window; centroid_off(t_w, S) -> arcsec."""
    out = {"T": np.nan, "T2": np.nan, "t_w": np.nan, "s0": False, "s1": False, "s2": False, "s3": False, "s4": False,
           "hyp": "", "A": np.nan, "var_amp": np.nan}
    s = stat(w)
    out["T"], out["T2"] = s["T"], s["T2"]
    if not (np.isfinite(s["T2"]) and s["T2"] > thr):
        return out
    if w0_restricted is not None:
        s0 = stat(w0_restricted)
        if np.isfinite(s0["T2"]) and s0["T2"] > thr:
            return out  # R-135's amendment: no credit during a real event
    out["s0"] = True
    # the detection window: the best window that satisfies the shape cut
    S, t0 = w.S.values, w.t0.values
    best, t_w = -np.inf, np.nan
    for j in range(len(w)):
        nb = (np.abs(t0 - t0[j] - WIN) < 1e-6) | (np.abs(t0 - t0[j] + WIN) < 1e-6)
        if nb.any() and S[nb].max() > 3.0 and S[j] > best:
            best, t_w = S[j], t0[j]
    out["t_w"] = t_w
    if ctrl_S_in_window(t_w) > CTRL_S_MAX:
        return out
    out["s1"] = True
    if not centroid_off(t_w, best) <= CENTROID_MAX:
        return out
    out["s2"] = True
    amp, _ = season_dev(r, e, jd, t_w)
    out["var_amp"] = amp
    if amp > VAR_AMP:
        return out
    out["s3"] = True
    cm = colour_mag(r, e, jd, k0, g, eg, jdg, k0g, t_w, lens)
    out["hyp"], out["A"] = cm["hyp"], cm["A"]
    out["s4"] = bool(cm["pass"])
    return out
