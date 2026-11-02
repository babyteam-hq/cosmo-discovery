"""Capture-recapture across gravitational-wave search catalogues (R-126, cosmo-deep-6).

The searches run on the same strain are treated as overlapping, imperfect registries of one population of real signals.
From the pattern of which lists found which candidate above a common threshold, estimate the number of real signals the
union of lists could have found, per observing run and per chirp-mass bin, with no injections.

The referee's caveat is built in: the lists share the strain and (partly) the templates, so their detections are
strongly positively dependent through the signal-to-noise ratio. Closed-population estimators that assume independence
(Lincoln-Petersen, the sample-coverage estimate, the independence log-linear model) are biased low under that
dependence, and log-linear interaction terms cannot identify a dependence that runs through one shared variable. So the
headline estimator conditions on the shared variable: within the observed network SNR (and chirp mass) each list's
catchability is modelled as a logistic function, the lists are taken as conditionally independent given SNR and chirp
mass, and the total is the Horvitz-Thompson sum over the observed candidates of 1 / P(found by at least one list).
Residual dependence given SNR is measured (Mantel-Haenszel odds ratio across SNR strata); where it is positive the
estimate is a lower bound, and Chao's lower bound is reported beside it.

Usage:
  venv/bin/python gw_capture_recapture.py            the data: all runs, both thresholds, chirp-mass bins
  venv/bin/python gw_capture_recapture.py --mock     the mock populations with known truth
  venv/bin/python gw_capture_recapture.py --vt       the per-pipeline VT ratios from the GWTC-2.1 / GWTC-3 tables
"""
from __future__ import annotations

import itertools
import math
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import optimize, stats, special

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "gw")
RNG = np.random.default_rng(126)

RUNS = {
    # run: (csv, lists, the lists that are LVK pipelines)
    "O1O2": ("lists_O1O2.csv", ["pycbc", "gstlal", "cwb", "ogc", "ias"], ["pycbc", "gstlal", "cwb"]),
    "O3": ("lists_O3.csv", ["cwb", "gstlal", "mbta", "pycbc", "pycbcbbh", "ogc", "ias"], ["cwb", "gstlal", "mbta", "pycbc", "pycbcbbh"]),
    "O4a": ("lists_O4a.csv", ["cwb", "gstlal", "mbta", "pycbc"], ["cwb", "gstlal", "mbta", "pycbc"]),
}
LABEL = {"cwb": "cWB", "gstlal": "GstLAL", "mbta": "MBTA", "pycbc": "PyCBC", "pycbcbbh": "PyCBC-BBH", "ogc": "4-OGC", "ias": "IAS"}


# ----------------------------------------------------------------------------------------------------------------
# Capture histories
# ----------------------------------------------------------------------------------------------------------------
def load(run: str) -> pd.DataFrame:
    csv, lists, _ = RUNS[run]
    df = pd.read_csv(os.path.join(DATA, csv))
    snrs = df[[c + "_snr" for c in lists if c + "_snr" in df.columns]]
    df["snr"] = snrs.max(axis=1)
    if "ias_snr" in df.columns:
        df["snr"] = df["snr"].where(df["snr"].notna(), df["ias_snr"])
    if "ogc_snr" in df.columns:
        df["snr"] = df["snr"].where(df["snr"].notna(), df["ogc_snr"])
    df["inst"] = df["inst"].fillna("").astype(str)
    return df


def histories(df: pd.DataFrame, lists: list[str], crit: str, far_thr: float = 1.0, pastro_thr: float = 0.5) -> np.ndarray:
    Y = np.zeros((len(df), len(lists)), dtype=bool)
    for j, L in enumerate(lists):
        if crit == "far":
            col = L + "_far"
            if col in df.columns:
                Y[:, j] = df[col].to_numpy() < far_thr
        else:
            col = L + "_pastro"
            if col in df.columns:
                Y[:, j] = df[col].to_numpy() > pastro_thr
    return Y


def select(df: pd.DataFrame, hl: bool, bbh: bool, mc_range=None, snr_range=None) -> pd.Series:
    m = pd.Series(True, index=df.index)
    if hl:
        m &= df["inst"].str.contains("H") & df["inst"].str.contains("L")
    if bbh:
        # BBH by the catalogue's template / PE chirp mass (> 5 Msun, i.e. both components above ~3 Msun in the usual
        # banks) or, where a p_BBH exists, by it
        mc = df["mchirp"]
        pb = df["pbbh"] if "pbbh" in df.columns else pd.Series(np.nan, index=df.index)
        m &= ((mc > 5.0) & ~(pb < 0.5)) | ((pb >= 0.5) & ~(mc <= 5.0))
    if mc_range is not None:
        m &= (df["mchirp"] >= mc_range[0]) & (df["mchirp"] < mc_range[1])
    if snr_range is not None:
        m &= (df["snr"] >= snr_range[0]) & (df["snr"] < snr_range[1])
    return m


# ----------------------------------------------------------------------------------------------------------------
# Estimators on a capture-history matrix Y (n x K, every row has >= 1 capture)
# ----------------------------------------------------------------------------------------------------------------
def chapman(Y: np.ndarray, i: int, j: int) -> tuple[float, float]:
    n1, n2 = Y[:, i].sum(), Y[:, j].sum()
    m = (Y[:, i] & Y[:, j]).sum()
    N = (n1 + 1) * (n2 + 1) / (m + 1) - 1
    var = (n1 + 1) * (n2 + 1) * (n1 - m) * (n2 - m) / ((m + 1) ** 2 * (m + 2))
    return N, math.sqrt(max(var, 0))


def chao_lower(Y: np.ndarray) -> tuple[float, float]:
    """Chao (1987) lower bound for N from the frequency counts of units seen by exactly 1 and 2 lists."""
    k = Y.sum(axis=1)
    D = len(k)
    f1, f2 = int((k == 1).sum()), int((k == 2).sum())
    if f2 > 0:
        N = D + f1 * f1 / (2 * f2)
        var = f2 * ((f1 / f2) ** 4 / 4 + (f1 / f2) ** 3 + (f1 / f2) ** 2 / 2)
    else:
        N = D + f1 * (f1 - 1) / 2
        var = f1 * (f1 - 1) / 2 + f1 * (2 * f1 - 1) ** 2 / 4 - f1**4 / (4 * max(N, 1))
    return N, math.sqrt(max(var, 0))


def sample_coverage(Y: np.ndarray) -> float:
    """Chao & Tsay (1998) sample-coverage estimate under independence: D / C, C = 1 - f1 / sum_i n_i."""
    k = Y.sum(axis=1)
    f1 = (k == 1).sum()
    C = 1 - f1 / Y.sum()
    return len(k) / C if C > 0 else float("inf")


def loglinear(Y: np.ndarray, order: int) -> tuple[float, float, float]:
    """Poisson log-linear model on the 2^K - 1 observed cells with main effects (order 1) or all pairwise
    interactions (order 2). Returns (N, se, AIC); N = D + fitted count of the unobserved all-zero cell."""
    import statsmodels.api as sm

    K = Y.shape[1]
    cells = np.array(list(itertools.product([0, 1], repeat=K)))[1:]
    counts = np.array([(Y == c.astype(bool)).all(axis=1).sum() for c in cells])
    cols = [cells[:, i] for i in range(K)]
    zero = [np.zeros(1) for _ in range(K)]
    if order >= 2:
        for i, j in itertools.combinations(range(K), 2):
            cols.append(cells[:, i] * cells[:, j])
            zero.append(np.zeros(1))
    X = sm.add_constant(np.column_stack(cols), has_constant="add")
    X0 = sm.add_constant(np.column_stack(zero), has_constant="add")
    try:
        res = sm.GLM(counts, X, family=sm.families.Poisson()).fit(maxiter=200)
        eta0 = float((X0 @ res.params)[0])
        n0 = math.exp(min(eta0, 20))
        se0 = math.sqrt(max(float((X0 @ np.asarray(res.cov_params()) @ X0.T)[0, 0]), 0))
        return len(Y) + n0, n0 * se0, res.aic
    except Exception:
        return float("nan"), float("nan"), float("nan")


def mantel_haenszel(Y: np.ndarray, strata: np.ndarray, i: int, j: int) -> tuple[float, int]:
    """Pooled odds ratio of 'found by i' vs 'found by j' across strata, conditioning on being in the union (lists i
    and j only, among units seen by any list). > 1: positive residual dependence given the stratum variable."""
    num = den = 0.0
    used = 0
    for s in np.unique(strata):
        m = strata == s
        a = (Y[m, i] & Y[m, j]).sum()
        b = (Y[m, i] & ~Y[m, j]).sum()
        c = (~Y[m, i] & Y[m, j]).sum()
        d = (~Y[m, i] & ~Y[m, j]).sum()
        n = a + b + c + d
        if n < 2 or (a + b) == 0 or (c + d) == 0 or (a + c) == 0 or (b + d) == 0:
            continue
        num += a * d / n
        den += b * c / n
        used += 1
    return (num / den if den > 0 else float("inf")), used


def fit_conditional(Y: np.ndarray, X: np.ndarray, w: np.ndarray | None = None) -> np.ndarray:
    """Huggins-Alho conditional likelihood: p_ij = expit(X_i . beta_j) per list j, conditioned on >= 1 capture.
    X is n x d (with a constant column); w weights each candidate's term (its probability of being a real signal).
    Returns beta (K x d)."""
    n, K = Y.shape
    d = X.shape[1]
    if w is None:
        w = np.ones(n)

    def nll(theta):
        B = theta.reshape(K, d)
        eta = X @ B.T  # n x K
        logp = -np.logaddexp(0, -eta)
        logq = -np.logaddexp(0, eta)
        ll_hist = (Y * logp + (~Y) * logq).sum(axis=1)
        log_pi = np.log1p(-np.exp(logq.sum(axis=1)) + 1e-300)
        pen = 1e-3 * (theta**2).sum()  # a weak ridge keeps separated lists finite
        return -(w * (ll_hist - log_pi)).sum() + pen

    theta0 = np.zeros(K * d)
    theta0[::d] = 0.0
    res = optimize.minimize(nll, theta0, method="L-BFGS-B")
    return res.x.reshape(K, d)


def horvitz_thompson(Y: np.ndarray, X: np.ndarray, B: np.ndarray) -> tuple[float, np.ndarray]:
    eta = X @ B.T
    q = special.expit(-eta)  # 1 - p
    pi = 1 - np.prod(q, axis=1)
    pi = np.clip(pi, 1e-6, 1)
    return float((1 / pi).sum()), pi


def design(snr: np.ndarray, mc: np.ndarray | None, snr0: float = 9.0) -> np.ndarray:
    cols = [np.ones_like(snr), np.clip(snr, 6, 40) - snr0]
    if mc is not None:
        cols.append(np.log(np.clip(mc, 1, 500)) - math.log(20))
    return np.column_stack(cols)


def ht_estimate(Y: np.ndarray, snr: np.ndarray, mc: np.ndarray | None, nboot: int = 300, w=None, masks: dict | None = None) -> dict:
    """The conditional-independence (given SNR, chirp mass) estimate with a bootstrap over candidates. w: each
    candidate's weight (its probability of being real); the total is sum of w / pi. masks: name -> boolean mask of the
    candidates whose w / pi are summed (a chirp-mass or SNR bin); the fit always uses every candidate. Returns
    {"N", "lo", "hi", "pi", "B", "bins": {name: (N, lo, hi, n)}}."""
    n = len(Y)
    if w is None:
        w = np.ones(n)
    masks = masks or {}
    X = design(snr, mc)
    B = fit_conditional(Y, X, w)
    N, pi = horvitz_thompson(Y, X, B)
    Ntot = float((w / pi).sum())
    boots = []
    bboots = {k: [] for k in masks}
    for _ in range(nboot):
        idx = RNG.integers(0, n, n)
        try:
            Bb = fit_conditional(Y[idx], X[idx], w[idx])
            _, pib = horvitz_thompson(Y[idx], X[idx], Bb)
        except Exception:
            continue
        boots.append(float((w[idx] / pib).sum()))
        for k, m in masks.items():
            mb = m[idx]
            bboots[k].append(float((w[idx][mb] / pib[mb]).sum()))
    boots = np.array(boots)
    lo, hi = (np.percentile(boots, [16, 84]) if len(boots) > 10 else (float("nan"), float("nan")))
    bins = {}
    for k, m in masks.items():
        bb = np.array(bboots[k])
        blo, bhi = (np.percentile(bb, [16, 84]) if len(bb) > 10 else (float("nan"), float("nan")))
        bins[k] = (float((w[m] / pi[m]).sum()), blo, bhi, int(m.sum()))
    return {"N": Ntot, "lo": lo, "hi": hi, "pi": pi, "B": B, "bins": bins}


# ----------------------------------------------------------------------------------------------------------------
# The analysis of one sample
# ----------------------------------------------------------------------------------------------------------------
def analyse(df: pd.DataFrame, lists: list[str], lvk: list[str], crit: str, label: str, mc_bins=None, nboot=300) -> list[str]:
    out = []
    Y = histories(df, lists, crit)
    keep = Y.any(axis=1)
    if keep.sum() < 5:
        return [f"**{label}**, threshold {crit}: {int(keep.sum())} candidates in the union; not analysed (GWTC-1 publishes no per-pipeline p_astro)"]
    df = df[keep].reset_index(drop=True)
    Y = Y[keep]
    snr = df["snr"].to_numpy(float)
    snr = np.where(np.isnan(snr), np.nanmedian(snr), snr)
    mc = df["mchirp"].to_numpy(float)
    mc = np.where(np.isnan(mc), np.nanmedian(mc), mc)
    n = len(df)
    k = Y.sum(axis=1)
    jl = [lists.index(L) for L in lvk]
    n_lvk = int(Y[:, jl].any(axis=1).sum())
    pcols = [L + "_pastro" for L in lists if L + "_pastro" in df.columns]
    w = df[pcols].max(axis=1).to_numpy(float) if pcols else np.ones(n)
    w = np.where(np.isnan(w), 1.0, np.clip(w, 0, 1))
    w_lvk = float(w[Y[:, jl].any(axis=1)].sum())
    out.append(f"**{label}**, threshold {'FAR < 1/yr' if crit == 'far' else 'p_astro > 0.5'}: union {n} (expected real, sum of max p_astro: {w.sum():.1f}), LVK union {n_lvk} ({w_lvk:.1f}); "
               + ", ".join(f"{LABEL[L]} {int(Y[:, j].sum())}" for j, L in enumerate(lists))
               + f"; seen by exactly one list {int((k == 1).sum())}, two {int((k == 2).sum())}, all {int((k == len(lists)).sum())}")
    # naive estimators
    chaoN, chaoS = chao_lower(Y)
    sc = sample_coverage(Y)
    ll1 = loglinear(Y, 1)
    ll2 = loglinear(Y, 2) if len(lists) <= 7 else (float("nan"),) * 3
    pairs = []
    for i, j in itertools.combinations(range(len(lists)), 2):
        N, s = chapman(Y, i, j)
        pairs.append((N, s, LABEL[lists[i]], LABEL[lists[j]]))
    pairs.sort()
    out.append(f"  pooled (dependence ignored): Chao lower bound {chaoN:.1f} ± {chaoS:.1f}; sample coverage {sc:.1f}; "
               f"log-linear independence {ll1[0]:.1f} ± {ll1[1]:.1f} (AIC {ll1[2]:.0f}), pairwise {ll2[0]:.1f} ± {ll2[1]:.1f} (AIC {ll2[2]:.0f}); "
               f"Lincoln-Petersen pairs {pairs[0][0]:.1f} ({pairs[0][2]}-{pairs[0][3]}) to {pairs[-1][0]:.1f} ({pairs[-1][2]}-{pairs[-1][3]})")
    # stratified: residual dependence given SNR
    strata = np.digitize(snr, [8.5, 9.5, 10.5, 12.0])
    ors = []
    for i, j in itertools.combinations(range(len(lists)), 2):
        o, used = mantel_haenszel(Y, strata, i, j)
        if used >= 2 and np.isfinite(o) and o > 0:
            ors.append((o, LABEL[lists[i]], LABEL[lists[j]]))
    if ors:
        med = np.median([o for o, _, _ in ors])
        ors.sort()
        out.append(f"  residual pairwise dependence given SNR (Mantel-Haenszel odds ratio over 5 SNR strata; 1 = none): median {med:.2f}, "
                   f"range {ors[0][0]:.2f} ({ors[0][1]}-{ors[0][2]}) to {ors[-1][0]:.2f} ({ors[-1][1]}-{ors[-1][2]}), {len(ors)} pairs")
    # the conditional-independence estimate
    bins = [(0, 8.5), (8.5, 9.5), (9.5, 10.5), (10.5, 12), (12, 99)]
    masks = {f"SNR {lo:g}-{hi:g}": (snr >= lo) & (snr < hi) for lo, hi in bins}
    if mc_bins:
        masks.update({f"Mc {lo:g}-{hi:g}": (mc >= lo) & (mc < hi) for lo, hi in mc_bins})
    masks = {k: m for k, m in masks.items() if m.sum() >= 1}
    res = ht_estimate(Y, snr, mc, nboot=nboot, masks=masks)
    out.append(f"  **conditional on SNR and chirp mass (Horvitz-Thompson): N = {res['N']:.1f} (68 % bootstrap {res['lo']:.1f}-{res['hi']:.1f})**; "
               f"union completeness {n / res['N']:.2f}, LVK completeness {n_lvk / res['N']:.2f}")
    res_s = ht_estimate(Y, snr, None, nboot=nboot)
    out.append(f"  conditional on SNR alone: N = {res_s['N']:.1f} ({res_s['lo']:.1f}-{res_s['hi']:.1f})")
    resw = ht_estimate(Y, snr, mc, nboot=nboot, w=w, masks=masks)
    out.append(f"  **weighted by p_astro (expected real signals): N_w = {resw['N']:.1f} ({resw['lo']:.1f}-{resw['hi']:.1f})**; union holds {w.sum():.1f} expected real "
               f"({w.sum() / resw['N']:.2f}), LVK {w_lvk:.1f} ({w_lvk / resw['N']:.2f})")
    # per-list catchability at SNR 8, 9, 10, 12 (chirp mass 20) and the implied miss probability of the union
    B = res["B"]
    rows = []
    for s in (8.0, 9.0, 10.0, 12.0):
        x = design(np.array([s]), np.array([20.0]))
        p = special.expit(x @ B.T)[0]
        rows.append(f"SNR {s:.0f}: " + ", ".join(f"{LABEL[L]} {p[j]:.2f}" for j, L in enumerate(lists)) + f", any {1 - np.prod(1 - p):.2f}")
    out.append("  catchability at chirp mass 20: " + "; ".join(rows))
    # SNR bins: observed vs estimated
    cells = []
    for lo, hi in bins:
        k = f"SNR {lo:g}-{hi:g}"
        if k not in masks:
            continue
        m = masks[k]
        r, rw = res["bins"][k], resw["bins"][k]
        cells.append(f"{k}: {int(m.sum())} seen ({w[m].sum():.1f} real), {int(Y[m][:, jl].any(axis=1).sum())} LVK, N {r[0]:.1f} ({r[1]:.1f}-{r[2]:.1f}), N_w {rw[0]:.1f} ({rw[1]:.1f}-{rw[2]:.1f})")
    out.append("  by SNR bin: " + "; ".join(cells))
    # chirp-mass bins: HT sum within bins (fit on all), Chao within bins
    if mc_bins:
        cells = []
        for lo, hi in mc_bins:
            k = f"Mc {lo:g}-{hi:g}"
            if k not in masks or masks[k].sum() < 3:
                continue
            m = masks[k]
            r, rw = res["bins"][k], resw["bins"][k]
            cN, cS = chao_lower(Y[m])
            cells.append(f"{k}: {int(m.sum())} seen ({w[m].sum():.1f} real), {int(Y[m][:, jl].any(axis=1).sum())} LVK ({w[m & Y[:, jl].any(axis=1)].sum():.1f}), "
                         + ", ".join(f"{LABEL[L]} {int(Y[m, j].sum())}" for j, L in enumerate(lists))
                         + f"; N_HT {r[0]:.1f} ({r[1]:.1f}-{r[2]:.1f}), N_w {rw[0]:.1f} ({rw[1]:.1f}-{rw[2]:.1f}), Chao {cN:.1f} ± {cS:.1f}")
        out.append("  by chirp-mass bin: " + " | ".join(cells))
    return out


def run_data(nboot: int = 300) -> list[str]:
    out = []
    mc_bins = [(0, 5), (5, 15), (15, 30), (30, 60), (60, 1000)]
    for run, (csv, lists, lvk) in RUNS.items():
        df = load(run)
        for crit in ("far", "pastro"):
            if run == "O3":
                # the seven-list population: Hanford-Livingston coincident BBH (what IAS and 4-OGC searched)
                m = select(df, hl=True, bbh=True)
                out += analyse(df[m].reset_index(drop=True), lists, lvk, crit, f"{run} HL BBH, 7 lists", mc_bins=mc_bins[1:], nboot=nboot)
                for sub in ("O3a", "O3b"):
                    ms = m & (df["run"] == sub)
                    out += analyse(df[ms].reset_index(drop=True), lists, lvk, crit, f"{sub} HL BBH, 7 lists", nboot=nboot)
                # the LVK-only lists on every candidate (all masses, all networks)
                out += analyse(df, lvk, lvk, crit, f"{run} all candidates, 5 LVK lists", mc_bins=mc_bins, nboot=nboot)
            elif run == "O1O2":
                m = select(df, hl=True, bbh=True)
                out += analyse(df[m].reset_index(drop=True), lists, lvk, crit, f"{run} HL BBH, 5 lists", nboot=nboot)
                out += analyse(df, lvk, lvk, crit, f"{run} all candidates, 3 LVK lists", nboot=nboot)
            else:
                out += analyse(df, lists, lvk, crit, f"{run} all candidates, 4 LVK lists", mc_bins=mc_bins, nboot=nboot)
                m = select(df, hl=False, bbh=True)
                out += analyse(df[m].reset_index(drop=True), lists, lvk, crit, f"{run} BBH, 4 LVK lists", mc_bins=mc_bins[1:], nboot=nboot)
    return out


# ----------------------------------------------------------------------------------------------------------------
# Mocks with known truth
# ----------------------------------------------------------------------------------------------------------------
def mock_population(n_true: int, K: int, thresholds, sigma_shared: float, sigma_list, miss_prob, mc_dep=None):
    """Signals with true network SNR ~ rho^-4 above 5. Every list sees base = rho + a shared noise term (the common
    strain realisation) plus its own noise, and reports the signal when that exceeds its threshold and it was not
    independently dropped (miss_prob: vetoes, gaps). mc_dep: (K,) extra threshold per e-fold of lower chirp mass (a
    list that needs louder low-mass signals). Returns base, log chirp mass, the capture matrix and each list's
    reported SNR (NaN where not found)."""
    u = RNG.random(n_true)
    rho = 5.0 * (1 - u) ** (-1 / 3)  # pdf ∝ rho^-4 above 5
    rho = rho[rho < 60]
    n = len(rho)
    logmc = RNG.normal(math.log(20), 0.6, n)
    base = rho + RNG.normal(0, sigma_shared, n)
    Y = np.zeros((n, K), bool)
    snr_rep = np.full((n, K), np.nan)
    for j in range(K):
        own = RNG.normal(0, sigma_list[j], n)
        thr = thresholds[j] + (mc_dep[j] * (math.log(20) - logmc) if mc_dep is not None else 0)
        stat = base + own
        det = (stat > thr) & (RNG.random(n) > miss_prob[j])
        Y[:, j] = det
        snr_rep[det, j] = stat[det]
    return base, logmc, Y, snr_rep


MOCK_BINS = [(8.0, 8.5), (8.5, 9.5), (9.5, 10.5), (10.5, 99)]


def run_mock(nboot: int = 100, reps: int = 3) -> list[str]:
    out = ["Mock populations (truth known). Lists see the same base SNR (true SNR plus a shared noise term, the common strain) "
           "plus their own noise; thresholds differ per list; each list also drops a random fraction (vetoes, gaps). The "
           "estimators see only the union of the lists and each candidate's reported SNR (the maximum over the lists that "
           "found it, as in the data). Truth per bin: the number of signals, found or not, whose base SNR lies in the bin; "
           "'truth >= 8' is the population a union of lists with thresholds near 8 could have found."]
    configs = [
        ("A: 5 lists, strong shared noise (sigma_shared 0.6, own 0.3), thresholds 8.0-9.0, 5 % random misses", 5, [8.0, 8.3, 8.6, 9.0, 8.4], 0.6, [0.3] * 5, [0.05] * 5, None),
        ("B: as A with no shared noise (conditional independence given the base holds exactly)", 5, [8.0, 8.3, 8.6, 9.0, 8.4], 0.0, [0.5] * 5, [0.05] * 5, None),
        ("C: as A with 20 % random misses per list", 5, [8.0, 8.3, 8.6, 9.0, 8.4], 0.6, [0.3] * 5, [0.2] * 5, None),
        ("D: as A, one list needs 2 more SNR per e-fold of lower chirp mass (a cWB-like list)", 5, [8.0, 8.3, 8.6, 9.0, 8.4], 0.6, [0.3] * 5, [0.05] * 5, [0, 0, 0, 2.0, 0]),
        ("E: 7 lists as A, thresholds 8.0-9.3 (the O3 configuration)", 7, [8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.3], 0.6, [0.3] * 7, [0.05] * 7, None),
        ("F: as A with a 3000-signal population (catalogue-sized samples)", 5, [8.0, 8.3, 8.6, 9.0, 8.4], 0.6, [0.3] * 5, [0.05] * 5, None),
    ]
    for label, K, thr, ss, sl, mp, mcd in configs:
        for rep in range(reps if label[0] != "F" else max(reps, 5)):
            base, logmc, Y, snr_rep = mock_population(3000 if label[0] == "F" else 20000, K, thr, ss, sl, mp, mcd)
            seen = Y.any(axis=1)
            snr = np.nanmax(np.where(np.isnan(snr_rep), -np.inf, snr_rep), axis=1)
            Ys, snrs, mcs = Y[seen], snr[seen], np.exp(logmc[seen])
            n_union = int(seen.sum())
            chaoN, chaoS = chao_lower(Ys)
            ll1 = loglinear(Ys, 1)
            ll2 = loglinear(Ys, 2)
            res = ht_estimate(Ys, snrs, mcs, nboot=nboot)
            # the oracle: the same estimator conditioned on the true base SNR (the shared variable itself) instead of
            # the selection-biased reported maximum; truth per bin of base among the seen vs all
            bases = base[seen]
            reso = ht_estimate(Ys, bases, mcs, nboot=max(nboot // 5, 5))
            ocells = []
            for lo, hi in MOCK_BINS:
                mo = (bases >= lo) & (bases < hi)
                t = int(((base >= lo) & (base < hi)).sum())
                ocells.append(f"{lo:g}-{hi:g}: seen {int(mo.sum())}, HT|base {(1 / reso['pi'][mo]).sum():.0f}, truth {t}")
            strata = np.digitize(snrs, [8.5, 9.5, 10.5, 12.0])
            ors = [mantel_haenszel(Ys, strata, i, j)[0] for i, j in itertools.combinations(range(K), 2)]
            ors = [o for o in ors if np.isfinite(o) and o > 0]
            truth8 = int((base >= 8.0).sum())
            cells = []
            for lo, hi in MOCK_BINS:
                m = (snrs >= lo) & (snrs < hi)
                t = int(((base >= lo) & (base < hi)).sum())
                cells.append(f"SNR {lo:g}-{hi:g}: seen {int(m.sum())}, HT {(1 / res['pi'][m]).sum():.0f}, truth {t}")
            out.append(f"- {label} (rep {rep + 1}): union {n_union}, truth (base >= 8) {truth8}, truth (base >= 7.5) {int((base >= 7.5).sum())}; "
                       f"Chao {chaoN:.0f} ± {chaoS:.0f}; log-linear independence {ll1[0]:.0f}, pairwise {ll2[0]:.0f}; "
                       f"HT|SNR,Mc {res['N']:.0f} ({res['lo']:.0f}-{res['hi']:.0f}); MH odds ratio median {np.median(ors) if ors else float('nan'):.2f}; "
                       + "; ".join(cells) + " | oracle (base SNR as covariate): " + "; ".join(ocells))
    out.append("The null of a nested list: list 2 is a random 60 % subsample of list 1 (independent catchability given list 1): " + nested_null())
    return out


def nested_null() -> str:
    base, logmc, Y, snr_rep = mock_population(20000, 1, [8.0], 0.0, [0.3], [0.0], None)
    seen = Y[:, 0]
    n1 = int(seen.sum())
    sub = seen & (RNG.random(len(base)) < 0.6)
    Y2 = np.column_stack([seen, sub])[seen]
    N, s = chapman(Y2, 0, 1)
    return f"list 1 {n1}, list 2 {int(sub.sum())}, Lincoln-Petersen {N:.0f} ± {s:.0f} (truth {n1}); a list nested in another returns the larger list, as it must"


# ----------------------------------------------------------------------------------------------------------------
# VT ratios from the catalogue papers (injection-based) for comparison with the catchabilities
# ----------------------------------------------------------------------------------------------------------------
def run_vt() -> list[str]:
    sys.path.insert(0, HERE)
    import gw_lists as g

    out = []
    vt3 = g.parse_macros(os.path.join(g.SRC, "2111.03606", "vt_macros.tex"))
    keys = ["THIRTYFIVETHIRTYFIVE", "THIRTYFIVETWENTY", "TWENTYTWENTY", "TWENTYTEN", "TENFIVE", "FIVEFIVE", "THIRTYFIVEONEPOINTFIVE", "ONEPOINTFIVEONEPOINTFIVE"]
    pipes = [("cwb", "CWBALLSKY"), ("gstlal", "GSTLALALLSKY"), ("mbta", "MBTAALLSKY"), ("pycbc", "PYCBCALLSKY"), ("pycbcbbh", "PYCBCHIGHMASS")]
    out.append("GWTC-3 Table (O3b, FAR < 1/yr... the table's threshold): sensitive hypervolume per pipeline / 'any pipeline', fixed-mass injections")
    for key in keys:
        bins = vt3["VTBINS"].get(key, key)
        anyv = g.tex_number(vt3["ANYCBCVT"].get(key))
        parts = []
        for p, pre in pipes:
            v = g.tex_number(vt3[pre + "VT"].get(key))
            parts.append(f"{LABEL[p]} {v / anyv:.2f}" if not math.isnan(v) and anyv else f"{LABEL[p]} -")
        out.append(f"  masses {bins.replace('ensuremath', '').strip('{}')}: " + ", ".join(parts) + f" (any = {anyv:g} Gpc^3 yr)")
    out.append("GWTC-3 'all O3' table:")
    for key in keys:
        bins = vt3["VTBINS"].get(key, key)
        anyv = g.tex_number(vt3["ANYCBCALLOTHREEVT"].get(key))
        parts = []
        for p, pre in pipes:
            v = g.tex_number(vt3.get(pre + "ALLOTHREEVT", {}).get(key))
            parts.append(f"{LABEL[p]} {v / anyv:.2f}" if not math.isnan(v) and anyv else f"{LABEL[p]} -")
        out.append(f"  masses {bins.replace('ensuremath', '').strip('{}')}: " + ", ".join(parts) + f" (any = {anyv:g})")
    vt21 = g.parse_macros(os.path.join(g.SRC, "2108.01045", "injection_vt_rate_macros.tex"))
    out.append("GWTC-2.1 (O3a) BBH populations: " + "; ".join(
        f"{k}: " + ", ".join(f"{LABEL[p]} {g.tex_number(vt21.get(pre + 'VT' + k, {}).get('', vt21.get(pre + 'VT' + k)) if False else '')}" for p, pre in pipes[1:])
        for k in ("BBHINJ", "BBHPOP")) if False else "GWTC-2.1 (O3a) BBH populations: see the macro dump below")
    for name in ("GSTLALVTBBHINJ", "MBTAVTBBHINJ", "PYCBCVTBBHINJ", "PYCBCBBHVTBBHINJ", "ANYCBCVTBBHINJ",
                 "GSTLALVTBBHPOP", "MBTAVTBBHPOP", "PYCBCVTBBHPOP", "PYCBCBBHVTBBHPOP", "ANYCBCVTBBHPOP"):
        out.append(f"  {name}: {macro_scalar(os.path.join(g.SRC, '2108.01045', 'injection_vt_rate_macros.tex'), name)}")
    return out


def macro_scalar(path: str, name: str) -> str:
    import re

    text = open(path, encoding="utf-8", errors="replace").read()
    m = re.search(r"\\(?:DeclareRobustCommand|newcommand)\{\\" + name + r"\}\{(.*)\}\s*$", text, re.M)
    return m.group(1) if m else "?"


# ----------------------------------------------------------------------------------------------------------------
# Injection-free relative catchability vs the injection-based relative VT, per pipeline and chirp-mass bin
# ----------------------------------------------------------------------------------------------------------------
VT_POINTS = {  # GWTC-3 'all O3' table, VT_i / VT_any at p_astro > 0.5 (printed by --vt); keyed by the injection chirp mass
    30.5: {"cwb": 0.49, "gstlal": 0.79, "mbta": 0.66, "pycbc": 0.62, "pycbcbbh": 0.82},
    22.9: {"cwb": 0.42, "gstlal": 0.77, "mbta": 0.61, "pycbc": 0.61, "pycbcbbh": 0.83},
    17.4: {"cwb": 0.33, "gstlal": 0.79, "mbta": 0.68, "pycbc": 0.67, "pycbcbbh": 0.84},
    12.2: {"cwb": 0.31, "gstlal": 0.80, "mbta": 0.71, "pycbc": 0.73, "pycbcbbh": 0.85},
    6.1: {"cwb": 0.12, "gstlal": 0.81, "mbta": 0.74, "pycbc": 0.87, "pycbcbbh": 0.84},
    4.35: {"cwb": 0.07, "gstlal": 0.82, "mbta": 0.62, "pycbc": 0.87, "pycbcbbh": 0.68},
}


def run_compare() -> list[str]:
    """For the O3 LVK lists: each pipeline's share of the union (weighted by each candidate's max p_astro, i.e. the
    expected real signals), per chirp-mass bin, with a binomial error, against the injection-based VT_i / VT_any of the
    GWTC-3 table at the nearest injection chirp mass. Both at p_astro > 0.5."""
    out = []
    df = load("O3")
    lists = RUNS["O3"][2]
    for crit in ("pastro", "far"):
        Y = histories(df, lists, crit)
        keep = Y.any(axis=1)
        d = df[keep].reset_index(drop=True)
        Y = Y[keep]
        pcols = [L + "_pastro" for L in lists if L + "_pastro" in d.columns]
        w = d[pcols].max(axis=1).to_numpy(float)
        w = np.where(np.isnan(w), 1.0, np.clip(w, 0, 1))
        mc = d["mchirp"].to_numpy(float)
        out.append(f"O3, five LVK lists, threshold {crit}: pipeline share of the union's expected real signals (sum of max p_astro) "
                   f"vs the injection-based VT_i / VT_any (GWTC-3 all-O3 table, p_astro > 0.5) at the nearest injection chirp mass")
        bins = [(0, 5, 4.35), (5, 15, 6.1), (15, 30, 22.9), (30, 60, 30.5), (60, 1000, 30.5), (5, 1000, None)]
        for lo, hi, key in bins:
            m = (mc >= lo) & (mc < hi)
            W = w[m].sum()
            if W < 1:
                continue
            parts = []
            for j, L in enumerate(lists):
                share = w[m & Y[:, j]].sum() / W
                err = max(math.sqrt(share * (1 - share) / W), 1.0 / (W + 2))  # a floor at one count for shares of 0 or 1
                if key is not None:
                    vt = VT_POINTS[key][L]
                    parts.append(f"{LABEL[L]} {share:.2f} ± {err:.2f} (VT {vt:.2f}, {(share - vt) / err:+.1f}σ)")
                else:
                    # population-weighted VT ratio: weight the injection points by the data's expected real signals per bin
                    parts.append(f"{LABEL[L]} {share:.2f} ± {err:.2f}")
            out.append(f"  Mc {lo:g}-{hi:g}: {int(m.sum())} candidates, {W:.1f} expected real: " + "; ".join(parts))
    return out


if __name__ == "__main__":
    if "--mock" in sys.argv:
        nb, reps = 100, 3
        for a in sys.argv[1:]:
            if a.startswith("--nboot="):
                nb = int(a.split("=")[1])
            if a.startswith("--reps="):
                reps = int(a.split("=")[1])
        lines = run_mock(nboot=nb, reps=reps)
    elif "--vt" in sys.argv:
        lines = run_vt()
    elif "--compare" in sys.argv:
        lines = run_compare()
    else:
        nb = 300
        for a in sys.argv[1:]:
            if a.startswith("--nboot="):
                nb = int(a.split("=")[1])
        lines = run_data(nboot=nb)
    print("\n".join(lines))