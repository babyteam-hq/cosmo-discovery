"""cosmo-ga-rev (R-171): inference from the injections (`garev_inject.py` -> data/r155/garev_runs.npz).

GA-like recovery: 30-70 members, largest separation 0.5-2x the GA's 742 Mpc, delta_mem within +-0.3 (variant +-0.5) of the
GA rebuild's own value under the same completeness model. L(A_true) = share of sites giving one.

Count-level comparison with the DESI measurement (8 absorbers in the GA's CHMS against 8.86 expected, W > 0.3; 8 against
8.40 at W > 0.6): under a hypothesis the hull's expected count is mu = N_exp (1 + bg + r A_true f_rec), with
  bg     the chance part of a GA-like hull under that hypothesis: pooled DESI N_obs / N_exp - 1 in the recovered hulls at
         the injection sites (real DESI data; the injected absorbers are not in it);
  f_rec  each recovery's share of its hull's DESI expectation inside the injected region, rotated back onto the GA;
  r      delta_DESI / delta_SDSS for the same structure (1 by default; scanned).
N ~ Poisson(mu x lognormal clustering scatter); the scatter's width is what the GA's shift null shows beyond Poisson
(sd 0.341 at N_exp 8.86 against Poisson 0.336). Reported: P(N <= N_obs) under "real" (the posterior mixture over A_true > 0
with a flat prior on the grid, variant prior 1/A) and under "chance" (A_true = 0), and their likelihood ratio at N_obs.
Errors: bootstrap over injection sites.

Run from the repo root: python -u research/cosmo-discovery/garev_infer.py
"""
import sys

import numpy as np
from scipy.stats import norm, poisson

sys.path.insert(0, "research/cosmo-discovery")
from r155_common import DATA  # noqa: E402

GA_SEP = 742.0
OBS = {0.3: (8, 8.86), 0.6: (8, 8.40)}
SIG_C = np.sqrt(max(0.0, 0.341 ** 2 - 1 / 8.86))
_eps = np.random.default_rng(5).normal(-SIG_C ** 2 / 2, SIG_C, 400) if SIG_C > 0 else np.zeros(1)


def cdf_pmf(nobs, mu):
    m = np.maximum(np.asarray(mu, float)[..., None] * np.exp(_eps), 1e-9)
    return poisson.cdf(nobs, m).mean(-1), poisson.pmf(nobs, m).mean(-1)


def grid_weights(A, prior):
    a = np.unique(A[A > 0])
    edges = np.r_[a[0] - (a[1] - a[0]) / 2, (a[1:] + a[:-1]) / 2, a[-1] + (a[-1] - a[-2]) / 2]
    w = np.diff(edges)
    if prior == "1/A":
        w = w / a
    return dict(zip(a, w))


def analyse(R, run, dwin=0.3, prior="flat", r=1.0, cut=0.3, sites=None):
    P, model = run
    nsite = int(R["k"].max()) + 1
    sel = (R["P"] == P) & (R["model"] == model)
    if sites is not None:  # bootstrap: rows of the resampled sites, with multiplicity
        idx = np.concatenate([np.nonzero(sel & (R["k"] == s))[0] for s in sites])
    else:
        idx = np.nonzero(sel)[0]
    dGA = float(R[f"dGA_{model}"])
    like = ((R["size"][idx] >= 30) & (R["size"][idx] <= 70) & (R["sep"][idx] >= 0.5 * GA_SEP) &
            (R["sep"][idx] <= 2 * GA_SEP) & (np.abs(R["dmem"][idx] - dGA) <= dwin))
    A = R["A"][idx]
    ko, ke, kf = ("bg_obs", "bg_exp", "f_rec") if cut == 0.3 else ("bg6_obs", "bg6_exp", "f_rec6")
    nobs, nexp = OBS[cut]
    out = dict(L={float(a): float(like[A == a].sum() / (A == a).sum()) for a in np.unique(A)},
               n_like={float(a): int(like[A == a].sum()) for a in np.unique(A)})
    ch = like & (A == 0)
    re = like & (A > 0)
    if ch.sum():
        bg_c = R[ko][idx][ch].sum() / R[ke][idx][ch].sum() - 1
    else:
        bg_c = np.nan
    gw = grid_weights(R["A"], prior)
    w = np.array([gw.get(float(a), 0.0) for a in A[re]])
    bo, be = R[ko][idx][re], R[ke][idx][re]
    bg_r = (w * bo).sum() / (w * be).sum() - 1 if re.sum() else np.nan
    f = R[kf][idx][re]
    dpred = bg_r + r * A[re] * f
    pc_r, pm_r = cdf_pmf(nobs, nexp * (1 + dpred))
    P_real = (w * pc_r).sum() / w.sum() if re.sum() else np.nan
    pm_real = (w * pm_r).sum() / w.sum() if re.sum() else np.nan
    P_ch, pm_ch = cdf_pmf(nobs, nexp * (1 + bg_c)) if np.isfinite(bg_c) else (np.nan, np.nan)
    o = np.argsort(dpred)
    cw = np.cumsum(w[o]) / w.sum() if re.sum() else np.zeros(0)
    q = (lambda p: float(dpred[o][np.searchsorted(cw, p)])) if re.sum() else (lambda p: np.nan)
    o2 = np.argsort(A[re])
    cw2 = np.cumsum(w[o2]) / w.sum() if re.sum() else np.zeros(0)
    qa = (lambda p: float(A[re][o2][np.searchsorted(cw2, p)])) if re.sum() else (lambda p: np.nan)
    out.update(n_chance=int(ch.sum()), n_real=int(re.sum()), bg_chance=float(bg_c), bg_real=float(bg_r),
               f_rec_med=float(np.median(f)) if re.sum() else np.nan,
               A_true_post=(qa(0.16), qa(0.5), qa(0.84)), delta_pred=(q(0.16), q(0.5), q(0.84)),
               P_real=float(P_real), P_chance=float(P_ch), LR_real_vs_chance=float(pm_real / pm_ch),
               merged_med=float(np.median(R["merged"][idx][re])) if re.sum() else np.nan,
               ninj_frac=float(np.median(R["ninj"][idx][re] / R["size"][idx][re])) if re.sum() else np.nan)
    out["sigma_real"] = float(norm.isf(P_real)) if P_real > 0 else np.inf
    out["nsite"] = nsite
    return out


def neff(cen, deg=20.0):
    keep = []
    for i, c in enumerate(cen):
        if all(np.degrees(np.arccos(np.clip(c @ cen[j], -1, 1))) > deg for j in keep):
            keep.append(i)
    return len(keep)


def main():
    R = dict(np.load(f"{DATA}/garev_runs.npz", allow_pickle=True))
    cen = np.load(f"{DATA}/garev_site_centres.npy")
    grp = np.load(f"{DATA}/garev_groups.npy")
    print(f"rows {len(R['k'])}; delta_GA u {float(R['dGA_u']):.3f}, g {float(R['dGA_g']):.3f}; clustering scatter "
          f"beyond Poisson {SIG_C:.3f}; sites {int(R['k'].max()) + 1}, independent (> 20 deg apart, greedy) {neff(cen)}")
    # chance check on the real slice's own groups (the GA excluded)
    dGA = float(R["dGA_u"])
    g_like = ((grp[:, 1] >= 30) & (grp[:, 1] <= 70) & (grp[:, 2] >= 0.5 * GA_SEP) & (grp[:, 2] <= 2 * GA_SEP) &
              np.isfinite(grp[:, 4]) & (grp[:, 5] > 0))
    notga = ~((grp[:, 1] == 46) & (np.abs(grp[:, 3] - dGA) < 1e-6))
    for win in (0.3, 0.5, 9):
        m = g_like & notga & (np.abs(grp[:, 3] - dGA) <= win)
        if m.sum():
            print(f"slice groups, GA-like within {win}: {m.sum()}, DESI pooled delta {grp[m, 4].sum() / grp[m, 5].sum() - 1:+.3f} "
                  f"({grp[m, 4].sum():.0f} / {grp[m, 5].sum():.1f})")
    rng = np.random.default_rng(171)
    # block bootstrap: sites grouped by the nearest of the greedy independent centres; blocks resampled
    keep = []
    for i, c in enumerate(cen):
        if all(np.degrees(np.arccos(np.clip(c @ cen[j], -1, 1))) > 20 for j in keep):
            keep.append(i)
    blk = np.argmax(cen @ cen[keep].T, 1)
    blocks = [np.nonzero(blk == b)[0] for b in range(len(keep))]
    boots = [np.concatenate([blocks[j] for j in rng.integers(len(blocks), size=len(blocks))]) for _ in range(100)]
    for run in (("H", "u"), ("S", "u"), ("L", "u"), ("H", "g")):
        for kw in (dict(), dict(dwin=0.5), dict(prior="1/A"), dict(cut=0.6)):
            o = analyse(R, run, **kw)
            bs = np.array([analyse(R, run, sites=b, **kw)["P_real"] for b in boots[:100]])
            bs = bs[np.isfinite(bs)]
            tag = f"{run[0]}/{run[1]} " + (" ".join(f"{k}={v}" for k, v in kw.items()) or "primary")
            print(f"\n== {tag}")
            print("  L(A):", {a: round(v, 3) for a, v in o["L"].items()}, " n:", o["n_like"])
            print(f"  GA-like: chance {o['n_chance']}, real {o['n_real']}; bg chance {o['bg_chance']:+.3f}, "
                  f"bg real {o['bg_real']:+.3f}; f_rec median {o['f_rec_med']:.2f}; injected share of members "
                  f"{o['ninj_frac']:.2f}; merged-from-real median {o['merged_med']:.0f}")
            print(f"  A_true posterior 16/50/84: {o['A_true_post']}; predicted delta in the hull 16/50/84: "
                  f"{tuple(round(x, 3) for x in o['delta_pred'])} (observed {OBS[kw.get('cut', 0.3)][0] / OBS[kw.get('cut', 0.3)][1] - 1:+.3f})")
            print(f"  P(N <= obs | real) {o['P_real']:.4f} ({o['sigma_real']:.2f} sigma), bootstrap 16-84 % "
                  f"{np.percentile(bs, 16):.4f}-{np.percentile(bs, 84):.4f}; P(N <= obs | chance) {o['P_chance']:.3f}; "
                  f"likelihood ratio real/chance at N_obs {o['LR_real_vs_chance']:.3f}")
        # relative bias scan
        for r in (0.25, 0.5, 0.75, 1.0, 1.25):
            o = analyse(R, run, r=r)
            print(f"  r = {r}: P_real {o['P_real']:.4f} ({o['sigma_real']:.2f} sigma), LR {o['LR_real_vs_chance']:.3f}")


if __name__ == "__main__":
    main()
