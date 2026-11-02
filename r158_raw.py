"""R-161 (cosmo-deep-25b): the R-158 Euclid Q1 result restated raw (no PSF or dilution correction), with the worst-case correction
bound from the data. Works only from R-158's landed per-field raw numbers (results-R158.md, R1 and R2); nothing is
refetched. Dipoles: the same per-field, Euclid-only, no-intercept fit as r158_stats.py, on raw A and raw sigma.
Usage: python3 r161_raw.py"""
import numpy as np

AX = {"G (SGP)": (12.86, -27.13), "H (HST best fit)": (78.0, 47.0), "S (SDSS z>0.15)": (71.0, 61.0)}
CEN = {"EDF-N": (269.73, 66.02), "EDF-F": (52.93, -28.09), "EDF-S": (61.24, -48.42), "GOODS-S": (53.12, -27.81)}
# R-158 R1, CNN, 0.3 < z < 1.2: raw A and its sigma (larger of binomial and tile bootstrap), per field, in %
RAW = {"EDF-N": (0.16, 0.33), "EDF-F": (0.01, 0.45), "EDF-S": (0.80, 0.28)}
ALL = (0.44, 0.19)                      # three fields together, 291,652 classified
Z3 = 3.0                                # 99.7 % two-sided, as R-158's verdicts

def cosang(a, b):
    a, b = np.radians(a), np.radians(b)
    return np.sin(a[1]) * np.sin(b[1]) + np.cos(a[1]) * np.cos(b[1]) * np.cos(a[0] - b[0])

def dipole(ax):
    F = list(RAW); c = np.array([cosang(CEN[f], ax) for f in F])
    y = np.array([RAW[f][0] for f in F]); s = np.array([RAW[f][1] for f in F])
    w = c ** 2 / s ** 2; Ad = (c * y / s ** 2).sum() / w.sum(); sd = 1 / np.sqrt(w.sum())
    proj = (c / s ** 2).sum() / w.sum()  # a uniform bias b enters the no-intercept fit as proj * b
    return Ad, sd, dict(zip(F, c.round(3))), proj

# Data bounds on the additive (chiral) bias b. The pipeline is exactly mirror-odd (classifier, prep, selection), so a
# net sign can only come from the images; the instrument's chiral part is the PSF. Null 4b re-convolved 3,000 real
# Euclid discs with the VIS PSF and with its mirror: dA = -0.10 +- 0.08 %, PSF term dA/2.
dA, sdA = -0.10, 0.08
b4b = abs(dA / 2) + Z3 * sdA / 2         # 99.7 % bound on |PSF term| from real galaxies
B = 3 * b4b                              # allowance: the re-convolution is not the native PSF (factor 3)

# Real-galaxy misclassification: our CNN against the published labels on JADES F200W (181 of 182 agree)
k, n = 1, 182
def cp_upper(k, n, alpha=0.00135):
    """One-sided Clopper-Pearson upper limit (alpha = the 3 sigma tail), by bisection on the binomial CDF."""
    from math import comb
    cdf = lambda p: sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k + 1))
    lo, hi = k / n, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if cdf(mid) > alpha else (lo, mid)
    return lo
eps_jades_up = cp_upper(k, n)

def need_eps(Araw, s, claim, b=B):
    """Smallest misclassification rate that lets a true asymmetry `claim` hide behind the measured raw A:
    A_raw = (1 - 2 eps) A_true + b, with the most hiding b = -B, at the raw 99.7 % upper limit."""
    d = (Araw + Z3 * s + b) / claim
    return d, (1 - d) / 2

def row(name, Araw, s, claims, proj=1.0):
    out = []
    for cname, cv in claims:
        d, e = need_eps(Araw, s, cv, B * abs(proj))
        out.append(f"| {name} | {cname} = {cv:g} % | {Araw:+.2f} +- {s:.2f} % | {Araw + Z3 * s:+.2f} % | {d:.3f} | {e:.2f} |")
    return out

if __name__ == "__main__":
    print(f"bias bound: null 4b PSF term {dA/2:+.2f} +- {sdA/2:.2f} %, 99.7 % |b| <= {b4b:.2f} %, allowance B = {B:.2f} %")
    print(f"JADES real-galaxy eps (1/182 disagree): 99.7 % upper {eps_jades_up:.3f}")
    print(f"cos(GOODS-S, SGP) = {cosang(CEN['GOODS-S'], AX['G (SGP)']):.3f}")
    D = {}
    for a, ax in AX.items():
        Ad, sd, c, p = dipole(ax); D[a] = (Ad, sd, p)
        print(f"raw dipole axis {a}: A_d = {Ad:+.2f} +- {sd:.2f} %  cos {c}  uniform-bias projection {p:+.2f}")
    print("\n| test | published | ours, raw | raw 99.7 % upper | largest 1 - 2 eps that hides it (bias at the bound) | eps needed |")
    print("|---|---|---|---|---|---|")
    L = row("EDF-F A", *RAW["EDF-F"], [("P2, all z", 20.2), ("P2 at 0.3<z<1.2", 9.3)])
    L += row("Euclid, three fields, A", *ALL, [("P2 at 0.3<z<1.2", 9.3)])
    L += row("dipole, axis G", *D["G (SGP)"][:2], [("implied", 24.8), ("at matched z", 11.4)], D["G (SGP)"][2])
    L += row("dipole, axis H", *D["H (HST best fit)"][:2], [("implied", 10.5), ("proposal", 5.0)], D["H (HST best fit)"][2])
    print("\n".join(L))
    # the most confident 5 % and the largest quartile (R2, post-verdict diagnostics), Euclid all fields
    for nm, (a, s) in {"most confident 5 % (|f| > 12.9)": (0.56, 0.64), "largest R50/FWHM quartile": (0.10, 0.34)}.items():
        d, e = need_eps(a, s, 9.3)
        print(f"{nm}: {a:+.2f} +- {s:.2f} %: hides 9.3 % only if 1-2eps <= {d:.3f} (eps >= {e:.2f})")
