"""R-143: statistics beside report.md: the transfer-stability rule (T_B/T_A, T_E/T_A per group), w̄_gs against zero, and
the change of the T_C-corrected signal when the fixed star mask is added (jackknife of the difference).
Usage: python -u spherex_r143_extra.py"""
import os, sys
import numpy as np
from scipy.stats import chi2
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from spherex_transfer import OUT, jkerr, NS, NG, group_edges
J = np.load(os.path.join(OUT, "result.npz"))
lamc = np.array([np.sqrt(a * b) for a, b in group_edges()])

def ivm_ratio(num, den, numj, denj):
    v, vj = num / den, numj / denj
    e = jkerr(vj); w = np.where(np.isfinite(e) & (e > 0), 1 / e ** 2, 0)
    return (np.nan_to_num(v) * w).sum(0) / w.sum(0), jkerr((np.nan_to_num(vj) * w).sum(1) / w.sum(0))

TA, TAe = ivm_ratio(J["recA"], J["truA"], J["recA_jk"], J["truA_jk"])
TB, TBe = ivm_ratio(J["recB"] / 5, J["truA"], J["recB_jk"] / 5, J["truA_jk"])
TE, TEe = ivm_ratio(J["recE"], J["truE"], J["recE_jk"], J["truE_jk"])
TC, _ = ivm_ratio(J["recC"], J["truC"], J["recC_jk"], J["truC_jk"])
TD, _ = ivm_ratio(J["recD"], J["truD"], J["recD_jk"], J["truD_jk"])
print("| group µm | T_B/T_A | T_E/T_A | stable (both within 15 %) | T_D/T_C |")
print("|---|---|---|---|---|")
for g in range(NG):
    st = abs(TB[g] / TA[g] - 1) <= 0.15 and abs(TE[g] / TA[g] - 1) <= 0.15
    print(f"| {lamc[g]:.3f} | {TB[g]/TA[g]:.2f} | {TE[g]/TA[g]:.2f} | {'yes' if st else 'no'} | {TD[g]/TC[g]:.2f} |")
w, wj = J["wgs"], J["wgs_jk"]
c = float(((w / jkerr(wj)) ** 2).sum())
print(f"\nw̄_gs against zero over 32 slices: χ² = {c:.1f}/32 (p = {chi2.sf(c, 32):.2f}); inverse-variance mean "
      f"{(w / jkerr(wj)**2).sum() / (1 / jkerr(wj)**2).sum():+.4f} ± {1/np.sqrt((1/jkerr(wj)**2).sum()):.4f}")
# signal with and without the star mask, T_C-corrected via the same T: difference of (real_sm - real)/(b wm T_C)
b, wmb = J["b"][:, None], J["wmb"][:, None]
bj, wmbj = J["b_jk"][:, :, None], J["wmb_jk"][:, :, None]
tC, tCj = J["recC"] / J["truC"], J["recC_jk"] / J["truC_jk"]
d = (J["real_sm"] - J["real"]) / (b * wmb * tC)
dj = (J["real_sm_jk"] - J["real_jk"]) / (bj * wmbj * tCj)
tonw = 2.99792458e14 / lamc * 1e-11
e = jkerr(dj)
print("\nStar mask minus R-138 masks, T_C-corrected signal (nW m^-2 sr^-1 per unit z), inverse-variance mean over slices, "
      "jackknife of the difference; χ² of the per-slice differences against zero:")
for g in range(NG):
    ok = np.isfinite(e[:, g]) & (e[:, g] > 0)
    wv = 1 / e[ok, g] ** 2
    mu = (d[ok, g] * wv).sum() / wv.sum() * tonw[g]
    se = jkerr((dj[:, ok, g] * wv).sum(1) / wv.sum()) * tonw[g]
    cc = float(((d[ok, g] / e[ok, g]) ** 2).sum())
    print(f"{lamc[g]:.3f}: {mu:+.2f} ± {se:.2f}; χ² {cc:.0f}/{ok.sum()} (p={chi2.sf(cc, ok.sum()):.2f})")
