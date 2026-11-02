"""R-156: refit after the pre-registered covariance failed (results F4).

The pre-registered C (size shuffles + analytic cosmic variance + a rank-1 mock-leak term) gives chi2/dof 4.8-6.7 and
the leak it subtracts is the mean of 3 mocks whose per-mock values scatter about zero (mock noise, not a leak). Here:
the 200-region jackknife covariance (Hartlap), which the mocks validate (jackknife error vs mock-to-mock scatter),
T_rel subtracted with its error, no mock-mean leak subtraction; the b_s leak enters as a systematic, the analytic
O(d/r) wide-angle estimate |leak| <= f b_s (d/r) |xi_2-like| (about 0.1 in A_D units). Everything else as in run.json.

usage: python -u r156_refit.py
"""
import json
import numpy as np
from scipy.stats import chi2 as C2
import r156_measure as M
import r156_run as RU
import r156_theory as T

D = 'data/r156/'


def main():
    run = json.load(open(D + 'run.json'))
    s, sel, k, u, rel = M.prepare(False)
    sh = M.Shells('data')
    kid = sh.kidx
    gj = sel[kid]
    zj = s['Z'][kid]
    reg = M.jk_regions(s['RA'][kid], s['DEC'][kid], np.ones(len(kid), bool))
    order = np.argsort(reg, kind='stable')
    bounds = np.searchsorted(reg[order], np.arange(M.NJK + 1))
    kj = np.where(gj, k[kid], 0.0)
    uj = np.where(gj, u[kid], 0.0)
    out = {}
    combo = []
    for bi, (lo, hi) in enumerate(M.BINS):
        R = run['bins'][str(bi)]
        B = gj & (zj > lo) & (zj < hi)
        w = np.where(B, uj, 0.0)
        num, den = RU.region_sums(sh, order, bounds, w * kj, w)
        xi, jk = sh.multipoles(num, den, 1)
        tD = np.array(R['tD'])
        tRel = np.array(R['tRel'])
        relerr = np.abs(tRel) * np.hypot(RU.C_FRAC, R['b_s_err'] / max(abs(R['b_s']), 1e-9))
        r, Hc, f, Dz = T.background(R['zeff'])
        res = {}
        for name, sl, nd in (('20-150', RU.FIT, 13), ('30-150', RU.FIT30, 12), ('20-60', slice(1, 5), 4)):
            C = M.jk_cov(jk)[sl, sl] + np.outer(relerr[sl], relerr[sl])
            y = (xi - tRel)[sl]
            A, Ae, chi2, chi0 = M.gls(y, tD[sl], C, nd)
            res[name] = {'A_D': A, 'err': Ae, 'chi2_fit': chi2, 'dof': nd - 1, 'p_fit': 1 - C2.cdf(chi2, nd - 1),
                         'chi2_zero': chi0, 'p_zero': 1 - C2.cdf(chi0, nd)}
        # raw dipole against zero in eta = 1 template units (no eta needed)
        tD1 = tD / max(abs(R['eta']), 1e-9) * np.sign(R['eta']) if abs(R['eta']) > 0.05 else None
        A = res['20-150']['A_D']
        eta_err = abs(A) * np.hypot(R['eta_err_stat'], R['eta_err_sys']) / max(abs(R['eta']), 1e-9)
        leak_sys = 0.1 * abs(A) + 0.1  # analytic wide-angle estimate, in A_D units (results F4)
        bfid = 1.34 / Dz
        b_sys = abs(A) * abs(1 - (R['b'] + 0.6 * f) / (bfid + 0.6 * f)) / 2  # b from 20-40 vs DESI fiducial 1.34/D
        res['A_D_err_b'] = b_sys
        res['A_D_err_eta'] = eta_err
        res['A_D_err_leak'] = leak_sys
        res['xi1'] = xi.tolist()
        res['xi1_err_jk'] = np.sqrt(np.diag(M.jk_cov(jk))).tolist()
        out[bi] = res
        if abs(R['eta']) > 0.05:
            tot = np.sqrt(res['20-150']['err'] ** 2 + eta_err ** 2 + leak_sys ** 2 + b_sys ** 2)
            combo.append((A, tot, res['20-150']['err']))
        print('bin %.2f-%.2f eta %.3f: A_D = %.2f +- %.2f (stat, jk) +- %.2f (eta) +- %.2f (leak) +- %.2f (b) | chi2 fit %.1f/%d '
              '(p %.3f) | vs zero chi2 %.1f/13 (p %.2g) | d>=30: %.2f +- %.2f | 20-60: %.2f +- %.2f' % (
                  lo, hi, R['eta'], A, res['20-150']['err'], eta_err, leak_sys, b_sys, res['20-150']['chi2_fit'],
                  res['20-150']['dof'], res['20-150']['p_fit'], res['20-150']['chi2_zero'], res['20-150']['p_zero'],
                  res['30-150']['A_D'], res['30-150']['err'], res['20-60']['A_D'], res['20-60']['err']), flush=True)
    a = np.array([c[0] for c in combo]); e = np.array([c[1] for c in combo]); es = np.array([c[2] for c in combo])
    wgt = 1 / e ** 2
    Ac = (a * wgt).sum() / wgt.sum()
    out['combined'] = {'A_D': Ac, 'err': 1 / np.sqrt(wgt.sum()),
                       'err_stat_only': 1 / np.sqrt((1 / es ** 2).sum()), 'bins_used': len(combo)}
    print('combined (bins with |eta| > 0.05): A_D = %.2f +- %.2f (total); %.1f sigma from 0; %.1f sigma from 1' % (
        Ac, 1 / np.sqrt(wgt.sum()), Ac * np.sqrt(wgt.sum()), (Ac - 1) * np.sqrt(wgt.sum())))
    json.dump(out, open(D + 'refit.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
