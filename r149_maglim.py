"""R-149: the magnitude model that puts a BGS / KiDS-bright selection on the UniverseMachine mock.

At each virtual redshift z_v, BGS galaxies with |z - z_v| < 0.02 give m_r - DM(z) = a_T + b_T (log M* - 10.5) + N(0, s_T)
for each type T (DESI QUIESCENT), fit by maximum likelihood with the truncation at r = 19.5 (only galaxies brighter than the
limit are in BGS). Writes data/r149/maglim.json: {z_v: {"Q": [a, b, s], "SF": [a, b, s]}}.

Usage: python -I r149_maglim.py <data dir>
"""
import json, os, sys
import numpy as np, fitsio
from scipy.optimize import minimize
from scipy.stats import norm
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r147_leak import dcom

ZV = (0.10, 0.15, 0.20, 0.25, 0.28, 0.35)


def dmod(z):
    return 5 * np.log10(dcom(z) * (1 + z) * 1e6 / 10)


def fit(x, y, lim):
    """Truncated linear regression y = a + b x + N(0, s), observed only where y < lim (lim per object)."""
    def nll(p):
        a, b, ls = p; s = np.exp(ls); mu = a + b * x
        return -np.sum(norm.logpdf(y, mu, s) - norm.logcdf(lim, mu, s))
    p0 = np.polyfit(x, y, 1)[::-1]
    r = minimize(nll, [p0[0], p0[1], np.log(0.4)], method='Nelder-Mead', options=dict(maxiter=4000, xatol=1e-5, fatol=1e-3))
    return [float(r.x[0]), float(r.x[1]), float(np.exp(r.x[2]))]


def main():
    d = sys.argv[1]
    g = fitsio.read(os.path.join(d, 'bgs-groups-dr1-v1.0.fits'), columns=['Z', 'LOGMSTAR', 'APP_MAG_R', 'QUIESCENT',
                                                                           'Z_ASSIGNED_FLAG'])
    ok = np.isfinite(g['LOGMSTAR']) & ((g['Z_ASSIGNED_FLAG'] == 0) | (g['Z_ASSIGNED_FLAG'] == -1)) & (g['APP_MAG_R'] < 19.5)
    g = g[ok]; lm = g['LOGMSTAR'] + 2 * np.log10(1 / 0.7); q = g['QUIESCENT'].astype(bool)
    out = {}
    for zv in ZV:
        zc = min(zv, 0.30)                  # beyond the BGS lens range, extrapolate the z = 0.30 fit with DM(z_v)
        s = (np.abs(g['Z'] - zc) < 0.02) & (lm > 8.0) & (lm < 12.0)
        y = g['APP_MAG_R'][s] - dmod(g['Z'][s]); lim = 19.5 - dmod(g['Z'][s]); x = lm[s] - 10.5
        out[f'{zv:.2f}'] = {}
        for tn, t in (('Q', q[s]), ('SF', ~q[s])):
            p = fit(x[t], y[t], lim[t]); out[f'{zv:.2f}'][tn] = p
            print(f'z_v {zv:.2f} (fit at {zc:.2f}) {tn}: N = {t.sum()}, m - DM = {p[0]:.3f} + {p[1]:.3f} (log M* - 10.5), '
                  f'scatter {p[2]:.3f}', flush=True)
    json.dump(out, open(os.path.join(d, 'maglim.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
