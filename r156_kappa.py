"""R-156: standardised size residual kappa_hat per galaxy (the convergence estimator) and its response.

kappa_hat = lambda - <lambda>(z, colour, type)            ("zc", the pre-registered default; no magnitude term)
kappa_hat = lambda - <lambda>(z, colour, type, M_r)       ("zcm", the alternative)

<lambda> is a least-squares fit, per Tractor type (REX, EXP, DEV, SER), of a cubic polynomial in z, quadratics in
g-r and r-z and their cross terms with z (and, for zcm, a quadratic in the apparent r magnitude, which at fixed z is
the absolute magnitude up to the k-correction captured by the colour terms), iterated with 4-sigma clipping.
Residuals are clipped (winsorised) at +-3 sigma_tot, symmetric, so the estimator is not driven by outliers.
Weights: w = 1 / (sigma_int^2 + sigma_lambda^2), sigma_int the type's clipped residual scatter.

The response eta = d<kappa_hat>/d kappa (a true convergence kappa scales R by 1+kappa and flux by 1+2kappa, and the
sample is cut at r < 19.5) is measured on the data by `response()`: shift lambda by kappa and r by -2.5log10(1+2kappa)
for every galaxy, drop those that cross r = 19.5 using the parent counts just above the limit (the decrement is
the loss through the bright edge of the observed distribution, d ln N / d r at 19.5), refit nothing (the relation is
fixed), and difference the weighted mean kappa_hat.
"""
import numpy as np

TYPES = ['REX', 'EXP', 'DEV', 'SER']


def design(z, gr, rz, rmag, with_mag):
    zc = (z - 0.3) / 0.1
    cols = [np.ones_like(z), zc, zc ** 2, zc ** 3, gr, gr ** 2, rz, rz ** 2, gr * rz, gr * zc, rz * zc, gr * zc ** 2]
    if with_mag:
        m = rmag - 18.5
        cols += [m, m ** 2, m * zc, m * gr]
    return np.vstack(cols).T


def fit_relation(s, sel, with_mag=False):
    """Return per-type coefficients and intrinsic scatter, fitted on galaxies in `sel`."""
    z, gr, rz = s['Z'], s['GMAG'] - s['RMAG'], s['RMAG'] - s['ZMAG']
    out = {}
    for t in TYPES:
        m = sel & (s['TYPE'] == t)
        X = design(z[m], gr[m], rz[m], s['RMAG'][m], with_mag)
        y = s['LAMBDA'][m]
        keep = np.ones(m.sum(), bool)
        for _ in range(4):
            c, *_ = np.linalg.lstsq(X[keep], y[keep], rcond=None)
            res = y - X @ c
            sig = 1.4826 * np.median(np.abs(res[keep] - np.median(res[keep])))
            keep = np.abs(res) < 4 * sig
        out[t] = (c, np.std(res[keep]), sig)
    return out


def kappa_hat(s, rel, with_mag=False, lam=None, rmag=None):
    z, gr, rz = s['Z'], s['GMAG'] - s['RMAG'], s['RMAG'] - s['ZMAG']
    lam = s['LAMBDA'] if lam is None else lam
    rmag = s['RMAG'] if rmag is None else rmag
    k = np.zeros(len(z))
    w = np.zeros(len(z))
    for t in TYPES:
        m = s['GOOD'] & (s['TYPE'] == t)
        c, sd, sig = rel[t]
        X = design(z[m], gr[m], rz[m], rmag[m], with_mag)
        res = lam[m] - X @ c
        st = np.sqrt(sd ** 2 + s['SIG_LAMBDA'][m] ** 2)
        k[m] = np.clip(res, -3 * st, 3 * st)
        w[m] = 1.0 / st ** 2
    return k, w
