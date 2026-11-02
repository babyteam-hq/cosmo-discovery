"""R-156: multipoles, templates, jackknife covariance, A_D fits and nulls from the stored shell sums.

usage: python r156_measure.py <pairs-tag> [--nulls-only] [--dmin 20]
Reads data/r156/sample.npz and data/r156/pairs_<tag>.npz; writes data/r156/measure_<tag>.json and prints the table.
"""
import sys, json
import numpy as np
from scipy.special import eval_legendre
import r156_kappa as K
import r156_theory as T

D = 'data/r156/'
BINS = [(0.15, 0.25), (0.25, 0.35), (0.35, 0.45)]
NJK = 200


def mu_weights(nmu, ell):
    edges = np.linspace(-1, 1, nmu + 1)
    x = np.linspace(-1, 1, 20001)
    w = np.zeros(nmu)
    for i in range(nmu):
        m = (x >= edges[i]) & (x <= edges[i + 1])
        w[i] = np.trapezoid(eval_legendre(ell, x[m]), x[m])
    return (2 * ell + 1) / 2 * w


def jk_regions(ra, dec, sel, n=NJK):
    import healpy as hp
    pix = hp.ang2pix(16, ra, dec, lonlat=True, nest=True)
    up, cnt = np.unique(pix[sel], return_counts=True)
    cum = np.cumsum(cnt) / cnt.sum()
    grp = np.minimum((cum * n).astype(int), n - 1)
    lut = dict(zip(up, grp))
    return np.array([lut.get(p, -1) for p in pix])


class Shells:
    def __init__(self, tag):
        p = np.load(D + 'pairs_%s.npz' % tag)
        self.Q, self.R, self.kidx, self.zk = p['Q'], p['R'], p['kidx'], p['zk']
        de = p['dedges']
        self.dcen = 0.5 * (de[1:] + de[:-1])
        self.nmu = int(p['nmu'])
        self.L = {l: mu_weights(self.nmu, l) for l in range(4)}

    def sums(self, mark, u, sel, reg):
        """Per-region numerator and denominator arrays (NJK, nd, nmu)."""
        m = sel
        wq = (u * mark)[m]
        num = np.zeros((NJK,) + self.Q.shape[1:])
        den = np.zeros_like(num)
        r = reg[m]
        for g in range(NJK):
            mm = r == g
            if mm.any():
                num[g] = np.einsum('j,jdm->dm', wq[mm], self.Q[m][mm], optimize=True)
                den[g] = np.einsum('j,jdm->dm', u[m][mm], self.R[m][mm], optimize=True)
        return num, den

    def multipoles(self, num, den, ell):
        tot = (num.sum(0) / den.sum(0)) @ self.L[ell]
        jk = np.array([((num.sum(0) - num[g]) / (den.sum(0) - den[g])) @ self.L[ell] for g in range(NJK)])
        return tot, jk


def jk_cov(jk):
    n = len(jk)
    d = jk - jk.mean(0)
    return (n - 1) / n * d.T @ d


def gls(y, t, C, nd_used, nsamp=NJK):
    Ci = np.linalg.inv(C) * (nsamp - nd_used - 2) / (nsamp - 1)  # Hartlap
    var = 1 / (t @ Ci @ t)
    a = var * (t @ Ci @ y)
    chi2 = (y - a * t) @ Ci @ (y - a * t)
    chi0 = y @ Ci @ y
    return a, np.sqrt(var), chi2, chi0


def eta_of_z(s, k, u, sel, zgrid, mlo=19.3):
    """eta(z) = 1 + 2.1715 d<kappa_hat>_u/dm_lim, slope over m_cut in [mlo, 19.5], z slices of 0.02."""
    r = s['RMAG']
    out = []
    for z0 in zgrid:
        m = sel & (np.abs(s['Z'] - z0) < 0.01)
        cuts = np.linspace(mlo, 19.5, 5)
        mu = [np.average(k[m & (r < c)], weights=u[m & (r < c)]) for c in cuts]
        out.append(1 + 2.1715 * np.polyfit(cuts, mu, 1)[0])
    return np.array(out)


def prepare(with_mag=False):
    s = dict(np.load(D + 'sample.npz'))
    z = s['Z']
    sel = s['GOOD'] & (z > 0.15) & (z < 0.45)
    rel = K.fit_relation(s, sel, with_mag)
    k, u = K.kappa_hat(s, rel, with_mag)
    # zero weighted mean in z slices of 0.01 per cap
    for cap in (0, 1):
        for z0 in np.arange(0.15, 0.45, 0.01):
            m = sel & (s['CAP'] == cap) & (z >= z0) & (z < z0 + 0.01)
            if m.any():
                k[m] -= np.average(k[m], weights=u[m])
    return s, sel, k, u, rel


def standardise(x, s, sel):
    """Residual of x against a cubic in z (per cap), divided by its scatter: a unit-variance mark."""
    y = np.zeros(len(x))
    for cap in (0, 1):
        m = sel & (s['CAP'] == cap)
        c = np.polyfit(s['Z'][m] - 0.3, x[m], 3)
        y[m] = x[m] - np.polyval(c, s['Z'][m] - 0.3)
    y[sel] /= np.std(y[sel])
    return y


def eta_full(s, k, u, sel, zc, a_resp, clipfrac):
    """eta(z) through every cut that magnification moves (method-check point 2):
    a_resp * clipfrac (Tractor response x winsorisation) + flux-limit edge + sigma_lambda < 0.5 edge.
    The flux edge is 2.1715 d<kappa_hat>/dm_lim (data, m_cut 19.3-19.5); the sigma_lambda edge counts the galaxies a
    magnification kappa0 brings under 0.5 (sigma_lambda scales as 1/SNR ~ 1/(1+kappa) at fixed surface brightness)."""
    flux = eta_of_z(s, k, u, sel, zc) - 1
    # sigma_lambda edge: candidates are GOOD except for the sigma cut
    base = s['GOOD_NOSIG'] & (s['Z'] > 0.15) & (s['Z'] < 0.45)
    k0 = 0.02
    sl = []
    for z0 in zc:
        m = sel & (np.abs(s['Z'] - z0) < 0.01)
        enter = base & (~s['GOOD']) & (s['SIG_LAMBDA_RAW'] > 0.5) & (s['SIG_LAMBDA_RAW'] < 0.5 * (1 + k0)) \
            & (np.abs(s['Z'] - z0) < 0.01)
        # their kappa_hat: lambda relative to the same relation, weight as the sample's
        n_in = u[m].sum()
        add = (s['KHAT_ALL'][enter] * s['U_ALL'][enter]).sum()
        wadd = s['U_ALL'][enter].sum()
        mean0 = np.average(k[m], weights=u[m])
        sl.append(((add + mean0 * n_in) / (n_in + wadd) - mean0) / k0)
    return a_resp * clipfrac + flux + np.array(sl), flux, np.array(sl)


ZC = np.arange(0.16, 0.45, 0.02)


def bin_template(sh, sel_j, u_j, eta_z, b, zj, kind='dip'):
    """Pair-weighted average over the bin's size galaxies of eta(z_j) t(z_j, d) (t: dipole or octupole)."""
    zgrid = np.arange(0.15, 0.4501, 0.01)
    if kind == 'dip':
        tz = np.array([T.template(z, b, sh.dcen) for z in zgrid])
    else:
        tz = np.array([T.template_oct(z, sh.dcen) for z in zgrid])
    ez = np.interp(zgrid, ZC, eta_z)
    iz = np.clip(np.round((zj - 0.15) / 0.01).astype(int), 0, len(zgrid) - 1)
    Rj = sh.R[sel_j].sum(2)  # (n, nd)
    w = u_j[:, None] * Rj
    out = np.zeros(len(sh.dcen))
    for g in np.unique(iz):
        m = iz == g
        out += (w[m] * (ez[g] * tz[g])[None, :]).sum(0)
    return out / w.sum(0)
