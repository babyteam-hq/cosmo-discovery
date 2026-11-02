"""R-156: linear-theory Doppler-magnification dipole template and its Gaussian covariance.

Template (Bonvin et al. 2017, arXiv:1610.05946, eq. 28, distant observer):
    xi_1(d) = f (H/c) (1 - c/(H r)) (b + 3f/5) nu_1(d),  nu_1(d) = (1/2pi^2) Int dk k P(k) j1(kd)
with H the conformal Hubble rate aH(z) in h/Mpc (c = 1), r the comoving distance. xi_1 is the coefficient of
P1(cos beta) in <Delta(x) kappa(x+d)>, beta the angle of d to the line of sight; it is negative (a galaxy behind an
overdensity is demagnified).

Covariance of the estimator xi_hat_1(d) = 3 <kappa_j D_j(d)> (D_j the shell-averaged mu-weighted contrast around the
kappa galaxy j), derived here for per-galaxy size noise sigma_k (uncorrelated between galaxies):
    C_noise(d,d') = 9 sigma_k^2 / N_k [ (b^2/3 + 2bf/5 + f^2/7) S(d,d') + delta_dd' / (3 nbar V_shell(d)) ]
    S(d,d') = (1/2pi^2) Int dk k^2 P(k) j1bar(kd) j1bar(kd')
The cosmic-variance term is eq. 13 of 2511.17858 (Andrianomena et al. 2019). Note: 2511.17858 eq. 14 normalises
the size-noise term by the number of grid cells V/l_p^3 instead of the number of galaxies N_k; at BGS densities
(nbar l_p^3 ~ 0.03-1) that makes its forecast optimistic by sqrt(1/(nbar l_p^3)) in S/N.
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import spherical_jn

C_KMS = 299792.458
_cache = {}


def cosmo():
    if 'c' in _cache:
        return _cache['c']
    import camb
    p = camb.CAMBparams()
    p.set_cosmology(H0=67.66, ombh2=0.02242, omch2=0.11933, mnu=0.06)
    p.InitPower.set_params(As=2.105e-9, ns=0.9665)
    p.set_matter_power(redshifts=[0.0], kmax=20.0)
    p.NonLinear = camb.model.NonLinear_none
    res = camb.get_results(p)
    kh, _, pk = res.get_matter_power_spectrum(minkh=1e-4, maxkh=10, npoints=1200)
    _cache['c'] = (res, kh, pk[0], p.omegam if hasattr(p, 'omegam') else (0.02242 + 0.11933 + 0.0006) / 0.6766 ** 2)
    return _cache['c']


def background(z):
    """Return r [Mpc/h], conformal H/c [h/Mpc], f, D (D(0)=1)."""
    res, kh, pk, om = cosmo()
    h = 0.6766
    r = res.comoving_radial_distance(z) * h
    Hz = res.hubble_parameter(z)  # km/s/Mpc
    Hc = Hz / (1 + z) / C_KMS / h  # conformal, h/Mpc
    omz = om * (1 + z) ** 3 / (Hz / 67.66) ** 2
    f = omz ** 0.55
    a = np.linspace(1e-3, 1, 4000)
    zz = 1 / a - 1
    E = res.hubble_parameter(zz) / 67.66
    # D(a) proportional to H(a) Int da / (a H)^3
    integ = np.concatenate([[0], np.cumsum(0.5 * (1 / (a[1:] * E[1:]) ** 3 + 1 / (a[:-1] * E[:-1]) ** 3) * np.diff(a))])
    Da = E * integ
    Da /= Da[-1]
    D = np.interp(1 / (1 + z), a, Da)
    return r, Hc, f, D


def j1bar(k, d, dd):
    """j1 averaged over a separation bin [d-dd/2, d+dd/2] with weight d^2."""
    x = np.linspace(d - dd / 2, d + dd / 2, 9)
    w = x ** 2
    return (spherical_jn(1, np.outer(k, x)) * w).sum(1) / w.sum()


def template(z, b, dcen, dd=10.0, f=None, prefactor=False):
    """xi_1(d) at redshift z for bias b (GR, A_D = 1). Returns array over dcen."""
    res, kh, pk0, om = cosmo()
    r, Hc, f0, D = background(z)
    f = f0 if f is None else f
    k = kh
    P = pk0 * D ** 2
    damp = np.exp(-(k / 2.0) ** 2)
    nu1 = np.array([np.trapezoid(k * P * j1bar(k, d, dd) * damp, k) for d in dcen]) / (2 * np.pi ** 2)
    pre = f * Hc * (1 - 1 / (Hc * r))
    if prefactor:
        return pre, nu1
    return pre * (b + 3 * f / 5) * nu1


def covariance(z, b, dcen, Nk, nbar, V, sigk, dd=10.0, f=None):
    res, kh, pk0, om = cosmo()
    r, Hc, f0, D = background(z)
    f = f0 if f is None else f
    k = kh
    P = pk0 * D ** 2
    damp = np.exp(-(k / 2.0) ** 2)
    J = np.array([j1bar(k, d, dd) for d in dcen])  # (nd, nk)
    S = np.einsum('ik,jk,k->ij', J, J, k ** 2 * P * damp * np.gradient(k)) / (2 * np.pi ** 2)
    Vsh = 4 * np.pi * dcen ** 2 * dd
    Cn = 9 * sigk ** 2 / Nk * ((b ** 2 / 3 + 2 * b * f / 5 + f ** 2 / 7) * S + np.diag(1 / (3 * nbar * Vsh)))
    S2 = np.einsum('ik,jk,k->ij', J, J, P ** 2 * damp * np.gradient(k)) / np.pi ** 2
    Ccv = 9 / V * (1 - 1 / (Hc * r)) ** 2 * (b ** 2 / 5 + 2 * b * f / 7 + f ** 2 / 9) * f ** 2 * Hc ** 2 * S2
    return Cn, Ccv


def snr(z, b, dcen, Nk, nbar, V, sigk, eta=1.0, dd=10.0):
    t = eta * template(z, b, dcen, dd)
    Cn, Ccv = covariance(z, b, dcen, Nk, nbar, V, sigk, dd)
    C = Cn + Ccv * eta ** 2
    return np.sqrt(t @ np.linalg.solve(C, t)), np.sqrt(np.sum(t ** 2 / np.diag(C))), t, C


def j3bar(k, d, dd):
    x = np.linspace(d - dd / 2, d + dd / 2, 9)
    w = x ** 2
    return (spherical_jn(3, np.outer(k, x)) * w).sum(1) / w.sum()


def template_oct(z, dcen, dd=10.0):
    """Octupole of <Delta kappa_v> (Bonvin et al. 2017 eq. 29): xi_3 = -f (H/c)(1 - c/(Hr)) (2f/5) nu_3(d),
    nu_3 = (1/2pi^2) Int dk k P(k) j3(kd). Independent of b; positive at z < 1."""
    res, kh, pk0, om = cosmo()
    r, Hc, f, D = background(z)
    P = pk0 * D ** 2
    damp = np.exp(-(kh / 2.0) ** 2)
    nu3 = np.array([np.trapezoid(kh * P * j3bar(kh, d, dd) * damp, kh) for d in dcen]) / (2 * np.pi ** 2)
    return -f * Hc * (1 - 1 / (Hc * r)) * (2 * f / 5) * nu3


def xi_matter_mono(z, dcen, dd=10.0):
    res, kh, pk0, om = cosmo()
    r, Hc, f, D = background(z)
    P = pk0 * D ** 2
    damp = np.exp(-(kh / 2.0) ** 2)
    out = []
    for d in dcen:
        x = np.linspace(d - dd / 2, d + dd / 2, 9)
        w = x ** 2
        j0 = (spherical_jn(0, np.outer(kh, x)) * w).sum(1) / w.sum()
        out.append(np.trapezoid(kh ** 2 * P * j0 * damp, kh) / (2 * np.pi ** 2))
    return np.array(out)


def c_rel(z, s, dlnn_dr):
    """Coefficient C(z) of v/c in the number counts' Doppler term (derivation in r156_plan.md, 'Model'):
    C = -Hdot/H^2 - 2/(rH) - 10 s - (dln nbar/dr)/H, after the luminosity-function evolution E is eliminated with
    the measured nbar(z) and s (E = dln nbar/dr + 5s(1/r + H))."""
    r, Hc, f, D = background(z)
    res = cosmo()[0]
    Hz = res.hubble_parameter(z)
    om = cosmo()[3]
    omz = om * (1 + z) ** 3 / (Hz / 67.66) ** 2
    hdot = 1 - 1.5 * omz
    return -hdot - 2 / (r * Hc) - 10 * s - dlnn_dr / Hc
