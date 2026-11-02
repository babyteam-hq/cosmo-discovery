"""R-148 (cosmo-deep-18) step 0: the model for the Lya forest x Compton-y cross-correlation.

  python -u lya_y_model.py            # writes data/lya/model.npz and prints <bP_e>(z) and xi(theta) per z bin

Model (stated before any cross-correlation, results-R148.md F1):
  xi_Fy(theta; z) = b_F(z) * (sigma_T / m_e c^2) * a(z) * <bP_e>(z) * (1/chi^2) * int l dl/2pi J0(l theta) B_l P_m(l/chi, z)
for a forest pixel at comoving distance chi(z). Limber; only k_par ~ 0 enters (y is a projection), so the forest's
Kaiser factor (1 + beta mu^2) is 1. <bP_e> is the proper (physical) bias-weighted mean electron pressure.
  - <bP_e>(z): halo model, Battaglia et al. 2012 (ApJ 758, 75) AGN-feedback pressure profile (P_th; P_e = 0.518 P_th
    for X = 0.76), Tinker 2008 M200m mass function and Tinker 2010 bias (colossus), M200c from 1e11 to 1e16 M_sun,
    profile integrated to 4 R200c. The lower limit 1e10 is printed as the bracket.
  - b_F(z): DESI DR1 Lya auto-correlation, b_F(2.33) = -0.1078 (arXiv:2404.03001), scaled as (1+z)^2.9.
  - P_m: CAMB halofit, Planck 2018 cosmology.
  - B_l: the y map's beam (ACT DR6 NILC ilc_beam.txt; Planck NILC 10' Gaussian).
"""
import numpy as np, sys, os
from scipy.special import j0
from scipy.integrate import simpson
import camb
from colossus.cosmology import cosmology
from colossus.lss import mass_function, bias
from colossus.halo import mass_defs

HERE = os.path.dirname(os.path.abspath(__file__))
cosmo = cosmology.setCosmology('planck18')
h = cosmo.h
SIGT = 6.6524587e-29          # m^2
MEC2_KEV = 510.99895          # keV
MPC_M = 3.0856775814913673e22
KEV_J = 1.602176634e-16
G_SI = 6.674e-11
MSUN = 1.98847e30

ZBINS = [(2.0, 2.5), (2.5, 3.0), (3.0, 3.6)]
THETA_EDGES = np.array([1, 2, 3, 5, 7.5, 10, 15, 20, 30, 45, 60.])  # arcmin
THETA_MID = np.sqrt(THETA_EDGES[1:] * THETA_EDGES[:-1])


def b_F(z):
    return -0.1078 * ((1 + z) / 3.33) ** 2.9


def battaglia_Pth(r_over_R, M200c, z):
    """Battaglia+2012 AGN profile, P_th in keV/cm^3 (proper); M200c in M_sun."""
    m = M200c / 1e14
    P0 = 18.1 * m ** 0.154 * (1 + z) ** -0.758
    xc = 0.497 * m ** -0.00865 * (1 + z) ** 0.731
    beta = 4.35 * m ** 0.0393 * (1 + z) ** 0.415
    x = r_over_R / xc
    rho_c = cosmo.rho_c(z) * 1e9 * h ** 2 * MSUN / MPC_M ** 3   # kg/m^3 (colossus: M_sun h^2 / kpc^3)
    R200 = (3 * M200c * MSUN / (4 * np.pi * 200 * rho_c)) ** (1 / 3)  # m
    fb = cosmo.Ob0 / cosmo.Om0
    P200 = G_SI * M200c * MSUN * 200 * rho_c * fb / (2 * R200)        # Pa = J/m^3
    P200_keVcm3 = P200 / KEV_J / 1e6
    return P0 * x ** -0.3 * (1 + x) ** -beta * P200_keVcm3, R200


def bPe(z, Mmin=1e11, Mmax=1e16, xmax=4.0):
    """<b P_e> in meV/cm^3 (proper)."""
    lnM = np.linspace(np.log(Mmin), np.log(Mmax), 120)
    M200c = np.exp(lnM)
    # convert to M200m (h-units) for Tinker MF/bias
    Mh = M200c * h
    c = 5.0 * np.ones_like(Mh)
    M200m, R200m, c200m = mass_defs.changeMassDefinition(Mh, c, z, '200c', '200m')
    dndlnM = mass_function.massFunction(M200m, z, mdef='200m', model='tinker08', q_out='dndlnM')  # (h/Mpc)^3
    b = bias.haloBias(M200m, z, mdef='200m', model='tinker10')
    # dlnM200m/dlnM200c ~ 1 (ratio nearly constant); number density comoving -> proper
    n_proper = dndlnM * h ** 3 * (1 + z) ** 3  # per Mpc^3 proper per lnM
    x = np.linspace(1e-3, xmax, 2000)
    E = np.empty_like(M200c)
    for k, M in enumerate(M200c):
        P, R = battaglia_Pth(x, M, z)
        Pe = 0.518 * P                      # keV/cm^3
        Rcm = R * 100
        E[k] = simpson(4 * np.pi * (x * Rcm) ** 2 * Pe * Rcm, x=x)  # keV
    E_Mpc = E / (MPC_M * 100) ** 3          # keV per (cm^3/Mpc^3) -> keV/cm^3 * Mpc^3
    return simpson(n_proper * b * E_Mpc, x=lnM) * 1e6  # meV/cm^3


def camb_pk(zs):
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.66, ombh2=0.02242, omch2=0.11933, mnu=0.06)
    pars.InitPower.set_params(As=2.105e-9, ns=0.9665)
    pars.set_matter_power(redshifts=sorted(set(zs), reverse=True), kmax=50)
    pars.NonLinear = camb.model.NonLinear_both
    res = camb.get_results(pars)
    return camb.get_matter_power_interpolator(pars, nonlinear=True, hubble_units=False, k_hunit=False, kmax=50)


def beam_act(ell):
    f = os.path.join(HERE, 'data/act_y/ilc_beam.txt')
    if os.path.exists(f):
        bl = np.loadtxt(f)
        if bl.ndim == 2:
            return np.interp(ell, bl[:, 0], bl[:, 1] / bl[0, 1], right=0)
        return np.interp(ell, np.arange(len(bl)), bl / bl[0], right=0)
    s = np.radians(1.6 / 60) / np.sqrt(8 * np.log(2))
    return np.exp(-0.5 * ell * (ell + 1) * s ** 2)


def beam_gauss(ell, fwhm_arcmin):
    s = np.radians(fwhm_arcmin / 60) / np.sqrt(8 * np.log(2))
    return np.exp(-0.5 * ell * (ell + 1) * s ** 2)


def xi_template(zgrid, theta_arcmin, PK, bpe_of_z, beam):
    """xi_Fy(theta, z) per unit nothing: returns array [nz, ntheta] (dimensionless y units)."""
    ell = np.logspace(1, 5, 4000)
    th = np.radians(np.asarray(theta_arcmin) / 60)
    out = np.empty((len(zgrid), len(th)))
    for i, z in enumerate(zgrid):
        chi = cosmo.comovingDistance(0, z) / h   # Mpc
        a = 1 / (1 + z)
        Pk = PK.P(z, ell / chi)                  # Mpc^3
        # W(chi) = sigma_T/mec2 * a * <bP_e>; <bP_e> meV/cm^3 -> keV/m^3 ; integral over dchi in m
        W = SIGT / MEC2_KEV * a * bpe_of_z(z) * 1e-6 * 1e6 * MPC_M   # per Mpc of chi, times Mpc^-2 from 1/chi^2 below
        # units: W [1/Mpc]; P/chi^2 [Mpc]; l dl dimensionless -> xi dimensionless
        integrand = ell[None, :] / (2 * np.pi) * j0(ell[None, :] * th[:, None]) * beam(ell)[None, :] * Pk[None, :] / chi ** 2
        out[i] = b_F(z) * W * np.trapezoid(integrand, ell, axis=1)
    return out


def main():
    zs = np.round(np.arange(1.9, 3.85, 0.05), 3)
    bp = np.array([bPe(z) for z in zs])
    bp10 = np.array([bPe(z, Mmin=1e10) for z in zs])
    for z, v, v10 in zip(zs[::4], bp[::4], bp10[::4]):
        print(f"z={z:.2f}  <bP_e> = {v:.4f} meV/cm^3 (M>1e11)  {v10:.4f} (M>1e10)")
    print('check z=0:', bPe(0.0), 'z=0.5:', bPe(0.5), 'z=1:', bPe(1.0))
    PK = camb_pk(list(zs) + [0.0])
    bfun = lambda z: np.interp(z, zs, bp)
    tmpl_act = xi_template(zs, THETA_MID, PK, bfun, beam_act)
    tmpl_pl = xi_template(zs, THETA_MID, PK, bfun, lambda l: beam_gauss(l, 10.0))
    # an unbeamed reference at theta=0 bins too
    np.savez(os.path.join(HERE, 'data/lya/model.npz'), z=zs, bPe=bp, bPe_1e10=bp10, theta=THETA_MID,
             theta_edges=THETA_EDGES, xi_act=tmpl_act, xi_planck=tmpl_pl)
    for zb in ZBINS:
        sel = (zs >= zb[0]) & (zs < zb[1])
        print(zb, 'xi_act(theta) [1e-9]:', np.round(tmpl_act[sel].mean(0) * 1e9, 3))
        print(zb, 'xi_planck(theta) [1e-9]:', np.round(tmpl_pl[sel].mean(0) * 1e9, 3))


if __name__ == '__main__':
    main()
