"""R-153: physics and line model for the ALP decay-line search in SPHEREx images of dwarf spheroidals.

- Line intensity: Gamma = g^2 m^3 / (64 pi) (natural units); each decay gives two photons of m/2, so the energy
  intensity is I_line = Gamma D_Omega / (4 pi), D_Omega the line-of-sight integral of rho (GeV cm^-2 sr^-1).
  Stimulated decay (an enhancement from ambient photons) is not included: it raises the rate, so leaving it out is the
  conservative choice.
- D(<theta): log-log interpolation through Pace & Strigari 2019 (arXiv:1802.06811, Table 4: D(0.1), D(0.2), D(0.5 deg))
  and Song, Hiroshima & Murase 2024 (arXiv:2401.15606, Table 2: D(theta_Delta), the truncation angle); flat beyond
  theta_Delta (the halo truncated there). The ALT shape continues the 0.5 deg -> theta_Delta slope past theta_Delta.
- Line profile: each pixel's LVF passband is a Gaussian in lambda with FWHM = the pipeline bandwidth (lambda / R),
  centred on the pixel wavelength; the dwarf's velocity dispersion (< 12 km/s, 4e-5 in lambda) is negligible against
  it and its systemic velocity shifts the line by v/c.
"""
import numpy as np

HBAR_GEV_S = 6.582119569e-25
GEV_J = 1.602176634e-10
C_UM = 2.99792458e14
SQDEG_SR = (np.pi / 180) ** 2

# name: ra, dec, v_sys (km/s), Plummer r_h (arcmin, McConnachie 2012), Pace&Strigari log D(0.1, 0.2, 0.5 deg) with
# (lo, hi) errors on D(0.5), Song+2024 (theta_Delta deg, log D, -err, +err)
DWARFS = {
    "draco":    dict(ra=260.0517, dec=57.9153, v=-291.0, rh=10.0, ps=(17.42, 17.92, 18.54), ps_err=(0.14, 0.11), song=(1.41, 18.72, 0.42, 0.30)),
    "umi":      dict(ra=227.2854, dec=67.2225, v=-246.9, rh=8.2, ps=(17.34, 17.75, 18.20), ps_err=(0.08, 0.14), song=(1.19, 18.62, 0.34, 0.24)),
    "sculptor": dict(ra=15.0392, dec=-33.7092, v=111.4, rh=11.3, ps=(17.28, 17.71, 18.20), ps_err=(0.08, 0.09), song=(1.78, 18.67, 0.49, 0.33)),
    "fornax":   dict(ra=39.9971, dec=-34.4492, v=55.3, rh=16.6, ps=(17.14, 17.54, 17.97), ps_err=(0.04, 0.05), song=(2.27, 18.46, 0.38, 0.31)),
    "sextans":  dict(ra=153.2625, dec=-1.6147, v=224.2, rh=27.8, ps=(16.90, 17.36, 17.90), ps_err=(0.09, 0.11), song=(1.28, 18.17, 0.47, 0.46)),
    "carina":   dict(ra=100.4029, dec=-50.9661, v=222.9, rh=8.2, ps=(16.95, 17.39, 17.88), ps_err=(0.15, 0.18), song=(0.70, 17.69, 0.48, 0.56)),
    "leo1":     dict(ra=152.1171, dec=12.3064, v=282.9, rh=3.4, ps=(17.06, 17.50, 18.01), ps_err=(0.28, 0.20), song=(0.44, 17.69, 0.31, 0.23)),
    "leo2":     dict(ra=168.3700, dec=22.1517, v=78.0, rh=2.6, ps=(16.97, 17.30, 17.64), ps_err=(0.33, 0.50), song=(0.20, 17.22, 0.24, 0.22)),
    "bootes1":  dict(ra=210.0250, dec=14.5000, v=101.8, rh=12.6, ps=(17.05, 17.52, 18.11), ps_err=(0.30, 0.25), song=None),
}


def d_within(name, theta, shape="fid"):
    """Decay D-factor within theta (deg), GeV cm^-2."""
    p = DWARFS[name]
    th = [0.1, 0.2, 0.5]
    lg = list(p["ps"])
    if p["song"] is not None and p["song"][0] > 0.5:
        th.append(p["song"][0])
        lg.append(p["song"][1])
    th, lg = np.log10(th), np.array(lg)
    x = np.log10(np.maximum(np.asarray(theta, float), 1e-4))
    y = np.interp(x, th, lg)
    lo = x < th[0]
    y = np.where(lo, lg[0] + (lg[1] - lg[0]) / (th[1] - th[0]) * (x - th[0]), y)   # inner slope
    hi = x > th[-1]
    if shape == "alt":
        sl = (lg[-1] - lg[-2]) / (th[-1] - th[-2])
        y = np.where(hi, lg[-1] + sl * (x - th[-1]), y)
    else:
        y = np.where(hi, lg[-1], y)
    return 10 ** y


def ring_domega(name, edges, shape="fid"):
    """Mean D_Omega (GeV cm^-2 sr^-1) in each ring between the edges (deg)."""
    e = np.asarray(edges, float)
    dw = d_within(name, e, shape)
    dw[e == 0] = 0.0
    area = np.pi * np.diff(e ** 2) * SQDEG_SR
    return np.diff(dw) / area


def plummer_ring(name, edges):
    """Mean of a unit-total Plummer surface density in each ring (per sr): stellar-light radial weight."""
    a = DWARFS[name]["rh"] / 60.0
    e = np.asarray(edges, float)
    frac = e ** 2 / (e ** 2 + a ** 2)
    area = np.pi * np.diff(e ** 2) * SQDEG_SR
    return np.diff(frac) / area


def gamma_per_g2(m_ev):
    """Decay rate per g^2 (s^-1 per GeV^-2)."""
    m = m_ev * 1e-9
    return m ** 3 / (64 * np.pi) / HBAR_GEV_S


def line_intensity(m_ev, g, domega):
    """Line energy intensity I_line, W m^-2 sr^-1, for coupling g (GeV^-1) and D_Omega (GeV cm^-2 sr^-1)."""
    return gamma_per_g2(m_ev) * g ** 2 * domega / (4 * np.pi) * GEV_J * 1e4


def line_profile_nuinu(m_ev, lam_b, bw_b, v_kms=0.0, shape="gauss"):
    """nu I_nu (in units of I_line) seen by pixels of central wavelength lam_b and FWHM bw_b (um)."""
    lam0 = 1.23984198 / m_ev * 2 * (1 + v_kms / 2.99792458e5)  # photon energy m/2
    if shape == "gauss":
        s = bw_b / 2.35482
        phi = np.exp(-0.5 * ((lam_b - lam0) / s) ** 2) / (np.sqrt(2 * np.pi) * s)
    else:  # top hat of full width bw
        phi = (np.abs(lam_b - lam0) < bw_b / 2) / bw_b
    # I_nu = I_line phi_lambda lam0^2 / c ;  nu I_nu at lam_b = I_nu c / lam_b
    return phi * lam0 ** 2 / lam_b


def lam_line(m_ev):
    return 2 * 1.23984198 / m_ev
