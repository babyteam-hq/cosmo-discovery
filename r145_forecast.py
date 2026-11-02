"""R-145 step 0: shape-noise forecast of the ETG/LTG lensing-RAR ratio with DES Y3 shear around DESI DR1 BGS lenses.

For a lens table (r145_isolate.py output), each lens's baryonic mass (Brouwer et al. 2021: M* (1 + f_cold), log f_cold =
-0.69 log M* + 6.63) maps the fixed g_bar bins (Brouwer's 15 log bins over 1e-15..5e-12 m/s^2; the statistic uses the 7
with g_bar > ~1e-13) to physical annuli. Source pairs per annulus from DES Y3 metacal n_eff = 5.59 arcmin^-2 and
sigma_e = 0.261 per component (Amon et al. 2022 / Gatti et al. 2021), with a smooth n(z) of mean 0.63 (Myles et al. 2021)
and sources at z_s > z_l + 0.1. Expected signal: g_obs from the McGaugh et al. 2016 RAR, a0 = 1.2e-10; g_obs = 4 G DeltaSigma
(the SIS conversion of Brouwer et al.). Shape noise only (no LSS covariance, no boost dilution): an optimistic bound.

Statistic: Delta = mean over the 7 bins of log10(g_obs,ETG / g_obs,LTG). Prints sigma(Delta) and the separation of
Delta = 0 from Brouwer's 0.14 dex (Sersic split) and 0.19 dex (colour split), both for g_bar > 1e-13.

Usage: python -I r145_forecast.py <lens table> [<lens table> ...]
"""
import sys
import numpy as np, fitsio
from scipy.integrate import cumulative_trapezoid

G = 6.674e-11; MSUN = 1.989e30; MPC = 3.0857e22; C = 299792.458
H0, OM = 70.0, 0.3
_zg = np.linspace(0, 3, 3001)
_dc = cumulative_trapezoid(C / H0 / np.sqrt(OM * (1 + _zg) ** 3 + 1 - OM), _zg, initial=0)
def dcom(z): return np.interp(z, _zg, _dc)

NEFF, SIGE = 5.59, 0.261
zs = np.linspace(0.01, 2.5, 500)
pz = zs ** 2 * np.exp(-(zs / 0.36) ** 1.5); pz /= np.trapezoid(pz, zs)
SIGCRIT_FAC = (C * 1e3) ** 2 / (4 * np.pi * G) * MPC / MSUN  # Sigma_crit [Msun/Mpc^2] = FAC * D_s / (D_l D_ls) [1/Mpc]
EDGES = -15 + (np.log10(5e-12) + 15) / 15 * np.arange(16)
USE = EDGES[:-1] >= -13.1  # the 7 bins with g_bar > ~1e-13


def inv_sigcrit_and_frac(zl):
    """<Sigma_crit^-1> (Mpc^2/Msun, physical) over sources at z_s > z_l + 0.1, and the fraction of sources kept."""
    m = zs > zl + 0.1
    dl, ds = dcom(zl), dcom(zs[m])
    inv = (ds - dl) * dl / ds / (1 + zl) / SIGCRIT_FAC  # physical: D_l D_ls / D_s with D_A = D_c/(1+z), flat
    frac = np.trapezoid(pz[m], zs[m])
    return np.trapezoid(pz[m] * inv, zs[m]) / max(frac, 1e-9), frac


def rar(gb, a0=1.2e-10): return gb / (1 - np.exp(-np.sqrt(gb / a0)))


def forecast(t):
    logm = t['LOGM']
    mb = 10 ** logm * (1 + 10 ** (-0.69 * logm + 6.63))
    w = np.zeros(15)  # sum over pairs of Sigma_crit^-2 weights, per bin
    for zl, m in zip(t['Z'], mb):
        inv, frac = inv_sigcrit_and_frac(zl)
        r = np.sqrt(G * m * MSUN / 10 ** EDGES) / MPC  # Mpc physical, decreasing with g
        da = dcom(zl) / (1 + zl)
        area = np.pi * (r[:-1] ** 2 - r[1:] ** 2) / da ** 2 * (180 / np.pi * 60) ** 2  # arcmin^2
        w += area * NEFF * frac * inv ** 2
    gb = 10 ** (0.5 * (EDGES[:-1] + EDGES[1:]))
    ds = rar(gb) / (4 * G) * MPC ** 2 / MSUN  # Msun / Mpc^2
    sig_ds = SIGE / np.sqrt(np.maximum(w, 1e-30))
    sig_log = sig_ds / ds / np.log(10)
    return sig_log


def main():
    for path in sys.argv[1:]:
        t = fitsio.read(path)
        etg_c = t['QUIESCENT'] == 1
        print(f'{path}: N = {t.size}, quiescent (ETG by colour/spectrum) {etg_c.sum()}, star-forming {(~etg_c).sum()}')
        if etg_c.sum() < 2 or (~etg_c).sum() < 2:
            print('  too few lenses'); continue
        sE, sL = forecast(t[etg_c]), forecast(t[~etg_c])
        sD = np.sqrt(np.sum(sE[USE] ** 2 + sL[USE] ** 2)) / USE.sum()
        print('  per-bin sigma(log g_obs) ETG:', np.round(sE[USE], 3), ' LTG:', np.round(sL[USE], 3))
        print(f'  sigma(Delta) = {sD:.3f} dex; Brouwer Sersic 0.14 dex at {0.14 / sD:.1f} sigma, colour 0.19 dex at '
              f'{0.19 / sD:.1f} sigma (shape noise only)')


if __name__ == '__main__':
    main()
