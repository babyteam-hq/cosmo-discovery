#!/usr/bin/env python
"""R-123 (cosmo-deep-3): which dataset pays for the DESI DR2 w0wa preference.

Shapley attribution of delta-chi2(LCDM - w0waCDM) over every coalition of datasets, and a
specification curve of the significance across frozen analysis choices (`results-R123.md`).

Background-only pipeline: flat (w0, wa) CDM with photons, 2.0293 massless + 1 massive neutrino
(CAMB's convention, N_eff = 3.044), distances by direct quadrature; r_d, r_* and z_* from CAMB,
tabulated once on an (omega_b, omega_c, sum m_nu) grid (`data/r123/rd_table.npz`).

Data (gitignored, `data/r123/`): CobayaSampler/bao_data (DESI DR1 and DR2 BAO), CobayaSampler/sn_data
(Pantheon+, Union3, DES-SN5YR, DES-Dovekie). The compressed CMB is DESI DR2's (theta_*, omega_b,
omega_bc) prior (arXiv:2503.14738, App. A, eqs. 35-36); the alternative is the Planck 2018 distance
priors (R, l_A, omega_b) of Chen, Huang & Wang 2019 (arXiv:1808.05724, Table 1, wCDM row).

Run with data/r123/venv/bin/python. Subcommands: validate, reproduce, shapley, bootstrap, speccurve.
"""
import argparse
import itertools
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy import integrate, interpolate, optimize, stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "r123")
OUT = os.path.join(DATA, "out")
os.makedirs(OUT, exist_ok=True)

C_KMS = 299792.458
TCMB = 2.7255
OMEGA_G_H2 = 2.47282e-5 * (TCMB / 2.7255) ** 4  # photons
NEFF = 3.044
# one massive species carries N_eff/3, two massless carry the rest (CAMB default with mnu set)
NU_MASSLESS = NEFF - NEFF / 3.0
NU_MASSIVE_G = NEFF / 3.0
KB_EV = 8.617333262e-5
# CAMB's convention: every species at the standard T_nu = (4/11)^(1/3) T_cmb; the N_eff excess is a
# degeneracy factor (massless 2.0293, the massive one N_eff/3 = 1.0147)
TNU0 = (4.0 / 11.0) ** (1.0 / 3.0) * TCMB  # K
NU_REL = 7.0 / 8.0 * (4.0 / 11.0) ** (4.0 / 3.0)  # one standard massless species relative to photons

# ----------------------------------------------------------------------------- neutrino density


def _fermi_integral(y):
    f = lambda x: x * x * np.sqrt(x * x + y * y) / (np.exp(x) + 1.0)
    return integrate.quad(f, 0, 60, limit=200)[0]


_F0 = _fermi_integral(0.0)
_YGRID = np.concatenate([[0.0], np.logspace(-3, 4, 400)])
_FGRID = np.array([_fermi_integral(y) / _F0 for y in _YGRID])
_F_INTERP = interpolate.interp1d(np.log1p(_YGRID), np.log(_FGRID), kind="cubic")


def nu_density_factor(y):
    """rho_nu(massive) / rho_nu(massless) at the same temperature, y = m c^2 / (k T_nu)."""
    y = np.asarray(y, dtype=float)
    out = np.exp(_F_INTERP(np.log1p(np.minimum(y, _YGRID[-1]))))
    big = y > _YGRID[-1]
    if np.any(big):
        out = np.where(big, _FGRID[-1] * y / _YGRID[-1], out)
    return out


def _camb_mass_scale():
    """CAMB converts sum m_nu to omega_nu h^2 with neutrino_mass_fac = 94.07 (times (N_eff/3)^0.75), a
    0.3 per cent lower density per eV than the exact Fermi-Dirac integral at T_nu = (4/11)^(1/3) T_cmb.
    Rescale the mass so today's density matches CAMB's; the relativistic limit is unaffected."""
    m = 0.06
    camb_omnuh2 = m / 94.07 * (NEFF / 3.0) ** 0.75
    mine = OMEGA_G_H2 * NU_REL * NU_MASSIVE_G * nu_density_factor(m / (KB_EV * TNU0))
    return camb_omnuh2 / mine


MASS_SCALE = float(_camb_mass_scale())


# ----------------------------------------------------------------------------- r_d, r_*, z_* table

RD_TABLE = os.path.join(DATA, "rd_table.npz")
OB_GRID = np.linspace(0.005, 0.06, 23)
OC_GRID = np.linspace(0.01, 0.50, 50)
MNU_GRID = np.array([0.0, 0.03, 0.06, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5])


def build_rd_table():
    import camb

    rd = np.zeros((len(OB_GRID), len(OC_GRID), len(MNU_GRID)))
    rs = np.zeros_like(rd)
    zs = np.zeros_like(rd)
    t0 = time.time()
    for i, ob in enumerate(OB_GRID):
        for j, oc in enumerate(OC_GRID):
            for k, mnu in enumerate(MNU_GRID):
                p = camb.set_params(H0=67.0, ombh2=ob, omch2=oc, mnu=mnu, nnu=NEFF, tau=0.055,
                                    As=2.1e-9, ns=0.965)
                d = camb.get_background(p).get_derived_params()
                rd[i, j, k] = d["rdrag"]
                rs[i, j, k] = d["rstar"]
                zs[i, j, k] = d["zstar"]
    np.savez(RD_TABLE, ob=OB_GRID, oc=OC_GRID, mnu=MNU_GRID, rd=rd, rs=rs, zs=zs)
    print(f"rd table built in {time.time() - t0:.0f} s")


class RdTable:
    def __init__(self):
        if not os.path.exists(RD_TABLE):
            build_rd_table()
        t = np.load(RD_TABLE)
        pts = (t["ob"], t["oc"], t["mnu"])
        kw = dict(method="cubic", bounds_error=False, fill_value=None)
        self.rd = interpolate.RegularGridInterpolator(pts, t["rd"], **kw)
        self.rs = interpolate.RegularGridInterpolator(pts, t["rs"], **kw)
        self.zs = interpolate.RegularGridInterpolator(pts, t["zs"], **kw)

    def __call__(self, ob, oc, mnu):
        x = np.array([[ob, oc, mnu]])
        return float(self.rd(x)[0]), float(self.rs(x)[0]), float(self.zs(x)[0])


_RD = None


def rd_table():
    global _RD
    if _RD is None:
        _RD = RdTable()
    return _RD


# ----------------------------------------------------------------------------- background

U_GRID = np.linspace(0.0, math.log(1300.0), 4000)  # u = ln(1+z)
Z_GRID = np.expm1(U_GRID)
A_GRID = 1.0 / (1.0 + Z_GRID)


class Background:
    """Flat w0waCDM. Parameters: ob = omega_b h^2, oc = omega_c h^2, h, w0, wa, mnu (eV, sum)."""

    def __init__(self, ob, oc, h, w0=-1.0, wa=0.0, mnu=0.06):
        self.ob, self.oc, self.h, self.w0, self.wa, self.mnu = ob, oc, h, w0, wa, mnu
        h2 = h * h
        og = OMEGA_G_H2 / h2
        onu_massless = og * NU_REL * NU_MASSLESS
        onu_massive0 = og * NU_REL * NU_MASSIVE_G  # massless-equivalent density of the massive species
        self.y0 = MASS_SCALE * mnu / (KB_EV * TNU0) if mnu > 0 else 0.0  # m c^2 / k T_nu today
        ocb = (ob + oc) / h2
        a = A_GRID
        fnu = nu_density_factor(self.y0 * a) if mnu > 0 else 1.0
        rho_nu = onu_massive0 * a ** -4 * fnu
        rho_r = (og + onu_massless) * a ** -4
        rho_m = ocb * a ** -3
        onu0 = onu_massive0 * (nu_density_factor(self.y0) if mnu > 0 else 1.0)
        self.omega_m = ocb + onu0  # includes massive neutrinos today (CAMB's omegam)
        self.omega_de = 1.0 - ocb - og - onu_massless - onu0
        fde = a ** (-3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * (1.0 - a))
        e2 = rho_r + rho_nu + rho_m + self.omega_de * fde
        self.E = np.sqrt(e2)
        # comoving distance in Mpc: D_C = (c/H0) int dz/E = (c/H0) int (1+z)/E du
        integrand = (1.0 + Z_GRID) / self.E
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(U_GRID))])
        self.dc_grid = (C_KMS / (100.0 * h)) * cum
        self.omnuh2 = onu0 * h2

    def E_of_z(self, z):
        return np.interp(np.log1p(z), U_GRID, self.E)

    def DM(self, z):
        return np.interp(np.log1p(z), U_GRID, self.dc_grid)

    def DH(self, z):
        return C_KMS / (100.0 * self.h * self.E_of_z(z))

    def DV(self, z):
        return (z * self.DM(z) ** 2 * self.DH(z)) ** (1.0 / 3.0)

    def DL(self, zhd, zhel):
        return (1.0 + zhel) * self.DM(zhd)

    def sound_horizons(self):
        return rd_table()(self.ob, self.oc, self.mnu)  # rd, rstar, zstar

    def theta_star(self):
        rd, rs, zs = self.sound_horizons()
        return rs / self.DM(zs)

    def distance_priors(self):
        rd, rs, zs = self.sound_horizons()
        dm = self.DM(zs)
        R = math.sqrt(self.omega_m) * 100.0 * self.h * dm / C_KMS
        lA = math.pi * dm / rs
        return R, lA


# ----------------------------------------------------------------------------- likelihoods


def read_cov(path):
    with open(path) as f:
        n = int(f.readline().split()[0])
    vals = np.loadtxt(path, skiprows=1)
    return vals.reshape(n, n)


class Gaussian:
    """chi2 = (d - m)^T Cinv (d - m); subclasses define model(bg)."""

    name = "gaussian"

    def chi2(self, bg):
        r = self.data - self.model(bg)
        return float(r @ self.cinv @ r)

    def simulate(self, bg, rng):
        new = self.__class__.__new__(self.__class__)
        new.__dict__.update(self.__dict__)
        new.data = self.model(bg) + self.chol @ rng.standard_normal(len(self.data))
        return new

    def finish(self):
        self.cinv = np.linalg.inv(self.cov)
        self.chol = np.linalg.cholesky(self.cov)


class BAO(Gaussian):
    def __init__(self, release="dr2"):
        self.name = "DESI-DR2" if release == "dr2" else "DESI-DR1"
        if release == "dr2":
            base = os.path.join(DATA, "bao_data", "desi_bao_dr2", "desi_gaussian_bao_ALL_GCcomb")
        else:
            base = os.path.join(DATA, "bao_data", "desi_2024_gaussian_bao_ALL_GCcomb")
        rows = [l.split() for l in open(base + "_mean.txt") if not l.startswith("#") and l.strip()]
        self.z = np.array([float(r[0]) for r in rows])
        self.data = np.array([float(r[1]) for r in rows])
        self.kind = [r[2] for r in rows]
        self.cov = np.loadtxt(base + "_cov.txt")
        self.finish()

    def model(self, bg):
        rd = bg.sound_horizons()[0]
        out = np.empty(len(self.z))
        for i, (z, k) in enumerate(zip(self.z, self.kind)):
            if k == "DM_over_rs":
                out[i] = bg.DM(z) / rd
            elif k == "DH_over_rs":
                out[i] = bg.DH(z) / rd
            else:
                out[i] = bg.DV(z) / rd
        return out


class SN(Gaussian):
    """SN distance moduli with the absolute magnitude marginalised analytically (as cobaya does)."""

    def __init__(self, sample, zmin=None):
        self.name = sample
        d = os.path.join(DATA, "sn_data")
        if sample == "PantheonPlus":
            cols = open(os.path.join(d, "PantheonPlus", "Pantheon+SH0ES.dat")).readline().split()
            arr = np.genfromtxt(os.path.join(d, "PantheonPlus", "Pantheon+SH0ES.dat"), skip_header=1,
                                usecols=[cols.index("zHD"), cols.index("zHEL"), cols.index("m_b_corr")])
            zhd, zhel, mag = arr.T
            cov = read_cov(os.path.join(d, "PantheonPlus", "Pantheon+SH0ES_STAT+SYS.cov"))
            keep = zhd > (0.01 if zmin is None else zmin)
        elif sample == "DESY5":
            arr = np.genfromtxt(os.path.join(d, "DESY5", "DES-SN5YR_HD.csv"), delimiter=",", names=True)
            zhd, zhel, mag, err = arr["zHD"], arr["zHEL"], arr["MU"], arr["MUERR_FINAL"]
            cov = read_cov(os.path.join(d, "DESY5", "covsys_000.txt")) + np.diag(err ** 2)
            keep = zhd > (0.0 if zmin is None else zmin)
        elif sample == "DESDovekie":
            arr = np.genfromtxt(os.path.join(d, "DES-Dovekie", "DES-Dovekie_HD.csv"), delimiter=",", names=True)
            zhd, zhel, mag = arr["zHD"], arr["zHEL"], arr["MU"]
            npz = np.load(os.path.join(d, "DES-Dovekie", "covtot_inv_000.npz"))
            n = int(npz["nsn"][0])  # upper triangle of the inverse total covariance (cobaya's reader)
            cinv = np.zeros((n, n))
            cinv[np.triu_indices(n)] = npz["cov"].astype(float)
            il = np.tril_indices(n, -1)
            cinv[il] = cinv.T[il]
            cov = np.linalg.inv(cinv)
            keep = zhd > (0.0 if zmin is None else zmin)
        elif sample == "Union3":
            rows = [l.split() for l in open(os.path.join(d, "Union3", "lcparam_full.txt"))
                    if l.strip() and not l.startswith("#")]
            zhd = np.array([float(r[1]) for r in rows])  # columns: name zcmb zhel dz mb ...
            zhel = zhd.copy()
            mag = np.array([float(r[4]) for r in rows])
            cov = read_cov(os.path.join(d, "Union3", "mag_covmat.txt"))
            keep = zhd > (0.0 if zmin is None else zmin)
        else:
            raise ValueError(sample)
        self.zhd, self.zhel, self.data = zhd[keep], zhel[keep], mag[keep]
        self.cov = cov[np.ix_(keep, keep)]
        self.finish()
        # marginalise the absolute magnitude (flat prior): project the constant mode out of Cinv
        cinv = self.cinv
        u = cinv @ np.ones(len(self.data))
        self.cinv = cinv - np.outer(u, u) / u.sum()

    def model(self, bg):
        return 5.0 * np.log10(bg.DL(self.zhd, self.zhel)) + 25.0


class CMBCompressed(Gaussian):
    """DESI DR2's (theta_*, omega_b, omega_bc) Gaussian prior, arXiv:2503.14738 eqs. 35-36."""

    name = "CMB-c"

    def __init__(self):
        self.data = np.array([0.01041, 0.02223, 0.14208])
        self.cov = 1e-9 * np.array([[0.006621, 0.12444, -1.1929],
                                    [0.12444, 21.344, -94.001],
                                    [-1.1929, -94.001, 1488.4]])
        self.finish()

    def model(self, bg):
        return np.array([bg.theta_star(), bg.ob, bg.ob + bg.oc])


class CMBDistancePriors(Gaussian):
    """Planck 2018 TT,TE,EE+lowE distance priors (R, l_A, omega_b), Chen+2019 Table 1, wCDM row."""

    name = "CMB-RlA"

    def __init__(self):
        self.data = np.array([1.7493, 301.462, 0.02239])
        sig = np.array([0.00465, 0.0895, 0.00015])
        corr = np.array([[1.0, 0.47, -0.66], [0.47, 1.0, -0.34], [-0.66, -0.34, 1.0]])
        self.cov = corr * np.outer(sig, sig)
        self.finish()

    def model(self, bg):
        R, lA = bg.distance_priors()
        return np.array([R, lA, bg.ob])


_LIKES = {}


def like(key):
    """Cached likelihood objects by key: 'dr2','dr1','cmbc','cmbrla','PantheonPlus','Union3','DESY5',
    'DESDovekie', and '<sn>@z0.1' for a low-z cut."""
    if key not in _LIKES:
        if key == "dr2":
            _LIKES[key] = BAO("dr2")
        elif key == "dr1":
            _LIKES[key] = BAO("dr1")
        elif key == "cmbc":
            _LIKES[key] = CMBCompressed()
        elif key == "cmbrla":
            _LIKES[key] = CMBDistancePriors()
        elif "@z" in key:
            s, z = key.split("@z")
            _LIKES[key] = SN(s, zmin=float(z))
        else:
            _LIKES[key] = SN(key)
    return _LIKES[key]


# ----------------------------------------------------------------------------- fitting

# DESI's priors: w0 in [-3, 1], wa in [-3, 2], w0 + wa < 0 (arXiv:2503.14738 Table 3)
PRIORS = {
    "desi": dict(w0=(-3.0, 1.0), wa=(-3.0, 2.0), sumcut=True),
    "desi-nocut": dict(w0=(-3.0, 1.0), wa=(-3.0, 2.0), sumcut=False),
    "narrow": dict(w0=(-2.0, 0.0), wa=(-2.0, 1.0), sumcut=True),
}
BOUNDS = dict(ob=(0.005, 0.06), oc=(0.01, 0.5), h=(0.3, 1.2), mnu=(0.0, 0.5))


class Problem:
    """A coalition of likelihoods and a model ('lcdm', 'w0wa', 'wcdm'), with options."""

    def __init__(self, keys, model="w0wa", prior="desi", free_mnu=False, mnu=0.06, likes=None):
        self.keys = list(keys)
        self.likes = likes if likes is not None else [like(k) for k in keys]
        self.model = model
        self.prior = PRIORS[prior]
        self.free_mnu = free_mnu
        self.mnu = mnu
        self.bbn = not any(k.startswith("cmb") for k in self.keys)
        self.names = ["ob", "oc", "h"]
        if model in ("w0wa", "wcdm"):
            self.names.append("w0")
        if model == "w0wa":
            self.names.append("wa")
        if free_mnu:
            self.names.append("mnu")
        self.bounds = [BOUNDS.get(n) or self.prior[n] for n in self.names]

    def background(self, x):
        p = dict(zip(self.names, x))
        return Background(p["ob"], p["oc"], p["h"], p.get("w0", -1.0), p.get("wa", 0.0),
                          p.get("mnu", self.mnu))

    def chi2(self, x, parts=False):
        p = dict(zip(self.names, x))
        pen = 0.0
        for n, (lo, hi) in zip(self.names, self.bounds):
            v = p[n]
            if v < lo:
                pen += 1e6 * (lo - v) ** 2
            if v > hi:
                pen += 1e6 * (v - hi) ** 2
        if self.prior["sumcut"] and "wa" in p and p["w0"] + p["wa"] > 0:
            pen += 1e6 * (p["w0"] + p["wa"]) ** 2
        xc = [min(max(p[n], lo), hi) for n, (lo, hi) in zip(self.names, self.bounds)]
        bg = self.background(xc)
        vals = [lk.chi2(bg) for lk in self.likes]
        if parts:
            return vals
        if self.bbn:  # DESI's BBN prior on omega_b when no CMB is in the coalition (Schoeneberg 2024)
            pen += ((p["ob"] - 0.02218) / 0.00055) ** 2
        return sum(vals) + pen

    def starts(self, n, rng):
        base = dict(ob=0.0222, oc=0.119, h=0.68, w0=-0.75, wa=-0.8, mnu=0.06)
        out = [np.array([base[k] for k in self.names])]
        if "w0" in self.names:
            out.append(np.array([base[k] if k not in ("w0", "wa") else {"w0": -1.0, "wa": 0.0}[k]
                                 for k in self.names]))
        while len(out) < n:
            x = []
            for k, (lo, hi) in zip(self.names, self.bounds):
                if k in ("w0", "wa"):
                    x.append(rng.uniform(lo, hi))
                elif k == "h":
                    x.append(rng.uniform(0.55, 0.8))
                elif k == "ob":
                    x.append(rng.uniform(0.019, 0.026))
                elif k == "oc":
                    x.append(rng.uniform(0.09, 0.16))
                else:
                    x.append(rng.uniform(0.0, 0.3))
            x = np.array(x)
            if "wa" in self.names and self.prior["sumcut"]:
                i0, ia = self.names.index("w0"), self.names.index("wa")
                if x[i0] + x[ia] > 0:
                    x[ia] = -abs(x[i0]) - rng.uniform(0, 0.5)
            out.append(x)
        return out

    def fit(self, nstart=6, seed=0, x0=None):
        rng = np.random.default_rng(seed)
        starts = self.starts(nstart, rng) if nstart > 0 else []
        if x0 is not None:
            starts = [np.asarray(x0)] + starts
        best = None
        for s in starts:
            r = optimize.minimize(self.chi2, s, method="Nelder-Mead",
                                  options=dict(xatol=1e-5, fatol=1e-6, maxiter=6000, maxfev=6000,
                                               adaptive=True))
            r = optimize.minimize(self.chi2, r.x, method="Powell",
                                  options=dict(xtol=1e-5, ftol=1e-8, maxiter=20000))
            if best is None or r.fun < best.fun:
                best = r
        return float(best.fun), dict(zip(self.names, [float(v) for v in best.x]))


def significance(dchi2, dof):
    """Gaussian-equivalent significance of a delta-chi2 with `dof` extra parameters."""
    if dchi2 <= 0:
        return 0.0
    p = stats.chi2.sf(dchi2, dof)
    return float(stats.norm.isf(p / 2.0))


def fit_pair(keys, prior="desi", free_mnu=False, model="w0wa", nstart=6, likes=None, x0=None):
    """chi2_min for LCDM and the extended model on a coalition; returns dict."""
    p0 = Problem(keys, "lcdm", prior, free_mnu, likes=likes)
    p1 = Problem(keys, model, prior, free_mnu, likes=likes)
    c0, x0_ = p0.fit(nstart=(max(2, nstart // 2) if nstart > 0 else 0), x0=(x0 or {}).get("lcdm"))
    c1, x1_ = p1.fit(nstart=nstart, x0=(x0 or {}).get(model))
    d = c0 - c1
    dof = 2 if model == "w0wa" else 1
    return dict(keys=list(keys), chi2_lcdm=c0, chi2_ext=c1, dchi2=d, sigma=significance(d, dof),
                x_lcdm=x0_, x_ext=x1_,
                parts_lcdm=p0.chi2([x0_[n] for n in p0.names], parts=True),
                parts_ext=p1.chi2([x1_[n] for n in p1.names], parts=True))


# ----------------------------------------------------------------------------- Shapley


def shapley(values, players):
    """values: dict frozenset(players) -> v(S); exact Shapley value per player."""
    n = len(players)
    phi = {}
    for i in players:
        tot = 0.0
        others = [p for p in players if p != i]
        for r in range(n):
            w = math.factorial(r) * math.factorial(n - r - 1) / math.factorial(n)
            for S in itertools.combinations(others, r):
                S = frozenset(S)
                tot += w * (values[S | {i}] - values[S])
        phi[i] = tot
    return phi


def all_coalitions(players):
    for r in range(len(players) + 1):
        for S in itertools.combinations(players, r):
            yield frozenset(S)


def _fit_coalition(args):
    keys, prior, free_mnu, nstart = args
    if not keys:
        return dict(keys=[], dchi2=0.0, chi2_lcdm=0.0, chi2_ext=0.0, sigma=0.0)
    t = time.time()
    r = fit_pair(keys, prior=prior, free_mnu=free_mnu, nstart=nstart)
    r["seconds"] = time.time() - t
    return r


def run_shapley(players, prior="desi", free_mnu=False, nstart=6, workers=4, tag="shapley"):
    coalitions = [tuple(sorted(S, key=players.index)) for S in all_coalitions(players)]
    args = [(c, prior, free_mnu, nstart) for c in coalitions]
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        results = list(ex.map(_fit_coalition, args))
    values = {frozenset(r["keys"]): r["dchi2"] for r in results}
    phi = shapley(values, players)
    out = dict(players=players, prior=prior, free_mnu=free_mnu, nstart=nstart,
               coalitions=results, shapley=phi, grand=values[frozenset(players)],
               seconds=time.time() - t0)
    with open(os.path.join(OUT, f"{tag}.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    return out


def print_shapley(out):
    players = out["players"]
    print(f"\nGrand coalition delta-chi2 = {out['grand']:.2f}  (sum of Shapley = "
          f"{sum(out['shapley'].values()):.2f}); {out['seconds']:.0f} s")
    for p in players:
        print(f"  {p:14s} phi = {out['shapley'][p]:6.2f}  ({100 * out['shapley'][p] / out['grand']:5.1f} %)")
    print("\ncoalition values:")
    for r in sorted(out["coalitions"], key=lambda r: (len(r["keys"]), r["keys"])):
        if r["keys"]:
            print(f"  {'+'.join(r['keys']):55s} dchi2 = {r['dchi2']:6.2f}  ({r['sigma']:.1f} sigma)"
                  f"  w0 = {r['x_ext']['w0']:6.3f} wa = {r['x_ext']['wa']:6.3f}")


# ----------------------------------------------------------------------------- bootstrap


def _bootstrap_one(args):
    players, prior, free_mnu, nstart, seed, truth, x0 = args
    rng = np.random.default_rng(seed)
    bg = Background(**truth)
    sims = {k: like(k).simulate(bg, rng) for k in players}
    values = {}
    for S in all_coalitions(players):
        if not S:
            values[S] = 0.0
            continue
        keys = [p for p in players if p in S]
        r = fit_pair(keys, prior=prior, free_mnu=free_mnu, nstart=nstart,
                     likes=[sims[k] for k in keys], x0=x0.get("+".join(keys)))
        values[S] = r["dchi2"]
    phi = shapley(values, players)
    return dict(seed=seed, shapley=phi, grand=values[frozenset(players)],
                values={"+".join(p for p in players if p in S): v for S, v in values.items()})


def run_bootstrap(players, truth, nboot=100, prior="desi", free_mnu=False, nstart=2, workers=4,
                  tag="bootstrap", x0=None):
    args = [(players, prior, free_mnu, nstart, 1000 + i, truth, x0 or {}) for i in range(nboot)]
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(_bootstrap_one, args))
    out = dict(players=players, truth=truth, nboot=nboot, prior=prior, results=res,
               seconds=time.time() - t0)
    with open(os.path.join(OUT, f"{tag}.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    return out


# ----------------------------------------------------------------------------- specification curve

SPEC_FACTORS = dict(
    sn=["none", "PantheonPlus", "Union3", "DESY5", "DESDovekie"],
    sncut=["released", "z0.1"],
    bao=["dr2", "dr1"],
    cmb=["none", "cmbc", "cmbrla"],
    de=["w0wa", "wcdm"],
    prior=["desi", "desi-nocut", "narrow"],
    nu=["fixed", "free"],
)


def spec_grid():
    specs = []
    for combo in itertools.product(*SPEC_FACTORS.values()):
        s = dict(zip(SPEC_FACTORS.keys(), combo))
        if s["sn"] == "none" and s["sncut"] != "released":
            continue  # collapses
        if s["de"] == "wcdm" and s["prior"] == "desi-nocut":
            continue  # the w0 + wa cut is moot for wCDM
        specs.append(s)
    return specs


def spec_keys(s):
    keys = [s["bao"]]
    if s["cmb"] != "none":
        keys.append(s["cmb"])
    if s["sn"] != "none":
        keys.append(s["sn"] + ("@z0.1" if s["sncut"] == "z0.1" else ""))
    return keys


# small steps: with the 20x larger ones the Hessian is not positive definite (the chi2 is not quadratic over
# 0.03 in w0 and 0.1 in wa); at these steps sigma(w0), sigma(wa) for DESI+CMB_c+DESY5 are 0.059, 0.24, DESI's values
HESS_STEPS = dict(ob=1.5e-5, oc=1e-4, h=2.5e-4, w0=1.5e-3, wa=5e-3, mnu=1e-3)


def param_shift_sigma(prob, x_ext, dof):
    """Distance of (w0, wa) = (-1, 0) from the best fit in the Gaussian (Hessian) metric: the full-parameter
    Hessian of chi2 at the best fit, inverted (pseudo-inverse: unconstrained directions such as h in an SN-only
    fit carry no curvature), the (w0, wa) block taken as the marginal covariance."""
    names = prob.names
    x0 = np.array([x_ext[n] for n in names])
    h = np.array([HESS_STEPS[n] for n in names])
    n = len(names)
    f0 = prob.chi2(x0)
    H = np.zeros((n, n))
    for a in range(n):
        ea = np.zeros(n); ea[a] = h[a]
        H[a, a] = (prob.chi2(x0 + ea) - 2 * f0 + prob.chi2(x0 - ea)) / h[a] ** 2
        for b in range(a + 1, n):
            eb = np.zeros(n); eb[b] = h[b]
            H[a, b] = H[b, a] = (prob.chi2(x0 + ea + eb) - prob.chi2(x0 + ea - eb) - prob.chi2(x0 - ea + eb)
                                 + prob.chi2(x0 - ea - eb)) / (4 * h[a] * h[b])
    F = 0.5 * H
    # scale to unit steps for a well-conditioned pseudo-inverse
    S = np.diag(h)
    Fs = S @ F @ S
    Cs = np.linalg.pinv(Fs, rcond=1e-8, hermitian=True)
    C = S @ Cs @ S
    idx = [names.index(k) for k in ("w0", "wa") if k in names]
    Cw = C[np.ix_(idx, idx)]
    d = np.array([-1.0, 0.0][:len(idx)]) - x0[idx]
    try:
        d2 = float(d @ np.linalg.solve(Cw, d))
    except Exception:
        d2 = float("nan")
    return significance(d2, dof) if np.isfinite(d2) and d2 > 0 else 0.0, d2


def _run_spec(s):
    t = time.time()
    keys = spec_keys(s)
    try:
        r = fit_pair(keys, prior=s["prior"], free_mnu=(s["nu"] == "free"), model=s["de"], nstart=3)
        dof = 2 if s["de"] == "w0wa" else 1
        p1 = Problem(keys, s["de"], s["prior"], s["nu"] == "free")
        ps_sigma, d2 = param_shift_sigma(p1, r["x_ext"], dof)
        out = dict(spec=s, keys=keys, dchi2=r["dchi2"], sigma_dchi2=r["sigma"], sigma_shift=ps_sigma, d2_shift=d2,
                   w0=r["x_ext"].get("w0"), wa=r["x_ext"].get("wa"), seconds=time.time() - t)
    except Exception as e:
        out = dict(spec=s, keys=keys, error=repr(e), seconds=time.time() - t)
    return out


def run_speccurve(workers=4, tag="speccurve", limit=None):
    specs = spec_grid()
    if limit:
        specs = specs[:limit]
    print(f"{len(specs)} specifications", flush=True)
    t0 = time.time()
    results = []
    with ProcessPoolExecutor(workers) as ex:
        for i, r in enumerate(ex.map(_run_spec, specs, chunksize=2)):
            results.append(r)
            if (i + 1) % 25 == 0:
                print(f"  {i + 1}/{len(specs)} done, {time.time() - t0:.0f} s", flush=True)
                json.dump(results, open(os.path.join(OUT, f"{tag}.json"), "w"), indent=0, default=float)
    json.dump(results, open(os.path.join(OUT, f"{tag}.json"), "w"), indent=0, default=float)
    return results


def anova(results, metric="sigma_dchi2", subset=None):
    """Main-effects decomposition of the variance of `metric` over the grid (type-I-free: the grid is nearly
    balanced; we report sum of squares of factor-level means around the grand mean, over the total)."""
    rows = [r for r in results if "error" not in r and r.get(metric) is not None and np.isfinite(r[metric])
            and (subset is None or subset(r["spec"]))]
    y = np.array([r[metric] for r in rows])
    gm = y.mean()
    sst = float(((y - gm) ** 2).sum())
    out = dict(n=len(rows), mean=float(gm), std=float(y.std()), min=float(y.min()), max=float(y.max()), factors={})
    for f, levels in SPEC_FACTORS.items():
        ss = 0.0
        means = {}
        for lv in levels:
            sel = np.array([r["spec"][f] == lv for r in rows])
            if sel.sum() == 0:
                continue
            m = y[sel].mean()
            means[lv] = (float(m), int(sel.sum()))
            ss += sel.sum() * (m - gm) ** 2
        out["factors"][f] = dict(share=float(ss / sst) if sst > 0 else 0.0, means=means)
    # the two-way SN x CMB interaction, the obvious one
    inter = {}
    for a in SPEC_FACTORS["sn"]:
        for b in SPEC_FACTORS["cmb"]:
            sel = np.array([r["spec"]["sn"] == a and r["spec"]["cmb"] == b for r in rows])
            if sel.sum():
                inter[f"{a}|{b}"] = (float(y[sel].mean()), float(y[sel].min()), float(y[sel].max()), int(sel.sum()))
    out["sn_x_cmb"] = inter
    return out


def cmd_speccurve(args):
    if args.analyse_only:
        results = json.load(open(os.path.join(OUT, f"{args.tag}.json")))
    else:
        results = run_speccurve(workers=args.workers, tag=args.tag, limit=args.limit)
    for metric in ("sigma_dchi2", "sigma_shift"):
        a = anova(results, metric)
        print(f"\n== {metric}: n = {a['n']}, mean {a['mean']:.2f}, sd {a['std']:.2f}, range {a['min']:.2f}-{a['max']:.2f}")
        for f, d in sorted(a["factors"].items(), key=lambda kv: -kv[1]["share"]):
            lv = ", ".join(f"{k} {v[0]:.2f}" for k, v in d["means"].items())
            print(f"  {f:6s} share {100 * d['share']:5.1f} %   {lv}")
        print("  SN x CMB cell means (min-max):")
        for k, v in a["sn_x_cmb"].items():
            print(f"    {k:22s} {v[0]:.2f} ({v[1]:.2f}-{v[2]:.2f}, n={v[3]})")
        json.dump(a, open(os.path.join(OUT, f"{args.tag}_anova_{metric}.json"), "w"), indent=1)
    errs = [r for r in results if "error" in r]
    print(f"\n{len(errs)} failed specs", errs[:3])
    bad = [r for r in results if "error" not in r and not np.isfinite(r.get("sigma_shift", 0.0))]
    print(f"{len(bad)} specs with a non-finite shift metric (best fit on a prior edge):",
          [(r["spec"]["sn"], r["spec"]["cmb"], r["spec"]["de"], r["spec"]["prior"], round(r["wa"] or 0, 2)) for r in bad[:8]])
    # the balanced sub-grid with an SN sample (the low-z cut factor is defined there)
    a = anova(results, "sigma_dchi2", subset=lambda s: s["sn"] != "none")
    print(f"\n== sigma_dchi2 on the SN != none sub-grid: n = {a['n']}, mean {a['mean']:.2f}, sd {a['std']:.2f}")
    for f, d in sorted(a["factors"].items(), key=lambda kv: -kv[1]["share"]):
        print(f"  {f:6s} share {100 * d['share']:5.1f} %   " + ", ".join(f"{k} {v[0]:.2f}" for k, v in d["means"].items()))
    # the DESI DR2 + CMB_c + SN headline slice: what each other choice does to it
    print("\n== headline slice (dr2, cmbc, released cut, w0wa, desi prior, fixed nu) and one-factor changes:")
    base = dict(sncut="released", bao="dr2", cmb="cmbc", de="w0wa", prior="desi", nu="fixed")
    for sn in SPEC_FACTORS["sn"]:
        line = []
        for r in results:
            s = r["spec"]
            if s["sn"] != sn or "error" in r:
                continue
            diff = [k for k, v in base.items() if s[k] != v]
            if len(diff) == 0:
                line.insert(0, f"base {r['sigma_dchi2']:.2f}")
            elif len(diff) == 1:
                line.append(f"{diff[0]}={s[diff[0]]} {r['sigma_dchi2']:.2f}")
        print(f"  {sn:13s} " + "; ".join(line))
    # the ranked curve
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        rows = sorted([r for r in results if "error" not in r], key=lambda r: r["sigma_dchi2"])
        cols = {"none": "0.5", "PantheonPlus": "tab:blue", "Union3": "tab:green", "DESY5": "tab:red", "DESDovekie": "tab:orange"}
        fig, ax = plt.subplots(figsize=(9, 4.2))
        for i, r in enumerate(rows):
            ax.plot(i, r["sigma_dchi2"], ".", ms=4, color=cols[r["spec"]["sn"]])
            if r["spec"]["sncut"] == "z0.1":
                ax.plot(i, -0.15, "|", ms=4, color="k")
            if r["spec"]["de"] == "wcdm":
                ax.plot(i, -0.3, "|", ms=4, color="tab:purple")
        for k, c in cols.items():
            ax.plot([], [], ".", color=c, label=k)
        ax.plot([], [], "|", color="k", label="SN z > 0.1 cut")
        ax.plot([], [], "|", color="tab:purple", label="wCDM (1 dof)")
        ax.set_xlabel("specification, ranked")
        ax.set_ylabel("significance vs LCDM (sigma)")
        ax.set_title("DESI BAO (+CMB prior) (+SN): 540 analysis specifications, compressed likelihoods")
        ax.legend(fontsize=7, ncol=4, loc="upper left")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(HERE, "desi_speccurve.png"), dpi=130)
        print("figure: desi_speccurve.png")
    except Exception as e:
        print("no figure:", e)


# ----------------------------------------------------------------------------- CLI


def cmd_validate(args):
    """Compare the background against CAMB at a few points."""
    import camb

    for (ob, oc, h, w0, wa, mnu) in [(0.0224, 0.12, 0.675, -1.0, 0.0, 0.06), (0.0222, 0.118, 0.64, -0.7, -1.0, 0.06),
                                     (0.0230, 0.13, 0.70, -1.2, 0.5, 0.06), (0.0224, 0.12, 0.675, -0.8, -0.5, 0.3),
                                     (0.0224, 0.12, 0.675, -1.0, 0.0, 0.0)]:
        bg = Background(ob, oc, h, w0, wa, mnu)
        p = camb.set_params(H0=100 * h, ombh2=ob, omch2=oc, mnu=mnu, nnu=NEFF, tau=0.055, As=2.1e-9, ns=0.965,
                            dark_energy_model="ppf", w=w0, wa=wa)
        r = camb.get_background(p)
        d = r.get_derived_params()
        zs = [0.1, 0.5, 1.0, 2.33, d["zstar"]]
        dm_camb = r.angular_diameter_distance(zs) * (1 + np.array(zs))
        dm_me = bg.DM(np.array(zs))
        hz = r.hubble_parameter(np.array(zs[:4])) / (100 * h)
        print(f"ob={ob} oc={oc} h={h} w0={w0} wa={wa} mnu={mnu}")
        print("  DM rel err:", np.array2string((dm_me / dm_camb - 1), precision=2))
        print("  E  rel err:", np.array2string((bg.E_of_z(np.array(zs[:4])) / hz - 1), precision=2))
        rd, rs, zst = bg.sound_horizons()
        print(f"  rd {rd:.3f} vs {d['rdrag']:.3f}; rs {rs:.3f} vs {d['rstar']:.3f}; z* {zst:.2f} vs {d['zstar']:.2f};"
              f" theta* {bg.theta_star():.7f} vs {d['thetastar'] / 100:.7f}  ({bg.theta_star() / (d['thetastar'] / 100) - 1:.1e});"
              f" omnuh2 {bg.omnuh2:.5f} vs {r.get_Omega('nu') * h * h:.5f}; Om {bg.omega_m:.4f} vs {r.get_Omega('baryon') + r.get_Omega('cdm') + r.get_Omega('nu'):.4f}")
    t = time.time()
    for i in range(200):
        Background(0.0224, 0.12, 0.675, -0.9, -0.3).theta_star()
    print(f"background eval: {(time.time() - t) / 200 * 1e3:.2f} ms")
    for k in ["dr2", "dr1", "cmbc", "cmbrla", "PantheonPlus", "Union3", "DESY5", "DESDovekie"]:
        lk = like(k)
        bg = Background(0.02223, 0.1198, 0.6766, -1.0, 0.0)
        print(f"  {k:12s} n = {len(lk.data):5d}  chi2(Planck LCDM) = {lk.chi2(bg):9.2f}")


def cmd_reproduce(args):
    """DESI DR2 Table 6 rows that the compressed pipeline can reproduce."""
    rows = [
        (["dr2"], -4.7, 1.7), (["dr2", "cmbc"], -8.0, 2.4), (["dr2", "PantheonPlus"], -4.9, 1.7),
        (["dr2", "Union3"], -10.1, 2.7), (["dr2", "DESY5"], -13.6, 3.3),
        (["dr2", "cmbc", "PantheonPlus"], None, None), (["dr2", "cmbc", "Union3"], None, None),
        (["dr2", "cmbc", "DESY5"], None, None), (["dr2", "cmbc", "DESY5@z0.1"], None, None),
    ]
    out = []
    for keys, pub, pubsig in rows:
        t = time.time()
        r = fit_pair(keys, nstart=args.nstart)
        r["published_dchi2"] = pub
        r["published_sigma"] = pubsig
        out.append(r)
        x = r["x_ext"]
        print(f"{'+'.join(keys):40s} dchi2 = {-r['dchi2']:6.1f} ({r['sigma']:.1f}s)  published {pub} ({pubsig}s)"
              f"  w0 = {x['w0']:.3f} wa = {x['wa']:.3f} h = {x['h']:.3f} Om = {Background(**{k: x[k] for k in ('ob', 'oc', 'h')}).omega_m:.3f}"
              f"  [{time.time() - t:.0f} s]")
    with open(os.path.join(OUT, "reproduce.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)


def cmd_shapley(args):
    players = args.players.split(",")
    out = run_shapley(players, prior=args.prior, free_mnu=args.free_mnu, nstart=args.nstart,
                      workers=args.workers, tag=args.tag)
    print_shapley(out)


def cmd_bootstrap(args):
    base = json.load(open(os.path.join(OUT, f"{args.base}.json")))
    players = base["players"]
    grand = [r for r in base["coalitions"] if set(r["keys"]) == set(players)][0]
    x = grand["x_ext"] if args.truth == "w0wa" else grand["x_lcdm"]
    truth = dict(ob=x["ob"], oc=x["oc"], h=x["h"], w0=x.get("w0", -1.0), wa=x.get("wa", 0.0))
    x0 = {"+".join(r["keys"]): dict(lcdm=[r["x_lcdm"][n] for n in ("ob", "oc", "h")],
                                    w0wa=[r["x_ext"][n] for n in ("ob", "oc", "h", "w0", "wa")])
          for r in base["coalitions"] if r["keys"]}
    out = run_bootstrap(players, truth, nboot=args.nboot, prior=base["prior"], nstart=args.nstart,
                        workers=args.workers, tag=args.tag, x0=x0)
    phis = np.array([[r["shapley"][p] for p in players] for r in out["results"]])
    grands = np.array([r["grand"] for r in out["results"]])
    print(f"truth = {args.truth}: grand dchi2 {grands.mean():.2f} +- {grands.std():.2f} (data: {base['grand']:.2f})")
    for i, p in enumerate(players):
        sh = phis[:, i] / grands
        print(f"  {p:14s} phi = {phis[:, i].mean():6.2f} +- {phis[:, i].std():5.2f}   share {np.median(sh):.2f}"
              f" [{np.percentile(sh, 16):.2f}, {np.percentile(sh, 84):.2f}]   data: {base['shapley'][p]:.2f}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("validate")
    p = sub.add_parser("reproduce")
    p.add_argument("--nstart", type=int, default=8)
    p = sub.add_parser("shapley")
    p.add_argument("--players", default="dr2,cmbc,PantheonPlus,Union3,DESY5")
    p.add_argument("--prior", default="desi")
    p.add_argument("--free-mnu", action="store_true")
    p.add_argument("--nstart", type=int, default=6)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--tag", default="shapley")
    p = sub.add_parser("bootstrap")
    p.add_argument("--base", default="shapley")
    p.add_argument("--truth", default="w0wa")
    p.add_argument("--nboot", type=int, default=100)
    p.add_argument("--nstart", type=int, default=2)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--tag", default="bootstrap")
    p = sub.add_parser("speccurve")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--tag", default="speccurve")
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--analyse-only", action="store_true")
    args = ap.parse_args()
    {"validate": cmd_validate, "reproduce": cmd_reproduce, "shapley": cmd_shapley,
     "bootstrap": cmd_bootstrap, "speccurve": cmd_speccurve}[args.cmd](args)


if __name__ == "__main__":
    main()
