"""R-153: the ALP decay-line search, one dwarf (or fake dwarf) at a time, and the joint limit.

Chain per target (the fixed plan, results-R153.md "Method check"):
  entries  = alp_cells.offsub (same-exposure OFF at 1.0-1.75 deg, same ln-lambda bin) on images reduced with the fixed
             detector pattern subtracted (alp_reduce.py --pattern), cells with theta < THETA_ON;
  clean    = alp_cells.robust_clean (5 robust sigma about each cell's continuum; variance rescaled per detector);
  model    = per-cell continuum (Legendre degree 4 in ln lambda per detector)
             + A x line template (D_Omega(cell) - D_Omega(OFF)) x passband response at the line
             + beta x stellar nuisance (the cell's own fitted continuum x the same passband response: any narrow
               spectral feature proportional to the local starlight, of free amplitude, per mass);
  A        = weighted least squares per mass; its error and null distribution from 1000 random sign flips of whole
             exposures' residual contributions (the noise as it is, correlated within an exposure).
The per-exposure contributions are kept, so flips, jackknife and injections are sums over exposures.
A is in units of (g / 1e-11 GeV^-1)^2 for the nominal D-factor.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import alp_cells as C  # noqa: E402
import alp_model as M  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REDP = os.path.join(HERE, "data", "alp", "redp")
OUT = os.path.join(HERE, "data", "alp", "search")
THETA_ON = 0.6
DEG = 4
NFLIP = 1000


def mass_grid():
    """Masses whose line falls in 0.75-4.95 um: steps of 1/4 FWHM (R ~ 41 below 3.8 um, ~ 130 above)."""
    lam = [0.76]
    while lam[-1] < 4.95:
        lam.append(lam[-1] * np.exp(0.006 if lam[-1] < 3.82 else 0.002))
    lam = np.array(lam[:-1])
    return 2 * 1.23984198 / lam


def prepare(name, centre_name=None, theta_on=THETA_ON, red=REDP, inject=None, shape="fid", prof="gauss"):
    """Entries, templates and stellar nuisance for one target. inject = (masses_inj, A_inj): add A_inj x template
    (at mass m_inj) to the data before cleaning (the whole chain then runs on it)."""
    d = C.load(name, red)
    ent = C.offsub(d)
    ent = C.subset(ent, ent["theta"] < theta_on)
    masses = mass_grid()
    if inject is not None:
        mi, ai = inject
        Ti = C.templates(name, ent, np.atleast_1d(mi), shape, prof, centre_name)
        ent["y"] = ent["y"] + (Ti * np.atleast_1d(ai)[None, :]).sum(1)
    ent, summ = C.robust_clean(ent, deg=DEG)
    T = C.templates(name, ent, masses, shape, prof, centre_name)
    cont = ent["y"] - C.continuum_residuals(ent, DEG)
    v = M.DWARFS[centre_name or name]["v"] if (centre_name or name) in M.DWARFS else 0.0
    N = np.empty_like(T)
    for j, m in enumerate(masses):
        N[:, j] = np.nan_to_num(cont) * M.line_profile_nuinu(m, ent["lam"], ent["bw"], v, prof) * ent["bw"] / ent["lam"]
    return ent, masses, T, N, summ


def decompose(ent, T, N, deg=DEG):
    """Per-exposure contributions to the normal equations after projecting each cell off its continuum basis.
    Returns dict with TT, TN, NN (per mass) and Ty_e, Ny_e (exposure x mass), plus the exposure ids."""
    w = 1.0 / ent["var"]
    exps, einv = np.unique(ent["e"], return_inverse=True)
    ne, k = len(exps), T.shape[1]
    out = {key: np.zeros(k) for key in ("TT", "TN", "NN")}
    Ty_e = np.zeros((ne, k))
    Ny_e = np.zeros((ne, k))
    order = np.argsort(ent["cell"], kind="stable")
    cells, st = np.unique(ent["cell"][order], return_index=True)
    en = np.r_[st[1:], len(order)]
    for c, a, b in zip(cells, st, en):
        idx = order[a:b]
        if len(idx) < 3 * (deg + 1) + 5:
            continue
        B = C.basis(ent["det"][idx], ent["lam"][idx], deg)
        ww = w[idx]
        BW = B * ww[:, None]
        F = np.linalg.pinv(B.T @ BW)

        def proj(X):
            return X - B @ (F @ (BW.T @ X))
        y = proj(ent["y"][idx])
        Tc = proj(T[idx])
        Nc = proj(N[idx])
        out["TT"] += (Tc * Tc * ww[:, None]).sum(0)
        out["TN"] += (Tc * Nc * ww[:, None]).sum(0)
        out["NN"] += (Nc * Nc * ww[:, None]).sum(0)
        np.add.at(Ty_e, einv[idx], Tc * (ww * y)[:, None])
        np.add.at(Ny_e, einv[idx], Nc * (ww * y)[:, None])
    out["Ty_e"], out["Ny_e"], out["exps"] = Ty_e, Ny_e, exps
    return out


def solve(dec, signs=None, stellar=True):
    """A per mass (and per sign-flip realisation if signs is (n_real, n_exp))."""
    Ty_e, Ny_e = dec["Ty_e"], dec["Ny_e"]
    if signs is None:
        Ty, Ny = Ty_e.sum(0), Ny_e.sum(0)
    else:
        Ty, Ny = signs @ Ty_e, signs @ Ny_e
    TT, TN, NN = dec["TT"], dec["TN"], dec["NN"]
    with np.errstate(invalid="ignore", divide="ignore"):
        if not stellar:
            return Ty / TT
        det = TT * NN - TN ** 2
        return (NN * Ty - TN * Ny) / det


def white_sigma(dec, stellar=True):
    with np.errstate(invalid="ignore", divide="ignore"):
        if not stellar:
            return 1 / np.sqrt(dec["TT"])
        return np.sqrt(dec["NN"] / (dec["TT"] * dec["NN"] - dec["TN"] ** 2))


def run_target(name, centre_name=None, tag=None, seed=153, stellar=True, **kw):
    ent, masses, T, N, summ = prepare(name, centre_name, **kw)
    dec = decompose(ent, T, N)
    A = solve(dec, stellar=stellar)
    rng = np.random.default_rng(seed)
    signs = rng.choice([-1.0, 1.0], size=(NFLIP, len(dec["exps"])))
    Af = solve(dec, signs, stellar=stellar)
    sF = np.nanstd(Af, 0)
    sW = white_sigma(dec, stellar)
    os.makedirs(OUT, exist_ok=True)
    tag = tag or name
    np.savez(os.path.join(OUT, f"{tag}.npz"), masses=masses, A=A, sF=sF, sW=sW, Af=Af.astype(np.float32),
             nexp=len(dec["exps"]), nent=len(ent["y"]), summ=str(summ), TT=dec["TT"], TN=dec["TN"], NN=dec["NN"])
    return masses, A, sF, sW, Af, summ


def profile_limit(Ahat, sig, logd_sig, cl_dchi2=2.71, gmax=400.0, n=4001):
    """Joint 95% upper limit on G = (g / 1e-11)^2 for one mass: dwarfs d with Ahat_d +- sig_d measuring G x 10^delta_d,
    delta_d ~ N(0, logd_sig_d) (the D-factor's log10 uncertainty), profiled. Returns (G95, G_best, TS of G_best)."""
    Gs = np.linspace(0, gmax, n)
    dl = np.linspace(-4, 4, 161)
    chi = np.zeros(n)
    for a, s, sd in zip(Ahat, sig, logd_sig):
        if not np.isfinite(a) or not np.isfinite(s):
            continue
        delta = dl * sd
        pred = Gs[:, None] * 10 ** delta[None, :]
        c = (a - pred) ** 2 / s ** 2 + dl[None, :] ** 2
        chi += c.min(1)
    i0 = np.argmin(chi)
    above = np.where((np.arange(n) > i0) & (chi > chi[i0] + cl_dchi2))[0]
    G95 = Gs[above[0]] if len(above) else np.inf
    return G95, Gs[i0], chi[0] - chi[i0]


if __name__ == "__main__":
    name = sys.argv[1]
    centre = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] != "-" else None
    masses, A, sF, sW, Af, summ = run_target(name, centre)
    print(name, "clean", summ)
    for i in range(0, len(masses), 12):
        print("m=%.3f A=%8.3f sF=%7.3f sW=%7.3f ratio=%.1f g95exp=%.2e" % (masses[i], A[i], sF[i], sW[i], sF[i] / sW[i],
              1e-11 * np.sqrt(1.645 * sF[i])))
