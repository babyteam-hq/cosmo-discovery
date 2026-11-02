"""R-147: type-dependent photo-z isolation leakage and the ESD / g_obs bias it predicts (cosmo-deep-16b).

Emulates KiDS-bright photometric isolation (Brouwer et al. 2021, 3 Mpc/h70; Mistele et al. 2024, 4 Mpc/h70; no neighbour with
> 10 % of the lens stellar mass) on the DESI DR1 BGS group VAC, finds which photo-z-isolated lenses have a real > 10 % neighbour
(projected R, |dv| < 500 km/s, catalogue redshifts), and adds those neighbours' group halos (NFW, group-VAC M_HALO) to the lens
profile in Brouwer's 15 g_bar bins. Everything is fixed in results-R147.md ("Fixed before computing").

Writes one npz per run (seed, error model) with per-lens arrays; r147_table.py does the type splits, mass matching, bootstrap
and the offset statistic.

Usage: python -I r147_leak.py <data dir> <out dir> --seed 1 [--model lorentz|gauss|none] [--specnb]
Units: physical Mpc, M_sun (h70), H0 = 70, Om = 0.3.
"""
import argparse, os, time
import numpy as np, fitsio
from scipy.spatial import cKDTree
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import RegularGridInterpolator

C = 299792.458; H0, OM = 70.0, 0.3
G_SI = 6.674e-11; MSUN = 1.989e30; MPC = 3.0857e22
SIGMA_UNIT = MSUN / (MPC / 1e6) ** 2           # kg m^-2 per (M_sun pc^-2)
RHO_CRIT0 = 2.775e11 * 0.7 ** 2                 # M_sun Mpc^-3 (h70 = 1)
_zg = np.linspace(0, 1.5, 3001)
_E = np.sqrt(OM * (1 + _zg) ** 3 + 1 - OM)
_dc = cumulative_trapezoid(C / H0 / _E, _zg, initial=0)
def dcom(z): return np.interp(z, _zg, _dc)
def dang(z): return dcom(z) / (1 + z)
def ez(z): return np.sqrt(OM * (1 + z) ** 3 + 1 - OM)

GBAR_EDGES = np.logspace(-15, np.log10(5e-12), 16)
GBAR = np.sqrt(GBAR_EDGES[1:] * GBAR_EDGES[:-1])                   # bin centres, m s^-2
RISOL = (3.0, 4.0)
RFRAC = (0.3, 0.5, 1.0, 2.0, 3.0, 4.0)


def g_m16(gbar, gdag=1.2e-10):
    return gbar / (1 - np.exp(-np.sqrt(gbar / gdag)))


DSIG_OWN = g_m16(GBAR) / (4 * G_SI) / SIGMA_UNIT                     # M_sun pc^-2, the shared own profile


# ---- NFW, untruncated, off-centre ESD table in units of rho_s r_s -------------------------------------------------------
def nfw_f(x):
    x = np.asarray(x, float); out = np.empty_like(x)
    lo = x < 1 - 1e-6; hi = x > 1 + 1e-6; mid = ~(lo | hi)
    a = x[lo]; out[lo] = (1 - 2 / np.sqrt(1 - a ** 2) * np.arctanh(np.sqrt((1 - a) / (1 + a)))) / (a ** 2 - 1)
    b = x[hi]; out[hi] = (1 - 2 / np.sqrt(b ** 2 - 1) * np.arctan(np.sqrt((b - 1) / (1 + b)))) / (b ** 2 - 1)
    out[mid] = 1 / 3
    return 2 * out                                                    # Sigma / (rho_s r_s)


def build_table():
    xs = np.logspace(-4, 2.7, 300); ys = np.logspace(-3, 2.7, 200)
    th = (np.arange(1024) + 0.5) / 1024 * np.pi
    G = np.empty((xs.size, ys.size))
    for j, y in enumerate(ys):
        r = np.sqrt(xs[:, None] ** 2 + y ** 2 + 2 * xs[:, None] * y * np.cos(th)[None, :])
        ring = nfw_f(np.maximum(r, 1e-7)).mean(axis=1)
        inner = ring[0] * xs[0] ** 2 / 2
        mean_in = 2 / xs ** 2 * (inner + cumulative_trapezoid(ring * xs ** 2, np.log(xs), initial=0))
        G[:, j] = mean_in - ring
    # centred, analytic mean within x for the check and for y -> 0
    return xs, ys, G


class NFW:
    def __init__(self):
        xs, ys, G = build_table()
        self.ip = RegularGridInterpolator((np.log(xs), np.log(ys)), G, bounds_error=False, fill_value=None)
        self.xs, self.ys = xs, ys
        # centred: mean Sigma within x minus Sigma(x)
        xx = np.logspace(-6, 3, 4000); s = nfw_f(xx)
        m = 2 / xx ** 2 * (s[0] * xx[0] ** 2 / 2 + cumulative_trapezoid(s * xx ** 2, np.log(xx), initial=0))
        self.cen = lambda x: np.interp(np.log(x), np.log(xx), m - s)

    @staticmethod
    def params(m200m, z):
        """Duffy et al. 2008 (full sample, M_200m): c = 10.14 (M / 2e12 h^-1)^-0.081 (1+z)^-1.01; physical r_s, rho_s."""
        c = 10.14 * (m200m * 0.7 / 2e12) ** -0.081 * (1 + z) ** -1.01
        rhom = OM * RHO_CRIT0 * (1 + z) ** 3
        r200 = (3 * m200m / (4 * np.pi * 200 * rhom)) ** (1 / 3)
        rs = r200 / c
        rhos = 200 / 3 * rhom * c ** 3 / (np.log(1 + c) - c / (1 + c))
        return rs, rhos

    def dsig_off(self, R, m, z, d):
        """ESD (M_sun pc^-2) about a point at projected distance d (Mpc) from an NFW halo of mass m, at radii R (Mpc)."""
        rs, rhos = self.params(m, z)
        x = R / rs; y = np.maximum(d / rs, self.ys[0])
        x = np.clip(x, self.xs[0], self.xs[-1]); y = np.minimum(y, self.ys[-1])
        x, y = np.broadcast_arrays(x, y)
        return rhos * rs * self.ip(np.c_[np.log(x).ravel(), np.log(y).ravel()]).reshape(x.shape) / 1e12

    def dsig_cen(self, R, m, z):
        rs, rhos = self.params(m, z)
        return rhos * rs * self.cen(np.clip(R / rs, 1e-6, 999)) / 1e12


def sigcrit_inv2(zl):
    """Relative Sigma_crit^-2 for a KiDS-like source n(z) ~ z^2 exp[-(z/0.5)^1.5], z_s > z_l + 0.2."""
    zs = np.linspace(0.01, 3, 600); nz = zs ** 2 * np.exp(-(zs / 0.5) ** 1.5)
    out = np.empty_like(zl)
    zgrid = np.linspace(0.0, 0.6, 121); val = np.empty_like(zgrid)
    for k, z in enumerate(zgrid):
        s = zs > z + 0.2
        dls = (dcom(zs[s]) - dcom(z)) / (1 + zs[s]); ds = dang(zs[s])
        val[k] = dang(z) * np.trapezoid(nz[s] * dls / ds, zs[s]) / np.trapezoid(nz, zs)
    return np.interp(zl, zgrid, val) ** 2


def unitvec(ra, dec):
    th = np.radians(90 - dec); ph = np.radians(ra)
    return np.c_[np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]


def pairs(tree, rmax, idx, X, keep):
    """(i in idx, j) with |X_i - X_j| < rmax that pass keep(i, j), in chunks; returns i, j arrays (global indices)."""
    I, J = [], []
    for s in range(0, idx.size, 20000):
        ii = idx[s:s + 20000]
        lists = tree.query_ball_point(X[ii], rmax, workers=4)
        n = np.fromiter((len(l) for l in lists), int, ii.size)
        i = np.repeat(ii, n); j = np.fromiter((j for l in lists for j in l), int, n.sum())
        k = keep(i, j); I.append(i[k]); J.append(j[k])
    return np.concatenate(I), np.concatenate(J)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('data'); ap.add_argument('out')
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--model', default='lorentz')
    ap.add_argument('--specnb', action='store_true'); ap.add_argument('--nmax', type=int, default=0)
    a = ap.parse_args(); rng = np.random.default_rng(a.seed); t0 = time.time()
    os.makedirs(a.out, exist_ok=True)
    g = fitsio.read(os.path.join(a.data, 'bgs-groups-dr1-v1.0.fits'),
                    columns=['TARGETID', 'RA', 'DEC', 'Z', 'LOGMSTAR', 'Z_ASSIGNED_FLAG', 'QUIESCENT', 'G_R', 'APP_MAG_R',
                             'IGRP', 'IS_SAT', 'M_HALO', 'N_SAT'])
    ok = np.isfinite(g['LOGMSTAR']) & np.isfinite(g['Z']) & (g['Z'] > 0.002)
    if a.specnb: ok &= (g['Z_ASSIGNED_FLAG'] == 0) | (g['Z_ASSIGNED_FLAG'] == -1)
    g = g[ok]; N = g.size
    z = g['Z']; logm = g['LOGMSTAR'] + 2 * np.log10(1 / 0.7)
    spec = (g['Z_ASSIGNED_FLAG'] == 0) | (g['Z_ASSIGNED_FLAG'] == -1)
    u = unitvec(g['RA'], g['DEC'])
    # ---- emulated photo-z view
    if a.model == 'lorentz':
        nu = 2 * 2.613 - 1; sc = 0.0149 * np.sqrt(2 * 2.613 / nu)
        dz = sc * rng.standard_t(nu, N)
    elif a.model == 'gauss':
        dz = 0.02 * rng.standard_normal(N)
    else:
        dz = np.zeros(N)
    zp = np.clip(z + dz * (1 + z), 0.002, 1.4)
    dl, dlp = dcom(z) * (1 + z), dcom(zp) * (1 + zp)
    logmp = logm + 2 * np.log10(dlp / dl) + (0.12 * rng.standard_normal(N) if a.model != 'none' else 0)
    # ---- lens candidates
    cand = np.where(spec & (z > 0.05) & (z < 0.3) & (zp > 0.1) & (zp < 0.5) & (logmp > 8.5) & (logmp < 11.0))[0]
    if a.nmax: cand = np.sort(rng.choice(cand, a.nmax, replace=False))
    print(f'{N} galaxies, {cand.size} lens candidates ({time.time() - t0:.0f}s)', flush=True)
    # ---- photo-z isolation: 3D comoving positions from photo-z, physical distance = comoving / (1 + zp_lens)
    Xp = u * dcom(zp)[:, None]
    tree_p = cKDTree(Xp)
    I, J = pairs(tree_p, max(RISOL) * 1.5, cand, Xp,
                 lambda i, j: (i != j) & (logmp[j] > logmp[i] - 1)
                 & (np.linalg.norm(Xp[i] - Xp[j], axis=1) < max(RISOL) * (1 + zp[i])))
    d3 = np.linalg.norm(Xp[I] - Xp[J], axis=1) / (1 + zp[I])
    pos = np.searchsorted(cand, I)                                   # cand is sorted
    iso = {}
    for R in RISOL:
        bad = np.zeros(cand.size, bool); bad[pos[d3 < R]] = True; iso[R] = ~bad
    print(f'photo-z isolated: ' + ', '.join(f'{R} Mpc {iso[R].mean():.3f}' for R in RISOL) + f' ({time.time() - t0:.0f}s)',
          flush=True)
    del I, J, d3, pos, tree_p, Xp
    # ---- truth view for the photo-z-isolated lenses at either radius
    lens = cand[iso[RISOL[0]] | iso[RISOL[1]]]
    X = u * dcom(z)[:, None]; tree = cKDTree(X)
    rmax = np.sqrt((max(RISOL) * 1.31) ** 2 + (500 * 1.31 / H0 / 0.9) ** 2) + 0.2
    def tkeep(i, j):
        k = (i != j) & (logm[j] > logm[i] - 1) & (C * np.abs(z[j] - z[i]) / (1 + z[i]) < 500)
        rp = np.arccos(np.clip(np.sum(u[i] * u[j], axis=1), -1, 1)) * dang(z[i])
        return k & (rp < max(RISOL))
    I, J = pairs(tree, rmax, lens, X, tkeep)
    cosang = np.clip(np.sum(u[I] * u[J], axis=1), -1, 1)
    rp = np.arccos(cosang) * dang(z[I]); dv = C * np.abs(z[J] - z[I]) / (1 + z[I])
    keep = (dv < 500) & (rp < max(RISOL)); I, J, rp = I[keep], J[keep], rp[keep]
    lpos = np.searchsorted(lens, I)
    rpmin = np.full(lens.size, np.inf); np.minimum.at(rpmin, lpos, rp)
    print(f'truth pairs {I.size} ({time.time() - t0:.0f}s)', flush=True)
    # ---- halos: group table, draws
    ig = g['IGRP']; ug, ginv = np.unique(ig, return_inverse=True)
    cen = np.full(ug.size, -1); cs = np.where(~g['IS_SAT'])[0]; cen[ginv[cs]] = cs
    # groups whose central was cut (spec-only variant): use the most massive member present
    miss = np.where(cen < 0)[0]
    if miss.size:
        o = np.lexsort((-logm, ginv)); first = np.r_[0, np.where(np.diff(ginv[o]))[0] + 1]
        best = o[first]; cen[ginv[best][np.isin(ginv[best], miss)]] = best[np.isin(ginv[best], miss)]
    lmh = np.log10(g['M_HALO'][cen] / 0.7)                           # per group, h70
    shift = rng.normal(0, 0.15); lmh_draw = lmh + rng.normal(0, 0.2, ug.size) + shift
    # isolated-central halo mass vs stellar mass (catalogue medians, 0.1 dex)
    iso_c = (~g['IS_SAT']) & (g['N_SAT'] == 0)
    mb = np.arange(8.0, 12.01, 0.1); mid = (mb[1:] + mb[:-1]) / 2
    med = np.array([np.median(np.log10(g['M_HALO'][iso_c & (logm >= lo) & (logm < hi)] / 0.7))
                    if np.sum(iso_c & (logm >= lo) & (logm < hi)) > 20 else np.nan for lo, hi in zip(mb[:-1], mb[1:])])
    good = np.isfinite(med); miso = lambda lm: np.interp(lm, mid[good], med[good]) + shift
    nfw = NFW()
    # per lens: radii of the g_bar bins from the photo-z mass (what KiDS would use)
    mp = 10 ** logmp[lens]; mgal = mp * (1 + 10 ** (-0.69 * logmp[lens] + 6.63))
    Rk = np.sqrt(G_SI * mgal[:, None] * MSUN / GBAR[None, :]) / MPC  # physical Mpc, (nlens, 15)
    w = (Rk / dang(z[lens])[:, None]) ** 2 * sigcrit_inv2(z[lens])[:, None]
    out = dict(lens=lens, z=z[lens], zp=zp[lens], logm=logm[lens], logmp=logmp[lens], spec=spec[lens],
               q=g['QUIESCENT'][lens], gr=g['G_R'][lens], mag=g['APP_MAG_R'][lens], tid=g['TARGETID'][lens],
               rpmin=rpmin, w=w, Rk=Rk, model=a.model, seed=a.seed, shift=shift)
    for R in RISOL:
        sel = rp < R
        li, jj = lpos[sel], J[sel]
        hh = ginv[jj]
        key = np.unique(li.astype(np.int64) * (ug.size + 1) + hh)
        li, hh = key // (ug.size + 1), key % (ug.size + 1)
        L = lens[li]; own = ginv[L] == hh; sat = g['IS_SAT'][L]
        cidx = cen[hh]
        d = np.arccos(np.clip(np.sum(u[L] * u[cidx], axis=1), -1, 1)) * dang(z[L])
        m = 10 ** lmh_draw[hh]; zz = z[L]
        R2 = Rk[li]
        dsig = np.zeros_like(R2); comp = np.zeros(li.size, int)
        a_ = ~own; b_ = own & sat; c_ = own & ~sat
        for s_, k in ((a_, 0), (b_, 1)):
            if s_.any():
                dsig[s_] = nfw.dsig_off(R2[s_], m[s_][:, None], zz[s_][:, None], d[s_][:, None]); comp[s_] = k
        if c_.any():
            m0 = 10 ** miso(logm[L[c_]])
            dsig[c_] = np.maximum(nfw.dsig_cen(R2[c_], m[c_][:, None], zz[c_][:, None])
                                  - nfw.dsig_cen(R2[c_], m0[:, None], zz[c_][:, None]), 0); comp[c_] = 2
        tot = np.zeros((lens.size, GBAR.size)); parts = np.zeros((3, lens.size, GBAR.size))
        np.add.at(tot, li, dsig)
        for k in range(3): np.add.at(parts[k], li[comp == k], dsig[comp == k])
        isoR = np.zeros(lens.size, bool); isoR[np.searchsorted(lens, cand[iso[R]])] = True
        out[f'iso{R:g}'] = isoR; out[f'ds{R:g}'] = tot; out[f'parts{R:g}'] = parts
        print(f'R_isol {R}: {li.size} lens-halo pairs; components a/b/c {np.bincount(comp, minlength=3)}', flush=True)
    out['cand_n'] = cand.size
    for R in RISOL: out[f'cand_iso{R:g}'] = iso[R].sum()
    fn = os.path.join(a.out, f'r147_{a.model}{"_spec" if a.specnb else ""}_s{a.seed}.npz')
    np.savez_compressed(fn, **out); print(f'wrote {fn} ({time.time() - t0:.0f}s)', flush=True)


if __name__ == '__main__':
    main()
