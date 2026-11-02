"""R-153: the cell-based line fit for one dwarf (or fake dwarf).

Data (alp_reduce.py): per exposure e, per sky cell c (0.05 deg, gnomonic about the centre), per ln(lambda) bin b: N, sum of
IMAGE - ZODI, sum of VARIANCE, sum lambda, sum bandwidth.

1. Same-exposure OFF: in each (e, b), the mean over cells with 1.0 <= theta < 1.75 deg is subtracted from every cell's
   mean (it removes the exposure's offset and zodi-model error at that wavelength).
2. Per-cell continuum: every ON cell (theta < THETA_ON) has its own smooth spectrum, a Legendre polynomial of degree
   DEG in ln(lambda) per detector, fitted across all the exposures that saw it (each at a different wavelength, the LVF
   moving it). It removes fixed sky structure (cirrus, stars, galaxies, the dwarf's own smooth starlight) that the
   LVF would otherwise map into spectral structure in a single exposure.
3. Line: template t = (D_Omega(cell) - D_Omega(OFF)) x the pixel passband response at the line, per entry; the
   amplitude A (in units of g^2 = 1 GeV^-2 scaled, see fit()) by weighted least squares after projecting each cell's
   template and data off its continuum basis. A stellar-feature nuisance can be added (fit(..., star=...)).
"""
import glob
import os

import numpy as np

import alp_model as M

HERE = os.path.dirname(os.path.abspath(__file__))
RED = os.path.join(HERE, "data", "alp", "red")
L0, DL, NB = np.log(0.74), 0.002, 958
CELL, NC = 0.05, 80
C_UM = 2.99792458e14
DET_EDGES = {1: (0.7455, 1.1167), 2: (1.1031, 1.6476), 3: (1.6317, 2.4237), 4: (2.4070, 3.8255),
             5: (3.8020, 4.4234), 6: (4.4104, 5.0089)}


def cell_theta(sub=5):
    """Cell-centre angular radius and the sub-sampled positions (for D_Omega averaging), deg."""
    c = (np.arange(NC) - NC / 2 + 0.5) * CELL
    xi, eta = np.meshgrid(c, c)            # cell = cy * NC + cx
    off = (np.arange(sub) - (sub - 1) / 2) / sub * CELL
    ox, oy = np.meshgrid(off, off)
    sx = xi.ravel()[:, None] + ox.ravel()[None, :]
    sy = eta.ravel()[:, None] + oy.ravel()[None, :]
    return np.hypot(xi, eta).ravel(), np.hypot(sx, sy)


def domega_theta(name, theta, shape="fid"):
    """Local D_Omega (GeV cm^-2 sr^-1) at radius theta (deg), from dD(<theta)/dOmega."""
    th = np.maximum(np.asarray(theta, float), 1e-3)
    e1, e2 = th * 0.98, th * 1.02
    return (M.d_within(name, e2, shape) - M.d_within(name, e1, shape)) / (np.pi * (e2 ** 2 - e1 ** 2) * M.SQDEG_SR)


def load(name, red=RED):
    cf = os.path.join(red, f"{name}_cells.npz")
    files = sorted(glob.glob(os.path.join(red, name, "*.npz")))
    if os.path.exists(cf):
        d = dict(np.load(cf))
        if d["nexp"] == len(files):
            return d
    parts = {k: [] for k in ("e", "det", "cell", "lb", "n", "s", "v", "lam", "bw", "z")}
    for i, f in enumerate(files):
        z = np.load(f)
        k = len(z["n"])
        parts["e"].append(np.full(k, i, np.int32))
        parts["det"].append(np.full(k, int(z["det"]), np.int8))
        for key in ("cell", "lb", "n", "s", "v", "lam", "bw", "z"):
            parts[key].append(z[key])
    d = {k: np.concatenate(v) for k, v in parts.items()}
    d["nexp"] = len(files)
    np.savez(cf, **d)
    return d


def offsub(d, off=(1.0, 1.75), nmin=100):
    """Per (exposure, bin): subtract the OFF-annulus mean. Returns entries (cells < off[0]) with y, var in nu I_nu nW."""
    th_c, _ = cell_theta(1)
    th = th_c[d["cell"]]
    isoff = (th >= off[0]) & (th < off[1])
    key = d["e"].astype(np.int64) * NB + d["lb"]
    uk, inv = np.unique(key, return_inverse=True)
    no = np.bincount(inv, weights=np.where(isoff, d["n"], 0), minlength=len(uk))
    so = np.bincount(inv, weights=np.where(isoff, d["s"], 0), minlength=len(uk))
    vo = np.bincount(inv, weights=np.where(isoff, d["v"], 0), minlength=len(uk))
    mo = so / np.maximum(no, 1)
    keep = (th < off[0]) & (no[inv] >= nmin) & (d["n"] >= 4)
    n = d["n"][keep].astype(float)
    lam = d["lam"][keep] / n
    conv = 1e-20 * C_UM / lam * 1e9
    y = (d["s"][keep] / n - mo[inv][keep]) * conv
    var = (d["v"][keep] / n ** 2 + vo[inv][keep] / np.maximum(no[inv][keep], 1) ** 2) * conv ** 2
    return dict(e=d["e"][keep], det=d["det"][keep], cell=d["cell"][keep], lb=d["lb"][keep], lam=lam,
                bw=d["bw"][keep] / n, y=y, var=var, n=n, zr=(d["z"][keep] / n) * conv,
                doff_frac=None, theta=th[keep])


def basis(det, lam, deg):
    """Legendre polynomials in ln(lambda) per detector (block-diagonal across detectors)."""
    cols = []
    for k in range(1, 7):
        sel = det == k
        if not sel.any():
            continue
        lo, hi = np.log(DET_EDGES[k][0]), np.log(DET_EDGES[k][1])
        x = (2 * (np.log(lam) - lo) / (hi - lo) - 1).clip(-1.2, 1.2)
        for j in range(deg + 1):
            c = np.zeros(len(lam))
            c[sel] = np.polynomial.legendre.legval(x[sel], [0] * j + [1])
            cols.append(c)
    return np.array(cols).T


def templates(name, ent, masses, shape="fid", prof="gauss", centre_name=None):
    """Line template per entry and mass: (D_Omega(cell) - D_Omega(OFF)) x passband response, nu I_nu nW per
    (g / 1e-11 GeV^-1)^2. centre_name: take D from another dwarf (fake-dwarf null uses the real dwarf's profile)."""
    dn = centre_name or name
    _, sub = cell_theta(5)
    dom_cell = domega_theta(dn, sub, shape).mean(1)
    th_c, _ = cell_theta(1)
    offm = (th_c >= 1.0) & (th_c < 1.75)
    dom_off = dom_cell[offm].mean()
    dc = dom_cell[ent["cell"]] - dom_off
    v = M.DWARFS[dn]["v"] if name in M.DWARFS else 0.0
    T = np.empty((len(ent["y"]), len(masses)))
    for j, m in enumerate(masses):
        I1 = M.line_intensity(m, 1e-11, 1.0) * 1e9            # nW m^-2 sr^-1 per unit D_Omega at g = 1e-11
        T[:, j] = I1 * dc * M.line_profile_nuinu(m, ent["lam"], ent["bw"], v, prof)
    return T


def fit(ent, T, deg=4, wscale=None, extra=None, return_proj=False):
    """Line amplitude per mass (column of T) after projecting each cell off its continuum basis (and the optional
    `extra` global regressors, columns of an (n_entries, k) array, fitted jointly). Returns A, sigma_A (white weights),
    and chi2/dof of the continuum-only residuals."""
    w = 1.0 / ent["var"]
    if wscale is not None:
        w = w / wscale
    order = np.argsort(ent["cell"], kind="stable")
    cells, starts = np.unique(ent["cell"][order], return_index=True)
    ends = np.r_[starts[1:], len(order)]
    k = T.shape[1]
    ke = 0 if extra is None else extra.shape[1]
    G = np.zeros((k + ke, k + ke))   # we only need diagonal blocks per mass with the extra columns
    num = np.zeros(k)
    den = np.zeros(k)
    cross_te = np.zeros((k, ke))
    cross_ty = np.zeros(k)
    ee = np.zeros((ke, ke))
    ey = np.zeros(ke)
    chi2, dof = 0.0, 0
    Ty = np.zeros(k)
    TT = np.zeros(k)
    for c, a, b in zip(cells, starts, ends):
        idx = order[a:b]
        if len(idx) < 3 * (deg + 1) + 5:
            continue
        B = basis(ent["det"][idx], ent["lam"][idx], deg)
        ww = w[idx]
        BW = B * ww[:, None]
        try:
            F = np.linalg.inv(B.T @ BW + 1e-9 * np.eye(B.shape[1]))
        except np.linalg.LinAlgError:
            continue
        def proj(X):
            return X - B @ (F @ (BW.T @ X))
        y = proj(ent["y"][idx])
        Tc = proj(T[idx])
        chi2 += (y * y * ww).sum()
        dof += len(idx) - B.shape[1]
        Ty += (Tc * (ww * y)[:, None]).sum(0)
        TT += (Tc * Tc * ww[:, None]).sum(0)
        if ke:
            Ec = proj(extra[idx])
            cross_te += (Tc * ww[:, None]).T @ Ec
            ee += (Ec * ww[:, None]).T @ Ec
            ey += (Ec * (ww * y)[:, None]).sum(0)
    if ke:
        A = np.empty(k)
        sA = np.empty(k)
        for j in range(k):
            Mx = np.zeros((ke + 1, ke + 1))
            Mx[0, 0] = TT[j]
            Mx[0, 1:] = Mx[1:, 0] = cross_te[j]
            Mx[1:, 1:] = ee
            C = np.linalg.inv(Mx)
            p = C @ np.r_[Ty[j], ey]
            A[j], sA[j] = p[0], np.sqrt(C[0, 0])
    else:
        good = TT > 1e-12
        A = np.where(good, Ty / np.where(good, TT, 1), np.nan)
        sA = np.where(good, 1 / np.sqrt(np.where(good, TT, 1)), np.nan)
    return A, sA, chi2 / max(dof, 1)


def offsub_plane(d, off=(1.0, 1.75), nmin=100, order=1, n_iter=4):
    """As offsub, but the same-exposure background is a per-bin offset (from the OFF annulus) plus a polynomial of
    `order` in the gnomonic offsets (xi, eta), one per exposure, fitted to the OFF cells by alternating the offsets and
    the plane. The plane takes out the exposure's zodi-model gradient, which a one-sided OFF would leak into ON."""
    c = (np.arange(NC) - NC / 2 + 0.5) * CELL
    cx, cy = d["cell"] % NC, d["cell"] // NC
    xi, eta = c[cx], c[cy]
    th = np.hypot(xi, eta)
    isoff = (th >= off[0]) & (th < off[1])
    mean = d["s"] / np.maximum(d["n"], 1)
    terms = [xi, eta] if order >= 1 else []
    if order >= 2:
        terms += [xi * xi, xi * eta, eta * eta]
    X = np.array(terms).T if terms else np.zeros((len(xi), 0))
    key = d["e"].astype(np.int64) * NB + d["lb"]
    uk, inv = np.unique(key, return_inverse=True)
    no = np.bincount(inv, weights=np.where(isoff, d["n"], 0), minlength=len(uk))
    vo = np.bincount(inv, weights=np.where(isoff, d["v"], 0), minlength=len(uk))
    plane = np.zeros(len(mean))
    wo = np.where(isoff, d["n"], 0).astype(float)
    order_e = np.argsort(d["e"], kind="stable")
    es, st = np.unique(d["e"][order_e], return_index=True)
    en = np.r_[st[1:], len(order_e)]
    for _ in range(n_iter):
        a = np.bincount(inv, weights=wo * (mean - plane), minlength=len(uk)) / np.maximum(no, 1)
        r = mean - a[inv]
        if X.shape[1] == 0:
            break
        for e, s0, s1 in zip(es, st, en):
            idx = order_e[s0:s1]
            m = wo[idx] > 0
            if m.sum() < 50:
                continue
            Xi = X[idx][m]
            Wi = wo[idx][m]
            coef = np.linalg.lstsq(Xi * np.sqrt(Wi)[:, None], r[idx][m] * np.sqrt(Wi), rcond=None)[0]
            plane[idx] = X[idx] @ coef
    a = np.bincount(inv, weights=wo * (mean - plane), minlength=len(uk)) / np.maximum(no, 1)
    keep = (th < off[0]) & (no[inv] >= nmin) & (d["n"] >= 4)
    n = d["n"][keep].astype(float)
    lam = d["lam"][keep] / n
    conv = 1e-20 * C_UM / lam * 1e9
    y = (mean[keep] - a[inv][keep] - plane[keep]) * conv
    var = (d["v"][keep] / n ** 2 + vo[inv][keep] / np.maximum(no[inv][keep], 1) ** 2) * conv ** 2
    return dict(e=d["e"][keep], det=d["det"][keep], cell=d["cell"][keep], lb=d["lb"][keep], lam=lam,
                bw=d["bw"][keep] / n, y=y, var=var, n=n, zr=(d["z"][keep] / n) * conv, theta=th[keep])


def subset(ent, m):
    return {k: (v[m] if isinstance(v, np.ndarray) and len(v) == len(m) else v) for k, v in ent.items()}


def continuum_residuals(ent, deg=4, keep=None):
    """Residual of each entry about its cell's continuum (fitted on the `keep` entries only)."""
    keep = np.ones(len(ent["y"]), bool) if keep is None else keep
    w = 1.0 / ent["var"]
    res = np.full(len(w), np.nan)
    order = np.argsort(ent["cell"], kind="stable")
    cells, st = np.unique(ent["cell"][order], return_index=True)
    en = np.r_[st[1:], len(order)]
    for c, a, b in zip(cells, st, en):
        idx = order[a:b]
        k = keep[idx]
        if k.sum() < 3 * (deg + 1) + 5:
            continue
        B = basis(ent["det"][idx], ent["lam"][idx], deg)
        ww = w[idx] * k
        F = np.linalg.pinv(B.T @ (B * ww[:, None]))
        y = ent["y"][idx]
        res[idx] = y - B @ (F @ (B.T @ (ww * y)))
    return res


def robust_clean(ent, deg=4, nsig=5.0, n_iter=3):
    """Clip entries whose continuum residual exceeds nsig robust sigma (per detector, in units of the pipeline sigma),
    refitting the continuum without them; then rescale each detector's variance by its robust chi2 of the kept entries.
    The clip is symmetric and set on entries whose own noise is ~100x a line at the limit, so it cannot bias the line.
    Returns the cleaned entries and a summary per detector (kept fraction, variance scale)."""
    keep = np.ones(len(ent["y"]), bool)
    for _ in range(n_iter):
        res = continuum_residuals(ent, deg, keep)
        z = res / np.sqrt(ent["var"])
        newkeep = np.isfinite(z)
        for k in range(1, 7):
            m = (ent["det"] == k) & np.isfinite(z)
            if not m.any():
                continue
            zz = z[m & keep] if (m & keep).any() else z[m]
            s = 1.4826 * np.median(np.abs(zz - np.median(zz)))
            newkeep[m] &= np.abs(z[m]) < nsig * s
        keep = newkeep
    res = continuum_residuals(ent, deg, keep)
    z = res / np.sqrt(ent["var"])
    out = subset(ent, keep)
    out["var"] = out["var"].copy()
    summary = {}
    for k in range(1, 7):
        m = out["det"] == k
        if not m.any():
            continue
        sc = np.mean(z[keep][m] ** 2)
        out["var"][m] *= sc
        summary[k] = (round(float(np.mean(keep[ent["det"] == k])), 4), round(float(sc), 2))
    return out, summary
