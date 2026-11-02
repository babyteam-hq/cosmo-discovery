"""Time delays from the ZTF scene-model light curves (R-127, cosmo-deep-7).

Light curves: per image, total flux = reference PSF-fit flux (the ZTF reference catalogue, calibrated with the
reference zero point) + the scene-model difference flux (ztf_lens_phot.py); same-night epochs averaged; magnitudes
with a 0.01 mag error floor; only the primary (field, ccd, quadrant) per filter, so one reference per image.

Estimator (PyCS-style curve shifting, simplified): for a pair (X, Y) and a trial delay tau, Y's epochs are shifted
back by tau and both curves are fitted jointly by one cubic B-spline (knots every `knot` days, a second-difference
ridge on the coefficients so season gaps do not run wild) for the intrinsic variability plus, for Y, a polynomial in
time (degree `mldeg`) for microlensing and the magnitude offset. chi2(tau) on a grid; the delay is the minimum.
Uncertainty: mock curves built from the best model at several true delays around the minimum, with the measured
noise (the formal errors plus the excess residual variance as exponentially correlated noise, 10 d), re-fitted the
same way; the error is the quadrature sum of the mock bias and scatter (the PyCS convention). Convention:
dt_XY = t_Y - t_X > 0 means Y lags X.

Usage (from research/cosmo-discovery):
    python ztf_lens_delays.py <target> X Y [--band zr] [--tmin -100 --tmax 100 --step 1] [--knot 30] [--mldeg 2]
                              [--nmock 40] [--plot]
Writes data/ztf_lens/<target>/delay_<X><Y>_<band>.json and (with --plot) a PNG beside it.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.interpolate import BSpline

from ztf_lens_targets import TARGETS

DATA = Path(__file__).resolve().parent / "data" / "ztf_lens"


def light_curves(target: str, band: str, floor: float = 0.01, snr_min: float = 2.0,
                 bin_days: float = 1.0) -> dict[str, pd.DataFrame]:
    """Per image: magnitudes from reference + difference flux, binned in flux over `bin_days` (1 = nightly) before
    the conversion so faint images are not truncated by a per-night S/N cut; primary reference set only."""
    t = TARGETS[target]
    ph = pd.read_csv(DATA / target / "phot.csv")
    ph = ph[(ph.filtercode == band) & (ph.infobits == 0) & (ph.seeing < 4.0) & (ph.chi2red < 3)]
    key = ph.groupby(["field", "ccdid", "qid"]).size().idxmax()
    ph = ph[(ph.field == key[0]) & (ph.ccdid == key[1]) & (ph.qid == key[2])]
    stem = f"{key[0]:06d}_{band}_c{key[1]:02d}_q{key[2]}"
    cat = fits.open(DATA / target / "ref" / f"{stem}_refpsfcat.fits")
    zp = float(cat[0].header["MAGZP"])
    c = cat[1].data
    out = {}
    scene = ref_scene(target, stem, t["images"], ph)  # scene fit on the reference cutout: the fallback reference flux
    # the mean difference-image PSF is narrower than the reference coadd's, so the scene magnitudes run ~0.3 mag
    # faint; the offset measured on the images the catalogue does hold corrects the fallback
    offs = []
    for name, (r0, d0) in t["images"].items():
        d = np.hypot((c["ra"] - r0) * np.cos(np.radians(d0)), c["dec"] - d0) * 3600
        i = int(np.argmin(d))
        if d[i] <= 1.5 and np.isfinite(scene.get(name, np.nan)):
            offs.append(scene[name] - (float(c["mag"][i]) + zp))
    scene_corr = float(np.median(offs)) if offs else 0.0
    for name, (r0, d0) in t["images"].items():
        d = np.hypot((c["ra"] - r0) * np.cos(np.radians(d0)), c["dec"] - d0) * 3600
        i = int(np.argmin(d))
        if d[i] > 1.5:  # no catalogue source (a faint image, or one the reference pipeline merged with a galaxy)
            mref = scene.get(name, np.nan) - scene_corr
            fref = 10 ** (-0.4 * (mref - 23.9))
            eref = np.nan
            ref_src = f"scene (corrected by {scene_corr:+.2f})"
        else:
            mref = float(c["mag"][i]) + zp
            fref = 10 ** (-0.4 * (mref - 23.9))
            eref = float(c["sigmag"][i])
            ref_src = "refpsfcat"
        f = ph[f"f_{name}"].to_numpy() + fref
        e = ph[f"e_{name}"].to_numpy()
        df = pd.DataFrame({"jd": ph.jd.to_numpy(), "f": f, "e": e, "maglim": ph.maglim.to_numpy()})
        df = df[np.isfinite(df.f) & np.isfinite(df.e) & (df.e > 0)]
        # nightly bins (inverse-variance mean in flux)
        df["night"] = np.floor((df.jd - 0.5) / bin_days).astype(int)
        g = df.groupby("night")
        w = 1 / df.e**2
        fb = g.apply(lambda x: np.average(x.f, weights=1 / x.e**2), include_groups=False)
        eb = g.apply(lambda x: 1 / np.sqrt((1 / x.e**2).sum()), include_groups=False)
        tb = g.apply(lambda x: np.average(x.jd, weights=1 / x.e**2), include_groups=False)
        lc = pd.DataFrame({"jd": tb, "f": fb, "e": eb}).reset_index(drop=True)
        lc = lc[(lc.f > 0) & (lc.f / lc.e > snr_min)]
        lc["mag"] = 23.9 - 2.5 * np.log10(lc.f)
        lc["err"] = np.sqrt((1.0857 * lc.e / lc.f) ** 2 + floor**2)
        lc.attrs = {"mref": mref, "fref": fref, "eref": eref, "refset": stem, "n_epochs": int(len(df)),
                    "ref_source": ref_src, "mref_scene": scene.get(name, np.nan)}
        out[name] = lc.reset_index(drop=True)
    return out




def ref_scene(target: str, stem: str, images: dict, ph: pd.DataFrame) -> dict[str, float]:
    """Magnitudes from a scene fit on the reference-image cutout (the band's mean difference-image PSF stamp as the
    PSF, one flux per image plus a constant background), with the reference zero point. The fallback when the
    reference PSF-fit catalogue holds no source at an image's position."""
    from astropy.wcs import WCS

    from ztf_lens_phot import place_psf

    rp = DATA / target / "ref" / f"{stem}_refimg.fits"
    if not rp.exists():
        return {}
    with fits.open(rp) as h:
        img = np.array(h[0].data, dtype=float)
        hd = h[0].header
    psfs = []
    for st in ph.sort_values("seeing").stem.head(40):
        pp = DATA / target / "cut" / f"{st}_psf.fits"
        if pp.exists():
            psfs.append(fits.getdata(pp).astype(float))
    if not psfs:
        return {}
    psf = np.mean(psfs, axis=0)
    psf /= psf.sum()
    wcs = WCS(hd)
    names = list(images)
    px = wcs.all_world2pix(np.array([images[n] for n in names]), 0)
    cols = [place_psf(psf, img.shape, x, y) for x, y in px] + [np.ones_like(img)]
    m = np.isfinite(img)
    A = np.stack([c[m] for c in cols], axis=1)
    p, *_ = np.linalg.lstsq(A, img[m], rcond=None)
    zp = float(hd["MAGZP"])
    return {n: (zp - 2.5 * np.log10(f) if f > 0 else np.nan) for n, f in zip(names, p[:-1])}


def design(t: np.ndarray, knots: np.ndarray, k: int = 3) -> np.ndarray:
    return BSpline.design_matrix(t, knots, k).toarray()


class Pair:
    """A pair of light curves with the spline basis fixed over the whole tau range (X's design matrix cached)."""

    def __init__(self, tx, mx, ex, ty, my, ey, tmin, tmax, knot, mldeg, ridge, mlknot=0.0, min_frac=0.25):
        self.tx, self.mx, self.ex, self.ty, self.my, self.ey = tx, mx, ex, ty, my, ey
        lo = min(tx.min(), ty.min() - tmax) - 1
        hi = max(tx.max(), ty.max() - tmin) + 1
        nk = max(int((hi - lo) / knot), 4)
        inner = np.linspace(lo, hi, nk + 1)
        self.knots = np.concatenate([[lo] * 3, inner, [hi] * 3])
        self.Bx = design(tx, self.knots)
        self.nb = self.Bx.shape[1]
        self.t0, self.span = 0.5 * (lo + hi), 0.5 * (hi - lo)
        if mlknot:  # microlensing as a slow spline in Y's own time (ridge-regularised), else a polynomial
            mlo, mhi = ty.min() - 1, ty.max() + 1
            mk = max(int((mhi - mlo) / mlknot), 1)
            mknots = np.concatenate([[mlo] * 3, np.linspace(mlo, mhi, mk + 1), [mhi] * 3])
            self.P = design(ty, mknots)
        else:
            self.P = np.stack([((ty - self.t0) / self.span) ** p for p in range(mldeg + 1)], axis=1)
        self.mldeg = self.P.shape[1] - 1
        mldeg = self.mldeg
        D = np.zeros((self.nb - 2, self.nb + mldeg + 1))
        for i in range(self.nb - 2):
            D[i, i:i + 3] = [1, -2, 1]
        self.Dr = np.sqrt(ridge) * D
        if mlknot and self.P.shape[1] > 2:
            Dm = np.zeros((self.P.shape[1] - 2, self.nb + mldeg + 1))
            for i in range(self.P.shape[1] - 2):
                Dm[i, self.nb + i:self.nb + i + 3] = [1, -2, 1]
            self.Dr = np.vstack([self.Dr, np.sqrt(ridge) * Dm])
        self.w = 1 / np.concatenate([ex, ey])
        self.nx = len(tx)
        self.knot = knot
        self.min_frac = min_frac

    def fit(self, tau, mx=None, my=None):
        """Joint spline + microlensing fit at delay tau. Returns the statistic (chi2 per point over the overlap of
        the two shifted curves, each widened by one knot spacing), the model on the observed epochs, the
        coefficients and (chi2 over all points, chi2 over the overlap, number of overlap points). Scoring only the
        overlap keeps a long trial delay from winning because little is left to compare."""
        mx = self.mx if mx is None else mx
        my = self.my if my is None else my
        tys = self.ty - tau
        By = design(np.clip(tys, self.knots[0], self.knots[-1]), self.knots)
        n = self.nx + len(self.ty)
        A = np.zeros((n, self.nb + self.mldeg + 1))
        A[: self.nx, : self.nb] = self.Bx
        A[self.nx:, : self.nb] = By
        A[self.nx:, self.nb:] = self.P
        y = np.concatenate([mx, my])
        Aw2 = np.vstack([A * self.w[:, None], self.Dr])
        yw2 = np.concatenate([y * self.w, np.zeros(self.Dr.shape[0])])
        c, *_ = np.linalg.lstsq(Aw2, yw2, rcond=None)
        model = A @ c
        r2 = ((y - model) * self.w) ** 2
        lo, hi = max(self.tx.min(), tys.min()) - self.knot, min(self.tx.max(), tys.max()) + self.knot
        t_all = np.concatenate([self.tx, tys])
        inov = (t_all >= lo) & (t_all <= hi)
        n_ov = int(inov.sum())
        chi2_ov = float(r2[inov].sum())
        # a trial delay must keep at least min_frac of the points in the overlap (and 20 of them): with a short
        # overlap the ridge-regularised spline fits each curve alone, and chi2 per point drops for no good reason
        stat = chi2_ov / n_ov if n_ov >= max(20, self.min_frac * n) else np.inf
        return stat, model, c, (float(r2.sum()), chi2_ov, n_ov)

    def scan(self, taus, mx=None, my=None):
        return np.array([self.fit(tau, mx, my)[0] for tau in taus])


def best_delay(taus, chi):
    i = int(np.argmin(chi))
    # parabolic refinement
    if 0 < i < len(taus) - 1:
        a, b, c = chi[i - 1], chi[i], chi[i + 1]
        den = a - 2 * b + c
        if den > 0:
            return float(taus[i] + 0.5 * (a - c) / den * (taus[1] - taus[0]))
    return float(taus[i])


def corr_noise(t, sigma, ell, rng):
    """Exponentially correlated Gaussian noise with std sigma and correlation length ell (days)."""
    n = len(t)
    if n == 0 or sigma <= 0:
        return np.zeros(n)
    order = np.argsort(t)
    ts = t[order]
    x = np.empty(n)
    x[0] = rng.normal()
    for i in range(1, n):
        r = np.exp(-abs(ts[i] - ts[i - 1]) / ell)
        x[i] = r * x[i - 1] + np.sqrt(1 - r * r) * rng.normal()
    out = np.empty(n)
    out[order] = x * sigma
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("target")
    ap.add_argument("X")
    ap.add_argument("Y")
    ap.add_argument("--band", default="zr")
    ap.add_argument("--tmin", type=float, default=-100)
    ap.add_argument("--tmax", type=float, default=100)
    ap.add_argument("--step", type=float, default=1.0)
    ap.add_argument("--knot", type=float, default=30)
    ap.add_argument("--mldeg", type=int, default=2)
    ap.add_argument("--mlknot", type=float, default=0.0, help="microlensing spline knot spacing (d); 0 = polynomial")
    ap.add_argument("--ridge", type=float, default=1.0)
    ap.add_argument("--nmock", type=int, default=40)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--bin", type=float, default=1.0, help="flux-binning window in days before the magnitude conversion")
    ap.add_argument("--snr", type=float, default=2.0, help="drop bins below this total-flux S/N")
    ap.add_argument("--inject", type=float, default=None, help="injection-recovery at this true delay (d) instead of a fit")
    ap.add_argument("--inject-tol", type=float, default=40.0)
    ap.add_argument("--mock-window", type=float, default=30, help="half-width (d) of the mocks' tau scan around the minimum")
    ap.add_argument("--min-frac", type=float, default=0.25, help="least fraction of all points inside the overlap")
    ap.add_argument("--min-overlap", type=float, default=150, help="days of overlap below which a tau is not scanned")
    a = ap.parse_args()
    lcs = light_curves(a.target, a.band, snr_min=a.snr, bin_days=a.bin)
    X, Y = lcs[a.X], lcs[a.Y]
    tx, mx, ex = X.jd.to_numpy(), X.mag.to_numpy(), X.err.to_numpy()
    ty, my, ey = Y.jd.to_numpy(), Y.mag.to_numpy(), Y.err.to_numpy()
    taus = np.arange(a.tmin, a.tmax + a.step / 2, a.step)
    # overlap check: the data must overlap once shifted
    ov = np.array([min(tx.max(), ty.max() - tau) - max(tx.min(), ty.min() - tau) for tau in taus])
    ok = ov > a.min_overlap
    if ok.sum() < 5:
        print(f"{a.target} {a.X}{a.Y} {a.band}: no tau in [{a.tmin},{a.tmax}] leaves {a.min_overlap} d of overlap "
              f"(baseline {max(tx.max(), ty.max()) - min(tx.min(), ty.min()):.0f} d)")
        return
    taus = taus[ok]
    pair = Pair(tx, mx, ex, ty, my, ey, taus[0], taus[-1], a.knot, a.mldeg, a.ridge, a.mlknot, a.min_frac)
    stat = pair.scan(taus)
    fin = np.isfinite(stat)
    if fin.sum() < 5:
        print(f"{a.target} {a.X}{a.Y} {a.band}: no tau keeps {a.min_frac:.0%} of the points in the overlap")
        return
    taus, stat = taus[fin], stat[fin]
    tau0 = best_delay(taus, stat)
    stat0, model0, _c, (chi2_all0, chi2_ov0, n_ov0) = pair.fit(tau0)
    chi = (stat - stat0) * n_ov0  # delta chi2 equivalent, on the overlap at the best delay
    chi0 = 0.0
    resid = np.concatenate([mx, my]) - model0
    err = np.concatenate([ex, ey])
    excess = np.sqrt(max(np.var(resid) - np.mean(err**2), 0.0))
    ndof = len(resid) - (int((max(tx.max(), ty.max() - tau0) - min(tx.min(), ty.min() - tau0)) / a.knot) + a.mldeg + 4)
    # mocks: true delays spread around tau0 over +-2 steps (PyCS: several true delays), refit on the same grid
    rng = np.random.default_rng(a.seed)
    if a.inject is not None:
        # injection-recovery: mocks from the model fitted at the injected delay (the data's own variability), with the
        # formal noise and the excess as correlated noise, scanned over the FULL grid; how often is the injected delay
        # found, and where does the estimator go when it is not?
        _, model_i, _, _ = pair.fit(a.inject)
        rec = []
        for k in range(a.nmock):
            mxm = model_i[: len(tx)] + rng.normal(0, ex) + corr_noise(tx, excess, 10.0, rng)
            mym = model_i[len(tx):] + rng.normal(0, ey) + corr_noise(ty, excess, 10.0, rng)
            sm = pair.scan(taus, mxm, mym)
            rec.append(best_delay(taus, sm))
        rec = np.array(rec)
        tol = a.inject_tol
        hit = np.abs(rec - a.inject) < tol
        res = {"target": a.target, "pair": f"{a.X}{a.Y}", "band": a.band, "inject": a.inject, "tol": tol,
               "nmock": a.nmock, "hit_fraction": float(hit.mean()), "recovered": rec.tolist(),
               "excess_mag": float(excess), "data_tau": tau0, "grid": [float(taus[0]), float(taus[-1]), a.step]}
        out = DATA / a.target / f"inject_{a.X}{a.Y}_{a.band}_{int(a.inject)}.json"
        out.write_text(json.dumps(res, indent=1))
        miss = np.sort(rec[~hit])
        print(f"{a.target} {a.X}{a.Y} {a.band}: injected {a.inject} d, found within {tol} d in {hit.sum()}/{a.nmock}; "
              f"hits {np.median(rec[hit]) if hit.any() else float('nan'):.1f} +- {np.std(rec[hit]) if hit.any() else float('nan'):.1f}; "
              f"misses at {np.round(miss).astype(int).tolist()}; data minimum {tau0:.1f}")
        return
    mocks = []
    width = max(2 * a.step, 0.02 * abs(tau0) if tau0 else 0)
    # the mocks are scanned on a window around tau0 (the full grid is for the data's profile and its other minima)
    half = max(a.mock_window, 4 * width)
    mtaus = taus[(taus >= tau0 - half) & (taus <= tau0 + half)]
    for k in range(a.nmock):
        tau_true = tau0 + rng.uniform(-width, width)
        _, model_t, _, _ = pair.fit(tau_true)
        # the fitted model at tau_true evaluated on the observed epochs (Y's microlensing included)
        mxm = model_t[: len(tx)] + rng.normal(0, ex) + corr_noise(tx, excess, 10.0, rng)
        mym = model_t[len(tx):] + rng.normal(0, ey) + corr_noise(ty, excess, 10.0, rng)
        chim = pair.scan(mtaus, mxm, mym)
        mocks.append((tau_true, best_delay(mtaus, chim)))
    mocks = np.array(mocks)
    d = mocks[:, 1] - mocks[:, 0]
    bias, scat = float(np.mean(d)), float(np.std(d))
    tot = float(np.hypot(bias, scat))
    # the chi2 profile's own interval: delta chi2 < 1 scaled by chi2red
    chired = stat0  # chi2 per point over the overlap at the best delay
    inside = taus[chi - chi0 < max(chired, 1.0)]
    # quoted error: the larger of the mock error and the half-width of the (chi2red-scaled) delta chi2 < 1 interval,
    # so a broad or multi-modal minimum is not reported with the mocks' local scatter alone
    halfw = 0.5 * (inside.max() - inside.min()) if len(inside) else 0.0
    quoted = float(max(tot, halfw))
    res = {
        "target": a.target, "pair": f"{a.X}{a.Y}", "band": a.band, "tau": tau0, "err": quoted, "err_mock": tot,
        "err_chi2_halfwidth": float(halfw), "bias": bias,
        "scatter": scat, "n_X": int(len(tx)), "n_Y": int(len(ty)), "chi2": chi2_all0, "ndof": int(ndof),
        "chi2_overlap": chi2_ov0, "n_overlap": n_ov0, "stat_min": stat0,
        "excess_mag": float(excess), "grid": [float(taus[0]), float(taus[-1]), a.step],
        "chi2_interval": [float(inside.min()), float(inside.max())] if len(inside) else None,
        "knot": a.knot, "mldeg": a.mldeg, "mlknot": a.mlknot, "bin": a.bin, "snr": a.snr, "min_frac": a.min_frac, "ridge": a.ridge, "nmock": a.nmock,
        "mref": {a.X: X.attrs["mref"], a.Y: Y.attrs["mref"]}, "refset": X.attrs["refset"],
        "ref_source": {a.X: X.attrs["ref_source"], a.Y: Y.attrs["ref_source"]},
        "mref_scene": {a.X: X.attrs["mref_scene"], a.Y: Y.attrs["mref_scene"]},
        "taus": taus.tolist(), "chi2_profile": chi.tolist(), "stat_profile": stat.tolist(),
        # secondary minima: local minima within delta chi2 < 25 of the best
        "other_minima": [float(taus[i]) for i in range(1, len(taus) - 1)
                         if chi[i] < chi[i - 1] and chi[i] < chi[i + 1] and 0 < chi[i] - chi0 < 25
                         and abs(taus[i] - tau0) > 3 * a.step],
    }
    out = DATA / a.target / f"delay_{a.X}{a.Y}_{a.band}.json"
    out.write_text(json.dumps(res, indent=1))
    print(f"{a.target} {a.X}{a.Y} {a.band}: tau = {tau0:.2f} +- {quoted:.2f} d (mock {tot:.2f}: bias {bias:+.2f}, scatter {scat:.2f}; "
          f"chi2 {chi2_all0:.0f}/{ndof}, overlap {chi2_ov0:.0f}/{n_ov0}, excess {excess:.3f} mag; N {len(tx)}/{len(ty)}; "
          f"chi2 interval {res['chi2_interval']}; other minima {res['other_minima'][:6]})")
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(2, 1, figsize=(11, 8))
        ax[0].errorbar(tx, mx, ex, fmt=".", ms=3, label=a.X)
        ml = pair.P @ _c[pair.nb:]
        ax[0].errorbar(ty - tau0, my - ml, ey, fmt=".", ms=3, label=f"{a.Y} shifted by {tau0:.1f} d, microlensing removed")
        ax[0].invert_yaxis()
        ax[0].set_xlabel("JD")
        ax[0].legend()
        ax[1].plot(taus, chi - chi0)
        ax[1].set_xlabel("tau (d)")
        ax[1].set_ylabel("delta chi2")
        ax[1].set_ylim(-5, min(np.max(chi - chi0), 200))
        fig.savefig(out.with_suffix(".png"), dpi=100)


if __name__ == "__main__":
    main()
