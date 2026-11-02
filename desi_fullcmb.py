#!/usr/bin/env python
"""R-123: the full-CMB check of the headline coalitions.

Planck 2018 plik_lite TT,TE,EE + Commander low-l TT + SimAll low-l EE (+ the PR3 lensing reconstruction), through
cobaya's clik likelihoods (clipy) and CAMB (PPF dark energy, DESI's choice), combined with DESI DR2 BAO and the SN
samples from `desi_shapley.py`. Parameters (omega_b, omega_c, 100 theta_MC, tau, logA, ns, A_planck, w0, wa) are
whitened with the Planck 2018 w0wa+BAO covariance (planck_supp_data_and_covmats), then the profile chi2 is minimised
by L-BFGS-B (finite differences of 0.1 sigma) and polished with Nelder-Mead, from a warm start. One fit per call;
results as JSON in data/r123/out/fullcmb/.

    data/r123/venv/bin/python desi_fullcmb.py fit --coalition dr2,cmb,DESY5 --model w0wa [--no-lensing]
    data/r123/venv/bin/python desi_fullcmb.py table
"""
import argparse
import glob
import json
import os
import sys
import time

import numpy as np
from scipy import optimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import desi_shapley as ds  # noqa: E402

PK = os.path.join(ds.DATA, "packages")
OUTDIR = os.path.join(ds.OUT, "fullcmb")
os.makedirs(OUTDIR, exist_ok=True)
COVMAT = os.path.join(PK, "data", "planck_supp_data_and_covmats", "covmats",
                      "base_w_wa_plikHM_TTTEEE_lowl_lowE_BAO.covmat")

NAMES = ["ombh2", "omch2", "theta", "tau", "logA", "ns", "A_planck", "w0", "wa"]
COVNAMES = dict(ombh2="omegabh2", omch2="omegach2", theta="theta", tau="tau", logA="logA", ns="ns",
                A_planck="calPlanck", w0="w", wa="wa")
P0 = dict(ombh2=0.02237, omch2=0.1200, theta=1.04092, tau=0.0544, logA=3.044, ns=0.9649, A_planck=1.0, w0=-1.0, wa=0.0)
BOUNDS = dict(ombh2=(0.017, 0.027), omch2=(0.09, 0.16), theta=(1.03, 1.05), tau=(0.01, 0.2), logA=(2.8, 3.3),
              ns=(0.9, 1.03), A_planck=(0.98, 1.02), w0=(-3.0, 1.0), wa=(-3.0, 2.0))


def planck_cov(free):
    names = open(COVMAT).readline().split()[1:]
    C = np.loadtxt(COVMAT)
    idx = [names.index(COVNAMES[n]) for n in free]
    return C[np.ix_(idx, idx)]


def build_cmb_model(lensing=True):
    from cobaya.model import get_model

    likes = {"planck_2018_highl_plik.TTTEEE_lite": None, "planck_2018_lowl.TT_clik": None,
             "planck_2018_lowl.EE_clik": None}
    if lensing:
        likes["planck_2018_lensing.clik"] = None
    params = {"ombh2": {"prior": {"min": 0.01, "max": 0.03}}, "omch2": {"prior": {"min": 0.05, "max": 0.2}},
              "cosmomc_theta": {"prior": {"min": 0.0100, "max": 0.0108}}, "tau": {"prior": {"min": 0.005, "max": 0.3}},
              "logA": {"prior": {"min": 2.5, "max": 3.5}}, "As": {"value": "lambda logA: 1e-10*np.exp(logA)"},
              "ns": {"prior": {"min": 0.8, "max": 1.1}}, "w": {"prior": {"min": -3.0, "max": 1.0}},
              "wa": {"prior": {"min": -3.0, "max": 2.0}}, "mnu": 0.06,
              "A_planck": {"prior": {"min": 0.9, "max": 1.1}}, "H0": None}
    info = {"likelihood": likes,
            "theory": {"camb": {"extra_args": {"lens_potential_accuracy": 1, "dark_energy_model": "ppf",
                                               "num_massive_neutrinos": 1, "nnu": 3.044}}},
            "params": params, "packages_path": PK}
    m = get_model(info)
    return m


class FullProblem:
    def __init__(self, keys, model="w0wa", lensing=True):
        self.keys = keys
        self.lensing = lensing
        self.model = model
        self.cmb = build_cmb_model(lensing) if "cmb" in keys else None
        if self.cmb is not None:
            self.ih0 = list(self.cmb.parameterization.derived_params()).index("H0")
        self.others = [ds.like(k) for k in keys if k != "cmb"]
        self.free = ["ombh2", "omch2", "theta"]
        if self.cmb is not None:
            self.free += ["tau", "logA", "ns", "A_planck"]
        if model in ("w0wa", "wcdm"):
            self.free.append("w0")
        if model == "w0wa":
            self.free.append("wa")
        C = planck_cov(self.free)
        if self.cmb is None:  # BAO/SN-only coalition: a diagonal, wide scale
            C = np.diag(np.diag(C)) * 25
        self.L = np.linalg.cholesky(C)  # p = p_start + L y
        self.neval = 0
        self.best = (np.inf, None, None)
        self.pstart = dict(P0)

    def unpack(self, y):
        p = dict(self.pstart)
        dp = self.L @ np.asarray(y)
        for n, v in zip(self.free, dp):
            p[n] = self.pstart[n] + v
        return p

    def h0_of(self, p):
        """H0 for a theta-parametrised point when no CMB model is loaded (background-only coalitions)."""
        import camb
        pars = camb.set_params(cosmomc_theta=p["theta"] / 100, ombh2=p["ombh2"], omch2=p["omch2"], mnu=0.06,
                               nnu=3.044, tau=0.055, As=2.1e-9, ns=0.965, dark_energy_model="ppf",
                               w=p["w0"], wa=p["wa"])
        return pars.H0

    def chi2_parts(self, p):
        parts = {}
        if self.cmb is not None:
            vals = dict(ombh2=p["ombh2"], omch2=p["omch2"], cosmomc_theta=p["theta"] / 100, tau=p["tau"],
                        logA=p["logA"], ns=p["ns"], w=p["w0"], wa=p["wa"], A_planck=p["A_planck"])
            ll, der = self.cmb.loglikes(vals, return_derived=True)
            for name, l in zip(self.cmb.likelihood, np.atleast_1d(ll)):
                parts[name] = -2.0 * float(l)
            parts["A_planck prior"] = ((p["A_planck"] - 1.0) / 0.0025) ** 2
            h0 = float(der[self.ih0])
        else:
            h0 = self.h0_of(p)
        p["H0"] = h0
        if self.others:
            bg = ds.Background(p["ombh2"], p["omch2"], h0 / 100.0, p["w0"], p["wa"], 0.06)
            for lk in self.others:
                parts[lk.name] = lk.chi2(bg)
        if self.cmb is None:
            parts["BBN prior"] = ((p["ombh2"] - 0.02218) / 0.00055) ** 2
        return parts

    def chi2(self, y):
        p = self.unpack(y)
        pen = 0.0
        for n in self.free:
            lo, hi = BOUNDS[n]
            sig = np.sqrt((self.L @ self.L.T)[self.free.index(n), self.free.index(n)])
            if p[n] < lo:
                pen += 1e4 * ((lo - p[n]) / sig) ** 2
                p[n] = lo
            if p[n] > hi:
                pen += 1e4 * ((p[n] - hi) / sig) ** 2
                p[n] = hi
        if self.model == "w0wa" and p["w0"] + p["wa"] > 0:
            pen += 1e4 * (p["w0"] + p["wa"]) ** 2
        try:
            parts = self.chi2_parts(p)
            tot = sum(parts.values())
        except Exception as e:  # CAMB failure far from the best fit
            print("eval failed:", repr(e)[:200], {k: round(v, 5) for k, v in p.items()}, flush=True)
            tot, parts = 1e8, {}
        if not np.isfinite(tot):
            tot = 1e8
        self.neval += 1
        if tot + pen < self.best[0]:
            self.best = (tot + pen, p, parts)
        if self.neval % 50 == 0:
            print(f"  eval {self.neval}: chi2 {tot:.3f} best {self.best[0]:.3f}", flush=True)
        return tot + pen

    def y_from(self, p):
        dp = np.array([p[n] - self.pstart[n] for n in self.free])
        return np.linalg.solve(self.L, dp)

    def fit(self, start):
        self.pstart = dict(P0)
        self.pstart.update({k: v for k, v in start.items() if k in self.pstart})
        y0 = np.zeros(len(self.free))
        t = time.time()
        # 1.7-2.5 s per evaluation: cap each pass (the chi2 is flat to 0.01 after ~300 evaluations in practice)
        r = optimize.minimize(self.chi2, y0, method="L-BFGS-B",
                              options=dict(eps=0.1, maxiter=40, maxfun=450, ftol=1e-12, gtol=1e-4, maxcor=30))
        print(f"L-BFGS-B pass 1: {r.fun:.3f} after {self.neval} evals, {time.time() - t:.0f} s: {r.message}", flush=True)
        r2 = optimize.minimize(self.chi2, self.y_from(self.best[1]), method="L-BFGS-B",
                               options=dict(eps=0.03, maxiter=15, maxfun=200, ftol=1e-12, gtol=1e-4, maxcor=30))
        print(f"L-BFGS-B pass 2: {r2.fun:.3f} after {self.neval} evals, {time.time() - t:.0f} s: {r2.message}", flush=True)
        n = len(self.free)
        simplex = np.vstack([self.y_from(self.best[1])] * (n + 1))
        for i in range(n):
            simplex[i + 1, i] += 0.15
        r3 = optimize.minimize(self.chi2, simplex[0], method="Nelder-Mead",
                               options=dict(xatol=0.005, fatol=0.001, maxfev=250, initial_simplex=simplex,
                                            adaptive=True))
        print(f"Nelder-Mead polish: {r3.fun:.3f} after {self.neval} evals, {time.time() - t:.0f} s", flush=True)
        return self.best


def warm_start(keys, model):
    """Start from the compressed-pipeline best fit of the same coalition (the CMB in its (theta, omega_b, omega_bc)
    form); theta_MC at Planck's value."""
    ckeys = [("cmbc" if k == "cmb" else k) for k in keys]
    p = dict(P0)
    try:
        r = ds.fit_pair(ckeys, nstart=4)
        x = r["x_ext"] if model != "lcdm" else r["x_lcdm"]
        p.update(ombh2=x["ob"], omch2=x["oc"], w0=x.get("w0", -1.0), wa=x.get("wa", 0.0))
    except Exception as e:
        print("warm start failed:", e, flush=True)
    if "cmb" not in keys:
        p["ombh2"] = 0.02218
    return p


def cmd_fit(args):
    keys = args.coalition.split(",")
    tag = "+".join(keys) + ("" if args.lensing else "-nolens") + "_" + args.model
    out_path = os.path.join(OUTDIR, tag + ".json")
    if os.path.exists(out_path) and not args.force:
        print("exists:", out_path)
        return
    t0 = time.time()
    prob = FullProblem(keys, args.model, args.lensing)
    start = warm_start(keys, args.model)
    if args.model == "lcdm":
        start["w0"], start["wa"] = -1.0, 0.0
    print("start:", {k: round(v, 4) for k, v in start.items()}, flush=True)
    chi2, p, parts = prob.fit(start)
    res = dict(keys=keys, model=args.model, lensing=args.lensing, chi2=chi2, params=p, parts=parts,
               nevals=prob.neval, seconds=time.time() - t0)
    json.dump(res, open(out_path, "w"), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))


def cmd_table(args):
    rows = {}
    for f in sorted(glob.glob(os.path.join(OUTDIR, "*.json"))):
        r = json.load(open(f))
        key = "+".join(r["keys"]) + ("" if r["lensing"] else "-nolens")
        rows.setdefault(key, {})[r["model"]] = r
    for key, d in rows.items():
        if "lcdm" in d and "w0wa" in d:
            dc = d["lcdm"]["chi2"] - d["w0wa"]["chi2"]
            p = d["w0wa"]["params"]
            parts = {k: d["lcdm"]["parts"].get(k, 0) - d["w0wa"]["parts"].get(k, 0) for k in d["w0wa"]["parts"]}
            print(f"{key:35s} dchi2 = {dc:6.2f} ({ds.significance(dc, 2):.1f} sigma)  w0 = {p['w0']:.3f} wa = {p['wa']:.3f}"
                  f" H0 = {p['H0']:.2f}  chi2 {d['lcdm']['chi2']:.2f} -> {d['w0wa']['chi2']:.2f}; per-part gain: "
                  + ", ".join(f"{k}: {v:.1f}" for k, v in parts.items()))
        else:
            print(f"{key:35s} have: {list(d)}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("fit")
    p.add_argument("--coalition", required=True)
    p.add_argument("--model", default="w0wa")
    p.add_argument("--no-lensing", dest="lensing", action="store_false")
    p.add_argument("--force", action="store_true")
    sub.add_parser("table")
    args = ap.parse_args()
    {"fit": cmd_fit, "table": cmd_table}[args.cmd](args)


if __name__ == "__main__":
    main()
