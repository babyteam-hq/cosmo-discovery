#!/usr/bin/env python
"""R-123 revision pass after the first finding check (`results-R123.md`, row 5).

(a) The specification-curve variance decomposition on a balanced sub-grid (4 SN samples x 2 cuts x 2 BAO x 3 CMB x
    2 models x 2 priors [DESI, narrow] x 2 neutrino = 384 specifications, all present), with the significance and,
    for CPL alone, the raw delta-chi2 as the response; main effects are orthogonal on a balanced grid.
(b) The expected effect of the SN z > 0.1 cut under the fitted w0wa signal: parametric realisations of
    DESI DR2 + CMB_c + DES-SN5YR (and Pantheon+) around the coalition's own w0wa best fit, each fitted with and without
    the cut, against the data's drop.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import desi_shapley as ds  # noqa: E402


def balanced_anova():
    res = json.load(open(os.path.join(ds.OUT, "speccurve.json")))
    keep = lambda s: s["sn"] != "none" and s["prior"] in ("desi", "narrow")
    rows = [r for r in res if "error" not in r and keep(r["spec"])]
    print(f"balanced sub-grid: {len(rows)} specifications")
    factors = dict(ds.SPEC_FACTORS)
    factors["sn"] = [x for x in factors["sn"] if x != "none"]
    factors["prior"] = ["desi", "narrow"]
    out = {}
    for label, sel, metric in [("sigma, both models", lambda s: True, "sigma_dchi2"),
                               ("delta-chi2, CPL only", lambda s: s["de"] == "w0wa", "dchi2"),
                               ("sigma, CPL only", lambda s: s["de"] == "w0wa", "sigma_dchi2")]:
        rr = [r for r in rows if sel(r["spec"])]
        y = np.array([r[metric] for r in rr])
        gm = y.mean()
        sst = ((y - gm) ** 2).sum()
        print(f"\n== {label}: n = {len(rr)}, mean {gm:.2f}, sd {y.std():.2f}, range {y.min():.2f}-{y.max():.2f}")
        shares = {}
        for f, levels in factors.items():
            if f == "de" and "CPL only" in label:
                continue
            ss = 0.0
            means = {}
            for lv in levels:
                m = np.array([r["spec"][f] == lv for r in rr])
                means[lv] = float(y[m].mean())
                ss += m.sum() * (y[m].mean() - gm) ** 2
            shares[f] = float(ss / sst)
        # the cut x sample interaction
        ss_int = 0.0
        for a in factors["sn"]:
            for b in factors["sncut"]:
                m = np.array([r["spec"]["sn"] == a and r["spec"]["sncut"] == b for r in rr])
                ma = np.array([r["spec"]["sn"] == a for r in rr])
                mb = np.array([r["spec"]["sncut"] == b for r in rr])
                ss_int += m.sum() * (y[m].mean() - y[ma].mean() - y[mb].mean() + gm) ** 2
        shares["sn x sncut"] = float(ss_int / sst)
        for f, s in sorted(shares.items(), key=lambda kv: -kv[1]):
            lv = ""
            if f in factors:
                lv = "   " + ", ".join(f"{k} {np.mean([r[metric] for r in rr if r['spec'][f] == k]):.2f}" for k in factors[f])
            print(f"  {f:10s} {100 * s:5.1f} %{lv}")
        print(f"  main effects {100 * sum(v for k, v in shares.items() if k in factors):.0f} %")
        out[label] = shares
    json.dump(out, open(os.path.join(ds.OUT, "speccurve_balanced_anova.json"), "w"), indent=1)


def cut_baseline(sample, nreal=12):
    keys = ["dr2", "cmbc", sample]
    data_full = ds.fit_pair(keys, nstart=3)
    data_cut = ds.fit_pair(["dr2", "cmbc", sample + "@z0.1"], nstart=3)
    print(f"\n== {sample}: data dchi2 {data_full['dchi2']:.2f} ({data_full['sigma']:.1f}s) -> cut {data_cut['dchi2']:.2f}"
          f" ({data_cut['sigma']:.1f}s)")
    x = data_full["x_ext"]
    truth = dict(ob=x["ob"], oc=x["oc"], h=x["h"], w0=x["w0"], wa=x["wa"])
    bg = ds.Background(**truth)
    full, cut = [], []
    x0 = dict(lcdm=[data_full["x_lcdm"][n] for n in ("ob", "oc", "h")],
              w0wa=[data_full["x_ext"][n] for n in ("ob", "oc", "h", "w0", "wa")])
    for i in range(nreal):
        rng = np.random.default_rng(5000 + i)
        sims = {k: ds.like(k).simulate(bg, rng) for k in keys}
        # the cut sample shares the realisation: same noise draw, rows removed
        sn_full = sims[sample]
        sn_cut = ds.like(sample + "@z0.1")
        import copy
        sn_cut = copy.copy(sn_cut)
        mask = sn_full.zhd > 0.1
        sn_cut.data = sn_full.data[mask]
        rf = ds.fit_pair(keys, nstart=1, likes=[sims["dr2"], sims["cmbc"], sn_full], x0=x0)
        rc = ds.fit_pair(keys, nstart=1, likes=[sims["dr2"], sims["cmbc"], sn_cut], x0=x0)
        full.append(rf["dchi2"])
        cut.append(rc["dchi2"])
        print(f"  real {i}: full {rf['dchi2']:.1f} cut {rc['dchi2']:.1f}", flush=True)
    full, cut = np.array(full), np.array(cut)
    ratio = cut / full
    print(f"  expected under the fitted w0wa: full {full.mean():.1f} +- {full.std():.1f}, cut {cut.mean():.1f} +- {cut.std():.1f},"
          f" ratio {ratio.mean():.2f} +- {ratio.std():.2f}; data ratio {data_cut['dchi2'] / data_full['dchi2']:.2f};"
          f" P(ratio <= data) = {(ratio <= data_cut['dchi2'] / data_full['dchi2']).mean():.2f}")
    return dict(sample=sample, data_full=data_full["dchi2"], data_cut=data_cut["dchi2"], truth=truth,
                full=full.tolist(), cut=cut.tolist())


if __name__ == "__main__":
    balanced_anova()
    out = [cut_baseline("DESY5"), cut_baseline("PantheonPlus"), cut_baseline("Union3")]
    json.dump(out, open(os.path.join(ds.OUT, "cut_baseline.json"), "w"), indent=1)
