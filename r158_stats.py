"""R-158 statistics from classified_*.csv: A per field and redshift bin (binomial and tile bootstrap), dilution
correction (eps from the chiral injections plus the armless bound), the signed PSF term, classifier B and agreement,
trends (R50/FWHM, PA, tile quadrant), nulls, the JADES per-galaxy test and the dipole fits. Writes stats.json and
prints markdown tables. Usage: r158_stats.py"""
import os, json, numpy as np, pandas as pd
from r158_common import DATA
TH = json.load(open(os.path.join(DATA, "thresholds_sims.json")))
BINS = [(0.3, 0.6), (0.6, 0.9), (0.9, 1.2), (0.3, 1.2)]
AX = {"H": (78.0, 47.0), "S": (71.0, 61.0), "G": (12.86, -27.13)}
CEN = {"COSMOS": (150.12, 2.20), "EDF-N": (269.73, 66.02), "EDF-F": (52.93, -28.09), "EDF-S": (61.24, -48.42), "GOODS-S": (53.12, -27.81)}
rng = np.random.default_rng(158)

def cosang(a, b):
    a, b = np.radians(a), np.radians(b)
    return np.sin(a[1]) * np.sin(b[1]) + np.cos(a[1]) * np.cos(b[1]) * np.cos(a[0] - b[0])

P4B = {k: json.load(open(os.path.join(DATA, f"psf4b_{k}.json"))) for k in ("cosmos", "euclid_disc") if os.path.exists(os.path.join(DATA, f"psf4b_{k}.json"))}

def measure(df, col, T, n_sample, simk=None):
    k = df[col] != 0; s = df[col][k].values; N = len(s)
    if N < 10: return None
    A = s.mean(); sb = np.sqrt(max(1 - A * A, 1e-9) / N)
    tiles = df.tile[k].values; ut = np.unique(tiles)
    if len(ut) > 5:
        idx = {t: np.where(tiles == t)[0] for t in ut}; bs = []
        for _ in range(500):
            pick = rng.choice(ut, len(ut)); ii = np.concatenate([idx[t] for t in pick]); bs.append(s[ii].mean())
        sboot = float(np.std(bs))
    else:
        sboot = sb
    ec = T["chiral_err"]; w = min(1.0, 0.05 * n_sample / N); eps = ec * (1 - w) + 0.5 * w
    d = 1 - 2 * eps
    # signed PSF term: the armless share of the classified (<= w) carries the achiral sign mean
    psf, psf_e = w * T["achiral_A"], w * T["achiral_A_err"]
    if col == "clsA" and simk in P4B:   # null 4b: PSF term on real galaxies = dA/2 (signed), added to the armless term
        psf += P4B[simk]["dA"] / 2; psf_e = np.hypot(psf_e, P4B[simk]["dA_err"] / 2)
    sig = max(sb, sboot)
    Ac = (A - psf) / d if d > 0 else np.nan
    sc = np.hypot(sig, psf_e) / d if d > 0 else np.nan
    return {"N": int(N), "n_cw": int((s > 0).sum()), "n_ccw": int((s < 0).sum()), "A": float(A), "sig_binom": float(sb),
            "sig_boot": float(sboot), "w": float(w), "eps": float(eps), "dil": float(d), "psf": float(psf),
            "A_corr": float(Ac), "sig_corr": float(sc), "not_constraining_floor": bool(d < 0.3)}

def verdict(A, s, claim):
    lo, hi = A - 3 * s, A + 3 * s
    if not (lo <= claim <= hi): return "excluded"
    if lo <= 0 <= hi: return "not constraining"
    return "consistent"

def main():
    out = {}; rows = []
    sets = {"cosmos": ("cosmos", "COSMOS"), "euclid_disc": ("euclid_disc", None)}
    for name, (simk, fld) in sets.items():
        f = os.path.join(DATA, f"classified_{name}.csv")
        if not os.path.exists(f): continue
        df = pd.read_csv(f)
        fields = [fld] if fld else ["EDF-N", "EDF-F", "EDF-S", "all-Euclid"]
        for F in fields:
            dF = df if F in ("COSMOS", "all-Euclid") else df[df.field == F]
            for lo, hi in BINS:
                d = dF[(dF.z > lo) & (dF.z < hi)]
                for col, key in (("clsA", "f"), ("clsB", "S")):
                    r = measure(d, col, TH[simk][key], len(d), simk)
                    if r is None: continue
                    r.update(dataset=name, field=F, zlo=lo, zhi=hi, clf="CNN" if col == "clsA" else "Fourier", n_sample=len(d))
                    rows.append(r)
            d = dF[(dF.z > 0.3) & (dF.z < 1.2)]
            both = (d.clsA != 0) & (d.clsB != 0)
            agree = float((d.clsA[both] == d.clsB[both]).mean()) if both.sum() else np.nan
            out[f"agree_{F}"] = {"n_both": int(both.sum()), "agree": agree}
            # trends (CNN, 0.3-1.2)
            k = d.clsA != 0; dd = d[k]
            q = np.quantile(dd.r50_fwhm, [0, .25, .5, .75, 1])
            out[f"trend_size_{F}"] = [[float(q[i]), float(q[i + 1]), float(dd.clsA[(dd.r50_fwhm >= q[i]) & (dd.r50_fwhm <= q[i + 1])].mean()),
                                       int(((dd.r50_fwhm >= q[i]) & (dd.r50_fwhm <= q[i + 1])).sum())] for i in range(4)]
            out[f"trend_pa_{F}"] = [[a, a + 45, float(dd.clsA[(dd.pa >= a) & (dd.pa < a + 45)].mean()), int(((dd.pa >= a) & (dd.pa < a + 45)).sum())] for a in (0, 45, 90, 135)]
            med = dd.groupby("tile")[["ra", "dec"]].transform("median")
            quad = (dd.ra.values > med.ra.values).astype(int) * 2 + (dd.dec.values > med.dec.values).astype(int)
            out[f"trend_quad_{F}"] = [[int(qq), float(dd.clsA.values[quad == qq].mean()), int((quad == qq).sum())] for qq in range(4)]
            out[f"mirror_{F}"] = float(np.max(np.abs(d.f + d.f_mirror)))
    # nulls
    for name in ("cosmos_ell", "cosmos_star", "euclid_ell", "euclid_stars"):
        f = os.path.join(DATA, f"classified_{name}.csv")
        if not os.path.exists(f): continue
        df = pd.read_csv(f)
        for col in ("clsA", "clsB"):
            k = df[col] != 0
            out[f"null_{name}_{col}"] = {"N": int(len(df)), "Nc": int(k.sum()), "frac": float(k.mean()),
                                         "A": float(df[col][k].mean()) if k.sum() else np.nan, "sig": float(1 / np.sqrt(max(k.sum(), 1)))}
    # JADES
    f = os.path.join(DATA, "classified_jades.csv")
    if os.path.exists(f):
        df = pd.read_csv(f); sh = np.where(df.shamir == "cw", 1, -1)
        for col in ("clsA", "clsB"):
            k = df[col] != 0
            ag = (df[col][k].values == sh[k.values])
            out[f"jades_{col}"] = {"N": int(len(df)), "Nc": int(k.sum()), "A_ours": float(df[col][k].mean()),
                                   "A_shamir_same": float(sh[k.values].mean()), "agree": float(ag.mean()), "n_agree": int(ag.sum()),
                                   "sign_all_f": float(np.sign(df.f).mean())}
        mz = (df.z > 0.3) & (df.z < 1.2)
        k = (df.clsA != 0) & mz
        out["jades_z03_12"] = {"N": int(mz.sum()), "Nc": int(k.sum()), "A_ours": float(df.clsA[k].mean()) if k.sum() else np.nan,
                               "A_shamir": float(sh[mz.values].mean())}
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(DATA, "stats_rows.csv"), index=False)
    # dipole fits, CNN, corrected, per z bin: Euclid fields only (primary), COSMOS check, all four with offset
    fits_ = {}
    for lo, hi in BINS:
        sub = R[(R.clf == "CNN") & (R.zlo == lo) & (R.zhi == hi) & R.field.isin(["EDF-N", "EDF-F", "EDF-S", "COSMOS"])]
        for ax, a in AX.items():
            c = np.array([cosang(CEN[F], a) for F in sub.field]); y = sub.A_corr.values; s = sub.sig_corr.values
            e = sub.field.str.startswith("EDF").values
            if e.sum() >= 2:
                wgt = c[e] ** 2 / s[e] ** 2; Ad = (c[e] * y[e] / s[e] ** 2).sum() / wgt.sum(); sAd = 1 / np.sqrt(wgt.sum())
                chi2 = float((((y[e] - Ad * c[e]) / s[e]) ** 2).sum())
            else:
                Ad = sAd = chi2 = np.nan
            # all four with HST-minus-Euclid offset
            Xd = np.c_[c, (~e).astype(float)]; Wt = np.diag(1 / s ** 2)
            try:
                cov = np.linalg.inv(Xd.T @ Wt @ Xd); bet = cov @ Xd.T @ Wt @ y
                Ad4, sAd4, off, soff = bet[0], np.sqrt(cov[0, 0]), bet[1], np.sqrt(cov[1, 1])
            except np.linalg.LinAlgError:
                Ad4 = sAd4 = off = soff = np.nan
            fits_[f"{ax}_{lo}_{hi}"] = {"Ad_euclid": float(Ad), "sig": float(sAd), "chi2_euclid": chi2,
                                        "Ad_all_offset": float(Ad4), "sig_all": float(sAd4), "offset": float(off), "sig_offset": float(soff)}
    out["dipole"] = fits_
    # free-axis look-elsewhere on the Euclid fields (1-deg grid of axes; shuffles of per-field A within errors)
    sub = R[(R.clf == "CNN") & (R.zlo == 0.3) & (R.zhi == 1.2) & R.field.isin(["EDF-N", "EDF-F", "EDF-S", "COSMOS"])]
    if len(sub) >= 3:
        ra_g, de_g = np.meshgrid(np.arange(0, 360, 2.0), np.arange(-89, 90, 2.0))
        axes = np.c_[ra_g.ravel(), de_g.ravel()]
        C_ = np.array([[cosang(CEN[F], a) for F in sub.field] for a in axes])
        def best(y, s):
            wgt = (C_ ** 2 / s ** 2).sum(1); Ad = (C_ * y / s ** 2).sum(1) / wgt
            z = np.abs(Ad) * np.sqrt(wgt); i = np.argmax(z); return float(z[i]), axes[i], float(Ad[i])
        z0, a0, Ad0 = best(sub.A_corr.values, sub.sig_corr.values)
        null = [best(rng.normal(0, sub.sig_corr.values), sub.sig_corr.values)[0] for _ in range(300)]
        out["free_axis"] = {"z_max": z0, "axis": a0.tolist(), "Ad": Ad0, "p_global": float(np.mean(np.array(null) >= z0))}
    json.dump(out, open(os.path.join(DATA, "stats.json"), "w"), indent=1, default=float)
    pd.set_option("display.width", 250)
    print(R[["dataset", "field", "zlo", "zhi", "clf", "n_sample", "N", "n_cw", "n_ccw", "A", "sig_binom", "sig_boot", "eps", "A_corr", "sig_corr"]].round(4).to_string())
    print(json.dumps({k: v for k, v in out.items()}, indent=1, default=float))

if __name__ == "__main__":
    main()
