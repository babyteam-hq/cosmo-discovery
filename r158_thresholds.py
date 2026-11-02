"""R-158: thresholds from the simulations (sign-free with respect to the targets). t = the |score| above which 5 % of
achiral simulated galaxies (dataset PSF, R50 and S/N) pass. Reports the chiral-injection pass rate and error at t, the
achiral sign mean (PSF term), and the sign mapping of each classifier to the model handedness (brightest third)."""
import json, sys, os, numpy as np
from r158_common import DATA
out = {}
for ds in sys.argv[1:]:
    d = json.load(open(os.path.join(DATA, f"sims_{ds}.json")))
    o = {}
    for k in ("f", "S"):
        a = np.array(d["psf"][k]); m = np.array(d["mirror_psf"][k])
        c = np.array(d["chiral"][k]); h = np.array(d["chiral"]["hand"]); snr = np.array(d["chiral"]["snr"])
        br = snr > np.percentile(snr, 67)
        agree = np.mean(np.sign(c[br]) == h[br]); sg = 1.0 if agree > 0.5 else -1.0
        t = float(np.quantile(np.abs(a), 0.95)); kc = np.abs(c) > t; ka = np.abs(a) > t; km = np.abs(m) > t
        o[k] = {"map_to_hand": sg, "bright_agreement": float(max(agree, 1 - agree)), "t": t,
                "achiral_pass": float(ka.mean()), "achiral_A": float(np.mean(np.sign(a[ka]))), "achiral_A_err": float(1 / np.sqrt(ka.sum())),
                "mirrorpsf_A": float(np.mean(np.sign(m[km]))), "achiral_mean": float(a.mean()), "achiral_mean_err": float(a.std() / np.sqrt(len(a))),
                "mirrorpsf_mean": float(m.mean()),
                "chiral_pass": float(kc.mean()), "chiral_err": float(np.mean(sg * np.sign(c[kc]) != h[kc]))}
        # injections by true handedness: subsample the chiral set to A_in, recover
        rng = np.random.default_rng(3); rec = {}
        for Ain in (0.05, 0.20):
            pos, neg = np.where(h > 0)[0], np.where(h < 0)[0]
            nn = min(len(pos), len(neg)); npos = int(round(nn * (1 + Ain) / (1 + Ain) )); nneg = int(round(nn * (1 - Ain) / (1 + Ain)))
            idx = np.r_[rng.choice(pos, npos, replace=False), rng.choice(neg, nneg, replace=False)]
            kk = idx[np.abs(c[idx]) > t]; Aobs = np.mean(sg * np.sign(c[kk])); e = o[k]["chiral_err"]
            rec[str(Ain)] = {"A_true": float(np.mean(h[idx])), "A_obs": float(Aobs), "A_corr": float(Aobs / (1 - 2 * e)), "err": float(1 / np.sqrt(len(kk)) / (1 - 2 * e))}
        o[k]["inject"] = rec
    out[ds] = o; print(ds, json.dumps(o, indent=1))
json.dump(out, open(os.path.join(DATA, "thresholds_sims.json"), "w"), indent=1)
