"""R-158: thresholds and error rates on the held-out GZ1 set, before any target galaxy is classified.
For each dataset, the held-out raw stamps are degraded with Gaussian noise so that prep's S/N follows that dataset's
S/N distribution (S/N from prep, a mirror-invariant quantity). t_A (CNN) is the smallest |f| with accuracy >= 0.90 above
it; t_B (Fourier) the smallest |S| with accuracy >= 0.85. eps = 1 - accuracy above threshold. Also fixes B's sign.
Usage: r158_validate.py cosmos euclid jades"""
import os, sys, json, numpy as np, torch
from r158_common import DATA
import r158_classify as C
torch.set_num_threads(4)
d = np.load(os.path.join(DATA, "gz1", "train.npz")); p = np.load(os.path.join(DATA, "gz1", "prepped.npz"))
vi = np.load(os.path.join(DATA, "val_idx.npz"))["vi"]
raw = d["x"][p["keep"][vi]]; y = (p["y"][vi] > 0)
net = C.make_net(); net.load_state_dict(torch.load(os.path.join(DATA, "net.pt"))); net.eval()
def run(stamps, rng=None, snr_t=None):
    X, ok, snr = [], [], []
    for j, s in enumerate(stamps):
        s = s.astype(np.float64)
        if snr_t is not None:
            _, i0 = C.prep(s)
            if "snr" not in i0: ok.append(False); X.append(np.zeros((64, 64), np.float32)); snr.append(np.nan); continue
            b = s[np.hypot(*np.mgrid[-47.5:48.5, -47.5:48.5]) > 0.42 * 96]
            sig = 1.4826 * np.median(np.abs(b - np.median(b)))
            new = sig * i0["snr"] / max(snr_t[j], 1e-3)
            if new > sig: s = s + rng.normal(0, np.sqrt(new ** 2 - sig ** 2), s.shape)
        x, info = C.prep(s)
        ok.append(x is not None); X.append(x if x is not None else np.zeros((64, 64), np.float32)); snr.append(info.get("snr", np.nan))
    X = np.array(X); ok = np.array(ok)
    with torch.no_grad():
        f = torch.cat([net(torch.tensor(X[b:b + 256])[:, None]) for b in range(0, len(X), 256)]).numpy()
    S = np.array([C.fft_handedness(x) for x in X])
    return f, S, ok, np.array(snr)
def thresh(score, y, ok, acc):
    a = np.abs(score); order = np.sort(a[ok])
    for t in np.quantile(order, np.linspace(0, 0.98, 99)):
        k = ok & (a > t)
        if k.sum() > 30 and ((score[k] > 0) == y[k]).mean() >= acc:
            return float(t), float(1 - ((score[k] > 0) == y[k]).mean()), float(k.sum() / ok.sum())
    return None, None, None
out = {}
f, S, ok, snr0 = run(raw)
signB = 1.0 if ((S[ok] > 0) == y[ok]).mean() > 0.5 else -1.0
out["clean"] = {"acc_A_all": float(((f[ok] > 0) == y[ok]).mean()), "acc_B_all": float(((signB * S[ok] > 0) == y[ok]).mean()),
                "signB": signB, "snr_median": float(np.nanmedian(snr0)),
                "A": thresh(f, y, ok, 0.90), "B": thresh(signB * S, y, ok, 0.85)}
print("clean", out["clean"], flush=True)
rng = np.random.default_rng(5)
for ds in sys.argv[1:]:
    from r158_common import prep_subset
    t = prep_subset(ds)["snr"]
    snr_t = rng.choice(t[np.isfinite(t)], len(raw))
    f, S, ok, snr1 = run(raw, rng, snr_t)
    o = {"snr_target_median": float(np.median(snr_t)), "snr_degraded_median": float(np.nanmedian(snr1)),
         "acc_A_all": float(((f[ok] > 0) == y[ok]).mean()), "acc_B_all": float(((signB * S[ok] > 0) == y[ok]).mean()),
         "A": thresh(f, y, ok, 0.90), "B": thresh(signB * S, y, ok, 0.85)}
    # error of the CNN at the clean threshold too
    tA = out["clean"]["A"][0]; k = ok & (np.abs(f) > tA)
    o["A_at_clean_t"] = [tA, float(1 - ((f[k] > 0) == y[k]).mean()), float(k.sum() / ok.sum())]
    out[ds] = o; print(ds, o, flush=True)
json.dump(out, open(os.path.join(DATA, "thresholds.json"), "w"), indent=1)
