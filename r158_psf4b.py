"""R-158 null 4b: real sample galaxies re-convolved with the dataset's empirical PSF and with its mirror image.
dA = A(PSF) - A(mirror PSF) on the same galaxies measures the response of A to a PSF chirality as large as the PSF's own
(a PSF fully chiral at its own scale); the PSF term already in the data is bounded by |dA|/2.
Usage: r158_psf4b.py cosmos|euclid_disc N"""
import os, sys, json, glob, numpy as np, torch
from scipy.signal import fftconvolve
from r158_common import DATA, cosmos_classes
import r158_classify as C
ds, N = sys.argv[1], int(sys.argv[2]); rng = np.random.default_rng(4)
from r158_psf import psf_from, stars
st_, crop_ = stars(ds); psf, _ = psf_from(st_, crop_)
T = json.load(open(os.path.join(DATA, "thresholds_sims.json")))[ds]["f"]["t"]
raw = []
if ds == "cosmos":
    cls, _ = cosmos_classes()
    for f in sorted(glob.glob(os.path.join(DATA, "cosmos_stamps", "*.npz"))):
        d = np.load(f, allow_pickle=True)
        for r, s in zip(d["rows"], d["stamps"]):
            if s is not None and cls[r] == "disc": raw.append(s if d["det"] < 0 else s[:, ::-1])
else:
    for f in sorted(glob.glob(os.path.join(DATA, "euclid_s3", "*.npz"))):
        d = np.load(f, allow_pickle=True)
        disc = set(d["object_id"][d["cls"] == "disc"].tolist())
        raw += [s for i, s in zip(d["raw_id"], d["raw"]) if i in disc]   # already sky parity
idx = rng.choice(len(raw), min(N, len(raw)), replace=False)
net = C.make_net(); net.load_state_dict(torch.load(os.path.join(DATA, "net.pt"))); net.eval()
res = {}
for tag, P in (("psf", psf), ("mirror", psf[:, ::-1])):
    X = []
    for i in idx:
        x, _ = C.prep(fftconvolve(np.asarray(raw[i], float), P, mode="same"))
        X.append(x if x is not None else np.zeros((64, 64), np.float32))
    with torch.no_grad():
        f = torch.cat([net(torch.tensor(np.array(X[b:b + 256]))[:, None]) for b in range(0, len(X), 256)]).numpy()
    res[tag] = f
k = (np.abs(res["psf"]) > T) | (np.abs(res["mirror"]) > T)
a1 = np.sign(res["psf"][np.abs(res["psf"]) > T]).mean(); a2 = np.sign(res["mirror"][np.abs(res["mirror"]) > T]).mean()
flip = np.mean(np.sign(res["psf"][k]) != np.sign(res["mirror"][k]))
o = {"n": int(len(idx)), "A_psf": float(a1), "A_mirror": float(a2), "dA": float(a1 - a2),
     "dA_err": float(np.sqrt(2 * flip / k.sum()) if k.sum() else np.nan), "sign_flips": float(flip), "mean_df": float(np.mean(res["psf"] - res["mirror"])),
     "mean_df_err": float(np.std(res["psf"] - res["mirror"]) / np.sqrt(len(idx)))}
json.dump(o, open(os.path.join(DATA, f"psf4b_{ds}.json"), "w")); print(ds, o)
