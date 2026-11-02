"""R-158: prep every stamp of a dataset once (sky parity from det CD, then r158_classify.prep), save x and info.
Signs are not computed here. Usage: r158_prepall.py cosmos|euclid_disc|euclid_ell|euclid_stars|jades"""
import os, sys, glob, numpy as np
from concurrent.futures import ProcessPoolExecutor
from r158_common import DATA
import r158_classify as C
ds = sys.argv[1]
def files():
    if ds == "cosmos": return sorted(glob.glob(os.path.join(DATA, "cosmos_stamps", "*.npz")))
    if ds.startswith("euclid_"): return sorted(glob.glob(os.path.join(DATA, "euclid_stamps", ds[7:], "*.npz")))
    if ds == "jades": return [os.path.join(DATA, "jades", "stamps.npz")]
def work(f):
    d = np.load(f, allow_pickle=True)
    ids = d["rows"] if "rows" in d else (d["object_id"] if "object_id" in d else d["id"])
    det = np.broadcast_to(d["det"], (len(ids),))
    X = np.zeros((len(ids), 64, 64), np.float32); ok = np.zeros(len(ids), bool)
    r50 = np.full(len(ids), np.nan); snr = np.full(len(ids), np.nan); edge = np.zeros(len(ids), bool); px = np.zeros(len(ids), int)
    for j, s in enumerate(d["stamps"]):
        if s is None or not np.isfinite(det[j]): continue
        s = np.asarray(s)
        if s.ndim != 2 or s.shape[0] != s.shape[1] or s.shape[0] < 21 or (np.abs(s) > 0).mean() < 0.5: continue
        if det[j] > 0: s = s[:, ::-1]
        x, info = C.prep(s)
        if x is None: continue
        X[j], ok[j] = x, True; r50[j], snr[j], edge[j], px[j] = info["r50_pix"], info["snr"], info["edge"], s.shape[0]
    return ids, X, ok, r50, snr, edge, px
out = {k: [] for k in ("ids", "X", "ok", "r50", "snr", "edge", "px")}
with ProcessPoolExecutor(4) as ex:
    for res in ex.map(work, files()):
        for k, v in zip(out, res): out[k].append(v)
out = {k: np.concatenate(v) for k, v in out.items()}
m = out["ok"]
np.savez(os.path.join(DATA, f"prep_{ds}.npz"), **{k: v[m] for k, v in out.items()})
print(ds, "stamps", len(m), "prepped", m.sum())
