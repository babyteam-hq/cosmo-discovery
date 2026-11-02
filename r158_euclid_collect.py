"""R-158: gather the per-tile Euclid S3 outputs into prep_euclid_{disc,ell,stars}.npz (the r158_measure input)."""
import os, glob, numpy as np
from r158_common import DATA
acc = {c: {k: [] for k in ("ids", "X", "r50", "snr", "px")} for c in ("disc", "ell", "star")}
for f in sorted(glob.glob(os.path.join(DATA, "euclid_s3", "*.npz"))):
    d = np.load(f, allow_pickle=True)
    for c in acc:
        m = d["cls"] == c
        if not m.any(): continue
        acc[c]["ids"].append(d["object_id"][m]); acc[c]["X"].append(d["X"][m])
        for k in ("r50", "snr", "px"): acc[c][k].append(d[k][m])
for c, v in acc.items():
    v = {k: np.concatenate(x) for k, x in v.items()}
    _, u = np.unique(v["ids"], return_index=True)   # objects in two tiles' overlap: keep one
    v = {k: x[u] for k, x in v.items()}
    np.savez(os.path.join(DATA, f"prep_euclid_{'stars' if c == 'star' else c}.npz"), **v); print(c, len(v["ids"]))
