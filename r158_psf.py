"""R-158: empirical PSFs (median of re-centred star stamps) shared by r158_sims.py and r158_psf4b.py."""
import os, glob, numpy as np
from scipy import ndimage as ndi
from r158_common import DATA

def psf_from(stamps, crop):
    acc = []
    for s in stamps:
        s = np.asarray(s, float)
        if s is None or s.ndim != 2: continue
        c = (s.shape[0] - 1) // 2
        b = np.median(np.r_[s[:3].ravel(), s[-3:].ravel()]); s = s - b
        cut = s[c - crop:c + crop + 1, c - crop:c + crop + 1]
        if cut.shape != (2 * crop + 1, 2 * crop + 1) or cut.max() <= 0: continue
        yy, xx = np.unravel_index(np.argmax(cut), cut.shape)
        if abs(yy - crop) > 2 or abs(xx - crop) > 2: continue
        w = np.clip(cut, 0, None); my, mx = ndi.center_of_mass(w[crop - 3:crop + 4, crop - 3:crop + 4])
        cut = ndi.shift(cut, (3 - my, 3 - mx), order=3, mode="constant")
        if cut.sum() <= 0: continue
        acc.append(cut / cut.sum())
    p = np.median(np.array(acc), 0); p = np.clip(p, 0, None); return p / p.sum(), len(acc)
def stars(ds):
    import glob
    if ds == "cosmos":
        from astropy.table import Table
        z = Table.read(os.path.join(DATA, "cosmos_zest_z.fits"))
        out = []
        for f in sorted(glob.glob(os.path.join(DATA, "cosmos_stamps", "*.npz")))[::3]:
            d = np.load(f, allow_pickle=True)
            for r, s, in zip(d["rows"], d["stamps"]):
                if s is not None and z["acs_class_star"][r] > 0.95 and 19 < z["acs_mag_auto"][r] < 21.5:
                    out.append(s if d["det"] < 0 else s[:, ::-1])
        return out, 12
    from astropy.table import Table
    sid = set(np.array(Table.read(os.path.join(DATA, "euclid_targets_stars800.fits"))["object_id"]).tolist())
    out = []
    for f in sorted(glob.glob(os.path.join(DATA, "euclid_s3", "*.npz"))):
        d = np.load(f, allow_pickle=True)
        out += [s for i, s in zip(d["raw_id"], d["raw"]) if int(i) in sid]
    return out, 7
