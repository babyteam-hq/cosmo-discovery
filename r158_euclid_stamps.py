"""R-158: Euclid Q1 VIS stamps via IRSA's IBE cutout service on the MER BGSUB VIS mosaics (0.1"/px).
Input: a FITS table with object_id, ra, dec, tileid, size_arcsec. Output: one npz per chunk of 2000 objects."""
import os, sys, io, numpy as np, subprocess, time
from concurrent.futures import ThreadPoolExecutor
from astropy.io import fits
from astropy.table import Table
from r158_common import DATA
inp, tag = sys.argv[1], sys.argv[2]
nthr = int(sys.argv[3]) if len(sys.argv) > 3 else 32
t = Table.read(inp)
uris = {}
for l in open(os.path.join(DATA, "vis_tiles.csv")).read().split("\n")[1:]:
    if l: u = l.split(",")[0]; uris[int(u.split("TILE")[1].split("-")[0])] = u
out = os.path.join(DATA, "euclid_stamps", tag); os.makedirs(out, exist_ok=True)
def one(i):
    u = uris.get(int(t["tileid"][i]))
    if u is None: return None
    url = f"https://irsa.ipac.caltech.edu/{u}?center={t['ra'][i]:.7f},{t['dec'][i]:.7f}&size={float(t['size_arcsec'][i]):.1f}arcsec"
    for k in range(5):
        r = subprocess.run(["curl", "-sS", "--fail", "--max-time", "60", url], capture_output=True)
        if r.returncode == 0 and len(r.stdout) > 2880:
            try:
                h = fits.open(io.BytesIO(r.stdout)); a = h[0].data.astype(np.float32); hd = h[0].header
                cd = np.array([[hd.get("CD1_1", 0), hd.get("CD1_2", 0)], [hd.get("CD2_1", 0), hd.get("CD2_2", 0)]])
                if not cd.any(): cd = np.diag([hd.get("CDELT1", -1), hd.get("CDELT2", 1)]) @ np.array([[hd.get("PC1_1", 1), hd.get("PC1_2", 0)], [hd.get("PC2_1", 0), hd.get("PC2_2", 1)]])
                return a, float(np.linalg.det(cd))
            except Exception: pass
        time.sleep(1 + 2 * k)
    return None
CH = 2000
t0 = time.time()
with ThreadPoolExecutor(nthr) as ex:
    for c0 in range(0, len(t), CH):
        of = os.path.join(out, f"{c0 // CH:04d}.npz")
        if os.path.exists(of): continue
        idx = list(range(c0, min(len(t), c0 + CH)))
        res = list(ex.map(one, idx))
        st = np.empty(len(idx), dtype=object); det = np.full(len(idx), np.nan)
        for j, r in enumerate(res):
            if r is not None: st[j], det[j] = r
        np.savez(of, object_id=np.array(t["object_id"][idx]), stamps=st, det=det, allow_pickle=True)
        print(c0 + len(idx), len(t), sum(r is not None for r in res), round(time.time() - t0), flush=True)
print("done", round(time.time() - t0))
