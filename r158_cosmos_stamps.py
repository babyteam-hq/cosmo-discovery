"""R-158: COSMOS ACS F814W stamps (IRSA acs_2.0, 0.03"/px). Stream each tile, cut stamps, delete the tile.
Targets: the pre-registered disc sample (0.25<z<1.3, type 2, r50 > 9.5 px, elongation < 2.5, I < 24), ellipticals
(type 1) and stars (acs_class_star > 0.95, 19 < I < 23) as nulls. Stamp half-size = clip(4 r50, 32, 150) px."""
import os, sys, numpy as np, pandas as pd, subprocess, time
from concurrent.futures import ThreadPoolExecutor
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from r158_common import DATA
z = Table.read(os.path.join(DATA, "cosmos_zest_z.fits")).to_pandas()
z["acstile"] = z.acstile.str.decode("ascii").str.slice(0, 17)
z["row"] = np.arange(len(z))
disc = (z.lp_zbest > 0.25) & (z.lp_zbest < 1.3) & (z.type == 2) & (z.r50 > 9.5) & (z.acs_elongation < 2.5) & (z.junkflag == 0)
ell = (z.type == 1) & (z.r50 > 9.5) & (z.junkflag == 0)
star = (z.acs_class_star > 0.95) & (z.acs_mag_auto > 19) & (z.acs_mag_auto < 23)
z["cls"] = np.where(disc, "disc", np.where(ell, "ell", np.where(star, "star", "")))
z = z[z.cls != ""]
print(z.cls.value_counts().to_dict(), flush=True)
base = "https://irsa.ipac.caltech.edu/data/COSMOS/images/acs_2.0/I/"
names = [l.split('"')[1] for l in subprocess.run(["curl", "-sS", base], capture_output=True, text=True).stdout.split("href=")[1:] if "_sci_20.fits" in l.split('"')[1]]
bytile = {n[:17]: n for n in names}
tmp = os.path.join(DATA, "acs"); out = os.path.join(DATA, "cosmos_stamps"); os.makedirs(out, exist_ok=True)
def do(tile):
    of = os.path.join(out, tile + ".npz")
    if os.path.exists(of): return tile, 0
    g = z[z.acstile == tile]
    f = os.path.join(tmp, tile + ".fits")
    for k in range(5):
        r = subprocess.run(["curl", "-sS", "--fail", "--max-time", "900", "-o", f, base + bytile[tile]], capture_output=True)
        if r.returncode == 0: break
        time.sleep(2 ** k)
    with fits.open(f, memmap=True) as h:
        a = h[0].data; hd = h[0].header; w = WCS(hd)
        x, y = w.all_world2pix(g.ra.values, g.dec.values, 0)
        st, keep = [], []
        for xi, yi, r50 in zip(x, y, np.nan_to_num(g.r50.values, nan=10.0)):
            if not (np.isfinite(xi) and np.isfinite(yi)): st.append(None); continue
            hs = int(np.clip(np.ceil(4 * r50), 32, 150)); xi, yi = int(round(xi)), int(round(yi))
            if xi - hs < 0 or yi - hs < 0 or xi + hs + 1 > a.shape[1] or yi + hs + 1 > a.shape[0]:
                st.append(None); continue
            s = np.array(a[yi - hs:yi + hs + 1, xi - hs:xi + hs + 1], dtype=np.float32)
            st.append(s)
        det = np.linalg.det(w.pixel_scale_matrix)
    np.savez(of, rows=g.row.values, stamps=np.array(st, dtype=object), det=det, allow_pickle=True)
    os.remove(f)
    return tile, len(g)
tiles = sorted(set(z.acstile) & set(bytile))
t0 = time.time()
with ThreadPoolExecutor(int(sys.argv[1]) if len(sys.argv) > 1 else 10) as ex:
    for j, (t, n) in enumerate(ex.map(do, tiles)):
        if j % 20 == 0: print(j, len(tiles), t, n, round(time.time() - t0), flush=True)
print("done", round(time.time() - t0))
