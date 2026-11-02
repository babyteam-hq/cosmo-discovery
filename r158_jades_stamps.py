"""R-158: stamps (301 px = 9") from the JADES DR2 F200W SCI array at Shamir's 263 positions."""
import os, numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from r158_common import DATA, HERE
od = os.path.join(DATA, "jades")
hdr = fits.Header.fromtextfile(os.path.join(od, "sci_header.txt"))
nx, ny = hdr["NAXIS1"], hdr["NAXIS2"]
a = np.memmap(os.path.join(od, "sci.raw"), dtype=">f4", mode="r", shape=(ny, nx))
w = WCS(hdr)
s = pd.read_csv(os.path.join(HERE, "r158_shamir_jades.csv"), comment="#")
x, y = w.all_world2pix(s.ra, s.dec, 0)
hs = 150; st = []
for xi, yi in zip(np.round(x).astype(int), np.round(y).astype(int)):
    if hs <= xi < nx - hs and hs <= yi < ny - hs:
        st.append(np.array(a[yi - hs:yi + hs + 1, xi - hs:xi + hs + 1], dtype=np.float32))
    else:
        st.append(None)
det = np.linalg.det(w.pixel_scale_matrix)
np.savez(os.path.join(od, "stamps.npz"), id=s.id.values, stamps=np.array(st, dtype=object), det=det, allow_pickle=True)
print("stamps", sum(v is not None for v in st), "det", det, "nonzero", sum(v is not None and np.abs(v).sum() > 0 for v in st))
