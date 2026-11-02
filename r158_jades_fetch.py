"""R-158: JADES GOODS-S DR2 F200W SCI extension (25000 x 27700 float32, 2.77 GB) by parallel HTTP ranges from MAST,
saved as a raw big-endian array plus its header. Then 263 stamps at Shamir's published positions (2502.18781 T1-2)."""
import os, subprocess, numpy as np, time
from concurrent.futures import ThreadPoolExecutor
from astropy.io import fits
from r158_common import DATA
URL = "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HLSP/jades/dr2/goods-s/images/hlsp_jades_jwst_nircam_goods-s-deep_f200w_v2.0_drz.fits"
od = os.path.join(DATA, "jades"); os.makedirs(od, exist_ok=True)
head = open(os.path.join(od, "head.bin"), "rb").read()
# primary header, then SCI header: find the 2nd END card, data start at next 2880 boundary
pos, ends = 0, []
while len(ends) < 2:
    card = head[pos:pos + 80]
    if card.startswith(b"END" + b" " * 77): ends.append(pos)
    pos += 80
start = ((ends[1] + 80 + 2879) // 2880) * 2880
hdr = fits.Header.fromstring(head[(((ends[0] + 80 + 2879) // 2880) * 2880):ends[1] + 80].decode("latin1"))
nx, ny = hdr["NAXIS1"], hdr["NAXIS2"]; nbytes = nx * ny * 4
hdr.tofile(os.path.join(od, "sci_header.txt"), sep="\n", endcard=False, padding=False, overwrite=True)
raw = os.path.join(od, "sci.raw")
NCH = 32; ch = (nbytes + NCH - 1) // NCH
def get(k):
    f = raw + f".{k:02d}"; lo = start + k * ch; hi = min(start + nbytes, lo + ch) - 1
    if os.path.exists(f) and os.path.getsize(f) == hi - lo + 1: return k
    for t in range(6):
        r = subprocess.run(["curl", "-sS", "--fail", "--max-time", "1800", "-r", f"{lo}-{hi}", "-o", f, URL], capture_output=True)
        if r.returncode == 0 and os.path.getsize(f) == hi - lo + 1: return k
        time.sleep(2 ** t)
    raise RuntimeError(f"chunk {k}")
t0 = time.time()
with ThreadPoolExecutor(12) as ex:
    for k in ex.map(get, range(NCH)): print("chunk", k, round(time.time() - t0), flush=True)
with open(raw, "wb") as o:
    for k in range(NCH):
        o.write(open(raw + f".{k:02d}", "rb").read()); os.remove(raw + f".{k:02d}")
print("done", nx, ny, round(time.time() - t0), flush=True)
