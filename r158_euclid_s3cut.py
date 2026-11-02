"""R-158: Euclid Q1 VIS stamps from the MER BGSUB mosaics on IRSA's S3 bucket, read in row stripes by HTTP range
(nothing written to disk but the outputs). IRSA's IBE cutout service reset connections under load (~10/s at best), so
this replaced r158_euclid_stamps.py. Each stamp is flipped to sky parity and prepped (sign-free) in the worker.
Saves per tile: object_id, class, prepped 64x64 (float16), prep info; raw stamps for objects with object_id % 100 == 0
(for the PSF re-convolution test) and for all stars (PSF). Targets: euclid_targets_{disc,ell2000,stars800}.fits."""
import os, sys, time, subprocess, numpy as np
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, vstack
from r158_common import DATA
import r158_classify as C
B = "https://nasa-irsa-euclid-q1.s3.amazonaws.com/"
od = os.path.join(DATA, "euclid_s3"); os.makedirs(od, exist_ok=True)

def get(url, lo, hi):
    for k in range(6):
        r = subprocess.run(["curl", "-sS", "--fail", "--max-time", "300", "-r", f"{lo}-{hi}", url], capture_output=True)
        if r.returncode == 0 and len(r.stdout) == hi - lo + 1:
            return r.stdout
        time.sleep(2 ** k)
    raise RuntimeError(f"range {url} {lo}")

def key(tile):
    r = subprocess.run(["curl", "-sS", B + f"?list-type=2&prefix=q1/MER/{tile}/VIS/"], capture_output=True, text=True).stdout
    ks = [k.split("<")[0] for k in r.split("<Key>")[1:] if "BGSUB-MOSAIC-VIS" in k]
    return ks[0] if ks else None

def do_tile(args):
    tile, rows = args
    of = os.path.join(od, f"{tile}.npz")
    if os.path.exists(of):
        return tile, -1
    k = key(tile)
    if k is None:
        return tile, 0
    url = B + k
    h = get(url, 0, 2880 * 12 - 1)
    for i in range(0, len(h), 80):
        if h[i:i + 80].startswith(b"END" + b" " * 77):
            hb = h[:i + 80]; start = ((i + 80 + 2879) // 2880) * 2880; break
    hdr = fits.Header.fromstring(hb.decode("latin1"))
    nx, ny = hdr["NAXIS1"], hdr["NAXIS2"]; w = WCS(hdr); det = np.linalg.det(w.pixel_scale_matrix)
    x, y = w.all_world2pix(rows["ra"], rows["dec"], 0); x = np.round(x).astype(int); y = np.round(y).astype(int)
    hs = np.clip(np.round(rows["size_arcsec"] / 0.1 / 2).astype(int), 30, 100)
    ok = (x - hs >= 0) & (y - hs >= 0) & (x + hs < nx) & (y + hs < ny)
    SH = 1200; stripes = []
    for a in range(0, ny, SH):
        sel = np.where(ok & (y >= a) & (y < a + SH))[0]
        if len(sel):
            stripes.append((a, sel))
    res = {k_: [] for k_ in ("object_id", "cls", "X", "r50", "snr", "px", "edge", "raw_id", "raw")}

    def fetch(st):
        a, sel = st; lo = max(0, a - 100); hi = min(ny, a + SH + 100)
        buf = get(url, start + lo * nx * 4, start + hi * nx * 4 - 1)
        return lo, np.frombuffer(buf, dtype=">f4").reshape(hi - lo, nx), sel
    with ThreadPoolExecutor(4) as ex:
        for lo, arr, sel in ex.map(fetch, stripes):
            for j in sel:
                s = np.array(arr[y[j] - hs[j] - lo:y[j] + hs[j] + 1 - lo, x[j] - hs[j]:x[j] + hs[j] + 1], dtype=np.float32)
                if det > 0:
                    s = s[:, ::-1]
                if (s != 0).mean() < 0.5:
                    continue
                oid = int(rows["object_id"][j]); c = rows["cls"][j]
                if c == "star" or oid % 100 == 0:
                    res["raw_id"].append(oid); res["raw"].append(s)
                xx, info = C.prep(s)
                if xx is None:
                    continue
                res["object_id"].append(oid); res["cls"].append(c); res["X"].append(xx.astype(np.float16))
                res["r50"].append(info["r50_pix"]); res["snr"].append(info["snr"]); res["px"].append(s.shape[0])
                res["edge"].append(info["edge"])
    np.savez(of, object_id=np.array(res["object_id"], np.int64), cls=np.array(res["cls"]),
             X=np.array(res["X"], dtype=np.float16).reshape(-1, 64, 64), r50=np.array(res["r50"]), snr=np.array(res["snr"]),
             px=np.array(res["px"]), edge=np.array(res["edge"]), det=det, raw_id=np.array(res["raw_id"], np.int64),
             raw=np.array(res["raw"] + [None], dtype=object)[:-1], allow_pickle=True)
    return tile, len(res["object_id"])

if __name__ == "__main__":
    ts = []
    for tag, fn in (("disc", "euclid_targets_disc.fits"), ("ell", "euclid_targets_ell2000.fits"), ("star", "euclid_targets_stars800.fits")):
        t = Table.read(os.path.join(DATA, fn))["object_id", "ra", "dec", "tileid", "size_arcsec"]; t["cls"] = tag; ts.append(t)
    T = vstack(ts)
    tiles = np.unique(T["tileid"])
    jobs = [(int(t), {c: np.array(T[c][T["tileid"] == t]) for c in T.colnames}) for t in tiles]
    jobs.sort(key=lambda j: -len(j[1]["ra"]))
    t0 = time.time(); n = 0
    with ProcessPoolExecutor(int(sys.argv[1]) if len(sys.argv) > 1 else 4) as ex:
        for j, (tile, k) in enumerate(ex.map(do_tile, jobs)):
            n += max(k, 0)
            if j % 10 == 0:
                print(j, len(jobs), tile, k, n, round(time.time() - t0), flush=True)
    print("done", n, round(time.time() - t0))
