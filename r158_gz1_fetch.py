"""R-158: training set. Galaxy Zoo 1 (Lintott+2011 table 2) spirals with P_CW or P_ACW > 0.8 (NVOTE >= 15),
r-band Legacy Surveys DR10 FITS cutouts (96 px at 0.262"/px), stored as one npz. Labels are GZ1's sky-view CW/ACW."""
import os, sys, io, numpy as np, pandas as pd, subprocess, time
from concurrent.futures import ThreadPoolExecutor
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from r158_common import DATA
N = int(sys.argv[1]) if len(sys.argv) > 1 else 7000
d = pd.read_csv(os.path.join(DATA, "gz1", "GalaxyZoo1_DR_table2.csv"))
cw = d[(d.P_CW > 0.8) & (d.NVOTE >= 15)].sample(N, random_state=1)
acw = d[(d.P_ACW > 0.8) & (d.NVOTE >= 15)].sample(N, random_state=2)
s = pd.concat([cw.assign(lab=1), acw.assign(lab=-1)]).sample(frac=1, random_state=3).reset_index(drop=True)
c = SkyCoord(s.RA.values, s.DEC.values, unit=(u.hourangle, u.deg))
s["ra"], s["dec"] = c.ra.deg, c.dec.deg
out = os.path.join(DATA, "gz1", "stamps"); os.makedirs(out, exist_ok=True)
def get(i):
    f = os.path.join(out, f"{i}.fits")
    if os.path.exists(f) and os.path.getsize(f) > 1000: return
    url = f"https://www.legacysurvey.org/viewer/cutout.fits?ra={s.ra[i]:.6f}&dec={s.dec[i]:.6f}&layer=ls-dr10&pixscale=0.262&size=96&bands=r"
    for k in range(4):
        r = subprocess.run(["curl", "-sS", "--fail", "--max-time", "60", "-o", f, url], capture_output=True)
        if r.returncode == 0: return
        time.sleep(1 + k)
t0 = time.time()
with ThreadPoolExecutor(48) as ex:
    for j, _ in enumerate(ex.map(get, range(len(s)))):
        if j % 1000 == 0: print(j, round(time.time() - t0), flush=True)
imgs, labs, keep = [], [], []
for i in range(len(s)):
    f = os.path.join(out, f"{i}.fits")
    try:
        a = fits.getdata(f).astype(np.float32)
        if a.ndim == 3: a = a[0]
        if a.shape != (96, 96) or not np.isfinite(a).all(): continue
        imgs.append(a); labs.append(s.lab[i]); keep.append(i)
    except Exception: pass
np.savez_compressed(os.path.join(DATA, "gz1", "train.npz"), x=np.array(imgs), y=np.array(labs), objid=s.OBJID.values[keep])
print("done", len(imgs), "s", round(time.time() - t0))
