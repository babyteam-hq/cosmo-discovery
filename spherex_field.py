"""R-143: rebuild R-138's field (data/spherex/region32.npy) from its stated rule, and the QR2 D1-D3 image list over it.

Field (results-R138.md, "Fixed before any cross-correlation"): NSIDE-32 RING pixels where the DESI DR1 NGC BGS_ANY, LRG and
ELG_LOPnotqso random counts are each >= 0.8 of their 90th percentile, and |b| > 40 deg. Images: every QR2 wide-survey
D1-D3 image whose centre falls in the field or one NSIDE-32 pixel around it, in a fixed random order (seed 138).

Usage: python -u spherex_field.py region | list
"""
import csv
import io
import os
import sys
import urllib.parse
import urllib.request

import numpy as np
import healpy as hp

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")
SP = os.path.join(D, "spherex")
TAP = "https://irsa.ipac.caltech.edu/TAP/sync"
S3 = "https://nasa-irsa-spherex.s3.amazonaws.com/"


def tap(q, tries=4):
    data = urllib.parse.urlencode({"QUERY": q, "FORMAT": "csv"}).encode()
    for k in range(tries):
        try:
            with urllib.request.urlopen(TAP, data=data, timeout=1800) as r:
                return list(csv.DictReader(io.StringIO(r.read().decode())))
        except Exception as e:  # IRSA answers 504 under load: retry
            print("TAP retry", k, repr(e)[:100], flush=True)
            if k == tries - 1:
                raise


def region():
    fp = np.load(os.path.join(D, "desi_dr1", "footprint32.npz"))
    ok = np.ones(hp.nside2npix(32), bool)
    for t in ("BGS_ANY", "LRG", "ELG_LOPnotqso"):
        c = fp[t].astype(float)
        ok &= c >= 0.8 * np.percentile(c[c > 0], 90)
    th, ph = hp.pix2ang(32, np.arange(hp.nside2npix(32)))
    vec = hp.ang2vec(th, ph)
    rot = hp.Rotator(coord=["C", "G"])
    l, b = rot(th, ph)
    ok &= np.abs(90 - np.degrees(l)) > 40
    reg = np.where(ok)[0]
    os.makedirs(SP, exist_ok=True)
    np.save(os.path.join(SP, "region32.npy"), reg)
    ra, dec = hp.pix2ang(32, reg, lonlat=True)
    print("field pixels", len(reg), "deg2", len(reg) * hp.nside2pixarea(32, True),
          "RA", ra.min(), ra.max(), "Dec", dec.min(), dec.max())


def image_list():
    reg = np.load(os.path.join(SP, "region32.npy"))
    pad = set(reg.tolist())
    for p in reg:
        pad.update(hp.get_all_neighbours(32, p).tolist())
    pad.discard(-1)
    padm = np.zeros(hp.nside2npix(32), bool)
    padm[list(pad)] = True
    ra, dec = hp.pix2ang(32, reg, lonlat=True)
    rows = []
    for det in (1, 2, 3):
        cache = os.path.join(SP, f"list_D{det}.csv")
        if os.path.exists(cache):
            rows += [(r[0], int(r[1]), float(r[2]), r[3]) for r in csv.reader(open(cache))]
            continue
        n_before = len(rows)
        q = (f"SELECT obs_id, s_ra, s_dec, t_min, t_max FROM spherex.obscore WHERE obs_collection='spherex_qr2' "
             f"AND energy_bandpassname='SPHEREx-D{det}' AND s_ra BETWEEN {ra.min()-4} AND {ra.max()+4} "
             f"AND s_dec BETWEEN {dec.min()-4} AND {dec.max()+4}")
        ob = tap(q)
        q2 = f"SELECT uri FROM spherex.artifact WHERE uri LIKE '%/qr2/level2/%D{det}_spx_l2b%'"
        art = tap(q2)
        uri = {}
        for a in art:
            base = a["uri"].rsplit("/", 1)[1]          # level2_<obs>D<det>_spx_<ver>.fits
            key = base[len("level2_"):].split(f"D{det}_spx")[0]
            uri[key] = a["uri"].replace("ibe/data/spherex/", S3)
        n0 = 0
        for o in ob:
            pix = hp.ang2pix(32, float(o["s_ra"]), float(o["s_dec"]), lonlat=True)
            if not padm[pix] or o["obs_id"] not in uri:
                n0 += o["obs_id"] not in uri
                continue
            rows.append((o["obs_id"], det, 0.5 * (float(o["t_min"]) + float(o["t_max"])), uri[o["obs_id"]]))
        with open(cache, "w") as f:
            csv.writer(f).writerows(rows[n_before:])
        print(f"D{det}: obscore {len(ob)}, artifacts {len(art)}, in field {sum(r[1]==det for r in rows)}, no uri {n0}",
              flush=True)
    rng = np.random.default_rng(138)
    order = rng.permutation(len(rows))
    with open(os.path.join(SP, "list_d123.csv"), "w") as f:
        f.write("obs_id,detector,t_mid,url\n")
        for i in order:
            f.write("%s,%d,%.6f,%s\n" % rows[i])
    print("images", len(rows))


if __name__ == "__main__":
    {"region": region, "list": image_list}[sys.argv[1]]()
