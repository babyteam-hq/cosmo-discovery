"""R-158: Euclid Q1 catalogue, one sync TAP query per MER tile (the all-sky async join did not finish in 45 min).
Galaxies: VIS < 23, Sersic VIS radius > 0.4" (cut to 0.5" later), with Sersic shape and PHZ.
Stars (null): point_like_prob > 0.95, 19 < VIS < 23, up to 40 per tile, on tiles with tileid % 8 == 0."""
import os, numpy as np, io
from concurrent.futures import ThreadPoolExecutor
from astropy.table import Table, vstack
from r158_common import *
tiles = sorted({int(l.split("TILE")[1].split("-")[0]) for l in open(os.path.join(DATA, "vis_tiles.csv")).read().split("\n")[1:] if l})
od = os.path.join(DATA, "euclid_cat"); os.makedirs(od, exist_ok=True)
G = """SELECT m.object_id, m.ra, m.dec, m.tileid, m.flux_detection_total, m.point_like_prob, m.spurious_prob, m.det_quality_flag,
 mo.sersic_sersic_vis_radius, mo.sersic_sersic_vis_axis_ratio, mo.sersic_sersic_vis_index, mo.sersic_angle,
 p.phz_median, p.phz_70_int1, p.phz_70_int2
 FROM euclid_q1_mer_catalogue m JOIN euclid_q1_mer_morphology mo ON m.object_id=mo.object_id
 JOIN euclid_q1_phz_photo_z p ON m.object_id=p.object_id
 WHERE m.tileid={t} AND m.flux_detection_total > 2.29 AND mo.sersic_sersic_vis_radius > 0.4"""
S = """SELECT TOP 40 m.object_id, m.ra, m.dec, m.tileid, m.flux_detection_total, m.point_like_prob
 FROM euclid_q1_mer_catalogue m WHERE m.tileid={t} AND m.point_like_prob > 0.95 AND m.flux_detection_total > 2.29
 AND m.flux_detection_total < 91.2 AND m.spurious_prob < 0.1"""
def one(t):
    for kind, q in (("gal", G),) + ((("star", S),) if t % 8 == 0 else ()):
        f = os.path.join(od, f"{kind}_{t}.csv")
        if os.path.exists(f) and os.path.getsize(f) > 0: continue
        try:
            r = tap(q.format(t=t)); open(f, "wb").write(r)
        except Exception as e:
            print("fail", t, kind, e, flush=True)
    return t
t0 = time.time()
with ThreadPoolExecutor(16) as ex:
    for j, t in enumerate(ex.map(one, tiles)):
        if j % 25 == 0: print(j, len(tiles), round(time.time() - t0), flush=True)
for kind in ("gal", "star"):
    ts = [Table.read(os.path.join(od, f"{kind}_{t}.csv"), format="ascii.csv") for t in tiles if os.path.exists(os.path.join(od, f"{kind}_{t}.csv")) and os.path.getsize(os.path.join(od, f"{kind}_{t}.csv")) > 0]
    ts = [x for x in ts if len(x)]
    vstack(ts).write(os.path.join(DATA, f"euclid_{kind}.fits"), overwrite=True)
print("done", round(time.time() - t0))
