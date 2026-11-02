"""R-153: list the SPHEREx QR2 L2 images (all detectors) whose footprint covers each dwarf spheroidal.

An image is kept when its centre lies within 2.3 deg of the dwarf (the 3.5 x 3.5 deg detector reaches 1.75-2.5 deg from
its centre), and its footprint, checked later from its WCS, holds the ON aperture and OFF annulus.
Usage: python -u alp_list.py [fake]   -> data/alp/list_<dwarf>.csv (obs_id,detector,t_mid,url,sep_deg)
"""
import csv
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spherex_field import tap, S3  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data", "alp")
# J2000 centres (McConnachie 2012)
DWARFS = {"draco": (260.0517, 57.9153), "umi": (227.2854, 67.2225), "sculptor": (15.0392, -33.7092),
          "fornax": (39.9971, -34.4492), "sextans": (153.2625, -1.6147), "carina": (100.4029, -50.9661),
          "leo1": (152.1171, 12.3064), "leo2": (168.3700, 22.1517), "bootes1": (210.0250, 14.5000)}
# fake dwarfs (blank fields, no dark-matter overdensity): each real dwarf's ecliptic latitude, moved ~7 deg in ecliptic
# longitude (results-R153.md, the fixed plan)
FAKES = {"fk_draco": (249.6930, 62.1368), "fk_umi": (234.4373, 60.9164), "fk_sculptor": (7.2734, -37.2183),
         "fk_fornax": (47.3905, -31.5759)}
RMAX = 2.3


def sep(ra1, de1, ra2, de2):
    r1, d1, r2, d2 = map(np.radians, (ra1, de1, ra2, de2))
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    uri = {}
    weeks_done = set()

    def load_week(week):
        # the artifact table, one survey week at a time (the full table drops the connection)
        art = tap(f"SELECT uri FROM spherex.artifact WHERE uri LIKE '%/qr2/level2/{week}/%spx_l2b%'")
        for a in art:
            base = a["uri"].rsplit("/", 1)[1]
            m = base[len("level2_"):].split("_spx")[0]          # <obs_id>D<det>
            key, det = m[:m.rindex("D")], int(m[m.rindex("D") + 1:])
            uri[(key, det)] = a["uri"].replace("ibe/data/spherex/", S3)
        weeks_done.add(week)
    targets = FAKES if sys.argv[1:] == ["fake"] else DWARFS
    for name, (ra, de) in targets.items():
        dra = RMAX / max(np.cos(np.radians(de)), 0.05)
        q = (f"SELECT obs_id, energy_bandpassname, s_ra, s_dec, t_min, t_max FROM spherex.obscore "
             f"WHERE obs_collection='spherex_qr2' AND s_dec BETWEEN {de-RMAX} AND {de+RMAX} "
             f"AND s_ra BETWEEN {ra-dra} AND {ra+dra}")
        ob = tap(q)
        rows = []
        for w in sorted({o["obs_id"].rsplit("_", 2)[0] for o in ob} - weeks_done):
            load_week(w)
        for o in ob:
            det = int(o["energy_bandpassname"][-1])
            s = sep(ra, de, float(o["s_ra"]), float(o["s_dec"]))
            if s > RMAX or (o["obs_id"], det) not in uri:
                continue
            rows.append((o["obs_id"], det, 0.5 * (float(o["t_min"]) + float(o["t_max"])), uri[(o["obs_id"], det)], s))
        with open(os.path.join(OUT, f"list_{name}.csv"), "w") as f:
            f.write("obs_id,detector,t_mid,url,sep_deg\n")
            for r in sorted(rows):
                f.write("%s,%d,%.6f,%s,%.3f\n" % r)
        per = {d: sum(r[1] == d for r in rows) for d in range(1, 7)}
        print(name, "obscore", len(ob), "kept", len(rows), per, f"({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
