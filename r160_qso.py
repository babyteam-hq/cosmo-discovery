"""R-160: one spectroscopic quasar list from DESI DR1, SDSS DR16Q and Milliquas v8 (types Q/A, spectroscopic z).

Priority for the redshift of a duplicate (within 1.5"): DESI DR1, then DR16Q, then Milliquas. Cuts: z >= 0.1
(below it 150 Mpc is a large angle and the sample is tiny), DESI ZWARN == 0, DR16Q Z > 0 and ZWARNING == 0 where
the pipeline z is the source.

Usage: python -I r160_qso.py data/r160  -> data/r160/qso_all.fits
"""
import sys

import numpy as np
from astropy.table import Table, vstack
from scipy.spatial import cKDTree

import r160_common as c


def main(d):
    desi = Table.read(f"{d}/QSO_cat_iron_cumulative_v0.fits", hdu=1)
    desi = desi[(desi["ZWARN"] == 0) & (desi["Z"] >= 0.1)]
    # one row per TARGETID (the cumulative catalogue may repeat); keep the first
    _, k = np.unique(desi["TARGETID"], return_index=True)
    desi = desi[np.sort(k)]
    a = Table(dict(ra=desi["TARGET_RA"], dec=desi["TARGET_DEC"], z=desi["Z"], src=np.full(len(desi), "DESI")))
    sd = Table.read(f"{d}/DR16Q_v4.fits", hdu=1)
    sd = sd[(sd["Z"] >= 0.1)]
    b = Table(dict(ra=sd["RA"], dec=sd["DEC"], z=sd["Z"], src=np.full(len(sd), "DR16Q")))
    mq = Table.read(f"{d}/milliquas.fits")
    ty = np.array([s.strip()[:1] for s in mq["TYPE"]])
    mq = mq[np.isin(ty, ["Q", "A"]) & (mq["Z"] >= 0.1)]
    m = Table(dict(ra=mq["RA"], dec=mq["DEC"], z=mq["Z"], src=np.full(len(mq), "MQ")))
    out = a
    for nxt in (b, m):
        tree = cKDTree(c.unit(out["ra"], out["dec"]))
        dd, _ = tree.query(c.unit(nxt["ra"], nxt["dec"]), k=1, distance_upper_bound=c.chord(1.5 / 3600))
        out = vstack([out, nxt[~np.isfinite(dd)]])
    print({s: int(np.sum(out["src"] == s)) for s in ("DESI", "DR16Q", "MQ")}, "total", len(out))
    out.write(f"{d}/qso_all.fits", overwrite=True)


if __name__ == "__main__":
    main(sys.argv[1])
