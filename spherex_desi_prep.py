"""R-138: compact DESI DR1 NGC clustering catalogues (data and one random file per tracer) to RA, DEC, Z, WEIGHT npz,
and the footprint map (NSIDE 32 random counts per tracer) used to fix the field.

Usage: python -u spherex_desi_prep.py <desi_dir>
"""
import os
import sys

import numpy as np
from astropy.io import fits
import healpy as hp

TRACERS = ["BGS_ANY", "LRG", "ELG_LOPnotqso", "QSO"]


def compact(path, out):
    if os.path.exists(out):
        return np.load(out)
    with fits.open(path, memmap=True) as h:
        d = h[1].data
        cols = dict(ra=np.asarray(d["RA"], np.float64), dec=np.asarray(d["DEC"], np.float64),
                    z=np.asarray(d["Z"], np.float32), w=np.asarray(d["WEIGHT"], np.float32))
    np.savez(out, **cols)
    return np.load(out)


def main():
    dd = sys.argv[1]
    fp = {}
    for t in TRACERS:
        g = compact(f"{dd}/{t}_NGC_clustering.dat.fits", f"{dd}/{t}_dat.npz")
        r = compact(f"{dd}/{t}_NGC_0_clustering.ran.fits", f"{dd}/{t}_ran.npz")
        print(t, "data", len(g["ra"]), "randoms", len(r["ra"]), flush=True)
        fp[t] = np.bincount(hp.ang2pix(32, r["ra"], r["dec"], lonlat=True), minlength=hp.nside2npix(32))
    np.savez(f"{dd}/footprint32.npz", **fp)


if __name__ == "__main__":
    main()
