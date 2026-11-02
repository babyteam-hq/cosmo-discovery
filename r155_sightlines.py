"""R-155: DESI DR1 sightlines and absorbers, SDSS sightlines removed.

Sightlines: DESI DR1 main-dark quasars (QSO_cat_iron_cumulative_v0, the afterburner catalogue the Mg II VAC searched;
every VAC absorber's TARGETID is main-dark). One row per TARGETID: TSNR2_QSO summed over the catalogue's per-tile
rows (that sum equals the VAC's healpix-coadd TSNR2_QSO for 98.9 % of VAC quasars), redshift from the row with the
highest TSNR2_QSO, as are the
Legacy Surveys g and r fluxes (nanomaggies; the quasar's brightness at the Mg II slice, ~5000 A, which sets detection
completeness together with the exposure depth TSNR2). SDSS removal: within 1.5" of a DR16Q quasar or of a quasar in the Z&M DR7/DR12 searched lists.

Absorbers: MgII-Absorbers-DR1 (the PI file, z_abs > z_qso + 5000 km/s, is already separate), attached to sightlines
by TARGETID.

Writes data/r155/desi_sightlines.npz and desi_absorbers.npz (and desi_absorbers_all.npz, every VAC row with an
'sdss' flag, for the EW calibration against Z&M).
"""
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, "research/cosmo-discovery")
from r155_common import DATA, sky_match  # noqa: E402


def main():
    q = fits.getdata(f"{DATA}/QSO_cat_iron_cumulative_v0.fits", 1)
    q = q[(q["SURVEY"] == "main") & (q["PROGRAM"] == "dark")]
    tid = np.asarray(q["TARGETID"])
    u, inv = np.unique(tid, return_inverse=True)
    tsnr = np.bincount(inv, weights=q["TSNR2_QSO"])
    o = np.lexsort((-q["TSNR2_QSO"], tid))
    first = o[np.r_[True, tid[o][1:] != tid[o][:-1]]]  # best row per TARGETID, in TARGETID order (= u order)
    ra, dec, zq, fg, fr = (np.asarray(q[c][first], float) for c in ("TARGET_RA", "TARGET_DEC", "Z", "FLUX_G", "FLUX_R"))
    print(f"DESI main-dark quasars: {len(u)} unique TARGETIDs")

    sd = fits.open(f"{DATA}/DR16Q_v4.fits", memmap=True)[1].data
    sra, sdec = np.asarray(sd["RA"], float), np.asarray(sd["DEC"], float)
    zs = [fits.getdata(f"{DATA}/Searched_Quasars_{x}_107.fits", 1) for x in ("SDSS_DR7", "BOSS_DR12")]
    sra = np.concatenate([sra] + [np.asarray(z["RA"], float) for z in zs])
    sdec = np.concatenate([sdec] + [np.asarray(z["DEC"], float) for z in zs])
    i1, _ = sky_match(ra, dec, sra, sdec, 1.5)
    sdss = np.zeros(len(u), bool)
    sdss[i1] = True
    print(f"in DR16Q or Z&M searched lists (1.5\"): {sdss.sum()}; independent sightlines: {(~sdss).sum()}")
    np.savez(f"{DATA}/desi_sightlines.npz", tid=u, ra=ra, dec=dec, zq=zq, tsnr2=tsnr, sdss=sdss, flux_g=fg, flux_r=fr)

    m = fits.getdata(f"{DATA}/MgII-Absorbers-DR1.fits", 1)
    j = np.searchsorted(u, m["TARGETID"])
    ok = u[np.clip(j, 0, len(u) - 1)] == m["TARGETID"]
    print(f"VAC absorbers {len(m)}, attached to a sightline {ok.sum()}")
    a = dict(sl=j[ok], z=np.asarray(m["Z_MGII"][ok], float), w=np.asarray(m["EW_2796"][ok], float),
             w3=np.asarray(m["EW_2803"][ok], float), ra=np.asarray(m["RA"][ok], float),
             dec=np.asarray(m["DEC"][ok], float), zq=np.asarray(m["Z_QSO"][ok], float))
    np.savez(f"{DATA}/desi_absorbers_all.npz", sdss=sdss[a["sl"]], **a)


if __name__ == "__main__":
    main()
