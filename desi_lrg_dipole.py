"""DESI DR1 LRG number-count dipole (sweep-2.md row 3): data counts (WEIGHT) over random counts per nside-32 pixel, both caps,
fitted with a monopole and a dipole (weighted least squares, Poisson weights from the random counts), in redshift shells,
against the kinematic expectation D = [2 + x(1 + alpha)] beta ~ 0.005 toward the CMB apex (l, b) = (264, 48). The randoms
carry the DESI selection (footprint, completeness); imaging systematics are in WEIGHT_SYS. The DR1 footprint is two caps,
so the dipole is partly degenerate with the mask; the fit includes the full dipole and the result is read against that.
Published: the QSO count dipole (Wu & Xia 2025, arXiv:2608.30914) and redshift dipoles; no LRG count dipole found.
usage: python desi_lrg_dipole.py data/desi_lss
"""
import sys, numpy as np, healpy as hp
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
NS = 32

def main(lss):
    dat = [fits.getdata(f"{lss}/LRG_{c}_clustering.dat.fits", 1) for c in ('NGC', 'SGC')]
    ran = [fits.getdata(f"{lss}/LRG_{c}_0_clustering.ran.fits", 1) for c in ('NGC', 'SGC')]
    ra = np.concatenate([d['RA'] for d in dat]); dec = np.concatenate([d['DEC'] for d in dat]); z = np.concatenate([d['Z'] for d in dat]); w = np.concatenate([d['WEIGHT'] for d in dat])
    rra = np.concatenate([r['RA'] for r in ran]); rdec = np.concatenate([r['DEC'] for r in ran]); rz = np.concatenate([r['Z'] for r in ran]); rw = np.concatenate([r['WEIGHT'] for r in ran])
    npix = hp.nside2npix(NS); vec = np.array(hp.pix2vec(NS, np.arange(npix))).T
    apex = SkyCoord(l=264.021*u.deg, b=48.253*u.deg, frame='galactic').icrs; av = hp.ang2vec(apex.ra.deg, apex.dec.deg, lonlat=True)
    print(f"LRG: {len(z)} galaxies, {len(rz)} randoms (one file per cap)")
    for lo, hi in [(0.4, 1.1), (0.4, 0.6), (0.6, 0.8), (0.8, 1.1)]:
        sd = (z > lo) & (z < hi); sr = (rz > lo) & (rz < hi)
        Nd = np.bincount(hp.ang2pix(NS, ra[sd], dec[sd], lonlat=True), weights=w[sd], minlength=npix); Nr = np.bincount(hp.ang2pix(NS, rra[sr], rdec[sr], lonlat=True), weights=rw[sr], minlength=npix)
        m = Nr > 0.5*np.median(Nr[Nr > 0]); ratio = Nd[m]/Nr[m]; nbar = Nd[m].sum()/Nr[m].sum(); r = ratio/nbar - 1
        wt = Nr[m]*nbar  # Poisson: var(Nd) ~ Nd ~ nbar Nr, var(r) ~ 1/(nbar Nr)
        X = np.column_stack([np.ones(m.sum()), vec[m]]); sw = np.sqrt(wt); bb, *_ = np.linalg.lstsq(X*sw[:, None], r*sw, rcond=None); res = r - X@bb
        chi = np.sum(res**2*wt)/(m.sum()-4); cov = np.linalg.inv((X*wt[:, None]).T@X)*chi; ee = np.sqrt(np.diag(cov))
        D = bb[1:]; amp = np.linalg.norm(D); gd = SkyCoord(*[x[0] for x in hp.vec2ang(D, lonlat=True)]*u.deg, frame='icrs').galactic
        proj = D@av; eproj = np.sqrt(av@cov[1:, 1:]@av)
        print(f"  {lo}<z<{hi}: {sd.sum():7d} galaxies, {m.sum():4d} pixels; dipole {amp:.4f} toward (l,b)=({gd.l.deg:.0f},{gd.b.deg:.0f}), components ±{ee[1:].mean():.4f}; projection on the CMB apex {proj:+.4f} ± {eproj:.4f} (kinematic ~ +0.005); chi2/dof of the pixel ratios {chi:.2f}")

main(sys.argv[1] if len(sys.argv) > 1 else "data/desi_lss")
