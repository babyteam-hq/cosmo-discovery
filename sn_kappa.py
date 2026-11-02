"""Pantheon+ Hubble residuals vs the Planck PR4 CMB-lensing convergence at the SN position (sweep-2.md row 17): lensing
magnification makes a SN behind an overdensity brighter by dm = -(5/ln10) kappa_SN = -2.17 kappa_SN, and kappa_SN (to the SN
redshift) is the low-redshift part of the CMB convergence (to z = 1100), so the residual should anticorrelate with the CMB
kappa with a coefficient -2.17 x cov(kappa_SN, kappa_CMB)/var(kappa_CMB), of order -0.1 to -0.4 at z ~ 0.3 on degree scales
and growing with z. Maps: PR4 2018-like minimum-variance klm minus its mean field (Carron+2022, GitHub release "Data"),
rendered at nside 512 and smoothed with Gaussian beams of FWHM 30', 1 deg and 2 deg; the PR4 analysis mask. Fit: full
STAT+SYS covariance GLS, intercept free, kappa column alone and kappa x z (the growth with redshift), for z > 0.05 and in
redshift halves. Also the ZTF sample at z < 0.1 as a null control (negligible lensing).
usage: python sn_kappa.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/pr4/PR4_variations
"""
import sys, numpy as np, pandas as pd, healpy as hp
from scipy.integrate import quad
from scipy.stats import chi2 as chi2d
from astropy.coordinates import SkyCoord
import astropy.units as u
C = 299792.458; OM = 0.334; NS = 512

def dl(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3+1-OM); return (1+z)*C/70*np.array([quad(e, 0, zi)[0] for zi in z])

def main(dat, covf, pr4):
    klm = hp.read_alm(f"{pr4}/PR42018like_klm_dat_MV.fits") - hp.read_alm(f"{pr4}/PR42018like_klm_mf_MV.fits")
    mask = hp.ud_grade(hp.read_map(f"{pr4}/mask.fits.gz"), NS)
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f: n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    m = (d.IS_CALIBRATOR == 0) & (d.zHD > 0.01); idx = np.where(m)[0]
    z = d.zHD.values[idx]; y = d.m_b_corr.values[idx] - (5*np.log10(dl(z))+25)
    g = SkyCoord(d.RA.values[idx]*u.deg, d.DEC.values[idx]*u.deg, frame='icrs').galactic
    pix = hp.ang2pix(NS, g.l.deg, g.b.deg, lonlat=True); inmask = mask[pix] > 0.5
    Lfull = cov[np.ix_(idx, idx)]
    def gls(y_, cols, sel):
        Lc = np.linalg.cholesky(Lfull[np.ix_(sel, sel)]); X = np.column_stack([np.ones(sel.sum())] + [c[sel] for c in cols]); Xw = np.linalg.solve(Lc, X); yw = np.linalg.solve(Lc, y_[sel])
        b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); e = np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))); r = yw - Xw@b
        X0 = np.linalg.solve(Lc, np.ones((sel.sum(), 1))); b0, *_ = np.linalg.lstsq(X0, yw, rcond=None); r0 = yw - X0@b0
        return b[1:], e[1:], r0@r0 - r@r
    print(f"Pantheon+: {len(idx)} SNe, {inmask.sum()} inside the PR4 lensing mask")
    for fwhm in (30, 60, 120):
        kmap = hp.alm2map(hp.smoothalm(klm.copy(), fwhm=np.radians(fwhm/60), inplace=False), NS, lmax=2048, verbose=False) if False else hp.alm2map(hp.smoothalm(klm.copy(), fwhm=np.radians(fwhm/60)), NS, lmax=2048)
        kap = kmap[pix]
        print(f"kappa smoothed to FWHM {fwhm}': rms over the masked sky {np.std(kmap[mask > 0.5]):.4f}, at the SNe {np.std(kap[inmask]):.4f}")
        for lbl, sel in [("z>0.05", inmask & (z > 0.05)), ("0.05<z<0.2", inmask & (z > 0.05) & (z < 0.2)), ("z>0.2", inmask & (z > 0.2)), ("z>0.4", inmask & (z > 0.4))]:
            b, e, dc = gls(y, [kap], sel); b2, e2, dc2 = gls(y, [kap*z], sel)
            print(f"  {lbl:10s} n={sel.sum():4d} slope on kappa {b[0]:+6.2f}±{e[0]:4.2f} (dchi2 {dc:4.1f}, p={chi2d.sf(dc,1):.2f}) | on kappa x z {b2[0]:+6.2f}±{e2[0]:4.2f} | magnification would give ~ -2.17 x (low-z share of the CMB kernel, ~0.05-0.2 at these z)")

main(*sys.argv[1:4])
