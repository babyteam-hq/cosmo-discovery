"""Is the GNILC-SFD slope about GNILC's coarser, variable resolution rather than its dust? Replace SFD at the SN with
SFD averaged over a disk of radius R (R = 5, 10, 20, 40 arcmin; 24-point rings at R/2 and R plus the centre) and GNILC
likewise, and refit the Pantheon+ residual slope on (GNILC - SFD_R), (GNILC_R - SFD_R), and on the small-scale part of
SFD alone (SFD - SFD_R). If the slope lives in the large-scale calibration difference it survives smoothing both maps;
if it is SFD's small-scale structure, (SFD - SFD_R) carries it with a negative sign.
usage: python dust_resolution.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/dustmaps
"""
import sys, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.stats import chi2 as chi2d
C = 299792.458; OM = 0.334

def dl(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3+1-OM); return (1+z)*C/70*np.array([quad(e, 0, zi)[0] for zi in z])

def main(dat, covf, data_dir):
    from dustmaps.config import config; config['data_dir'] = data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.planck import PlanckGNILCQuery
    from astropy.coordinates import SkyCoord; import astropy.units as u
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f: n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    m = (d.IS_CALIBRATOR == 0) & (d.zHD > 0.01); idx = np.where(m)[0]
    z = d.zHD.values[idx]; y = d.m_b_corr.values[idx] - (5*np.log10(dl(z))+25)
    c = SkyCoord(d.RA.values[idx]*u.deg, d.DEC.values[idx]*u.deg, frame='icrs')
    sfd, gn = SFDQuery(), PlanckGNILCQuery()
    def disk(q, R):
        pts = [q(c)]
        for r in (R/2, R):
            for th in np.linspace(0, 2*np.pi, 24, endpoint=False): pts.append(q(c.directional_offset_by(th*u.rad, r*u.arcmin)))
        return np.mean(pts, axis=0)
    s0, g0 = sfd(c), gn(c)
    L = np.linalg.cholesky(cov[np.ix_(idx, idx)]); w = lambda X: np.linalg.solve(L, X)
    def fit(cols):
        X = np.column_stack([np.ones(len(y))]+cols); Xw, yw = w(X), w(y); b, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        e = np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))); r = yw-Xw@b; return b, e, r@r
    _, _, c0 = fit([])
    def line(name, x):
        b, e, c1 = fit([x]); print(f"  {name:34s} slope {b[1]:6.2f}±{e[1]:4.2f} dchi2={c0-c1:5.2f} p={chi2d.sf(c0-c1,1):.2g}  rms x {np.std(x):.4f}")
    line("GNILC - SFD (point values)", g0 - s0)
    for R in (5, 10, 20, 40):
        sR, gR = disk(sfd, R), disk(gn, R)
        print(f"R = {R}':")
        line(f"GNILC - SFD_{R}", g0 - sR); line(f"GNILC_{R} - SFD_{R}", gR - sR); line(f"SFD - SFD_{R} (small-scale SFD)", s0 - sR)
        line(f"GNILC - GNILC_{R} (small-scale GNILC)", g0 - gR)
main(*sys.argv[1:4])
