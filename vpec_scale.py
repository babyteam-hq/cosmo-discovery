"""Is the amplitude of Pantheon+'s 2M++ peculiar-velocity correction right?

Residual y = m_b_corr - mu(zCMB). The correction shifts mu by d = mu(zHD) - mu(zCMB). Fit y = M + beta*d with the full
STAT+SYS covariance; beta = 1 means the reconstructed flow has the right amplitude, beta < 1 it over-predicts the
velocities, beta > 1 under-predicts. Then beta by distance shell, and the dipole left over after scaling.
usage: python vpec_scale.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov
"""
import sys, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.stats import chi2 as chi2d
C = 299792.458; OM = 0.334

def dl(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3 + 1-OM)
    return (1+z)*C/70*np.array([quad(e, 0, zi)[0] for zi in z])

def main(dat, covf):
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f:
        n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    m = (d.IS_CALIBRATOR == 0) & (d.zCMB > 0.01)
    idx = np.where(m)[0]
    zc, zh = d.zCMB.values[idx], d.zHD.values[idx]
    y = d.m_b_corr.values[idx] - (5*np.log10(dl(zc)) + 25)
    corr = (5*np.log10(dl(zh)) + 25) - (5*np.log10(dl(zc)) + 25)   # mag shift the VPEC correction applies
    vpec = d.VPEC.values[idx]
    L = np.linalg.cholesky(cov[np.ix_(idx, idx)]); w = lambda X: np.linalg.solve(L, X)
    def fit(cols):
        X = np.column_stack([np.ones(len(y))] + cols); Xw, yw = w(X), w(y)
        b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); r = yw - Xw@b
        return b, np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))), r@r
    _, _, c0 = fit([])
    b, e, c1 = fit([corr]); print(f"N={len(y)} chi2(no corr)={c0:.1f} chi2(beta)={c1:.1f}  beta={b[1]:.3f} ± {e[1]:.3f}  (beta=1 fixed: dchi2 vs best = {((b[1]-1)/e[1])**2:.2f})")
    print(f"rms VPEC = {np.sqrt(np.mean(vpec**2)):.0f} km/s; rms corr = {np.sqrt(np.mean(corr**2)):.4f} mag")
    print("by zCMB shell (beta per shell, others held at the global beta):")
    shells = [(0.01, 0.02), (0.02, 0.03), (0.03, 0.045), (0.045, 0.065), (0.065, 0.1), (0.1, 0.2)]
    cols = [np.where((zc >= lo) & (zc < hi), corr, 0.0) for lo, hi in shells] + [np.where(zc >= 0.2, corr, 0.0)]
    b, e, c = fit(cols)
    for (lo, hi), bi, ei in zip(shells, b[1:], e[1:]):
        s = (zc >= lo) & (zc < hi)
        print(f"  z[{lo:.3f},{hi:.3f}) n={s.sum():4d} rmsV={np.sqrt(np.mean(vpec[s]**2)):4.0f} km/s  beta={bi:6.2f} ± {ei:4.2f}  ({(bi-1)/ei:+.1f} sigma from 1)")
    print(f"  shells vs single beta: dchi2={c1-c:.2f} for {len(shells)} extra dof, p={chi2d.sf(c1-c, len(shells)):.3f}")
    # does the flow model's *direction* match? split the correction into its projection on the CMB dipole axis vs the rest
    ra, dec = np.radians(d.RA.values[idx]), np.radians(d.DEC.values[idx])
    nhat = np.stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)], -1)
    print("VPEC vs sky: mean VPEC on each axis (equatorial) =", np.round(nhat.T @ vpec / len(vpec), 1))

main(*sys.argv[1:])
