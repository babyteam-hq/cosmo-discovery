"""CF4 follow-up 2: the cf4_robust.py test on the DES-SN5YR (Dovekie release) low-z SNe, 0.01<z<0.1, with DES's own
STAT+SYS covariance (the release ships the inverse covariance, packed upper-triangular, HD order). Caveats: the DES
low-z sample is CfA3/4 and Foundation SNe that are also in Pantheon+ (independent photometry pipeline, calibration
(Dovekie) and bias corrections, not independent objects), and the release has no SN coordinates for them, so RA/Dec
come from Pantheon+ by CID. zCMB and VPEC (2M++) from the DES metadata.
usage: python cf4_des.py data/des5yr_repo/4_DISTANCES_COVMAT data/Pantheon+SH0ES.dat data/new_CF4pp_mean_std_grids.npz
"""
import sys, os, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.ndimage import map_coordinates
from astropy.coordinates import SkyCoord
import astropy.units as u
C = 299792.458; OM = 0.334; N = 128; L = 1000.0; DELTA = L/N
def dc_h(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3 + 1-OM)
    return C/100*np.array([quad(e, 0, zi)[0] for zi in np.atleast_1d(z)])
def mu(z): return 5*np.log10((1+z)*dc_h(z)/0.7) + 25
def zhd(zcmb, v): return (1+zcmb)/(1+v/C) - 1

def main(desdir, pdat, grid):
    g = np.load(grid); V = g['v_mean']
    rd = lambda f: pd.read_csv(os.path.join(desdir, f), comment='#', sep=r'\s+').drop(columns=['VARNAMES:'], errors='ignore')
    hd = rd('DES-Dovekie_HD.csv'); md = rd('DES-Dovekie_Metadata.csv').set_index('CID')
    z = np.load(os.path.join(desdir, 'STAT+SYS.npz')); n = int(z['nsn'][0]); inv = np.zeros((n, n)); inv[np.triu_indices(n)] = z['cov']
    inv = inv + inv.T - np.diag(np.diag(inv)); cov = np.linalg.inv(inv)
    p = pd.read_csv(pdat, sep=r'\s+').drop_duplicates('CID').set_index('CID')
    hd['zCMB'] = md.loc[hd.CID, 'zCMB'].values; hd['VPEC'] = md.loc[hd.CID, 'VPEC'].values
    hd['RA'] = p.RA.reindex(hd.CID).values; hd['DEC'] = p.DEC.reindex(hd.CID).values
    base = (hd.zCMB > 0.01) & (hd.zCMB < 0.1)
    print(f"DES low-z 0.01<z<0.1: {base.sum()} SNe; with Pantheon+ coordinates: {(base & hd.RA.notna()).sum()}; "
          f"surveys {hd[base].IDSURVEY.value_counts().to_dict()}; diag cov median {np.median(np.sqrt(np.diag(cov))):.3f}")
    base &= hd.RA.notna()
    def run(mask, label):
        idx = np.where(mask)[0]
        zc = hd.zCMB.values[idx]; vp2 = hd.VPEC.values[idx]; r = dc_h(zc)
        sc = SkyCoord(ra=hd.RA.values[idx]*u.deg, dec=hd.DEC.values[idx]*u.deg, distance=r*u.Mpc, frame='fk5', equinox='J2000').supergalactic
        xyz = np.stack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value]); rhat = xyz/np.linalg.norm(xyz, axis=0)
        ci = (xyz + L/2)/DELTA - 0.5
        vcf = np.sum(np.stack([map_coordinates(V[k], ci, order=1, mode='nearest') for k in range(3)])*rhat, 0)
        y = hd.MU.values[idx] - mu(zc)
        Lc = np.linalg.cholesky(cov[np.ix_(idx, idx)]); w = lambda X: np.linalg.solve(Lc, X)
        def beta(v):
            cc = mu(zhd(zc, v)) - mu(zc)
            X = np.column_stack([np.ones(len(y)), cc]); Xw, yw = w(X), w(y)
            b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); e = np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw)))
            return f"{b[1]:5.2f}±{e[1]:4.2f}"
        A = np.column_stack([r, rhat.T])
        def split(v):
            coef, *_ = np.linalg.lstsq(A, v, rcond=None); return A@np.r_[0, coef[1:]], v - A@coef
        dcf, rcf = split(vcf); d2, r2 = split(vp2)
        print(f"{label:34s} n={len(idx):4d} | CF4++ full {beta(vcf)} dipole {beta(dcf)} rest {beta(rcf)} (rms rest {np.sqrt(np.mean(rcf**2)):3.0f}) | 2M++ full {beta(vp2)} dipole {beta(d2)} rest {beta(r2)} (rms rest {np.sqrt(np.mean(r2**2)):3.0f})")
    run(base, "DES low-z, all 0.01<z<0.1")
    run(base & (hd.zCMB < 0.05), "z<0.05")
    run(base & (hd.zCMB >= 0.05), "z>=0.05")
    for sid, cnt in hd[base].IDSURVEY.value_counts().items():
        if cnt >= 15: run(base & (hd.IDSURVEY != sid), f"  without survey {sid} ({cnt})")

main(*sys.argv[1:4])
