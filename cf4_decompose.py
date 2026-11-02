"""Decompose the CF4++ radial velocity at Pantheon+ positions into monopole (a*r: an H0 offset), dipole (bulk flow) and
the rest, then redo the beta fit on the rest. Also the grid's own bulk flow in spheres, to compare with Watkins+2023.
usage: python cf4_decompose.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/new_CF4pp_mean_std_grids.npz
"""
import sys, numpy as np, pandas as pd
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

def main(dat, covf, grid, zmax=0.1):
    g = np.load(grid)
    # 1. the grid's own monopole and bulk flow in spheres (voxel centres)
    ax = (np.arange(N) + 0.5)*DELTA - L/2
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij'); R = np.sqrt(X**2 + Y**2 + Z**2)
    V = g['v_mean']; VR = (V[0]*X + V[1]*Y + V[2]*Z)/np.maximum(R, 1e-9)
    print("CF4++ grid in spheres: R(Mpc/h), n_vox, bulk |B| km/s (and direction SGX,SGY,SGZ), monopole slope a=<vr/r> km/s per Mpc/h, rms vr")
    for Rmax in [50, 100, 150, 200, 250, 300]:
        s = R < Rmax
        B = np.array([V[k][s].mean() for k in range(3)])
        a = np.mean(VR[s]/np.maximum(R[s], 1e-9))
        print(f"  R<{Rmax:3d}: n={s.sum():7d} |B|={np.linalg.norm(B):5.0f} dir={np.round(B/np.linalg.norm(B),2)} a={a:+.3f} (dH0={a*100:+.1f}% of H0) rms vr={np.sqrt(np.mean(VR[s]**2)):4.0f}")
    # 2. at SN positions
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f:
        n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    m = (d.IS_CALIBRATOR == 0) & (d.zCMB > 0.01) & (d.zCMB < zmax); idx = np.where(m)[0]
    zc = d.zCMB.values[idx]; zh2 = d.zHD.values[idx]; vp2 = d.VPEC.values[idx]
    r = dc_h(zc)
    sc = SkyCoord(ra=d.RA.values[idx]*u.deg, dec=d.DEC.values[idx]*u.deg, distance=r*u.Mpc, frame='fk5', equinox='J2000').supergalactic
    xyz = np.stack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value]); rhat = xyz/np.linalg.norm(xyz, axis=0)
    ci = (xyz + L/2)/DELTA - 0.5
    vcf = np.sum(np.stack([map_coordinates(V[k], ci, order=1, mode='nearest') for k in range(3)])*rhat, 0)
    # vr = a*r + B.rhat + rest  (plain least squares at SN positions)
    A = np.column_stack([r, rhat.T]); coef, *_ = np.linalg.lstsq(A, vcf, rcond=None)
    rest = vcf - A@coef
    print(f"\nat {len(idx)} SN positions (0.01<z<{zmax}): a={coef[0]:+.3f} km/s per Mpc/h (dH0={coef[0]*100:+.1f}%), bulk |B|={np.linalg.norm(coef[1:]):.0f} km/s, "
          f"rms: total={np.sqrt(np.mean(vcf**2)):.0f}, after monopole+dipole={np.sqrt(np.mean(rest**2)):.0f} km/s; corr(rest, 2M++)={np.corrcoef(rest, vp2)[0,1]:.3f}")
    A2 = np.column_stack([r]); c2m, *_ = np.linalg.lstsq(A2, vcf, rcond=None); nomono = vcf - A2@c2m
    A3 = np.column_stack([r, rhat.T]); c3, *_ = np.linalg.lstsq(A3, vp2, rcond=None)
    print(f"for reference 2M++ VPEC at the same positions: a={c3[0]:+.3f}, |B|={np.linalg.norm(c3[1:]):.0f} km/s, rms rest={np.sqrt(np.mean((vp2-A3@c3)**2)):.0f}")
    y = d.m_b_corr.values[idx] - mu(zc)
    Lc = np.linalg.cholesky(cov[np.ix_(idx, idx)]); w = lambda X: np.linalg.solve(Lc, X)
    def fit(cols, yy=y):
        X = np.column_stack([np.ones(len(yy))] + cols); Xw, yw = w(X), w(yy)
        b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); rr = yw - Xw@b
        return b, np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))), rr@rr
    _, _, c0 = fit([])
    print(f"\nbeta fits (chi2 none={c0:.1f}):")
    for name, v in [('2M++', vp2), ('CF4++ full', vcf), ('CF4++ minus monopole', nomono), ('CF4++ minus mono+dipole', rest), ('CF4++ dipole only', A@np.r_[0, coef[1:]]), ('CF4++ monopole only', A@np.r_[coef[0], 0, 0, 0])]:
        cc = mu(zhd(zc, v)) - mu(zc)
        _, _, cf = fit([], y - cc); b, e, cb = fit([cc])
        print(f"  {name:26s} beta=1 chi2={cf:.1f} ({c0-cf:+.1f}) | beta={b[1]:.2f}±{e[1]:.2f} chi2={cb:.1f}")

main(*sys.argv[1:4])
