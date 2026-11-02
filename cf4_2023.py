"""Same test on the 2023 CF4 grid (Courtois+2023, CF4_new_64-z008_velocity.fits, x52 for km/s, axes (SGZ,SGY,SGX)).
The box size isn't in the header, so scan L and axis order, keep the convention that best matches 2M++ VPEC, and report beta,
plus the dipole / rest split under it.
usage: python cf4_2023.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/CF4_new_64-z008_velocity.fits
"""
import sys, itertools, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.ndimage import map_coordinates
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
C = 299792.458; OM = 0.334; N = 64
def dc_h(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3 + 1-OM)
    return C/100*np.array([quad(e, 0, zi)[0] for zi in np.atleast_1d(z)])
def mu(z): return 5*np.log10((1+z)*dc_h(z)/0.7) + 25
def zhd(zcmb, v): return (1+zcmb)/(1+v/C) - 1

def main(dat, covf, fitsf, zmax=0.08):
    V = fits.getdata(fitsf).astype(float)*52.0   # (3, a, b, c)
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f:
        n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    m = (d.IS_CALIBRATOR == 0) & (d.zCMB > 0.01) & (d.zCMB < zmax); idx = np.where(m)[0]
    zc = d.zCMB.values[idx]; vp2 = d.VPEC.values[idx]; r = dc_h(zc)
    sc = SkyCoord(ra=d.RA.values[idx]*u.deg, dec=d.DEC.values[idx]*u.deg, distance=r*u.Mpc, frame='fk5', equinox='J2000').supergalactic
    xyz0 = np.stack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value]); rhat = xyz0/np.linalg.norm(xyz0, axis=0)
    y = d.m_b_corr.values[idx] - mu(zc)
    Lc = np.linalg.cholesky(cov[np.ix_(idx, idx)]); w = lambda X: np.linalg.solve(Lc, X)
    def fit(cols, yy=y):
        X = np.column_stack([np.ones(len(yy))] + cols); Xw, yw = w(X), w(yy)
        b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); rr = yw - Xw@b
        return b, np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))), rr@rr
    rows = []
    for L in [400.0, 480.0, 500.0, 640.0, 1000.0]:
        for perm in itertools.permutations(range(3)):       # which grid array axis each of (x,y,z) indexes
            for vperm in [(0, 1, 2), (2, 1, 0)]:            # which component is vx,vy,vz
                ci = (xyz0[list(perm)] + L/2)/(L/N) - 0.5
                comps = [map_coordinates(V[vperm[k]], ci, order=1, mode='nearest') for k in range(3)]
                # comps are velocity components ordered as the grid's position axes under perm; map back to x,y,z
                vxyz = np.empty((3, len(idx)))
                for k in range(3): vxyz[perm[k]] = comps[k]
                vr = np.sum(vxyz*rhat, 0)
                rows.append((np.corrcoef(vr, vp2)[0, 1], L, perm, vperm, vr))
    rows.sort(key=lambda t: -t[0])
    print("best conventions (corr with 2M++ VPEC):")
    for c, L, perm, vperm, vr in rows[:6]:
        print(f"  L={L:6.0f} pos-perm={perm} v-perm={vperm} corr={c:+.3f} rms={np.sqrt(np.mean(vr**2)):.0f}")
    c, L, perm, vperm, vr = rows[0]
    _, _, c0 = fit([])
    A = np.column_stack([r, rhat.T]); coef, *_ = np.linalg.lstsq(A, vr, rcond=None); rest = vr - A@coef; dip = A@np.r_[0, coef[1:]]
    print(f"\nCF4-2023 at {len(idx)} SN positions (z<{zmax}): bulk |B|={np.linalg.norm(coef[1:]):.0f} km/s, rms total={np.sqrt(np.mean(vr**2)):.0f}, rest={np.sqrt(np.mean(rest**2)):.0f}; chi2 none={c0:.1f}")
    for name, v in [('2M++', vp2), ('CF4-2023 full', vr), ('CF4-2023 dipole only', dip), ('CF4-2023 minus mono+dipole', rest)]:
        cc = mu(zhd(zc, v)) - mu(zc); _, _, cf = fit([], y - cc); b, e, cb = fit([cc])
        print(f"  {name:28s} beta=1 chi2={cf:.1f} ({c0-cf:+.1f}) | beta={b[1]:.2f}±{e[1]:.2f} chi2={cb:.1f}")

main(*sys.argv[1:4])
