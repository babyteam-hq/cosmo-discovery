"""Grid-convention check for the CF4++ velocity grid: for every signed permutation of the three velocity components
(6 orders x 8 signs = 48), the correlation of the field's radial velocity at the Pantheon+ SN positions with the 2M++
VPEC, the fitted beta of the SN residuals on the whole field, and the Galactic direction of the field's mean 3-vector
at the SNe. The documented convention (example script: v_mean[:, ix, iy, iz] = (vx, vy, vz) in supergalactic cartesian)
is the identity. If another convention matched the SNe or 2M++ clearly better, finding.md's rejection of the non-dipole
structure would be a convention artefact. usage: python cf4_axes.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/new_CF4pp_mean_std_grids.npz
"""
import sys, itertools, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.ndimage import map_coordinates
from astropy.coordinates import SkyCoord, CartesianRepresentation
import astropy.units as u
C = 299792.458; OM = 0.334; N = 128; L = 1000.0; DELTA = L/N
def dc_h(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3 + 1-OM)
    return C/100*np.array([quad(e, 0, zi)[0] for zi in np.atleast_1d(z)])
def mu(z): return 5*np.log10((1+z)*dc_h(z)/0.7) + 25
def zhd(zcmb, v): return (1+zcmb)/(1+v/C) - 1

def main(dat, covf, grid):
    V = np.load(grid)['v_mean']; d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f: n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    m = (d.IS_CALIBRATOR == 0) & (d.zCMB > 0.01) & (d.zCMB < 0.1); idx = np.where(m)[0]
    zc = d.zCMB.values[idx]; vp2 = d.VPEC.values[idx]; r = dc_h(zc); y = d.m_b_corr.values[idx] - mu(zc)
    sc = SkyCoord(ra=d.RA.values[idx]*u.deg, dec=d.DEC.values[idx]*u.deg, distance=r*u.Mpc, frame='fk5', equinox='J2000').supergalactic
    xyz = np.stack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value]); rhat = xyz/np.linalg.norm(xyz, axis=0)
    ci = (xyz + L/2)/DELTA - 0.5
    v3 = np.stack([map_coordinates(V[k], ci, order=1, mode='nearest') for k in range(3)])
    Lc = np.linalg.cholesky(cov[np.ix_(idx, idx)]); w = lambda X: np.linalg.solve(Lc, X)
    def beta(v):
        cc = mu(zhd(zc, v)) - mu(zc); X = np.column_stack([np.ones(len(y)), cc]); Xw, yw = w(X), w(y)
        b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); e = np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))); rr = yw - Xw@b
        X0w = w(np.ones((len(y), 1))); b0, *_ = np.linalg.lstsq(X0w, yw, rcond=None); r0 = yw - X0w@b0
        return b[1], e[1], r0@r0 - rr@rr
    rows = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product([1, -1], repeat=3):
            vv = np.stack([signs[k]*v3[perm[k]] for k in range(3)]); vr = np.sum(vv*rhat, 0)
            b, e, dc = beta(vr); corr = np.corrcoef(vr, vp2)[0, 1]; mean = vv.mean(axis=1)
            g = SkyCoord(CartesianRepresentation(*mean), frame='supergalactic').galactic
            rows.append((dc, perm, signs, corr, b, e, np.linalg.norm(mean), g.l.deg, g.b.deg))
    rows.sort(key=lambda t: -t[0])
    print("signed permutations of (vx,vy,vz), sorted by dchi2 of the SN-residual fit on the whole field:")
    for dc, perm, signs, corr, b, e, amp, l, bb in rows[:8] + [t for t in rows if t[1] == (0, 1, 2) and t[2] == (1, 1, 1)]:
        tag = "DOCUMENTED" if perm == (0, 1, 2) and signs == (1, 1, 1) else ""
        print(f"  perm={perm} signs={signs} corr(vr,2M++)={corr:5.2f} beta={b:5.2f}±{e:4.2f} dchi2={dc:6.1f} mean|v|={amp:4.0f} toward l={l:5.1f} b={bb:5.1f} {tag}")
main(*sys.argv[1:4])
