"""What bulk flow do the Pantheon+ SNe themselves measure, and how does it compare with the CF4++ and 2M++ fields' mean
bulk flow over the same positions? Free 3-vector B (plus monopole and intercept) fitted to residuals by GLS, in shells.
usage: python bulkflow_free.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/new_CF4pp_mean_std_grids.npz
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
R_gal = np.array([[-0.0548755604,-0.8734370902,-0.4838350155],[0.4941094279,-0.4448296300,0.7469822445],[-0.8676661490,-0.1980763734,0.4559837762]])
def lb(v):
    g = R_gal @ v; return np.degrees(np.arctan2(g[1], g[0])) % 360, np.degrees(np.arcsin(g[2]/np.linalg.norm(g)))

def main(dat, covf, grid):
    g = np.load(grid); V = g['v_mean']
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f:
        n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    print("Pantheon+: free bulk flow from the residuals (equatorial 3-vector, GLS, intercept and monopole free) vs the fields' mean 3-vector over the same SNe")
    for lo, hi in [(0.01, 0.1), (0.01, 0.03), (0.03, 0.06), (0.06, 0.1), (0.01, 0.05), (0.05, 0.1)]:
        m = (d.IS_CALIBRATOR == 0) & (d.zCMB > lo) & (d.zCMB < hi); idx = np.where(m)[0]
        zc = d.zCMB.values[idx]; r = dc_h(zc)
        ra, dec = np.radians(d.RA.values[idx]), np.radians(d.DEC.values[idx])
        nh = np.stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)], -1)
        y = d.m_b_corr.values[idx] - mu(zc)
        # dmu from a radial velocity v at low z: dmu ~ -(5/ln10) v/(c z) (sign: an SN moving away looks farther than it is -> z too high -> mu_model too big -> y negative). use exact mapping via derivative
        f = (5/np.log(10))/(C*zc)   # mag per km/s
        X = np.column_stack([np.ones(len(y)), r*f, -f*nh[:, 0], -f*nh[:, 1], -f*nh[:, 2]])
        Lc = np.linalg.cholesky(cov[np.ix_(idx, idx)]); Xw, yw = np.linalg.solve(Lc, X), np.linalg.solve(Lc, y)
        b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); e = np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw)))
        B = b[2:]; eB = e[2:]
        # fields' mean 3-vectors at the SN positions
        sc = SkyCoord(ra=d.RA.values[idx]*u.deg, dec=d.DEC.values[idx]*u.deg, distance=r*u.Mpc, frame='fk5', equinox='J2000')
        sg = sc.supergalactic; xyz = np.stack([sg.cartesian.x.value, sg.cartesian.y.value, sg.cartesian.z.value])
        ci = (xyz + L/2)/DELTA - 0.5
        vsg = np.stack([map_coordinates(V[k], ci, order=1, mode='nearest') for k in range(3)])   # supergalactic components
        # convert the mean SG vector to equatorial via astropy
        msg = vsg.mean(1)
        from astropy.coordinates import CartesianRepresentation, Supergalactic
        vec = SkyCoord(Supergalactic(CartesianRepresentation(msg[0], msg[1], msg[2], unit=u.Mpc))).transform_to('fk5')
        Bcf = np.array([vec.cartesian.x.value, vec.cartesian.y.value, vec.cartesian.z.value])
        # 2M++ : only the radial component is given (VPEC); its best-fit 3-vector at these positions
        vp = d.VPEC.values[idx]; A = np.column_stack([r, nh]); c2, *_ = np.linalg.lstsq(A, vp, rcond=None); B2 = c2[1:]
        cos = lambda a, bb: a@bb/np.linalg.norm(a)/np.linalg.norm(bb)
        print(f"z[{lo},{hi}) n={len(idx):4d}: SNe B={np.linalg.norm(B):4.0f}±{np.mean(eB):3.0f} km/s (l,b)=({lb(B)[0]:.0f},{lb(B)[1]:.0f}) "
              f"| CF4++ mean B={np.linalg.norm(Bcf):4.0f} (l,b)=({lb(Bcf)[0]:.0f},{lb(Bcf)[1]:.0f}) cos={cos(B,Bcf):+.2f} "
              f"| 2M++ fit B={np.linalg.norm(B2):4.0f} (l,b)=({lb(B2)[0]:.0f},{lb(B2)[1]:.0f}) cos={cos(B,B2):+.2f}")

main(*sys.argv[1:4])
