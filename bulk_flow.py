"""For each low-z SN sample (Pantheon+ 0.01<z<0.1 with full covariance; DES-SN5YR Dovekie low-z with its covariance; ZTF
DR2 with diagonal errors, alpha/beta free), fit a free bulk-flow vector B (3 params, plus the magnitude monopole) to the
Hubble residuals y = m - mu(z_CMB): v_r = B.rhat enters through the linearised distance-modulus shift per km/s at each
SN. Report |B|, its Galactic direction and errors, next to the bulk flow the CF4++ field carries at the same positions
(mean of its 3-vector there, and the bulk vector fitted to its radial projection the way cf4_robust.py splits it) and
the 2M++ bulk vector fitted to the catalogue's radial VPEC (Pantheon+ and DES only). Frames: fits in supergalactic
cartesian, reported in Galactic (l, b).
usage: python bulk_flow.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/des5yr_repo/4_DISTANCES_COVMAT data/ztf/ztfsniadr2_lite/tables/snia_data.csv data/new_CF4pp_mean_std_grids.npz
"""
import sys, os, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.ndimage import map_coordinates
from astropy.coordinates import SkyCoord, CartesianRepresentation
import astropy.units as u
C = 299792.458; OM = 0.334; N = 128; L = 1000.0; DELTA = L/N; AL, BE = 0.15, 3.1
def dc_h(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3 + 1-OM)
    return C/100*np.array([quad(e, 0, zi)[0] for zi in np.atleast_1d(z)])
def mu(z): return 5*np.log10((1+z)*dc_h(z)/0.7) + 25
def zhd(zcmb, v): return (1+zcmb)/(1+v/C) - 1

def sg_to_gal(vec):
    s = SkyCoord(CartesianRepresentation(vec[0], vec[1], vec[2]), frame='supergalactic').galactic
    return s.l.deg, s.b.deg

def describe(B, covB=None):
    amp = np.linalg.norm(B); l, b = sg_to_gal(B)
    if covB is None: return f"|B|={amp:4.0f} km/s toward l={l:5.1f} b={b:5.1f}"
    n = B/amp; e_amp = np.sqrt(n @ covB @ n); e_dir = np.degrees(np.sqrt(max(np.trace(covB) - n @ covB @ n, 0)/2)/amp)
    return f"|B|={amp:4.0f}±{e_amp:3.0f} km/s toward l={l:5.1f} b={b:5.1f} (±{e_dir:3.0f} deg); components SG {np.round(B).astype(int)} ± {np.round(np.sqrt(np.diag(covB))).astype(int)}"

def sample(label, ra, dec, zc, y, whiten, base_cols, vpec2, V):
    r = dc_h(zc)
    sc = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, distance=r*u.Mpc, frame='fk5', equinox='J2000').supergalactic
    xyz = np.stack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value]); rhat = xyz/np.linalg.norm(xyz, axis=0)
    gain = mu(zhd(zc, 1.0)) - mu(zc)          # dmu per km/s of radial velocity
    X = np.column_stack(base_cols + [gain*rhat[k] for k in range(3)]); Xw, yw = whiten(X), whiten(y)
    b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); cov = np.linalg.inv(Xw.T@Xw); k = len(base_cols)
    B = b[k:]; covB = cov[k:, k:]
    X0w = whiten(np.column_stack(base_cols)); b0, *_ = np.linalg.lstsq(X0w, yw, rcond=None)
    dchi2 = (yw - X0w@b0)@(yw - X0w@b0) - (yw - Xw@b)@(yw - Xw@b)
    print(f"{label}: n={len(y)}, z median {np.median(zc):.3f}, depth median {np.median(r):.0f} Mpc")
    print(f"  SN residual bulk flow (free vector, 3 params): {describe(B, covB)}; dchi2={dchi2:.1f} for 3 dof")
    # the fields at these positions
    ci = (xyz + L/2)/DELTA - 0.5
    v3 = np.stack([map_coordinates(V[kk], ci, order=1, mode='nearest') for kk in range(3)])
    print(f"  CF4++ mean 3-vector at the SNe:            {describe(v3.mean(axis=1))}")
    A = rhat.T
    coef, *_ = np.linalg.lstsq(A, np.sum(v3*rhat, 0), rcond=None); print(f"  CF4++ bulk fitted to its radial projection: {describe(coef)}")
    if vpec2 is not None:
        coef2, *_ = np.linalg.lstsq(A, vpec2, rcond=None); print(f"  2M++ bulk fitted to the catalogue VPEC:     {describe(coef2)}")
    return B, covB

def main(pdat, pcov, desdir, ztfcsv, grid):
    V = np.load(grid)['v_mean']
    # Pantheon+
    d = pd.read_csv(pdat, sep=r'\s+')
    with open(pcov) as f: n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    m = (d.IS_CALIBRATOR == 0) & (d.zCMB > 0.01) & (d.zCMB < 0.1); idx = np.where(m)[0]
    zc = d.zCMB.values[idx]; y = d.m_b_corr.values[idx] - mu(zc)
    Lc = np.linalg.cholesky(cov[np.ix_(idx, idx)]); w = lambda X: np.linalg.solve(Lc, X)
    sample("Pantheon+ 0.01<z<0.1", d.RA.values[idx], d.DEC.values[idx], zc, y, w, [np.ones(len(idx))], d.VPEC.values[idx], V)
    # DES Dovekie low-z
    rd = lambda f: pd.read_csv(os.path.join(desdir, f), comment='#', sep=r'\s+').drop(columns=['VARNAMES:'], errors='ignore')
    hd = rd('DES-Dovekie_HD.csv'); md = rd('DES-Dovekie_Metadata.csv').set_index('CID')
    z = np.load(os.path.join(desdir, 'STAT+SYS.npz')); nn = int(z['nsn'][0]); inv = np.zeros((nn, nn)); inv[np.triu_indices(nn)] = z['cov']
    inv = inv + inv.T - np.diag(np.diag(inv)); dcov = np.linalg.inv(inv)
    p = d.drop_duplicates('CID').set_index('CID')
    hd['zCMB'] = md.loc[hd.CID, 'zCMB'].values; hd['VPEC'] = md.loc[hd.CID, 'VPEC'].values
    hd['RA'] = p.RA.reindex(hd.CID).values; hd['DEC'] = p.DEC.reindex(hd.CID).values
    bm = (hd.zCMB > 0.01) & (hd.zCMB < 0.1) & hd.RA.notna(); idx = np.where(bm)[0]
    zc = hd.zCMB.values[idx]; y = hd.MU.values[idx] - mu(zc)
    Lc = np.linalg.cholesky(dcov[np.ix_(idx, idx)]); w = lambda X: np.linalg.solve(Lc, X)
    sample("DES Dovekie low-z 0.01<z<0.1", hd.RA.values[idx], hd.DEC.values[idx], zc, y, w, [np.ones(len(idx))], hd.VPEC.values[idx], V)
    # ZTF DR2
    t = pd.read_csv(ztfcsv)
    sel = (t.sn_type == 'snia-cosmo') & (t.lccoverage_flag == 1) & (t.fitquality_flag == 1) & (t.source == 'z_gal') \
        & (t.x1.abs() < 3) & (t.c.abs() < 0.3) & (t.redshift > 0.015) & (t.redshift < 0.1) & t.x0.notna()
    t = t[sel].reset_index(drop=True)
    ct = SkyCoord(t.ra.values*u.deg, t.dec.values*u.deg, frame='icrs'); apex = SkyCoord(l=264.021*u.deg, b=48.253*u.deg, frame='galactic')
    zc = (1+t.redshift.values)/(1 - 369.82*np.cos(ct.galactic.separation(apex).rad)/C) - 1
    mb = -2.5*np.log10(t.x0.values); x1 = t.x1.values; cc = t.c.values
    var = (2.5/np.log(10)*t.x0_err.values/t.x0.values)**2 + AL**2*t.x1_err.values**2 + BE**2*t.c_err.values**2 \
        + 2*(-AL)*(-2.5/np.log(10)/t.x0.values)*t.cov_x0_x1.values + 2*BE*(-2.5/np.log(10)/t.x0.values)*t.cov_x0_c.values \
        + 2*(-AL)*BE*t.cov_x1_c.values + (5/np.log(10)*300/(C*zc))**2 + (5/np.log(10)*t.redshift_err.values/zc)**2 + 0.078**2
    wz = 1/np.sqrt(var); w = lambda X: X*wz[:, None] if X.ndim == 2 else X*wz
    y = mb - mu(zc)
    for zmax in (0.1, 0.06):
        s = zc < zmax
        wzs = wz[s]; ws = lambda X, wzs=wzs: X*wzs[:, None] if X.ndim == 2 else X*wzs
        sample(f"ZTF DR2 0.015<z<{zmax}", t.ra.values[s], t.dec.values[s], zc[s], y[s], ws, [np.ones(s.sum()), -x1[s], cc[s]], None, V)

main(*sys.argv[1:6])
