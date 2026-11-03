"""R-168 (cosmo-cf4-rev): shared loaders for the calibrated-null revision of the CF4 supernova test.

Samples (each a dict: name, ra, dec, zc, y, Lc (Cholesky of the covariance) or wdiag (1/sigma), nuis (nuisance columns),
plus sample-specific extras), the fields under test (CF4++ posterior mean, 2M++ Carrick+2015 grid) read at the SN
positions, the dipole/rest split and the joint GLS fit. Data in data/ (gitignored):
  Pantheon+SH0ES.dat, Pantheon+SH0ES_STATONLY.cov, pplus_sys/<GROUP>.cov (PantheonPlusSH0ES/DataRelease, sytematic_groupings)
  des5yr_repo/4_DISTANCES_COVMAT (DES-SN5YR Dovekie), ztf/ztfsniadr2_lite/tables/snia_data.csv (ZTF SN Ia DR2 lite)
  new_CF4pp_mean_std_grids.npz (CF4++), twompp/twompp_velocity.npy (2M++), cf4_table2.tsv (CF4 galaxies, VizieR J/ApJ/944/94)
"""
import os, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.ndimage import map_coordinates
from astropy.coordinates import SkyCoord
import astropy.units as u

C = 299792.458; OM = 0.334
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
_zg = np.linspace(0, 0.2, 2001)
_dg = C/100*np.array([quad(lambda x: 1/np.sqrt(OM*(1+x)**3 + 1-OM), 0, z)[0] for z in _zg])
dc_h = interp1d(_zg, _dg, kind='cubic')                  # comoving distance, Mpc/h
def mu(z): return 5*np.log10((1+z)*dc_h(z)/0.7) + 25
def zhd(zcmb, v): return (1+zcmb)/(1+v/C) - 1
def corr_mag(zc, v): return mu(zhd(zc, v)) - mu(zc)     # distance-modulus correction of a radial velocity v (km/s)
SYS_GROUPS = ['BS21', 'C11', 'CALIB', 'CSP_RECAL', 'G10', 'HSTCALSPEC', 'MASSSTEP', 'MWCOLORLAW', 'MWEBV', 'SCATTERMOD',
              'SPECEFF', 'VPEC', 'ZSHIFT']

def readcov(f):
    with open(f) as fh:
        n = int(fh.readline()); return np.loadtxt(fh).reshape(n, n)

def sg_xyz(ra, dec, r):
    sc = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, distance=r*u.Mpc, frame='fk5', equinox='J2000').supergalactic
    return np.stack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value])

def gal_xyz(ra, dec, r):
    sc = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, distance=r*u.Mpc, frame='fk5', equinox='J2000').galactic
    return np.stack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value])

def pantheon(zmax=0.1, drop_vpec_sys=True):
    """Pantheon+ 0.01 < zCMB < zmax, calibrators out. Covariance, when drop_vpec_sys: STAT+SYS minus the VPEC systematic
    (the 2M++-based peculiar-velocity systematic). The grouping files are not STATONLY + one systematic: each carries a common
    baseline that differs from STATONLY (min eigenvalue of G - STATONLY is -0.09 even for negligible systematics), so the VPEC
    term is taken as G_VPEC minus the mean of the seven groupings whose systematic is negligible (CSP_RECAL, HSTCALSPEC,
    MASSSTEP, MWCOLORLAW, MWEBV, SPECEFF, ZSHIFT: spread 2 % among them). At 0.01<z<0.1 that term is 38 % of SYS (Frobenius),
    0.019 mag on the diagonal. The diagonal 250 km/s velocity-noise term inside STAT (m_b_corr_err_VPEC) is kept: it is
    field-agnostic noise, not a 2M++ term."""
    d = pd.read_csv(os.path.join(D, 'Pantheon+SH0ES.dat'), sep=r'\s+')
    m = (d.IS_CALIBRATOR == 0) & (d.zCMB > 0.01) & (d.zCMB < zmax); idx = np.where(m)[0]
    ix = np.ix_(idx, idx)
    C_ = readcov(os.path.join(D, 'Pantheon+SH0ES_STAT+SYS.cov'))[ix]
    if drop_vpec_sys:
        base = ['CSP_RECAL', 'HSTCALSPEC', 'MASSSTEP', 'MWCOLORLAW', 'MWEBV', 'SPECEFF', 'ZSHIFT']
        ref = np.mean([readcov(os.path.join(D, 'pplus_sys', g + '.cov'))[ix] for g in base], 0)
        C_ = C_ - (readcov(os.path.join(D, 'pplus_sys', 'VPEC.cov'))[ix] - ref)
    zc = d.zCMB.values[idx]
    ev = d.m_b_corr_err_VPEC.values[idx]
    return dict(name='Pantheon+', ra=d.RA.values[idx], dec=d.DEC.values[idx], zc=zc, y=d.m_b_corr.values[idx] - mu(zc),
                cov=C_, Lc=np.linalg.cholesky(C_), nuis=[np.ones(len(idx))], cid=d.CID.values[idx], vpec=d.VPEC.values[idx],
                err_vpec=ev, host_ra=d.HOST_RA.values[idx], host_dec=d.HOST_DEC.values[idx], zhd=d.zHD.values[idx])

def des():
    dd = os.path.join(D, 'des5yr_repo', '4_DISTANCES_COVMAT')
    rd = lambda f: pd.read_csv(os.path.join(dd, f), comment='#', sep=r'\s+').drop(columns=['VARNAMES:'], errors='ignore')
    hd = rd('DES-Dovekie_HD.csv'); md = rd('DES-Dovekie_Metadata.csv').set_index('CID')
    z = np.load(os.path.join(dd, 'STAT+SYS.npz')); n = int(z['nsn'][0]); inv = np.zeros((n, n)); inv[np.triu_indices(n)] = z['cov']
    inv = inv + inv.T - np.diag(np.diag(inv)); cov = np.linalg.inv(inv)
    p = pd.read_csv(os.path.join(D, 'Pantheon+SH0ES.dat'), sep=r'\s+').drop_duplicates('CID').set_index('CID')
    hd['zCMB'] = md.loc[hd.CID, 'zCMB'].values
    hd['RA'] = p.RA.reindex(hd.CID).values; hd['DEC'] = p.DEC.reindex(hd.CID).values
    idx = np.where((hd.zCMB > 0.01) & (hd.zCMB < 0.1) & hd.RA.notna())[0]
    zc = hd.zCMB.values[idx]; C_ = cov[np.ix_(idx, idx)]
    return dict(name='DES', ra=hd.RA.values[idx], dec=hd.DEC.values[idx], zc=zc, y=hd.MU.values[idx] - mu(zc), cov=C_,
                Lc=np.linalg.cholesky(C_), nuis=[np.ones(len(idx))], cid=hd.CID.values[idx])

AL, BE = 0.15, 3.1
def ztf(zmax=0.1):
    """ZTF SN Ia DR2 lite, as ztf_tests.py: snia-cosmo, good coverage and fit, host spectroscopic z, |x1|<3, |c|<0.3,
    0.015 < z < zmax; diagonal errors with sigma_int set for chi2/dof = 1 on the base model; alpha, beta free (nuisance)."""
    t = pd.read_csv(os.path.join(D, 'ztf', 'ztfsniadr2_lite', 'tables', 'snia_data.csv'))
    sel = (t.sn_type == 'snia-cosmo') & (t.lccoverage_flag == 1) & (t.fitquality_flag == 1) & (t.source == 'z_gal') \
        & (t.x1.abs() < 3) & (t.c.abs() < 0.3) & (t.redshift > 0.015) & (t.redshift < 0.1) & t.x0.notna()
    t = t[sel].reset_index(drop=True)
    ct = SkyCoord(t.ra.values*u.deg, t.dec.values*u.deg, frame='icrs'); apex = SkyCoord(l=264.021*u.deg, b=48.253*u.deg, frame='galactic')
    zc = (1+t.redshift.values)/(1 - 369.82*np.cos(ct.galactic.separation(apex).rad)/C) - 1
    keep = zc < zmax; t = t[keep].reset_index(drop=True); zc = zc[keep]
    mb = -2.5*np.log10(t.x0.values); x1 = t.x1.values; cc = t.c.values
    var = (2.5/np.log(10)*t.x0_err.values/t.x0.values)**2 + AL**2*t.x1_err.values**2 + BE**2*t.c_err.values**2 \
        + 2*(-AL)*(-2.5/np.log(10)/t.x0.values)*t.cov_x0_x1.values + 2*BE*(-2.5/np.log(10)/t.x0.values)*t.cov_x0_c.values \
        + 2*(-AL)*BE*t.cov_x1_c.values + (5/np.log(10)*300/(C*zc))**2 + (5/np.log(10)*t.redshift_err.values/zc)**2
    y = mb - mu(zc); nuis = [np.ones(len(y)), -x1, cc]
    lo, hi = 0.0, 0.5
    for _ in range(40):
        s = (lo+hi)/2; w = 1/np.sqrt(var + s**2); X = np.column_stack(nuis)*w[:, None]
        b, *_ = np.linalg.lstsq(X, y*w, rcond=None); r = y*w - X@b
        lo, hi = (s, hi) if r@r > len(y) - 3 else (lo, s)
    sig = np.sqrt(var + ((lo+hi)/2)**2)
    return dict(name='ZTF', ra=t.ra.values, dec=t.dec.values, zc=zc, y=y, sig=sig, nuis=nuis, sig_int=(lo+hi)/2)

# ---------------------------------------------------------------- fields at SN positions (radial velocity, km/s, CMB frame)
_F = {}
def cf4pp_at(ra, dec, zc, key='v_mean'):
    if 'cf4pp' not in _F: _F['cf4pp'] = np.load(os.path.join(D, 'new_CF4pp_mean_std_grids.npz'))[key]
    V = _F['cf4pp']; N, L = 128, 1000.0
    xyz = sg_xyz(ra, dec, dc_h(zc)); rhat = xyz/np.linalg.norm(xyz, axis=0); ci = (xyz + L/2)/(L/N) - 0.5
    return np.sum(np.stack([map_coordinates(V[k], ci, order=1, mode='nearest') for k in range(3)])*rhat, 0)

VEXT_2MPP = np.array([89., -131., 17.])  # Carrick+2015 table 3, Galactic Cartesian, km/s (already added inside the grid)
def twompp_at(ra, dec, zc, order='kji'):
    """Carrick+2015 grid (257^3, Galactic Cartesian, cell centres -200..200 Mpc/h, 'i running fastest' => array [c,k,j,i]),
    real-space, beta*=0.43 with V_ext added. Outside the grid (r > 200 Mpc/h on an axis) the field is V_ext alone."""
    if 'twompp' not in _F: _F['twompp'] = np.load(os.path.join(D, 'twompp', 'twompp_velocity.npy'), mmap_mode='r')
    V = _F['twompp']
    xyz = gal_xyz(ra, dec, dc_h(zc)); rhat = xyz/np.linalg.norm(xyz, axis=0)
    ci = (xyz + 200.0)/(400.0/256)
    if order == 'kji': ci = ci[::-1]
    inside = np.all((xyz > -200) & (xyz < 200), axis=0)
    v = np.stack([map_coordinates(np.asarray(V[k]), ci, order=1, mode='nearest') for k in range(3)])
    v[:, ~inside] = VEXT_2MPP[:, None]
    return np.sum(v*rhat, 0)

# ---------------------------------------------------------------- statistic
def split(ra, dec, zc, v):
    """radial velocity v = a r + B.rhat + rest (least squares at the SN positions): returns (dipole part, rest, B)."""
    r = dc_h(zc); rad, de = np.radians(ra), np.radians(dec)
    nh = np.stack([np.cos(de)*np.cos(rad), np.cos(de)*np.sin(rad), np.sin(de)], -1)
    A = np.column_stack([r, nh]); coef, *_ = np.linalg.lstsq(A, v, rcond=None)
    return A@np.r_[0, coef[1:]], v - A@coef, coef[1:]

def whitener(s):
    if 'Lc' in s:
        from scipy.linalg import solve_triangular
        return lambda X: solve_triangular(s['Lc'], X, lower=True)
    w = 1/s['sig']; return lambda X: X*w[:, None] if X.ndim == 2 else X*w

def gls(s, cols, y=None, W=None):
    """GLS of y on [nuisance columns, cols]; returns (coefficients of cols, their errors, chi2, full cov of cols)."""
    y = s['y'] if y is None else y; W = W or whitener(s)
    X = np.column_stack(s['nuis'] + list(cols)); Xw, yw = W(X), W(y)
    b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); rr = yw - Xw@b; cv = np.linalg.inv(Xw.T@Xw); k = len(s['nuis'])
    return b[k:], np.sqrt(np.diag(cv))[k:], rr@rr, cv[k:, k:]

def beta_fits(s, v, y=None, W=None):
    """the statistic: joint GLS of beta_dipole and beta_rest (both columns at once), and the whole-field beta."""
    dip, rest, B = split(s['ra'], s['dec'], s['zc'], v)
    cd, cr = corr_mag(s['zc'], dip), corr_mag(s['zc'], rest)
    b, e, chi, cv = gls(s, [cd, cr], y, W)
    bw, ew, _, _ = gls(s, [corr_mag(s['zc'], v)], y, W)
    # predictive statistic: whitened correlation of the residual with c_rest, both after the nuisance and dipole columns
    W = W or whitener(s); y = s['y'] if y is None else y
    X0 = W(np.column_stack(s['nuis'] + [cd])); Q, _ = np.linalg.qr(X0)
    ry = W(y) - Q@(Q.T@W(y)); rc = W(cr) - Q@(Q.T@W(cr))
    return dict(bd=b[0], ed=e[0], br=b[1], er=e[1], rho=cv[0, 1]/np.sqrt(cv[0, 0]*cv[1, 1]), bw=bw[0], ew=ew[0], B=B,
                corr_rest=float(ry@rc/np.sqrt((ry@ry)*(rc@rc))))
