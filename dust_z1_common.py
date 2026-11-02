"""Shared pieces of the R-150 halo-dust reddening stack (DESI DR1 ELG/LRG foregrounds x DESI DR1 z > 2 quasars).

Everything here is fixed before any stack (results-R150.md, M1): the quasar cuts, the colour residuals, the extinction
curves and the per-quasar A_V estimator, the bins, and the Menard et al. 2010 reference profile.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, 'data', 'dust')

# --- bins (M1) -----------------------------------------------------------------------------------------------------
RP_EDGES = np.geomspace(20., 3000., 9)          # proper kpc (h = 0.6766, Planck 2018)
ZBINS = {'ELG': [(0.8, 1.1), (1.1, 1.3), (1.3, 1.6)], 'LRG': [(0.4, 0.6), (0.6, 0.8), (0.8, 1.1)]}
DZ_BEHIND = 0.5                                   # pair kept only if z_q > z_f + 0.5
ZQ_MIN = 2.0
COLOURS = ('g-r', 'r-z', 'z-W1')
LAM = {'g': 4800., 'r': 6400., 'z': 9200., 'W1': 34000., 'W2': 46000.}   # effective wavelengths, Angstrom
WINSOR = 1.0                                      # colour residuals clipped at +-1 mag
H = 0.6766


def cosmo():
    from astropy.cosmology import Planck18
    return Planck18


def da_kpc(z):
    """Angular-diameter distance in proper kpc per radian."""
    return cosmo().angular_diameter_distance(z).to_value('kpc')


# --- extinction curves ---------------------------------------------------------------------------------------------
def k_curve(lam_rest_A, law='SMC'):
    """A_lambda / A_V at rest-frame wavelength (Angstrom). SMC bar (Gordon+2003), primary; MW O'Donnell 1994 R_V = 3.1."""
    import astropy.units as u
    from dust_extinction.averages import G03_SMCBar
    from dust_extinction.parameter_averages import O94
    lam = np.atleast_1d(np.asarray(lam_rest_A, float))
    m = G03_SMCBar() if law == 'SMC' else O94(Rv=3.1)
    x = 1e4 / lam
    lo, hi = (0.3, 10.0) if law == 'SMC' else (0.3, 10.0)
    out = np.zeros_like(lam)
    ok = (x >= lo) & (x <= hi)
    out[ok] = m(lam[ok] * u.AA)
    out[x < lo] = m(1e4 / lo * u.AA) * (x[x < lo] / lo) ** 1.7    # IR power law below the curve's range
    return out


def response(zf, law='SMC'):
    """E(colour)_obs / A_V,rest for dust at redshift zf, for COLOURS."""
    k = {b: k_curve(LAM[b] / (1 + zf), law)[0] for b in ('g', 'r', 'z', 'W1')}
    return np.array([k['g'] - k['r'], k['r'] - k['z'], k['z'] - k['W1']])


def band_k(zf, law='SMC'):
    return {b: k_curve(LAM[b] / (1 + zf), law)[0] for b in LAM}


# --- Menard et al. 2010 reference ------------------------------------------------------------------------------------
def menard_av_rest(rp_kpc, law='MW'):
    """Menard+2010 eq. 28, <E(g-i)>(theta) = 1.4e-3 (theta/1')^-0.84 (SDSS i < 21 galaxies, z_eff = 0.36), converted
    to rest-frame A_V at the galaxy with the stated law; theta -> proper r_p (kpc, input) at z = 0.36 with their
    Omega_m = 0.3 and h = 0.6766 (their 100 h^-1 kpc = 148 kpc)."""
    from astropy.cosmology import FlatLambdaCDM
    z0 = 0.36
    c = FlatLambdaCDM(H0=100 * H, Om0=0.3)
    kpc_per_arcmin = c.angular_diameter_distance(z0).to_value('kpc') * np.pi / 180 / 60
    theta = np.asarray(rp_kpc) / kpc_per_arcmin
    egi = 1.4e-3 * theta ** -0.84
    kg, ki = k_curve(4686. / (1 + z0), law)[0], k_curve(7480. / (1 + z0), law)[0]
    return egi / (kg - ki)


def menard_amplitude(law='MW', rp0=100.):
    """rest-frame A_V of the Menard profile at rp0 proper kpc."""
    return float(menard_av_rest(rp0, law))


# --- quasars ---------------------------------------------------------------------------------------------------------
def mags(q):
    """Galactic-extinction-corrected AB magnitudes from the Legacy Surveys DR9 fluxes in the DESI catalogue."""
    out = {}
    for b, f, t in (('g', 'flux_g', 'mw_transmission_g'), ('r', 'flux_r', 'mw_transmission_r'),
                    ('z', 'flux_z', 'mw_transmission_z'), ('W1', 'flux_w1', 'mw_transmission_w1'),
                    ('W2', 'flux_w2', 'mw_transmission_w2')):
        fl = q[f] / q[t]
        with np.errstate(divide='ignore', invalid='ignore'):
            out[b] = np.where(fl > 0, 22.5 - 2.5 * np.log10(np.clip(fl, 1e-30, None)), np.nan)
    return out


def load_qso(zmin=None):
    q = dict(np.load(os.path.join(D, 'qso.npz')))
    m = mags(q)
    good = np.ones(len(q['z']), bool)
    for b in ('g', 'r', 'z', 'W1'):
        good &= np.isfinite(m[b])
    good &= (q['flux_ivar_w1'] > 0) & (q['flux_w1'] * np.sqrt(q['flux_ivar_w1']) > 3)   # W1 at > 3 sigma
    good &= (m['r'] > 16.5) & (m['r'] < 23)
    good &= (q['desi_target'] & 4) > 0                                                       # main-survey QSO target
    if zmin is not None:
        good &= q['z'] > zmin
    q = {k: v[good] for k, v in q.items()}
    m = {k: v[good] for k, v in m.items()}
    q['mag'] = m
    return q


def colour_residuals(q, zstep=0.02, nmag=4):
    """Colour residuals against the median at the quasar's redshift (0.02 bins), imaging system (N/S) and r-mag quartile;
    winsorized at +-WINSOR. Returns (n, 3) array and the r-mag quartile index."""
    m = q['mag']
    col = np.c_[m['g'] - m['r'], m['r'] - m['z'], m['z'] - m['W1']]
    res = np.full_like(col, np.nan)
    zb = np.floor(q['z'] / zstep).astype(int)
    qr = np.zeros(len(zb), int)
    for north in (False, True):
        s = q['north'] == north
        edges = np.percentile(m['r'][s], np.linspace(0, 100, nmag + 1))
        qr[s] = np.clip(np.searchsorted(edges, m['r'][s]) - 1, 0, nmag - 1)
    key = (zb * 2 + q['north'].astype(int)) * nmag + qr
    order = np.argsort(key, kind='stable')
    ks = key[order]
    bounds = np.flatnonzero(np.diff(ks)) + 1
    for idx in np.split(order, bounds):
        if len(idx) < 30:
            continue
        res[idx] = col[idx] - np.median(col[idx], axis=0)
    res = np.clip(res, -WINSOR, WINSOR)
    return res, qr


def av_estimator(res, qr, zf, law='SMC'):
    """Per-quasar A_V,rest estimate for dust at zf: GLS projection of the colour residuals on the reddening vector,
    covariance measured from the data in each r-mag quartile. Returns (A_hat, weight = 1/var)."""
    R = response(zf, law)
    ah = np.full(len(res), np.nan)
    w = np.zeros(len(res))
    ok = np.all(np.isfinite(res), axis=1)
    for k in np.unique(qr):
        s = ok & (qr == k)
        C = np.cov(res[s].T)
        Ci = np.linalg.inv(C)
        f = Ci @ R / (R @ Ci @ R)
        ah[s] = res[s] @ f
        w[s] = R @ Ci @ R
    return ah, w


# --- pairs -----------------------------------------------------------------------------------------------------------
def unitvec(ra, dec):
    r, d = np.radians(ra), np.radians(dec)
    return np.c_[np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)]


def pairs(qra, qdec, qz, fra, fdec, fz, zlo, zhi, rmax_kpc=RP_EDGES[-1], dz=DZ_BEHIND, behind=False):
    """All (quasar, foreground) pairs with zlo <= z_f < zhi, r_p (proper, at z_f) < rmax, and z_q > z_f + dz
    (behind=True: the null with the 'foreground' behind the quasar, z_f > z_q + 0.1 instead).
    Returns qi, fi, rp_kpc."""
    from scipy.spatial import cKDTree
    s = np.flatnonzero((fz >= zlo) & (fz < zhi))
    tree = cKDTree(unitvec(fra[s], fdec[s]))
    th_max = rmax_kpc / da_kpc(max(zlo, 0.3))
    chord = 2 * np.sin(th_max / 2)
    qv = unitvec(qra, qdec)
    lists = tree.query_ball_point(qv, chord)
    n = np.fromiter((len(l) for l in lists), int, len(lists))
    qi = np.repeat(np.arange(len(qra)), n)
    fi = s[np.concatenate([np.asarray(l, int) for l in lists])] if n.sum() else np.zeros(0, int)
    cosang = np.sum(qv[qi] * unitvec(fra[fi], fdec[fi]), axis=1)
    th = np.arccos(np.clip(cosang, -1, 1))
    zg = np.linspace(zlo, zhi, 64)
    rp = th * np.interp(fz[fi], zg, da_kpc(zg))
    keep = (rp < rmax_kpc) & (rp >= RP_EDGES[0])
    keep &= (fz[fi] > qz[qi] + 0.1) if behind else (qz[qi] > fz[fi] + dz)
    return qi[keep], fi[keep], rp[keep]


# --- secondary: spectroscopic continuum slope from the Redrock fit ----------------------------------------------------
SLOPE_WIN = ((6000., 6800.), (8800., 9600.))      # observed-frame windows (redward of Lya to z = 3.5)


def redrock_slope(q):
    """-2.5 log10 of the Redrock QSO-HIZ v1.0 model's mean flux density in the blue window over the red window (observed
    frame), from each quasar's COEFF at Z_RR (the DR1/iron template set). NaN where Z_RR and the catalogue z disagree."""
    import fitsio
    T = fitsio.read(os.path.join(D, 'ext', 'rr', 'rrtemplate-QSO-HIZ-v1.0.fits'))
    hdr = fitsio.read_header(os.path.join(D, 'ext', 'rr', 'rrtemplate-QSO-HIZ-v1.0.fits'))
    lam = 10 ** (hdr['CRVAL1'] + hdr['CDELT1'] * np.arange(T.shape[1]))
    c = q['coeff'][:, :4]
    out = np.full(len(c), np.nan)
    ok = (np.abs(q['z_rr'] - q['z']) < 0.01) & (q['ncoeff'] == 4)
    flux = []
    for (lo, hi) in SLOPE_WIN:
        grid = np.linspace(lo, hi, 41)
        # template basis integrated over the window at each quasar's rest frame
        vals = np.zeros((len(c), 4))
        for k in range(4):
            vals[:, k] = np.mean(np.interp(grid[None, :] / (1 + q['z_rr'][:, None]), lam, T[k]), axis=1)
        flux.append(np.sum(vals * c, axis=1))
    with np.errstate(invalid='ignore', divide='ignore'):
        good = ok & (flux[0] > 0) & (flux[1] > 0)
        out[good] = -2.5 * np.log10(flux[0][good] / flux[1][good])
    return out


def slope_residual(q, slope, zstep=0.02):
    res = np.full(len(slope), np.nan)
    key = np.floor(q['z'] / zstep).astype(int) * 2 + q['north'].astype(int)
    order = np.argsort(key, kind='stable'); ks = key[order]
    for idx in np.split(order, np.flatnonzero(np.diff(ks)) + 1):
        f = idx[np.isfinite(slope[idx])]
        if len(f) >= 30:
            res[f] = slope[f] - np.median(slope[f])
    return np.clip(res, -WINSOR, WINSOR)


def slope_response(zf, law='SMC'):
    lb, lr = np.mean(SLOPE_WIN[0]), np.mean(SLOPE_WIN[1])
    return k_curve(lb / (1 + zf), law)[0] - k_curve(lr / (1 + zf), law)[0]
