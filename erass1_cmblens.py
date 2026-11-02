"""eRASS1 cluster masses confronted with Planck PR3 CMB lensing (R-121, cosmo-deep-1).

Model under test: the eRASS1 primary catalogue's M500 (Bulbul+2024, weak-lensing-calibrated count-rate-mass relation,
selection-corrected posteriors). Data: the Planck 2018 lensing convergence reconstruction from the tSZ-deprojected SMICA map
(COM_Lensing-Szdeproj_4096_R3.00, TT), the MV reconstruction (COM_Lensing_4096_R3.00) as the cross-check. The stacked,
band-limited kappa profile at the clusters is compared with the NFW (Duffy+2008 c-M) plus two-halo (linear P(k) x Tinker+2010
bias) prediction from the catalogue masses, filtered identically; the amplitude A and a mass rescaling b_M (M_true = b_M x
M_eRASS1) are fitted with a covariance from stacks of the catalogue rotated in Galactic longitude (null stacks).

usage:
  python erass1_cmblens.py forecast                       # Fisher S/N of the stack before touching the map
  python erass1_cmblens.py profiles [--rec szdeproj|mv] [--nside 1024] [--lmin 8] [--lmax 2048]
                                                          # per-cluster annulus means -> data/erass1_cmblens/<rec>_profiles.npz
  python erass1_cmblens.py nulls [--rec ...] [--nrand 300] [--seed 1]   # rotated-catalogue stacks (same selection and filter)
  python erass1_cmblens.py fit [--rec ...]                # stacks, covariance, A and b_M fits, splits, trials
  python erass1_cmblens.py ystack [--rec ...] [--nrand 50] # Planck MILCA y in the same annuli (tSZ diagnostic)
All data under data/ (gitignored): data/erass1/erass1cl_primary_v3.2.fits, data/pr3lens/COM_Lensing_Szdeproj_4096_R3.00/,
data/pr3lens/COM_Lensing_4096_R3.00/. Run with python -I.
"""
import argparse, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
OUT = os.path.join(DATA, 'erass1_cmblens')
CAT = os.path.join(DATA, 'erass1', 'erass1cl_primary_v3.2.fits')
REC = {'szdeproj': os.path.join(DATA, 'pr3lens', 'COM_Lensing_Szdeproj_4096_R3.00', 'TT'),
       'mv': os.path.join(DATA, 'pr3lens', 'COM_Lensing_4096_R3.00', 'MV'),
       'pp': os.path.join(DATA, 'pr3lens', 'COM_Lensing_4096_R3.00', 'PP'),
       'tt': os.path.join(DATA, 'pr3lens', 'COM_Lensing_4096_R3.00', 'TT')}
MASK = {'szdeproj': os.path.join(DATA, 'pr3lens', 'COM_Lensing_Szdeproj_4096_R3.00', 'mask.fits.gz'),
        'mv': os.path.join(DATA, 'pr3lens', 'COM_Lensing_4096_R3.00', 'mask.fits.gz'),
        'pp': os.path.join(DATA, 'pr3lens', 'COM_Lensing_4096_R3.00', 'mask.fits.gz'),
        'tt': os.path.join(DATA, 'pr3lens', 'COM_Lensing_4096_R3.00', 'mask.fits.gz')}

# eROSITA cosmology (Ghirardini+2024 best fit; the catalogue masses are computed in it)
H0, OM, OB = 67.77, 0.29, 0.049
ARCMIN = np.pi / 180 / 60
# annuli (arcmin) for the stacked profile
EDGES = np.array([0, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120, 150, 180.0])
THETA_MAX = 3.0  # deg, radius of the disc read per cluster


def cosmo():
    from astropy.cosmology import FlatLambdaCDM
    return FlatLambdaCDM(H0=H0, Om0=OM, Ob0=OB, Tcmb0=2.7255)


# ----------------------------------------------------------------------------------------------------------- catalogue
def load_cat(pcont_max=0.5):
    from astropy.io import fits
    d = fits.open(CAT)[1].data
    ok = (d['PCONT'] < pcont_max) & np.isfinite(d['M500']) & (d['M500'] > 0) & (d['BEST_Z'] > 0.01)
    c = {k: np.asarray(d[k][ok], dtype=float) for k in
         ['RA', 'DEC', 'BEST_Z', 'M500', 'M500_L', 'M500_H', 'R500', 'EXT_LIKE', 'DET_LIKE_0', 'PCONT', 'CR500',
          'CR500_L', 'CR500_H', 'EXP', 'L500', 'KT']}
    c['M500'] *= 1e13; c['M500_L'] *= 1e13; c['M500_H'] *= 1e13
    c['ztype'] = np.asarray(d['BEST_Z_TYPE'][ok])
    c['psz'] = np.array(['PSZ' in s for s in d['MATCH_NAME'][ok]])
    c['act'] = np.array([('ACT' in s) or ('SPT' in s) for s in d['MATCH_NAME'][ok]])
    c['name'] = np.asarray(d['NAME'][ok])
    cosmo_ids = set(fits.open(os.path.join(DATA, 'erass1', 'erass1cl_cosmology_v1.1.fits'))[1].data['DETUID'])
    c['cosmo'] = np.array([u in cosmo_ids for u in d['DETUID'][ok]])  # the 5,259-cluster cosmology sample (~95 % pure)
    # posterior mean mass from the PDF (the catalogue M500 is the median)
    pdf = np.asarray(d['M500_PDF'][ok], dtype=float); grid = np.asarray(d['M500_PDF_array'][ok], dtype=float)
    w = pdf / np.maximum(pdf.sum(1, keepdims=True), 1e-300)
    c['M500_mean'] = (w * grid).sum(1)
    bad = ~np.isfinite(c['M500_mean']) | (c['M500_mean'] <= 0)
    c['M500_mean'][bad] = c['M500'][bad]
    return c


def radec_to_gal(ra, dec):
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    g = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs').galactic
    return g.l.deg, g.b.deg


# ------------------------------------------------------------------------------------------------------- halo model
class HaloModel:
    """NFW + two-halo convergence for a CMB source, per (M500c, z); flat-sky Hankel transform to multipoles."""

    def __init__(self, lmax=3000, cscale=1.0):
        import camb
        self.cscale = cscale
        from scipy.interpolate import RectBivariateSpline
        self.co = cosmo()
        self.lmax = lmax
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=H0, ombh2=OB * (H0 / 100) ** 2, omch2=(OM - OB) * (H0 / 100) ** 2, mnu=0.06, tau=0.054)
        pars.InitPower.set_params(As=2.1e-9, ns=0.965)
        pars.set_matter_power(redshifts=list(np.linspace(0, 1.6, 17))[::-1], kmax=50.0)
        pars.NonLinear = camb.model.NonLinear_none
        res = camb.get_results(pars)
        s8 = res.get_sigma8_0()
        self.s8_scale = (0.88 / s8) ** 2  # eRASS1 best-fit sigma8 = 0.88 sets the two-halo amplitude; see finding
        kh, zs, pk = res.get_matter_power_spectrum(minkh=1e-4, maxkh=50, npoints=600)
        self.h = H0 / 100
        self.k = kh * self.h  # 1/Mpc
        self.zs = np.array(zs)
        self.pk = np.array(pk) / self.h ** 3 * self.s8_scale  # Mpc^3, rows = zs ascending? camb returns ascending z
        order = np.argsort(self.zs); self.zs = self.zs[order]; self.pk = self.pk[order]
        self._pk_spl = RectBivariateSpline(self.zs, np.log(self.k), np.log(self.pk))
        # sigma(M) at z=0 for the bias (top hat, mean matter density)
        self.rho_m0 = (self.co.Om0 * self.co.critical_density0).to('Msun/Mpc3').value
        self.rho_c0 = self.co.critical_density0.to('Msun/Mpc3').value
        self.theta = np.logspace(np.log10(0.05 * ARCMIN), np.log10(6 * np.pi / 180), 900)  # rad
        self.ell = np.arange(0, lmax + 1)
        self._cache = {}

    def plin(self, z, k):
        return np.exp(self._pk_spl(z, np.log(k), grid=False))

    def growth(self, z):
        return np.sqrt(self.plin(z, 0.05) / self.plin(0.0, 0.05))

    def sigma_m(self, M, z):
        R = (3 * M / (4 * np.pi * self.rho_m0)) ** (1 / 3)
        k = self.k; x = k * R
        W = 3 * (np.sin(x) - x * np.cos(x)) / x ** 3
        P = self.plin(z, k)
        return np.sqrt(np.trapezoid(k ** 2 * P * W ** 2, k) / (2 * np.pi ** 2))

    def bias(self, M200c, z):
        """Tinker+2010 bias at Delta = 200 (mean) applied to M200c (a few-% approximation at these masses)."""
        nu = 1.686 / self.sigma_m(M200c, z)
        y = np.log10(200.0)
        A = 1 + 0.24 * y * np.exp(-(4 / y) ** 4); a = 0.44 * y - 0.88
        B = 0.183; b = 1.5
        C = 0.019 + 0.107 * y + 0.19 * np.exp(-(4 / y) ** 4); c = 2.4
        return 1 - A * nu ** a / (nu ** a + 1.686 ** a) + B * nu ** b + C * nu ** c

    def m500_to_nfw(self, M500c, z):
        """Duffy+2008 c200c-M200c (full sample, NFW); solve M200c from M500c; returns (M200c, c200c, r_s [Mpc], rho_s)."""
        from scipy.optimize import brentq
        rho_c = self.rho_c0 * self.co.efunc(z) ** 2
        mfun = lambda x: np.log(1 + x) - x / (1 + x)

        def resid(lnM200):
            M200 = np.exp(lnM200)
            c200 = self.cscale * 5.71 * (M200 * self.h / 2e12) ** -0.084 * (1 + z) ** -0.47
            r200 = (3 * M200 / (4 * np.pi * 200 * rho_c)) ** (1 / 3)
            # find r500: M(<r) = M200 m(c r/r200)/m(c); 500 rho_c 4/3 pi r^3 = M(<r)
            f = lambda r: M200 * mfun(c200 * r / r200) / mfun(c200) - 500 * rho_c * 4 / 3 * np.pi * r ** 3
            r500 = brentq(f, 0.05 * r200, r200)
            return np.log(M200 * mfun(c200 * r500 / r200) / mfun(c200)) - np.log(M500c)

        lnM200 = brentq(resid, np.log(M500c), np.log(3 * M500c))
        M200 = np.exp(lnM200)
        c200 = self.cscale * 5.71 * (M200 * self.h / 2e12) ** -0.084 * (1 + z) ** -0.47
        r200 = (3 * M200 / (4 * np.pi * 200 * rho_c)) ** (1 / 3)
        rs = r200 / c200
        rho_s = M200 / (4 * np.pi * rs ** 3 * mfun(c200))
        return M200, c200, rs, rho_s

    def sigma_crit(self, z):
        from astropy.constants import c, G
        import astropy.units as u
        Dl = self.co.angular_diameter_distance(z); Ds = self.co.angular_diameter_distance(1100.0)
        Dls = self.co.angular_diameter_distance_z1z2(z, 1100.0)
        return (c ** 2 / (4 * np.pi * G) * Ds / (Dl * Dls)).to('Msun/Mpc2').value

    @staticmethod
    def nfw_sigma(R, rs, rho_s):
        x = R / rs
        f = np.empty_like(x)
        lo = x < 1 - 1e-6; hi = x > 1 + 1e-6; eq = ~(lo | hi)
        xl = x[lo]; xh = x[hi]
        f[lo] = (1 - 2 / np.sqrt(1 - xl ** 2) * np.arctanh(np.sqrt((1 - xl) / (1 + xl)))) / (xl ** 2 - 1)
        f[hi] = (1 - 2 / np.sqrt(xh ** 2 - 1) * np.arctan(np.sqrt((xh - 1) / (1 + xh)))) / (xh ** 2 - 1)
        f[eq] = 1 / 3
        return 2 * rs * rho_s * f

    def kappa_theta(self, M500c, z, b_M=1.0):
        """1-halo and 2-halo kappa(theta) on self.theta; b_M rescales the mass before the profile is built."""
        M = b_M * M500c
        M200, c200, rs, rho_s = self.m500_to_nfw(M, z)
        Dl = self.co.angular_diameter_distance(z).value  # Mpc
        scrit = self.sigma_crit(z)
        R = self.theta * Dl
        k1 = self.nfw_sigma(R, rs, rho_s) / scrit
        # two-halo: Sigma_com(R_c) = rho_m0 b int k dk/2pi P(k,z) J0(k R_c); physical Sigma = (1+z)^2 Sigma_com
        from scipy.special import j0
        b = self.bias(M200, z)
        kk = np.logspace(-3.5, 1.3, 1200)
        P = self.plin(z, kk)
        Rc = R * (1 + z)
        integ = np.trapezoid(kk[None, :] * P[None, :] * j0(kk[None, :] * Rc[:, None]), kk, axis=1) / (2 * np.pi)
        k2 = (1 + z) ** 2 * self.rho_m0 * b * integ / scrit
        return k1, k2, dict(M200=M200, c200=c200, bias=b, r500_arcmin=(3 * M / (4 * np.pi * 500 * self.rho_c0 * self.co.efunc(z) ** 2)) ** (1 / 3) / Dl / ARCMIN)

    def to_ell(self, ktheta):
        """flat-sky Hankel transform kappa(theta) -> kappa(ell) on self.ell."""
        from scipy.special import j0
        th = self.theta
        # integrate 2 pi int kappa(theta) J0(l theta) theta dtheta on the log grid
        out = np.empty(len(self.ell))
        chunk = 400
        for i in range(0, len(self.ell), chunk):
            l = self.ell[i:i + chunk]
            J = j0(l[:, None] * th[None, :])
            out[i:i + chunk] = 2 * np.pi * np.trapezoid(ktheta[None, :] * J * th[None, :], th, axis=1)
        return out

    @staticmethod
    def annulus_means_from_ell(kell, ell, filt, edges_arcmin):
        """back-transform a filtered kappa(ell) to the mean in each annulus (area-weighted)."""
        from scipy.special import j1
        e = edges_arcmin * ARCMIN
        l = ell.astype(float); f = kell * filt
        means = np.empty(len(e) - 1)
        for i in range(len(e) - 1):
            a, b = e[i], e[i + 1]
            # mean over annulus of int l dl/2pi f(l) J0(l theta) = int l dl/2pi f(l) [2 (b J1(lb) - a J1(la)) / (l (b^2-a^2))]
            with np.errstate(invalid='ignore', divide='ignore'):
                W = 2 * (b * j1(l * b) - a * j1(l * a)) / (l * (b ** 2 - a ** 2))
            W[l == 0] = 1.0
            means[i] = np.trapezoid(l * f * W, l) / (2 * np.pi)
        return means

    def cluster_kell(self, M500c, z, b_M=1.0):
        key = (round(float(M500c) / 1e12, 3), round(float(z), 4), round(float(b_M), 4))
        if key in self._cache:
            return self._cache[key]
        k1, k2, info = self.kappa_theta(M500c, z, b_M)
        r = (self.to_ell(k1), self.to_ell(k2), info)
        self._cache[key] = r
        return r


def model_grid(hm, b_M=1.0, nM=22, nz=20, tag=''):
    """kappa(ell) 1h and 2h on a (log M500, z) grid, cached to disk; bilinear interpolation for every cluster."""
    os.makedirs(OUT, exist_ok=True)
    fn = os.path.join(OUT, f'modelgrid_bM{b_M:.3f}{tag}.npz')
    if os.path.exists(fn):
        g = np.load(fn)
        return g['lnM'], g['z'], g['k1'], g['k2'], g['bias']
    lnM = np.linspace(np.log(3e12), np.log(3e15), nM)
    zg = np.concatenate([np.linspace(0.015, 0.1, 6), np.linspace(0.13, 1.5, nz - 6)])
    k1 = np.zeros((nM, len(zg), len(hm.ell))); k2 = np.zeros_like(k1); bias = np.zeros((nM, len(zg)))
    t0 = time.time()
    for i, lm in enumerate(lnM):
        for j, z in enumerate(zg):
            a, b, info = hm.cluster_kell(np.exp(lm), z, b_M)
            k1[i, j] = a; k2[i, j] = b; bias[i, j] = info['bias']
        print(f'  grid row {i + 1}/{nM} ({time.time() - t0:.0f}s)', flush=True)
    np.savez(fn, lnM=lnM, z=zg, k1=k1, k2=k2, bias=bias)
    return lnM, zg, k1, k2, bias


def interp_grid(lnM, zg, arr, lm, z):
    """bilinear interpolation of arr[nM, nz, ...] at points (lm, z); linear in ln M and z."""
    i = np.clip(np.searchsorted(lnM, lm) - 1, 0, len(lnM) - 2)
    j = np.clip(np.searchsorted(zg, z) - 1, 0, len(zg) - 2)
    fi = np.clip((lm - lnM[i]) / (lnM[i + 1] - lnM[i]), 0, 1)
    fj = np.clip((z - zg[j]) / (zg[j + 1] - zg[j]), 0, 1)
    # work in log for the 1h (positive); arr may be ~0 at large ell -> use linear
    out = ((1 - fi) * (1 - fj))[..., None] * arr[i, j] + (fi * (1 - fj))[..., None] * arr[i + 1, j] + \
          ((1 - fi) * fj)[..., None] * arr[i, j + 1] + (fi * fj)[..., None] * arr[i + 1, j + 1]
    return out


# ------------------------------------------------------------------------------------------------------- forecast
def forecast(args):
    c = load_cat()
    hm = HaloModel(lmax=3000)
    n = np.loadtxt(os.path.join(REC['szdeproj'], 'nlkk.dat'))
    nl = n[:, 1]; cn = n[:, 2]
    L = np.arange(len(nl))
    lnM, zg, k1, k2, bias = model_grid(hm)
    lm = np.log(c['M500']); z = c['BEST_Z']
    K1 = interp_grid(lnM, zg, k1, lm, z); K2 = interp_grid(lnM, zg, k2, lm, z)
    ell = hm.ell
    cnL = np.interp(ell, L, cn)
    def sn2(K, lmax):
        sel = (ell >= 8) & (ell <= lmax)
        return np.sum(ell[sel] * np.gradient(ell)[sel] / (2 * np.pi) * K[:, sel] ** 2 / cnL[sel], axis=1)
    print(f'clusters {len(z)}, median M500 {np.median(c["M500"]):.2e}, median z {np.median(z):.3f}')
    for lmax in [1000, 2048, 3000]:
        s1 = sn2(K1, lmax); s12 = sn2(K1 + K2, lmax)
        print(f'lmax {lmax}: Fisher S/N 1-halo {np.sqrt(s1.sum()):.1f}, 1h+2h {np.sqrt(s12.sum()):.1f} '
              f'(per-cluster median S/N {np.sqrt(np.median(s1)):.3f})')
    for lo, hi in [(0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.5), (0.5, 2)]:
        s = (z >= lo) & (z < hi)
        print(f'  z {lo}-{hi}: n {s.sum()}, S/N 1h {np.sqrt(sn2(K1[s], 2048).sum()):.1f}')
    for lo, hi in [(0, 1e14), (1e14, 2e14), (2e14, 4e14), (4e14, 1e16)]:
        s = (c['M500'] >= lo) & (c['M500'] < hi)
        print(f'  M {lo:.0e}-{hi:.0e}: n {s.sum()}, S/N 1h {np.sqrt(sn2(K1[s], 2048).sum()):.1f}')
    # what the 1-halo kappa(ell) looks like for a typical cluster
    i = np.argmin(np.abs(c['M500'] - 2e14) + 10 * np.abs(z - 0.25))
    print('typical cluster M500 %.2e z %.2f: kappa_ell at L=10,100,500,1000,2000:' % (c['M500'][i], z[i]),
          [f'{K1[i, l]:.2e}' for l in [10, 100, 500, 1000, 2000]], 'two-halo:', [f'{K2[i, l]:.2e}' for l in [10, 100, 500, 1000]])


# ------------------------------------------------------------------------------------------------------- map side
def filt_ell(lmin, lmax, ell, nside):
    import healpy as hp
    f = np.zeros(len(ell))
    sel = (ell >= lmin) & (ell <= lmax)
    f[sel] = 1.0
    # cosine taper over the last 10 % of the band and the first 4 multipoles
    tap = (ell > 0.9 * lmax) & sel
    f[tap] = 0.5 * (1 + np.cos(np.pi * (ell[tap] - 0.9 * lmax) / (0.1 * lmax)))
    low = (ell >= lmin) & (ell < lmin + 4)
    f[low] *= 0.5 * (1 - np.cos(np.pi * (ell[low] - lmin + 1) / 5))
    pw = hp.pixwin(nside, lmax=len(ell) - 1)
    return f, f * pw[:len(ell)]


def build_map(rec, nside, lmin, lmax):
    import healpy as hp
    alm = hp.read_alm(os.path.join(REC[rec], 'dat_klm.fits'))
    lm = hp.Alm.getlmax(len(alm))
    ell = np.arange(lm + 1)
    f, fpw = filt_ell(lmin, lmax, ell, nside)
    m = hp.alm2map(hp.almxfl(alm, f), nside, lmax=lm, verbose=False) if 'verbose' in hp.alm2map.__code__.co_varnames else hp.alm2map(hp.almxfl(alm, f), nside, lmax=lm)
    mask = hp.read_map(MASK[rec])
    mask = hp.ud_grade(mask, nside)
    return m, mask, ell, fpw


def selection(c, mask, nside):
    """clusters in the lensing footprint: unmasked within 30' (>= 0.99) and within 2 deg (>= 0.8)."""
    import healpy as hp
    l, b = radec_to_gal(c['RA'], c['DEC'])
    vec = hp.ang2vec(l, b, lonlat=True)
    keep = np.zeros(len(l), bool); f30 = np.zeros(len(l)); f120 = np.zeros(len(l))
    for i in range(len(l)):
        p = hp.query_disc(nside, vec[i], 0.5 * np.pi / 180)
        f30[i] = mask[p].mean()
        p = hp.query_disc(nside, vec[i], 2.0 * np.pi / 180)
        f120[i] = mask[p].mean()
    keep = (f30 >= 0.99) & (f120 >= 0.8)
    return keep, l, b, f30, f120


def annulus_profiles(m, mask, nside, l, b):
    """mean of the masked map in each annulus around each (l, b); NaN where an annulus has no unmasked pixel."""
    import healpy as hp
    vec = hp.ang2vec(l, b, lonlat=True)
    e = EDGES * ARCMIN
    prof = np.full((len(l), len(e) - 1), np.nan); npx = np.zeros((len(l), len(e) - 1))
    for i in range(len(l)):
        p = hp.query_disc(nside, vec[i], THETA_MAX * np.pi / 180)
        pv = np.array(hp.pix2vec(nside, p))
        cosd = np.clip(vec[i] @ pv, -1, 1)
        d = np.arccos(cosd)
        ok = mask[p] > 0.5
        idx = np.searchsorted(e, d[ok], side='right') - 1
        v = m[p][ok]
        good = (idx >= 0) & (idx < len(e) - 1)
        s = np.bincount(idx[good], weights=v[good], minlength=len(e) - 1)
        n = np.bincount(idx[good], minlength=len(e) - 1)
        with np.errstate(invalid='ignore'):
            prof[i] = s / n
        npx[i] = n
    return prof, npx


def profiles(args):
    import healpy as hp
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    c = load_cat()
    m, mask, ell, fpw = build_map(args.rec, args.nside, args.lmin, args.lmax)
    print(f'map built ({time.time() - t0:.0f}s): nside {args.nside}, band {args.lmin}-{args.lmax}, rms in mask {m[mask > 0.5].std():.3e}', flush=True)
    keep, l, b, f30, f120 = selection(c, mask, args.nside)
    print(f'selection: {keep.sum()} of {len(keep)} clusters in the lensing footprint ({time.time() - t0:.0f}s); '
          f'PSZ matches kept {c["psz"][keep].sum()} of {c["psz"].sum()}', flush=True)
    prof, npx = annulus_profiles(m, mask, args.nside, l[keep], b[keep])
    print(f'profiles done ({time.time() - t0:.0f}s)', flush=True)
    # mean-field stack as a diagnostic
    mf = hp.read_alm(os.path.join(REC[args.rec], 'mf_klm.fits'))
    f, _ = filt_ell(args.lmin, args.lmax, np.arange(hp.Alm.getlmax(len(mf)) + 1), args.nside)
    mfmap = hp.alm2map(hp.almxfl(mf, f), args.nside)
    mfprof, _ = annulus_profiles(mfmap, mask, args.nside, l[keep], b[keep])
    np.savez(os.path.join(OUT, f'{args.rec}{args.tag}_profiles.npz'), prof=prof, npx=npx, mfprof=mfprof, keep=keep, l=l, b=b,
             f30=f30, f120=f120, lmin=args.lmin, lmax=args.lmax, nside=args.nside, fpw=fpw, ell=ell,
             map_rms=m[mask > 0.5].std())
    print('saved', time.time() - t0)


def nulls(args):
    """rotated-catalogue stacks: the whole catalogue turned by a random Galactic longitude (and a random reflection in
    latitude), the same footprint selection, the same annuli; one mean profile per realisation (plus per-mass-bin means)."""
    import healpy as hp
    rng = np.random.default_rng(args.seed)
    c = load_cat()
    P = np.load(os.path.join(OUT, f'{args.rec}{args.tag}_profiles.npz'))
    m, mask, ell, fpw = build_map(args.rec, int(P['nside']), int(P['lmin']), int(P['lmax']))
    nside = int(P['nside'])
    l0, b0 = P['l'], P['b']
    z = c['BEST_Z']; M = c['M500']
    out = []; outw = []; nk = []
    t0 = time.time()
    fn = os.path.join(OUT, f'{args.rec}{args.tag}_nulls_seed{args.seed}.npz')
    for r in range(args.nrand):
        dl = rng.uniform(25, 335)
        refl = rng.integers(2) == 1
        l = (l0 + dl) % 360; b = -b0 if refl else b0
        keep, _, _, f30, f120 = selection({'RA': None, 'DEC': None}, mask, nside) if False else sel_lb(l, b, mask, nside)
        prof, npx = annulus_profiles(m, mask, nside, l[keep], b[keep])
        # store the per-cluster profiles compactly: the realisation's mean, and means by mass tercile and z bin
        rec = {'mean': np.nanmean(prof, 0), 'n': keep.sum()}
        Mk = M[keep]; zk = z[keep]
        for name, s in subsets(Mk, zk, c, keep).items():
            rec[name] = np.nanmean(prof[s], 0) if s.sum() > 0 else np.full(prof.shape[1], np.nan)
        out.append(rec); nk.append(keep.sum())
        if (r + 1) % 10 == 0 or r == 0:
            print(f'null {r + 1}/{args.nrand}: kept {keep.sum()}, {time.time() - t0:.0f}s', flush=True)
            np.savez(fn, recs=np.array(out, dtype=object), allow_pickle=True)
    np.savez(fn, recs=np.array(out, dtype=object))
    print('saved', fn)


def sel_lb(l, b, mask, nside):
    import healpy as hp
    vec = hp.ang2vec(l, b, lonlat=True)
    f30 = np.zeros(len(l)); f120 = np.zeros(len(l))
    for i in range(len(l)):
        f30[i] = mask[hp.query_disc(nside, vec[i], 0.5 * np.pi / 180)].mean()
        f120[i] = mask[hp.query_disc(nside, vec[i], 2.0 * np.pi / 180)].mean()
    keep = (f30 >= 0.99) & (f120 >= 0.8)
    return keep, l, b, f30, f120


def resid_fixed_mz(q, M, z):
    """log10 q minus its median in 6 x 6 cells of (log M500, z); NaN where a cell has <= 5 clusters."""
    lq = np.log10(np.maximum(q, 1e-9)); r = np.full_like(lq, np.nan)
    mb = np.digitize(np.log10(M), [13.8, 14.0, 14.2, 14.4, 14.6]); zb = np.digitize(z, [0.1, 0.2, 0.3, 0.45, 0.6])
    for i in np.unique(mb):
        for j in np.unique(zb):
            sij = (mb == i) & (zb == j)
            if sij.sum() > 5: r[sij] = lq[sij] - np.median(lq[sij])
    return r


SPLIT_PAIRS = [('M<4e14 DETres>0', 'M<4e14 DETres<=0'), ('cosmo M<4e14 DETres>0', 'cosmo M<4e14 DETres<=0'),
               ('cosmo M<4e14 compres>0', 'cosmo M<4e14 compres<=0'), ('cosmo M<4e14 extres>0', 'cosmo M<4e14 extres<=0'),
               ('PCONT<0.1 M<4e14 DETres>0', 'PCONT<0.1 M<4e14 DETres<=0'), ('EXT>=12 M<4e14 DETres>0', 'EXT>=12 M<4e14 DETres<=0')]


def subsets(Mk, zk, c, keep):
    """named subsets of the kept clusters (boolean over the kept array)."""
    s = {}
    s['M<1e14'] = Mk < 1e14; s['M1-2e14'] = (Mk >= 1e14) & (Mk < 2e14); s['M2-4e14'] = (Mk >= 2e14) & (Mk < 4e14)
    s['M>4e14'] = Mk >= 4e14
    s['z<0.1'] = zk < 0.1; s['z0.1-0.3'] = (zk >= 0.1) & (zk < 0.3); s['z0.3-0.6'] = (zk >= 0.3) & (zk < 0.6); s['z>0.6'] = zk >= 0.6
    s['z0.1-0.6'] = (zk >= 0.1) & (zk < 0.6)
    ext = c['EXT_LIKE'][keep]; pc = c['PCONT'][keep]; det = c['DET_LIKE_0'][keep]
    s['EXT_LIKE>=12'] = ext >= 12; s['EXT_LIKE<12'] = ext < 12
    s['PCONT<0.1'] = pc < 0.1; s['DET_LIKE>=40'] = det >= 40; s['DET_LIKE<40'] = det < 40
    s['nopsz'] = ~c['psz'][keep]; s['specz'] = np.array([t != 'photo_z' for t in c['ztype'][keep]])
    s['M>1e14 z0.1-0.6'] = (Mk >= 1e14) & (zk >= 0.1) & (zk < 0.6)
    s['M<1e14 z0.1-0.6'] = (Mk < 1e14) & (zk >= 0.1) & (zk < 0.6)
    # the referee's check (2026-11-01): purity. The cosmology sample, PCONT < 0.1 and EXT_LIKE >= 12 versions of the
    # fixed-(M, z) brightness splits, each with its own null stacks
    cosmo = c['cosmo'][keep]; exp_ = c['EXP'][keep]; low4 = Mk < 4e14
    s['cosmo'] = cosmo; s['cosmo M<4e14'] = cosmo & low4; s['PCONT<0.1 M<4e14'] = (pc < 0.1) & low4
    s['M<4e14'] = low4
    def halves(base, label, q):
        # residual medians computed inside the restricted sample, so the halves are matched in (M, z) within it
        r = np.full(len(Mk), np.nan); r[base] = resid_fixed_mz(q[base], Mk[base], zk[base])
        s[f'{label}>0'] = base & (r > 0); s[f'{label}<=0'] = base & (r <= 0)
    halves(low4, 'M<4e14 DETres', det)
    halves(cosmo & low4, 'cosmo M<4e14 DETres', det)
    halves(cosmo & low4, 'cosmo M<4e14 compres', det / np.maximum(exp_, 1e-3))
    halves(cosmo & low4, 'cosmo M<4e14 extres', ext / np.maximum(det, 1e-3))
    halves((pc < 0.1) & low4, 'PCONT<0.1 M<4e14 DETres', det)
    halves((ext >= 12) & low4, 'EXT>=12 M<4e14 DETres', det)
    return s


# ------------------------------------------------------------------------------------------------------- fitting
def stacked_model(hm, lnMg, zg, k1g, k2g, lm, z, fpw, ell, edges, bias_scale=1.0):
    """mean filtered annulus means of the model over clusters (lm = ln M500, z arrays)."""
    K1 = interp_grid(lnMg, zg, k1g, lm, z); K2 = interp_grid(lnMg, zg, k2g, lm, z)
    k1 = K1.mean(0); k2 = K2.mean(0)
    f = np.interp(ell, np.arange(len(fpw)), fpw)
    m1 = HaloModel.annulus_means_from_ell(k1, ell, f, edges)
    m2 = HaloModel.annulus_means_from_ell(k2, ell, f, edges) * bias_scale
    return m1, m2


def gls(d, C, X):
    Ci = np.linalg.inv(C)
    F = X.T @ Ci @ X; cov = np.linalg.inv(F)
    p = cov @ X.T @ Ci @ d
    chi2 = (d - X @ p) @ Ci @ (d - X @ p)
    return p, cov, chi2


def fit(args):
    c = load_cat()
    P = np.load(os.path.join(OUT, f'{args.rec}{args.tag}_profiles.npz'))
    prof, keep = P['prof'], P['keep']
    fpw, ell = P['fpw'], P['ell']
    z = c['BEST_Z'][keep]; M = c['M500'][keep]; Mmean = c['M500_mean'][keep]
    nb = prof.shape[1]
    # null stacks
    import glob
    recs = []
    for fn in sorted(glob.glob(os.path.join(OUT, f'{args.rec}{args.tag}_nulls_seed*.npz'))):
        recs += list(np.load(fn, allow_pickle=True)['recs'])
    print(f'{args.rec}: {keep.sum()} clusters, {len(recs)} null stacks, band {int(P["lmin"])}-{int(P["lmax"])}, nside {int(P["nside"])}')
    hm = HaloModel(lmax=3000)
    lnMg, zg, k1g, k2g, biasg = model_grid(hm)
    ellm = hm.ell
    mid = 0.5 * (EDGES[1:] + EDGES[:-1])

    def cov_of(name):
        A = np.array([r[name] for r in recs if name in r and np.all(np.isfinite(r[name]))])
        C = np.cov(A.T)
        hart = (len(A) - nb - 2) / (len(A) - 1)
        return C / hart, A

    def do_fit(sel, name, verbose=True):
        d = np.nanmean(prof[sel], 0)
        C, A = cov_of(name)
        sig = np.sqrt(np.diag(C))
        m1, m2 = stacked_model(hm, lnMg, zg, k1g, k2g, np.log(M[sel]), z[sel], fpw, ellm, EDGES)
        # (1) one amplitude on 1h+2h; (2) 1h amplitude with 2h fixed; (3) both free; (4) null chi2
        p1, c1, x1 = gls(d, C, (m1 + m2)[:, None])
        p2, c2, x2 = gls(d - m2, C, m1[:, None])
        p3, c3, x3 = gls(d, C, np.vstack([m1, m2]).T)
        x0 = d @ np.linalg.inv(C) @ d
        # null distribution of the single amplitude
        Ci = np.linalg.inv(C); X = (m1 + m2)[:, None]
        w = np.linalg.inv(X.T @ Ci @ X) @ X.T @ Ci
        Anull = (A @ w.T)[:, 0]
        # b_M grid: rescale masses, 2h bias follows the mass (recomputed from the grid at b_M M)
        bgrid = np.linspace(0.3, 2.0, 35); chi = []
        for bm in bgrid:
            mm1, mm2 = stacked_model(hm, lnMg, zg, k1g, k2g, np.log(bm * M[sel]), z[sel], fpw, ellm, EDGES)
            r = d - (mm1 + mm2); chi.append(r @ Ci @ r)
        chi = np.array(chi); ib = np.argmin(chi)
        # parabolic 68 % interval
        ok = chi - chi[ib] < 1
        blo, bhi = bgrid[ok].min(), bgrid[ok].max()
        # null-mean debiasing (a residual mean field or mask effect shows up as a nonzero mean over the rotated stacks)
        dn = A.mean(0)
        p1d, _, _ = gls(d - dn, C, (m1 + m2)[:, None])
        # delete-one-region jackknife over nside-2 Galactic HEALPix regions, and a bootstrap over clusters, for sigma_A
        import healpy as hp
        reg = hp.ang2pix(2, P['l'][keep][sel], P['b'][keep][sel], lonlat=True)
        ureg = np.unique(reg); Aj = []
        Xw = w
        for rr in ureg:
            dd = np.nanmean(prof[sel][reg != rr], 0); Aj.append((Xw @ dd)[0])
        Aj = np.array(Aj); nj = len(Aj)
        sA_jk = np.sqrt((nj - 1) / nj * np.sum((Aj - Aj.mean()) ** 2))
        rng = np.random.default_rng(7); ps = prof[sel]; Ab = []
        for _ in range(400):
            ii = rng.integers(0, len(ps), len(ps)); Ab.append((Xw @ np.nanmean(ps[ii], 0))[0])
        sA_bs = np.std(Ab)
        res = dict(name=name, n=int(sel.sum()), d=d, sig=sig, m1=m1, m2=m2, A=p1[0], sA=np.sqrt(c1[0, 0]), chi2_A=x1,
                   A_debiased=p1d[0], null_mean_profile=dn, sA_jk=sA_jk, njk=nj, sA_bs=sA_bs,
                   A1h=p2[0], sA1h=np.sqrt(c2[0, 0]), A1=p3[0], A2=p3[1], cov3=c3, chi2_3=x3, chi2_null=x0,
                   bM=bgrid[ib], bM_lo=blo, bM_hi=bhi, chi2_bM=chi[ib], Anull_mean=Anull.mean(), Anull_std=Anull.std(),
                   pnull=(np.abs(Anull - Anull.mean()) >= abs(p1[0] - Anull.mean())).mean(), nnull=len(A),
                   medM=np.median(M[sel]), medz=np.median(z[sel]))
        if verbose:
            print(f'\n[{name}] n={sel.sum()}  median M500 {res["medM"]:.2e}  median z {res["medz"]:.3f}  ({len(A)} nulls)')
            print('  theta(arcmin) | data x1e4 | sigma x1e4 | model 1h x1e4 | 2h x1e4')
            for i in range(nb):
                print(f'  {mid[i]:6.1f} | {d[i] * 1e4:8.3f} | {sig[i] * 1e4:6.3f} | {m1[i] * 1e4:7.3f} | {m2[i] * 1e4:7.3f}')
            print(f'  null chi2 {x0:.1f}/{nb}; A(1h+2h) = {p1[0]:.3f} +- {np.sqrt(c1[0, 0]):.3f} (chi2 {x1:.1f}/{nb - 1}); '
                  f'A_1h (2h fixed) = {p2[0]:.3f} +- {np.sqrt(c2[0, 0]):.3f}; free: A_1h {p3[0]:.3f} +- {np.sqrt(c3[0, 0]):.3f}, '
                  f'A_2h {p3[1]:.3f} +- {np.sqrt(c3[1, 1]):.3f} (corr {c3[0, 1] / np.sqrt(c3[0, 0] * c3[1, 1]):.2f}, chi2 {x3:.1f}/{nb - 2})')
            print(f'  null-mean profile x1e4: {np.round(dn * 1e4, 3)}; A after subtracting it: {p1d[0]:.3f}; '
                  f'sigma_A: GLS {np.sqrt(c1[0, 0]):.3f}, null spread {Anull.std():.3f}, jackknife ({nj} regions) {sA_jk:.3f}, bootstrap {sA_bs:.3f}')
            print(f'  b_M = {bgrid[ib]:.2f} [{blo:.2f}, {bhi:.2f}] (chi2 {chi[ib]:.1f}); null A: {Anull.mean():.3f} +- {Anull.std():.3f}, '
                  f'empirical two-sided p of the data A = {res["pnull"]:.4f}; detection {p1[0] / np.sqrt(c1[0, 0]):.1f} sigma (GLS), '
                  f'{(p1[0] - Anull.mean()) / Anull.std():.1f} sigma (null spread)')
        return res

    allsel = np.ones(len(z), bool)
    results = [do_fit(allsel, 'mean')]
    S = subsets(M, z, c, keep)
    for name in ['z0.1-0.6', 'M>1e14 z0.1-0.6', 'M<1e14 z0.1-0.6', 'M<1e14', 'M1-2e14', 'M2-4e14', 'M>4e14', 'z<0.1', 'z0.1-0.3',
                 'z0.3-0.6', 'z>0.6', 'EXT_LIKE>=12', 'EXT_LIKE<12', 'PCONT<0.1', 'DET_LIKE>=40', 'DET_LIKE<40', 'nopsz', 'specz',
                 'cosmo', 'cosmo M<4e14', 'PCONT<0.1 M<4e14', 'M<4e14'] + [n for pair in SPLIT_PAIRS for n in pair]:
        if name in recs[0] and name in S and S[name].sum() > 50:
            results.append(do_fit(S[name], name))
    # split differences with their own null: A_above - A_below per null realisation (both halves stacked in the same rotation)
    print('\nsplit differences (A_above - A_below), null from the paired rotated stacks:')
    for a, b in SPLIT_PAIRS:
        ra = [r for r in results if r['name'] == a]; rb = [r for r in results if r['name'] == b]
        if not ra or not rb or a not in recs[0]: continue
        ra, rb = ra[0], rb[0]
        Ca, Aa = cov_of(a); Cb, Ab = cov_of(b)
        Xa = (ra['m1'] + ra['m2'])[:, None]; Xb = (rb['m1'] + rb['m2'])[:, None]
        wa = np.linalg.inv(Xa.T @ np.linalg.inv(Ca) @ Xa) @ Xa.T @ np.linalg.inv(Ca)
        wb = np.linalg.inv(Xb.T @ np.linalg.inv(Cb) @ Xb) @ Xb.T @ np.linalg.inv(Cb)
        n = min(len(Aa), len(Ab)); dnull = (Aa[:n] @ wa.T)[:, 0] - (Ab[:n] @ wb.T)[:, 0]
        dA = ra['A'] - rb['A']
        print(f'  {a} minus {b}: dA = {dA:+.3f}; null {dnull.mean():+.3f} +- {dnull.std():.3f}; '
              f'{(dA - dnull.mean()) / dnull.std():+.1f} sigma; empirical two-sided p {(np.abs(dnull - dnull.mean()) >= abs(dA - dnull.mean())).mean():.4f} ({n} nulls)')
    # extra subsets (defined after the nulls were run): covariance = the full-sample null covariance scaled by n_mean / n_sub
    # (noise-dominated stacks), checked against the jackknife and bootstrap sigma_A printed with each fit
    Cm, Am = cov_of('mean'); nm = len(z)
    det = c['DET_LIKE_0'][keep]; ext = c['EXT_LIKE'][keep]; exp_ = c['EXP'][keep]; cr = c['CR500'][keep]
    crsn = cr / np.maximum(0.5 * (c['CR500_H'][keep] - c['CR500_L'][keep]), 1e-9)
    def resid(q):
        lq = np.log10(np.maximum(q, 1e-9)); r = np.zeros_like(lq)
        mb = np.digitize(np.log10(M), [13.8, 14.0, 14.2, 14.4, 14.6]); zb = np.digitize(z, [0.1, 0.2, 0.3, 0.45, 0.6])
        for i in np.unique(mb):
            for j in np.unique(zb):
                sij = (mb == i) & (zb == j)
                if sij.sum() > 5: r[sij] = lq[sij] - np.median(lq[sij])
                else: r[sij] = np.nan
        return r
    rdet, rext, rexp, rcr = resid(det), resid(ext), resid(exp_), resid(crsn)
    rdetx = resid(det / np.maximum(exp_, 1e-3))      # detection likelihood per unit exposure: count rate x compactness
    rextdet = resid(ext / np.maximum(det, 1e-3))     # extent per detection likelihood: extendedness
    low4 = M < 4e14
    extra = {'M<4e14': low4, 'M<4e14 z0.1-0.6': low4 & (z >= 0.1) & (z < 0.6),
             'M<4e14 DETres>0': low4 & (rdet > 0), 'M<4e14 DETres<=0': low4 & (rdet <= 0),
             'M<4e14 EXTres>0': low4 & (rext > 0), 'M<4e14 EXTres<=0': low4 & (rext <= 0),
             'M<4e14 EXPres>0': low4 & (rexp > 0), 'M<4e14 EXPres<=0': low4 & (rexp <= 0),
             'M<4e14 CRSNres>0': low4 & (rcr > 0), 'M<4e14 CRSNres<=0': low4 & (rcr <= 0),
             'M<4e14 DET/EXPres>0': low4 & (rdetx > 0), 'M<4e14 DET/EXPres<=0': low4 & (rdetx <= 0),
             'M<4e14 EXT/DETres>0': low4 & (rextdet > 0), 'M<4e14 EXT/DETres<=0': low4 & (rextdet <= 0),
             'M<4e14 DET>=40': low4 & (det >= 40), 'M<4e14 DET<40': low4 & (det < 40),
             'M2-4e14 DET>=40': (M >= 2e14) & low4 & (det >= 40), 'M2-4e14 DET<40': (M >= 2e14) & low4 & (det < 40),
             'M1-2e14 DET>=40': (M >= 1e14) & (M < 2e14) & (det >= 40), 'M1-2e14 DET<40': (M >= 1e14) & (M < 2e14) & (det < 40),
             'M>4e14 z>0.5': (M >= 4e14) & (z > 0.5), 'M>4e14 z<=0.5': (M >= 4e14) & (z <= 0.5)}
    print('\nextra subsets (scaled full-sample covariance):')
    print('  name | n | A +- s(scaled) | s_jk | s_bs | A_1h(2h fixed) | b_M [68%]')
    for name, sel in extra.items():
        if sel.sum() < 50: continue
        d = np.nanmean(prof[sel], 0); C = Cm * nm / sel.sum(); Ci = np.linalg.inv(C)
        m1, m2 = stacked_model(hm, lnMg, zg, k1g, k2g, np.log(M[sel]), z[sel], fpw, ellm, EDGES)
        p1, c1, x1 = gls(d, C, (m1 + m2)[:, None]); p2, c2, x2 = gls(d - m2, C, m1[:, None])
        X = (m1 + m2)[:, None]; w = np.linalg.inv(X.T @ Ci @ X) @ X.T @ Ci
        import healpy as hp
        reg = hp.ang2pix(2, P['l'][keep][sel], P['b'][keep][sel], lonlat=True); ureg = np.unique(reg)
        Aj = np.array([(w @ np.nanmean(prof[sel][reg != rr], 0))[0] for rr in ureg]); nj = len(Aj)
        sjk = np.sqrt((nj - 1) / nj * np.sum((Aj - Aj.mean()) ** 2))
        rng = np.random.default_rng(3); ps = prof[sel]
        sbs = np.std([(w @ np.nanmean(ps[rng.integers(0, len(ps), len(ps))], 0))[0] for _ in range(300)])
        bgrid = np.linspace(0.3, 2.0, 35); chi = []
        for bm in bgrid:
            mm1, mm2 = stacked_model(hm, lnMg, zg, k1g, k2g, np.log(bm * M[sel]), z[sel], fpw, ellm, EDGES)
            r = d - (mm1 + mm2); chi.append(r @ Ci @ r)
        chi = np.array(chi); ib = np.argmin(chi); ok = chi - chi[ib] < 1
        print(f'  {name:20s} | {sel.sum():5d} | {p1[0]:.3f} +- {np.sqrt(c1[0, 0]):.3f} | {sjk:.3f} | {sbs:.3f} | '
              f'{p2[0]:.3f} +- {np.sqrt(c2[0, 0]):.3f} | {bgrid[ib]:.2f} [{bgrid[ok].min():.2f}, {bgrid[ok].max():.2f}]')
        results.append(dict(name=name, n=int(sel.sum()), A=p1[0], sA=np.sqrt(c1[0, 0]), A1h=p2[0], sA1h=np.sqrt(c2[0, 0]),
                            bM=bgrid[ib], bM_lo=bgrid[ok].min(), bM_hi=bgrid[ok].max(), pnull=np.nan, sA_jk=sjk, sA_bs=sbs, d=d, m1=m1, m2=m2))
    # mean-field diagnostic
    mfp = np.nanmean(P['mfprof'], 0)
    print('\nmean-field stack x1e4 (should be ~0):', np.round(mfp * 1e4, 3))
    # posterior-mean masses instead of medians
    m1, m2 = stacked_model(hm, lnMg, zg, k1g, k2g, np.log(Mmean), z, fpw, ellm, EDGES)
    m1b, m2b = results[0]['m1'], results[0]['m2']
    print(f'model with posterior-mean masses vs medians: 1h ratio {np.sum(m1) / np.sum(m1b):.3f}, 2h ratio {np.sum(m2) / np.sum(m2b):.3f}')
    np.savez(os.path.join(OUT, f'{args.rec}{args.tag}_fit.npz'), results=np.array(results, dtype=object))
    print('\nsummary: name | n | A +- s | A_1h(2h fixed) | b_M [68%] | null p | GLS sigma')
    for r in results:
        print(f'  {r["name"]:18s} | {r["n"]:5d} | {r["A"]:.3f} +- {r["sA"]:.3f} | {r["A1h"]:.3f} +- {r["sA1h"]:.3f} | '
              f'{r["bM"]:.2f} [{r["bM_lo"]:.2f}, {r["bM_hi"]:.2f}] | {r["pnull"]:.4f} | {r["A"] / r["sA"]:.1f}')


def ystack(args):
    """the Planck MILCA Compton-y map in the same annuli at the same clusters (tSZ-leakage diagnostic for the MV
    reconstruction; the fallback measurement's first step), with rotated-catalogue nulls."""
    import healpy as hp
    from astropy.io import fits
    P = np.load(os.path.join(OUT, f'{args.rec}{args.tag}_profiles.npz'))
    keep = P['keep']; l0 = P['l']; b0 = P['b']
    y = hp.read_map(os.path.join(DATA, 'planck_y', 'milca_ymaps.fits'), field=0)
    nside = hp.get_nside(y)
    # the NILC release's own mask column (point sources and Galaxy), RING ordering like the maps
    mask = hp.read_map(os.path.join(DATA, 'planck_y', 'nilc_ymaps.fits'), field=3)
    mask = mask * (y > -1e30)
    ymap = np.where(mask > 0.5, y, 0.0)
    prof, npx = annulus_profiles(ymap, mask, nside, l0[keep], b0[keep])
    rng = np.random.default_rng(11); nulls = []
    for r in range(args.nrand):
        dl = rng.uniform(25, 335); refl = rng.integers(2) == 1
        l = (l0 + dl) % 360; b = -b0 if refl else b0
        k, _, _, _, _ = sel_lb(l, b, hp.ud_grade(mask, 512) if False else mask, nside)
        pr, _ = annulus_profiles(ymap, mask, nside, l[k], b[k])
        nulls.append(np.nanmean(pr, 0))
        print(f'y null {r + 1}/{args.nrand}', flush=True)
    np.savez(os.path.join(OUT, f'{args.rec}{args.tag}_yprofiles.npz'), prof=prof, npx=npx, nulls=np.array(nulls))
    d = np.nanmean(prof, 0); n = np.array(nulls)
    print('y stack x1e6 per annulus:', np.round(d * 1e6, 3))
    print('null spread x1e6       :', np.round(n.std(0) * 1e6, 3), ' null mean:', np.round(n.mean(0) * 1e6, 3))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['forecast', 'profiles', 'nulls', 'fit', 'ystack'])
    ap.add_argument('--rec', default='szdeproj', choices=list(REC))
    ap.add_argument('--nside', type=int, default=1024)
    ap.add_argument('--lmin', type=int, default=8)
    ap.add_argument('--lmax', type=int, default=2048)
    ap.add_argument('--nrand', type=int, default=300)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--tag', default='', help='suffix for the output files (e.g. _l1000 for a band variant)')
    args = ap.parse_args()
    {'forecast': forecast, 'profiles': profiles, 'nulls': nulls, 'fit': fit, 'ystack': ystack}[args.cmd](args)


if __name__ == '__main__':
    main()
