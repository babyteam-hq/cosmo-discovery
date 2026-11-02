"""The stacked thermal SZ signal of the eRASS1 clusters against their weak-lensing-calibrated masses (R-131, cosmo-deep-11).

Maps: the Planck 2015 MILCA and NILC Compton-y maps (PLA, nside 2048, Galactic, 10' beam; the NILC file's MASK column) and the
ACT DR6 NILC Compton-y map (NERSC dr6_nilc/ymaps_20230220: ACT+Planck ILC, CAR 0.5', equatorial, 1.6' beam, with the wide
GAL070 mask, the inpainted-regions and subtracted-source masks, the cluster mask and the CIB-deprojected variants). Catalogue:
eRASS1 primary v3.2 and cosmology v1.1 through R-121's `erass1_cmblens.load_cat` (M500 and R500 in the eRASS1 cosmology).

Per cluster the template is the Arnaud+2010 universal pressure profile at the catalogue R500, truncated at 5 R500, normalised to
the spherical Y_500 of the reference relation (Planck 2013 XX with 1-b = 0.8), beam-convolved, averaged over the lognormal
mass posterior; the data are masked annulus means around the X-ray positions; each bin's statistic is the GLS amplitude of the
bin-averaged template with a free constant, A = Y_500(observed) / Y_500(reference); covariance and null from rotated-map stacks.

usage (python -I; all data under data/, gitignored; logs under logs/):
  python -I erass1_ysz.py prep  --map milca|nilc|act|act_cib   # footprint, selection, stacking matrices, randoms
  python -I erass1_ysz.py stack --map ...  [--nrot 200]           # data, random and neighbour-model stacks, rotation nulls
  python -I erass1_ysz.py fit   --map ...                         # per-bin amplitudes, the Y-M fit, b_M
"""
import argparse, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import erass1_cmblens as r121  # noqa: E402

DATA = os.path.join(HERE, 'data')
OUT = os.path.join(DATA, 'erass1_ysz')
ARCMIN = np.pi / 180 / 60
MAPS = {
    'milca': dict(kind='planck', file=os.path.join(DATA, 'planck_y', 'milca_ymaps.fits'), nside=2048, fwhm=10.0,
                  edges=[0, 2.5, 5, 7.5, 10, 15, 20, 30, 40, 60.0], lmax=4096, frame='gal'),
    'nilc': dict(kind='planck', file=os.path.join(DATA, 'planck_y', 'nilc_ymaps.fits'), nside=2048, fwhm=10.0,
                 edges=[0, 2.5, 5, 7.5, 10, 15, 20, 30, 40, 60.0], lmax=4096, frame='gal'),
    'act': dict(kind='act', file=os.path.join(DATA, 'act_y', 'ilc_actplanck_ymap.fits'), nside=4096, fwhm=None,
                edges=[0, 1, 2, 3, 4, 6, 8, 10, 15, 20, 30.0], lmax=19999, frame='eq'),
    'act_cib': dict(kind='act', file=os.path.join(DATA, 'act_y', 'ilc_actplanck_ymap_deproj_cib_1.7_10.7.fits'), nside=4096,
                    fwhm=None, edges=[0, 1, 2, 3, 4, 6, 8, 10, 15, 20, 30.0], lmax=19999, frame='eq'),
}
NSIDE_COARSE = 512          # the 1-degree footprint check
MASS_EDGES = np.array([1, 1.4, 2, 2.8, 4, 8]) * 1e14
Z_EDGES = [0.1, 0.3, 0.8]
# reference relation (Planck 2013 XX, arXiv:1303.5080): E^-2/3 D_A^2 Y_500 / 1e-4 Mpc^2 = 10^-0.19 [(1-b) M / 6e14]^1.79
REF_LOGA, REF_B, REF_MPIV, REF_1MB = -0.19, 1.79, 6e14, 0.8
GH_X, GH_W = np.polynomial.hermite.hermgauss(5)


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def yref(M, z, co):
    """reference Y_500 in Mpc^2 (D_A^2 Y_500, the sr-to-Mpc^2 form) at an eRASS1 mass M (the relation at 0.8 M)."""
    return 1e-4 * 10 ** REF_LOGA * (REF_1MB * M / REF_MPIV) ** REF_B * co.efunc(z) ** (2 / 3)


def r500_mpc(M, z, co):
    rho_c = co.critical_density(z).to('Msun/Mpc3').value
    return (3 * M / (4 * np.pi * 500 * rho_c)) ** (1 / 3)


# --------------------------------------------------------------------------------------------------------- the template
class UPP:
    """Arnaud+2010 universal pressure profile shape; u(x) the projected (truncated at 5 R500) profile, S the spherical
    integral within R500, so that y(theta) = Y500[Mpc^2] / (R500[Mpc]^2 S) u(theta / theta500)."""
    P0, c500, gam, alp, bet = 8.403, 1.177, 0.3081, 1.0510, 5.4905

    def __init__(self):
        self.x = np.logspace(-3, np.log10(5.0), 600)
        t = np.linspace(0, 1, 3000) ** 2
        u = np.zeros_like(self.x)
        for i, xx in enumerate(self.x):
            L = np.sqrt(max(25 - xx ** 2, 0))
            l = L * t
            u[i] = 2 * np.trapezoid(self.p(np.sqrt(xx ** 2 + l ** 2)), l)
        self.u = u
        xs = np.logspace(-5, 0, 4000)
        self.S = np.trapezoid(self.p(xs) * 4 * np.pi * xs ** 2, xs)
        # cylindrical within 5 R500 over spherical within R500 (Arnaud+2010: 1.796)
        self.cyl5_over_sph = np.trapezoid(u * 2 * np.pi * self.x, self.x) / self.S

    @classmethod
    def p(cls, x):
        return cls.P0 / ((cls.c500 * x) ** cls.gam * (1 + (cls.c500 * x) ** cls.alp) ** ((cls.bet - cls.gam) / cls.alp))

    def profile(self, theta, theta500):
        """u(theta/theta500) on theta (rad), zero beyond 5 theta500."""
        xx = theta / theta500
        out = np.interp(np.log(np.clip(xx, self.x[0], 5.0)), np.log(self.x), self.u)
        out[xx >= 5.0] = 0.0
        return out


def beam_ell(mapname, ell):
    m = MAPS[mapname]
    if m['kind'] == 'planck':
        sig = m['fwhm'] * ARCMIN / np.sqrt(8 * np.log(2))
        return np.exp(-0.5 * ell * (ell + 1) * sig ** 2)
    b = np.loadtxt(os.path.join(DATA, 'act_y', 'ilc_beam.txt'))
    return np.interp(ell, b[:, 0], b[:, 1], right=0.0)


class TemplateGrid:
    """beam-convolved annulus means and real-space profiles of the unit-normalised template on a theta500 grid."""

    def __init__(self, mapname, pixwin=True):
        import healpy as hp
        from scipy.special import j0, j1
        m = MAPS[mapname]
        self.edges = np.array(m['edges']) * ARCMIN
        self.nb = len(self.edges) - 1
        self.upp = UPP()
        self.theta = np.logspace(np.log10(0.02 * ARCMIN), np.log10(4 * np.pi / 180), 1400)
        self.ell = np.arange(0, m['lmax'] + 1).astype(float)
        f = beam_ell(mapname, self.ell)
        if pixwin and m['kind'] == 'planck':
            pw = hp.pixwin(m['nside'], lmax=int(self.ell[-1]))
            f = f * np.concatenate([pw, np.zeros(len(self.ell) - len(pw))])[:len(self.ell)]
        self.filt = f
        self.t500 = np.logspace(np.log10(0.3 * ARCMIN), np.log10(90 * ARCMIN), 70)
        # annulus window functions W_k(ell)
        l = self.ell; e = self.edges
        W = np.empty((self.nb, len(l)))
        with np.errstate(invalid='ignore', divide='ignore'):
            for k in range(self.nb):
                a, b = e[k], e[k + 1]
                W[k] = 2 * (b * j1(l * b) - a * j1(l * a)) / (l * (b ** 2 - a ** 2))
        W[:, 0] = 1.0
        self.W = W
        self.ann = np.zeros((len(self.t500), self.nb))      # annulus means for y(theta) = u(theta/theta500)
        self.real = np.zeros((len(self.t500), len(self.theta)))  # the convolved real-space profile
        th = self.theta
        t0 = time.time()
        chunk = 1500
        for i, ts in enumerate(self.t500):
            y = self.upp.profile(th, ts)
            fl = np.empty(len(l))
            for a in range(0, len(l), chunk):
                J = j0(l[a:a + chunk, None] * th[None, :])
                fl[a:a + chunk] = 2 * np.pi * np.trapezoid(y[None, :] * J * th[None, :], th, axis=1)
            fl *= f
            self.ann[i] = np.trapezoid(l[None, :] * fl[None, :] * W, l, axis=1) / (2 * np.pi)
            # real space, on a coarser theta set (every 2nd point), interpolated back
            sub = th[::2]
            r = np.empty(len(sub))
            for a in range(0, len(l), chunk):
                J = j0(l[a:a + chunk, None] * sub[None, :])
                r_part = np.trapezoid(l[a:a + chunk, None] * fl[a:a + chunk, None] * J, l[a:a + chunk], axis=0) / (2 * np.pi)
                r = r + r_part if a else r_part
            self.real[i] = np.interp(th, sub, r)
        log(f'template grid {mapname}: {len(self.t500)} theta500 values in {time.time() - t0:.0f}s; '
            f'Y_cyl(5R500)/Y_sph(R500) = {self.upp.cyl5_over_sph:.3f}; S = {self.upp.S:.4f}')

    def annulus_means_of(self, y_theta):
        """beam-filtered annulus means of an arbitrary radial profile y(theta) given on self.theta (the two-halo shape)."""
        from scipy.special import j0
        l = self.ell; th = self.theta
        fl = np.empty(len(l)); chunk = 1500
        for a in range(0, len(l), chunk):
            J = j0(l[a:a + chunk, None] * th[None, :])
            fl[a:a + chunk] = 2 * np.pi * np.trapezoid(y_theta[None, :] * J * th[None, :], th, axis=1)
        fl *= self.filt
        return np.trapezoid(l[None, :] * fl[None, :] * self.W, l, axis=1) / (2 * np.pi)

    def annulus_means(self, theta500):
        """(n, nb) annulus means of u(theta/theta500) convolved, for an array of theta500 (rad)."""
        lt = np.log(np.clip(theta500, self.t500[0], self.t500[-1]))
        i = np.clip(np.searchsorted(self.t500, np.exp(lt)) - 1, 0, len(self.t500) - 2)
        fi = np.clip((lt - np.log(self.t500[i])) / (np.log(self.t500[i + 1]) - np.log(self.t500[i])), 0, 1)
        return (1 - fi)[:, None] * self.ann[i] + fi[:, None] * self.ann[i + 1]

    def real_profile(self, theta500, theta):
        lt = np.log(np.clip(theta500, self.t500[0], self.t500[-1]))
        i = int(np.clip(np.searchsorted(self.t500, np.exp(lt)) - 1, 0, len(self.t500) - 2))
        fi = np.clip((lt - np.log(self.t500[i])) / (np.log(self.t500[i + 1]) - np.log(self.t500[i])), 0, 1)
        prof = (1 - fi) * self.real[i] + fi * self.real[i + 1]
        return np.interp(theta, self.theta, prof, right=0.0)


def posterior_templates(tg, M, z, sig, co):
    """per-cluster template annulus means (n, nb) in y units, and the reference Y (Mpc^2, E^-2/3 removed), both averaged
    over the lognormal mass posterior; plus the median-mass versions."""
    n = len(M)
    nodes = np.log(M)[:, None] + np.sqrt(2) * sig[:, None] * GH_X[None, :]
    w = GH_W / np.sqrt(np.pi)
    Da = co.angular_diameter_distance(z).value
    T = np.zeros((n, tg.nb)); Yt = np.zeros(n)
    for j in range(len(w)):
        Mj = np.exp(nodes[:, j])
        R = r500_mpc(Mj, z, co); Y = yref(Mj, z, co)
        T += w[j] * (Y / (R ** 2 * tg.upp.S))[:, None] * tg.annulus_means(R / Da)
        Yt += w[j] * Y / co.efunc(z) ** (2 / 3)
    R = r500_mpc(M, z, co); Y = yref(M, z, co)
    Tmed = (Y / (R ** 2 * tg.upp.S))[:, None] * tg.annulus_means(R / Da)
    return T, Yt, Tmed, Y / co.efunc(z) ** (2 / 3)


# --------------------------------------------------------------------------------------------------------- map sources
class PlanckSource:
    def __init__(self, mapname):
        import healpy as hp
        m = MAPS[mapname]
        self.nside = m['nside']
        self.y = hp.read_map(m['file'], field=0).astype(np.float64)
        # the release's masks product (COM_CompMap_Compton-SZMap-masks_2048_R2.01, NESTED): M1-M4 the apodized Galactic masks
        # leaving 60, 50, 40, 30 % of the sky, M5 the point-source mask; clean = GAL60 > 0.99 and no point source
        mf = os.path.join(DATA, 'planck_y', 'masks_R2.01.fits')
        gal = hp.reorder(hp.read_map(mf, field=0, nest=True), n2r=True) > 0.99
        ps = hp.reorder(hp.read_map(mf, field=4, nest=True), n2r=True) > 0.5
        self.b = gal & ps & np.isfinite(self.y) & (self.y > -1e30)
        log(f'{mapname}: clean fraction of the sky {self.b.mean():.3f} (GAL60 x point sources)')
        self.y[~self.b] = 0.0
        self.frame = 'gal'

    def at(self, theta, phi):
        import healpy as hp
        p = hp.ang2pix(self.nside, theta, phi)
        return self.y[p], self.b[p]


class ACTSource:
    """the CAR map sampled at the nearest 0.5' pixel; clean = wide mask >= 0.99 and not inpainted and not a subtracted source."""

    def __init__(self, mapname):
        from astropy.io import fits
        m = MAPS[mapname]
        h = fits.open(m['file'], memmap=True)[0]
        self.hdr = h.header
        self.y = np.asarray(h.data, dtype=np.float32)
        self.b = act_clean_mask()
        bad = ~np.isfinite(self.y)
        self.y[bad] = 0.0; self.b[bad] = False
        self.frame = 'eq'
        self.ny, self.nx = self.y.shape
        self.crpix1, self.crpix2 = self.hdr['CRPIX1'], self.hdr['CRPIX2']
        self.cdelt1, self.cdelt2 = self.hdr['CDELT1'], self.hdr['CDELT2']
        self.crval1, self.crval2 = self.hdr['CRVAL1'], self.hdr['CRVAL2']

    def at(self, theta, phi):
        ra = np.degrees(phi); dec = 90 - np.degrees(theta)
        ra = (ra - self.crval1 + 180) % 360 - 180
        ix = np.rint((ra) / self.cdelt1 + self.crpix1 - 1).astype(np.int64)
        iy = np.rint((dec - self.crval2) / self.cdelt2 + self.crpix2 - 1).astype(np.int64)
        ok = (ix >= 0) & (ix < self.nx) & (iy >= 0) & (iy < self.ny)
        v = np.zeros(len(ra)); b = np.zeros(len(ra), bool)
        v[ok] = self.y[iy[ok], ix[ok]]; b[ok] = self.b[iy[ok], ix[ok]]
        return v, b



def act_bool(fn, thresh, above=True):
    """a CAR mask FITS (float64, 3.6 GB) to a bool array, read in row chunks through the memmap; cached as .npy."""
    from astropy.io import fits
    d = os.path.join(DATA, 'act_y')
    cache = os.path.join(d, fn.replace('.fits', '_bool.npy'))
    if os.path.exists(cache):
        return np.load(cache)
    h = fits.open(os.path.join(d, fn), memmap=True)[0]
    ny, nx = h.data.shape
    out = np.zeros((ny, nx), bool)
    for a in range(0, ny, 500):
        blk = np.asarray(h.data[a:a + 500])
        out[a:a + 500] = (blk >= thresh) if above else (blk > thresh)
    np.save(cache, out)
    return out


def act_clean_mask():
    """clean = wide mask >= 0.99 and not inpainted and not a subtracted source (cached)."""
    d = os.path.join(DATA, 'act_y')
    cache = os.path.join(d, 'clean_bool.npy')
    if os.path.exists(cache):
        return np.load(cache)
    b = act_bool('wide_mask.fits', 0.99)
    for fn in ['inpainted_mask.fits', 'subtracted_sources_mask.fits']:   # 1 = untouched, 0 = inpainted / subtracted
        b &= act_bool(fn, 0.5, above=False)
    np.save(cache, b)
    return b


def source(mapname):
    return PlanckSource(mapname) if MAPS[mapname]['kind'] == 'planck' else ACTSource(mapname)


def cluster_lonlat(c, frame):
    if frame == 'gal':
        return r121.radec_to_gal(c['RA'], c['DEC'])
    return c['RA'], c['DEC']


# ------------------------------------------------------------------------------------------------- stacking geometry
class Geometry:
    """the pixel sets: for each object the nside-work pixels of its annuli (rows (i, k)), its 15' disc, and its 1-degree disc at
    the coarse nside; all columns re-indexed to the touched pixels so a rotation only samples those."""

    def __init__(self, lon, lat, mapname, label):
        import healpy as hp
        from scipy.sparse import csr_matrix
        m = MAPS[mapname]
        self.nside = m['nside']; e = np.array(m['edges']) * ARCMIN; self.nb = len(e) - 1
        vec = hp.ang2vec(lon, lat, lonlat=True)
        n = len(lon); t0 = time.time()
        rows, cols, r15, c15, r60, c60 = [], [], [], [], [], []
        for i in range(n):
            p = hp.query_disc(self.nside, vec[i], e[-1])
            pv = np.array(hp.pix2vec(self.nside, p))
            d = np.arccos(np.clip(vec[i] @ pv, -1, 1))
            k = np.searchsorted(e, d, side='right') - 1
            ok = (k >= 0) & (k < self.nb)
            rows.append(i * self.nb + k[ok]); cols.append(p[ok])
            r15.append(np.full((d <= 15 * ARCMIN).sum(), i)); c15.append(p[d <= 15 * ARCMIN])
            q = hp.query_disc(NSIDE_COARSE, vec[i], np.pi / 180)
            r60.append(np.full(len(q), i)); c60.append(q)
            if (i + 1) % 2000 == 0:
                log(f'  {label}: geometry {i + 1}/{n} ({time.time() - t0:.0f}s)')
        cols = np.concatenate(cols); c15 = np.concatenate(c15); c60 = np.concatenate(c60)
        self.touched = np.unique(np.concatenate([cols, c15])); self.touched_c = np.unique(c60)
        ci = np.searchsorted(self.touched, cols); ci15 = np.searchsorted(self.touched, c15); ci60 = np.searchsorted(self.touched_c, c60)
        self.S = csr_matrix((np.ones(len(ci), np.float32), (np.concatenate(rows), ci)), shape=(n * self.nb, len(self.touched)))
        self.S15 = csr_matrix((np.ones(len(ci15), np.float32), (np.concatenate(r15), ci15)), shape=(n, len(self.touched)))
        self.S60 = csr_matrix((np.ones(len(ci60), np.float32), (np.concatenate(r60), ci60)), shape=(n, len(self.touched_c)))
        self.n15 = np.asarray(self.S15.sum(1)).ravel(); self.n60 = np.asarray(self.S60.sum(1)).ravel()
        self.n = n
        self.th, self.ph = hp.pix2ang(self.nside, self.touched)
        self.thc, self.phc = hp.pix2ang(NSIDE_COARSE, self.touched_c)

    def sample(self, src, dphi=0.0, reflect=False):
        th = np.pi - self.th if reflect else self.th; thc = np.pi - self.thc if reflect else self.thc
        v, b = src.at(th, (self.ph - dphi) % (2 * np.pi))
        _, bc = src.at(thc, (self.phc - dphi) % (2 * np.pi))
        return v, b, bc

    def profiles(self, v, b):
        bf = b.astype(np.float32)
        with np.errstate(invalid='ignore', divide='ignore'):
            return ((self.S @ (v * bf)) / (self.S @ bf)).reshape(self.n, self.nb)

    def selection(self, b, bc):
        f15 = (self.S15 @ b.astype(np.float32)) / self.n15
        f60 = (self.S60 @ bc.astype(np.float32)) / self.n60
        return (f15 >= 0.99) & (f60 >= 0.8)


def subsets(c, keep, insz=None, common=None):
    """boolean subsets over the kept clusters (the pre-registered bins and the informative rows)."""
    M = c['M500'][keep]; z = c['BEST_Z'][keep]; cosmo = c['cosmo'][keep]
    inz = (z >= Z_EDGES[0]) & (z < Z_EDGES[-1]); inM = (M >= MASS_EDGES[0]) & (M < MASS_EDGES[-1])
    s = {}
    for base, lab in [(cosmo, 'cosmo'), (np.ones(len(M), bool), 'all')]:
        s[f'{lab}'] = base & inz & inM
        for i in range(len(MASS_EDGES) - 1):
            mb = (M >= MASS_EDGES[i]) & (M < MASS_EDGES[i + 1])
            s[f'{lab} M{i}'] = base & inz & mb
            for j in range(len(Z_EDGES) - 1):
                s[f'{lab} M{i} z{j}'] = base & mb & (z >= Z_EDGES[j]) & (z < Z_EDGES[j + 1])
        for j in range(len(Z_EDGES) - 1):
            s[f'{lab} z{j}'] = base & inM & (z >= Z_EDGES[j]) & (z < Z_EDGES[j + 1])
    s['all z<0.1'] = (z < 0.1) & inM; s['all M<1e14'] = inz & (M < 1e14); s['all M>8e14'] = inz & (M >= 8e14)
    s['cosmo M<1e14'] = cosmo & inz & (M < 1e14)
    s['all PCONT<0.1'] = (c['PCONT'][keep] < 0.1) & inz & inM
    s['all nopsz'] = ~c['psz'][keep] & inz & inM; s['all psz'] = c['psz'][keep] & inz & inM
    if insz is not None:
        s['cosmo insz'] = s['cosmo'] & insz; s['cosmo nosz'] = s['cosmo'] & ~insz
        s['all insz'] = s['all'] & insz; s['all nosz'] = s['all'] & ~insz
    if common is not None:   # the clusters kept by both the Planck and the ACT footprints (defined after the stacks; no own rotations)
        for lab in ['cosmo', 'all']:
            s[f'{lab} common'] = s[lab] & common
            for i in range(len(MASS_EDGES) - 1):
                s[f'{lab} common M{i}'] = s[f'{lab} M{i}'] & common
    return s


# ------------------------------------------------------------------------------------------------------------- prep
def prep(args):
    import healpy as hp
    from scipy.sparse import save_npz
    os.makedirs(OUT, exist_ok=True)
    m = MAPS[args.map]; tag = args.map
    c = r121.load_cat()
    lon, lat = cluster_lonlat(c, m['frame'])
    src = source(args.map)
    log(f'{tag}: map read; catalogue {len(lon)} clusters (PCONT < 0.5, M500 > 0, z > 0.01)')
    G = Geometry(lon, lat, args.map, 'clusters')
    v, b, bc = G.sample(src)
    keep = G.selection(b, bc)
    M = c['M500']; z = c['BEST_Z']
    inrange = (z >= 0.1) & (z < 0.8) & (M >= 1e14) & (M < 8e14)
    log(f'{tag}: kept {keep.sum()} of {len(keep)} (clean within 15\', >= 80 % within 1 deg); in the mass-z range {int((keep & inrange).sum())} '
        f'of {inrange.sum()}; cosmology sample in range kept {int((keep & inrange & c["cosmo"]).sum())} of {int((inrange & c["cosmo"]).sum())}')
    for i in range(len(MASS_EDGES) - 1):
        mb = (M >= MASS_EDGES[i]) & (M < MASS_EDGES[i + 1]) & (z >= 0.1) & (z < 0.8)
        log(f'  M500 {MASS_EDGES[i]:.1e}-{MASS_EDGES[i + 1]:.1e}: all {int((mb & keep).sum())}/{mb.sum()}, cosmo {int((mb & keep & c["cosmo"]).sum())}/{int((mb & c["cosmo"]).sum())}')
    insz = None
    if m['kind'] == 'act':
        from astropy.io import fits
        cm = act_bool('cluster_mask.fits', 0.5, above=False)   # 1 = clean, 0 = an SZ cluster (confirmed or SNR > 6)
        th = np.radians(90 - lat); ph = np.radians(lon)
        ra = (np.degrees(ph) - src.crval1 + 180) % 360 - 180
        ix = np.rint(ra / src.cdelt1 + src.crpix1 - 1).astype(int); iy = np.rint((lat - src.crval2) / src.cdelt2 + src.crpix2 - 1).astype(int)
        ok = (ix >= 0) & (ix < src.nx) & (iy >= 0) & (iy < src.ny)
        insz = np.zeros(len(lon), bool); insz[ok] = ~cm[iy[ok], ix[ok]]
        del cm
        log(f'{tag}: clusters inside the release cluster mask (SZ SNR > 6): {int((insz & keep & inrange).sum())} of the kept in-range '
            f'({int((insz & keep & inrange & c["cosmo"]).sum())} in the cosmology sample)')
    # mask-matched randoms over the kept clusters' latitude range
    rng = np.random.default_rng(5)
    lo, hi = lat[keep].min() - 1, lat[keep].max() + 1
    got = []; tried = 0
    while len(got) < args.nrand:
        nn = 4 * (args.nrand - len(got)) + 100
        rl = rng.uniform(0, 360, nn); rb = np.degrees(np.arcsin(rng.uniform(np.sin(np.radians(lo)), np.sin(np.radians(hi)), nn)))
        tried += nn
        vec = hp.ang2vec(rl, rb, lonlat=True)
        _, b1 = src.at(np.radians(90 - rb), np.radians(rl))
        cand = np.where(b1)[0]
        for i in cand:
            p = hp.query_disc(m['nside'], vec[i], 15 * ARCMIN)
            _, bb = src.at(*hp.pix2ang(m['nside'], p))
            if bb.mean() < 0.99: continue
            q = hp.query_disc(NSIDE_COARSE, vec[i], np.pi / 180)
            _, bq = src.at(*hp.pix2ang(NSIDE_COARSE, q))
            if bq.mean() < 0.8: continue
            got.append((rl[i], rb[i]))
            if len(got) >= args.nrand: break
    got = np.array(got)
    log(f'{tag}: {len(got)} randoms of {tried} tried in lat {lo:.1f}..{hi:.1f}')
    GR = Geometry(got[:, 0], got[:, 1], args.map, 'randoms')
    vr, br, _ = GR.sample(src)
    rprof = GR.profiles(vr, br)
    rstack = np.nanmean(rprof, 0); rerr = np.nanstd(rprof, 0) / np.sqrt(len(got))
    log(f'{tag}: random stack x1e7 {np.round(rstack * 1e7, 2)} +- {np.round(rerr * 1e7, 2)}')
    import pickle
    with open(os.path.join(OUT, f'geom_{tag}.pkl'), 'wb') as f:
        pickle.dump(dict(G=G, keep=keep, lon=lon, lat=lat, insz=insz, rstack=rstack, rerr=rerr, rpos=got), f, protocol=4)
    log('prep done')


# ------------------------------------------------------------------------------------------------------------ stack
def neighbour_model(G, src, tg, c, lon, lat, co, mapname):
    """the model y map of every catalogue cluster (its convolved template at its reference Y_500, median mass), painted on the
    touched pixels only; returns (n_touched,) y values."""
    import healpy as hp
    m = MAPS[mapname]; nside = m['nside']
    idx = np.full(12 * nside ** 2, -1, np.int64); idx[G.touched] = np.arange(len(G.touched))
    model = np.zeros(len(G.touched))
    M = c['M500']; z = c['BEST_Z']
    Da = co.angular_diameter_distance(z).value; R = r500_mpc(M, z, co); Y = yref(M, z, co)
    t500 = R / Da; amp = Y / (R ** 2 * tg.upp.S)
    sig_b = (m['fwhm'] if m['fwhm'] else 1.6) * ARCMIN / np.sqrt(8 * np.log(2))
    vec = hp.ang2vec(lon, lat, lonlat=True)
    for i in range(len(M)):
        rad = 5 * t500[i] + 4 * sig_b
        p = hp.query_disc(nside, vec[i], rad)
        j = idx[p]; ok = j >= 0
        if not ok.any(): continue
        pv = np.array(hp.pix2vec(nside, p[ok]))
        d = np.arccos(np.clip(vec[i] @ pv, -1, 1))
        model[j[ok]] += amp[i] * tg.real_profile(t500[i], d)
    return model


def stack(args):
    import pickle
    m = MAPS[args.map]; tag = args.map
    c = r121.load_cat(); co = r121.cosmo()
    with open(os.path.join(OUT, f'geom_{tag}.pkl'), 'rb') as f:
        P = pickle.load(f)
    G, keep, lon, lat = P['G'], P['keep'], P['lon'], P['lat']
    src = source(args.map)
    tg = TemplateGrid(args.map)
    v, b, bc = G.sample(src)
    prof = G.profiles(v, b)
    log(f'{tag}: data stack x1e7 (all kept): {np.round(np.nanmean(prof[keep], 0) * 1e7, 2)}')
    # neighbour model: the catalogue's own clusters painted, stacked, minus each cluster's own template
    model = neighbour_model(G, src, tg, c, lon, lat, co, args.map)
    mprof = G.profiles(model, b)
    M = c['M500']; z = c['BEST_Z']
    Da = co.angular_diameter_distance(z).value; R = r500_mpc(M, z, co); Y = yref(M, z, co)
    own = (Y / (R ** 2 * tg.upp.S))[:, None] * tg.annulus_means(R / Da)
    neigh = mprof - own
    log(f'{tag}: model stack (own template, pixel-painted) x1e7 {np.round(np.nanmean(mprof[keep], 0) * 1e7, 2)}; analytic own '
        f'{np.round(np.nanmean(own[keep], 0) * 1e7, 2)}; neighbours {np.round(np.nanmean(neigh[keep], 0) * 1e7, 2)}')
    # rotation nulls: per-subset mean profiles and counts per rotation
    S = subsets(c, keep, None if P['insz'] is None else P['insz'][keep])
    names = sorted(S)
    rng = np.random.default_rng(args.seed)
    means = []; ns = []
    t0 = time.time()
    for r in range(args.nrot):
        dphi = np.radians(rng.uniform(25, 335)); refl = bool(rng.integers(2)) if m['frame'] == 'gal' else False
        vr, br, bcr = G.sample(src, dphi, refl)
        kr = G.selection(br, bcr)
        pr = G.profiles(vr, br)
        sel = {n: (np.zeros(len(keep), bool)) for n in names}
        for n in names:
            sel[n][keep] = S[n]
        means.append(np.array([np.nanmean(pr[sel[n] & kr], 0) if (sel[n] & kr).sum() else np.full(G.nb, np.nan) for n in names]))
        ns.append(np.array([(sel[n] & kr).sum() for n in names]))
        if (r + 1) % 10 == 0 or r == 0:
            log(f'  rotation {r + 1}/{args.nrot}: kept {kr.sum()} ({time.time() - t0:.0f}s)')
            np.savez(os.path.join(OUT, f'stack_{tag}.npz'), prof=prof, neigh=neigh, mprof=mprof, own=own, rot_means=np.array(means),
                     rot_ns=np.array(ns), names=np.array(names), rstack=P['rstack'], rerr=P['rerr'])
    np.savez(os.path.join(OUT, f'stack_{tag}.npz'), prof=prof, neigh=neigh, mprof=mprof, own=own, rot_means=np.array(means),
             rot_ns=np.array(ns), names=np.array(names), rstack=P['rstack'], rerr=P['rerr'])
    log('stack done')


# -------------------------------------------------------------------------------------------------------------- fit
def fit(args):
    import pickle
    m = MAPS[args.map]; tag = args.map
    c = r121.load_cat(); co = r121.cosmo()
    with open(os.path.join(OUT, f'geom_{tag}.pkl'), 'rb') as f:
        P = pickle.load(f)
    keep = P['keep']
    D = np.load(os.path.join(OUT, f'stack_{tag}.npz'), allow_pickle=True)
    prof, neigh, rot_means, rot_ns, names = D['prof'], D['neigh'], D['rot_means'], D['rot_ns'], list(D['names'])
    rstack = D['rstack']
    nb = prof.shape[1]; nrot = len(rot_means)
    tg = TemplateGrid(args.map)
    M = c['M500']; z = c['BEST_Z']
    sig = 0.5 * (np.log(c['M500_H']) - np.log(c['M500_L'])); sig = np.where(np.isfinite(sig) & (sig > 0), sig, np.nanmedian(sig))
    T, Yt, Tmed, Ytmed = posterior_templates(tg, M, z, sig, co)
    hm = r121.HaloModel(lmax=3000)   # the two-halo shape: the projected linear correlation at the subset's median (M, z), R-121's model
    h2cache = {}

    def twohalo(Mm, zm):
        key = (round(np.log10(Mm), 2), round(zm, 2))
        if key not in h2cache:
            _, k2, _ = hm.kappa_theta(Mm, zm)
            prof = np.interp(tg.theta, hm.theta, k2, right=0.0)
            h2 = tg.annulus_means_of(prof)
            h2cache[key] = h2 / h2[0]   # unit central annulus; the amplitude is free
        return h2cache[key]
    other = 'milca' if m['kind'] == 'act' else 'act'
    common = None
    if os.path.exists(os.path.join(OUT, f'geom_{other}.pkl')):
        with open(os.path.join(OUT, f'geom_{other}.pkl'), 'rb') as f:
            common = pickle.load(f)['keep'][keep]
        log(f'{tag}: common footprint with {other}: {int(common.sum())} of the kept clusters')
    S = subsets(c, keep, None if P['insz'] is None else P['insz'][keep], common)
    edges = np.array(m['edges'])
    # per-cluster aperture: the annulus area fractions inside 2 theta500 (catalogue mass), so sum_k w_ik prof_ik is the cluster's
    # cylindrical y integral within 2 R500 (per unit solid angle of the aperture), compared with the same functional of its template
    Da_all = co.angular_diameter_distance(z).value
    th2 = 2 * r500_mpc(M, z, co) / Da_all / ARCMIN
    frac = np.clip((th2[:, None] ** 2 - edges[None, :-1] ** 2) / (edges[None, 1:] ** 2 - edges[None, :-1] ** 2), 0, 1)
    wper = frac * (edges[1:] ** 2 - edges[:-1] ** 2)[None, :]   # area-weighted, in arcmin^2
    if m['kind'] == 'planck':
        ring = (edges[:-1] >= 20) & (edges[1:] <= 40)
    else:
        ring = (edges[:-1] >= 10) & (edges[1:] <= 20)
    wring = np.zeros(nb); wring[ring] = (edges[1:] ** 2 - edges[:-1] ** 2)[ring]; wring /= wring.sum()
    mid = 0.5 * (edges[1:] + edges[:-1])
    log(f'{tag}: {keep.sum()} kept clusters, {nrot} rotations, {nb} annuli; sigma_lnM median {np.median(sig[keep]):.3f}; '
        f'posterior/median template ratio (all kept, central annulus) {np.nanmean(T[keep], 0)[0] / np.nanmean(Tmed[keep], 0)[0]:.3f}')
    # aperture functional: y within the first K annuli (area weights) minus the annuli covering 20-40' (Planck) / 10-20' (ACT)
    area = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    if m['kind'] == 'planck':
        inner = edges[1:] <= 10.0; outer = (edges[:-1] >= 20) & (edges[1:] <= 40)
    else:
        inner = edges[1:] <= 5.0; outer = (edges[:-1] >= 10) & (edges[1:] <= 20)
    wap = np.zeros(nb); wap[inner] = area[inner] / area[inner].sum(); wap[outer] -= area[outer] / area[outer].sum()
    results = {}

    def one(name, verbose=True):
        s = np.zeros(len(keep), bool); s[keep] = S[name]
        n = int(s.sum())
        if n < 20: return None
        base = name if name in names else ('cosmo' if name.startswith('cosmo') else 'all')   # fallback: the parent's rotations, rescaled
        j = names.index(base)
        A = rot_means[:, j]; nr = rot_ns[:, j]; ok = np.all(np.isfinite(A), 1) & (nr > 10)
        A = A[ok]; nr = nr[ok]
        Ar = A.mean(0) + (A - A.mean(0)) * np.sqrt(nr / n)[:, None]   # rescaled to the sample's size
        C = np.cov(Ar.T); hart = (len(Ar) - nb - 2) / (len(Ar) - 1); C = C / hart
        Ci = np.linalg.inv(C); sg = np.sqrt(np.diag(C))
        d0 = np.nanmean(prof[s], 0); dd = d0 - rstack
        t = np.nanmean(T[s], 0); tmed = np.nanmean(Tmed[s], 0); ng = np.nanmean(neigh[s], 0)
        X = np.vstack([t, np.ones(nb)]).T
        p, cov, chi2 = r121.gls(dd, C, X)
        # neighbours subtracted with their amplitude tied to A: solve d = A (t + ng) + const
        pn, covn, chi2n = r121.gls(dd, C, np.vstack([t + ng, np.ones(nb)]).T)
        p1, cov1, chi21 = r121.gls(dd, C, t[:, None])                       # no constant
        h2 = twohalo(np.median(M[s]), np.median(z[s]))
        p2h, cov2h, chi2h = r121.gls(dd, C, np.vstack([t, h2, np.ones(nb)]).T)   # template + two-halo shape + constant
        p2hn, cov2hn, chi2hn = r121.gls(dd, C, np.vstack([t, h2]).T)              # template + two-halo shape
        pm, covm, _ = r121.gls(dd, C, np.vstack([tmed, np.ones(nb)]).T)       # median-mass template
        praw, _, _ = r121.gls(d0, C, X)
        # the inner annuli only (< 20' Planck, < 10' ACT)
        inn = edges[1:] <= (20.0 if m['kind'] == 'planck' else 10.0)
        pi_, covi, _ = r121.gls(dd[inn], C[np.ix_(inn, inn)], X[inn])
        # null distribution of A from the rotations (same weights)
        Wt = np.linalg.inv(X.T @ Ci @ X) @ X.T @ Ci
        Anull = (Ar @ Wt.T)[:, 0]
        # aperture statistic
        Aap = (wap @ dd) / (wap @ t); sap = np.sqrt(wap @ C @ wap) / (wap @ t)
        # per-cluster 2 theta500 aperture (mean over the clusters of the cluster's own aperture integral), raw and ring-subtracted
        Pi = prof[s] - rstack[None, :]; Wp = wper[s]
        ok_i = np.all(np.isfinite(Pi), 1)
        num = np.nanmean((Wp[ok_i] * Pi[ok_i]).sum(1)); den = np.nanmean((Wp[ok_i] * T[s][ok_i]).sum(1))
        numr = np.nanmean((Wp[ok_i] * (Pi[ok_i] - (Pi[ok_i] @ wring)[:, None])).sum(1))
        denr = np.nanmean((Wp[ok_i] * (T[s][ok_i] - (T[s][ok_i] @ wring)[:, None])).sum(1))
        wbar = np.nanmean(Wp[ok_i], 0)   # the mean weight vector, for the error through the rotation covariance
        A2r = num / den; sA2r = np.sqrt(wbar @ C @ wbar) / den
        wbr = wbar - wbar.sum() * wring
        A2s = numr / denr; sA2s = np.sqrt(wbr @ C @ wbr) / denr
        Ym = np.nanmean(Yt[s]); Ymmed = np.nanmean(Ytmed[s])
        res = dict(name=name, n=n, d=dd, draw=d0, sig=sg, t=t, ng=ng, A=p[0], sA=np.sqrt(cov[0, 0]), const=p[1], chi2=chi2,
                   A_neigh=pn[0], sA_neigh=np.sqrt(covn[0, 0]), A_noconst=p1[0], sA_noconst=np.sqrt(cov1[0, 0]), A_med=pm[0],
                   A_r2=A2r, sA_r2=sA2r, A_r2s=A2s, sA_r2s=sA2s, th2_med=np.median(th2[s]),
                   A_2h=p2h[0], sA_2h=np.sqrt(cov2h[0, 0]), A2h_amp=p2h[1] * 1e7, chi2_2h=chi2h, A_2hn=p2hn[0], sA_2hn=np.sqrt(cov2hn[0, 0]), chi2_2hn=chi2hn,
                   A_raw=praw[0], A_inner=pi_[0], sA_inner=np.sqrt(covi[0, 0]), Anull=Anull, A_ap=Aap, sA_ap=sap, Yref=Ym, Yref_med=Ymmed,
                   medM=np.median(M[s]), meanM=np.mean(M[s]), medz=np.median(z[s]), meankT=np.nanmean(c['KT'][s]), nrot=len(Ar),
                   chi2_null=dd @ Ci @ dd, nb=nb, neigh_frac=(wap @ ng) / (wap @ t), Wt=Wt)
        if verbose:
            print(f'\n[{name}] n={n}  median M500 {res["medM"]:.2e}  mean {res["meanM"]:.2e}  z {res["medz"]:.3f}  kT {res["meankT"]:.2f} keV  ({len(Ar)} rotations)')
            print('  theta  | data-rand x1e7 | sigma | template x1e7 | neighbours x1e7')
            for i in range(nb):
                print(f'  {mid[i]:6.1f} | {dd[i] * 1e7:8.3f} | {sg[i] * 1e7:6.3f} | {t[i] * 1e7:8.3f} | {ng[i] * 1e7:7.3f}')
            print(f'  A = {p[0]:.3f} +- {np.sqrt(cov[0, 0]):.3f} (const {p[1] * 1e7:+.3f}e-7; chi2 {chi2:.1f}/{nb - 2}; null chi2 {res["chi2_null"]:.1f}/{nb}); '
                  f'{p[0] / np.sqrt(cov[0, 0]):.1f} sigma; rotation null {Anull.mean():+.3f} +- {Anull.std():.3f}, '
                  f'empirical p {(Anull - Anull.mean() >= p[0] - Anull.mean()).mean():.4f}')
            print(f'  variants: neighbours tied {pn[0]:.3f} +- {np.sqrt(covn[0, 0]):.3f} (neighbour/own aperture {res["neigh_frac"]:.3f}); no constant {p1[0]:.3f} +- '
                  f'{np.sqrt(cov1[0, 0]):.3f}; median-mass template {pm[0]:.3f}; raw (no random subtraction) {praw[0]:.3f}; inner annuli only '
                  f'{pi_[0]:.3f} +- {np.sqrt(covi[0, 0]):.3f}; fixed aperture {Aap:.3f} +- {sap:.3f}')
            print(f'  two-halo shape: template + 2h + const A = {p2h[0]:.3f} +- {np.sqrt(cov2h[0, 0]):.3f} (2h central x1e7 {p2h[1] * 1e7:+.3f}, const {p2h[2] * 1e7:+.3f}e-7, chi2 {chi2h:.1f}/{nb - 3}); '
                  f'template + 2h A = {p2hn[0]:.3f} +- {np.sqrt(cov2hn[0, 0]):.3f} (chi2 {chi2hn:.1f}/{nb - 2})')
            print(f'  2 theta500 aperture (median {np.median(th2[s]):.1f}\'): raw A = {A2r:.3f} +- {sA2r:.3f}; ring-subtracted A = {A2s:.3f} +- {sA2s:.3f}')
            print(f'  Y_500 E^-2/3 D_A^2 at the bin: reference {Ym:.3e} Mpc^2 (median-mass {Ymmed:.3e}); measured {p[0] * Ym:.3e} +- {np.sqrt(cov[0, 0]) * Ym:.3e}')
        return res

    for name in sorted(S):
        r = one(name, verbose=True)
        if r is not None: results[name] = r
    # the Y-M fit over the five mass bins, per sample, with the joint covariance of the amplitudes from the rotations
    out = {}
    for lab, akey, skey in [('cosmo', 'A', 'sA'), ('all', 'A', 'sA'), ('cosmo', 'A_2h', 'sA_2h'), ('all', 'A_2h', 'sA_2h'),
                            ('cosmo', 'A_inner', 'sA_inner'), ('cosmo', 'A_ap', 'sA_ap'), ('cosmo', 'A_neigh', 'sA_neigh'),
                            ('cosmo', 'A_r2', 'sA_r2'), ('cosmo', 'A_r2s', 'sA_r2s'), ('all', 'A_r2', 'sA_r2'), ('all', 'A_r2s', 'sA_r2s'),
                            ('cosmo common', 'A', 'sA'), ('cosmo common', 'A_r2', 'sA_r2'), ('cosmo common', 'A_r2s', 'sA_r2s')]:
        bins = [f'{lab} M{i}' for i in range(len(MASS_EDGES) - 1)]
        if not all(b in results for b in bins): continue
        rs = [results[b] for b in bins]
        nmin = min(len(r['Anull']) for r in rs)
        An = np.array([r['Anull'][:nmin] for r in rs]).T   # (nrot, 5)
        Y = np.array([r[akey] * r['Yref'] for r in rs]); Yr = np.array([r['Yref'] for r in rs])
        CA = np.cov(An.T) * (nmin - 1) / (nmin - len(rs) - 2)
        # the variants' errors: the pre-registered statistic's null covariance rescaled to the variant's GLS error per bin
        if akey != 'A':
            sc = np.array([r[skey] / r['sA'] for r in rs]); CA = CA * np.outer(sc, sc)
        CY = CA * np.outer(Yr, Yr)
        # effective mass of each bin: the reference relation inverted at the bin's mean reference Y (E^-2/3 removed)
        Meff = np.array([REF_MPIV / REF_1MB * (r['Yref'] / (1e-4 * 10 ** REF_LOGA)) ** (1 / REF_B) for r in rs])
        Bgrid = np.arange(0.3, 3.5001, 0.01); chi = np.zeros(len(Bgrid)); amp = np.zeros(len(Bgrid))
        CYi = np.linalg.inv(CY)
        for i, B in enumerate(Bgrid):
            x = (Meff / 3e14) ** B
            a = (x @ CYi @ Y) / (x @ CYi @ x); amp[i] = a
            r_ = Y - a * x; chi[i] = r_ @ CYi @ r_
        ib = np.argmin(chi)
        def interval(dchi):
            ok = chi - chi[ib] < dchi
            return Bgrid[ok].min(), Bgrid[ok].max()
        b68 = interval(1.0); b95 = interval(3.84)
        x = (Meff / 3e14) ** Bgrid[ib]; sa = np.sqrt(1 / (x @ CYi @ x))
        # the reference relation's slope and amplitude at the pivot, for comparison
        yref_piv = 1e-4 * 10 ** REF_LOGA * (REF_1MB * 3e14 / REF_MPIV) ** REF_B
        # one amplitude with B fixed at the reference slope: the mass-scale ratio
        xr = (Meff / 3e14) ** REF_B; ar = (xr @ CYi @ Y) / (xr @ CYi @ xr); sar = np.sqrt(1 / (xr @ CYi @ xr))
        Aratio = ar / yref_piv; sAratio = sar / yref_piv
        bM = lambda Ar_, omb: REF_1MB * Ar_ ** (1 / REF_B) / omb
        print(f'\n== Y-M fit, {lab} sample, statistic {akey}, five mass bins (joint covariance from {nmin} rotations):')
        print('  bin | n | M_eff | A +- s | Y_500 E^-2/3 D_A^2 [1e-5 Mpc^2] +- s')
        for r, me in zip(rs, Meff):
            print(f'  {r["name"]:10s} | {r["n"]:5d} | {me:.2e} | {r[akey]:.3f} +- {r[skey]:.3f} | {r[akey] * r["Yref"] * 1e5:.3f} +- {r[skey] * r["Yref"] * 1e5:.3f}')
        print(f'  correlation of the bin amplitudes (null): off-diagonal max |r| {np.max(np.abs(np.corrcoef(An.T) - np.eye(len(rs)))):.2f}')
        print(f'  free fit: B = {Bgrid[ib]:.2f} 68 % [{b68[0]:.2f}, {b68[1]:.2f}] 95 % [{b95[0]:.2f}, {b95[1]:.2f}]; '
              f'log10 Y(3e14) = {np.log10(amp[ib]):.3f} +- {sa / amp[ib] / np.log(10):.3f} (Y = {amp[ib] * 1e5:.3f} +- {sa * 1e5:.3f} x1e-5 Mpc^2); '
              f'chi2 {chi[ib]:.1f}/{len(rs) - 2}; reference at the pivot {yref_piv * 1e5:.3f} (slope 1.79)')
        print(f'  slope fixed at 1.79: A_ratio = {Aratio:.3f} +- {sAratio:.3f} (chi2 {np.min([((Y - a_ * xr) @ CYi @ (Y - a_ * xr)) for a_ in [ar]]):.1f}/{len(rs) - 1}); '
              f'b_M at 1-b = 0.8: {bM(Aratio, 0.8):.3f} [{bM(Aratio - sAratio, 0.8):.3f}, {bM(Aratio + sAratio, 0.8):.3f}]; '
              f'at 1-b = 0.7: {bM(Aratio, 0.7):.3f}; 0.9: {bM(Aratio, 0.9):.3f}; 0.844 (Aymerich+2025): {bM(Aratio, 0.844):.3f}')
        out[f'{lab} {akey}'] = dict(bins=bins, Meff=Meff, Y=Y, CY=CY, B=Bgrid[ib], b68=b68, b95=b95, amp=amp[ib], samp=sa, Aratio=Aratio, sAratio=sAratio, chi=chi, Bgrid=Bgrid)
    np.savez(os.path.join(OUT, f'fit_{tag}.npz'), results=np.array([results], dtype=object), ym=np.array([out], dtype=object))
    print('\nsummary: sample | n | A +- s | sigma | null std | fixed-aperture A | A_inner | A_2h+const | A_neigh | chi2 | Y_meas x1e5 | 2th500 raw | 2th500 ring-sub')
    for name, r in results.items():
        print(f'  {name:16s} | {r["n"]:5d} | {r["A"]:.3f} +- {r["sA"]:.3f} | {r["A"] / r["sA"]:5.1f} | {r["Anull"].std():.3f} | '
              f'{r["A_ap"]:.3f} +- {r["sA_ap"]:.3f} | {r["A_inner"]:.3f} +- {r["sA_inner"]:.3f} | {r["A_2h"]:.3f} +- {r["sA_2h"]:.3f} | '
              f'{r["A_neigh"]:.3f} | {r["chi2"]:.1f}/{r["chi2_2h"]:.1f} | {r["A"] * r["Yref"] * 1e5:.3f} | {r["A_r2"]:.3f} +- {r["sA_r2"]:.3f} | {r["A_r2s"]:.3f} +- {r["sA_r2s"]:.3f}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['prep', 'stack', 'fit'])
    ap.add_argument('--map', default='milca', choices=list(MAPS))
    ap.add_argument('--nrand', type=int, default=10000)
    ap.add_argument('--nrot', type=int, default=200)
    ap.add_argument('--seed', type=int, default=1)
    args = ap.parse_args()
    {'prep': prep, 'stack': stack, 'fit': fit}[args.cmd](args)


if __name__ == '__main__':
    main()
