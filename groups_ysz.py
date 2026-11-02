"""The one-halo thermal SZ signal of group catalogues at 10^12.5-10^14.5 M_sun on the ACT DR6 y map (R-132, cosmo-deep-12).

Pre-registration: `results-R132.md`. Maps as R-131 (`erass1_ysz.py`): the ACT DR6 NILC ACT+Planck ILC Compton-y map (CAR, 0.5',
1.6' beam, with the wide GAL070, inpainted and subtracted-source masks), its CIB-deprojected variant, and the Planck 2015 MILCA
y map (10', reprojected to a 2' CAR grid). Catalogues: the DESI DR9 Yang et al. 2021 groups (DESI DR1 VAC gfinder v1.0), the
Tully 2015 2MRS groups, the Saulder et al. 2016 SDSS DR12 groups. Templates: R-131's Arnaud+2010 UPP at the group's R500 from
its catalogue halo mass (converted to M500c), normalised to the Planck 2013 XX reference Y500, beam-convolved, averaged over a
lognormal mass posterior. Statistic: annulus means of the masked map around each group, from FFT convolutions of the map with the
annulus kernels in declination strips (so RA-rotation nulls and randoms are shifted lookups); per bin the GLS amplitude A of the
bin template with the neighbour (two-halo) model and a constant; Y500 = A <Y_ref>.

usage (python -I; data under data/, gitignored; logs under logs/):
  python -I groups_ysz.py cats                        # catalogues -> data/groups_ysz/cat_<name>.npz (M500c, radio flag)
  python -I groups_ysz.py sample --map act|act_cib|milca   # annulus means at every group, rotations, randoms, footprint
  python -I groups_ysz.py paint  --map act|milca --cat desi|tully|saulder   # neighbour model map and its samples
  python -I groups_ysz.py fit    --map act|act_cib|milca --cat desi|tully|saulder
"""
import argparse, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import erass1_ysz as r131  # noqa: E402  (UPP, yref, r500_mpc, GH nodes, log)

DATA = os.path.join(HERE, 'data')
OUT = os.path.join(DATA, 'groups_ysz')
ARCMIN = np.pi / 180 / 60
log = r131.log
H_YANG = 0.674          # Yang et al. 2021 use Planck18
EDGES_ACT = np.array([0, 1, 2, 3, 4, 6, 8, 10, 15, 20, 30, 45, 60.0])
EDGES_PLANCK = np.array([0, 2.5, 5, 7.5, 10, 15, 20, 30, 40, 60.0])
MAPS = {
    'act': dict(kind='act', file=os.path.join(DATA, 'act_y', 'ilc_actplanck_ymap.fits'), pix=0.5, edges=EDGES_ACT, lmax=20000),
    'act_cib': dict(kind='act', file=os.path.join(DATA, 'act_y', 'ilc_actplanck_ymap_deproj_cib_1.7_10.7.fits'), pix=0.5,
                    edges=EDGES_ACT, lmax=20000),
    'act_cibdb': dict(kind='act', file=os.path.join(DATA, 'act_y', 'ilc_actplanck_ymap_deproj_cib_cibdBeta_1.7_10.7.fits'), pix=0.5,
                      edges=EDGES_ACT, lmax=20000),
    'milca': dict(kind='planck', file=os.path.join(DATA, 'planck_y', 'COM_CompMap_Compton-SZMap-milca-ymaps_2048_R2.00.fits'),
                  pix=2.0, edges=EDGES_PLANCK, lmax=4096, fwhm=10.0),
}
DEC_LO, DEC_HI = -63.0, 23.0        # the ACT y map's declination range (CRPIX2 7561 at dec 0, 0.5' pixels, 10320 rows)
NROT = 23                           # RA rotations of 15 n degrees
NSUB_ROT = 150000                   # groups per catalogue carrying rotation samples
SIG_M = 0.25                        # lognormal mass scatter (dex), pre-registered
RADIO_FLUX_MJY, RADIO_RADIUS_ARCMIN = 5.0, 3.0

CATS = {
    'desi': dict(zlo=0.1, zhi=0.6, edges=[12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0, 14.5], headline_from=13.0,
                 zsplit=[0.1, 0.3, 0.6], richcheck=3, rmax=30.0, minrich=1),
    # Tully: 1,500 < <Vcmba> < 10,000 km/s; Nmb >= 2 (pre-registered); Saulder: Ntot >= 2
    'tully': dict(zlo=1500 / 299792.458, zhi=10000 / 299792.458, edges=[12.5, 13.0, 13.5, 14.0, 15.5], headline_from=12.5, zsplit=None,
                  richcheck=3, rmax=60.0, minrich=2),
    'saulder': dict(zlo=0.02, zhi=0.11, edges=[12.5, 13.0, 13.5, 14.0, 15.5], headline_from=12.5, zsplit=None, richcheck=3, rmax=60.0,
                    minrich=2),
}


def cosmo():
    from astropy.cosmology import FlatLambdaCDM
    return FlatLambdaCDM(H0=67.4, Om0=0.315)


# ------------------------------------------------------------------------------------------------------- mass conversion
def _nfw_m(x):
    return np.log1p(x) - x / (1 + x)


def c200c(M200c, z):
    """Dutton & Maccio 2014 NFW c200c(M, z), M in M_sun (h = 0.674)."""
    a = 0.520 + (0.905 - 0.520) * np.exp(-0.617 * z ** 1.21)
    b = -0.101 + 0.026 * z
    return 10 ** (a + b * np.log10(M200c * H_YANG / 1e12))


def convert_mass(M, z, delta_in, ref_in, co):
    """M at overdensity delta_in wrt ref_in ('mean' or 'crit') -> (M200c, M500c), NFW with the DM14 concentration; computed
    on a (log M, z) grid and interpolated (the direct version took 25 min on 23M groups)."""
    from scipy.interpolate import RegularGridInterpolator
    z = np.asarray(z, float); M = np.asarray(M, float)
    lg = np.linspace(10.5, 16.0, 111); zg = np.linspace(0.0, 1.6, 33)
    LG, ZG = np.meshgrid(lg, zg, indexing='ij')
    m2, m5 = _convert_mass_direct(10 ** LG.ravel(), ZG.ravel(), delta_in, ref_in, co)
    pts = np.c_[np.clip(np.log10(M), lg[0], lg[-1]), np.clip(z, zg[0], zg[-1])]
    f2 = RegularGridInterpolator((lg, zg), np.log10(m2).reshape(LG.shape)); f5 = RegularGridInterpolator((lg, zg), np.log10(m5).reshape(LG.shape))
    return 10 ** f2(pts), 10 ** f5(pts)


def _convert_mass_direct(M, z, delta_in, ref_in, co):
    z = np.asarray(z, float); M = np.asarray(M, float)
    rho_c = co.critical_density(z).to('Msun/Mpc3').value
    rho_in = rho_c * (co.Om(z) if ref_in == 'mean' else 1.0) * delta_in
    M200 = M.copy()
    for _ in range(12):
        c = c200c(M200, z)
        r200 = (3 * M200 / (4 * np.pi * 200 * rho_c)) ** (1 / 3)
        rs = r200 / c
        # r_in: mean density inside = rho_in  ->  M200 m(r/rs)/m(c) = (4/3) pi r^3 rho_in; solve for x = r/rs by bisection
        lo = np.full_like(M200, 0.05); hi = np.full_like(M200, 60.0)
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            f = M200 * _nfw_m(mid) / _nfw_m(c) - 4 / 3 * np.pi * (mid * rs) ** 3 * rho_in
            lo = np.where(f > 0, mid, lo); hi = np.where(f > 0, hi, mid)
        Min = M200 * _nfw_m(0.5 * (lo + hi)) / _nfw_m(c)
        M200 = M200 * (M / Min) ** 0.9
    c = c200c(M200, z); r200 = (3 * M200 / (4 * np.pi * 200 * rho_c)) ** (1 / 3); rs = r200 / c
    lo = np.full_like(M200, 0.05); hi = np.full_like(M200, 60.0)
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        f = M200 * _nfw_m(mid) / _nfw_m(c) - 4 / 3 * np.pi * (mid * rs) ** 3 * rho_c * 500
        lo = np.where(f > 0, mid, lo); hi = np.where(f > 0, hi, mid)
    M500 = M200 * _nfw_m(0.5 * (lo + hi)) / _nfw_m(c)
    return M200, M500


# ------------------------------------------------------------------------------------------------------------ catalogues
def read_tsv(fn):
    rows = [l.rstrip('\n').split('\t') for l in open(fn) if not l.startswith('#') and l.strip()]
    hdr = rows[0]; data = rows[3:]
    out = {}
    for j, k in enumerate(hdr):
        col = [r[j].strip() if j < len(r) else '' for r in data]
        try:
            out[k] = np.array([float(v) if v else np.nan for v in col])
        except ValueError:
            out[k] = np.array(col)
    return out


def radio_sources():
    """NVSS and SUMSS (HEASARC tdat) above the pre-registered flux: ra, dec (deg)."""
    import gzip
    cache = os.path.join(OUT, 'radio.npz')
    if os.path.exists(cache):
        d = np.load(cache); return d['ra'], d['dec']
    ras, decs = [], []
    for fn, icol, scale in [('heasarc_nvss.tdat.gz', 'flux_20_cm', 1.0), ('heasarc_sumss.tdat.gz', 'int_flux_36_cm', (1400 / 843) ** -0.7)]:
        cols = None; indata = False; n = 0
        with gzip.open(os.path.join(DATA, 'radio', fn), 'rt', errors='replace') as f:
            for line in f:
                if line.startswith('line[1]'):
                    cols = line.split('=', 1)[1].split(); continue
                if line.strip() == '<DATA>':
                    indata = True; continue
                if line.strip() == '<END>':
                    break
                if not indata: continue
                p = line.rstrip('\n').split('|')
                try:
                    ra = float(p[cols.index('ra')]); dec = float(p[cols.index('dec')])
                    s = float(p[cols.index(icol)]) * scale
                except (ValueError, IndexError):
                    continue
                if s >= RADIO_FLUX_MJY:
                    ras.append(ra); decs.append(dec); n += 1
        log(f'radio: {fn} columns {cols[:12]}...; {n} sources above {RADIO_FLUX_MJY} mJy (1.4 GHz)')
    ra = np.array(ras); dec = np.array(decs)
    os.makedirs(OUT, exist_ok=True); np.savez(cache, ra=ra, dec=dec)
    return ra, dec


def radio_flag(ra, dec):
    from scipy.spatial import cKDTree
    rra, rdec = radio_sources()
    def xyz(a, d):
        a = np.radians(a); d = np.radians(d)
        return np.c_[np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)]
    t = cKDTree(xyz(rra, rdec))
    r = 2 * np.sin(0.5 * RADIO_RADIUS_ARCMIN * ARCMIN)
    cnt = t.query_ball_point(xyz(ra, dec), r, return_length=True)
    return cnt > 0


def load_desi(co):
    from astropy.io import fits
    h = fits.open(os.path.join(DATA, 'groups', 'DESIDR9.y1.v1_group.fits'), memmap=True)[1]
    log(f'desi: {h.header["NAXIS2"]} groups in the VAC; columns {h.columns.names}')
    logm = np.asarray(h.data['GRP_LOGM'], dtype=np.float64) - np.log10(H_YANG)     # M_h in M_sun
    dec = np.asarray(h.data['GRP_DEC'], dtype=np.float64)
    keep = (logm >= 12.5) & (dec > DEC_LO - 1.5) & (dec < DEC_HI + 1.5)
    idx = np.where(keep)[0]
    log(f'desi: {len(idx)} with log M_h >= 12.5 in dec {DEC_LO - 1.5}..{DEC_HI + 1.5}')
    ra = np.asarray(h.data['GRP_RA'], dtype=np.float64)[idx]; dec = dec[idx]
    z = np.asarray(h.data['GRP_Z'], dtype=np.float64)[idx]; rich = np.asarray(h.data['RICH'], dtype=np.int64)[idx]
    logm = logm[idx]
    ok = (z > 0.005) & (z < 1.5) & np.isfinite(logm)
    ra, dec, z, rich, logm = ra[ok], dec[ok], z[ok], rich[ok], logm[ok]
    M200, M500 = convert_mass(10 ** logm, z, 180.0, 'mean', co)
    return dict(ra=ra, dec=dec, z=z, logMh=logm, logM500=np.log10(M500), logM200=np.log10(M200), rich=rich.astype(float))


def load_tully(co):
    t = read_tsv(os.path.join(DATA, 'groups', 'tully2015_table3.tsv'))
    V = t['<Vcmba>']; D = 10 ** (t['DM'] / 5 - 5)          # Mpc, the catalogue's distance modulus
    ok = np.isfinite(t['Mlum']) & (t['Mlum'] > 0) & np.isfinite(V) & (V > 0)
    z = V[ok] / 299792.458
    M200, M500 = convert_mass(t['Mlum'][ok], z, 200.0, 'crit', co)
    ra = t['_RA.icrs'] if '_RA.icrs' in t else t['_RAJ2000']; de = t['_DE.icrs'] if '_DE.icrs' in t else t['_DEJ2000']
    return dict(ra=ra[ok], dec=de[ok], z=z, logMh=np.log10(t['Mlum'][ok]), logM500=np.log10(M500),
                logM200=np.log10(M200), rich=t['Nmb'][ok], dist=D[ok], V=V[ok])


def load_saulder(co):
    t = read_tsv(os.path.join(DATA, 'groups', 'saulder2016_sdss_groups.tsv'))
    ok = np.isfinite(t['logMtot']) & np.isfinite(t['z']) & (t['z'] > 0)
    z = t['z'][ok]
    M200, M500 = convert_mass(10 ** t['logMtot'][ok], z, 200.0, 'crit', co)
    return dict(ra=t['RAJ2000'][ok], dec=t['DEJ2000'][ok], z=z, logMh=t['logMtot'][ok], logM500=np.log10(M500),
                logM200=np.log10(M200), rich=t['Ntot'][ok])


def cats(args):
    os.makedirs(OUT, exist_ok=True)
    co = cosmo()
    for name, loader in [('tully', load_tully), ('saulder', load_saulder), ('desi', load_desi)]:
        fn = os.path.join(OUT, f'cat_{name}.npz')
        if os.path.exists(fn) and not args.force:
            continue
        t0 = time.time(); c = loader(co)
        c['radio'] = radio_flag(c['ra'], c['dec'])
        np.savez(fn, **c)
        inrange = (c['dec'] > DEC_LO) & (c['dec'] < DEC_HI)
        log(f'{name}: {len(c["ra"])} groups, {inrange.sum()} in the ACT dec range, radio fraction {c["radio"][inrange].mean():.3f} '
            f'({time.time() - t0:.0f}s)')
        e = CATS[name]['edges']
        for i in range(len(e) - 1):
            m = inrange & (c['logMh'] >= e[i]) & (c['logMh'] < e[i + 1]) & (c['z'] >= CATS[name]['zlo']) & (c['z'] < CATS[name]['zhi'])
            log(f'  log M_h {e[i]}-{e[i + 1]}: {m.sum()} (radio {c["radio"][m].mean() if m.any() else 0:.3f}; '
                f'<log M500c> {c["logM500"][m].mean() if m.any() else 0:.2f}; <z> {c["z"][m].mean() if m.any() else 0:.3f})')


def load_cat(name):
    d = np.load(os.path.join(OUT, f'cat_{name}.npz'))
    return {k: d[k] for k in d.files}


# -------------------------------------------------------------------------------------------------------- the CAR grids
class Grid:
    """a CAR grid at pix arcmin: row 0 at DEC_LO, column ix = round(nx/2 - ra_wrapped/pix); RA periodic."""

    def __init__(self, pix):
        self.pix = pix
        self.nx = int(round(360 * 60 / pix)); self.ny = int(round((DEC_HI - DEC_LO) * 60 / pix))
        self.dec0 = DEC_LO                    # the centre of row 0 (the ACT map's row 0 is dec -63.0 exactly)

    def rowcol(self, ra, dec):
        raw = (np.asarray(ra) + 180) % 360 - 180
        ix = np.rint(self.nx / 2 - raw * 60 / self.pix).astype(np.int64) % self.nx
        iy = np.rint((np.asarray(dec) - self.dec0) * 60 / self.pix).astype(np.int64)
        return iy, ix

    def dec_of_row(self, iy):
        return self.dec0 + iy * self.pix / 60

    def ra_of_col(self, ix):
        return ((self.nx / 2 - ix) * self.pix / 60) % 360


def act_grid_offset():
    """the ACT CAR map's row of dec DEC_LO: CRPIX2 = 7561 at dec 0 -> row index 7560 is dec 0; row of DEC_LO = 7560 + DEC_LO*120."""
    return int(round(7560 + DEC_LO * 120))


def read_act(mapname):
    """the ACT y map on our Grid(0.5) (its own grid, shifted) and the clean mask."""
    from astropy.io import fits
    m = MAPS[mapname]
    h = fits.open(m['file'], memmap=True)[0]
    g = Grid(0.5)
    r0 = act_grid_offset()
    assert h.header['CRPIX1'] == 21601.0 and h.header['CRPIX2'] == 7561.0 and abs(h.header['CDELT1'] + 1 / 120) < 1e-9
    assert r0 >= 0 and r0 + g.ny <= h.header['NAXIS2'], (r0, g.ny)
    y = np.asarray(h.data[r0:r0 + g.ny], dtype=np.float32)
    # map column ix_map = CRPIX1-1 - ra/CDELT... : ra = -(ix_map - 21600)/120 -> ix_map = 21600 - ra*120 = our ix (nx/2 = 21600)
    b = np.load(os.path.join(DATA, 'act_y', 'wide_mask_bool.npy'), mmap_mode='r')[r0:r0 + g.ny]
    b = np.array(b)
    for fn in ['inpainted_mask_bool.npy', 'subtracted_sources_mask_bool.npy']:   # R-131's act_bool: True = untouched (> 0.5)
        b &= np.array(np.load(os.path.join(DATA, 'act_y', fn), mmap_mode='r')[r0:r0 + g.ny])
    bad = ~np.isfinite(y); y[bad] = 0; b[bad] = False
    y[~b] = 0
    log(f'{mapname}: {y.shape}, clean fraction {b.mean():.3f}')
    return g, y, b


def read_planck():
    """the MILCA y map and its GAL60 x point-source mask, reprojected (nearest pixel) onto Grid(2.0)."""
    import healpy as hp
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    m = MAPS['milca']
    cache = os.path.join(OUT, 'milca_car2.npz')
    g = Grid(2.0)
    if os.path.exists(cache):
        d = np.load(cache); return g, d['y'], d['b']
    y = hp.read_map(m['file'], field=0).astype(np.float64)
    mf = os.path.join(DATA, 'planck_y', 'COM_CompMap_Compton-SZMap-masks_2048_R2.01.fits')
    gal = hp.reorder(hp.read_map(mf, field=0, nest=True), n2r=True) > 0.99
    ps = hp.reorder(hp.read_map(mf, field=4, nest=True), n2r=True) > 0.5
    bh = gal & ps & np.isfinite(y)
    y[~bh] = 0
    iy, ix = np.meshgrid(np.arange(g.ny), np.arange(g.nx), indexing='ij')
    ra = g.ra_of_col(ix).ravel(); dec = g.dec_of_row(iy).ravel()
    c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs').galactic
    p = hp.ang2pix(2048, c.l.deg, c.b.deg, lonlat=True)
    yc = y[p].reshape(g.ny, g.nx).astype(np.float32); bc = bh[p].reshape(g.ny, g.nx)
    np.savez(cache, y=yc, b=bc)
    log(f'milca: reprojected to {yc.shape}, clean fraction {bc.mean():.3f}')
    return g, yc, bc


# ------------------------------------------------------------------------------------------------- annulus means by FFT
def strip_kernels(g, dec_c, edges_arcmin, pad, padx):
    """indicator kernels (nk, 2pad+1, 2padx+1) of the annuli for a stamp centred at declination dec_c, exact spherical
    distances on the CAR pixel grid; plus the 15' disc (footprint)."""
    dy = (np.arange(-pad, pad + 1) * g.pix / 60)[:, None]
    dx = (np.arange(-padx, padx + 1) * g.pix / 60)[None, :]
    d1 = np.radians(dec_c); d2 = np.radians(dec_c + dy); da = np.radians(dx)
    cosd = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(da)
    dist = np.degrees(np.arccos(np.clip(cosd, -1, 1))) * 60
    e = np.asarray(edges_arcmin)
    K = np.zeros((len(e), dist.shape[0], dist.shape[1]), np.float32)
    for k in range(len(e) - 1):
        K[k] = (dist >= e[k]) & (dist < e[k + 1])
    K[-1] = dist < 15.0
    return K


def annulus_maps(y, b, g, edges_arcmin, rmax_arcmin, rows, pad):
    """annulus means (nk-1, nrows, nx) and the 15' clean fraction (nrows, nx) for the map rows [rows[0], rows[1])."""
    import scipy.fft as sf
    r0, r1 = rows
    dec_c = g.dec_of_row(0.5 * (r0 + r1 - 1))
    padx = int(np.ceil(pad / np.cos(np.radians(max(abs(g.dec_of_row(r0)), abs(g.dec_of_row(r1 - 1))))) * 1.02)) + 2
    padx = min(padx, g.nx // 2 - 1)
    K = strip_kernels(g, dec_c, edges_arcmin, pad, padx)
    nr = r1 - r0 + 2 * pad
    ys = np.zeros((nr, g.nx), np.float32); bs = np.zeros((nr, g.nx), np.float32)
    a0, a1 = max(r0 - pad, 0), min(r1 + pad, g.ny)
    ys[a0 - (r0 - pad):a1 - (r0 - pad)] = y[a0:a1] * b[a0:a1]
    bs[a0 - (r0 - pad):a1 - (r0 - pad)] = b[a0:a1]
    Fy = sf.rfft2(ys, workers=4); Fb = sf.rfft2(bs, workers=4)
    out = np.zeros((K.shape[0] - 1, r1 - r0, g.nx), np.float32); frac = None
    kern = np.zeros((nr, g.nx), np.float32)
    for k in range(K.shape[0]):
        kern[:] = 0
        kk = K[k]
        # centred at (0, 0) with wrap-around
        kern[:pad + 1, :padx + 1] = kk[pad:, padx:]
        kern[:pad + 1, g.nx - padx:] = kk[pad:, :padx]
        kern[nr - pad:, :padx + 1] = kk[:pad, padx:]
        kern[nr - pad:, g.nx - padx:] = kk[:pad, :padx]
        Fk = sf.rfft2(kern, workers=4)
        num = sf.irfft2(Fy * Fk, s=(nr, g.nx), workers=4)[pad:pad + (r1 - r0)]
        den = sf.irfft2(Fb * Fk, s=(nr, g.nx), workers=4)[pad:pad + (r1 - r0)]
        if k < K.shape[0] - 1:
            with np.errstate(invalid='ignore', divide='ignore'):
                out[k] = np.where(den > 0.5, num / np.maximum(den, 0.5), np.nan)
        else:
            frac = den / kk.sum()
    return out, frac


def rotations(g):
    return [int(round(15 * n * 60 / g.pix)) for n in range(1, NROT + 1)]


def sample(args):
    """annulus means at every catalogue group (all three catalogues) on one map; rotation samples for a subsample; randoms."""
    os.makedirs(OUT, exist_ok=True)
    m = MAPS[args.map]
    if m['kind'] == 'act':
        g, y, b = read_act(args.map)
    else:
        g, y, b = read_planck()
    edges = m['edges']; nb = len(edges) - 1
    pad = int(np.ceil(edges[-1] / g.pix))
    strip = int(round(2 * 60 / g.pix))     # 2-degree strips
    cats_ = {}
    rng = np.random.default_rng(132)
    for name in ['desi', 'tully', 'saulder']:
        c = load_cat(name)
        iy, ix = g.rowcol(c['ra'], c['dec'])
        n = len(iy)
        sub = np.zeros(n, bool)
        cand = np.where((iy >= 0) & (iy < g.ny))[0]
        sub[rng.choice(cand, size=min(NSUB_ROT, len(cand)), replace=False)] = True
        cats_[name] = dict(iy=iy, ix=ix, prof=np.full((n, nb), np.nan, np.float32), f15=np.zeros(n, np.float32), sub=sub,
                           rot=np.full((sub.sum(), NROT, nb), np.nan, np.float32), subidx=np.cumsum(sub) - 1)
    # randoms: uniform on the sphere in the dec range, later cut on footprint
    nrand = 400000
    rra = rng.uniform(0, 360, nrand)
    rdec = np.degrees(np.arcsin(rng.uniform(np.sin(np.radians(DEC_LO + 0.3)), np.sin(np.radians(DEC_HI - 0.3)), nrand)))
    riy, rix = g.rowcol(rra, rdec)
    rprof = np.full((nrand, nb), np.nan, np.float32); rf15 = np.zeros(nrand, np.float32)
    shifts = rotations(g)
    t0 = time.time()
    for r0 in range(0, g.ny, strip):
        r1 = min(r0 + strip, g.ny)
        A, F = annulus_maps(y, b, g, edges, edges[-1], (r0, r1), pad)
        for name, d in cats_.items():
            s = np.where((d['iy'] >= r0) & (d['iy'] < r1))[0]
            if len(s) == 0: continue
            d['prof'][s] = A[:, d['iy'][s] - r0, d['ix'][s]].T
            d['f15'][s] = F[d['iy'][s] - r0, d['ix'][s]]
            ss = s[d['sub'][s]]
            for j, sh in enumerate(shifts):
                d['rot'][d['subidx'][ss], j] = A[:, d['iy'][ss] - r0, (d['ix'][ss] + sh) % g.nx].T
        s = np.where((riy >= r0) & (riy < r1))[0]
        rprof[s] = A[:, riy[s] - r0, rix[s]].T; rf15[s] = F[riy[s] - r0, rix[s]]
        log(f'{args.map}: rows {r0}-{r1} of {g.ny} ({time.time() - t0:.0f}s)')
    for name, d in cats_.items():
        np.savez(os.path.join(OUT, f'samp_{args.map}_{name}.npz'), prof=d['prof'], f15=d['f15'], sub=d['sub'], rot=d['rot'])
    np.savez(os.path.join(OUT, f'samp_{args.map}_randoms.npz'), prof=rprof, f15=rf15, ra=rra, dec=rdec)
    log('sample done')


def sample_randoms(args):
    """mask-matched randoms only (post-hoc, R-136: the 400k of `sample` leave the random profile's error comparable to the
    low-mass bins' errors): NRAND uniform positions in the dec range, annulus means and the 15' clean fraction."""
    m = MAPS[args.map]
    g, y, b = read_act(args.map) if m['kind'] == 'act' else read_planck()
    edges = m['edges']; nb = len(edges) - 1
    pad = int(np.ceil(edges[-1] / g.pix)); strip = int(round(2 * 60 / g.pix))
    rng = np.random.default_rng(136)
    nrand = args.nrand
    rra = rng.uniform(0, 360, nrand)
    rdec = np.degrees(np.arcsin(rng.uniform(np.sin(np.radians(DEC_LO + 0.3)), np.sin(np.radians(DEC_HI - 0.3)), nrand)))
    riy, rix = g.rowcol(rra, rdec)
    rprof = np.full((nrand, nb), np.nan, np.float32); rf15 = np.zeros(nrand, np.float32)
    t0 = time.time()
    for r0 in range(0, g.ny, strip):
        r1 = min(r0 + strip, g.ny)
        A, F = annulus_maps(y, b, g, edges, edges[-1], (r0, r1), pad)
        s = np.where((riy >= r0) & (riy < r1))[0]
        rprof[s] = A[:, riy[s] - r0, rix[s]].T; rf15[s] = F[riy[s] - r0, rix[s]]
        log(f'{args.map} randoms: rows {r0}-{r1} of {g.ny} ({time.time() - t0:.0f}s)')
    np.savez(os.path.join(OUT, f'rand_{args.map}.npz'), prof=rprof, f15=rf15, ra=rra, dec=rdec)
    log('randoms done')


# ------------------------------------------------------------------------------------------------------------ templates
class Templates:
    """R-131's TemplateGrid for arbitrary annulus edges and beam: annulus means and the real-space profile of the
    unit-normalised, beam-convolved UPP on a theta500 grid."""

    def __init__(self, mapname, edges_arcmin, xoff=0.0):
        """xoff (R-136, post-hoc): a 2D Gaussian centre offset of rms xoff * R500 per axis, applied as exp(-l^2 (xoff theta500)^2 / 2)."""
        from scipy.special import j0, j1
        m = MAPS[mapname]
        mk = 'act' if m['kind'] == 'act' else mapname          # every ACT variant has the same beam
        cache = os.path.join(OUT, f'templates_{mk}_{len(edges_arcmin)}' + (f'_x{xoff:.2f}' if xoff else '') + '.npz')
        self.edges = np.asarray(edges_arcmin) * ARCMIN; self.nb = len(self.edges) - 1
        self.upp = r131.UPP()
        self.t500 = np.logspace(np.log10(0.3 * ARCMIN), np.log10(120 * ARCMIN), 80)
        self.theta = np.logspace(np.log10(0.02 * ARCMIN), np.log10(6 * np.pi / 180), 1500)
        if os.path.exists(cache):
            d = np.load(cache); self.ann = d['ann']; self.real = d['real']; return
        ell = np.arange(0, m['lmax'] + 1).astype(float)
        if m['kind'] == 'act':
            bm = np.loadtxt(os.path.join(DATA, 'act_y', 'ilc_beam.txt'))
            f = np.interp(ell, bm[:, 0], bm[:, 1], right=0.0)
        else:
            import healpy as hp
            sig = m['fwhm'] * ARCMIN / np.sqrt(8 * np.log(2))
            f = np.exp(-0.5 * ell * (ell + 1) * sig ** 2) * hp.pixwin(2048, lmax=int(ell[-1]))[:len(ell)]
        l = ell; e = self.edges
        W = np.empty((self.nb, len(l)))
        with np.errstate(invalid='ignore', divide='ignore'):
            for k in range(self.nb):
                a, bb = e[k], e[k + 1]
                W[k] = 2 * (bb * j1(l * bb) - a * j1(l * a)) / (l * (bb ** 2 - a ** 2))
        W[:, 0] = 1.0
        self.ann = np.zeros((len(self.t500), self.nb)); self.real = np.zeros((len(self.t500), len(self.theta)))
        th = self.theta; t0 = time.time(); chunk = 1000
        for i, ts in enumerate(self.t500):
            yy = self.upp.profile(th, ts)
            fl = np.empty(len(l))
            for a in range(0, len(l), chunk):
                J = j0(l[a:a + chunk, None] * th[None, :])
                fl[a:a + chunk] = 2 * np.pi * np.trapezoid(yy[None, :] * J * th[None, :], th, axis=1)
            fl *= f * (np.exp(-0.5 * (l * xoff * ts) ** 2) if xoff else 1.0)
            self.ann[i] = np.trapezoid(l[None, :] * fl[None, :] * W, l, axis=1) / (2 * np.pi)
            sub = th[::3]; r = np.zeros(len(sub))
            for a in range(0, len(l), chunk):
                J = j0(l[a:a + chunk, None] * sub[None, :])
                r += np.trapezoid(l[a:a + chunk, None] * fl[a:a + chunk, None] * J, l[a:a + chunk], axis=0) / (2 * np.pi)
            self.real[i] = np.interp(th, sub, r)
        np.savez(cache, ann=self.ann, real=self.real)
        log(f'templates {mapname}: {len(self.t500)} theta500 values in {time.time() - t0:.0f}s')

    def _w(self, theta500):
        lt = np.log(np.clip(theta500, self.t500[0], self.t500[-1]))
        i = np.clip(np.searchsorted(self.t500, np.exp(lt)) - 1, 0, len(self.t500) - 2)
        fi = np.clip((lt - np.log(self.t500[i])) / (np.log(self.t500[i + 1]) - np.log(self.t500[i])), 0, 1)
        return i, fi

    def annulus_means(self, theta500):
        i, fi = self._w(np.atleast_1d(theta500))
        return (1 - fi)[:, None] * self.ann[i] + fi[:, None] * self.ann[i + 1]

    def real_profile(self, theta500, theta):
        i, fi = self._w(np.atleast_1d(theta500)); i = int(i[0]); fi = float(fi[0])
        prof = (1 - fi) * self.real[i] + fi * self.real[i + 1]
        return np.interp(theta, self.theta, prof, right=0.0)


def group_templates(tg, M500, z, co, sig=SIG_M):
    """per-group annulus means (n, nb) in y and the reference <Y_ref E^-2/3> (Mpc^2), averaged over the lognormal mass
    posterior (Gauss-Hermite, 5 nodes); the median-mass template too."""
    n = len(M500)
    nodes = np.log(M500)[:, None] + np.sqrt(2) * sig * np.log(10) * r131.GH_X[None, :]
    w = r131.GH_W / np.sqrt(np.pi)
    Da = co.angular_diameter_distance(z).value
    T = np.zeros((n, tg.nb)); Yt = np.zeros(n)
    for j in range(len(w)):
        Mj = np.exp(nodes[:, j])
        R = r131.r500_mpc(Mj, z, co); Y = r131.yref(Mj, z, co)
        T += w[j] * (Y / (R ** 2 * tg.upp.S))[:, None] * tg.annulus_means(R / Da)
        Yt += w[j] * Y / co.efunc(z) ** (2 / 3)
    return T, Yt


# -------------------------------------------------------------------------------------------------- the neighbour model
def paint(args):
    """the y map of every catalogue group (log M_h >= 12.5) at its reference template, at 1' (ACT) or 2' (Planck) pixels,
    then sampled at the catalogue's positions with the map's annuli; the group's own painted template is sampled the same
    way on an empty stamp, so T2h = sampled - own is the catalogued-neighbour term."""
    co = cosmo()
    m = MAPS[args.map]
    pix = 1.0 if m['kind'] == 'act' else 2.0
    g = Grid(pix)
    tg = Templates(args.map, m['edges'])
    c = load_cat(args.cat)
    iy, ix = g.rowcol(c['ra'], c['dec'])
    inmap = (iy >= -200) & (iy < g.ny + 200)
    M500 = 10 ** c['logM500']; z = c['z']
    if args.cat == 'desi':
        # 23M groups above 12.5: paint those above 13.0 at z < 0.9 (the 12.5-13.0 half-dex adds ~0.4 of the 13.0-13.5 half's
        # y, by n(M) M^1.8; B absorbs the amplitude, the shape is the same clustering); the deviation is in results-R132.md
        inmap &= (c['logMh'] >= 13.0) & (z < 0.9)
    Da = co.angular_diameter_distance(z).value
    R = r131.r500_mpc(M500, z, co); Y = r131.yref(M500, z, co)
    t500 = R / Da; amp = Y / (R ** 2 * tg.upp.S)
    fwhm = 1.6 if m['kind'] == 'act' else 10.0
    ymap = np.zeros((g.ny, g.nx), np.float32)
    idx = np.where(inmap)[0]
    order = idx[np.argsort(t500[idx])]
    t0 = time.time()
    # paint in groups of similar theta500 and declination row to vectorise: per object, a stamp
    for cnt, i in enumerate(order):
        rmax = 5 * t500[i] / ARCMIN + 1.5 * fwhm       # arcmin
        ry = int(np.ceil(rmax / pix)); cosd = max(np.cos(np.radians(c['dec'][i])), 0.2)
        rx = int(np.ceil(rmax / pix / cosd))
        ys = np.arange(iy[i] - ry, iy[i] + ry + 1); ok = (ys >= 0) & (ys < g.ny); ys = ys[ok]
        if len(ys) == 0: continue
        xs = (np.arange(ix[i] - rx, ix[i] + rx + 1)) % g.nx
        dec_c = np.radians(c['dec'][i]); d2 = np.radians(g.dec_of_row(ys))[:, None]
        da = np.radians((xs - ix[i] + rx) * 0 + (np.arange(-rx, rx + 1) * pix / 60))[None, :]
        # offsets from the exact centre: the group's position, not the pixel centre
        dra = np.radians(((g.ra_of_col(xs) - c['ra'][i] + 180) % 360 - 180))[None, :]
        cosd_ = np.sin(dec_c) * np.sin(d2) + np.cos(dec_c) * np.cos(d2) * np.cos(dra)
        dist = np.arccos(np.clip(cosd_, -1, 1))
        prof = amp[i] * tg.real_profile(t500[i], dist.ravel()).reshape(dist.shape)
        ymap[np.ix_(ys, xs)] += prof.astype(np.float32)
        if (cnt + 1) % 200000 == 0:
            log(f'  paint {args.cat}/{args.map}: {cnt + 1}/{len(order)} ({time.time() - t0:.0f}s)')
    np.save(os.path.join(OUT, f'paint_{args.map}_{args.cat}.npy'), ymap)
    log(f'paint {args.cat}/{args.map}: {len(order)} groups painted ({time.time() - t0:.0f}s); map mean y {ymap.mean():.3e}')
    # sample at the catalogue positions with the map's annuli (mask = ones), strips
    edges = m['edges']; nb = len(edges) - 1
    pad = int(np.ceil(edges[-1] / g.pix)); strip = int(round(2 * 60 / g.pix))
    ones = np.ones_like(ymap, dtype=bool)
    S = np.full((len(iy), nb), np.nan, np.float32)
    for r0 in range(0, g.ny, strip):
        r1 = min(r0 + strip, g.ny)
        A, _ = annulus_maps(ymap, ones, g, edges, edges[-1], (r0, r1), pad)
        s = np.where((iy >= r0) & (iy < r1))[0]
        S[s] = A[:, iy[s] - r0, ix[s]].T
    del ymap
    # the own template as painted and sampled on an empty stamp, for the groups in the analysis range (others: nan)
    sel = np.where((iy >= 0) & (iy < g.ny) & (z >= CATS[args.cat]['zlo']) & (z < CATS[args.cat]['zhi']))[0]
    # own is deterministic in (M, z), so a random subsample of the in-range groups gives each bin's mean own without noise
    # worth the name; the fit takes the bin mean of S over all groups minus the bin mean of own over the subsample
    rng = np.random.default_rng(7)
    if len(sel) > 300000:
        sel = np.sort(rng.choice(sel, size=300000, replace=False))
    own = np.full((len(iy), nb), np.nan, np.float32)
    e = np.asarray(edges)
    for cnt, i in enumerate(sel):
        rmax = max(5 * t500[i] / ARCMIN + 1.5 * fwhm, edges[-1] + 2 * pix)
        ry = int(np.ceil(rmax / pix)); cosd = max(np.cos(np.radians(c['dec'][i])), 0.2)
        rx = int(np.ceil(rmax / pix / cosd))
        ys = np.arange(iy[i] - ry, iy[i] + ry + 1); xs = np.arange(ix[i] - rx, ix[i] + rx + 1)
        dec_c = np.radians(c['dec'][i]); d2 = np.radians(g.dec_of_row(ys))[:, None]
        dra = np.radians(((g.ra_of_col(xs % g.nx) - c['ra'][i] + 180) % 360 - 180))[None, :]
        cosd_ = np.sin(dec_c) * np.sin(d2) + np.cos(dec_c) * np.cos(d2) * np.cos(dra)
        dist = np.arccos(np.clip(cosd_, -1, 1))
        prof = amp[i] * tg.real_profile(t500[i], dist.ravel()).reshape(dist.shape)
        # distances from the pixel centre the sampling used (iy, ix)
        dp = np.radians(g.dec_of_row(iy[i])); dra_p = np.radians(((g.ra_of_col(xs % g.nx) - g.ra_of_col(ix[i]) + 180) % 360 - 180))[None, :]
        cosp = np.sin(dp) * np.sin(d2) + np.cos(dp) * np.cos(d2) * np.cos(dra_p)
        distp = np.degrees(np.arccos(np.clip(cosp, -1, 1))) * 60
        k = np.searchsorted(e, distp.ravel(), side='right') - 1
        okk = (k >= 0) & (k < nb)
        num = np.bincount(k[okk], weights=prof.ravel()[okk], minlength=nb); den = np.bincount(k[okk], minlength=nb)
        with np.errstate(invalid='ignore', divide='ignore'):
            own[i] = num / den
        if (cnt + 1) % 200000 == 0:
            log(f'  own {args.cat}/{args.map}: {cnt + 1}/{len(sel)} ({time.time() - t0:.0f}s)')
    np.savez(os.path.join(OUT, f'neigh_{args.map}_{args.cat}.npz'), S=S, own=own)
    log(f'paint done: neighbour term saved for {len(sel)} groups ({time.time() - t0:.0f}s)')


# ------------------------------------------------------------------------------------------------------------------ fit
def gls(d, C, X):
    Ci = np.linalg.inv(C)
    F = X.T @ Ci @ X; cov = np.linalg.inv(F)
    th = cov @ X.T @ Ci @ d
    r = d - X @ th
    return th, cov, float(r @ Ci @ r)


def fit(args):
    co = cosmo()
    m = MAPS[args.map]; spec = CATS[args.cat]
    edges = m['edges']; nbf = int(np.sum(np.asarray(edges) <= spec['rmax'] + 1e-6)) - 1   # annuli used in the fit
    c = load_cat(args.cat)
    s = np.load(os.path.join(OUT, f'samp_{args.map}_{args.cat}.npz'))
    prof = s['prof'][:, :nbf]; f15 = s['f15']; sub = s['sub']; rot = s['rot'][:, :, :nbf]
    nf = os.path.join(OUT, f'neigh_{args.map}_{args.cat}.npz')
    if not os.path.exists(nf) and args.map.startswith('act_'):
        nf = os.path.join(OUT, f'neigh_act_{args.cat}.npz')    # the same neighbours, the same beam
    neigh = np.load(nf) if os.path.exists(nf) else None
    tg = Templates(args.map, edges, xoff=args.xoff)
    z = c['z']; M500 = 10 ** c['logM500']; logMh = c['logMh']
    minrich = args.minrich if args.minrich else spec['minrich']
    if args.common:      # R-136: the groups in both footprints (ACT 15'-clean and Planck 15'-clean, finite profiles)
        o = np.load(os.path.join(OUT, f'samp_{args.common}_{args.cat}.npz'))
        f15 = np.where((o['f15'] >= 0.99) & np.all(np.isfinite(o['prof']), axis=1), f15, 0.0)
        del o
    base = (c['rich'] >= minrich) & (z >= spec['zlo']) & (z < spec['zhi']) & (f15 >= 0.99) & np.all(np.isfinite(prof), axis=1) & (c['dec'] > DEC_LO + 0.5) & (c['dec'] < DEC_HI - 0.5)
    log(f'{args.cat}/{args.map}: {base.sum()} groups in range and footprint of {len(z)}; radio-flagged {c["radio"][base].mean():.3f}')
    # restrict every per-group array to the in-range groups
    bi = np.where(base)[0]
    subidx_full = np.cumsum(sub) - 1
    rot_rows = subidx_full[bi]          # row in rot for each in-range group (valid where sub)
    c = {k: v[bi] for k, v in c.items()}; prof = prof[bi]; sub = sub[bi]
    if neigh is not None:
        neigh = dict(S=neigh['S'][bi], own=neigh['own'][bi])
    z = c['z']; M500 = 10 ** c['logM500']; logMh = c['logMh']; base = np.ones(len(bi), bool)
    T, Yt = group_templates(tg, M500, z, co, sig=args.sig)
    Pk = tg.annulus_means(np.array([tg.t500[0]]))[0][:nbf]; Pk = Pk / Pk[0]    # the beam (theta500 0.3' UPP), unit central annulus
    T = T[:, :nbf]
    bins = {}
    e = spec['edges']
    def add(label, sel):
        bins[label] = sel
    for i in range(len(e) - 1):
        mb = base & (logMh >= e[i]) & (logMh < e[i + 1])
        add(f'M{i} all-z', mb)
        if spec['zsplit']:
            for j in range(len(spec['zsplit']) - 1):
                add(f'M{i} z{j}', mb & (z >= spec['zsplit'][j]) & (z < spec['zsplit'][j + 1]))
    rfn = os.path.join(OUT, f'rand_{args.map}.npz') if args.rsub else os.path.join(OUT, f'samp_{args.map}_randoms.npz')
    rand = np.load(rfn)
    rp = rand['prof'][:, :nbf]; rok = (rand['f15'] >= 0.99) & np.all(np.isfinite(rp), axis=1) & (rand['dec'] > DEC_LO + 0.5) & (rand['dec'] < DEC_HI - 0.5)
    if args.common and args.rsub:      # the randoms share positions across maps (seed 136, same count): intersect footprints
        ro = np.load(os.path.join(OUT, f'rand_{args.common}.npz'))
        assert np.allclose(ro['ra'], rand['ra'])
        rok &= (ro['f15'] >= 0.99) & np.all(np.isfinite(ro['prof']), axis=1)
        rp = rand['prof'][:, :nbf]
    rp = rp[rok]
    # R-136 deviation (results-R132.md): subtract the mask-matched random profile, add its covariance; error x sqrt(chi2/dof)
    # the radio-masked stacks take the randoms with the same radio condition (no NVSS/SUMSS > 5 mJy within 3'): the
    # condition alone moves the 0-3' annuli against 3-6' (R-136, results-R132.md D5)
    rradio = radio_flag(rand['ra'][rok], rand['dec'][rok]) if args.rsub else np.zeros(len(rp), bool)
    dRs, CRs = {}, {}
    for variant, rs in [('radio-masked', ~rradio), ('unmasked', np.ones(len(rp), bool))]:
        dRs[variant] = rp[rs].mean(0).astype(float) if args.rsub else np.zeros(nbf)
        CRs[variant] = np.cov(rp[rs].T.astype(float)) / rs.sum() if args.rsub else np.zeros((nbf, nbf))
    dR, CR = dRs['radio-masked'], CRs['radio-masked']
    if args.rsub:
        log(f'{args.cat}/{args.map}: subtracting the mean of {len(rp)} mask-matched randoms')
    out = []
    lines = []
    for label, sel in bins.items():
        row = dict(label=label)
        for variant in ['radio-masked', 'unmasked']:
            ss = sel & (~c['radio'] if variant == 'radio-masked' else True)
            n = int(ss.sum())
            dR, CR = dRs[variant], CRs[variant]
            if n < 30:
                row[variant] = None; continue
            d = np.nanmean(prof[ss], 0).astype(float) - dR
            C = np.cov(prof[ss].T.astype(float)) / n + CR
            T1 = T[ss].mean(0)
            X = [T1]
            if neigh is not None:
                Sm = np.nanmean(neigh['S'][ss][:, :nbf].astype(float), 0)
                ow = neigh['own'][ss][:, :nbf].astype(float); okn = np.all(np.isfinite(ow), axis=1)
                T2 = Sm - ow[okn].mean(0) if okn.sum() >= 20 else None
                if T2 is not None: X.append(T2)
            if args.psrc:
                X.append(Pk)      # R-136 post-hoc: a beam-shaped central source of free amplitude (radio/dust in the centre)
            X.append(np.ones(nbf))
            X = np.array(X).T
            th, cov, chi2 = gls(d, C, X)
            A, sA = th[0], np.sqrt(cov[0, 0])
            B = (th[1], np.sqrt(cov[1, 1])) if X.shape[1] == 3 + int(args.psrc) else (np.nan, np.nan)
            # one-halo only (no neighbour term) and fixed-aperture check
            th1, cov1, chi21 = gls(d, C, np.c_[T1, np.ones(nbf)])
            # rotation nulls on the subsample, scaled to the bin: the GLS A on each rotation with the same X and C of the subsample
            sr = ss & sub
            nsub = int(sr.sum()); infl = 1.0; Anull = np.array([])
            if nsub >= 200:
                ridx = rot_rows[sr]
                R = rot[ridx]                                   # (nsub, NROT, nbf)
                Csub = np.cov(prof[sr].T.astype(float)) / nsub
                Anull = np.array([gls(np.nanmean(R[:, j], 0).astype(float) - dR, Csub, X)[0][0] for j in range(NROT)])
                thsub, covsub, _ = gls(np.nanmean(prof[sr], 0).astype(float) - dR, Csub, X)
                infl = max(1.0, np.std(Anull, ddof=1) / np.sqrt(covsub[0, 0]))
            if args.rsub:
                infl = max(infl, np.sqrt(chi2 / max(nbf - X.shape[1], 1)))
            Yref = Yt[ss].mean()
            # fixed aperture: 0-3' disc minus 3-6' ring (ACT) or 0-10' minus 10-20' (Planck), in the template's units
            ee = np.asarray(edges[:nbf + 1]); area = ee[1:] ** 2 - ee[:-1] ** 2
            if m['kind'] == 'act':
                inn = ee[1:] <= 3.0; outr = (ee[:-1] >= 3.0) & (ee[1:] <= 6.0)
            else:
                inn = ee[1:] <= 10.0; outr = (ee[:-1] >= 10.0) & (ee[1:] <= 20.0)
            w = np.zeros(nbf); w[inn] = area[inn] / area[inn].sum(); w[outr] = -area[outr] / area[outr].sum()
            cap = w @ d; capT = w @ T1; scap = np.sqrt(w @ C @ w)
            Acap = cap / capT
            S = (th[-2], np.sqrt(cov[-2, -2])) if args.psrc else (np.nan, np.nan)
            row[variant] = dict(S=S, n=n, A=A, sA=sA * infl, sA_gls=sA, infl=infl, B=B, c=th[-1], chi2=chi2, dof=nbf - X.shape[1],
                                A1h=th1[0], sA1h=np.sqrt(cov1[0, 0]), Y=A * Yref, sY=sA * infl * Yref, Yref=Yref, Acap=Acap,
                                sAcap=scap / capT * infl, logM500=np.log10(M500[ss]).mean(), logMh=logMh[ss].mean(), z=z[ss].mean(),
                                d=d, T1=T1, C=C, nnull=len(Anull), null_mean=float(np.mean(Anull)) if len(Anull) else np.nan)
        out.append(row)
        r = row['radio-masked']; u = row['unmasked']
        if r is None: continue
        lines.append(f"| {label} | {r['n']} | {r['logM500']:.2f} | {r['z']:.3f} | {r['A']:.3f} ± {r['sA']:.3f} (infl {r['infl']:.2f}) "
                     f"| {r['Y'] * 1e8:.3f} ± {r['sY'] * 1e8:.3f} | B {r['B'][0]:.2f} ± {r['B'][1]:.2f} | A1h {r['A1h']:.3f} ± {r['sA1h']:.3f} "
                     f"| cap {r['Acap']:.3f} ± {r['sAcap']:.3f} | S {r['S'][0] * 1e8:.2f} ± {r['S'][1] * 1e8:.2f} | chi2 {r['chi2']:.1f}/{r['dof']} "
                     f"| unmasked {u['n']} A {u['A']:.3f} ± {u['sA']:.3f} |" if u else '| (no unmasked) |')
    # randoms
    d = rp.mean(0).astype(float); C = np.cov(rp.T.astype(float)) / len(rp)
    log(f'{args.cat}/{args.map}: randoms {len(rp)}: mean y x1e8 per annulus {np.round(d * 1e8, 3)} +- {np.round(np.sqrt(np.diag(C)) * 1e8, 3)}')
    hdr = ('| bin | N | <log M500c> | <z> | A (obs/ref) | Y500 E^-2/3 D_A^2 [1e-8 Mpc^2] | B (neighbours) | A one-halo only | A fixed-aperture | chi2/dof | unmasked |\n'
           '|---|---|---|---|---|---|---|---|---|---|---|')
    print(f'\n### {args.cat} on {args.map}\n' + hdr + '\n' + '\n'.join(lines))
    # the power law over the headline bins (all-z rows with log M_h >= headline_from), radio-masked
    rows = [(row, row['radio-masked']) for row in out if row['label'].endswith('all-z') and row['radio-masked'] is not None]
    hl = [(row, r) for row, r in rows if float(row['label'].split()[0][1:]) >= [i for i, v in enumerate(e[:-1]) if v >= spec['headline_from']][0]]
    for tag, rr in [('headline', hl), ('all bins', rows)]:
        if len(rr) < 3: continue
        lm = np.array([r['logM500'] for _, r in rr]); Y = np.array([r['Y'] for _, r in rr]); sY = np.array([r['sY'] for _, r in rr])
        from scipy.optimize import curve_fit
        f = lambda x, la, b: 10 ** la * 10 ** (b * (x - 13.5))
        try:
            p, pc = curve_fit(f, lm, Y, p0=[np.log10(max(Y[-1], 1e-9)) - 1.8 * (lm[-1] - 13.5), 1.8], sigma=sY, absolute_sigma=True, maxfev=20000)
            chi2 = np.sum(((Y - f(lm, *p)) / sY) ** 2)
            # the reference relation at the same bins
            print(f'power law ({tag}, {len(rr)} bins, radio-masked): log10 A_Y(10^13.5 M500c) = {p[0]:.3f} ± {np.sqrt(pc[0, 0]):.3f} '
                  f'(Y in Mpc^2, E^-2/3 D_A^2 Y500), B_Y = {p[1]:.2f} ± {np.sqrt(pc[1, 1]):.2f}, chi2 {chi2:.1f}/{len(rr) - 2}; '
                  f'reference (Planck 2013, 1-b = 0.8) at 10^13.5: log10 Y = {np.log10(r131.yref(10 ** 13.5, 0.0, co)):.3f}, slope 1.79')
        except Exception as ex:
            print('power law fit failed', ex)
    import pickle
    tag = (f'_common' if args.common else '') + ('_rsub' if args.rsub else '') + ('_psrc' if args.psrc else '') + (f'_x{args.xoff:.2f}' if args.xoff else '') + (f'_rich{args.minrich}' if args.minrich else '') + (f'_sig{args.sig}' if args.sig != SIG_M else '')
    with open(os.path.join(OUT, f'fit_{args.map}_{args.cat}{tag}.pkl'), 'wb') as fh:
        pickle.dump(out, fh)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['cats', 'sample', 'paint', 'fit', 'randoms'])
    ap.add_argument('--map', default='act'); ap.add_argument('--cat', default='desi'); ap.add_argument('--force', action='store_true')
    ap.add_argument('--nrand', type=int, default=3000000); ap.add_argument('--rsub', action='store_true')
    ap.add_argument('--psrc', action='store_true'); ap.add_argument('--xoff', type=float, default=0.0)
    ap.add_argument('--common', default='')
    ap.add_argument('--minrich', type=int, default=0); ap.add_argument('--sig', type=float, default=SIG_M)
    a = ap.parse_args()
    {'cats': cats, 'sample': sample, 'paint': paint, 'fit': fit, 'randoms': sample_randoms}[a.cmd](a)


if __name__ == '__main__':
    main()
