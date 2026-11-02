"""Loaders shared by the Cosmicflows-4 catalogue-level scripts (R-119): VizieR J/ApJ/944/94 tables fetched as TSV from
vizier.idia.ac.za (`data/cf4_table2.tsv`, `data/cf4_groups.tsv`), and dust maps at galaxy positions."""
import io, numpy as np, pandas as pd

H0_CF4 = 74.6  # the catalogue's own Hubble constant (Tully+2023)
LINKS = {'Sloan', 'Simbad', 'NED', 'LEDA'}

def vizier_tsv(path):
    raw = [l for l in open(path).read().split('\n') if not l.startswith('#')]
    i = [k for k, l in enumerate(raw) if l.startswith('recno')][0]
    d = pd.read_csv(io.StringIO('\n'.join([raw[i]] + [l for l in raw[i + 3:] if l.strip()])), sep='\t',
                    low_memory=False, dtype=str)
    d.columns = [c.strip() for c in d.columns]
    for c in d.columns:
        if c not in LINKS: d[c] = pd.to_numeric(d[c].str.strip(), errors='coerce')
    return d

def hubble_dm(vcmb, h0=H0_CF4):
    """distance modulus of the Hubble law at CMB-frame velocity vcmb (km/s), the catalogue's own convention (D = V/H0)."""
    return 5 * np.log10(vcmb / h0) + 25

def fp_survey(d):
    """SDSS (north) vs 6dFGSv (south) membership of FP galaxies by sky and depth: 6dFGSv is Dec < 0, |b| > 10, cz < ~16,500;
    SDSS reaches z = 0.1 and dips below Dec = 0 only near the equator. Returns 'sdss', '6df' or 'eq' (the ambiguous strip)."""
    s = np.where(d.DEJ2000 > 0, 'sdss', np.where((d.DEJ2000 < -3) & (d.Vcmb < 16500), '6df', 'eq'))
    s = np.where((d.DEJ2000 < 0) & (d.Vcmb >= 16500), 'sdss', s)  # too deep for 6dFGSv: SDSS southern stripes
    return s

def unit_eq(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], -1)

R_GAL = np.array([[-0.0548755604, -0.8734370902, -0.4838350155], [0.4941094279, -0.4448296300, 0.7469822445],
                  [-0.8676661490, -0.1980763734, 0.4559837762]])  # equatorial -> galactic

def gal_lb(v):
    g = v @ R_GAL.T
    return np.degrees(np.arctan2(g[:, 1], g[:, 0])) % 360, np.degrees(np.arcsin(np.clip(g[:, 2], -1, 1)))

_CACHE = {}

def sfd_ebv(l, b, datadir='data/dustmaps'):
    from astropy.coordinates import SkyCoord
    if 'sfd' not in _CACHE:
        from dustmaps.config import config; config['data_dir'] = datadir
        from dustmaps.sfd import SFDQuery; _CACHE['sfd'] = SFDQuery()
    return _CACHE['sfd'](SkyCoord(l, b, unit='deg', frame='galactic'))

def gnilc_ebv(l, b, path='data/dustmaps/planck/COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits'):
    """Planck 2016 GNILC E(B-V) = 1.49e4 x tau353 (Planck 2014 XI calibration), nearest nside-2048 pixel (RING)."""
    import healpy as hp
    from astropy.io import fits
    if 'gnilc' not in _CACHE: _CACHE['gnilc'] = np.asarray(fits.getdata(path, 1)['TAU353'], dtype=float).ravel()
    return 1.49e4 * _CACHE['gnilc'][hp.ang2pix(2048, l, b, lonlat=True)]

def zhou_ebv(l, b, path='data/zhou/desi_dust_gr_256.fits', min_stars=1):
    """Zhou+2024 DESI stellar reddening map (g-r based, nside 256 RING, EBV_GR on the SFD scale); NaN where the pixel has
    fewer than min_stars stars (inpainted pixels) or is outside the footprint. Also returns the map's own SFD sampling."""
    import healpy as hp
    from astropy.io import fits
    key = ('zhou', path, min_stars)
    if key not in _CACHE:
        t = fits.getdata(path, 1)
        m = np.full(hp.nside2npix(256), np.nan); s = np.full(hp.nside2npix(256), np.nan)
        ok = t['N_STAR'] >= min_stars
        m[t['HPXPIXEL'][ok]] = t['EBV_GR'][ok]; s[t['HPXPIXEL'][ok]] = t['EBV_SFD'][ok]; _CACHE[key] = (m, s)
    m, s = _CACHE[key]; p = hp.ang2pix(256, l, b, lonlat=True)
    return m[p], s[p]
