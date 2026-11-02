"""R-128 (cosmo-deep-8): is the de-magnified M_UV, stellar mass and size of z > 8 galaxies behind Abell 2744 independent of the
magnification, beyond what the selection predicts? The sample, statistics and forward model are fixed in results-R128.md.

Usage (from research/cosmo-discovery, with a Python with the scripts' dependencies):
    python jwst_mu_test.py prep              # build the per-field samples -> data/r128/<field>.fits, footprint + depth cells
    python jwst_mu_test.py fit [--snr 3]     # R1-R3, T4 on the UNCOVER A2744 sample with bootstrap, lens and photo-z draws
    python jwst_mu_test.py forward           # the forward model: expected slopes under a correct lens model, no contamination
    python jwst_mu_test.py counts            # C5: lensed-vs-blank de-magnified bright-end counts
    python jwst_mu_test.py all
Data: data/uncover (UNCOVER DR3 SUPER, DR4 magnifications, DR4 SPS, v2.0 lens maps), data/dja (grizli v7 + eazy-py).
"""
import sys, os, json, glob, argparse
import numpy as np
from astropy.table import Table, join
from astropy.io import fits
from astropy.wcs import WCS
from astropy.cosmology import FlatLambdaCDM
from scipy import stats
from scipy.special import erf

COSMO = FlatLambdaCDM(H0=67.7, Om0=0.31)
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, 'data')
OUT = os.path.join(D, 'r128'); os.makedirs(OUT, exist_ok=True)
RPSF = 0.08            # arcsec, half-light radius of the F444W-matched detection PSF
MU_MIN, MU_MAX = 1.2, 30.0
ZMIN, Z16MIN = 8.0, 7.0
CELL = 6.0 / 3600.     # deg, footprint cell
A2744 = (3.5875, -30.40)
rng = np.random.default_rng(128)


def dm(z):
    return 5 * np.log10(COSMO.luminosity_distance(z).to('pc').value / 10.)


def muv(mab, z, mu):
    """flat-f_nu K-correction; mu >= 1."""
    return mab - dm(z) + 2.5 * np.log10(1 + z) + 2.5 * np.log10(mu)


def split_normal(lo, mid, hi, size):
    """draws from a two-sided normal with sigma_lo = mid-lo, sigma_hi = hi-mid (arrays broadcast)."""
    u = rng.standard_normal(size)
    return np.where(u < 0, mid + u * (mid - lo), mid + u * (hi - mid))


# ----------------------------------------------------------------------------------------------------------------- lens maps
class LensMap:
    def __init__(self):
        self.h8 = fits.open(os.path.join(D, 'uncover/lens_v2.0/z08-magnif.fits.gz'))[0]
        self.h10 = fits.open(os.path.join(D, 'uncover/lens_v2.0/z010-magnif.fits.gz'))[0]
        self.wcs = WCS(self.h8.header)
        self.mu9 = 0.5 * (np.abs(self.h8.data) + np.abs(self.h10.data))   # |mu| at z ~ 9

    def at(self, ra, dec):
        x, y = self.wcs.all_world2pix(ra, dec, 0)
        xi, yi = np.round(x).astype(int), np.round(y).astype(int)
        ok = (xi >= 0) & (yi >= 0) & (xi < self.mu9.shape[1]) & (yi < self.mu9.shape[0])
        out = np.full(len(ra), np.nan)
        out[ok] = self.mu9[yi[ok], xi[ok]]
        return out

    def cell_pixels(self, ra0, dec0, n=40):
        """|mu| values of n pixels spread through the 6" cell at (ra0, dec0); NaN outside the map."""
        off = (rng.random((n, 2)) - 0.5) * CELL
        ra = ra0 + off[:, 0] / np.cos(np.deg2rad(dec0)); dec = dec0 + off[:, 1]
        return self.at(ra, dec)


# ----------------------------------------------------------------------------------------------------------------- samples
def uv_band(z):
    return np.where(z < 10.5, 'f150w', 'f200w')


def build_uncover():
    t = Table.read(os.path.join(D, 'uncover/UNCOVER_DR3_LW_SUPER_catalog.fits'))
    m = Table.read(os.path.join(D, 'uncover/uncover-DR4-super-magnifications.fits'))
    s = Table.read(os.path.join(D, 'uncover/UNCOVER_DR4_SPS_catalog.fits'))
    assert np.array_equal(t['id'], m['id']) and np.array_equal(t['id'], s['id'])
    ok = (t['use_phot'] == 1) & (t['flag_eazy'] == 1)
    z = np.array(t['z_phot']); z16 = np.array(t['z160']); z84 = np.array(t['z840'])
    band = uv_band(z)
    f = np.where(band == 'f150w', t['f_f150w'], t['f_f200w']).astype(float)
    e = np.where(band == 'f150w', t['e_f150w'], t['e_f200w']).astype(float)
    snr = np.where(e > 0, f / np.where(e > 0, e, 1), 0)
    base = ok & (z > ZMIN) & (z16 > Z16MIN) & (snr >= 3) & (f > 0)
    mu = np.array(m['mu'], float); mu16 = np.array(m['mu_low_68'], float); mu84 = np.array(m['mu_high_68'], float)
    # multiple images: drop anything within 0.5"
    mi = Table.read(os.path.join(D, 'uncover/lens_v2.0/A2744_multiple-images_v2.0.fits'))
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(np.array(t['ra'])[base] * u.deg, np.array(t['dec'])[base] * u.deg)
    cm = SkyCoord(np.array(mi['RA'], float) * u.deg, np.array(mi['Dec'], float) * u.deg)
    idx, sep, _ = c.match_to_catalog_sky(cm)
    near = np.zeros(len(t), bool); near[np.where(base)[0][sep.arcsec < 0.5]] = True
    sel = base & ~near
    mab = 28.9 - 2.5 * np.log10(f[sel])
    out = Table()
    out['id'] = t['id'][sel]; out['ra'] = t['ra'][sel]; out['dec'] = t['dec'][sel]
    out['z'] = z[sel]; out['z16'] = z16[sel]; out['z84'] = z84[sel]; out['band'] = band[sel]
    out['mab'] = mab; out['snr'] = snr[sel]
    out['mu'] = mu[sel]; out['mu16'] = mu16[sel]; out['mu84'] = mu84[sel]; out['mu_v11'] = np.array(t['mu'], float)[sel]
    out['muv'] = muv(mab, z[sel], mu[sel])
    out['lmass'] = np.array(s['mstar_50'], float)[sel]; out['lmass16'] = np.array(s['mstar_16'], float)[sel]; out['lmass84'] = np.array(s['mstar_84'], float)[sel]
    out['z_sps'] = np.array(s['z_50'], float)[sel]
    out['fr'] = np.array(t['flux_radius'], float)[sel]
    out['e_uv'] = e[sel]
    out['quad'] = (np.array(out['ra']) > A2744[0]).astype(int) + 2 * (np.array(out['dec']) > A2744[1]).astype(int)
    out['in_mu'] = (out['mu'] >= MU_MIN) & (out['mu'] <= MU_MAX)
    out.write(os.path.join(OUT, 'a2744_uncover.fits'), overwrite=True)
    print(f'UNCOVER A2744: {sel.sum()} with S/N>=3 ({(sel & (snr >= 5)).sum()} at >=5); {near.sum()} near multiple images dropped; '
          f'in 1.2<=mu<=30: {out["in_mu"].sum()} (S/N>=5: {(out["in_mu"] & (out["snr"] >= 5)).sum()})')
    # footprint + depth cells for the forward model: cells with >= 1 use_phot object; depth = median e_f150w of compact use_phot objects
    cells(t[ok], 'a2744_uncover', ecol='e_f150w', compact=np.array(t['iso_area'][ok] < 0.5), lens=LensMap())


def cells(t, name, ecol, compact, lens=None, ra='ra', dec='dec', zp=28.9):
    RA = np.array(t[ra], float); DEC = np.array(t[dec], float)
    dec0 = np.median(DEC)
    ix = np.floor((RA - RA.min()) * np.cos(np.deg2rad(dec0)) / CELL).astype(int)
    iy = np.floor((DEC - DEC.min()) / CELL).astype(int)
    key = ix * 100000 + iy
    uk, inv = np.unique(key, return_inverse=True)
    e = np.array(t[ecol], float)
    rows = []
    for k, u in enumerate(uk):
        w = inv == k
        wc = w & compact & (e > 0)
        med = np.median(e[wc]) if wc.sum() >= 3 else (np.median(e[w & (e > 0)]) if (w & (e > 0)).sum() else np.nan)
        cra = RA.min() + (u // 100000 + 0.5) * CELL / np.cos(np.deg2rad(dec0)); cdec = DEC.min() + (u % 100000 + 0.5) * CELL
        rows.append((cra, cdec, w.sum(), med))
    c = Table(rows=rows, names=['ra', 'dec', 'n', 'e_med'])
    c['mlim5'] = zp - 2.5 * np.log10(5 * c['e_med'])
    c['area'] = CELL ** 2 * 3600 ** 2 * np.ones(len(c))     # arcsec^2 (image plane)
    if lens is not None:
        c['mu9'] = lens.at(np.array(c['ra']), np.array(c['dec']))
        pix = np.array([lens.cell_pixels(r, d) for r, d in zip(c['ra'], c['dec'])])
        c['mu9_pix'] = pix
    c.write(os.path.join(OUT, f'cells_{name}.fits'), overwrite=True)
    good = np.isfinite(c['mlim5'])
    print(f'  cells {name}: {len(c)} cells = {len(c) * CELL ** 2 * 3600:.1f} arcmin^2; 5-sigma depth median {np.nanmedian(c["mlim5"]):.2f} '
          f'(10-90 %: {np.nanpercentile(c["mlim5"][good], 10):.2f}-{np.nanpercentile(c["mlim5"][good], 90):.2f})')


def build_dja():
    lens = LensMap()
    m = Table.read(os.path.join(D, 'uncover/uncover-DR4-super-magnifications.fits'))
    t3 = Table.read(os.path.join(D, 'uncover/UNCOVER_DR3_LW_SUPER_catalog.fits'))
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    cu = SkyCoord(np.array(t3['ra']) * u.deg, np.array(t3['dec']) * u.deg)
    for zf in sorted(glob.glob(os.path.join(D, 'dja/*.eazypy.zout.fits'))):
        field = os.path.basename(zf).split('-fix')[0]
        pf = glob.glob(os.path.join(D, f'dja/{field}-fix_phot_apcorr.fits'))
        if not pf:
            continue
        zt = Table.read(zf); pt = Table.read(pf[0])
        if not np.array_equal(zt['id'], pt['id']):
            zt = join(zt, pt, keys='id', join_type='inner', table_names=['z', 'p'], metadata_conflicts='silent')
            # columns from zt keep their names unless duplicated
            ra = np.array(zt['ra_z']); dec = np.array(zt['dec_z'])
        else:
            zt = join(zt, pt, keys='id', metadata_conflicts='silent', table_names=['z', 'p'])
            ra = np.array(zt['ra_z']); dec = np.array(zt['dec_z'])
        z = np.array(zt['z_phot'], float); z16 = np.array(zt['z160'], float); z84 = np.array(zt['z840'], float)
        band = uv_band(z)
        def col(b, suf):
            return np.array(zt[f'{b}_{suf}'], float) if f'{b}_{suf}' in zt.colnames else np.full(len(zt), np.nan)
        f = np.where(band == 'f150w', col('f150w', 'tot_1'), col('f200w', 'tot_1'))
        e = np.where(band == 'f150w', col('f150w', 'etot_1'), col('f200w', 'etot_1'))
        snr = np.where((e > 0) & np.isfinite(e), f / np.where(e > 0, e, 1), 0)
        flag = np.array(zt['flag'], int)
        ok = (flag < 4) & np.isfinite(z) & (np.array(zt['nusefilt'], int) >= 6)
        sel = ok & (z > ZMIN) & (z16 > Z16MIN) & (snr >= 3) & (f > 0)
        out = Table()
        out['id'] = zt['id'][sel]; out['ra'] = ra[sel]; out['dec'] = dec[sel]
        out['z'] = z[sel]; out['z16'] = z16[sel]; out['z84'] = z84[sel]; out['band'] = band[sel]
        out['mab'] = 23.9 - 2.5 * np.log10(f[sel]); out['snr'] = snr[sel]   # DJA fluxes in microJy
        out['lmass'] = np.log10(np.clip(np.array(zt['mass'], float)[sel], 1, None))
        out['fr'] = np.array(zt['flux_radius'], float)[sel] * 0.04   # DJA pixels are 0.04"
        mu = np.ones(sel.sum()); mu16 = np.ones(sel.sum()); mu84 = np.ones(sel.sum())
        if field.startswith('abell2744clu'):
            c = SkyCoord(out['ra'] * u.deg, out['dec'] * u.deg)
            idx, sep, _ = c.match_to_catalog_sky(cu)
            hit = sep.arcsec < 0.3
            mu[hit] = np.array(m['mu'], float)[idx[hit]]; mu16[hit] = np.array(m['mu_low_68'], float)[idx[hit]]; mu84[hit] = np.array(m['mu_high_68'], float)[idx[hit]]
            mu[~hit] = np.nan
            print(f'  {field}: {hit.sum()} of {len(out)} matched to UNCOVER DR4 within 0.3"')
        out['mu'] = mu; out['mu16'] = mu16; out['mu84'] = mu84
        out['muv'] = muv(np.array(out['mab']), np.array(out['z']), np.where(np.isfinite(mu), mu, 1))
        out.write(os.path.join(OUT, f'dja_{field}.fits'), overwrite=True)
        print(f'DJA {field}: {ok.sum()} ok, {sel.sum()} z>8 (S/N>=3), {(sel & (snr >= 5)).sum()} at S/N>=5')
        ecol = 'f150w_etot_1'
        compact = np.array(zt['flux_radius'], float) < 5
        cells(zt[ok], f'dja_{field}', ecol=ecol, compact=compact[ok], lens=lens if field.startswith('abell2744clu') else None, ra='ra_z', dec='dec_z', zp=23.9)


# ----------------------------------------------------------------------------------------------------------------- estimators
def ols(y, X):
    """y on columns of X (intercept added); returns coefficients."""
    A = np.column_stack([np.ones(len(y))] + [np.asarray(x, float) for x in X])
    beta, *_ = np.linalg.lstsq(A, np.asarray(y, float), rcond=None)
    return beta


def estimators(z, mab, mu, lmass, fr, zcut9=9.0, muv_cut=-20.0, mass_cut=9.5):
    """R1-R3 and T4 on one realisation. Returns dict."""
    lmu = np.log10(mu)
    M = muv(mab, z, mu)
    r = {}
    r['R1'] = ols(M, [lmu])[1]
    r['R1_rho'] = stats.spearmanr(lmu, M).statistic if len(M) > 3 else np.nan
    okm = np.isfinite(lmass)
    r['R2'] = ols(lmass[okm], [M[okm], lmu[okm]])[2] if okm.sum() > 5 else np.nan
    res = np.isfinite(fr) & (fr > 1.1 * RPSF)
    if res.sum() > 8:
        lr = 0.5 * np.log10((fr[res] ** 2 - RPSF ** 2) / mu[res])
        r['R3'] = ols(lr, [M[res], lmu[res]])[2]
        r['R3_n'] = int(res.sum())
    else:
        r['R3'] = np.nan; r['R3_n'] = int(res.sum())
    # T4: stratified shift of log mu for bright / massive z>9 sources vs the z>9 parent at the same observed mab
    p9 = z > zcut9
    for name, flag in [('T4_bright', p9 & (M < muv_cut)), ('T4_massive', p9 & okm & (lmass > mass_cut))]:
        r[name], r[name + '_n'] = strat_shift(lmu, mab, flag, p9)
    return r


def strat_shift(lmu, mab, flag, parent):
    """mean over flagged objects of (log mu - mean log mu of the parent in the same 0.5-mag bin of observed mab)."""
    if flag.sum() < 2:
        return np.nan, int(flag.sum())
    b = np.floor(mab / 0.5)
    d = []
    for i in np.where(flag)[0]:
        w = parent & (b == b[i])
        if w.sum() >= 2:
            d.append(lmu[i] - lmu[w].mean())
    return (np.mean(d) if d else np.nan), int(flag.sum())


def strat_perm(lmu, mab, flag, parent, nperm=10000):
    """permutation p-value (two-sided) of strat_shift: the flag is permuted within the mag bins of the parent."""
    obs, n = strat_shift(lmu, mab, flag, parent)
    if not np.isfinite(obs):
        return obs, n, np.nan
    b = np.floor(mab / 0.5); idx_par = np.where(parent)[0]
    null = []
    for _ in range(nperm):
        f2 = np.zeros_like(flag)
        for bb in np.unique(b[idx_par]):
            ii = idx_par[b[idx_par] == bb]
            k = (flag[ii]).sum()
            if k:
                f2[rng.choice(ii, k, replace=False)] = True
        null.append(strat_shift(lmu, mab, f2, parent)[0])
    null = np.array(null); null = null[np.isfinite(null)]
    p = np.mean(np.abs(null) >= abs(obs)) if len(null) else np.nan
    return obs, n, p


def fit(snr_min=5, tag='', table='a2744_uncover', nboot=1000, ndraw=1000, lensnoise=0.0):
    t = Table.read(os.path.join(OUT, f'{table}.fits'))
    w = np.array(t['in_mu']) & (np.array(t['snr']) >= snr_min) if 'in_mu' in t.colnames else (np.array(t['snr']) >= snr_min) & (t['mu'] >= MU_MIN) & (t['mu'] <= MU_MAX)
    t = t[w]
    z = np.array(t['z']); z16 = np.array(t['z16']); z84 = np.array(t['z84']); mab = np.array(t['mab']); mu = np.array(t['mu'])
    mu16 = np.array(t['mu16']); mu84 = np.array(t['mu84']); lmass = np.array(t['lmass'], float); fr = np.array(t['fr'], float)
    lm16 = np.array(t['lmass16'], float) if 'lmass16' in t.colnames else lmass; lm84 = np.array(t['lmass84'], float) if 'lmass84' in t.colnames else lmass
    n = len(t)
    print(f'\n=== fit {table} S/N>={snr_min}: N={n}, z 8-{z.max():.1f}, mu {mu.min():.2f}-{mu.max():.1f} (median {np.median(mu):.2f}), '
          f'median 68% half-width in log mu {np.median(0.5 * np.log10(mu84 / mu16)):.3f}, M_UV {np.min(muv(mab, z, mu)):.1f} to {np.max(muv(mab, z, mu)):.1f}')
    base = estimators(z, mab, mu, lmass, fr)
    res = {'N': n, 'base': base}
    keys = ['R1', 'R2', 'R3', 'T4_bright', 'T4_massive']
    # bootstrap
    B = {k: [] for k in keys}
    for _ in range(nboot):
        i = rng.integers(0, n, n)
        e = estimators(z[i], mab[i], mu[i], lmass[i], fr[i])
        for k in keys: B[k].append(e[k])
    # lens draws (split-normal in log mu from the 16/84 posterior percentiles), photo-z draws, mass draws, and all together
    L = {k: [] for k in keys}; P = {k: [] for k in keys}; A = {k: [] for k in keys}
    slo = np.sqrt(np.log10(mu / mu16) ** 2 + lensnoise ** 2); shi = np.sqrt(np.log10(mu84 / mu) ** 2 + lensnoise ** 2)
    for _ in range(ndraw):
        dmu = 10 ** split_normal(-slo, 0 * slo, shi, n)
        e = estimators(z, mab, mu * dmu, lmass - np.log10(dmu), fr)       # the catalogue mass was de-magnified with the best-fit mu
        for k in keys: L[k].append(e[k])
        zz = np.clip(split_normal(z16, z, z84, n), 7.0, 20)
        e = estimators(zz, mab, mu, lmass, fr)
        for k in keys: P[k].append(e[k])
        dmu = 10 ** split_normal(-slo, 0 * slo, shi, n); zz = np.clip(split_normal(z16, z, z84, n), 7.0, 20)
        mm = split_normal(lm16, lmass, lm84, n) - np.log10(dmu)
        i = rng.integers(0, n, n)
        e = estimators(zz[i], mab[i], (mu * dmu)[i], mm[i], fr[i])
        for k in keys: A[k].append(e[k])
    def ci(a):
        a = np.array(a, float); a = a[np.isfinite(a)]
        return (np.std(a), np.percentile(a, 16), np.percentile(a, 84), np.percentile(a, 2.5), np.percentile(a, 97.5)) if len(a) > 10 else (np.nan,) * 5
    for k in keys:
        b = base[k]; sb = ci(B[k]); sl = ci(L[k]); sp = ci(P[k]); sa = ci(A[k])
        res[k] = dict(value=b, boot_sd=sb[0], lens_sd=sl[0], pz_sd=sp[0], all_sd=sa[0], all_68=[sa[1], sa[2]], all_95=[sa[3], sa[4]])
        extra = f" (N={base[k + '_n']})" if k + '_n' in base else ''
        print(f'{k:10s} = {b:+.3f}  boot ±{sb[0]:.3f}  lens ±{sl[0]:.3f}  photo-z ±{sp[0]:.3f}  all ±{sa[0]:.3f}  68% [{sa[1]:+.3f}, {sa[2]:+.3f}]  95% [{sa[3]:+.3f}, {sa[4]:+.3f}]{extra}')
    print(f"R1 Spearman rho = {base['R1_rho']:+.3f}; R3 resolved N = {base['R3_n']}")
    # T4 permutation p-values
    lmu = np.log10(mu); M = muv(mab, z, mu); p9 = z > 9
    for name, flag in [('T4_bright', p9 & (M < -20)), ('T4_massive', p9 & np.isfinite(lmass) & (lmass > 9.5))]:
        obs, nn, p = strat_perm(lmu, mab, flag, p9, nperm=5000)
        res[name]['perm_p'] = p; res[name]['n'] = nn; res[name]['n_parent'] = int(p9.sum())
        print(f'{name}: shift {obs:+.3f} dex, N={nn} of {p9.sum()} z>9, permutation p = {p:.3f}')
    # jackknife over quadrants
    if 'quad' in t.colnames:
        q = np.array(t['quad'])
        for k in ['R1', 'R2', 'R3']:
            jk = []
            for qq in np.unique(q):
                w = q != qq
                jk.append(estimators(z[w], mab[w], mu[w], lmass[w], fr[w])[k])
            jk = np.array(jk, float); jk = jk[np.isfinite(jk)]
            sd = np.sqrt((len(jk) - 1) / len(jk) * np.sum((jk - jk.mean()) ** 2)) if len(jk) > 1 else np.nan
            res[k]['jack_sd'] = sd; res[k]['jack_values'] = jk.tolist()
            print(f'{k} jackknife over 4 quadrants: values {np.round(jk, 3).tolist()}, sd {sd:.3f}')
    # v1.1 lens model substituted
    if 'mu_v11' in t.colnames:
        mu11 = np.array(t['mu_v11'])
        e = estimators(z, mab, mu11, lmass + np.log10(mu / mu11), fr)
        res['v11'] = {k: e[k] for k in keys}
        print('v1.1 model instead of v2.0:', {k: round(float(e[k]), 3) for k in keys})
    # by redshift and by mu range
    for lab, w in [('z 8-10', (z <= 10)), ('z > 10', z > 10), ('mu < 3', mu < 3), ('mu >= 3', mu >= 3)]:
        if w.sum() > 8:
            e = estimators(z[w], mab[w], mu[w], lmass[w], fr[w])
            print(f'  subsample {lab}: N={w.sum()}, R1 {e["R1"]:+.3f}, R2 {e["R2"]:+.3f}')
    with open(os.path.join(OUT, f'fit_{table}_snr{snr_min}{tag}{"_ln%.2f" % lensnoise if lensnoise else ""}.json'), 'w') as f:
        json.dump(res, f, indent=1, default=float)
    return res


# ----------------------------------------------------------------------------------------------------------------- forward model
LFS = {}   # filled from results of the literature lookup; set in main


def lf_dpl(M, phistar, Mstar, alpha, beta):
    return phistar / (10 ** (0.4 * (alpha + 1) * (M - Mstar)) + 10 ** (0.4 * (beta + 1) * (M - Mstar)))


def lf_schechter(M, phistar, Mstar, alpha):
    x = 10 ** (-0.4 * (M - Mstar))
    return 0.4 * np.log(10) * phistar * x ** (alpha + 1) * np.exp(-x)


def draw_M(lf, n, Mlo=-24.5, Mhi=-14.5):
    grid = np.linspace(Mlo, Mhi, 2001); p = lf(grid); c = np.cumsum(p); c /= c[-1]
    return np.interp(rng.random(n), c, grid)


def forward(table='a2744_uncover', snr_min=5, nmock=1000, lf_name='dpl', lensnoise=0.0):
    """Mocks of the data's size: sources from the LF, uniform in the source plane, magnified by the v2.0 z=9 map, detected with the
    cell's own 5-sigma depth (0.3-mag soft edge); the lens-model noise of the data's posterior widths added to mu; then the same
    estimators. Returns the distribution of R1-R3 and T4."""
    t = Table.read(os.path.join(OUT, f'{table}.fits'))
    w = np.array(t['in_mu']) & (np.array(t['snr']) >= snr_min); t = t[w]; n_data = len(t)
    zdata = np.array(t['z']); wid = np.sqrt((0.5 * np.log10(np.array(t['mu84']) / np.array(t['mu16']))) ** 2 + lensnoise ** 2)   # the data's posterior half-widths (+ the version scatter)
    c = Table.read(os.path.join(OUT, f'cells_{table}.fits'))
    good = np.isfinite(c['mlim5']) & np.isfinite(c['mu9'])
    c = c[good]
    pix = np.array(c['mu9_pix']); pix = np.where(np.isfinite(pix) & (pix > 0), pix, np.nan)
    # source-plane weight of each cell: sum over its pixels of 1/mu
    wcell = np.nansum(1 / pix, axis=1); wcell /= wcell.sum()
    wpix = np.where(np.isfinite(pix), 1 / pix, 0); cdf = np.cumsum(wpix, axis=1); cdf /= cdf[:, -1:]
    pixf = np.where(np.isfinite(pix), pix, 1.0)
    mlim = np.array(c['mlim5'])
    lf = LFS[lf_name]
    z9 = 9.0; DM9 = dm(z9) - 2.5 * np.log10(1 + z9)
    out = {k: [] for k in ['R1', 'R2', 'R3', 'T4_bright', 'T4_massive', 'R1_true', 'frac_bright', 'R1_nonoise']}
    for _ in range(nmock):
        # draw many, keep the detected ones, then subsample to n_data
        N = 40000
        M = draw_M(lf, N)
        ic = rng.choice(len(c), N, p=wcell)
        # pixel within the cell, weight 1/mu (per-cell CDF over its pixels; NaN pixels get zero weight)
        ipix = np.sum(cdf[ic] < rng.random(N)[:, None], axis=1).clip(0, pix.shape[1] - 1)
        mu = pixf[ic, ipix]
        zz = rng.choice(zdata, N)            # the data's redshift distribution, which the LF-only model does not set
        DMz = dm(zz) - 2.5 * np.log10(1 + zz)
        m_obs = M + DMz - 2.5 * np.log10(mu)
        # noise: the data's typical photometric error at the 5-sigma limit; detection P = Phi((mlim - m)/0.3)
        pdet = 0.5 * (1 + erf((mlim[ic] - m_obs) / (0.3 * np.sqrt(2))))
        det = (rng.random(N) < pdet) & (mu >= MU_MIN) & (mu <= MU_MAX)
        if det.sum() < n_data:
            continue
        i = rng.choice(np.where(det)[0], n_data, replace=False)
        M, mu, zz, m_obs = M[i], mu[i], zz[i], m_obs[i]
        m_obs = m_obs + rng.standard_normal(n_data) * 0.2 / 1.0 * (10 ** (0.4 * (m_obs - mlim[ic[i]])))   # flux noise growing to 0.2 mag at the limit
        # mass and size from the intrinsic M
        lmass_true = 9.0 - 0.4 * (M + 20) + rng.standard_normal(n_data) * 0.4
        kpc = COSMO.arcsec_per_kpc_proper(z9).value
        re = 0.6 * 10 ** (-0.4 * 0.3 * (M + 21)) * 10 ** (rng.standard_normal(n_data) * 0.2) * kpc   # arcsec
        fr = np.sqrt(re ** 2 * mu + RPSF ** 2) + rng.standard_normal(n_data) * 0.01
        # the lens-model noise: the estimated mu differs from the true one by the data's posterior widths
        dlm = rng.standard_normal(n_data) * rng.choice(wid, n_data)
        mu_est = mu * 10 ** dlm
        lmass_est = lmass_true - dlm            # the catalogue de-magnifies with the estimated mu
        e = estimators(zz, m_obs, mu_est, lmass_est, fr)
        for k in ['R1', 'R2', 'R3', 'T4_bright', 'T4_massive']: out[k].append(e[k])
        e0 = estimators(zz, m_obs, mu, lmass_true, fr)
        out['R1_nonoise'].append(e0['R1'])
        out['R1_true'].append(ols(M, [np.log10(mu)])[1])
        out['frac_bright'].append(np.mean(muv(m_obs, zz, mu_est) < -20))
    summ = {}
    print(f'\n=== forward model ({lf_name}, lens noise {lensnoise} dex, N={n_data} per mock, {len(out["R1"])} mocks, cells {len(c)})')
    for k, v in out.items():
        v = np.array(v, float); v = v[np.isfinite(v)]
        if len(v):
            summ[k] = dict(mean=v.mean(), sd=v.std(), p16=np.percentile(v, 16), p84=np.percentile(v, 84), p2=np.percentile(v, 2.5), p97=np.percentile(v, 97.5), n=len(v))
            print(f'{k:12s}: mean {v.mean():+.3f} sd {v.std():.3f}  68% [{np.percentile(v, 16):+.3f}, {np.percentile(v, 84):+.3f}]  95% [{np.percentile(v, 2.5):+.3f}, {np.percentile(v, 97.5):+.3f}]')
    with open(os.path.join(OUT, f'forward_{table}_snr{snr_min}_{lf_name}{"_ln%.2f" % lensnoise if lensnoise else ""}.json'), 'w') as f:
        json.dump(dict(summary=summ, draws={k: list(map(float, v)) for k, v in out.items()}), f)
    return summ


# ----------------------------------------------------------------------------------------------------------------- counts
def eh98_pk(k, h=0.677, Om=0.31, Ob=0.049, ns=0.967, s8=0.81):
    """Eisenstein & Hu 1998 no-wiggle linear P(k) at z=0, k in h/Mpc, normalised to sigma_8."""
    Th = 2.7255 / 2.7; Omh2 = Om * h ** 2; Obh2 = Ob * h ** 2
    s = 44.5 * np.log(9.83 / Omh2) / np.sqrt(1 + 10 * Obh2 ** 0.75)
    alpha = 1 - 0.328 * np.log(431 * Omh2) * Ob / Om + 0.38 * np.log(22.3 * Omh2) * (Ob / Om) ** 2
    def T(k):
        kk = k * h  # 1/Mpc
        G = alpha + (1 - alpha) / (1 + (0.43 * kk * s) ** 4)
        q = kk * Th ** 2 / (Omh2 * G)
        L = np.log(2 * np.e + 1.8 * q); C = 14.2 + 731 / (1 + 62.5 * q)
        return L / (L + C * q ** 2)
    def P(k): return k ** ns * T(k) ** 2
    kk = np.logspace(-4, 2, 3000)
    W = 3 * (np.sin(kk * 8) - kk * 8 * np.cos(kk * 8)) / (kk * 8) ** 3
    s8sq = np.trapezoid(P(kk) * W ** 2 * kk ** 2, kk) / (2 * np.pi ** 2)
    return P(k) * s8 ** 2 / s8sq


def growth(z, Om=0.31):
    """linear growth D(z)/D(0) via Carroll et al. 1992 approximation."""
    def g(a):
        Oz = Om / (Om + (1 - Om) * a ** 3); Lz = 1 - Oz
        return 2.5 * Oz / (Oz ** (4 / 7) - Lz + (1 + Oz / 2) * (1 + Lz / 70))
    a = 1 / (1 + z)
    return a * g(a) / g(1.0)


def sigma_box(Lx, Ly, Lz, z, n=60):
    """rms linear matter fluctuation in a comoving box Lx x Ly x Lz (Mpc/h) at redshift z."""
    k1 = np.logspace(-3, 1.3, n)
    kx, ky, kz = np.meshgrid(k1, k1, k1, indexing='ij')
    k = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
    W = np.sinc(kx * Lx / 2 / np.pi) * np.sinc(ky * Ly / 2 / np.pi) * np.sinc(kz * Lz / 2 / np.pi)
    P = eh98_pk(k)
    # integrate over the positive octant in log k: d^3k = kx ky kz dlnkx dlnky dlnkz; factor 8 for octants / (2pi)^3
    dl = np.log(k1[1] / k1[0])
    s2 = 8 * np.sum(P * W ** 2 * kx * ky * kz) * dl ** 3 / (2 * np.pi) ** 3
    return np.sqrt(s2) * growth(z)


def counts(snr_min=5, bias=8.0):
    """C5: z_phot>8 sources brighter than M_UV=-20 (and -20.5) per source-plane comoving volume, A2744 (DJA, mu from UNCOVER) vs the
    blank fields together; z 8-10 and 10-13."""
    lens = LensMap()
    rows = []
    fields = sorted(glob.glob(os.path.join(OUT, 'dja_*.fits')))
    res = {}
    for zlo, zhi in [(8, 10), (10, 13), (8, 13)]:
        dV = (COSMO.comoving_volume(zhi) - COSMO.comoving_volume(zlo)).value / (4 * np.pi * (180 / np.pi) ** 2 * 3600)   # Mpc^3 per arcmin^2 (image plane at mu=1)
        for Mcut in [-20.0, -20.5, -21.0]:
            tot = {'blank': [0, 0.0, []], 'a2744': [0, 0.0, []]}
            per = {}
            for f in fields:
                field = os.path.basename(f)[4:-5]
                t = Table.read(f); c = Table.read(os.path.join(OUT, f'cells_dja_{field}.fits'))
                lensed = field.startswith('abell2744clu')
                z = np.array(t['z']); mu = np.array(t['mu'], float)
                if lensed:
                    okmu = np.isfinite(mu) & (mu >= MU_MIN) & (mu <= MU_MAX)
                    area = np.nansum(np.where(np.isfinite(c['mu9']) & (c['mu9'] >= MU_MIN) & (c['mu9'] <= MU_MAX), np.nanmean(1 / np.array(c['mu9_pix']), axis=1), 0)) * CELL ** 2 * 3600  # arcmin^2 source plane
                    # cells partly outside the mu window: drop (objects there dropped by okmu)
                else:
                    okmu = np.ones(len(t), bool); area = len(c) * CELL ** 2 * 3600
                M = np.array(t['muv'])
                w = okmu & (np.array(t['snr']) >= snr_min) & (z > zlo) & (z <= zhi) & (M < Mcut)
                key = 'a2744' if lensed else 'blank'
                tot[key][0] += int(w.sum()); tot[key][1] += area * dV; tot[key][2].append((field, int(w.sum()), area, area * dV))
                per[field] = (int(w.sum()), area, area * dV)
            nA, VA, _ = tot['a2744']; nB, VB, _ = tot['blank']
            # Poisson intervals (Gehrels) on the counts; cosmic variance: sigma_cv = b * sigma_DM(box) per field, added in quadrature per field
            def cv(fieldlist):
                s2 = 0
                for field, nn, area, vol in fieldlist:
                    side = np.sqrt(area) * 60 * COSMO.kpc_comoving_per_arcmin(0.5 * (zlo + zhi)).value / 1000 * 0.677 / 60  # Mpc/h per side
                    depth = (COSMO.comoving_distance(zhi) - COSMO.comoving_distance(zlo)).value * 0.677
                    s2 += (bias * sigma_box(side, side, depth, 0.5 * (zlo + zhi)) * nn) ** 2
                return np.sqrt(s2)
            cvA = cv(tot['a2744'][2]); cvB = cv(tot['blank'][2])
            ratio = (nA / VA) / (nB / VB) if nB > 0 and nA > 0 else np.nan
            # interval on the ratio by sampling: gamma for Poisson, lognormal-ish for cosmic variance (relative)
            samp = []
            for _ in range(20000):
                a = rng.gamma(nA + 0.5) * (1 + rng.standard_normal() * (cvA / max(nA, 1))) if nA else 0
                b = rng.gamma(nB + 0.5) * (1 + rng.standard_normal() * (cvB / max(nB, 1)))
                if a > 0 and b > 0: samp.append((a / VA) / (b / VB))
            samp = np.array(samp)
            lo, hi = (np.percentile(samp, [16, 84]) if len(samp) else (np.nan, np.nan)); lo2, hi2 = (np.percentile(samp, [2.5, 97.5]) if len(samp) else (np.nan, np.nan))
            print(f'z {zlo}-{zhi} M_UV<{Mcut}: A2744 N={nA} in {VA:.0f} Mpc^3 (source plane) = {nA / VA * 1e5:.2f}e-5; blank N={nB} in {VB:.0f} Mpc^3 = {nB / VB * 1e5:.2f}e-5; '
                  f'ratio {ratio:.2f} 68% [{lo:.2f}, {hi:.2f}] 95% [{lo2:.2f}, {hi2:.2f}]; cv sigma A {cvA:.1f} B {cvB:.1f} objects')
            for fl, (nn, area, vol) in per.items():
                print(f'     {fl:28s} N={nn:3d} area {area:6.1f} arcmin^2 vol {vol:9.0f} Mpc^3 -> {nn / vol * 1e5 if vol else np.nan:.2f}e-5')
            # field-to-field spread of the blank reference: per-field densities, their log scatter, A2744's rank among them
            dens = {fl: nn / vol for fl, (nn, area, vol) in per.items() if vol > 0}
            bl = np.array([v for k, v in dens.items() if not k.startswith('abell2744clu')]); cl_d = [v for k, v in dens.items() if k.startswith('abell2744clu')]
            if len(bl) and cl_d and np.all(bl > 0):
                lsd = np.std(np.log10(bl)); rank = int(np.sum(bl > cl_d[0]))
                print(f'     blank fields: log10 density scatter {lsd:.2f} dex over {len(bl)} fields (Poisson alone would give {np.sqrt(np.mean(1 / np.array([nn for k, (nn, a, v) in per.items() if not k.startswith("abell2744clu") and nn > 0]))) / np.log(10):.2f}); '
                      f'A2744 exceeded by {rank} of {len(bl)} blank fields; ratio to the median blank field {cl_d[0] / np.median(bl):.2f}')
                res_extra = dict(blank_log_sd=lsd, rank=rank, ratio_to_median=cl_d[0] / np.median(bl))
            else:
                res_extra = {}
            res[f'z{zlo}-{zhi}_M{Mcut}'] = dict(nA=nA, VA=VA, nB=nB, VB=VB, ratio=ratio, ci68=[lo, hi], ci95=[lo2, hi2], cvA=cvA, cvB=cvB, per=per, **res_extra)
    with open(os.path.join(OUT, f'counts_snr{snr_min}.json'), 'w') as f:
        json.dump(res, f, indent=1, default=float)
    return res


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('cmd'); ap.add_argument('--snr', type=float, default=5); ap.add_argument('--nmock', type=int, default=1000)
    ap.add_argument('--table', default='a2744_uncover'); ap.add_argument('--lf', default='dpl'); ap.add_argument('--lensnoise', type=float, default=0.0)
    a = ap.parse_args()
    if a.cmd in ('prep', 'all'):
        build_uncover(); build_dja()
    if a.cmd in ('fit', 'all'):
        fit(snr_min=a.snr, table=a.table, lensnoise=a.lensnoise)
    if a.cmd in ('forward', 'all'):
        forward(table=a.table, snr_min=a.snr, nmock=a.nmock, lf_name=a.lf, lensnoise=a.lensnoise)
    if a.cmd in ('counts', 'all'):
        counts(snr_min=a.snr)


# UV LFs at z ~ 9 (phi* per Mpc^3 per mag), from the papers' tables (results-R128.md row F0):
# Donnan et al. 2024 (arXiv:2403.03171) Table 3, z = 9 DPL: phi* = 23.5e-5, M* = -19.70, alpha = -2.00, beta = -3.81;
# Bouwens et al. 2021 (arXiv:2102.07775) Table 5, z ~ 9 Schechter: phi* = 0.021e-3, M* = -21.15 (fixed), alpha = -2.33.
LFS['dpl'] = lambda M: lf_dpl(M, 23.5e-5, -19.70, -2.00, -3.81)
LFS['schechter'] = lambda M: lf_schechter(M, 0.021e-3, -21.15, -2.33)

if __name__ == '__main__':
    main()
