"""CMB-lensing contamination budget from zero-signal templates (R-129, cosmo-deep-9).

The ACT DR6 lensing convergence (five map variants) and the Planck PR3 MV convergence crossed with 24 templates that cannot carry
lensing signal (Gaia DR3 star counts per G bin and colour, Legacy DR9 imaging depth / PSF / dust / HI / stellar density, the Quaia
selection-function templates), per Galactic-latitude band, with the null from the release's 400 simulated reconstructions (and
rotated maps as a check), and the contamination propagated into the kappa x galaxy amplitude of the DESI Legacy LRG bins and
Quaia with the DES statistic X_S = C^{kS} C^{gS} / C^{SS} (Omori+2019 Sec. V, Chang+2022 App. B1). Pseudo-C_L on a common
weight per ratio, decoupled with a binned MASTER matrix computed here (closed-form Wigner 3j, m = 0).

usage (python -I; data under data/, gitignored; the plan is fixed in results-R129.md):
  python -I kappa_templates.py prep                   # templates, galaxy maps, masks, weights at nside 512 -> data/ktemp/
  python -I kappa_templates.py coupling [--map M]     # binned MASTER matrices for every weight of map M (baseline weights shared by the ACT variants)
  python -I kappa_templates.py data --map M           # data band powers: k x S, k x g, g x S, S x S per weight
  python -I kappa_templates.py sims --map baseline    # k_sim x S and k_sim x g per weight, streamed (act_dr6_simstream.sh)
  python -I kappa_templates.py rotations --map M [--nrot 60]
  python -I kappa_templates.py analyse --map M [--null sims|rotations]   # PTEs and the propagation table -> logs/
"""
import argparse, glob, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
ACT = os.path.join(DATA, 'act_dr6')
OUT = os.path.join(DATA, 'ktemp')
NSIDE = 512
LMAX = 1100  # alm and coupling-matrix limit; the top band ends at 800, so the coupling across the edge is inside the matrix
EDGES = np.array([20, 50, 100, 150, 200, 300, 400, 600, 800])
NB = len(EDGES) - 1
LAT = {'all': (0, 91), 'b<30': (0, 30), 'b30-50': (30, 50), 'b>50': (50, 91)}
GAIA = ['g12_16', 'g16_18', 'g18_20', 'g20_21', 'g16_20_blue', 'g16_20_red']
LEGACY = ['STARDENS', 'EBV', 'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'PSFDEPTH_W1', 'GALDEPTH_Z', 'PSFSIZE_G', 'PSFSIZE_R',
          'PSFSIZE_Z']  # HI was in the fixed list; the DR1 LSS maps that carry it cover only the DR1 spectroscopic footprint
QUAIA_T = ['dust', 'stars', 'm10', 'mcs', 'unwise', 'unwisescan', 'mcsunwise']
TEMPLATES = [('gaia_' + g, 'sky') for g in GAIA] + [('ls_' + c, 'legacy') for c in LEGACY] + [('quaia_' + q, 'sky') for q in QUAIA_T]
GALAXIES = ['lrg1', 'lrg2', 'lrg3', 'lrg4', 'quaia']
ACT_VARIANTS = ['baseline', 'tonly', 'polonly', 'cibdeproj', 'galcut040']


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def weights_for(mapname):
    """the weights used with a map: (name, footprint group, latitude band, galaxy or None); the ACT foreground variants share
    the baseline mask, galcut040 and planck have their own."""
    W = []
    for fp in ('sky', 'legacy'):
        for lat in LAT:
            W.append((f'{fp}_{lat}', fp, lat, None))
    for g in GALAXIES:
        W.append((f'legacy_all_{g}', 'legacy', 'all', g))
    return W


def mask_tag(mapname):
    return mapname if mapname in ('galcut040', 'planck') else 'baseline'


# ------------------------------------------------------------------------------------------------------------- prep
def gal_lat_map():
    import healpy as hp
    th, ph = hp.pix2ang(NSIDE, np.arange(hp.nside2npix(NSIDE)))
    r = hp.Rotator(coord=['C', 'G'])
    thg, _ = r(th, ph)
    return np.abs(90 - np.degrees(thg)).astype(np.float32)


def act_mask512(variant):
    import healpy as hp
    fn = os.path.join(ACT, f'mask_{variant}_nside512.npy')
    if os.path.exists(fn):
        return np.load(fn)
    log(f'reading the nside-4096 {variant} mask')
    m = hp.read_map(os.path.join(ACT, f'mask_act_dr6_lensing_v1_healpix_nside_4096_{variant}.fits'), dtype=np.float32)
    m512 = hp.ud_grade(m, NSIDE).astype(np.float32)
    np.save(fn, m512)
    return m512


def planck_mask512():
    """the PR3 lensing mask (Galactic, nside 2048, binary) rotated to equatorial at nside 512 and apodized (1 deg Gaussian,
    squared), so its coupling is as benign as ACT's."""
    import healpy as hp
    fn = os.path.join(DATA, 'planck', 'mask_planck_eq_nside512.npy')
    if os.path.exists(fn):
        return np.load(fn)
    m = hp.read_map(os.path.join(DATA, 'planck', 'COM_Lensing_4096_R3.00', 'mask.fits.gz'))
    m = hp.ud_grade(m, NSIDE)
    th, ph = hp.pix2ang(NSIDE, np.arange(hp.nside2npix(NSIDE)))
    r = hp.Rotator(coord=['C', 'G'])
    thg, phg = r(th, ph)
    meq = m[hp.ang2pix(NSIDE, thg, phg)]
    b = (meq > 0.5).astype(float)
    s = hp.smoothing(b, fwhm=np.radians(1.0), lmax=1500)
    w = np.clip(s, 0, 1) ** 2 * (b > 0)
    np.save(fn, w.astype(np.float32))
    return w.astype(np.float32)


def planck_alm_eq():
    """the PR3 MV klm truncated to LMAX and rotated from Galactic to equatorial."""
    import healpy as hp
    fn = os.path.join(DATA, 'planck', 'MV_klm_eq_lmax850.npy')
    if os.path.exists(fn):
        return np.load(fn)
    a = hp.read_alm(os.path.join(DATA, 'planck', 'COM_Lensing_4096_R3.00', 'MV', 'dat_klm.fits'))
    a = truncate_alm(a, LMAX)
    hp.Rotator(coord=['G', 'C']).rotate_alm(a, inplace=True)
    np.save(fn, a)
    return a


def truncate_alm(alm, lmax):
    import healpy as hp
    lm = hp.Alm.getlmax(len(alm))
    if lm <= lmax:
        return alm
    l, m = hp.Alm.getlm(lm)
    keep = l <= lmax
    out = np.zeros(hp.Alm.getsize(lmax), dtype=np.complex128)
    out[hp.Alm.getidx(lmax, l[keep], m[keep])] = alm[keep]
    return out


def load_gaia(name):
    import healpy as hp
    fn = os.path.join(DATA, 'gaia', f'{name}.csv')
    d = np.loadtxt(fn, delimiter=',', skiprows=1, dtype=np.int64)
    m = np.zeros(hp.nside2npix(256))
    m[d[:, 0]] = d[:, 1]
    m = m[hp.ring2nest(256, np.arange(len(m)))]  # NESTED -> RING: ring pixel r holds nested pixel ring2nest(r)
    return m


def load_legacy():
    """the Legacy Surveys DR9 imaging properties per nside-256 pixel from the DESI dark-time target pixweight file (NESTED;
    the whole DR9 imaging footprint, unlike the DR1 LSS maps, which cover the DR1 spectroscopic footprint only)."""
    from astropy.io import fits
    import healpy as hp
    d = fits.open(os.path.join(DATA, 'desi_lrg', 'pixweight-1-dark.fits'))[1].data
    foot = (np.asarray(d['FRACAREA'], dtype=float) > 0.5)
    for c in ('PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z', 'GALDEPTH_Z'):
        foot &= np.asarray(d[c], dtype=float) > 0
    ring = hp.ring2nest(256, np.arange(len(foot)))  # NESTED -> RING
    maps = {}
    for c in LEGACY:
        v = np.asarray(d[c], dtype=np.float64).copy()
        v[~foot] = 0
        maps[c] = v[ring]
    return maps, foot[ring]


def load_quaia_templates():
    import healpy as hp
    d = os.path.join(DATA, 'quaia', 'templates', 'selection_function_template_maps')
    return {q: np.load(os.path.join(d, f'map_{q}_NSIDE64.npy')).astype(float) for q in QUAIA_T}


def quaia_galaxy():
    """delta_g = n / (nbar s) - 1 at nside 256 with the released selection function s (nside 64), footprint s > 0.5."""
    from astropy.io import fits
    import healpy as hp
    c = fits.open(os.path.join(DATA, 'quaia', 'quaia_G20.5.fits'))[1].data
    s64 = np.asarray(fits.open(os.path.join(DATA, 'quaia', 'selection_function_NSIDE64_G20.5.fits'))[1].data['T']).ravel()
    s = hp.ud_grade(s64, 256)
    n = np.bincount(hp.ang2pix(256, c['ra'], c['dec'], lonlat=True), minlength=hp.nside2npix(256)).astype(float)
    foot = s > 0.5
    nbar = n[foot].sum() / s[foot].sum()
    dg = np.zeros_like(n)
    dg[foot] = n[foot] / (nbar * s[foot]) - 1
    log(f'quaia: {len(c)} sources, footprint fsky {foot.mean():.3f}, nbar per pixel {nbar:.2f}')
    return dg, foot.astype(float)


def lrg_galaxy(i):
    """White et al. 2022's LRG overdensity and mask for photometric bin i (nside 256). They are in Galactic coordinates
    (a third of the mask would lie at |b| < 15 deg read as equatorial; none does read as Galactic), so they are resampled onto
    equatorial nside-512 pixels (nearest pixel)."""
    import healpy as hp
    d = hp.read_map(os.path.join(DATA, 'desi_lrg_photo', 'maps', f'lrg_s0{i}_del.hpx0256.fits'))
    m = hp.read_map(os.path.join(DATA, 'desi_lrg_photo', 'maps', f'lrg_s0{i}_msk.hpx0256.fits'))
    th, ph = hp.pix2ang(NSIDE, np.arange(hp.nside2npix(NSIDE)))
    thg, phg = hp.Rotator(coord=['C', 'G'])(th, ph)
    p = hp.ang2pix(256, thg, phg)
    return d[p], m[p]


def up(m):
    import healpy as hp
    return hp.ud_grade(np.asarray(m, dtype=np.float64), NSIDE)


def prep(args):
    import healpy as hp
    os.makedirs(OUT, exist_ok=True)
    lat = gal_lat_map()
    np.save(os.path.join(OUT, 'abslat.npy'), lat)
    for v in ('baseline', 'galcut040'):
        if os.path.exists(os.path.join(ACT, f'mask_act_dr6_lensing_v1_healpix_nside_4096_{v}.fits')):
            m = act_mask512(v)
            log(f'ACT {v} mask at 512: fsky(>0) {(m > 0).mean():.4f}, <W^2> {(m ** 2).mean():.4f}')
    if os.path.exists(os.path.join(DATA, 'planck', 'COM_Lensing_4096_R3.00', 'mask.fits.gz')):
        m = planck_mask512(); log(f'Planck mask at 512 (equatorial, apodized): fsky(>0) {(m > 0).mean():.4f}, <W^2> {(m ** 2).mean():.4f}')
    # templates at nside 512, standardised over (own footprint) x (ACT baseline footprint)
    act = act_mask512('baseline') > 0
    T = np.zeros((len(TEMPLATES), hp.nside2npix(NSIDE)), dtype=np.float32)
    legacy, lfoot = load_legacy()
    lfoot512 = up(lfoot) > 0.5
    np.save(os.path.join(OUT, 'legacy_foot.npy'), lfoot512)
    quaia_t = load_quaia_templates()
    for i, (name, fp) in enumerate(TEMPLATES):
        if name.startswith('gaia_'):
            if not os.path.exists(os.path.join(DATA, 'gaia', f'{name[5:]}.csv')):
                log(f'template {i:2d} {name:16s} MISSING (Gaia query not finished); zeros, rerun prep'); continue
            m = up(load_gaia(name[5:]))
        elif name.startswith('ls_'):
            m = up(legacy[name[3:]])
        else:
            m = up(quaia_t[name[6:]])
        foot = act & (lfoot512 if fp == 'legacy' else True)
        mu, sd = m[foot].mean(), m[foot].std()
        if not sd > 0:
            log(f'template {i:2d} {name:16s} CONSTANT on the footprint (mean {mu:.4g}); zeros'); continue
        T[i] = np.where(foot, (m - mu) / sd, 0).astype(np.float32)
        lo = foot & (lat < 25); hi = foot & (lat > 50)
        log(f'template {i:2d} {name:16s} footprint fsky {foot.mean():.4f}, mean {mu:.4g}, std {sd:.4g}; standardised mean at |b|<25 '
            f'{T[i][lo].mean():+.2f}, at |b|>50 {T[i][hi].mean():+.2f}')
    np.save(os.path.join(OUT, 'templates.npy'), T)
    # galaxies at nside 512
    G = np.zeros((len(GALAXIES), hp.nside2npix(NSIDE)), dtype=np.float32)
    GM = np.zeros_like(G)
    for i, g in enumerate(GALAXIES):
        if g == 'quaia':
            d, m = quaia_galaxy()
        else:
            d, m = lrg_galaxy(int(g[-1]))
        G[i] = up(d) if g == 'quaia' else d; GM[i] = up(m) if g == 'quaia' else m
        both = (GM[i] > 0) & act
        log(f'galaxy {g}: fsky {(GM[i] > 0).mean():.4f}, in ACT footprint {both.mean():.4f}, delta rms there {G[i][both].std():.3f}')
    np.save(os.path.join(OUT, 'galaxies.npy'), G); np.save(os.path.join(OUT, 'galaxy_masks.npy'), GM)
    log('prep done')


# ---------------------------------------------------------------------------------------------------- weights, coupling
def load_prep():
    P = {}
    P['lat'] = np.load(os.path.join(OUT, 'abslat.npy'))
    P['T'] = np.load(os.path.join(OUT, 'templates.npy'), mmap_mode='r')
    P['G'] = np.load(os.path.join(OUT, 'galaxies.npy'))
    P['GM'] = np.load(os.path.join(OUT, 'galaxy_masks.npy'))
    P['legacy_foot'] = np.load(os.path.join(OUT, 'legacy_foot.npy'))
    return P


def base_mask(mapname):
    return planck_mask512() if mapname == 'planck' else act_mask512(mask_tag(mapname))


def build_weight(P, base, fp, lat, gal):
    w = base.astype(np.float64).copy()
    if fp == 'legacy':
        w *= P['legacy_foot']
    lo, hi = LAT[lat]
    w *= (P['lat'] >= lo) & (P['lat'] < hi)
    if gal is not None:
        w *= P['GM'][GALAXIES.index(gal)]
    return w


def wigner3j00_sq_log(l1, l2, l3):
    """log of (l1 l2 l3; 0 0 0)^2 for integer arrays with l1 + l2 + l3 even and triangle satisfied."""
    from scipy.special import gammaln
    L = l1 + l2 + l3
    g = L // 2
    return (gammaln(L - 2 * l1 + 1) + gammaln(L - 2 * l2 + 1) + gammaln(L - 2 * l3 + 1) - gammaln(L + 2)
            + 2 * (gammaln(g + 1) - gammaln(g - l1 + 1) - gammaln(g - l2 + 1) - gammaln(g - l3 + 1)))


def coupling_matrix(w_l, lmax):
    """MASTER: M_ll' = (2l'+1)/(4pi) sum_L (2L+1) W_L (l l' L; 0 0 0)^2, l, l' <= lmax, from the weight's power spectrum W_L
    (needs W_L to 2 lmax)."""
    M = np.zeros((lmax + 1, lmax + 1))
    lp = np.arange(lmax + 1)
    for l in range(lmax + 1):
        Lmin = np.abs(l - lp); Lmax = l + lp
        Ls = np.arange(0, 2 * lmax + 1)
        LL, LP = np.meshgrid(Ls, lp, indexing='ij')  # (L, l')
        ok = (LL >= Lmin[None, :]) & (LL <= Lmax[None, :]) & (((l + LP + LL) % 2) == 0)
        val = np.zeros_like(LL, dtype=float)
        idx = np.where(ok)
        val[idx] = np.exp(wigner3j00_sq_log(np.full(len(idx[0]), l), LP[idx], LL[idx]))
        M[l] = (2 * lp + 1) / (4 * np.pi) * ((2 * Ls + 1) * w_l[:2 * lmax + 1])[:, None].__mul__(val).sum(0)
    return M


def bin_matrix(lmax):
    P = np.zeros((NB, lmax + 1))
    for b in range(NB):
        s = slice(EDGES[b], EDGES[b + 1])
        P[b, s] = 1.0 / (EDGES[b + 1] - EDGES[b])
    Q = (P > 0).astype(float).T
    return P, Q


def coupling(args):
    import healpy as hp
    P = load_prep()
    base = base_mask(args.map)
    tag = mask_tag(args.map)
    k, n = (int(x) for x in args.part.split('/'))
    for j, (name, fp, lat, gal) in enumerate(weights_for(args.map)):
        fn = os.path.join(OUT, f'coupling_{tag}_{name}.npz')
        if os.path.exists(fn) or j % n != k:
            continue
        w = build_weight(P, base, fp, lat, gal)
        t0 = time.time()
        w_l = hp.anafast(w, lmax=2 * LMAX)
        M = coupling_matrix(w_l, LMAX)
        Pb, Qb = bin_matrix(LMAX)
        Mbb = Pb @ M @ Qb
        np.savez(fn, Mbb=Mbb, Minv=np.linalg.inv(Mbb), fsky=(w ** 2).mean(), w_l=w_l)
        log(f'coupling {tag} {name}: <W^2> {(w ** 2).mean():.4f}, diag(Mbb) / <W^2> {np.round(np.diag(Mbb) / (w ** 2).mean(), 3)} ({time.time() - t0:.0f}s)')


# ----------------------------------------------------------------------------------------------------------- spectra
def read_kappa_alm(mapname):
    import healpy as hp
    if mapname == 'planck':
        return planck_alm_eq()
    a = hp.read_alm(os.path.join(ACT, f'kappa_alm_data_act_dr6_lensing_v1_{mapname}.fits')).astype(np.complex128)
    a[~np.isfinite(a)] = 0
    return truncate_alm(a, LMAX)


def kappa_map(alm, dphi=0.0):
    import healpy as hp
    a = alm
    if dphi != 0.0:
        l, m = hp.Alm.getlm(LMAX)
        a = alm * np.exp(-1j * m * dphi)
    return hp.alm2map(a, NSIDE, lmax=LMAX, verbose=False) if 'verbose' in hp.alm2map.__code__.co_varnames else hp.alm2map(a, NSIDE, lmax=LMAX)


def masked_alm(m, w):
    """alm of (m - <m>_w) w."""
    import healpy as hp
    s = w.sum()
    mu = (m * w).sum() / s if s > 0 else 0.0
    return hp.map2alm((m - mu) * w, lmax=LMAX, iter=0)


def binned(cl):
    Pb, _ = bin_matrix(LMAX)
    return Pb @ cl[:LMAX + 1]


class Fields:
    """per weight: the template and galaxy alms on that weight (cached on disk), the coupling inverse."""

    def __init__(self, mapname):
        self.P = load_prep()
        self.map = mapname
        self.tag = mask_tag(mapname)
        self.base = base_mask(mapname)
        self.W = weights_for(mapname)
        self.w = {}; self.talm = {}; self.galm = {}; self.Minv = {}
        for name, fp, lat, gal in self.W:
            self.w[name] = build_weight(self.P, self.base, fp, lat, gal)
            c = np.load(os.path.join(OUT, f'coupling_{self.tag}_{name}.npz'))
            self.Minv[name] = c['Minv']
            fn = os.path.join(OUT, f'talm_{self.tag}_{name}.npy')
            if os.path.exists(fn):
                self.talm[name] = np.load(fn)
            else:
                t0 = time.time()
                idx = [i for i, (tn, tfp) in enumerate(TEMPLATES) if fp == 'legacy' or tfp == 'sky']
                A = np.zeros((len(TEMPLATES), 0), dtype=np.complex64)
                arr = []
                for i in range(len(TEMPLATES)):
                    if i in idx:
                        arr.append(masked_alm(np.asarray(self.P['T'][i], dtype=np.float64), self.w[name]).astype(np.complex64))
                    else:
                        arr.append(None)
                n = next(a for a in arr if a is not None).shape[0]
                A = np.zeros((len(TEMPLATES), n), dtype=np.complex64)
                for i, a in enumerate(arr):
                    if a is not None:
                        A[i] = a
                np.save(fn, A); self.talm[name] = A
                log(f'template alms on {name}: {len(idx)} templates ({time.time() - t0:.0f}s)')
            if gal is not None:
                self.galm[name] = masked_alm(self.P['G'][GALAXIES.index(gal)].astype(np.float64), self.w[name])
        self.zero = [i for i in range(len(TEMPLATES)) if not np.any(self.P['T'][i])]
        self.tidx = {name: [i for i, (tn, tfp) in enumerate(TEMPLATES) if (fp == 'legacy' or tfp == 'sky') and i not in self.zero]
                     for name, fp, lat, gal in self.W}
        if self.zero:
            log(f'constant templates skipped: {[TEMPLATES[i][0] for i in self.zero]}')

    def cross_all(self, kmap):
        """decoupled band powers of kappa x template (per weight, per template; NaN where the template is not on that weight) and
        kappa x galaxy (per galaxy weight)."""
        import healpy as hp
        KS = np.full((len(self.W), len(TEMPLATES), NB), np.nan)
        KG = np.full((len(self.W), NB), np.nan)
        for j, (name, fp, lat, gal) in enumerate(self.W):
            ka = masked_alm(kmap, self.w[name])
            for i in self.tidx[name]:
                KS[j, i] = self.Minv[name] @ binned(hp.alm2cl(ka, self.talm[name][i].astype(np.complex128), lmax=LMAX))
            if gal is not None:
                KG[j] = self.Minv[name] @ binned(hp.alm2cl(ka, self.galm[name], lmax=LMAX))
        return KS, KG

    def template_spectra(self):
        """g x S and S x S on each galaxy weight (data only)."""
        import healpy as hp
        GS = np.full((len(self.W), len(TEMPLATES), NB), np.nan)
        SS = np.full((len(self.W), len(TEMPLATES), NB), np.nan)
        for j, (name, fp, lat, gal) in enumerate(self.W):
            for i in self.tidx[name]:
                ta = self.talm[name][i].astype(np.complex128)
                SS[j, i] = self.Minv[name] @ binned(hp.alm2cl(ta, lmax=LMAX))
                if gal is not None:
                    GS[j, i] = self.Minv[name] @ binned(hp.alm2cl(self.galm[name], ta, lmax=LMAX))
        return GS, SS


def data(args):
    F = Fields(args.map)
    t0 = time.time()
    k = kappa_map(read_kappa_alm(args.map))
    KS, KG = F.cross_all(k)
    GS, SS = F.template_spectra()
    np.savez(os.path.join(OUT, f'data_{args.map}.npz'), KS=KS, KG=KG, GS=GS, SS=SS, weights=np.array([w[0] for w in F.W]),
             templates=np.array([t[0] for t in TEMPLATES]), galaxies=np.array(GALAXIES), edges=EDGES)
    log(f'data {args.map}: done ({time.time() - t0:.0f}s)')
    for j, (name, fp, lat, gal) in enumerate(F.W):
        if gal is not None:
            log(f'  kappa x {gal} on {name} x1e6: {np.round(KG[j] * 1e6, 2)}')
    j = [w[0] for w in F.W].index('sky_all')
    for i, (tn, tfp) in enumerate(TEMPLATES):
        jj = j if tfp == 'sky' else [w[0] for w in F.W].index('legacy_all')
        log(f'  kappa x {tn:16s} (whole footprint) x1e5: {np.round(KS[jj, i] * 1e5, 2)}')


def sims(args):
    """the same on each simulated reconstruction streamed into data/act_dr6/sims by act_dr6_simstream.sh; each file deleted once
    used."""
    import healpy as hp
    F = Fields(args.map)
    fn = os.path.join(OUT, f'sims_{args.map}.npz')
    done = {}
    if os.path.exists(fn):
        g = np.load(fn)
        done = {int(i): (ks, kg) for i, ks, kg in zip(g['idx'], g['KS'], g['KG'])}
    simdir = os.path.join(ACT, 'sims')
    idle = 0
    while len(done) < args.nsims and idle < args.idle_max:
        files = sorted(glob.glob(os.path.join(simdir, f'kappa_alm_sim_act_dr6_lensing_v1_{args.map}_*.fits')))
        files = [f for f in files if int(f[-9:-5]) not in done]
        if not files:
            lg = os.path.join(ACT, f'simstream_{args.map}.log')
            if os.path.exists(lg) and 'ALLDONE' in open(lg).read():
                break
            time.sleep(10); idle += 1; continue
        idle = 0
        f = files[0]; i = int(f[-9:-5])
        t0 = time.time()
        try:
            a = hp.read_alm(f).astype(np.complex128)
        except Exception as e:  # noqa: BLE001
            log(f'sim {i}: unreadable ({e}); removing'); os.remove(f); continue
        a[~np.isfinite(a)] = 0
        k = kappa_map(truncate_alm(a, LMAX))
        KS, KG = F.cross_all(k)
        done[i] = (KS, KG)
        os.remove(f); open(os.path.join(simdir, f'done_{args.map}_{i:04d}'), 'w').close()
        idx = sorted(done)
        np.savez(fn, idx=np.array(idx), KS=np.array([done[j][0] for j in idx]), KG=np.array([done[j][1] for j in idx]))
        log(f'sim {i}: {len(done)} done ({time.time() - t0:.0f}s); kappa x lrg1 x1e6 {np.round(KG[-5] * 1e6, 2)}')
    log(f'sims: {len(done)} done')


def rotations(args):
    """the data alm turned in right ascension by 360/nrot * k (k = 1..nrot-1 would include 0; we use k = 1..nrot with step
    360/(nrot+1)), the weights and templates fixed: a null that keeps the map's own noise and the templates' structure."""
    F = Fields(args.map)
    alm = read_kappa_alm(args.map)
    fn = os.path.join(OUT, f'rotations_{args.map}.npz')
    KS = []; KG = []; ang = []
    if os.path.exists(fn):
        g = np.load(fn); KS = list(g['KS']); KG = list(g['KG']); ang = list(g['ang'])
    for k in range(1, args.nrot + 1):
        dphi = 2 * np.pi * k / (args.nrot + 1)
        if any(abs(a - dphi) < 1e-9 for a in ang):
            continue
        t0 = time.time()
        km = kappa_map(alm, dphi=dphi)
        ks, kg = F.cross_all(km)
        KS.append(ks); KG.append(kg); ang.append(dphi)
        np.savez(fn, KS=np.array(KS), KG=np.array(KG), ang=np.array(ang))
        log(f'rotation {k}/{args.nrot} ({np.degrees(dphi):.1f} deg, {time.time() - t0:.0f}s)')
    log(f'rotations {args.map}: {len(ang)} done')


# ----------------------------------------------------------------------------------------------------------- analysis
def rotation_factors(mapname):
    """a kappa map turned in right ascension covers only part of each weight W (the ACT footprint is a third of its declination
    band), so a rotation's cross band powers are those of a smaller area: their null spread is sqrt(f) of the data's, with
    f = sum W^2 [M_rot > 0] / sum W^2 and M the map's own mask. Returns f per (rotation, weight), cached."""
    import healpy as hp
    fn = os.path.join(OUT, f'rotfactors_{mapname}.npz')
    R = np.load(os.path.join(OUT, f'rotations_{mapname}.npz'))
    if os.path.exists(fn):
        g = np.load(fn)
        if len(g['f']) == len(R['ang']):
            return g['f']
    P = load_prep(); base = base_mask(mapname)
    W = [build_weight(P, base, fp, lat, gal) for name, fp, lat, gal in weights_for(mapname)]
    th, ph = hp.pix2ang(NSIDE, np.arange(hp.nside2npix(NSIDE)))
    f = np.zeros((len(R['ang']), len(W)))
    for k, dphi in enumerate(R['ang']):
        mrot = (base[hp.ang2pix(NSIDE, th, ph - dphi)] > 0)
        for j, w in enumerate(W):
            f[k, j] = (w ** 2 * mrot).sum() / (w ** 2).sum()
    np.savez(fn, f=f, ang=R['ang'])
    log(f'rotation overlap fractions for {mapname}: mean {f.mean(0).round(2)}')
    return f


def load_null(mapname, null):
    """the null band powers (KS: rotations/sims x weights x templates x bands; KG: x weights x bands), rotations rescaled by
    1/sqrt(f) so their spread is the full footprint's."""
    N = np.load(os.path.join(OUT, f'{null}_{mapname}.npz'))
    KS, KG = N['KS'].copy(), N['KG'].copy()
    if null == 'rotations':
        f = rotation_factors(mapname)
        KS /= np.sqrt(f)[:, :, None, None]; KG /= np.sqrt(f)[:, :, None]
    return KS, KG


def pte_chi2(d, nulls):
    """chi2 of d against the null sample (Hartlap-corrected covariance) and the empirical PTE (fraction of nulls with chi2 >=
    the data's, each null's chi2 computed against the others' covariance approximated by the full one)."""
    n = len(nulls)
    C = np.cov(nulls.T)
    Ci = np.linalg.pinv(C) * (n - len(d) - 2) / (n - 1)
    mu = nulls.mean(0)
    chi = float((d - mu) @ Ci @ (d - mu))
    chin = np.einsum('ij,jk,ik->i', nulls - mu, Ci, nulls - mu)
    return chi, float((chin >= chi).mean()), chin


def published_lrg(g, mapname):
    """Sailer et al. 2024's measured kappa x LRG band powers (Zenodo 12613408) for the LRG bins: (ell, C_L, error)."""
    import json
    fn = os.path.join(DATA, 'sailer', 'zenodo_LRGxPR4-DR6', 'data', 'lrg_cross_pr4+dr6.json')
    if not g.startswith('lrg') or not os.path.exists(fn):
        return None
    d = json.load(open(fn))
    k = 'PR4' if mapname == 'planck' else 'DR6'
    name = f'{k}_LRGz{g[-1]}'
    cl = np.array(d[f'cl_{name}']); cov = np.array(d[f'cov_{name}_{name}'])
    return np.array(d['ell']), cl, np.sqrt(np.diag(cov))


def analyse(args):
    D = np.load(os.path.join(OUT, f'data_{args.map}.npz'))
    weights = list(D['weights']); templates = list(D['templates'])
    NKS, NKG = load_null(args.map, args.null)
    nn = len(NKS)
    lines = [f'# {args.map}: kappa x template band powers against {nn} {args.null}', '',
             f'bands {EDGES.tolist()}; chi2 over {NB} bands, PTE empirical ({args.null}); "z" per band = (data - null mean) / null spread', '',
             '| template | footprint | weight | chi2 | PTE | z per band |', '|---|---|---|---|---|---|']
    res = {}
    for j, wname in enumerate(weights):
        if '_lrg' in wname or '_quaia' in wname:
            continue
        for i, tn in enumerate(templates):
            d = D['KS'][j, i]
            if not np.isfinite(d).all():
                continue
            nulls = NKS[:, j, i]
            if not np.isfinite(nulls).all():
                continue
            chi, pte, _ = pte_chi2(d, nulls)
            z = (d - nulls.mean(0)) / nulls.std(0)
            res[(wname, tn)] = (chi, pte, z)
            lines.append(f'| {tn} | {TEMPLATES[i][1]} | {wname} | {chi:.1f} | {pte:.4f} | {" ".join(f"{v:+.1f}" for v in z)} |')
    lines += ['', '## Propagation into kappa x galaxy amplitudes (DES statistic X_S = C^kS C^gS / C^SS on the common weight)', '']
    for g in GALAXIES:
        j = weights.index(f'legacy_all_{g}')
        ckg = D['KG'][j]; nkg = NKG[:, j]
        sig = nkg.std(0)
        wts = ckg / sig ** 2
        norm = (ckg ** 2 / sig ** 2).sum()
        sigA = norm ** -0.5
        snr = (ckg * nkg.mean(0) * 0 + ckg) @ (ckg / sig ** 2)  # placeholder to keep shape; detection below
        # detection of kappa x g itself: chi2 of ckg against the null (kappa_sim x g has no signal)
        chi, pte, _ = pte_chi2(ckg, nkg)
        lines += [f'### {g}: kappa x g band powers x1e6 {np.round(ckg * 1e6, 2)}, null spread x1e6 {np.round(sig * 1e6, 2)}; '
                  f'chi2 {chi:.1f} against the null (PTE {pte:.4f}); amplitude error sigma_A / A = {sigA:.3f}', '']
        pub = published_lrg(g, args.map)
        if pub is not None:
            ell_pub, cl_pub, err_pub = pub
            centres = 0.5 * (EDGES[1:] + EDGES[:-1])
            mine = np.interp(ell_pub, centres, ckg)
            lines += [f'published (Sailer et al. 2024, {"ACT DR6" if args.map != "planck" else "Planck PR4"}) C^kg x1e6 at L = {ell_pub.astype(int).tolist()}: '
                      f'{np.round(cl_pub * 1e6, 3).tolist()} ± {np.round(err_pub * 1e6, 3).tolist()}; this pipeline interpolated: {np.round(mine * 1e6, 3).tolist()}', '']
        lines += [
                  '| template | X_S per band / C^kg (%) | dA_S / A (%) | null spread of dA_S / A (%) | dA_S / sigma_A | 95 % limit on |dA_S| / A (%) |',
                  '|---|---|---|---|---|---|']
        tot2 = 0.0; mx = 0.0; rows = []; tot2_sig = 0.0
        for i, tn in enumerate(templates):
            cks = D['KS'][j, i]; cgs = D['GS'][j, i]; css = D['SS'][j, i]
            if not (np.isfinite(cks).all() and np.isfinite(css).all() and (css > 0).all()):
                continue
            X = cks * cgs / css
            dA = (X * wts).sum() / norm
            Xn = NKS[:, j, i] * cgs / css
            dAn = (Xn * wts).sum(1) / norm
            sd = dAn.std()
            lim = abs(dA) + 1.64 * sd
            if tn != 'quaia_unwise':  # the unWISE source counts are galaxies: a positive control, not contamination
                tot2 += dA ** 2; mx = max(mx, abs(dA))
                if tn in ('ls_EBV', 'ls_PSFDEPTH_W1', 'quaia_m10'):
                    tot2_sig += dA ** 2
            rows.append((abs(dA), f'| {tn} | {" ".join(f"{v:+.1f}" for v in 100 * X / ckg)} | {100 * dA:+.2f} | {100 * sd:.2f} | {dA / sigA:+.2f} | {100 * lim:.2f} |'))
        rows.sort(key=lambda r: -r[0])
        lines += [r[1] for r in rows]
        lines += ['', f'**{g}: quadrature sum over the 20 zero-signal templates |dA| / A = {100 * np.sqrt(tot2):.2f} % ({np.sqrt(tot2) / sigA:.2f} sigma_A); '
                  f'largest single template {100 * mx:.2f} % ({mx / sigA:.2f} sigma_A); the three templates with PTE < 1/400 (EBV, W1 depth, m10) in quadrature {100 * np.sqrt(tot2_sig):.2f} % ({np.sqrt(tot2_sig) / sigA:.2f} sigma_A)**', '']
    out = os.path.join(HERE, 'logs', f'analyse_{args.map}_{args.null}.md')
    open(out, 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    log(f'written {out}')


EXTRA = ['csfd_ebv', 'sfd_lss', 'ecl_lat']


def extra_templates():
    """post-hoc diagnostics (after unblinding; stated as such): Chiang 2023's CSFD E(B-V) (SFD with its extragalactic, large-
    scale-structure component removed; Zenodo 8207175, Galactic nside 2048), that removed LSS intensity component, and
    |ecliptic latitude|; all resampled to equatorial nside 512 and standardised on the ACT
    baseline x Legacy footprint."""
    import healpy as hp
    from astropy.io import fits
    fn = os.path.join(OUT, 'extra_templates.npy')
    if os.path.exists(fn):
        return np.load(fn)
    P = load_prep(); act = act_mask512('baseline') > 0; foot = act & P['legacy_foot']
    th, ph = hp.pix2ang(NSIDE, np.arange(hp.nside2npix(NSIDE)))
    thg, phg = hp.Rotator(coord=['C', 'G'])(th, ph)
    T = np.zeros((len(EXTRA), len(th)), dtype=np.float32)
    maps = {}
    for key, f in (('csfd_ebv', 'csfd_ebv.fits'), ('sfd_lss', 'lss_intensity.fits')):
        m = hp.read_map(os.path.join(DATA, 'csfd', f))
        m = hp.ud_grade(m, 1024)
        maps[key] = m[hp.ang2pix(1024, thg, phg)]
    the, _ = hp.Rotator(coord=['C', 'E'])(th, ph)
    maps['ecl_lat'] = np.abs(90 - np.degrees(the))
    for i, k in enumerate(EXTRA):
        m = maps[k]; good = foot & np.isfinite(m)
        mu, sd = m[good].mean(), m[good].std()
        T[i] = np.where(good, (m - mu) / sd, 0)
        log(f'extra template {k}: mean {mu:.4g} std {sd:.4g}')
    np.save(fn, T)
    return T


def extra(args):
    """kappa x the post-hoc templates on the legacy_all and legacy_b>50 weights, data and 60 rotations (the simulations were
    deleted after use), and the propagation into the LRG bins and Quaia on their weights."""
    import healpy as hp
    T = extra_templates()
    F = Fields.__new__(Fields)
    F.P = load_prep(); F.map = args.map; F.tag = mask_tag(args.map); F.base = base_mask(args.map)
    names = ['legacy_all', 'legacy_b>50'] + [f'legacy_all_{g}' for g in GALAXIES]
    F.W = [w for w in weights_for(args.map) if w[0] in names]
    F.w = {}; F.Minv = {}; F.talm = {}; F.galm = {}; F.tidx = {}
    for name, fp, lat, gal in F.W:
        F.w[name] = build_weight(F.P, F.base, fp, lat, gal)
        F.Minv[name] = np.load(os.path.join(OUT, f'coupling_{F.tag}_{name}.npz'))['Minv']
        F.talm[name] = np.array([masked_alm(T[i].astype(np.float64), F.w[name]) for i in range(len(EXTRA))])
        F.tidx[name] = list(range(len(EXTRA)))
        if gal is not None:
            F.galm[name] = masked_alm(F.P['G'][GALAXIES.index(gal)].astype(np.float64), F.w[name])
    alm = read_kappa_alm(args.map)
    KS, KG = F.cross_all(kappa_map(alm))
    GS, SS = F.template_spectra()
    RKS = []; RKG = []; ang = []
    for k in range(1, args.nrot + 1):
        dphi = 2 * np.pi * k / (args.nrot + 1)
        ks, kg = F.cross_all(kappa_map(alm, dphi=dphi)); RKS.append(ks); RKG.append(kg); ang.append(dphi)
        if k % 10 == 0:
            log(f'extra rotation {k}/{args.nrot}')
    np.savez(os.path.join(OUT, f'extra_{args.map}.npz'), KS=KS, KG=KG, GS=GS, SS=SS, RKS=np.array(RKS), RKG=np.array(RKG), ang=np.array(ang),
             weights=np.array([w[0] for w in F.W]), templates=np.array(EXTRA))
    # rotation overlap factors for these weights
    f_all = rotation_factors(args.map)
    allnames = [w[0] for w in weights_for(args.map)]
    f = f_all[:, [allnames.index(n) for n in [w[0] for w in F.W]]]
    RKS = np.array(RKS) / np.sqrt(f)[:, :, None, None]; RKG = np.array(RKG) / np.sqrt(f)[:, :, None]
    lines = [f'# post-hoc templates on {args.map}, null = {args.nrot} rotations', '', '| template | weight | chi2 | PTE | z per band |', '|---|---|---|---|---|']
    for j, (name, fp, lat, gal) in enumerate(F.W):
        if gal is not None:
            continue
        for i, tn in enumerate(EXTRA):
            chi, pte, _ = pte_chi2(KS[j, i], RKS[:, j, i]); z = (KS[j, i] - RKS[:, j, i].mean(0)) / RKS[:, j, i].std(0)
            lines.append(f'| {tn} | {name} | {chi:.1f} | {pte:.3f} | {" ".join(f"{v:+.1f}" for v in z)} |')
    lines += ['', '| sample | template | dA_S / A (%) | null spread (%) | dA_S / sigma_A |', '|---|---|---|---|---|']
    for j, (name, fp, lat, gal) in enumerate(F.W):
        if gal is None:
            continue
        ckg = KG[j]; sig = RKG[:, j].std(0); wts = ckg / sig ** 2; norm = (ckg ** 2 / sig ** 2).sum(); sigA = norm ** -0.5
        for i, tn in enumerate(EXTRA):
            X = KS[j, i] * GS[j, i] / SS[j, i]; dA = (X * wts).sum() / norm
            dAn = ((RKS[:, j, i] * GS[j, i] / SS[j, i]) * wts).sum(1) / norm
            lines.append(f'| {gal} | {tn} | {100 * dA:+.2f} | {100 * dAn.std():.2f} | {dA / sigA:+.2f} |')
    out = os.path.join(HERE, 'logs', f'extra_{args.map}.md'); open(out, 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))


def extra2(args):
    """referee's test (pass 1): does a template's own extragalactic content explain its kappa correlation? For S in (W1 depth,
    m10, SFD, CSFD, stellar density): C^{S,lss} with Chiang's LSS/CIB map on the legacy_all and legacy_b>50 weights, and the part
    of C^{kS} that the CIB content predicts, X = C^{k,lss} C^{S,lss} / C^{lss,lss}, against the measured C^{kS} and its
    400-simulation spread; also chi2 of the three detections over bands 2-8 (L >= 50, the ACT baseline's validated range)."""
    import healpy as hp
    T = extra_templates(); P = load_prep(); D = np.load(os.path.join(OUT, 'data_baseline.npz')); N = np.load(os.path.join(OUT, 'sims_baseline.npz'))
    E = np.load(os.path.join(OUT, 'extra_baseline.npz'))
    weights = list(D['weights']); templates = list(D['templates']); ew = list(E['weights']); et = list(E['templates'])
    base = base_mask('baseline')
    lines = ['# referee test: the CIB content of the templates', '', '| template | weight | C^{S,lss} z per band (vs rotations of the LSS map: none; sign and size only, x1e3 of sigma_S sigma_lss) | X_lss per band / C^{kS} (the part of kappa x S the CIB content predicts, %) | C^{kS} z per band (400 sims) | residual (C^{kS} - X_lss) z per band |',
             '|---|---|---|---|---|---|']
    for wname in ('legacy_all', 'legacy_b>50'):
        fp, lat = 'legacy', wname.split('_')[1]
        w = build_weight(P, base, fp, lat, None)
        Minv = np.load(os.path.join(OUT, f'coupling_baseline_{wname}.npz'))['Minv']
        lss = masked_alm(T[et.index('sfd_lss')].astype(np.float64), w)
        c_ll = Minv @ binned(hp.alm2cl(lss, lmax=LMAX))
        j = weights.index(wname); je = ew.index(wname)
        c_kl = E['KS'][je, et.index('sfd_lss')]
        for tn in ('ls_PSFDEPTH_W1', 'quaia_m10', 'ls_EBV', 'ls_STARDENS', 'ls_PSFDEPTH_Z', 'ls_PSFSIZE_Z'):
            i = templates.index(tn)
            ta = masked_alm(np.asarray(P['T'][i], dtype=np.float64), w)
            c_sl = Minv @ binned(hp.alm2cl(ta, lss, lmax=LMAX))
            X = c_kl * c_sl / c_ll
            cks = D['KS'][j, i]; nulls = N['KS'][:, j, i]; sd = nulls.std(0); mu = nulls.mean(0)
            lines.append(f'| {tn} | {wname} | {" ".join(f"{v * 1e3:+.2f}" for v in c_sl)} | {" ".join(f"{100 * v:+.0f}" for v in X / cks)} | {" ".join(f"{v:+.1f}" for v in (cks - mu) / sd)} | {" ".join(f"{v:+.1f}" for v in (cks - X - mu) / sd)} |')
        # CSFD too, from the extra run
        i = et.index('csfd_ebv')
        ta = masked_alm(T[i].astype(np.float64), w); c_sl = Minv @ binned(hp.alm2cl(ta, lss, lmax=LMAX))
        lines.append(f'| csfd_ebv | {wname} | {" ".join(f"{v * 1e3:+.2f}" for v in c_sl)} | — | (rotations) | — |')
    lines += ['', '## chi2 and PTE over bands 2-8 (L = 50-800) against the 400 simulations, whole footprint', '', '| template | chi2 (7 bands) | PTE |', '|---|---|---|']
    for tn, wname in (('ls_EBV', 'legacy_all'), ('ls_PSFDEPTH_W1', 'legacy_all'), ('quaia_m10', 'sky_all'), ('ls_EBV', 'legacy_b>50'), ('ls_PSFDEPTH_W1', 'legacy_b>50'), ('quaia_m10', 'sky_b>50')):
        j = weights.index(wname); i = templates.index(tn)
        chi, pte, _ = pte_chi2(D['KS'][j, i][1:], N['KS'][:, j, i][:, 1:])
        lines.append(f'| {tn} on {wname} | {chi:.1f} | {pte:.4f} |')
    out = os.path.join(HERE, 'logs', 'extra2_baseline.md'); open(out, 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))


def extra3(args):
    """post-gate addendum (referee pass 2): (a) the joint multi-template projection X_joint(b) = C^{kS}^T (C^{SS})^-1 C^{Sg} over
    the 20 zero-signal templates per band (the 20 x 20 template cross-spectra decoupled per pair; the inverse regularised by
    cutting eigenvalues below 1e-3 of the largest), propagated like a single template, with its null from the 400 simulations;
    (b) the noise-bias offset of the single-template statistic, E[C^{kS} C^{gS}] = C^{kg} C^{SS} / nu_b for a template independent
    of kappa and g, nu_b = f_eff sum(2L+1) with f_eff = <W^2>^2 / <W^4>, as a fraction of the amplitude; the single-template
    biases after subtracting it."""
    import healpy as hp
    D = np.load(os.path.join(OUT, 'data_baseline.npz')); N = np.load(os.path.join(OUT, 'sims_baseline.npz'))
    weights = list(D['weights']); templates = list(D['templates'])
    P = load_prep(); base = base_mask('baseline')
    lines = ['# post-gate addendum: the joint projection and the noise-bias offset (baseline map, 400 simulations)', '',
             '| sample | sigma_A / A | noise-bias offset per template, dA/A (%) | mean dA_S/A of the 17 undetected templates before / after the offset (%) | three detected, in quadrature, after the offset (%) | joint 20-template projection dA/A (%), null spread, / sigma_A | joint, the three detected only | joint, all 20 without SFD (CSFD in its place is not available on these weights) |',
             '|---|---|---|---|---|---|---|---|']
    keep = [i for i, tn in enumerate(templates) if tn not in ('quaia_unwise', 'quaia_mcs', 'quaia_mcsunwise')]
    det = [templates.index(t) for t in ('ls_EBV', 'ls_PSFDEPTH_W1', 'quaia_m10')]
    undet = [i for i in keep if i not in det]
    nosfd = [i for i in keep if templates[i] != 'ls_EBV']
    L = np.arange(LMAX + 1)
    for g in GALAXIES:
        name = f'legacy_all_{g}'; j = weights.index(name)
        w = build_weight(P, base, 'legacy', 'all', g)
        f_eff = (w ** 2).mean() ** 2 / (w ** 4).mean()
        nu = np.array([f_eff * ((2 * L[EDGES[b]:EDGES[b + 1]] + 1).sum()) for b in range(NB)])
        Minv = np.load(os.path.join(OUT, f'coupling_baseline_{name}.npz'))['Minv']
        A = np.load(os.path.join(OUT, f'talm_baseline_{name}.npy'))
        ckg = D['KG'][j]; nkg = N['KG'][:, j]; sig = nkg.std(0); wts = ckg / sig ** 2; norm = (ckg ** 2 / sig ** 2).sum(); sigA = norm ** -0.5
        # template cross-spectrum matrix per band
        n = len(keep); SS = np.zeros((NB, n, n))
        for a in range(n):
            for c in range(a, n):
                cl = Minv @ binned(hp.alm2cl(A[keep[a]].astype(np.complex128), A[keep[c]].astype(np.complex128), lmax=LMAX))
                SS[:, a, c] = cl; SS[:, c, a] = cl
        def joint(ks_rows, idx):
            # ks_rows: (ntemp_all, NB) kappa x S; returns dA/A of the joint projection over templates idx
            sel = [keep.index(i) for i in idx]
            X = np.zeros(NB)
            for b in range(NB):
                M = SS[b][np.ix_(sel, sel)]
                ev, U = np.linalg.eigh(M); cut = ev > 1e-3 * ev.max()
                Mi = (U[:, cut] / ev[cut]) @ U[:, cut].T
                X[b] = ks_rows[idx, b] @ Mi @ D['GS'][j, idx, b]
            return (X * wts).sum() / norm
        dj = joint(D['KS'][j], keep); djn = np.array([joint(N['KS'][k, j], keep) for k in range(len(N['KS']))])
        dj3 = joint(D['KS'][j], det); dj3n = np.array([joint(N['KS'][k, j], det) for k in range(len(N['KS']))])
        djns = joint(D['KS'][j], nosfd); djnsn = np.array([joint(N['KS'][k, j], nosfd) for k in range(len(N['KS']))])
        # single-template biases and the offset
        off_b = ckg / nu  # X offset per band for an independent template
        off = (off_b * wts).sum() / norm
        single = {}
        for i in keep:
            X = D['KS'][j, i] * D['GS'][j, i] / D['SS'][j, i]
            single[i] = (X * wts).sum() / norm
        mu_before = np.mean([single[i] for i in undet]); mu_after = mu_before - off
        q3 = np.sqrt(sum((single[i] - off) ** 2 for i in det))
        lines.append(f'| {g} | {sigA:.3f} | {100 * off:+.3f} | {100 * mu_before:+.3f} / {100 * mu_after:+.3f} | {100 * q3:.2f} | {100 * dj:+.2f} ({100 * djn.std():.2f}), {dj / sigA:+.2f} | {100 * dj3:+.2f} ({100 * dj3n.std():.2f}), {dj3 / sigA:+.2f} | {100 * djns:+.2f} ({100 * djnsn.std():.2f}), {djns / sigA:+.2f} |')
        log(f'{g}: offset {100 * off:+.3f} %, joint {100 * dj:+.2f} +- {100 * djn.std():.2f} %')
    out = os.path.join(HERE, 'logs', 'extra3_baseline.md'); open(out, 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))


def summary(args):
    """the pass criterion across maps: a template is 'detected' if its whole-footprint PTE against the 400 simulations on the
    baseline map is < 1/400 and no rotation exceeds the data on at least one of tonly, polonly, galcut040; 'marginal' if any
    (template, latitude band) has PTE < 0.01 on the baseline; counts against the chance expectation. Also the detections per
    map against its own rotations."""
    rows = []
    def ptes(mapname, null):
        D = np.load(os.path.join(OUT, f'data_{mapname}.npz')); NKS, _ = load_null(mapname, null)
        weights = list(D['weights']); templates = list(D['templates'])
        out = {}
        for j, wname in enumerate(weights):
            if '_lrg' in wname or '_quaia' in wname:
                continue
            for i, tn in enumerate(templates):
                d = D['KS'][j, i]; nulls = NKS[:, j, i]
                if not (np.isfinite(d).all() and np.isfinite(nulls).all()):
                    continue
                chi, pte, _ = pte_chi2(d, nulls)
                out[(wname, tn)] = (chi, pte, len(nulls))
        return out
    base = ptes('baseline', 'sims')
    nb = next(iter(base.values()))[2]
    variants = {}
    for v in ('tonly', 'polonly', 'galcut040', 'cibdeproj', 'planck', 'baseline'):
        if os.path.exists(os.path.join(OUT, f'data_{v}.npz')) and os.path.exists(os.path.join(OUT, f'rotations_{v}.npz')):
            variants[v] = ptes(v, 'rotations')
    lines = [f'# Summary: the pass criterion ({nb} simulations on the baseline; rotations on the variants)', '',
             '| template | whole-footprint PTE (baseline, sims) | PTE per latitude band (b<30, 30-50, >50) | rotations PTE: ' + ', '.join(variants) + ' | verdict |',
             '|---|---|---|---|---|']
    nmarg = 0; ntests = 0; detected = []
    for tn, fp in TEMPLATES:
        wall = f'{fp}_all'
        if (wall, tn) not in base:
            continue
        chi, pte, _ = base[(wall, tn)]
        lat = [base.get((f'{fp}_{b}', tn), (np.nan, np.nan, 0))[1] for b in ('b<30', 'b30-50', 'b>50')]
        ntests += 1 + sum(np.isfinite(lat)); nmarg += int(pte < 0.01) + sum(int(p < 0.01) for p in lat if np.isfinite(p))
        vp = [variants[v].get((wall, tn), (np.nan, np.nan, 0))[1] for v in variants]
        confirm = any(variants[v].get((wall, tn), (0, 1, 0))[1] < 1 / variants[v].get((wall, tn), (0, 1, 60))[2] + 1e-9 for v in ('tonly', 'polonly', 'galcut040') if v in variants)
        if pte < 1.0 / nb - 1e-9 and confirm:  # strictly: no simulation above the data
            verdict = 'DETECTED'; detected.append(tn)
        elif pte < 0.01:
            verdict = f'marginal ({int(round(pte * nb))} of {nb} simulations above the data' + ('; variants at their floor)' if confirm else '; not confirmed by a variant)')
        else:
            verdict = 'upper limit'
        if tn == 'quaia_unwise':
            verdict += ' (galaxy signal, not a zero-signal template)'
        lines.append(f'| {tn} | {pte:.4f} (chi2 {chi:.1f}) | ' + ' '.join(f'{p:.3f}' for p in lat) + ' | ' + ' '.join(f'{p:.3f}' for p in vp) + f' | {verdict} |')
    lines += ['', f'{ntests} (template, latitude band) tests on the baseline map; {nmarg} with PTE < 0.01 against {0.01 * ntests:.1f} expected by chance.',
              f'Detected (PTE < 1/{nb} on the whole footprint and no rotation above the data on tonly, polonly or galcut040): {detected if detected else "none"}.']
    out = os.path.join(HERE, 'logs', 'summary.md')
    open(out, 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('cmd', choices=['prep', 'coupling', 'data', 'sims', 'rotations', 'analyse', 'summary', 'extra', 'extra2', 'extra3'])
    p.add_argument('--map', default='baseline')
    p.add_argument('--nsims', type=int, default=400)
    p.add_argument('--idle-max', type=int, default=360)
    p.add_argument('--nrot', type=int, default=60)
    p.add_argument('--null', default='sims')
    p.add_argument('--part', default='0/1', help='coupling: compute weights j with j %% n == k, as k/n')
    args = p.parse_args()
    globals()[args.cmd](args)


if __name__ == '__main__':
    main()
