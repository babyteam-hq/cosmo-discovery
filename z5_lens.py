"""z ~ 5 Lyman-break galaxies (GOLDRUSH r-dropouts) and z ~ 4 (g-dropouts) crossed with CMB lensing (R-146, cosmo-deep-17).

The plan (samples, masks, multipoles, estimator, model, fit, nulls) is fixed in `results-R146.md` before the first
cross-spectrum. Pseudo-C_ell with NaMaster's mode-coupling matrix at nside 1024; ACT DR6 baseline, CIB-deprojected and curl
kappa, the release's 400 simulated baseline reconstructions for the covariance; Planck PR3 MV kappa with Gaussian C+N
realisations; the galaxy auto-spectra with shot noise subtracted, their covariance from correlated Gaussian kappa-galaxy sims.

usage (python -I; data under data/, gitignored; outputs and logs in data/z5/):
  python -I z5_lens.py region                 # the HSC region where kappa is stored (needs data/z5/cat_g.npy)
  python -I z5_lens.py simstore --n 400       # consumer of act_dr6_simstream.sh: ACT sim kappa on the region -> float16 store
  python -I z5_lens.py prep                   # galaxy maps, windows, workspaces, theory, magnification slopes, f_int
  python -I z5_lens.py data                   # data cross- and auto-spectra
  python -I z5_lens.py sims                   # stored ACT sims x every galaxy map; Planck Gaussian sims; joint g-kappa sims
  python -I z5_lens.py fit                    # b, A, significance, nulls -> data/z5/fit.md
"""
import argparse, glob, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
DATA = os.path.join(HERE, 'data')
ACT = os.path.join(DATA, 'act_dr6')
OUT = os.path.join(DATA, 'z5')
NSIDE = 1024
LMAXK = 3000          # kappa alm kept to here before the map (no truncation inside the fitted range)
LMAXN = 2000          # NaMaster lmax (fields, coupling); bins to 1500 used, the overflow dropped
EDGES = {'act': [100, 200, 300, 400, 550, 700, 850, 1000, 1300], 'planck': [100, 150, 200, 300, 400],
         'auto': [100, 200, 300, 400, 550, 700, 850, 1000]}


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


# ------------------------------------------------------------------------------------------------------------- region
def region():
    """nside-1024 pixels within ~2 deg of any g-dropout: where kappa is stored and where the kappa window is non-zero."""
    import healpy as hp
    g = np.load(os.path.join(OUT, 'cat_g.npy'))
    n64 = np.bincount(hp.ang2pix(64, g['ra'], g['dec'], lonlat=True), minlength=hp.nside2npix(64))
    r = n64 > 0
    for _ in range(2):
        nb = hp.get_all_neighbours(64, np.where(r)[0]).ravel()
        r[nb[nb >= 0]] = True
    r1024 = hp.ud_grade(r.astype(float), NSIDE) > 0.5
    pix = np.where(r1024)[0]
    np.save(os.path.join(OUT, 'region_pix.npy'), pix)
    log(f'region: {len(pix)} pixels, {len(pix) * hp.nside2pixarea(NSIDE, degrees=True):.0f} deg2')
    return pix


def region_pix():
    return np.load(os.path.join(OUT, 'region_pix.npy'))


def truncate(alm, lmax):
    import healpy as hp
    lm = hp.Alm.getlmax(len(alm))
    if lm <= lmax:
        return alm
    l, m = hp.Alm.getlm(lm)
    keep = l <= lmax
    out = np.zeros(hp.Alm.getsize(lmax), dtype=np.complex128)
    out[hp.Alm.getidx(lmax, l[keep], m[keep])] = alm[keep]
    return out


def alm_map(alm, lmax=LMAXK):
    import healpy as hp
    a = truncate(np.nan_to_num(alm.astype(np.complex128)), lmax)   # the ACT data alm hold NaN at L = 0, 1
    return hp.alm2map(a, NSIDE, lmax=lmax)


def simstore(args):
    """consume streamed ACT simulated reconstructions: kappa map at nside 1024 on the region, stored as float16."""
    import healpy as hp
    pix = region_pix()
    variant = args.map
    fn = os.path.join(OUT, f'actsims_{variant}.npy')
    store = np.load(fn) if os.path.exists(fn) else np.zeros((args.n, len(pix)), dtype=np.float16)
    donef = os.path.join(OUT, f'actsims_{variant}_done.npy')
    done = set(np.load(donef).tolist()) if os.path.exists(donef) else set()
    simdir = os.path.join(ACT, 'inputs' if variant == 'input' else 'sims')     # 'input': the lensing inputs (transfer)
    pat = 'input_kappa_alm_sim_*.fits' if variant == 'input' else f'kappa_alm_sim_act_dr6_lensing_v1_{variant}_*.fits'
    idle = 0
    while len(done) < args.n and idle < 900:
        files = sorted(glob.glob(os.path.join(simdir, pat)))
        if not files:
            lg = os.path.join(ACT, f'simstream_{variant}.log')
            if os.path.exists(lg) and 'ALLDONE' in open(lg).read():
                break
            time.sleep(2); idle += 1; continue
        idle = 0
        f = files[0]
        i = int(f[-9:-5])
        try:
            alm = hp.read_alm(f)
        except Exception as e:
            log(f'sim {i}: unreadable ({e}); waiting'); time.sleep(8)
            try:
                alm = hp.read_alm(f)
            except Exception:
                os.remove(f); continue
        k = alm_map(alm)[pix]
        store[i - 1] = k.astype(np.float16)     # 0.05 % precision; far outside the mask the sims reach |kappa| ~ 10
        os.remove(f); open(os.path.join(simdir, f'done_{variant}_{i:04d}'), 'w').close()
        done.add(i)
        if len(done) % 20 == 0 or len(done) == args.n:
            np.save(fn, store); np.save(donef, np.array(sorted(done)))
            log(f'{variant}: {len(done)} sims stored')
    np.save(fn, store); np.save(donef, np.array(sorted(done)))
    log(f'{variant}: done, {len(done)} sims')


# ------------------------------------------------------------------------------------------------------------- inputs
# Published inputs, fixed in results-R146.md before any cross-spectrum.
# Selection completeness C(z) (Harikane et al. 2022, Fig. 4, black curves, read off the figure at 0.1 in z)
CZ = {'r': ([4.19, 4.29, 4.39, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 5.10, 5.20, 5.29, 5.39],
            [0.0, 0.01, 0.10, 0.39, 0.52, 0.58, 0.58, 0.60, 0.57, 0.54, 0.21, 0.025, 0.0]),
      'g': ([3.03, 3.19, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50],
            [0.0, 0.01, 0.18, 0.56, 0.61, 0.66, 0.645, 0.62, 0.61, 0.62, 0.61, 0.56, 0.25, 0.04, 0.0])}
# foreground-interloper fraction vs magnitude, Wide layer (Harikane et al. 2022, Table 4; z-band for r-drops, i for g-drops)
FINT = {'r': ([23.0, 24.0, 24.5, 25.0, 25.5, 26.0], [0.45, 0.35, 0.15, 0.09, 0.08, 0.06]),
        'g': ([22.0, 23.0, 24.0, 24.5, 25.0, 25.5, 26.0], [0.60, 0.61, 0.20, 0.11, 0.08, 0.05, 0.04])}
FINT_ERR = 0.10       # prior width on the sample-mean f_int (the COSMOS2015 cross-check allows ~30 % at the faint end)
# interlopers' redshift (Gaussian, the 4000 A / Balmer break in the dropout band pair) and bias priors
ZINT = {'r': (0.75, 0.25), 'g': (0.45, 0.20)}
BINT = {'r': (1.5, 0.5), 'g': (1.3, 0.4)}
# UV luminosity functions (Harikane et al. 2022, Table 6): galaxy DPL (M*, log phi*, alpha, beta) + AGN DPL; and the
# galaxy-Schechter alternative (DPL+Schechter fit)
LF = {'r': dict(z=4.9, gal=(-21.54, -3.63, -2.01, -4.91), agn=(-27.32, -8.35, -1.92, -4.77), sch=(-21.09, -3.16, -1.76),
                agn_s=(-27.67, -8.71, -2.27, -5.92)),
      'g': dict(z=3.8, gal=(-20.99, -3.00, -1.86, -4.77), agn=(-25.69, -6.74, -1.41, -3.24), sch=(-20.49, -2.52, -1.59),
                agn_s=(-27.49, -8.00, -2.20, -5.05))}
BAND = {'r': 'z', 'g': 'i'}
MAGCUT = {'r': 99.0, 'g': 25.5}     # g: Miyatake et al.'s i < 25.5; r: every r-dropout (the catalogue's own z S/N > 5 limit)


def dpl(M, Ms, lphi, a, b):
    x = M - Ms
    return 10 ** lphi / (10 ** (0.4 * (a + 1) * x) + 10 ** (0.4 * (b + 1) * x))


def sch(M, Ms, lphi, a):
    x = M - Ms
    return 10 ** lphi * 10 ** (-0.4 * (a + 1) * x) * np.exp(-10 ** (-0.4 * x))


def lf_slope(m, tag, form='dpl', dz=0.0):
    """s(m) = d log10 Phi / dm of the UV LF (galaxies + AGN) at apparent magnitude m, flat f_nu K-correction."""
    import qso_z4_theory as T
    from astropy.cosmology import FlatLambdaCDM
    c = FlatLambdaCDM(H0=T.H0, Om0=0.315)
    L = LF[tag]; z = L['z'] + dz
    DM = c.distmod(z).value - 2.5 * np.log10(1 + z)
    M = np.asarray(m) - DM
    def phi(M):
        if form == 'dpl':
            return dpl(M, *L['gal']) + dpl(M, *L['agn'])
        return sch(M, *L['sch']) + dpl(M, *L['agn_s'])
    h = 0.01
    return (np.log10(phi(M + h)) - np.log10(phi(M - h))) / (2 * h)


def mags(cat, tag):
    b = BAND[tag]
    return cat[f'{b}_tot'] - cat[f'a_{b}']


def load_cat(tag):
    c = np.load(os.path.join(OUT, f'cat_{tag}.npy'))
    m = mags(c, tag)
    return c[np.isfinite(m) & (m < MAGCUT[tag])]


# ------------------------------------------------------------------------------------------------------------- maps
def coverage():
    """galaxy window: the fraction of each nside-1024 pixel's four nside-2048 sub-pixels holding >= 1 g-dropout (all
    magnitudes, ~4.8 per sub-pixel: an empty sub-pixel is a mask hole or an edge, P(0) ~ 1 % by chance); pixels below 0.5
    dropped. GOLDRUSH ships no randoms; the masked catalogues' holes are the only record of the masks. Field id per pixel
    from the g-dropouts in it."""
    import healpy as hp
    g = np.load(os.path.join(OUT, 'cat_g.npy'))
    p2 = hp.ang2pix(2048, g['ra'], g['dec'], lonlat=True, nest=True)
    occ = np.bincount(p2, minlength=hp.nside2npix(2048)) > 0
    f = occ.reshape(-1, 4).mean(1)                     # nest: the 4 children of a 1024 pixel are consecutive
    f = hp.reorder(f, n2r=True)
    w = np.where(f >= 0.5, f, 0.0)
    p1 = hp.ang2pix(NSIDE, g['ra'], g['dec'], lonlat=True)
    field = np.zeros(hp.nside2npix(NSIDE), dtype=np.int8)
    field[p1] = g['field']
    field[w == 0] = 0
    return w, field


def delta_map(N, w, field):
    d = np.zeros_like(w); nbar = {}
    for f in np.unique(field[w > 0]):
        s = (field == f) & (w > 0)
        nbar[int(f)] = N[s].sum() / w[s].sum()
        d[s] = N[s] / (nbar[int(f)] * w[s]) - 1
    return d, nbar


def act_window(variant='baseline'):
    import healpy as hp
    fn = os.path.join(OUT, f'actmask2_{variant}_ns{NSIDE}.npy')
    if os.path.exists(fn):
        return np.load(fn)
    m = hp.read_map(os.path.join(ACT, f'mask_act_dr6_lensing_v1_healpix_nside_4096_{variant}.fits'), dtype=np.float32)
    m2 = hp.ud_grade(m.astype(np.float64) ** 2, NSIDE)
    m4 = hp.ud_grade(m.astype(np.float64) ** 4, NSIDE)
    np.save(fn, m2); np.save(fn.replace('mask2', 'mask4'), m4)
    return m2


def planck_window():
    import healpy as hp
    fn = os.path.join(OUT, f'planckmask_ns{NSIDE}.npy')
    if os.path.exists(fn):
        return np.load(fn)
    m = hp.read_map(os.path.join(DATA, 'planck', 'COM_Lensing_4096_R3.00', 'mask.fits.gz'))
    th, ph = hp.pix2ang(NSIDE, np.arange(hp.nside2npix(NSIDE)))
    thg, phg = hp.Rotator(coord=['C', 'G'])(th, ph)
    meq = (m[hp.ang2pix(hp.get_nside(m), thg, phg)] > 0.5).astype(float)
    np.save(fn, meq)
    return meq


def kwindows():
    import healpy as hp
    reg = np.zeros(hp.nside2npix(NSIDE)); reg[region_pix()] = 1
    a = act_window('baseline') * reg
    p = planck_window() * reg
    return {'act': a, 'planck': p, 'pin': p * (a > 0), 'pout': p * (a <= 0)}


def planck_kappa():
    import healpy as hp
    fn = os.path.join(OUT, 'planck_kappa_eq.npy')
    if os.path.exists(fn):
        return np.load(fn)
    a = hp.read_alm(os.path.join(DATA, 'planck', 'COM_Lensing_4096_R3.00', 'MV', 'dat_klm.fits'))
    a = truncate(a, LMAXK)
    hp.Rotator(coord=['G', 'C']).rotate_alm(a, inplace=True)
    k = hp.alm2map(a, NSIDE, lmax=LMAXK)
    np.save(fn, k)
    return k


# ------------------------------------------------------------------------------------------------------------- feasibility
def feas(args):
    """step 0: counts, overlaps with ACT DR6 and Planck PR3, and the Gaussian forecast of the cross-power significance."""
    import healpy as hp
    import qso_z4_theory as T
    w, field = coverage()
    K = kwindows()
    pa = hp.nside2pixarea(NSIDE, degrees=True)
    th = T.Theory(lmax=2000)
    nka = np.loadtxt(os.path.join(ACT, 'N_L_kk_act_dr6_lensing_v1_baseline.txt'))
    nk_act = np.zeros(4097); nk_act[nka[:, 0].astype(int)] = nka[:, 1]
    nl = np.loadtxt(os.path.join(DATA, 'planck', 'COM_Lensing_4096_R3.00', 'MV', 'nlkk.dat'))
    nk_pl = np.zeros(4097); nk_pl[nl[:, 0].astype(int)] = nl[:, 1]
    log(f'galaxy window: {np.sum(w > 0) * pa:.1f} deg2 (effective {w.sum() * pa:.1f}); fields ' +
        ', '.join(f'W0{f}: {np.sum(w[field == f]) * pa:.1f}' for f in np.unique(field[w > 0])))
    for tag in ('r', 'g'):
        c = load_cat(tag)
        p = hp.ang2pix(NSIDE, c['ra'], c['dec'], lonlat=True)
        inw = w[p] > 0
        zc = np.linspace(2.5, 6.5, 401)
        nz = np.interp(zc, *CZ[tag], left=0, right=0) * th.res.comoving_radial_distance(zc) ** 2 / th.res.hubble_parameter(zc)
        sp = th.spectra(zc, nz)
        m = mags(c, tag)
        s = np.mean(lf_slope(m[inw], tag))
        f = np.mean(np.interp(m[inw], *FINT[tag]))
        for lens, nk, lmin, lmax in (('act', nk_act, 100, 1000), ('planck', nk_pl, 100, 400)):
            ov = (w > 0) & (K[lens] > 0)
            n_ov = np.sum(ov[p] & inw); area = np.sum(w[ov]) * pa
            nbar = n_ov / (area * (np.pi / 180) ** 2)
            sn = {b: T.forecast(sp, (1 - f) * b, s, nbar, area / 41253, nk, lmin, lmax) for b in (3, 4, 5, 6)}
            log(f'{tag}-drops x {lens}: N {n_ov} on {area:.0f} deg2 (n = {n_ov / area:.0f}/deg2); s = {s:.2f}, f_int = {f:.2f}; '
                'S/N ' + ', '.join(f'b={b}: {v:.1f}' for b, v in sn.items()))
        log(f'{tag}-drops: {len(c)} after the magnitude cut, {inw.sum()} in the galaxy window')


# ------------------------------------------------------------------------------------------------------------- prep
FIELDS = (2, 3, 4, 5)       # the HSC Wide fields under ACT (W06 is outside; W07 is 1 deg2 and dropped from the per-field test)
SAMPLES = ('r', 'g', 'r_ran', 'g_ran')


def bins(kind):
    import pymaster as nmt
    e = EDGES[kind] + [LMAXN + 1]          # the last (overflow) bin is dropped from every output
    return nmt.NmtBin.from_edges(np.array(e[:-1]), np.array(e[1:]))


def galaxy_maps():
    """delta maps (nside 1024) for every sample and the randoms; returns maps, window, field, nbar per field (per sr)."""
    import healpy as hp
    w, field = coverage()
    npix = hp.nside2npix(NSIDE); om = hp.nside2pixarea(NSIDE)
    out = {}; nb = {}
    rng = np.random.default_rng(146)
    for tag in ('r', 'g'):
        c = load_cat(tag)
        p = hp.ang2pix(NSIDE, c['ra'], c['dec'], lonlat=True)
        N = np.bincount(p, minlength=npix).astype(float)
        d, nbar = delta_map(N, w, field)
        out[tag] = d; nb[tag] = {f: v / om for f, v in nbar.items()}
        # randoms: each field's count, Poisson on the window
        Nr = np.zeros(npix)
        for f in nbar:
            s = (field == f) & (w > 0)
            ntot = int(N[s].sum())
            pr = rng.choice(np.where(s)[0], ntot, p=w[s] / w[s].sum())
            Nr += np.bincount(pr, minlength=npix)
        out[f'{tag}_ran'], _ = delta_map(Nr, w, field)
        log(f'{tag}: N = {int(N[w > 0].sum())} in the window; n per field (/deg2): ' +
            ', '.join(f'W0{f}: {v * (np.pi / 180) ** 2:.0f}' for f, v in nb[tag].items()))
    return out, w, field, nb


def theory_tag(tag, dz=0.0):
    """Limber pieces (unit amplitudes) for the LBG n(z) and an interloper grid (z_int means 0.2-1.6)."""
    import qso_z4_theory as T
    th = T.Theory(lmax=LMAXN)
    zc = np.linspace(2.5, 6.5, 401)
    nz = np.interp(zc - dz, *CZ[tag], left=0, right=0) * th.res.comoving_radial_distance(zc) ** 2 / th.res.hubble_parameter(zc)
    WgL, WmuL = th.kernels(zc, nz)
    out = dict(km=th.limber(th.Wk, WgL), kmu=th.limber(th.Wk, WmuL), gg=th.limber(WgL, WgL), gmu=th.limber(WgL, WmuL),
               mumu=th.limber(WmuL, WmuL), kk=th.limber(th.Wk, th.Wk), zmean=np.trapezoid(zc * nz, zc) / np.trapezoid(nz, zc))
    zi = np.round(np.arange(0.2, 1.61, 0.1), 2)
    kc, cc, cmu = [], [], []
    for z0 in zi:
        zz = np.linspace(0.01, 3.0, 300)
        n = np.exp(-0.5 * ((zz - z0) / ZINT[tag][1]) ** 2)
        Wc, _ = th.kernels(zz, n)
        kc.append(th.limber(th.Wk, Wc)); cc.append(th.limber(Wc, Wc)); cmu.append(th.limber(Wc, WmuL))
    out.update(zint=zi, kc=np.array(kc), cc=np.array(cc), cmu=np.array(cmu))
    # sigma8(z) of the fiducial (for A -> sigma8(z))
    import camb
    p = camb.set_params(H0=T.H0, ombh2=T.OMBH2, omch2=T.OMCH2, ns=T.NS, As=T.AS, mnu=0.06, WantTransfer=True, kmax=10,
                        redshifts=[out['zmean'], 0.0])
    r = camb.get_results(p)
    s8 = r.get_sigma8()      # ordered by decreasing z
    out['sigma8_z'] = s8[0]; out['sigma8_0'] = s8[1]
    return out


class Pipe:
    """galaxy fields and workspaces: cross (lens x galaxy window, also per field for ACT) and auto."""
    def __init__(self, samples=SAMPLES, need=('act', 'planck', 'pin', 'pout'), perfield=True, auto=True):
        import pymaster as nmt
        self.nmt = nmt
        g = np.load(os.path.join(OUT, 'maps.npz'))
        self.w = g['w']; self.field = g['field']
        self.samples = samples
        self.gf = {k: nmt.NmtField(self.w, [g[k]], lmax=LMAXN, n_iter=1) for k in samples}
        self.gfw = {}
        if perfield:
            for f in FIELDS:
                wf = self.w * (self.field == f)
                for k in ('r', 'g'):
                    if k in samples:
                        self.gfw[(k, f)] = nmt.NmtField(wf, [g[k]], lmax=LMAXN, n_iter=1)
        self.wk = {k: v for k, v in kwindows().items() if k in need}
        self.ws = {}
        for lens in self.wk:
            gws = [('all', self.w)] + ([(f, self.w * (self.field == f)) for f in FIELDS] if perfield and lens == 'act' else [])
            for name, wg in gws:
                self.ws[(lens, name)] = self._ws(f'ws_{lens}_{name}.fits', self.wk[lens], wg, 'act' if lens != 'planck' and lens != 'pin' and lens != 'pout' else 'planck')
        if auto:
            self.ws[('auto', 'all')] = self._ws('ws_auto_all.fits', self.w, self.w, 'auto')

    def _ws(self, name, w1, w2, kind):
        nmt = self.nmt
        fn = os.path.join(OUT, name)
        ws = nmt.NmtWorkspace()
        if os.path.exists(fn):
            ws.read_from(fn)
        else:
            f1 = nmt.NmtField(w1, None, lmax=LMAXN, spin=0); f2 = nmt.NmtField(w2, None, lmax=LMAXN, spin=0)
            ws.compute_coupling_matrix(f1, f2, bins(kind))
            ws.write_to(fn)
            log(f'workspace {name}')
        return ws

    def kfield(self, kmap, lens):
        return self.nmt.NmtField(self.wk[lens], [kmap], lmax=LMAXN, n_iter=1)

    def cross(self, kmap, lens, samples=None, perfield=False):
        kf = self.kfield(kmap, lens)
        out = {}
        for k in (samples or self.samples):
            out[k] = self.ws[(lens, 'all')].decouple_cell(self.nmt.compute_coupled_cell(kf, self.gf[k]))[0][:-1]
        if perfield:
            for (k, f), gf in self.gfw.items():
                out[f'{k}_W0{f}'] = self.ws[(lens, f)].decouple_cell(self.nmt.compute_coupled_cell(kf, gf))[0][:-1]
        return out

    def auto(self, k, nbar_map=None, gf=None):
        gf = gf or self.gf[k]
        cl = self.nmt.compute_coupled_cell(gf, gf)
        if nbar_map is not None:
            noise = np.mean(np.where(self.w > 0, self.w / np.where(nbar_map > 0, nbar_map, 1), 0.0))
            cl = cl - noise
        return self.ws[('auto', 'all')].decouple_cell(cl)[0][:-1]


def nbar_map(nb, field):
    m = np.zeros(len(field))
    for f, v in nb.items():
        m[field == f] = v
    return m


def prep(args):
    import healpy as hp
    maps, w, field, nb = galaxy_maps()
    np.savez(os.path.join(OUT, 'maps.npz'), w=w, field=field, **maps)
    np.save(os.path.join(OUT, 'nbar.npy'), nb, allow_pickle=True)
    for tag in ('r', 'g'):
        for dz, sfx in ((0.0, ''), (0.1, '_dz')):
            t = theory_tag(tag, dz)
            np.savez(os.path.join(OUT, f'theory_{tag}{sfx}.npz'), **t)
            log(f'theory {tag}{sfx}: <z> {t["zmean"]:.2f}, sigma8(z) {t["sigma8_z"]:.4f}; C^km(300) {t["km"][300]:.3e}, '
                f'C^kmu(300) {t["kmu"][300]:.3e}, C^gg(300) {t["gg"][300]:.3e}')
    # sample-mean f_int and s
    nuis = {}
    c_ = {}
    for tag in ('r', 'g'):
        c = load_cat(tag); p = hp.ang2pix(NSIDE, c['ra'], c['dec'], lonlat=True); m = mags(c, tag)[w[p] > 0]
        fm = np.interp(m, *FINT[tag]); sm = lf_slope(m, tag)
        nuis[tag] = dict(f=float(np.mean(fm)), s=float(np.mean(sm)), E1=float(np.mean((1 - fm) * (5 * sm - 2))),
                         E2=float(np.mean(5 * sm - 2)), E3=float(np.mean(1 - fm)))
        log(f'{tag}: f_int {nuis[tag]["f"]:.3f}, s {nuis[tag]["s"]:.3f}')
    np.save(os.path.join(OUT, 'nuis.npy'), nuis, allow_pickle=True)
    Pipe()
    log('prep done')


# ------------------------------------------------------------------------------------------------------------- data
def data(args):
    import healpy as hp
    P = Pipe()
    nb = np.load(os.path.join(OUT, 'nbar.npy'), allow_pickle=True).item()
    out = {}
    kact = alm_map(hp.read_alm(os.path.join(ACT, 'kappa_alm_data_act_dr6_lensing_v1_baseline.fits')))
    for k, v in P.cross(kact, 'act', perfield=True).items():
        out[f'{k}|act'] = v
    del kact
    kcib = alm_map(hp.read_alm(os.path.join(ACT, 'kappa_alm_data_act_dr6_lensing_v1_cibdeproj.fits')))
    for k, v in P.cross(kcib, 'act', samples=('r', 'g')).items():
        out[f'{k}|cib'] = v
    del kcib
    kc = alm_map(hp.read_alm(os.path.join(ACT, 'kappa_alm_data_act_dr6_lensing_v1_curl.fits')))
    for k, v in P.cross(kc, 'act').items():
        out[f'{k}|curl'] = v
    del kc
    kp = planck_kappa()
    for lens in ('planck', 'pin', 'pout'):
        for k, v in P.cross(kp, lens).items():
            out[f'{k}|{lens}'] = v
    for k in ('r', 'g', 'r_ran', 'g_ran'):
        out[f'{k}|auto'] = P.auto(k, nbar_map(nb[k.split('_')[0]], P.field))
    np.savez(os.path.join(OUT, 'data_cls.npz'), **out)
    for key, v in out.items():
        log(key, ' '.join(f'{x:+.2e}' for x in v))


# ------------------------------------------------------------------------------------------------------------- sims
def full_map(vals):
    import healpy as hp
    m = np.zeros(hp.nside2npix(NSIDE)); m[region_pix()] = vals.astype(np.float64)
    return m


def actsims(args):
    """the stored ACT simulated reconstructions x every galaxy map (and per field); the transfer from the stored inputs."""
    P = Pipe(need=('act',), auto=False)
    S = np.load(os.path.join(OUT, 'actsims_baseline.npy'), mmap_mode='r')
    done = np.load(os.path.join(OUT, 'actsims_baseline_done.npy'))
    store = {}
    for n, i in enumerate(done):
        for k, v in P.cross(full_map(S[i - 1]), 'act', perfield=True).items():
            store.setdefault(k, []).append(v)
        if (n + 1) % 50 == 0:
            log(f'actsims {n + 1}/{len(done)}')
    # transfer: rec_i x in_i(galaxy window) / in_i x in_i(galaxy window)
    I = np.load(os.path.join(OUT, 'actsims_input.npy'), mmap_mode='r')
    idone = np.load(os.path.join(OUT, 'actsims_input_done.npy'))
    nmt = P.nmt
    for i in idone:
        if i not in set(done.tolist()):
            continue
        kin = full_map(I[i - 1]); krec = full_map(S[i - 1])
        fin = nmt.NmtField(P.w, [kin], lmax=LMAXN, n_iter=1)
        ws = P.ws[('act', 'all')]
        store.setdefault('T_rec', []).append(ws.decouple_cell(nmt.compute_coupled_cell(P.kfield(krec, 'act'), fin))[0][:-1])
        store.setdefault('T_in', []).append(ws.decouple_cell(nmt.compute_coupled_cell(P.kfield(kin, 'act'), fin))[0][:-1])
    np.savez(os.path.join(OUT, 'sims_act.npz'), ids=done, **{k: np.array(v) for k, v in store.items()})
    log(f'actsims done: {len(done)} sims, {len(store.get("T_rec", []))} transfer pairs')


def gsims(args):
    """Planck PR3: Gaussian C+N kappa realisations x the galaxy maps (planck, pin, pout windows)."""
    import healpy as hp
    P = Pipe(need=('planck', 'pin', 'pout'), perfield=False, auto=False)
    nl = np.loadtxt(os.path.join(DATA, 'planck', 'COM_Lensing_4096_R3.00', 'MV', 'nlkk.dat'))
    cl = np.zeros(LMAXN + 1); cl[nl[:LMAXN + 1, 0].astype(int)] = nl[:LMAXN + 1, 2]
    rng = np.random.default_rng(1460 + args.seed)
    store = {}
    for i in range(args.n):
        np.random.seed(int(rng.integers(1 << 31)))
        kmap = hp.alm2map(hp.synalm(cl, lmax=LMAXN), NSIDE, lmax=LMAXN)
        for lens in ('planck', 'pin', 'pout'):
            for k, v in P.cross(kmap, lens).items():
                store.setdefault(f'{k}|{lens}', []).append(v)
        if (i + 1) % 50 == 0:
            log(f'gsims {i + 1}')
    np.savez(os.path.join(OUT, f'gsims_planck_{args.seed}.npz'), **{k: np.array(v) for k, v in store.items()})


def fid_spectra(tag, b, A=1.0, f=None, s=None, zint=None, bint=None, sfx=''):
    """unbinned model C^{kg} and C^{gg} (unit pixel window) at the given parameters."""
    t = np.load(os.path.join(OUT, f'theory_{tag}{sfx}.npz'))
    nuis = np.load(os.path.join(OUT, 'nuis.npy'), allow_pickle=True).item()[tag]
    f = nuis['f'] if f is None else f
    s = nuis['s'] if s is None else s
    zint = ZINT[tag][0] if zint is None else zint
    bint = BINT[tag][0] if bint is None else bint
    j = int(np.argmin(np.abs(t['zint'] - zint)))
    a = 5 * s - 2
    ckg = A ** 2 * (1 - f) * b * t['km'] + f * bint * t['kc'][j] + a * t['kmu']
    cgg = (A ** 2 * (1 - f) ** 2 * b ** 2 * t['gg'] + f ** 2 * bint ** 2 * t['cc'][j] + 2 * f * bint * a * t['cmu'][j]
           + 2 * A ** 2 * (1 - f) * b * a * t['gmu'] + a ** 2 * t['mumu'])
    return ckg, cgg, t['kk']


def jsims(args):
    """correlated Gaussian kappa-galaxy realisations with Poisson sampling: the auto and cross-auto covariance blocks."""
    import healpy as hp
    tag = args.tag
    P = Pipe(samples=(tag,), need=('act', 'planck', 'pout'), perfield=False)
    nb = np.load(os.path.join(OUT, 'nbar.npy'), allow_pickle=True).item()[tag]
    nbm = nbar_map(nb, P.field)
    om = hp.nside2pixarea(NSIDE)
    ckg, cgg, ckk = fid_spectra(tag, b={'r': 4.0, 'g': 3.0}[tag])
    L = LMAXN
    nka = np.loadtxt(os.path.join(ACT, 'N_L_kk_act_dr6_lensing_v1_baseline.txt'))
    na = np.zeros(L + 1); na[nka[:L + 1, 0].astype(int)] = nka[:L + 1, 1]
    nl = np.loadtxt(os.path.join(DATA, 'planck', 'COM_Lensing_4096_R3.00', 'MV', 'nlkk.dat'))
    npl = np.zeros(L + 1); npl[nl[:L + 1, 0].astype(int)] = nl[:L + 1, 1]
    cls = [ckk[:L + 1], cgg[:L + 1], ckg[:L + 1]]       # healpy synalm 'new=False' order: TT, EE, TE (kk, gg, kg)
    rng = np.random.default_rng(14600 + args.seed + (0 if tag == 'r' else 777))
    store = {}
    for i in range(args.n):
        np.random.seed(int(rng.integers(1 << 31)))
        ak, ag = hp.synalm([cls[0], cls[1], cls[2]], lmax=L, new=False)
        kt = hp.alm2map(ak, NSIDE, lmax=L)
        dg = hp.alm2map(ag, NSIDE, lmax=L, pixwin=True)
        lam = np.where(P.w > 0, nbm * om * P.w * np.clip(1 + dg, 0, None), 0.0)
        N = np.random.poisson(lam).astype(float)
        d, _ = delta_map(N, P.w, P.field)
        gf = P.nmt.NmtField(P.w, [d], lmax=LMAXN, n_iter=1)
        P.gf[tag] = gf
        ka = kt + hp.alm2map(hp.synalm(na, lmax=L), NSIDE, lmax=L)
        kp = kt + hp.alm2map(hp.synalm(npl, lmax=L), NSIDE, lmax=L)
        store.setdefault('act', []).append(P.cross(ka, 'act')[tag])
        store.setdefault('planck', []).append(P.cross(kp, 'planck')[tag])
        store.setdefault('pout', []).append(P.cross(kp, 'pout')[tag])
        store.setdefault('auto', []).append(P.auto(tag, nbm))
        if (i + 1) % 25 == 0:
            log(f'jsims {tag} {i + 1}')
            np.savez(os.path.join(OUT, f'jsims_{tag}_{args.seed}.npz'), **{k: np.array(v) for k, v in store.items()})
    np.savez(os.path.join(OUT, f'jsims_{tag}_{args.seed}.npz'), **{k: np.array(v) for k, v in store.items()})


# ------------------------------------------------------------------------------------------------------------- fit
PIECES = ('km', 'kmu', 'kk')
APIECES = ('gg', 'gmu', 'mumu')


def binned(tag, lens, sfx=''):
    """binned unit pieces for the cross with `lens` (pixel window once) and the auto (pixel window squared)."""
    import pymaster as nmt, healpy as hp
    t = np.load(os.path.join(OUT, f'theory_{tag}{sfx}.npz'))
    pw = hp.pixwin(NSIDE, lmax=LMAXN)
    def bw(name):
        w = nmt.NmtWorkspace(); w.read_from(os.path.join(OUT, name))
        return w.get_bandpower_windows()[0, :-1, 0, :]
    def b1(B, c, p):
        x = np.zeros(B.shape[1]); n = min(len(c), len(x)); x[:n] = c[:n]
        return B @ (x * p[:B.shape[1]])
    out = {}
    if lens != 'auto':
        B = bw(f'ws_{lens}_all.fits')
        out['km'] = b1(B, t['km'], pw); out['kmu'] = b1(B, t['kmu'], pw)
        out['kc'] = np.array([b1(B, c, pw) for c in t['kc']])
    else:
        B = bw('ws_auto_all.fits')
        for k in APIECES:
            out[k] = b1(B, t[k], pw ** 2)
        out['cc'] = np.array([b1(B, c, pw ** 2) for c in t['cc']])
        out['cmu'] = np.array([b1(B, c, pw ** 2) for c in t['cmu']])
    out['zint'] = t['zint']; out['sigma8_z'] = float(t['sigma8_z']); out['zmean'] = float(t['zmean'])
    return out


MODEL = 'amended'   # 'fixed': (5s-2) on the whole sample (the plan as fixed); 'amended': the LBG-weighted <(1-f)(5s-2)>


def mag_coef(nuis, df, ds):
    """the magnification coefficient multiplying C^{kappa mu} for a shift df of f(m) and ds of s(m)."""
    if MODEL == 'fixed':
        return 5 * (nuis['s'] + ds) - 2
    return nuis['E1'] - df * nuis['E2'] + 5 * ds * nuis['E3'] - 5 * df * ds


def draws(tag, n=300, seed=7):
    """nuisance draws: f, a (the magnification coefficient), z_int, b_int."""
    nuis = np.load(os.path.join(OUT, 'nuis.npy'), allow_pickle=True).item()[tag]
    rng = np.random.default_rng(seed)
    f = np.clip(rng.normal(nuis['f'], FINT_ERR, n), 0, 0.6)
    ds = rng.normal(0, 0.15, n)
    a = mag_coef(nuis, f - nuis['f'], ds)
    zi = np.clip(rng.normal(*ZINT[tag], n), 0.2, 1.6)
    bi = np.clip(rng.normal(*BINT[tag], n), 0.3, None)
    return f, a, zi, bi


def cross_model(T, b, A, f, a, zi, bi, transfer=1.0):
    j = int(np.argmin(np.abs(T['zint'] - zi)))
    return transfer * (A ** 2 * (1 - f) * b * T['km'] + f * bi * T['kc'][j] + a * T['kmu'])


def auto_model(T, b, A, f, a, zi, bi):
    j = int(np.argmin(np.abs(T['zint'] - zi)))
    return (A ** 2 * (1 - f) ** 2 * b ** 2 * T['gg'] + f ** 2 * bi ** 2 * T['cc'][j] + 2 * f * bi * a * T['cmu'][j]
            + 2 * A ** 2 * (1 - f) * b * a * T['gmu'] + a ** 2 * T['mumu'])


def hinv(C, nsim):
    nb = len(C)
    return (nsim - nb - 2) / (nsim - 1) * np.linalg.inv(C)


def post_1d(grid, lnL):
    p = np.exp(lnL - lnL.max()); p /= p.sum()
    c = np.cumsum(p)
    q = [grid[np.searchsorted(c, x)] for x in (0.16, 0.5, 0.84)]
    return q[1], q[1] - q[0], q[2] - q[1], grid[np.argmax(p)]


def fit_cross_b(d, Ci, Ts, nd, transfer=None, bgrid=np.linspace(-10, 30, 801)):
    """b at A = 1 from cross vectors d (concatenated legs), Ts list of binned theory per leg, nuisances marginalised."""
    f, s, zi, bi = nd
    L = np.zeros((len(f), len(bgrid)))
    for i in range(len(f)):
        u = np.concatenate([(transfer[k] if transfer is not None else 1.0) * (1 - f[i]) * T['km'] for k, T in enumerate(Ts)])
        v = np.concatenate([cross_model(T, 0.0, 1.0, f[i], s[i], zi[i], bi[i], transfer[k] if transfer is not None else 1.0)
                            for k, T in enumerate(Ts)])
        r = d[None, :] - v[None, :] - bgrid[:, None] * u[None, :]
        L[i] = -0.5 * np.einsum('ij,jk,ik->i', r, Ci, r)
    lnL = np.log(np.mean(np.exp(L - L.max()), 0)) + L.max()
    return post_1d(bgrid, lnL), bgrid, lnL


def _lnL_cross(d, Ci, Ts, nd, transfer, bgrid=np.linspace(-10, 30, 801)):
    return fit_cross_b(d, Ci, Ts, nd, transfer, bgrid)[2]


def central(tag):
    nuis = np.load(os.path.join(OUT, 'nuis.npy'), allow_pickle=True).item()[tag]
    return nuis['f'], mag_coef(nuis, 0.0, 0.0), ZINT[tag][0], BINT[tag][0]


def snr_cross(d, Ci, Ts, tag, transfer=None):
    """sqrt(chi2_null - chi2_best) with the nuisances at their centres and b free (A = 1); also chi2 and PTE."""
    from scipy.stats import chi2
    f, s, zi, bi = central(tag)
    tr = transfer if transfer is not None else [1.0] * len(Ts)
    u = np.concatenate([tr[k] * (1 - f) * T['km'] for k, T in enumerate(Ts)])
    v = np.concatenate([cross_model(T, 0.0, 1.0, f, s, zi, bi, tr[k]) for k, T in enumerate(Ts)])
    bh = (u @ Ci @ (d - v)) / (u @ Ci @ u)
    r = d - v - bh * u
    c2 = r @ Ci @ r; c0 = d @ Ci @ d
    # amplitude of the best-fit shape against zero (the significance of the cross power)
    m = v + bh * u
    amp = (m @ Ci @ d) / (m @ Ci @ m); samp = 1 / np.sqrt(m @ Ci @ m)
    return dict(b=bh, sb=1 / np.sqrt(u @ Ci @ u), chi2=c2, dof=len(d) - 1, pte=chi2.sf(c2, len(d) - 1), chi0=c0,
                pte0=chi2.sf(c0, len(d)), snr=np.sqrt(max(c0 - c2, 0)), amp=amp, samp=samp)


def fit_joint(dx, Tx, tr, da, Ta, Cfull_inv, nd, bgrid=np.linspace(0, 20, 201), Agrid=np.linspace(0.01, 3, 300)):
    """(b, A) from cross legs + auto, nuisances marginalised; returns marginal posteriors."""
    f, aa, zi, bi = nd
    d = np.concatenate([dx, da])
    Bg, Ag = np.meshgrid(bgrid, Agrid, indexing='ij')
    acc = None; mx = -np.inf
    Ls = []
    for i in range(len(f)):
        jx = [int(np.argmin(np.abs(T['zint'] - zi[i]))) for T in Tx]
        a = aa[i]
        # cross: tr*(A^2 (1-f) b km + f bi kc + a kmu); auto: A^2(1-f)^2 b^2 gg + const + 2 A^2 (1-f) b a gmu
        km = np.concatenate([tr[k] * (1 - f[i]) * T['km'] for k, T in enumerate(Tx)])
        cx = np.concatenate([tr[k] * (f[i] * bi[i] * T['kc'][jx[k]] + a * T['kmu']) for k, T in enumerate(Tx)])
        ja = int(np.argmin(np.abs(Ta['zint'] - zi[i])))
        gg = (1 - f[i]) ** 2 * Ta['gg']; gmu = 2 * (1 - f[i]) * a * Ta['gmu']
        ca = f[i] ** 2 * bi[i] ** 2 * Ta['cc'][ja] + 2 * f[i] * bi[i] * a * Ta['cmu'][ja] + a ** 2 * Ta['mumu']
        X = (Ag ** 2 * Bg)[..., None]
        mod = np.concatenate([X * km + cx, (Ag ** 2 * Bg ** 2)[..., None] * gg + X * gmu + ca], axis=-1)
        r = d - mod
        L = -0.5 * np.einsum('abi,ij,abj->ab', r, Cfull_inv, r)
        Ls.append(L)
    Ls = np.array(Ls); m = Ls.max()
    lnL = np.log(np.mean(np.exp(Ls - m), 0)) + m
    p = np.exp(lnL - lnL.max())
    pb = p.sum(1); pA = p.sum(0)
    return post_1d(bgrid, np.log(pb + 1e-300)), post_1d(Agrid, np.log(pA + 1e-300)), lnL, bgrid, Agrid


def fit(args):
    global MODEL
    from scipy.stats import chi2
    MODEL = args.model
    D = dict(np.load(os.path.join(OUT, 'data_cls.npz')))
    S = dict(np.load(os.path.join(OUT, 'sims_act.npz')))
    G = {}
    for fn in sorted(glob.glob(os.path.join(OUT, 'gsims_planck_*.npz'))):
        for k, v in np.load(fn).items():
            G[k] = np.concatenate([G[k], v]) if k in G else v
    nS = len(S['ids']); nG = len(G['r|planck'])
    lines = [f'# R-146 fit ({time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime())})', '',
             f'ACT simulated reconstructions: {nS}; Planck Gaussian sims: {nG}', '']
    if 'T_rec' in S:
        T = S['T_rec'].mean(0) / S['T_in'].mean(0)
        eT = S['T_rec'].std(0) / np.sqrt(len(S['T_rec'])) / np.abs(S['T_in'].mean(0))
        lines.append(f'ACT transfer ({len(S["T_rec"])} sims): ' + ' '.join(f'{t:.3f}±{e:.3f}' for t, e in zip(T, eT)))
    else:
        T = np.ones(len(EDGES['act']) - 1); lines.append('ACT transfer: not measured (1)')
    lines.append('')
    res = {}
    for tag in ('r', 'g'):
        Ta = binned(tag, 'act'); Tp = binned(tag, 'planck'); To = binned(tag, 'pout'); Ti = binned(tag, 'pin'); Tauto = binned(tag, 'auto')
        nd = draws(tag)
        f0, a0, _, _ = central(tag)
        nu = np.load(os.path.join(OUT, 'nuis.npy'), allow_pickle=True).item()[tag]
        lines.append(f'## {tag}-dropouts, model {MODEL} (<z> = {Ta["zmean"]:.2f}; f_int = {f0:.3f} ± {FINT_ERR}, s = {nu["s"]:.3f} ± 0.15; '
                     f'magnification coefficient {a0:.2f})')
        lines.append('')
        CA = np.cov(S[tag].T); CP = np.cov(G[f'{tag}|planck'].T); CO = np.cov(G[f'{tag}|pout'].T); CI = np.cov(G[f'{tag}|pin'].T)
        legs = {
            'ACT DR6': (D[f'{tag}|act'], hinv(CA, nS), [Ta], [T]),
            'Planck PR3 (all fields)': (D[f'{tag}|planck'], hinv(CP, nG), [Tp], None),
            'Planck PR3 inside ACT': (D[f'{tag}|pin'], hinv(CI, nG), [Ti], None),
            'Planck PR3 outside ACT (W06)': (D[f'{tag}|pout'], hinv(CO, nG), [To], None),
        }
        dC = np.concatenate([D[f'{tag}|act'], D[f'{tag}|pout']])
        CC = np.zeros((len(dC), len(dC))); na = len(CA); CC[:na, :na] = CA; CC[na:, na:] = CO
        legs['ACT ⊕ Planck outside ACT'] = (dC, np.linalg.inv(CC) * 0 + np.block([[hinv(CA, nS), np.zeros((na, len(CO)))], [np.zeros((len(CO), na)), hinv(CO, nG)]]), [Ta, To], [T, np.ones(len(CO))])
        lines.append('| lensing | cross S/N | χ²/dof (PTE) | χ²_null (PTE) | b (A = 1, nuisances marginalised) | b (central nuisances) |')
        lines.append('|---|---|---|---|---|---|')
        for name, (d, Ci, Ts, tr) in legs.items():
            q = snr_cross(d, Ci, Ts, tag, tr)
            (bm, lo, hi, bmax), _, _ = fit_cross_b(d, Ci, Ts, nd, tr)
            res[(tag, name)] = dict(q=q, b=(bm, lo, hi))
            lines.append(f'| {name} | {q["snr"]:.1f} | {q["chi2"]:.1f}/{q["dof"]} ({q["pte"]:.2f}) | {q["chi0"]:.1f} ({q["pte0"]:.1e}) | '
                         f'{bm:.2f} −{lo:.2f} +{hi:.2f} | {q["b"]:.2f} ± {q["sb"]:.2f} |')
        lines.append('')
        # nulls
        cu = D[f'{tag}|curl']; c2 = cu @ hinv(CA, nS) @ cu
        rn = D[f'{tag}_ran|act']; Cr = np.cov(S[f'{tag}_ran'].T); c2r = rn @ hinv(Cr, nS) @ rn
        rp = D[f'{tag}_ran|planck']; Crp = np.cov(G[f'{tag}_ran|planck'].T); c2rp = rp @ hinv(Crp, nG) @ rp
        lines.append(f'- curl × {tag}: χ² {c2:.1f}/{len(cu)} (PTE {chi2.sf(c2, len(cu)):.2f}) (baseline-sim covariance)')
        lines.append(f'- randoms × ACT: χ² {c2r:.1f}/{len(rn)} (PTE {chi2.sf(c2r, len(rn)):.2f}); randoms × Planck: χ² {c2rp:.1f}/{len(rp)} '
                     f'(PTE {chi2.sf(c2rp, len(rp)):.2f})')
        # ACT vs Planck on the same sky
        qa = res[(tag, 'ACT DR6')]['q']; qi = res[(tag, 'Planck PR3 inside ACT')]['q']
        dd = qa['b'] - qi['b']; sd = np.hypot(qa['sb'], qi['sb'])
        lines.append(f'- ACT vs Planck on the ACT sky (central nuisances): {qa["b"]:.2f} ± {qa["sb"]:.2f} vs {qi["b"]:.2f} ± {qi["sb"]:.2f}: '
                     f'Δ = {dd:.2f} ± {sd:.2f} ({dd / sd:.1f}σ)')
        # per field
        bf = []
        for fld in FIELDS:
            k = f'{tag}_W0{fld}'
            Ts_ = binned_field(tag, fld)
            q = snr_cross(D[f'{k}|act'], hinv(np.cov(S[k].T), nS), [Ts_], tag, [T])
            bf.append((q['b'], q['sb']))
        bf = np.array(bf); wmean = np.sum(bf[:, 0] / bf[:, 1] ** 2) / np.sum(1 / bf[:, 1] ** 2)
        c2f = np.sum(((bf[:, 0] - wmean) / bf[:, 1]) ** 2)
        ptef = chi2.sf(c2f, len(bf) - 1)
        lines.append('- per field (ACT, central nuisances): ' + ', '.join(f'W0{f}: {b:.2f} ± {e:.2f}' for f, (b, e) in zip(FIELDS, bf)) +
                     f'; χ² {c2f:.1f}/{len(bf) - 1} (PTE {ptef:.2f})')
        # CIB
        dcib = D[f'{tag}|cib']
        nb_ = np.loadtxt(os.path.join(ACT, 'N_L_kk_act_dr6_lensing_v1_baseline.txt'))[:, 1]
        nc_ = np.loadtxt(os.path.join(ACT, 'N_L_kk_act_dr6_lensing_v1_cibdeproj.txt'))[:, 1]
        e = EDGES['act']; ratio = np.array([np.mean(nc_[e[i]:e[i + 1]] / nb_[e[i]:e[i + 1]]) for i in range(len(e) - 1)])
        CAc = CA * np.sqrt(np.outer(ratio, ratio))
        qc = snr_cross(dcib, hinv(CAc, nS), [Ta], tag, [T])
        sdc = np.sqrt(max(qc['sb'] ** 2 - qa['sb'] ** 2, 0))
        lines.append(f'- CIB-deprojected ACT: b = {qc["b"]:.2f} ± {qc["sb"]:.2f} (baseline {qa["b"]:.2f}); Δ = {qc["b"] - qa["b"]:+.2f}, expected '
                     f'scatter {sdc:.2f}; S/N {qc["snr"]:.1f}')
        cib_extra = max(abs(qc['b'] - qa['b']) - sdc, 0)
        infl = 1.0
        if ptef < 0.05:
            infl = max(infl, np.sqrt(c2f / (len(bf) - 1)))
        if abs(dd / sd) > 2:
            infl = max(infl, abs(dd / sd))
        lines.append(f'- error inflation from the consistency rules: ×{infl:.2f}; CIB systematic added in quadrature to b: {cib_extra:.2f}')
        res[(tag, 'infl')] = infl; res[(tag, 'cib')] = cib_extra
        lines.append('')
        # joint cross + auto
        J = {}
        for fn in sorted(glob.glob(os.path.join(OUT, f'jsims_{tag}_*.npz'))):
            for k, v in np.load(fn).items():
                J[k] = np.concatenate([J[k], v]) if k in J else v
        if J:
            nJ = len(J['auto'])
            X = np.concatenate([J['act'], J['pout'], J['auto']], axis=1)
            R = np.corrcoef(X.T)
            sig = np.concatenate([np.sqrt(np.diag(CA)), np.sqrt(np.diag(CO)), np.std(J['auto'], 0)])
            Cf = R * np.outer(sig, sig)
            Cfi = hinv(Cf, min(nJ, nS, nG))
            da = D[f'{tag}|auto']
            lines.append(f'auto ({tag}, shot noise subtracted): ' + ' '.join(f'{x:.2e}' for x in da) + f'; Gaussian-sim σ: ' +
                         ' '.join(f'{x:.1e}' for x in np.std(J['auto'], 0)) + f' ({nJ} joint sims)')
            Jx = np.cov(J['act'].T)
            lines.append('ACT cross σ, release sims vs joint Gaussian sims: ' + ' '.join(f'{a:.2f}' for a in np.sqrt(np.diag(CA) / np.diag(Jx))))
            ra = da @ hinv(np.cov(J['auto'].T), nJ) @ da
            lines.append(f'auto χ² against zero: {ra:.1f}/{len(da)}')
            nd2 = draws(tag, n=150, seed=11)
            (bm, blo, bhi, _), (Am, Alo, Ahi, Amax), lnL, bg, Ag = fit_joint(
                np.concatenate([D[f'{tag}|act'], D[f'{tag}|pout']]), [Ta, To], [T, np.ones(len(CO))], da, Tauto, Cfi, nd2)
            s8 = Ta['sigma8_z']
            res[(tag, 'joint')] = dict(b=(bm, blo, bhi), A=(Am, Alo, Ahi), s8=s8)
            # best-fit chi2 of the joint
            i, j = np.unravel_index(np.argmax(lnL), lnL.shape)
            lines.append(f'**joint (cross ACT ⊕ Planck W06 + auto): b = {bm:.2f} −{blo:.2f} +{bhi:.2f}; A = {Am:.2f} −{Alo:.2f} +{Ahi:.2f} '
                         f'(peak {Amax:.2f}); σ8(z = {Ta["zmean"]:.2f}) = {Am * s8:.3f} −{Alo * s8:.3f} +{Ahi * s8:.3f} vs Planck ΛCDM '
                         f'{s8:.3f}** (grid peak b {bg[i]:.2f}, A {Ag[j]:.2f})')
            if infl > 1 or cib_extra > 0:
                lines.append(f'with the consistency inflation ×{infl:.2f} and the CIB term: A = {Am:.2f} −{Alo * infl:.2f} +{Ahi * infl:.2f}')
            # the two halves separately (A = 1): the cross measures (1-f) b A^2, the auto (1-f)^2 b^2 A^2; their ratio is A
            (bc, bclo, bchi, _), _, _ = fit_cross_b(np.concatenate([D[f'{tag}|act'], D[f'{tag}|pout']]),
                                                   legs['ACT ⊕ Planck outside ACT'][1], [Ta, To], nd, [T, np.ones(len(CO))])
            Cai = hinv(np.cov(J['auto'].T), nJ)
            bgrid = np.linspace(0, 20, 801)
            La = []
            for i in range(len(nd[0])):
                La.append([-0.5 * (lambda r: r @ Cai @ r)(da - auto_model(Tauto, b_, 1.0, nd[0][i], nd[1][i], nd[2][i], nd[3][i]))
                           for b_ in bgrid])
            La = np.array(La); mL = La.max(); lnLa = np.log(np.mean(np.exp(La - mL), 0)) + mL
            bam, balo, bahi, bamax = post_1d(bgrid, lnLa)
            f0_, a0_, z0_, bi0_ = central(tag)
            ra_ = da - auto_model(Tauto, bamax, 1.0, f0_, a0_, z0_, bi0_)
            c2a = ra_ @ Cai @ ra_
            lines.append(f'cross alone (ACT ⊕ W06, A = 1): b_× = {bc:.2f} −{bclo:.2f} +{bchi:.2f} [= b A²]; auto alone (A = 1): '
                         f'b_auto = {bam:.2f} −{balo:.2f} +{bahi:.2f} [= b A], best-fit χ² {c2a:.1f}/{len(da) - 1} (PTE {chi2.sf(c2a, len(da) - 1):.3f})')
            rng_ = np.random.default_rng(3)
            # ratio A = b_x / b_auto by sampling each half's posterior (independent to the extent the auto is shot-noise/systematics dominated)
            px = np.exp(_lnL_cross(np.concatenate([D[f'{tag}|act'], D[f'{tag}|pout']]), legs['ACT ⊕ Planck outside ACT'][1], [Ta, To], nd,
                                   [T, np.ones(len(CO))]))
            gx = np.linspace(-10, 30, 801); px /= px.sum(); pa_ = np.exp(lnLa - lnLa.max()); pa_ /= pa_.sum()
            xs = rng_.choice(gx, 20000, p=px); as_ = rng_.choice(bgrid, 20000, p=pa_)
            Ar = xs / np.clip(as_, 1e-3, None)
            q16, q50, q84 = np.percentile(Ar, [16, 50, 84])
            lines.append(f'ratio A = b_× / b_auto: {q50:.2f} −{q50 - q16:.2f} +{q84 - q50:.2f} (A = 1 is Planck ΛCDM; P(A > 1) = {np.mean(Ar > 1):.3f})')
            res[(tag, 'ratio')] = (q50, q50 - q16, q84 - q50, float(np.mean(Ar > 1)))
            res[(tag, 'halves')] = dict(bx=(bc, bclo, bchi), bauto=(bam, balo, bahi), chi2auto=c2a)
            lines.append('')
        np.save(os.path.join(OUT, f'fit_res_{MODEL}.npy'), res, allow_pickle=True)
    open(os.path.join(OUT, f'fit_{MODEL}.md'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def binned_field(tag, fld):
    import pymaster as nmt, healpy as hp
    t = np.load(os.path.join(OUT, f'theory_{tag}.npz'))
    pw = hp.pixwin(NSIDE, lmax=LMAXN)
    w = nmt.NmtWorkspace(); w.read_from(os.path.join(OUT, f'ws_act_{fld}.fits'))
    B = w.get_bandpower_windows()[0, :-1, 0, :]
    def b1(c):
        x = np.zeros(B.shape[1]); n = min(len(c), len(x)); x[:n] = c[:n]
        return B @ (x * pw[:B.shape[1]])
    return dict(km=b1(t['km']), kmu=b1(t['kmu']), kc=np.array([b1(c) for c in t['kc']]), zint=t['zint'])


def pvariant(args):
    """post hoc (after the fit): Planck PR3 with ACT's band edges to L = 1000 on all fields, for the z ~ 4 reproduction
    (Miyatake's stack reaches theta ~ 6-20', i.e. L up to ~1000+); data and Gaussian C+N sims, g and r."""
    import healpy as hp, pymaster as nmt
    g = np.load(os.path.join(OUT, 'maps.npz'))
    w = g['w']
    wk = kwindows()['planck']
    fn = os.path.join(OUT, 'ws_p1000_all.fits')
    ws = nmt.NmtWorkspace()
    if os.path.exists(fn):
        ws.read_from(fn)
    else:
        ws.compute_coupling_matrix(nmt.NmtField(wk, None, lmax=LMAXN, spin=0), nmt.NmtField(w, None, lmax=LMAXN, spin=0), bins('act'))
        ws.write_to(fn)
    gf = {k: nmt.NmtField(w, [g[k]], lmax=LMAXN, n_iter=1) for k in ('r', 'g')}
    def cross(kmap):
        kf = nmt.NmtField(wk, [kmap], lmax=LMAXN, n_iter=1)
        return {k: ws.decouple_cell(nmt.compute_coupled_cell(kf, f))[0][:-1] for k, f in gf.items()}
    out = {f'{k}|p1000': v for k, v in cross(planck_kappa()).items()}
    nl = np.loadtxt(os.path.join(DATA, 'planck', 'COM_Lensing_4096_R3.00', 'MV', 'nlkk.dat'))
    cl = np.zeros(LMAXN + 1); cl[nl[:LMAXN + 1, 0].astype(int)] = nl[:LMAXN + 1, 2]
    rng = np.random.default_rng(1469)
    st = {}
    for i in range(args.n):
        np.random.seed(int(rng.integers(1 << 31)))
        for k, v in cross(hp.alm2map(hp.synalm(cl, lmax=LMAXN), NSIDE, lmax=LMAXN)).items():
            st.setdefault(k, []).append(v)
        if (i + 1) % 50 == 0:
            log(f'pvariant sims {i + 1}')
    np.savez(os.path.join(OUT, 'pvariant.npz'), **out, **{f'sim_{k}': np.array(v) for k, v in st.items()})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd')
    ap.add_argument('--n', type=int, default=400)
    ap.add_argument('--map', default='baseline')
    ap.add_argument('--tag', default='r')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--model', default='amended')
    args = ap.parse_args()
    if args.cmd == 'region':
        region()
    else:
        globals()[args.cmd](args)


if __name__ == '__main__':
    main()
