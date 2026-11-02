"""Quasars at 3.5 < z < 5 crossed with CMB lensing: the linear bias with magnification (R-139, cosmo-deep-15).

The plan (samples, masks, multipoles, estimator, model, fit, nulls) is fixed in `results-R139.md` before the first cross-spectrum.
Pseudo-C_ell with NaMaster's mode-coupling matrix; ACT DR6 baseline (and curl) kappa, its simulated reconstructions for the
covariance and the input kappa for the transfer function; Planck PR3 MV kappa with Gaussian C+N realisations for its covariance.

usage (python -I; all data under data/, gitignored; outputs and logs in data/qso_z4/):
  python -I qso_z4_lens.py prep                  # quasar maps, masks, s from the counts, workspaces, binned theory
  python -I qso_z4_lens.py data                  # data cross-spectra (ACT baseline, ACT curl, Planck PR3)
  python -I qso_z4_lens.py sims --map baseline   # consumer of act_dr6_simstream.sh (+ the transfer from inputs)
  python -I qso_z4_lens.py gsims [--n 400]       # Planck PR3 Gaussian C+N realisations
  python -I qso_z4_lens.py fit                   # b, intervals, significance, likelihood ratios, nulls -> fit.md
"""
import argparse, glob, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
DATA = os.path.join(HERE, 'data')
ACT = os.path.join(DATA, 'act_dr6')
OUT = os.path.join(DATA, 'qso_z4')
NSIDE = 512
LMAXK = 1500
EDGES = {'act': [40, 70, 100, 150, 200, 300, 400, 550, 700, 850, 1000], 'planck': [8, 40, 70, 100, 150, 200, 300, 400]}
EDGES['pout'] = EDGES['planck']
EDGES['pin'] = EDGES['planck']
LENS = {'act600': 'ACT DR6 (L < 600)', 'act': 'ACT DR6', 'planck': 'Planck PR3', 'pout': 'Planck PR3 outside ACT', 'pin': 'Planck PR3 inside ACT'}
# samples: name -> (survey, zmin, zmax, randoms-null?)
SAMPLES = {'desi': ('desi', 3.5, 5.0), 'dr16': ('dr16', 3.5, 5.0), 'desi_mid': ('desi', 2.2, 3.5),
           'desi_ran': ('desi', 3.5, 5.0), 'dr16_ran': ('dr16', 3.5, 5.0), 'desi_ebv': ('desi', 3.5, 5.0),
           'dr16_ebv': ('dr16', 3.5, 5.0)}
NTILE_COMP = np.array([np.nan, 0.644, 0.687, 0.677, 0.668, 0.663, 0.659, 0.654])   # fixed in results-R139.md before any fit


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def gal_lat(nside):
    import healpy as hp
    th, ph = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    thg, _ = hp.Rotator(coord=['C', 'G'])(th, ph)
    return 90 - np.degrees(thg)


# ------------------------------------------------------------------------------------------------------------- maps
def desi_maps(zmin, zmax, randoms_null=False, seed=1):
    """(delta, weight) at nside 512: weight = random coverage fraction (NTILE completeness folded in), delta = N_w/(nbar f)-1
    per photometric region (N: DEC > 32.375 in the NGC; S: the rest)."""
    import healpy as hp, fitsio
    npix = hp.nside2npix(NSIDE)
    fn = os.path.join(OUT, 'desi_randoms_ns512.npz')
    if os.path.exists(fn):
        g = np.load(fn); R = g['R']; Rw = g['Rw']; reg = g['reg']; rra = g['rra']; rdec = g['rdec']
    else:
        R = np.zeros(npix); Rw = np.zeros(npix); keep_ra, keep_dec = [], []
        for cap in ('NGC', 'SGC'):
            r = fitsio.read(os.path.join(DATA, 'qso', f'QSO_{cap}_0_clustering.ran.fits'), columns=['RA', 'DEC', 'NTILE'])
            p = hp.ang2pix(NSIDE, r['RA'], r['DEC'], lonlat=True)
            nt = np.clip(r['NTILE'], 1, 7)
            R += np.bincount(p, minlength=npix); Rw += np.bincount(p, weights=NTILE_COMP[nt], minlength=npix)
            sub = np.random.default_rng(7).random(len(r)) < 0.02
            keep_ra.append(r['RA'][sub]); keep_dec.append(r['DEC'][sub])
        rra = np.concatenate(keep_ra); rdec = np.concatenate(keep_dec)
        th, ph = hp.pix2ang(NSIDE, np.arange(npix), lonlat=False)
        dec = 90 - np.degrees(th)
        reg = ((dec > 32.375) & (gal_lat(NSIDE) > 0)).astype(int)
        np.savez(fn, R=R, Rw=Rw, reg=reg, rra=rra, rdec=rdec)
    full = np.median(R[R > 0][R[R > 0] > np.percentile(R[R > 0], 50)])     # the interior random count per pixel
    f = Rw / (full * np.nanmean(NTILE_COMP[1:]))
    cover = R / full
    w = np.where((cover > 0.5) & (np.abs(gal_lat(NSIDE)) > 20), f, 0.0)
    if randoms_null:
        rng = np.random.default_rng(seed)
        d = np.load(os.path.join(DATA, 'qso', 'desi_full_slim.npy'))
        z = d['Z_not4clus']; nsel = int(((z > zmin) & (z < zmax)).sum())
        idx = rng.choice(len(rra), nsel, replace=False)
        p = hp.ang2pix(NSIDE, rra[idx], rdec[idx], lonlat=True); wt = np.ones(nsel)
    else:
        d = np.load(os.path.join(DATA, 'qso', 'desi_full_slim.npy'))
        z = d['Z_not4clus']; s = np.isfinite(z) & (z > zmin) & (z < zmax)
        p = hp.ang2pix(NSIDE, d['RA'][s], d['DEC'][s], lonlat=True)
        wt = d['WEIGHT_ZFAIL'][s] / np.clip(d['FRACZ_TILELOCID'][s], 0.05, None)
    N = np.bincount(p, weights=wt, minlength=npix)
    return to_delta(N, w, reg), w, (N[w > 0].sum())


def dr16_maps(zmin, zmax, randoms_null=False, seed=2):
    import healpy as hp, fitsio
    npix = hp.nside2npix(NSIDE)
    d = fitsio.read(os.path.join(DATA, 'qso', 'DR16Q_v4.fits'), columns=['RA', 'DEC', 'Z'])
    cnt = np.bincount(hp.ang2pix(64, d['RA'], d['DEC'], lonlat=True), minlength=hp.nside2npix(64))
    foot64 = (cnt > 0.3 * cnt[cnt > 0].mean()).astype(float)
    foot = hp.ud_grade(foot64, NSIDE) > 0.5
    w = np.where(foot & (np.abs(gal_lat(NSIDE)) > 20), 1.0, 0.0)
    reg = (gal_lat(NSIDE) > 0).astype(int)
    s = (d['Z'] > zmin) & (d['Z'] < zmax)
    if randoms_null:
        rng = np.random.default_rng(seed)
        ok = np.where(w > 0)[0]
        pr = rng.choice(ok, int(s.sum()), replace=True)
        N = np.bincount(pr, minlength=npix).astype(float)
    else:
        N = np.bincount(hp.ang2pix(NSIDE, d['RA'][s], d['DEC'][s], lonlat=True), minlength=npix).astype(float)
    return to_delta(N, w, reg), w, N[w > 0].sum()


def to_delta(N, w, reg):
    delta = np.zeros_like(w)
    for r in np.unique(reg):
        s = (reg == r) & (w > 0)
        if not s.any():
            continue
        nbar = N[s].sum() / w[s].sum()
        delta[s] = N[s] / (nbar * w[s]) - 1
    return delta


def ebv_cut(w, maxebv=0.15):
    import healpy as hp
    e = np.load(os.path.join(OUT, 'ebv_ns64.npy'))
    e512 = hp.ud_grade(e, NSIDE)
    return np.where(e512 < maxebv, w, 0.0)


# ------------------------------------------------------------------------------------------------------------- kappa
def act_window(variant='baseline'):
    import healpy as hp
    fn = os.path.join(OUT, f'actmask2_{variant}_ns512.npy')
    if os.path.exists(fn):
        return np.load(fn)
    src = os.path.join(ACT, 'mask_act_dr6_lensing_v1_healpix_nside_4096_baseline.fits' if variant == 'baseline' else 'mask_curl.fits')
    m = hp.read_map(src, dtype=np.float32)
    m512 = hp.ud_grade(m ** 2, NSIDE).astype(np.float64)
    np.save(fn, m512)
    return m512


def planck_window():
    import kappa_templates as kt
    return kt.planck_mask512().astype(np.float64)


def truncate(alm, lmax):
    import kappa_templates as kt
    return kt.truncate_alm(alm, lmax)


def alm_map(alm, lmax=LMAXK):
    import healpy as hp
    a = truncate(np.nan_to_num(alm.astype(np.complex128)), lmax)   # the ACT data alm hold NaN at L = 0, 1
    return hp.alm2map(a, NSIDE, lmax=lmax)


def read_act(path):
    import healpy as hp
    return alm_map(hp.read_alm(path))


def planck_kappa():
    import kappa_templates as kt
    return alm_map(kt.planck_alm_eq(), lmax=hp_lmax_planck())


def hp_lmax_planck():
    import kappa_templates as kt
    return kt.LMAX


# ------------------------------------------------------------------------------------------------------------- NaMaster
def bins(lens):
    import pymaster as nmt
    e = EDGES[lens] + [LMAXK + 1]          # the last (overflow) bin is dropped from every output
    return nmt.NmtBin.from_edges(np.array(e[:-1]), np.array(e[1:]))


class Pipe:
    """the quasar fields and the workspaces for each (sample, lens)."""
    def __init__(self, need=('act', 'planck', 'pout', 'pin')):
        import pymaster as nmt
        self.nmt = nmt
        g = np.load(os.path.join(OUT, 'maps.npz'))
        self.samples = [k for k in SAMPLES if f'{k}_d' in g]
        self.gf = {k: nmt.NmtField(g[f'{k}_w'], [g[f'{k}_d']], lmax=LMAXK, n_iter=1) for k in self.samples}
        self.wk = {}
        if 'act' in need:
            self.wk['act'] = act_window('baseline')
        if 'planck' in need:
            self.wk['planck'] = planck_window()
        if 'pout' in need:   # Planck on the sky ACT does not cover: independent of the ACT leg, so the two combine
            self.wk['pout'] = planck_window() * (act_window('baseline') <= 0)
        if 'pin' in need:    # Planck on the ACT footprint: the ACT-Planck consistency test on the same sky
            self.wk['pin'] = planck_window() * (act_window('baseline') > 0)
        self.ws = {}
        for lens in self.wk:
            for k in self.samples:
                fn = os.path.join(OUT, f'ws_{k}_{lens}.fits')
                w = nmt.NmtWorkspace()
                if os.path.exists(fn):
                    w.read_from(fn)
                else:
                    kf = nmt.NmtField(self.wk[lens], None, lmax=LMAXK, spin=0)
                    w.compute_coupling_matrix(kf, self.gf[k], bins(lens))
                    w.write_to(fn)
                self.ws[(k, lens)] = w

    def kfield(self, kmap, lens):
        return self.nmt.NmtField(self.wk[lens], [kmap], lmax=LMAXK, n_iter=1)

    def cross(self, kmap, lens, samples=None):
        kf = self.kfield(kmap, lens)
        out = {}
        for k in (samples or self.samples):
            cl = self.nmt.compute_coupled_cell(kf, self.gf[k])
            out[k] = self.ws[(k, lens)].decouple_cell(cl)[0][:-1]
        return out


# ------------------------------------------------------------------------------------------------------------- prep
def intrinsic_slope(r_ref, lo=18.0, hi=20.75, s_faint=0.25):
    """s(r) = d log10 n/dm of the intrinsic differential counts: a quadratic fit of log10 n(r) in 0.25-mag bins over [lo, hi]
    (the reference sample's range taken as complete) for r <= hi, and the faint-end value s_faint beyond (amended plan,
    results-R139.md: the selection is not a hard flux cut, so s is the sample-weighted intrinsic count slope)."""
    e = np.arange(lo, hi + 1e-6, 0.25); h, _ = np.histogram(r_ref, e); c = 0.5 * (e[1:] + e[:-1])
    ok = h > 0
    p = np.polyfit(c[ok], np.log10(h[ok]), 2, w=np.sqrt(h[ok]))
    def s_of(r):
        r = np.asarray(r)
        return np.where(r <= hi, np.polyval(np.polyder(p), np.clip(r, lo, hi)), s_faint)
    return s_of


def sample_slope(r, r_ref):
    """s_eff = <s(r_i)> over the sample; error from the faint-end value (0.12-0.40: z~4 faint-end LF slopes alpha = -1.3 to
    -2.0), the bright-side fit range (hi 20.5 / 21.0, lo 17.5 / 18.5), and a 0.05 floor for the completeness assumption."""
    base = sample = np.mean(intrinsic_slope(r_ref)(r))
    var = []
    for kw in [dict(s_faint=0.12), dict(s_faint=0.40), dict(hi=20.5), dict(hi=21.0), dict(lo=17.5), dict(lo=18.5)]:
        var.append(np.mean(intrinsic_slope(r_ref, **kw)(r)) - base)
    var = np.array(var)
    err = np.sqrt((0.5 * (var[1] - var[0])) ** 2 + np.max(np.abs(var[2:4])) ** 2 + np.max(np.abs(var[4:])) ** 2 + 0.05 ** 2)
    return base, err, var


def prep(args):
    import healpy as hp, fitsio
    import qso_z4_theory as T
    os.makedirs(OUT, exist_ok=True)
    # E(B-V) map at nside 64 from the DESI targets' SFD column (median per pixel; pixels without targets from DR16Q EXTINCTION_r/2.165)
    d = np.load(os.path.join(DATA, 'qso', 'desi_full_slim.npy'))
    d16 = fitsio.read(os.path.join(DATA, 'qso', 'DR16Q_v4.fits'), columns=['RA', 'DEC', 'Z', 'PSFMAG', 'EXTINCTION'])
    if not os.path.exists(os.path.join(OUT, 'ebv_ns64.npy')):
        full = fitsio.read(os.path.join(DATA, 'qso', 'QSO_full_HPmapcut.dat.fits'), columns=['EBV'])
        e = np.full(hp.nside2npix(64), np.nan)
        p = hp.ang2pix(64, d['RA'], d['DEC'], lonlat=True)
        o = np.argsort(p); ps = p[o]; ev = full['EBV'][o]
        u, st = np.unique(ps, return_index=True)
        for i, (pp, a) in enumerate(zip(u, st)):
            b_ = st[i + 1] if i + 1 < len(st) else len(ps)
            e[pp] = np.median(ev[a:b_])
        p2 = hp.ang2pix(64, d16['RA'], d16['DEC'], lonlat=True)
        e16 = np.bincount(p2, weights=d16['EXTINCTION'][:, 2] / 2.165, minlength=len(e)) / np.maximum(np.bincount(p2, minlength=len(e)), 1)
        e = np.where(np.isfinite(e), e, e16)
        np.save(os.path.join(OUT, 'ebv_ns64.npy'), e)
    maps = {}
    for k, (surv, z0, z1) in SAMPLES.items():
        fnc = desi_maps if surv == 'desi' else dr16_maps
        delta, w, n = fnc(z0, z1, randoms_null=k.endswith('_ran'))
        if k.endswith('_ebv'):
            w = ebv_cut(w)
            delta = np.where(w > 0, delta, 0.0)
        maps[f'{k}_d'] = delta; maps[f'{k}_w'] = w
        log(f'{k}: N = {n:.0f} (weighted), area {np.sum(w > 0) * hp.nside2pixarea(NSIDE, degrees=True):.0f} deg2, '
            f'fsky_eff {np.mean(w):.4f}')
    np.savez(os.path.join(OUT, 'maps.npz'), **maps)
    # magnification slopes from the samples' own counts
    zd = d['Z_not4clus']; sd = np.isfinite(zd) & (zd > 3.5) & (zd < 5)
    rd = 22.5 - 2.5 * np.log10(np.clip(d['FLUX_R'][sd] / d['MW_TRANSMISSION_R'][sd], 1e-4, None))
    s16 = (d16['Z'] > 3.5) & (d16['Z'] < 5)
    r16 = d16['PSFMAG'][s16, 2] - d16['EXTINCTION'][s16, 2]
    zm = (zd > 2.2) & (zd < 3.5) & np.isfinite(zd)
    rm = 22.5 - 2.5 * np.log10(np.clip(d['FLUX_R'][zm] / d['MW_TRANSMISSION_R'][zm], 1e-4, None))
    slopes = {}
    for name, m, ref in [('desi', rd, rd), ('dr16', r16, rd), ('desi_mid', rm, rm)]:
        s_, e_, var = sample_slope(m, ref)
        slopes[name] = (s_, e_)
        log(f's[{name}] = {s_:.3f} ± {e_:.3f} (5s-2 = {5 * s_ - 2:.2f} ± {5 * e_:.2f}); shifts (faint 0.12, 0.40, hi 20.5, 21.0, '
            f'lo 17.5, 18.5): ' + ' '.join(f'{v:+.3f}' for v in var) + f'; fraction beyond r 20.75: {np.mean(m > 20.75):.2f}')
    # theory per sample (n(z) from the sample)
    th = T.Theory(lmax=LMAXK)
    spec = {}
    for name, z, wz in [('desi', zd[sd], None), ('dr16', d16['Z'][s16], None), ('desi_mid', zd[zm], None)]:
        h, e_ = np.histogram(z, bins=np.arange(z.min() - 0.025, z.max() + 0.05, 0.05))
        zc = 0.5 * (e_[1:] + e_[:-1])
        sp = th.spectra(zc, h.astype(float))
        spec[name] = sp
        np.savez(os.path.join(OUT, f'theory_{name}.npz'), zc=zc, nz=h, **sp)
        log(f'theory {name}: mean z {z.mean():.3f}; C^km(100) {sp["km"][100]:.3e}, C^kmu(100) {sp["kmu"][100]:.3e}')
    np.save(os.path.join(OUT, 'slopes.npy'), slopes, allow_pickle=True)
    Pipe()   # builds the workspaces
    log('prep done')


# ------------------------------------------------------------------------------------------------------------- data
def data(args):
    P = Pipe()
    out = {}
    for lens, kmap in [('act', read_act(os.path.join(ACT, 'kappa_alm_data_act_dr6_lensing_v1_baseline.fits'))),
                       ('planck', planck_kappa())]:
        for k, v in P.cross(kmap, lens).items():
            out[f'{k}|{lens}'] = v
        if lens == 'planck':
            for lx in ('pout', 'pin'):
                for k, v in P.cross(kmap, lx).items():
                    out[f'{k}|{lx}'] = v
    # ACT curl on the baseline window (the curl mask is the baseline one; checked in prep log)
    kc = read_act(os.path.join(ACT, 'kappa_alm_data_act_dr6_lensing_v1_curl.fits'))
    for k, v in P.cross(kc, 'act').items():
        out[f'{k}|curl'] = v
    np.savez(os.path.join(OUT, 'data_cls.npz'), **out)
    for key, v in out.items():
        log(key, ' '.join(f'{x:+.2e}' for x in v))


# ------------------------------------------------------------------------------------------------------------- sims
def sims(args):
    """consume streamed simulated reconstructions (and inputs 1..N for the transfer): kappa_sim x every quasar map."""
    import healpy as hp
    P = Pipe(need=('act',))
    variant = args.map
    fn = os.path.join(OUT, f'sims_{variant}.npz')
    res = dict(np.load(fn)) if os.path.exists(fn) else {}
    done = set(int(x) for x in res.get('ids', []))
    store = {k: list(res[k]) for k in res if k != 'ids'}
    simdir = os.path.join(ACT, 'sims'); indir = os.path.join(ACT, 'inputs')
    idle = 0
    # transfer: window for the input in a quasar-like geometry: the DESI weight
    g = np.load(os.path.join(OUT, 'maps.npz'))
    import pymaster as nmt
    wg = g['desi_w']
    b = bins('act')
    tw = None
    while len(done) < args.nsims and idle < 600:
        files = sorted(glob.glob(os.path.join(simdir, f'kappa_alm_sim_act_dr6_lensing_v1_{variant}_*.fits')))
        if not files:
            if os.path.exists(os.path.join(ACT, f'simstream_{variant}.log')) and 'ALLDONE' in open(os.path.join(ACT, f'simstream_{variant}.log')).read():
                break
            time.sleep(2); idle += 1; continue
        idle = 0
        f = files[0]
        i = int(f[-9:-5])
        try:
            alm = hp.read_alm(f)
        except Exception as e:
            log(f'sim {i}: unreadable ({e}); waiting'); time.sleep(5)
            try:
                alm = hp.read_alm(f)
            except Exception:
                os.remove(f); continue
        kmap = alm_map(alm)
        cl = P.cross(kmap, 'act')
        for k, v in cl.items():
            store.setdefault(k, []).append(v)
        if variant == 'baseline':
            inf = os.path.join(indir, f'input_kappa_alm_sim_{i:04d}.fits')
            if i <= args.ntransfer:
                t0 = time.time()
                while not os.path.exists(inf) and time.time() - t0 < 600:
                    time.sleep(3)
            if os.path.exists(inf):
                time.sleep(1)
                kin = alm_map(hp.read_alm(inf))
                if tw is None:
                    fa = nmt.NmtField(P.wk['act'], None, lmax=LMAXK, spin=0); fg = nmt.NmtField(wg, None, lmax=LMAXK, spin=0)
                    tw = nmt.NmtWorkspace(); tw.compute_coupling_matrix(fa, fg, b)
                fin_g = nmt.NmtField(wg, [kin], lmax=LMAXK, n_iter=1)
                rec = tw.decouple_cell(nmt.compute_coupled_cell(P.kfield(kmap, 'act'), fin_g))[0][:-1]
                inn = tw.decouple_cell(nmt.compute_coupled_cell(P.kfield(kin, 'act'), fin_g))[0][:-1]
                store.setdefault('T_rec', []).append(rec); store.setdefault('T_in', []).append(inn)
                os.remove(inf); open(os.path.join(indir, f'done_{i:04d}'), 'w').close()
        os.remove(f); open(os.path.join(simdir, f'done_{variant}_{i:04d}'), 'w').close()
        done.add(i)
        if len(done) % 10 == 0:
            np.savez(fn, ids=np.array(sorted(done)), **{k: np.array(v) for k, v in store.items()})
            log(f'{variant}: {len(done)} sims')
    np.savez(fn, ids=np.array(sorted(done)), **{k: np.array(v) for k, v in store.items()})
    log(f'{variant}: done, {len(done)} sims')


def gsims(args):
    import healpy as hp
    P = Pipe(need=('planck', 'pout', 'pin'))
    nl = np.loadtxt(os.path.join(DATA, 'planck', 'COM_Lensing_4096_R3.00', 'MV', 'nlkk.dat'))
    lmax = hp_lmax_planck()
    cl = np.zeros(lmax + 1); cl[nl[:lmax + 1, 0].astype(int)] = nl[:lmax + 1, 2]
    rng = np.random.default_rng(139)
    store = {}
    for i in range(args.n):
        np.random.seed(int(rng.integers(1 << 31)))
        a = hp.synalm(cl, lmax=lmax)
        kmap = hp.alm2map(a, NSIDE, lmax=lmax)
        for lens in ('planck', 'pout', 'pin'):
            for k, v in P.cross(kmap, lens).items():
                store.setdefault(k if lens == 'planck' else f'{k}|{lens}', []).append(v)
        if (i + 1) % 50 == 0:
            log(f'gsims {i + 1}')
    np.savez(os.path.join(OUT, 'gsims_planck.npz'), **{k: np.array(v) for k, v in store.items()})


# ------------------------------------------------------------------------------------------------------------- fit
def binned_theory(name, lens, sample_key, ws=None):
    import pymaster as nmt
    t = np.load(os.path.join(OUT, f'theory_{name}.npz'))
    w = ws
    if w is None:
        w = nmt.NmtWorkspace(); w.read_from(os.path.join(OUT, f'ws_{sample_key}_{lens}.fits'))
    bpw = w.get_bandpower_windows()[0, :-1, 0, :]
    import healpy as hp
    pw = hp.pixwin(NSIDE, lmax=bpw.shape[1] - 1)
    out = []
    for c in ('km', 'kmu'):
        cl = np.zeros(bpw.shape[1]); n = min(len(t[c]), len(cl)); cl[:n] = t[c][:n]
        out.append(bpw @ (cl * pw))
    return np.array(out)    # (2, nbin)


def fit_one(d, C, Tm, s, es, nsim, bgrid=np.linspace(-10, 40, 5001)):
    """b free, (5s-2) Gaussian prior marginalised exactly. Returns dict."""
    nb = len(d)
    hart = (nsim - nb - 2) / (nsim - 1)
    a = 5 * s - 2; ea = 5 * es
    Ce = C + ea ** 2 * np.outer(Tm[1], Tm[1])
    Ci = hart * np.linalg.inv(Ce)
    r0 = d - a * Tm[1]
    chi = np.array([(r0 - b * Tm[0]) @ Ci @ (r0 - b * Tm[0]) for b in bgrid])
    F = Tm[0] @ Ci @ Tm[0]
    bhat = (Tm[0] @ Ci @ r0) / F
    sig = 1 / np.sqrt(F)
    chi_min = (r0 - bhat * Tm[0]) @ Ci @ (r0 - bhat * Tm[0])
    chi_null = d @ (hart * np.linalg.inv(C)) @ d   # no signal at all (b = 0 and no magnification)
    def lnL(b):
        return -0.5 * (r0 - b * Tm[0]) @ Ci @ (r0 - b * Tm[0])
    # amplitude significance with the shape fixed at the best fit: sqrt(chi_null - chi(best model)) against C alone
    model = bhat * Tm[0] + a * Tm[1]
    Cn = hart * np.linalg.inv(C)
    A = (model @ Cn @ d) / (model @ Cn @ model); sA = 1 / np.sqrt(model @ Cn @ model)
    from scipy.stats import chi2
    return dict(b=bhat, sb=sig, chi2=chi_min, dof=nb - 1, pte=chi2.sf(chi_min, nb - 1), chi_null=chi_null,
                pte_null=chi2.sf(chi_null, nb), snr=A / sA, lr13_4=np.exp(lnL(13) - lnL(4)), a=a, ea=ea)


def fit(args):
    import healpy as hp
    slopes = np.load(os.path.join(OUT, 'slopes.npy'), allow_pickle=True).item()
    D = np.load(os.path.join(OUT, 'data_cls.npz'))
    S = dict(np.load(os.path.join(OUT, 'sims_baseline.npz')))
    G = dict(np.load(os.path.join(OUT, 'gsims_planck.npz')))
    Sc = dict(np.load(os.path.join(OUT, 'sims_curl.npz'))) if os.path.exists(os.path.join(OUT, 'sims_curl.npz')) else None
    lines = [f'# R-139 fit ({time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime())})', '',
             f'ACT sims: {len(S["ids"])}; Planck Gaussian sims: {len(G["desi"])}', '']
    # transfer
    if 'T_rec' in S:
        T = S['T_rec'].mean(0) / S['T_in'].mean(0)
        eT = S['T_rec'].std(0) / np.sqrt(len(S['T_rec'])) / np.abs(S['T_in'].mean(0))
        lines.append(f'ACT transfer (rec x in / in x in, {len(S["T_rec"])} sims): ' + ' '.join(f'{t:.3f}±{e:.3f}' for t, e in zip(T, eT)))
    else:
        T = np.ones(len(EDGES['act']) - 1)
    lines.append('')
    s_of = {'desi': 'desi', 'dr16': 'dr16', 'desi_mid': 'desi_mid', 'desi_ebv': 'desi', 'dr16_ebv': 'dr16', 'desi_ran': 'desi', 'dr16_ran': 'dr16'}
    th_of = {'desi': 'desi', 'dr16': 'dr16', 'desi_mid': 'desi_mid', 'desi_ebv': 'desi', 'dr16_ebv': 'dr16', 'desi_ran': 'desi', 'dr16_ran': 'dr16'}
    results = {}
    lines.append('| sample | map | L range | b | 68 % | χ²/dof (PTE) | χ²_null (PTE) | S/N | L(13)/L(4) | 5s−2 |')
    lines.append('|---|---|---|---|---|---|---|---|---|---|')
    def simof(k, lens):
        return S[k] if lens == 'act' else (G[k] if lens == 'planck' else G[f'{k}|{lens}'])
    for lens in ('act', 'planck', 'pout', 'pin'):
        for k in ['desi', 'dr16', 'desi_mid', 'desi_ebv', 'dr16_ebv', 'desi_ran', 'dr16_ran']:
            for lcut in ([None, 600] if lens == 'act' else [None]):
                d = D[f'{k}|{lens}'].copy()
                sim = simof(k, lens).copy()
                Tm = binned_theory(th_of[k], lens, k)
                if lens == 'act':
                    d = d / T; sim = sim / T
                nb = len(d)
                if lcut:
                    nb = int(np.sum(np.array(EDGES['act'][1:]) <= lcut))
                d = d[:nb]; sim = sim[:, :nb]; Tm = Tm[:, :nb]
                C = np.cov(sim.T)
                s, es = slopes[s_of[k]]
                r = fit_one(d, C, Tm, s, es, len(sim))
                results[(k, lens, lcut)] = (r, d, C, Tm)
                rng_ = f'{EDGES[lens][0]}-{lcut or EDGES[lens][-1]}'
                lines.append(f'| {k} | {LENS[lens]} | {rng_} | {r["b"]:.2f} | ±{r["sb"]:.2f} | {r["chi2"]:.1f}/{r["dof"]} ({r["pte"]:.2f}) | '
                             f'{r["chi_null"]:.1f} ({r["pte_null"]:.1e}) | {r["snr"]:.1f} | {r["lr13_4"]:.2e} | {r["a"]:.2f}±{r["ea"]:.2f} |')
    # combined: joint DESI + DR16Q on both maps (ACT and Planck sims independent -> block-diagonal between maps; the
    # overlap of the two maps' footprints is shared sky, handled in the ACT-Planck consistency and stated as a caveat)
    lines.append('')
    lines.append('## Combined')
    for combo in [('desi',), ('dr16',), ('desi_mid',), ('desi', 'dr16')]:
        for maps_ in [('act',), ('planck',), ('pout',), ('act', 'pout'), ('act600', 'pout')]:
            ds, Ts, blocks = [], [], []
            for lens0 in maps_:
                lens, lcut = ('act', 600) if lens0 == 'act600' else (lens0, None)
                simcat = []
                for k in combo:
                    r, d, C, Tm = results[(k, lens, lcut)]
                    ds.append(d); Ts.append(Tm)
                    sim = (simof(k, lens) / (T if lens == 'act' else 1))[:, :len(d)]
                    simcat.append(sim)
                blocks.append(np.cov(np.hstack(simcat).T))
                nsim = len(simcat[0])
            from scipy.linalg import block_diag
            C = block_diag(*blocks)
            d = np.concatenate(ds); Tm = np.hstack(Ts)
            # magnification: one (5s-2) per survey; their errors independent; build the prior covariance
            nper = [len(x) for x in ds]
            labels = [k for lens in maps_ for k in combo]
            Cm = np.zeros((len(d), len(d)))
            mu = np.zeros(len(d))
            off = np.cumsum([0] + nper)
            for i, ki in enumerate(labels):
                si, ei = slopes[ki]
                mu[off[i]:off[i + 1]] = (5 * si - 2) * Tm[1, off[i]:off[i + 1]]
                for j, kj in enumerate(labels):
                    if ki == kj:
                        sj, ej = slopes[kj]
                        Cm[off[i]:off[i + 1], off[j]:off[j + 1]] = 25 * ei * ej * np.outer(Tm[1, off[i]:off[i + 1]], Tm[1, off[j]:off[j + 1]])
            hart = (nsim - len(d) - 2) / (nsim - 1)
            Ci = hart * np.linalg.inv(C + Cm)
            r0 = d - mu
            F = Tm[0] @ Ci @ Tm[0]; bhat = Tm[0] @ Ci @ r0 / F; sb = 1 / np.sqrt(F)
            chi = (r0 - bhat * Tm[0]) @ Ci @ (r0 - bhat * Tm[0])
            lnL = lambda b: -0.5 * (r0 - b * Tm[0]) @ Ci @ (r0 - b * Tm[0])
            model = bhat * Tm[0] + mu
            Cn = hart * np.linalg.inv(C)
            snr = (model @ Cn @ d) / np.sqrt(model @ Cn @ model)
            from scipy.stats import chi2
            lines.append(f'- {" + ".join(combo)} on {" + ".join(LENS[m] for m in maps_)}: b = {bhat:.2f} ± {sb:.2f}; χ² {chi:.1f}/{len(d) - 1} '
                         f'(PTE {chi2.sf(chi, len(d) - 1):.2f}); S/N {snr:.1f}; L(13)/L(4) = {np.exp(lnL(13) - lnL(4)):.2e}')
            results[('comb', combo, maps_)] = (bhat, sb, chi, snr, np.exp(lnL(13) - lnL(4)))
    # nulls: curl
    lines.append('')
    lines.append('## Nulls')
    for k in ['desi', 'dr16', 'desi_mid']:
        d = D[f'{k}|curl']
        sim = (Sc[k] if Sc is not None else S[k])
        C = np.cov(sim.T); hart = (len(sim) - len(d) - 2) / (len(sim) - 1)
        from scipy.stats import chi2
        c2 = d @ (hart * np.linalg.inv(C)) @ d
        Tm = binned_theory(k if k != 'desi_mid' else 'desi_mid', 'act', k)
        F = Tm[0] @ (hart * np.linalg.inv(C)) @ Tm[0]; bc = Tm[0] @ (hart * np.linalg.inv(C)) @ d / F
        lines.append(f'- curl × {k}: χ² {c2:.1f}/{len(d)} (PTE {chi2.sf(c2, len(d)):.2f}); fitted "b" {bc:.2f} ± {1 / np.sqrt(F):.2f} '
                     f'(covariance from {"curl" if Sc is not None else "baseline"} sims)')
    # ACT vs Planck on the same quasars (difference of b; independent-sim covariance, conservative)
    for k in ['desi', 'dr16', 'desi_mid']:
        ra = results[(k, 'act', None)][0]; rp = results[(k, 'pin', None)][0]
        diff = ra['b'] - rp['b']; sd = np.hypot(ra['sb'], rp['sb'])
        lines.append(f'- ACT − Planck on the ACT footprint, {k}: b {ra["b"]:.2f} ± {ra["sb"]:.2f} vs {rp["b"]:.2f} ± {rp["sb"]:.2f}, '
                     f'Δb = {diff:.2f} ± {sd:.2f} ({diff / sd:+.1f}σ; errors added as independent: conservative, the shared κ signal '
                     'and quasar map make the difference smaller)')
    open(os.path.join(OUT, 'fit.md'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def bright(args):
    """post hoc (not in the fixed plan): DR16Q at i_dered < 20.2, Shen et al. 2007's limit, on Planck PR3 (full, and outside
    and inside ACT) with its own Gaussian sims; slope from the same intrinsic-count method on DR16Q's i counts."""
    import healpy as hp, fitsio, pymaster as nmt
    import qso_z4_theory as T
    d = fitsio.read(os.path.join(DATA, 'qso', 'DR16Q_v4.fits'), columns=['RA', 'DEC', 'Z', 'PSFMAG', 'EXTINCTION'])
    i = d['PSFMAG'][:, 3] - d['EXTINCTION'][:, 3]
    s = (d['Z'] > 3.5) & (d['Z'] < 5)
    g = np.load(os.path.join(OUT, 'maps.npz'))
    w = g['dr16_w']; npix = len(w)
    sel = s & (i < 20.2)
    N = np.bincount(hp.ang2pix(NSIDE, d['RA'][sel], d['DEC'][sel], lonlat=True), minlength=npix).astype(float)
    reg = (gal_lat(NSIDE) > 0).astype(int)
    delta = to_delta(N, w, reg)
    log(f'bright: {sel.sum()} DR16Q quasars at 3.5<z<5, i<20.2; {N[w > 0].sum():.0f} in the footprint')
    # slope: the intrinsic i counts of all DR16Q 3.5<z<5 (complete at i < 20.2 by construction of the cut: the cut is the
    # hard limit here), s = cumulative slope at i = 20.2 over 0.25 mag, error Poisson + dm spread, floor 0.05
    vals = []
    for dm in (0.15, 0.25, 0.35, 0.5):
        n1 = np.sum(i[s] < 20.2); n0 = np.sum(i[s] < 20.2 - dm); vals.append((np.log10(n1) - np.log10(n0)) / dm)
    sb = vals[1]; esb = max(np.std(vals), 0.05)
    log(f'bright slope s = {sb:.3f} ± {esb:.3f} (by dm {", ".join(f"{v:.3f}" for v in vals)})')
    th = T.Theory(lmax=LMAXK)
    h, e_ = np.histogram(d['Z'][sel], bins=np.arange(3.475, 5.05, 0.05)); zc = 0.5 * (e_[1:] + e_[:-1])
    sp = th.spectra(zc, h.astype(float)); np.savez(os.path.join(OUT, 'theory_dr16b.npz'), zc=zc, nz=h, **sp)
    gf = nmt.NmtField(w, [delta], lmax=LMAXK, n_iter=1)
    nl = np.loadtxt(os.path.join(DATA, 'planck', 'COM_Lensing_4096_R3.00', 'MV', 'nlkk.dat'))
    lmax = hp_lmax_planck(); cl = np.zeros(lmax + 1); cl[nl[:lmax + 1, 0].astype(int)] = nl[:lmax + 1, 2]
    wins = {'planck': planck_window()}
    wins['pout'] = wins['planck'] * (act_window('baseline') <= 0); wins['pin'] = wins['planck'] * (act_window('baseline') > 0)
    ws = {}
    for lens, wk in wins.items():
        fn = os.path.join(OUT, f'ws_dr16b_{lens}.fits'); W = nmt.NmtWorkspace()
        if os.path.exists(fn):
            W.read_from(fn)
        else:
            W.compute_coupling_matrix(nmt.NmtField(wk, None, lmax=LMAXK, spin=0), gf, bins(lens)); W.write_to(fn)
        ws[lens] = W
    def cross(kmap):
        return {lens: ws[lens].decouple_cell(nmt.compute_coupled_cell(nmt.NmtField(wk, [kmap], lmax=LMAXK, n_iter=1), gf))[0][:-1]
                for lens, wk in wins.items()}
    dat = cross(planck_kappa())
    rng = np.random.default_rng(1391); sims_ = {l_: [] for l_ in wins}
    for k in range(args.n):
        np.random.seed(int(rng.integers(1 << 31)))
        km = hp.alm2map(hp.synalm(cl, lmax=lmax), NSIDE, lmax=lmax)
        for l_, v in cross(km).items():
            sims_[l_].append(v)
    lines = ['', '## Post hoc: DR16Q bright (i < 20.2, Shen et al. 2007 limit), Planck PR3', f'N = {N[w > 0].sum():.0f}; s = {sb:.3f} ± {esb:.3f}']
    for lens in wins:
        Tm = binned_theory('dr16b', lens, None, ws=ws[lens])
        sim = np.array(sims_[lens]); C = np.cov(sim.T)
        r = fit_one(dat[lens], C, Tm, sb, esb, len(sim))
        lines.append(f'- {LENS[lens]}: b = {r["b"]:.2f} ± {r["sb"]:.2f}; χ² {r["chi2"]:.1f}/{r["dof"]} (PTE {r["pte"]:.2f}); S/N {r["snr"]:.1f}; '
                     f'L(13)/L(4) = {r["lr13_4"]:.2e}')
    open(os.path.join(OUT, 'bright.md'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd')
    ap.add_argument('--map', default='baseline')
    ap.add_argument('--nsims', type=int, default=400)
    ap.add_argument('--ntransfer', type=int, default=60)
    ap.add_argument('--n', type=int, default=400)
    a = ap.parse_args()
    {'prep': prep, 'data': data, 'sims': sims, 'gsims': gsims, 'fit': fit, 'bright': bright}[a.cmd](a)


if __name__ == '__main__':
    main()
