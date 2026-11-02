"""eRASS1 cluster masses confronted with ACT DR6 CMB lensing (R-122, cosmo-deep-2).

The Planck PR3 version is `erass1_cmblens.py` (R-121); this script reuses its catalogue handling, halo model and fits, and
changes what the three referees required: the predicted kappa per cluster is averaged over the cluster's mass posterior, the
two-halo amplitude is free (with a halo-exclusion form as a check), the covariance and the null distribution come from the
release's 400 simulated reconstructions, and the null-mean correction from mask-matched random positions.

Map: ACT DR6 lensing v1 (Qu+2024 arXiv:2304.05202, Madhavacheril+2024 arXiv:2304.05203; sims Farren+2024 arXiv:2309.05659),
kappa alm in equatorial coordinates, mean-field subtracted and normalised, band 2-3000; the release mask (nside 4096, apodized,
reject mask < 0.99). Band used here: LMIN <= L <= LMAX with a cosine taper, synthesised at nside 1024 for the annuli inside 30'
and at nside 512 for the annuli outside (the pixel windows are in the model).

usage (python -I; all data under data/, gitignored):
  python -I erass1_act.py prep  [--variant baseline]     # selection, counts by mass bin, stacking matrices, randoms
  python -I erass1_act.py data  [--variant baseline]     # the data stack and the random stack on the data map
  python -I erass1_act.py sims  [--variant baseline]     # the same stack on each simulation in data/act_dr6/sims (streamed)
  python -I erass1_act.py rotations [--nrot 200]          # rotated-map nulls (map and mask turned in RA, the clusters fixed)
  python -I erass1_act.py fit   [--variant baseline] [--cscale 1.0] [--exclusion]
  python -I erass1_act.py planck                          # the Planck PR3 tSZ-deprojected stack on the same clusters
"""
import argparse, glob, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import erass1_cmblens as r121  # noqa: E402

DATA = os.path.join(HERE, 'data')
ACT = os.path.join(DATA, 'act_dr6')
OUT = os.path.join(DATA, 'erass1_act')
ARCMIN = np.pi / 180 / 60
EDGES_IN = np.array([0, 5, 10, 15, 20, 30.0])        # arcmin, nside 1024
EDGES_OUT = np.array([30, 50, 80, 110, 150.0])       # arcmin, nside 512
EDGES = np.concatenate([EDGES_IN, EDGES_OUT[1:]])
NB = len(EDGES) - 1
NSIDE_IN, NSIDE_OUT = 1024, 512
LMIN, LMAX = 20, 2048
MASS_BINS = [(0, 1e14), (1e14, 2e14), (2e14, 4e14), (4e14, 6e14), (6e14, 1e16)]

# the samples and the statistic, fixed before the first fit (R-122 brief, step 5): the one-halo amplitude with the
# two-halo amplitude free, and b_M with the two-halo amplitude profiled out, per sample
SAMPLES = ['all nosz', 'M<4e14 nosz', 'cosmo nosz', 'cosmo M<4e14 nosz']   # outside the release's SZ-cluster (inpainting) mask
INFO_SUBSETS = ['all', 'M<4e14', 'cosmo', 'cosmo M<4e14', 'insz', 'insz M<4e14', 'insz M>4e14', 'M<2e14 nosz', 'cosmo M<2e14 nosz',
                'M<1e14 nosz', 'M1-2e14 nosz', 'M2-4e14 nosz', 'M>4e14 nosz', 'M4-6e14 nosz', 'M>6e14 nosz', 'z<0.1 nosz',
                'z0.1-0.3 nosz', 'z0.3-0.6 nosz', 'z>0.6 nosz', 'PCONT<0.1 nosz', 'cosmo M>4e14 nosz',
                'M<1e14', 'M1-2e14', 'M2-4e14', 'M>4e14', 'z<0.1', 'z0.1-0.3', 'z0.3-0.6', 'z>0.6', 'PCONT<0.1', 'nopsz',
                'cosmo M>4e14', 'M4-6e14', 'M>6e14']


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def act_file(kind, variant):
    if kind == 'alm':
        return os.path.join(ACT, f'kappa_alm_data_act_dr6_lensing_v1_{variant}.fits')
    if kind == 'mask':
        # the three foreground variants share the baseline analysis mask; galcut040 has its own
        v = variant if variant in ('baseline', 'galcut040') else 'baseline'
        return os.path.join(ACT, f'mask_act_dr6_lensing_v1_healpix_nside_{4096}_{v}.fits')


def mask_tag(variant):
    return variant if variant in ('galcut040', 'planck') else 'baseline'


# ------------------------------------------------------------------------------------------------------------ masks
def binary_masks(variant):
    """the release mask degraded to the two working resolutions, binary (clean = every sub-pixel >= 0.99 -> mean >= 0.99)."""
    import healpy as hp
    fn = os.path.join(OUT, f'masks_{mask_tag(variant)}.npz')
    if os.path.exists(fn):
        g = np.load(fn)
        return g['b_in'], g['b_out']
    log('reading the nside-4096 mask')
    m = hp.read_map(act_file('mask', variant), dtype=np.float32)
    b_in = hp.ud_grade(m, NSIDE_IN) >= 0.99
    b_out = hp.ud_grade(m, NSIDE_OUT) >= 0.99
    log(f'mask: fsky(>=0.99) at 4096 {np.mean(m >= 0.99):.4f}, at {NSIDE_IN} {b_in.mean():.4f}, at {NSIDE_OUT} {b_out.mean():.4f}')
    os.makedirs(OUT, exist_ok=True)
    np.savez(fn, b_in=b_in, b_out=b_out)
    return b_in, b_out


def disc_pixels(nside, vec, radius_rad):
    import healpy as hp
    p = hp.query_disc(nside, vec, radius_rad)
    pv = np.array(hp.pix2vec(nside, p))
    d = np.arccos(np.clip(vec @ pv, -1, 1))
    return p, d


def footprint(lon, lat, b_in, b_out):
    """f30: clean fraction within 30' (nside 1024); f120: within 2 deg (nside 512); keep = f30 >= 0.99 & f120 >= 0.8."""
    import healpy as hp
    vec = hp.ang2vec(lon, lat, lonlat=True)
    f30 = np.zeros(len(lon)); f120 = np.zeros(len(lon))
    for i in range(len(lon)):
        f30[i] = b_in[hp.query_disc(NSIDE_IN, vec[i], 0.5 * np.pi / 180)].mean()
        f120[i] = b_out[hp.query_disc(NSIDE_OUT, vec[i], 2.0 * np.pi / 180)].mean()
    return f30, f120, (f30 >= 0.99) & (f120 >= 0.8)


def stacking_matrices(lon, lat):
    """sparse (n*nb_in x npix_in) and (n*nb_out x npix_out) matrices of ones: row (i, k) holds the pixels of annulus k of
    object i; the masked mean is (S @ (m * b)) / (S @ b). Also the 30' disc (nside_in) and 2 deg disc (nside_out) selection
    matrices."""
    import healpy as hp
    from scipy.sparse import csr_matrix
    vec = hp.ang2vec(lon, lat, lonlat=True)
    n = len(lon)
    rows_in, cols_in, rows_out, cols_out = [], [], [], []
    rows30, cols30, rows120, cols120 = [], [], [], []
    e_in = EDGES_IN * ARCMIN; e_out = EDGES_OUT * ARCMIN
    t0 = time.time()
    for i in range(n):
        p, d = disc_pixels(NSIDE_IN, vec[i], e_in[-1])
        k = np.searchsorted(e_in, d, side='right') - 1
        ok = (k >= 0) & (k < len(e_in) - 1)
        rows_in.append(i * (len(e_in) - 1) + k[ok]); cols_in.append(p[ok])
        p2, d2 = disc_pixels(NSIDE_OUT, vec[i], e_out[-1])
        k2 = np.searchsorted(e_out, d2, side='right') - 1
        ok2 = (k2 >= 0) & (k2 < len(e_out) - 1)
        rows_out.append(i * (len(e_out) - 1) + k2[ok2]); cols_out.append(p2[ok2])
        rows30.append(np.full(len(p), i)); cols30.append(p)
        p3 = hp.query_disc(NSIDE_OUT, vec[i], 2.0 * np.pi / 180)
        rows120.append(np.full(len(p3), i)); cols120.append(p3)
        if (i + 1) % 2000 == 0:
            log(f'  matrices {i + 1}/{n} ({time.time() - t0:.0f}s)')
    def mk(r, c, nrow, npix):
        r = np.concatenate(r); c = np.concatenate(c)
        return csr_matrix((np.ones(len(r), np.float32), (r, c)), shape=(nrow, npix))
    npix_in = 12 * NSIDE_IN ** 2; npix_out = 12 * NSIDE_OUT ** 2
    return (mk(rows_in, cols_in, n * (len(e_in) - 1), npix_in), mk(rows_out, cols_out, n * (len(e_out) - 1), npix_out),
            mk(rows30, cols30, n, npix_in), mk(rows120, cols120, n, npix_out))


def random_weights(b_in, b_out, nrand, seed, dec_range, ra_ranges=None):
    """mask-matched randoms: uniform positions over the footprint's Dec range that pass the same selection; their stacking
    weights accumulated into dense per-annulus weight maps (1 / n_clean_pixels per random, averaged over randoms), with the
    fixed data mask; the random stack of any map on this mask is W_in @ m_in, W_out @ m_out."""
    import healpy as hp
    rng = np.random.default_rng(seed)
    W_in = np.zeros((len(EDGES_IN) - 1, 12 * NSIDE_IN ** 2), np.float32)
    W_out = np.zeros((len(EDGES_OUT) - 1, 12 * NSIDE_OUT ** 2), np.float32)
    e_in = EDGES_IN * ARCMIN; e_out = EDGES_OUT * ARCMIN
    got = 0; tried = 0; f120s = []; pos = []
    t0 = time.time()
    while got < nrand:
        ra = rng.uniform(0, 360); dec = np.degrees(np.arcsin(rng.uniform(np.sin(np.radians(dec_range[0])), np.sin(np.radians(dec_range[1])))))
        if ra_ranges is not None and not any(lo <= ra < hi for lo, hi in ra_ranges):
            continue
        tried += 1
        vec = hp.ang2vec(ra, dec, lonlat=True)
        p30 = hp.query_disc(NSIDE_IN, vec, 0.5 * np.pi / 180)
        if b_in[p30].mean() < 0.99:
            continue
        p120 = hp.query_disc(NSIDE_OUT, vec, 2.0 * np.pi / 180)
        f120 = b_out[p120].mean()
        if f120 < 0.8:
            continue
        p, d = disc_pixels(NSIDE_IN, vec, e_in[-1]); k = np.searchsorted(e_in, d, side='right') - 1
        ok = (k >= 0) & (k < len(e_in) - 1) & b_in[p]
        nk = np.bincount(k[ok], minlength=len(e_in) - 1)
        W_in[k[ok], p[ok]] += (1.0 / np.maximum(nk, 1))[k[ok]]
        p, d = disc_pixels(NSIDE_OUT, vec, e_out[-1]); k = np.searchsorted(e_out, d, side='right') - 1
        ok = (k >= 0) & (k < len(e_out) - 1) & b_out[p]
        nk = np.bincount(k[ok], minlength=len(e_out) - 1)
        W_out[k[ok], p[ok]] += (1.0 / np.maximum(nk, 1))[k[ok]]
        got += 1; f120s.append(f120); pos.append((ra, dec))
        if got % 2000 == 0:
            log(f'  randoms {got}/{nrand} (tried {tried}, {time.time() - t0:.0f}s)')
    return W_in / got, W_out / got, np.array(f120s), np.array(pos), tried


def mass_bin_table(M, masks, labels):
    lines = []
    for lo, hi in MASS_BINS:
        s = (M >= lo) & (M < hi)
        lines.append(f'  M500 {lo:.0e}-{hi:.0e}: ' + ', '.join(f'{lab} {int((s & m).sum())}' for lab, m in zip(labels, masks)))
    return '\n'.join(lines)


# ------------------------------------------------------------------------------------------------------------- prep
def prep(args):
    import healpy as hp
    from scipy.sparse import save_npz
    os.makedirs(OUT, exist_ok=True)
    tag = mask_tag(args.variant)
    c = r121.load_cat()
    b_in, b_out = binary_masks(args.variant)
    lon, lat = c['RA'], c['DEC']
    log(f'catalogue: {len(lon)} clusters (PCONT < 0.5, M500 > 0, z > 0.01); footprint selection ({tag} mask)')
    f30, f120, keep = footprint(lon, lat, b_in, b_out)
    inreg = f120 >= 0.8
    near = f120 > 0.0
    M = c['M500']
    log(f'selection: {keep.sum()} kept of {inreg.sum()} with >= 80 % clean within 2 deg ({near.sum()} touch the footprint; '
        f'{len(lon)} in the catalogue); cosmology sample kept {int((keep & c["cosmo"]).sum())} of {int((inreg & c["cosmo"]).sum())}')
    centre_hole = inreg & (b_in[hp.ang2pix(NSIDE_IN, lon, lat, lonlat=True)] == 0)
    log('by mass bin: in the footprint region (f120 >= 0.8), kept (clean within 30\'), removed by holes, centre itself in a hole:\n' +
        mass_bin_table(M, [inreg, keep, inreg & ~keep, centre_hole], ['region', 'kept', 'hole-removed', 'centre-in-hole']))
    log('cosmology sample by mass bin:\n' + mass_bin_table(M, [inreg & c['cosmo'], keep & c['cosmo'], inreg & c['cosmo'] & ~keep],
                                                           ['region', 'kept', 'hole-removed']))
    for lo, hi in MASS_BINS:
        s = inreg & (M >= lo) & (M < hi)
        if s.sum():
            log(f'  hole-removed fraction M500 {lo:.0e}-{hi:.0e}: {1 - keep[s].mean():.3f} (n {s.sum()}), '
                f'PSZ2/ACT/SPT-matched among them {int((s & (c["psz"] | c["act"])).sum())}, of which kept {int((s & keep & (c["psz"] | c["act"])).sum())}')
    # the SZ-cluster mask the release provides separately (6' holes, nside 8192): which eRASS1 clusters sit in it
    szfn = os.path.join(ACT, 'sz_cluster_mask_radius_6.0_arcmin_nside_8192.fits')
    if os.path.exists(szfn) and os.path.getsize(szfn) > 7e8 and not os.path.exists(os.path.join(OUT, 'insz_baseline.npy')):
        try:
            sz = hp.read_map(szfn, dtype=np.float32)
            insz = sz[hp.ang2pix(8192, lon, lat, lonlat=True)] < 0.5
            del sz
            log(f'SZ-cluster mask (release other/): fsky(<0.5) ; eRASS1 clusters in region inside an SZ hole: {int((inreg & insz).sum())}, '
                f'of which kept by our selection {int((keep & insz).sum())}\n' +
                mass_bin_table(M, [inreg & insz, keep & insz], ['in-SZ-hole', 'kept-but-in-SZ-hole']))
            np.save(os.path.join(OUT, f'insz_{tag}.npy'), insz)
        except Exception as e:  # noqa: BLE001
            log('SZ mask not read:', e)
    log('building the stacking matrices for the kept clusters')
    S_in, S_out, S30, S120 = stacking_matrices(lon[keep], lat[keep])
    save_npz(os.path.join(OUT, f'S_in_{tag}.npz'), S_in); save_npz(os.path.join(OUT, f'S_out_{tag}.npz'), S_out)
    save_npz(os.path.join(OUT, f'S30_{tag}.npz'), S30); save_npz(os.path.join(OUT, f'S120_{tag}.npz'), S120)
    np.savez(os.path.join(OUT, f'selection_{tag}.npz'), keep=keep, f30=f30, f120=f120, inreg=inreg, lon=lon, lat=lat)
    log('mask-matched randoms')
    dec_range = (lat[keep].min() - 1, lat[keep].max() + 1)
    W_in, W_out, f120r, pos, tried = random_weights(b_in, b_out, args.nrand, 5, dec_range)
    log(f'randoms: {len(pos)} of {tried} tried in Dec {dec_range[0]:.1f}..{dec_range[1]:.1f}; mean f120 randoms {f120r.mean():.3f} '
        f'vs clusters {f120[keep].mean():.3f}; f30 of randoms = 1 by construction')
    np.savez(os.path.join(OUT, f'randoms_{tag}.npz'), W_in=W_in, W_out=W_out, f120=f120r, pos=pos, dec_range=dec_range)
    log('prep done')


# ------------------------------------------------------------------------------------------------------------- maps
def band_filter(ell, lmin, lmax):
    f = np.zeros(len(ell))
    sel = (ell >= lmin) & (ell <= lmax)
    f[sel] = 1.0
    tap = (ell > 0.9 * lmax) & sel
    f[tap] = 0.5 * (1 + np.cos(np.pi * (ell[tap] - 0.9 * lmax) / (0.1 * lmax)))
    low = (ell >= lmin) & (ell < lmin + 4)
    f[low] *= 0.5 * (1 - np.cos(np.pi * (ell[low] - lmin + 1) / 5))
    return f


def read_alm(fn):
    import healpy as hp
    alm = hp.read_alm(fn).astype(np.complex128)
    bad = ~np.isfinite(alm)
    if bad.any():
        alm[bad] = 0  # the release's L = 0, 1 entries are NaN; they are outside the band anyway
    return alm


def maps_from_alm(alm, lmin=LMIN, lmax=LMAX, dphi=0.0):
    """filtered maps at the two working resolutions; dphi (rad) turns the map in longitude (phi -> phi + dphi)."""
    import healpy as hp
    lm = hp.Alm.getlmax(len(alm))
    ell = np.arange(lm + 1)
    a = hp.almxfl(alm, band_filter(ell, lmin, lmax))
    if dphi != 0.0:
        l, m = hp.Alm.getlm(lm)
        a = a * np.exp(-1j * m * dphi)
    m_in = hp.alm2map(a, NSIDE_IN, lmax=lm)
    m_out = hp.alm2map(a, NSIDE_OUT, lmax=lm)
    return m_in, m_out


def rotate_binary(b, nside, dphi):
    """nearest-pixel rotation of a binary map by dphi in longitude (the map's value at phi - dphi)."""
    import healpy as hp
    th, ph = hp.pix2ang(nside, np.arange(len(b)))
    return b[hp.ang2pix(nside, th, ph - dphi)]


class Stacker:
    def __init__(self, tag):
        from scipy.sparse import load_npz
        self.S_in = load_npz(os.path.join(OUT, f'S_in_{tag}.npz')); self.S_out = load_npz(os.path.join(OUT, f'S_out_{tag}.npz'))
        self.S30 = load_npz(os.path.join(OUT, f'S30_{tag}.npz')); self.S120 = load_npz(os.path.join(OUT, f'S120_{tag}.npz'))
        self.n = self.S30.shape[0]
        self.n30 = np.asarray(self.S30.sum(1)).ravel(); self.n120 = np.asarray(self.S120.sum(1)).ravel()
        R = np.load(os.path.join(OUT, f'randoms_{tag}.npz'))
        self.W_in, self.W_out = R['W_in'], R['W_out']

    def profiles(self, m_in, m_out, b_in, b_out):
        nb_in = len(EDGES_IN) - 1; nb_out = len(EDGES_OUT) - 1
        with np.errstate(invalid='ignore', divide='ignore'):
            pi = (self.S_in @ (m_in * b_in)) / (self.S_in @ b_in.astype(np.float32))
            po = (self.S_out @ (m_out * b_out)) / (self.S_out @ b_out.astype(np.float32))
        return np.hstack([pi.reshape(self.n, nb_in), po.reshape(self.n, nb_out)])

    def random_stack(self, m_in, m_out):
        return np.concatenate([self.W_in @ m_in, self.W_out @ m_out])

    def selection(self, b_in, b_out):
        f30 = (self.S30 @ b_in.astype(np.float32)) / self.n30
        f120 = (self.S120 @ b_out.astype(np.float32)) / self.n120
        return (f30 >= 0.99) & (f120 >= 0.8)


def subset_masks(c, keep):
    """boolean subsets over the kept clusters, from R-121's definitions plus the mass bins above 4e14."""
    M = c['M500'][keep]; z = c['BEST_Z'][keep]
    s = r121.subsets(M, z, c, keep)
    s['all'] = np.ones(len(M), bool)
    s['cosmo M>4e14'] = s['cosmo'] & (M >= 4e14)
    s['M4-6e14'] = (M >= 4e14) & (M < 6e14); s['M>6e14'] = M >= 6e14
    # the release's SZ-cluster mask (6' holes, inpainted in the CMB maps before the reconstruction): the clusters inside it
    # are a separate sample; every primary sample has a version outside it ("nosz")
    fn = os.path.join(OUT, 'insz_baseline.npy')
    insz = np.load(fn)[keep] if os.path.exists(fn) else np.zeros(len(M), bool)
    s['insz'] = insz; s['insz M<4e14'] = insz & (M < 4e14); s['insz M>4e14'] = insz & (M >= 4e14)
    for name in ['all', 'M<4e14', 'cosmo', 'cosmo M<4e14', 'M<1e14', 'M1-2e14', 'M2-4e14', 'M>4e14', 'z<0.1', 'z0.1-0.3', 'z0.3-0.6',
                 'z>0.6', 'PCONT<0.1', 'cosmo M>4e14', 'M4-6e14', 'M>6e14']:
        s[name + ' nosz'] = s[name] & ~insz
    s['M<2e14 nosz'] = (M < 2e14) & ~insz; s['cosmo M<2e14 nosz'] = s['cosmo'] & (M < 2e14) & ~insz
    return s


def data(args):
    tag = mask_tag(args.variant)
    c = r121.load_cat()
    sel = np.load(os.path.join(OUT, f'selection_{tag}.npz')); keep = sel['keep']
    b_in, b_out = binary_masks(args.variant)
    st = Stacker(tag)
    t0 = time.time()
    m_in, m_out = maps_from_alm(read_alm(act_file('alm', args.variant)), args.lmin, args.lmax)
    log(f'{args.variant}: maps built ({time.time() - t0:.0f}s); rms in the clean footprint {m_in[b_in].std():.3e} / {m_out[b_out].std():.3e}')
    prof = st.profiles(m_in, m_out, b_in, b_out)
    rs = st.random_stack(m_in, m_out)
    log(f'data stack x1e4: {np.round(np.nanmean(prof, 0) * 1e4, 3)}')
    log(f'random stack x1e4: {np.round(rs * 1e4, 3)}')
    # the data map's mean in the clean footprint, for comparison with the randoms
    log(f'map mean in clean pixels x1e4: {m_in[b_in].mean() * 1e4:.3f} (nside {NSIDE_IN}), {m_out[b_out].mean() * 1e4:.3f} (nside {NSIDE_OUT})')
    np.savez(os.path.join(OUT, f'data_{args.variant}{args.tag}.npz'), prof=prof, random=rs, lmin=args.lmin, lmax=args.lmax,
             map_rms=(m_in[b_in].std(), m_out[b_out].std()))


def sims(args):
    """the same stack on each simulated reconstruction (streamed into data/act_dr6/sims by simstream.sh); per-subset means and
    the random stack per simulation, appended to sims_<variant>.npz; each file deleted once stacked."""
    tag = mask_tag(args.variant)
    c = r121.load_cat()
    sel = np.load(os.path.join(OUT, f'selection_{tag}.npz')); keep = sel['keep']
    b_in, b_out = binary_masks(args.variant)
    st = Stacker(tag)
    S = subset_masks(c, keep)
    names = sorted(S)
    fn = os.path.join(OUT, f'sims_{args.variant}{args.tag}.npz')
    done = {}
    if os.path.exists(fn):
        g = np.load(fn, allow_pickle=True)
        done = {int(i): (m, r) for i, m, r in zip(g['idx'], g['means'], g['randoms'])}
    simdir = os.path.join(ACT, 'sims')
    idle = 0
    while len(done) < args.nsims and idle < args.idle_max:
        files = sorted(glob.glob(os.path.join(simdir, f'kappa_alm_sim_act_dr6_lensing_v1_{args.variant}_*.fits')))
        files = [f for f in files if int(f[-9:-5]) not in done]
        if not files:
            if os.path.exists(os.path.join(ACT, 'simstream.log')) and 'ALLDONE' in open(os.path.join(ACT, 'simstream.log')).read():
                break
            time.sleep(10); idle += 1; continue
        idle = 0
        f = files[0]; i = int(f[-9:-5])
        t0 = time.time()
        try:
            alm = read_alm(f)
        except Exception as e:  # noqa: BLE001
            log(f'sim {i}: unreadable ({e}); removing'); os.remove(f); continue
        m_in, m_out = maps_from_alm(alm, args.lmin, args.lmax)
        prof = st.profiles(m_in, m_out, b_in, b_out)
        means = np.array([np.nanmean(prof[S[n]], 0) if S[n].sum() else np.full(NB, np.nan) for n in names])
        done[i] = (means, st.random_stack(m_in, m_out))
        os.remove(f); open(os.path.join(simdir, f'done_{i:04d}'), 'w').close()
        idx = sorted(done)
        np.savez(fn, idx=np.array(idx), means=np.array([done[j][0] for j in idx]), randoms=np.array([done[j][1] for j in idx]),
                 names=np.array(names))
        log(f'sim {i}: {len(done)} done ({time.time() - t0:.0f}s); all-sample stack x1e4 {np.round(means[names.index("all")] * 1e4, 3)}')
    log(f'sims: {len(done)} stacked')


def rotations(args):
    """rotated-map nulls: the map and the mask turned together by a random RA offset of 25-335 deg, the clusters fixed, the
    same footprint selection re-applied; per-subset means per rotation (R-121's rotated-catalogue null, in this frame)."""
    tag = mask_tag(args.variant)
    c = r121.load_cat()
    sel = np.load(os.path.join(OUT, f'selection_{tag}.npz')); keep = sel['keep']
    b_in, b_out = binary_masks(args.variant)
    st = Stacker(tag)
    S = subset_masks(c, keep); names = sorted(S)
    alm = read_alm(act_file('alm', args.variant))
    rng = np.random.default_rng(args.seed)
    fn = os.path.join(OUT, f'rotations_{args.variant}{args.tag}_seed{args.seed}.npz')
    means = []; ns = []; dphis = []
    for r in range(args.nrot):
        t0 = time.time()
        dphi = np.radians(rng.uniform(25, 335))
        m_in, m_out = maps_from_alm(alm, args.lmin, args.lmax, dphi=dphi)
        bi = rotate_binary(b_in, NSIDE_IN, dphi); bo = rotate_binary(b_out, NSIDE_OUT, dphi)
        kr = st.selection(bi, bo)
        prof = st.profiles(m_in, m_out, bi, bo)
        means.append(np.array([np.nanmean(prof[S[n] & kr], 0) if (S[n] & kr).sum() else np.full(NB, np.nan) for n in names]))
        ns.append(np.array([(S[n] & kr).sum() for n in names])); dphis.append(dphi)
        if (r + 1) % 5 == 0 or r == 0:
            np.savez(fn, means=np.array(means), ns=np.array(ns), dphi=np.array(dphis), names=np.array(names))
            log(f'rotation {r + 1}/{args.nrot}: kept {kr.sum()} of {len(kr)} ({time.time() - t0:.0f}s)')
    np.savez(fn, means=np.array(means), ns=np.array(ns), dphi=np.array(dphis), names=np.array(names))
    log('rotations done')


# ------------------------------------------------------------------------------------------------------------ model
class ExclusionModel(r121.HaloModel):
    """R-121's halo model with the two-halo term zeroed inside r_ex = ex_r200 x r200c (halo exclusion)."""

    def __init__(self, ex_r200=1.5, **kw):
        super().__init__(**kw)
        self.ex_r200 = ex_r200
        self._xi = {}

    def xi_lin(self, z):
        if z in self._xi:
            return self._xi[z]
        from scipy.special import spherical_jn
        r = np.logspace(-2, 2.7, 400)
        kk = np.logspace(-4, 1.5, 3000)
        P = self.plin(z, kk) * np.exp(-(kk / 10.0) ** 2)
        xi = np.array([np.trapezoid(kk ** 2 * P * spherical_jn(0, kk * ri), kk) / (2 * np.pi ** 2) for ri in r])
        self._xi[z] = (r, xi)
        return r, xi

    def kappa_theta(self, M500c, z, b_M=1.0):
        k1, k2, info = super().kappa_theta(M500c, z, b_M)
        M = b_M * M500c
        M200, c200, rs, rho_s = self.m500_to_nfw(M, z)
        r200 = c200 * rs
        rex = self.ex_r200 * r200 * (1 + z)  # comoving Mpc
        Dl = self.co.angular_diameter_distance(z).value
        Rc = self.theta * Dl * (1 + z)
        r, xi = self.xi_lin(z)
        lxi = lambda rr: np.interp(np.log(rr), np.log(r), xi)
        corr = np.zeros_like(Rc)
        inside = Rc < rex
        for i in np.where(inside)[0]:
            chimax = np.sqrt(rex ** 2 - Rc[i] ** 2)
            chi = np.linspace(0, chimax, 60)
            corr[i] = 2 * np.trapezoid(lxi(np.sqrt(Rc[i] ** 2 + chi ** 2)), chi)
        k2 = k2 - (1 + z) ** 2 * self.rho_m0 * info['bias'] * corr / self.sigma_crit(z)
        return k1, k2, info


class NonlinearModel(r121.HaloModel):
    """R-121's halo model with the two-halo term on the halofit (non-linear) P(k) instead of the linear one: at the 2-7 Mpc
    the 10-30' annuli probe, non-linear clustering exceeds linear theory by 1.5-2, and a free linear-shape two-halo amplitude
    absorbs that at the one-halo amplitude's expense (the shape systematic of the statistic)."""

    def __init__(self, **kw):
        import camb
        from scipy.interpolate import RectBivariateSpline
        super().__init__(**kw)
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=r121.H0, ombh2=r121.OB * (r121.H0 / 100) ** 2, omch2=(r121.OM - r121.OB) * (r121.H0 / 100) ** 2, mnu=0.06, tau=0.054)
        pars.InitPower.set_params(As=2.1e-9 * self.s8_scale, ns=0.965)   # sigma8 = 0.88 at the linear level, as R-121
        pars.set_matter_power(redshifts=list(np.linspace(0, 1.6, 17))[::-1], kmax=50.0)
        pars.NonLinear = camb.model.NonLinear_both
        res = camb.get_results(pars)
        kh, zs, pk = res.get_matter_power_spectrum(minkh=1e-4, maxkh=50, npoints=600)
        zs = np.array(zs); pk = np.array(pk) / self.h ** 3
        order = np.argsort(zs)
        self._pk_spl = RectBivariateSpline(zs[order], np.log(kh * self.h), np.log(pk[order]))
        self._lin_spl = RectBivariateSpline(self.zs, np.log(self.k), np.log(self.pk))

    def sigma_m(self, M, z):
        # the bias stays on the linear sigma(M)
        R = (3 * M / (4 * np.pi * self.rho_m0)) ** (1 / 3)
        k = self.k; x = k * R
        W = 3 * (np.sin(x) - x * np.cos(x)) / x ** 3
        P = np.exp(self._lin_spl(z, np.log(k), grid=False))
        return np.sqrt(np.trapezoid(k ** 2 * P * W ** 2, k) / (2 * np.pi ** 2))


def make_model(args):
    if args.nonlinear:
        hm = NonlinearModel(lmax=3000, cscale=args.cscale)
    elif args.exclusion:
        hm = ExclusionModel(ex_r200=args.ex_r200, lmax=3000, cscale=args.cscale)
    else:
        hm = r121.HaloModel(lmax=3000, cscale=args.cscale)
    gtag = ('_nl' if args.nonlinear else '') + ('_ex%.1f' % args.ex_r200 if args.exclusion else '') + ('_c%.2f' % args.cscale if args.cscale != 1.0 else '')
    return hm, gtag


def model_grid_here(hm, tag):
    old = r121.OUT; r121.OUT = OUT
    try:
        return r121.model_grid(hm, 1.0, tag=tag)
    finally:
        r121.OUT = old


GH_X, GH_W = np.polynomial.hermite.hermgauss(5)


def posterior_nodes(lnM, sig):
    """Gauss-Hermite nodes of a lognormal mass posterior: (n, 5) ln M values and the weights."""
    return lnM[:, None] + np.sqrt(2) * sig[:, None] * GH_X[None, :], GH_W / np.sqrt(np.pi)


def annulus_grid(hm, grid, filt_in, filt_out):
    """the filtered annulus means of every grid node (nM, nz, NB) for the 1h and 2h terms; the annulus mean is linear in
    kappa(L) and the grid interpolation is linear, so interpolating these equals interpolating kappa(L) first."""
    lnMg, zg, k1g, k2g, _ = grid
    ell = hm.ell
    def means(k):
        return np.concatenate([r121.HaloModel.annulus_means_from_ell(k, ell, filt_in, EDGES_IN),
                               r121.HaloModel.annulus_means_from_ell(k, ell, filt_out, EDGES_OUT)])
    A1 = np.zeros(k1g.shape[:2] + (NB,)); A2 = np.zeros_like(A1)
    for i in range(k1g.shape[0]):
        for j in range(k1g.shape[1]):
            A1[i, j] = means(k1g[i, j]); A2[i, j] = means(k2g[i, j])
    return lnMg, zg, A1, A2


def stacked_model_posterior(hm, agrid, lnM_nodes, w_nodes, z, filt_in=None, filt_out=None, bias_scale=1.0):
    """mean filtered annulus means of the model over clusters, each cluster's prediction averaged over its mass posterior
    (lnM_nodes: (n, nn) with weights w_nodes (nn,) or (n, nn)); agrid from annulus_grid."""
    lnMg, zg, A1, A2 = agrid
    n, nn = lnM_nodes.shape
    w = np.broadcast_to(np.asarray(w_nodes), (n, nn)) if np.ndim(w_nodes) == 1 else w_nodes
    m1 = np.zeros(NB); m2 = np.zeros(NB)
    for j in range(nn):
        m1 += (w[:, j, None] * r121.interp_grid(lnMg, zg, A1, lnM_nodes[:, j], z)).sum(0)
        m2 += (w[:, j, None] * r121.interp_grid(lnMg, zg, A2, lnM_nodes[:, j], z)).sum(0)
    return m1 / n, m2 / n * bias_scale


def filters(hm, lmin, lmax):
    import healpy as hp
    ell = hm.ell
    f = band_filter(ell, lmin, lmax)
    def pw(nside):
        p = hp.pixwin(nside, lmax=min(4 * nside, len(ell) - 1))
        return np.concatenate([p, np.zeros(len(ell) - len(p))])
    return f * pw(NSIDE_IN), f * pw(NSIDE_OUT)


def profile_chi2(d, Ci, m1, m2):
    """chi2 minimised over the two-halo amplitude, for a one-halo template m1 fixed (A_1h = 1): returns (chi2, A_2h)."""
    r = d - m1
    a2 = (m2 @ Ci @ r) / (m2 @ Ci @ m2)
    rr = r - a2 * m2
    return rr @ Ci @ rr, a2


def interval(bgrid, chi, dchi):
    """the b_M values where chi2 - min < dchi, by linear interpolation of the crossings (grid edge if not crossed)."""
    ib = np.argmin(chi); c0 = chi[ib]
    lo = bgrid[0]; hi = bgrid[-1]
    for i in range(ib, 0, -1):
        if chi[i - 1] - c0 >= dchi:
            lo = np.interp(dchi, [chi[i] - c0, chi[i - 1] - c0], [bgrid[i], bgrid[i - 1]]); break
    for i in range(ib, len(bgrid) - 1):
        if chi[i + 1] - c0 >= dchi:
            hi = np.interp(dchi, [chi[i] - c0, chi[i + 1] - c0], [bgrid[i], bgrid[i + 1]]); break
    return bgrid[ib], lo, hi


def fit(args):
    """stacks, covariance from the simulations, the fits per sample; every amplitude computed is printed (the trials record)."""
    tag = mask_tag(args.variant)
    c = r121.load_cat()
    sel = np.load(os.path.join(OUT, f'selection_{tag}.npz')); keep = sel['keep']
    D = np.load(os.path.join(OUT, f'data_{args.variant}{args.tag}.npz'))
    prof, rstack = D['prof'], D['random']
    lmin, lmax = int(D['lmin']), int(D['lmax'])
    simfn = os.path.join(OUT, f'sims_{args.simvariant or args.variant}{args.tag}.npz')
    rot = None
    rf = sorted(glob.glob(os.path.join(OUT, f'rotations_{args.variant}{args.tag}_seed*.npz')))
    if rf:
        R = [np.load(f, allow_pickle=True) for f in rf]
        rot = dict(means=np.concatenate([r['means'] for r in R]), names=list(R[0]['names']))
    if args.cov == 'sims':
        SM = np.load(simfn, allow_pickle=True)
        names = list(SM['names']); simmeans = SM['means']; simrand = SM['randoms']
    else:
        # covariance from the rotations: each rotation's stack rescaled to the sample's size (noise-dominated stacks scale as
        # 1/n), used where a variant has no simulations (polonly, tonly, cibdeproj, galcut040)
        names = rot['names']; simmeans = rot['means']; simrand = np.zeros((len(simmeans), NB))
    nsim = len(simmeans)
    S = subset_masks(c, keep)
    z = c['BEST_Z'][keep]; M = c['M500'][keep]
    sig = 0.5 * (np.log(c['M500_H'][keep]) - np.log(c['M500_L'][keep]))
    sig = np.where(np.isfinite(sig) & (sig > 0), sig, np.nanmedian(sig))
    log(f'{args.variant}{args.tag}: {keep.sum()} clusters, {nsim} simulations, {0 if rot is None else len(rot["means"])} rotations, '
        f'band {lmin}-{lmax}; quoted sigma_lnM median {np.median(sig):.3f}')
    hm, gtag = make_model(args)
    f_in, f_out = filters(hm, lmin, lmax)
    grid = annulus_grid(hm, model_grid_here(hm, gtag), f_in, f_out)
    mid = 0.5 * (EDGES[1:] + EDGES[:-1])
    # the catalogue's M500_PDF as weights on its own grid (the check on the lognormal)
    from astropy.io import fits as pf
    d = pf.open(r121.CAT)[1].data
    ok = (d['PCONT'] < 0.5) & np.isfinite(d['M500']) & (d['M500'] > 0) & (d['BEST_Z'] > 0.01)
    pdf = np.asarray(d['M500_PDF'][ok][keep], float); pgrid = np.asarray(d['M500_PDF_array'][ok][keep], float)  # Msun
    pdfw = pdf / np.maximum(pdf.sum(1, keepdims=True), 1e-300)
    bad = ~np.isfinite(pdfw).all(1) | ~(pgrid > 0).all(1)
    log(f'M500_PDF grid: {pgrid.shape[1]} points per cluster; {bad.sum()} clusters without a usable PDF (median mass used)')
    bgrid = np.arange(0.1, 2.001, 0.05)
    results = {}
    FULL = set(SAMPLES + ['all', 'M<4e14', 'cosmo', 'cosmo M<4e14', 'insz', 'insz M<4e14', 'M<2e14 nosz', 'cosmo M<2e14 nosz',
                          'M2-4e14 nosz', 'M1-2e14 nosz'])
    rot_ns = None
    if rot is not None:
        rot_ns = np.concatenate([r['ns'] for r in R])

    def covs(name):
        j = names.index(name)
        A = simmeans[:, j]
        ok = np.all(np.isfinite(A), 1)
        if args.cov == 'rotations':
            nr = rot_ns[:, j]; ok &= nr > 20
            A = A[ok]; A = A.mean(0) + (A - A.mean(0)) * np.sqrt(nr[ok] / S[name].sum())[:, None]
        else:
            A = A[ok]
        C = np.cov(A.T); hart = (len(A) - NB - 2) / (len(A) - 1)
        return C / hart, A

    def one(name, verbose=True):
        s = S[name]
        if s.sum() < 30 or name not in names:
            return None
        d0 = np.nanmean(prof[s], 0)
        C, A = covs(name)
        Ci = np.linalg.inv(C); sg = np.sqrt(np.diag(C))
        # three null-mean estimates: the mask-matched randoms on the data map; the simulations' mean at the cluster positions;
        # the rotations' mean (where run)
        nm_rand = rstack; nm_sim = A.mean(0)
        nm_rot = rot['means'][:, rot['names'].index(name)] if rot is not None and name in rot['names'] else None
        if nm_rot is not None:
            nm_rot = np.nanmean(nm_rot, 0)
        dd = d0 - nm_rand  # the estimator: random-corrected
        lnM = np.log(M[s]); zz = z[s]
        nodes, w = posterior_nodes(lnM, sig[s])
        m1, m2 = stacked_model_posterior(hm, grid, nodes, w, zz, f_in, f_out)
        m1_med, m2_med = stacked_model_posterior(hm, grid, lnM[:, None], np.ones(1), zz, f_in, f_out)
        # the PDF check
        pw = pdfw[s].copy(); pg = np.log(np.where(pgrid[s] > 0, pgrid[s], 1e13))
        pw[bad[s]] = 0; pw[bad[s], pg.shape[1] // 2] = 1
        m1_pdf, m2_pdf = stacked_model_posterior(hm, grid, pg, pw, zz, f_in, f_out)
        # fits: (a) A_1h, A_2h free (the statistic); (b) one amplitude on 1h+2h; (c) A_1h with the 2h fixed at sigma8 0.88/0.81/0.76
        p3, c3, x3 = r121.gls(dd, C, np.vstack([m1, m2]).T)
        p1, c1, x1 = r121.gls(dd, C, (m1 + m2)[:, None])
        fixed = {s8: r121.gls(dd - m2 * (s8 / 0.88) ** 2, C, m1[:, None])[0][0] for s8 in (0.88, 0.81, 0.76)}
        p3m, c3m, _ = r121.gls(dd, C, np.vstack([m1_med, m2_med]).T)
        p3p, c3p, _ = r121.gls(dd, C, np.vstack([m1_pdf, m2_pdf]).T)
        p3raw, _, _ = r121.gls(d0, C, np.vstack([m1, m2]).T)
        p3sim, _, _ = r121.gls(d0 - nm_sim, C, np.vstack([m1, m2]).T)
        p3rot = r121.gls(d0 - nm_rot, C, np.vstack([m1, m2]).T)[0] if nm_rot is not None else None
        x0 = dd @ Ci @ dd
        # the null distribution of the statistic from the simulations (random-corrected the same way with each sim's randoms)
        X = np.vstack([m1, m2]).T; Wt = np.linalg.inv(X.T @ Ci @ X) @ X.T @ Ci
        Anull = (A - simrand[:len(A)]) @ Wt.T
        # the rotations' null of the statistic: each rotation keeps a different number of clusters in the footprint, so its
        # amplitude is rescaled by sqrt(n_rot / n) to the data sample's size before the spread is taken
        Arot = None
        if nm_rot is not None:
            j = rot['names'].index(name)
            Rm = rot['means'][:, j]; nr = rot_ns[:, j]; okr = np.all(np.isfinite(Rm), 1) & (nr > 20)
            Arot = (Rm[okr] @ Wt.T) * np.sqrt(nr[okr] / s.sum())[:, None]
        # b_M with the two-halo amplitude profiled out, and with the two-halo fixed (information)
        chi = np.zeros(len(bgrid)); a2s = np.zeros(len(bgrid)); chif = np.zeros(len(bgrid))
        for i, bm in enumerate(bgrid):
            if name not in FULL:
                chi[i] = chif[i] = (bm - 1) ** 2 * 1e3; continue
            mm1, mm2 = stacked_model_posterior(hm, grid, nodes + np.log(bm), w, zz, f_in, f_out)
            chi[i], a2s[i] = profile_chi2(dd, Ci, mm1, mm2)
            r = dd - mm1 - mm2; chif[i] = r @ Ci @ r
        b68 = interval(bgrid, chi, 1.0); b95 = interval(bgrid, chi, 3.84)
        bf68 = interval(bgrid, chif, 1.0); bf95 = interval(bgrid, chif, 3.84)
        chi_at = lambda bm: np.interp(bm, bgrid, chi) - chi.min()
        res = dict(name=name, n=int(s.sum()), d=dd, draw=d0, sig=sg, m1=m1, m2=m2, A1=p3[0], A2=p3[1], sA1=np.sqrt(c3[0, 0]), sA2=np.sqrt(c3[1, 1]),
                   corr=c3[0, 1] / np.sqrt(c3[0, 0] * c3[1, 1]), chi2_3=x3, A=p1[0], sA=np.sqrt(c1[0, 0]), chi2_A=x1, chi2_null=x0,
                   A1_fixed=fixed, A1_med=p3m[0], A1_pdf=p3p[0], A2_med=p3m[1], A2_pdf=p3p[1], A1_raw=p3raw[0], A1_simcorr=p3sim[0],
                   A1_rotcorr=None if p3rot is None else p3rot[0], bM=b68[0], b68=b68[1:], b95=b95[1:], bMf=bf68[0], bf68=bf68[1:], bf95=bf95[1:],
                   dchi_05=chi_at(0.5), dchi_08=chi_at(0.8), dchi_1=chi_at(1.0), dchif_05=np.interp(0.5, bgrid, chif) - chif.min(),
                   Anull_mean=Anull.mean(0), Anull_std=Anull.std(0), nsim=len(A), Arot_mean=None if Arot is None else Arot.mean(0),
                   Arot_std=None if Arot is None else Arot.std(0), nrot=0 if Arot is None else len(Arot),
                   p_emp=(Anull[:, 0] - Anull[:, 0].mean() >= p3[0] - Anull[:, 0].mean()).mean(),
                   nm_rand=nm_rand, nm_sim=nm_sim, nm_rot=nm_rot, medM=np.median(M[s]), medz=np.median(zz), sig_med=np.median(sig[s]),
                   ratio_post_med_1h=m1.sum() / m1_med.sum(), ratio_pdf_med_1h=m1_pdf.sum() / m1_med.sum())
        if verbose:
            print(f'\n[{name}] n={s.sum()}  median M500 {res["medM"]:.2e}  z {res["medz"]:.3f}  sigma_lnM {res["sig_med"]:.2f}  ({len(A)} sims)')
            print('  theta  | data-rand x1e4 | raw x1e4 | sigma | 1h | 2h   (model: posterior-averaged)')
            for i in range(NB):
                print(f'  {mid[i]:6.1f} | {dd[i] * 1e4:8.3f} | {d0[i] * 1e4:8.3f} | {sg[i] * 1e4:6.3f} | {m1[i] * 1e4:7.3f} | {m2[i] * 1e4:7.3f}')
            print(f'  null-mean x1e4: randoms {np.round(nm_rand * 1e4, 3)}\n                  sims    {np.round(nm_sim * 1e4, 3)}' +
                  ('' if nm_rot is None else f'\n                  rotations {np.round(nm_rot * 1e4, 3)}'))
            print(f'  STATISTIC  A_1h = {p3[0]:.3f} +- {np.sqrt(c3[0, 0]):.3f}, A_2h = {p3[1]:.3f} +- {np.sqrt(c3[1, 1]):.3f} (corr {res["corr"]:.2f}; '
                  f'chi2 {x3:.1f}/{NB - 2}; null chi2 {x0:.1f}/{NB})')
            print(f'    b_M (2h profiled) = {b68[0]:.3f}  68 % [{b68[1]:.3f}, {b68[2]:.3f}]  95 % [{b95[1]:.3f}, {b95[2]:.3f}];  '
                  f'dchi2 at b_M = 0.5: {res["dchi_05"]:.1f}, 0.8: {res["dchi_08"]:.1f}, 1.0: {res["dchi_1"]:.1f}')
            print(f'    detection: A_1h / sigma {p3[0] / np.sqrt(c3[0, 0]):.1f}; single A(1h+2h) {p1[0]:.3f} +- {np.sqrt(c1[0, 0]):.3f} '
                  f'({p1[0] / np.sqrt(c1[0, 0]):.1f} sigma); sim null of A_1h {Anull[:, 0].mean():+.3f} +- {Anull[:, 0].std():.3f} '
                  f'(empirical p {res["p_emp"]:.4f}, {len(A)} sims)' + ('' if Arot is None else f'; rotation null {Arot[:, 0].mean():+.3f} +- {Arot[:, 0].std():.3f} ({len(Arot)})'))
            print(f'    information: A_1h with 2h fixed at sigma8 0.88 / 0.81 / 0.76: {fixed[0.88]:.3f} / {fixed[0.81]:.3f} / {fixed[0.76]:.3f}; '
                  f'b_M with 2h fixed {bf68[0]:.3f} [{bf68[1]:.3f}, {bf68[2]:.3f}] 95 % [{bf95[1]:.3f}, {bf95[2]:.3f}] (dchi2 at 0.5: {res["dchif_05"]:.1f})')
            print(f'    model variants: medians A_1h {p3m[0]:.3f} (1h prediction x{1 / res["ratio_post_med_1h"]:.3f} of the posterior one), '
                  f'M500_PDF A_1h {p3p[0]:.3f} (x{res["ratio_pdf_med_1h"]:.3f}); null-mean: raw {p3raw[0]:.3f}, sim-corrected {p3sim[0]:.3f}' +
                  ('' if p3rot is None else f', rotation-corrected {p3rot[0]:.3f}'))
        return res

    for name in SAMPLES + INFO_SUBSETS:
        r = one(name)
        if r is not None:
            results[name] = r
    np.savez(os.path.join(OUT, f'fit_{args.variant}{args.tag}{gtag}.npz'), results=np.array([results], dtype=object))
    print('\nsummary: sample | n | A_1h +- s | A_2h +- s | b_M 68 % | 95 % | A_1h/s | sim null std | rot null std')
    for name, r in results.items():
        print(f'  {name:16s} | {r["n"]:5d} | {r["A1"]:.3f} +- {r["sA1"]:.3f} | {r["A2"]:.3f} +- {r["sA2"]:.3f} | {r["bM"]:.2f} [{r["b68"][0]:.2f}, {r["b68"][1]:.2f}] | '
              f'[{r["b95"][0]:.2f}, {r["b95"][1]:.2f}] | {r["A1"] / r["sA1"]:.1f} | {r["Anull_std"][0]:.3f} | ' +
              ('-' if r['Arot_std'] is None else f'{r["Arot_std"][0]:.3f}'))


# ----------------------------------------------------------------------------------------------------------- Planck
def planck(args):
    """the Planck PR3 tSZ-deprojected stack on the clusters kept in the ACT footprint (Galactic frame, R-121's band 8-2048,
    the Planck lensing mask), with rotation nulls for its covariance; written as data_planck / rotations_planck files that `fit`
    reads with --variant planck (the sims file is the ACT one only for the subset names; covariance is from the rotations)."""
    import healpy as hp
    from scipy.sparse import save_npz
    pdir = os.path.join(DATA, 'pr3lens', 'COM_Lensing_Szdeproj_4096_R3.00')
    c = r121.load_cat()
    selA = np.load(os.path.join(OUT, 'selection_baseline.npz')); keepA = selA['keep']
    l, b = r121.radec_to_gal(c['RA'], c['DEC'])
    mask = hp.read_map(os.path.join(pdir, 'mask.fits.gz'))
    b_in = hp.ud_grade(mask, NSIDE_IN) >= 0.99; b_out = hp.ud_grade(mask, NSIDE_OUT) >= 0.99
    f30, f120, keepP = footprint(l, b, b_in, b_out)
    keep = keepA & keepP
    log(f'Planck: {keep.sum()} of the {keepA.sum()} ACT-footprint clusters pass the Planck footprint selection')
    S_in, S_out, S30, S120 = stacking_matrices(l[keep], b[keep])
    for nm, Sx in [('S_in', S_in), ('S_out', S_out), ('S30', S30), ('S120', S120)]:
        save_npz(os.path.join(OUT, f'{nm}_planck.npz'), Sx)
    np.savez(os.path.join(OUT, 'selection_planck.npz'), keep=keep, f30=f30, f120=f120, inreg=keepA, lon=l, lat=b)
    np.savez(os.path.join(OUT, 'masks_planck.npz'), b_in=b_in, b_out=b_out)
    W_in, W_out, f120r, pos, tried = random_weights(b_in, b_out, args.nrand, 9, (b[keep].min() - 1, b[keep].max() + 1))
    np.savez(os.path.join(OUT, 'randoms_planck.npz'), W_in=W_in, W_out=W_out, f120=f120r, pos=pos)
    alm = read_alm(os.path.join(pdir, 'TT', 'dat_klm.fits'))
    st = Stacker('planck')
    m_in, m_out = maps_from_alm(alm, 8, 2048)
    prof = st.profiles(m_in, m_out, b_in, b_out); rs = st.random_stack(m_in, m_out)
    np.savez(os.path.join(OUT, 'data_planck.npz'), prof=prof, random=rs, lmin=8, lmax=2048, map_rms=(m_in[b_in].std(), m_out[b_out].std()))
    log(f'Planck data stack x1e4: {np.round(np.nanmean(prof, 0) * 1e4, 3)}; randoms {np.round(rs * 1e4, 3)}')
    # rotations in Galactic longitude for the covariance (and a sims file in the same format from them, for `fit`)
    S = subset_masks(c, keep); names = sorted(S)
    rng = np.random.default_rng(3); means = []; rands = []
    for r in range(args.nrot):
        t0 = time.time()
        dphi = np.radians(rng.uniform(25, 335))
        mi, mo = maps_from_alm(alm, 8, 2048, dphi=dphi)
        bi = rotate_binary(b_in, NSIDE_IN, dphi); bo = rotate_binary(b_out, NSIDE_OUT, dphi)
        kr = st.selection(bi, bo); pr = st.profiles(mi, mo, bi, bo)
        means.append(np.array([np.nanmean(pr[S[n] & kr], 0) if (S[n] & kr).sum() else np.full(NB, np.nan) for n in names]))
        rands.append(np.zeros(NB))
        if (r + 1) % 5 == 0 or r == 0:
            np.savez(os.path.join(OUT, 'sims_planck.npz'), idx=np.arange(len(means)), means=np.array(means), randoms=np.array(rands), names=np.array(names))
            log(f'Planck rotation {r + 1}/{args.nrot}: kept {kr.sum()} ({time.time() - t0:.0f}s)')
    np.savez(os.path.join(OUT, 'sims_planck.npz'), idx=np.arange(len(means)), means=np.array(means), randoms=np.array(rands), names=np.array(names))
    log('planck done')


def transfer(args):
    """the release's recommended normalisation check: the reconstructed simulations against their input convergence fields,
    both filtered to the band and masked with the clean footprint; T(L) = C_L^{out,in} / C_L^{in,in} in bands (pseudo-C_L,
    the same mask on both so its coupling cancels at leading order), averaged over the first `nsims` simulations."""
    import healpy as hp
    b_in, _ = binary_masks(args.variant)
    w = b_in.astype(float)
    tdir = os.path.join(ACT, 'transfer'); os.makedirs(tdir, exist_ok=True)
    edges = np.array([20, 40, 80, 150, 250, 400, 600, 800, 1000, 1300, 1600, 2048])
    rows = []
    for i in range(1, args.nsims + 1):
        fo = os.path.join(tdir, f'out_{i:04d}.fits'); fi = os.path.join(tdir, f'in_{i:04d}.fits')
        for src, dst in [(f'maps/{args.variant}/simulations/kappa_alm_sim_act_dr6_lensing_v1_{args.variant}_{i:04d}.fits', fo),
                         (f'sim_inputs/kappa_alm/input_kappa_alm_sim_{i:04d}.fits', fi)]:
            if os.system(f'sh {os.path.join(ACT, "fetch.sh")} {src} {dst}') != 0:
                log(f'transfer: fetch failed for sim {i}'); continue
        try:
            ao = read_alm(fo); ai = read_alm(fi)
        except Exception as e:  # noqa: BLE001
            log(f'transfer: sim {i} unreadable ({e})'); continue
        lmo = hp.Alm.getlmax(len(ao)); lmi = hp.Alm.getlmax(len(ai))
        mo = hp.alm2map(hp.almxfl(ao, band_filter(np.arange(lmo + 1), args.lmin, args.lmax)), NSIDE_IN, lmax=lmo) * w
        mi = hp.alm2map(hp.almxfl(ai, band_filter(np.arange(lmi + 1), args.lmin, args.lmax)), NSIDE_IN, lmax=lmi) * w
        cx = hp.anafast(mo, mi, lmax=args.lmax); cii = hp.anafast(mi, lmax=args.lmax); coo = hp.anafast(mo, lmax=args.lmax)
        L = np.arange(args.lmax + 1)
        T = []; R = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            s = (L >= lo) & (L < hi)
            T.append(cx[s].sum() / cii[s].sum()); R.append(coo[s].sum() / cii[s].sum())
        rows.append((T, R))
        os.remove(fo); os.remove(fi)
        log(f'transfer sim {i}: T(L) by band {np.round(T, 3)}')
        np.savez(os.path.join(OUT, f'transfer_{args.variant}.npz'), T=np.array([r[0] for r in rows]), R=np.array([r[1] for r in rows]), edges=edges)
    T = np.array([r[0] for r in rows])
    log(f'transfer function, {len(T)} sims, bands {edges.tolist()}:\n  mean {np.round(T.mean(0), 3)}\n  error of the mean {np.round(T.std(0) / np.sqrt(len(T)), 3)}')


def grid(args):
    """compute and cache the halo-model grid for the given model options (so the variants can run in parallel)."""
    hm, gtag = make_model(args)
    model_grid_here(hm, gtag)
    log(f'grid {gtag or "baseline"} done')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['prep', 'data', 'sims', 'rotations', 'fit', 'planck', 'grid', 'transfer'])
    ap.add_argument('--variant', default='baseline')
    ap.add_argument('--simvariant', default='', help='use another variant\'s simulations for the covariance (default: own)')
    ap.add_argument('--lmin', type=int, default=LMIN)
    ap.add_argument('--lmax', type=int, default=LMAX)
    ap.add_argument('--nrand', type=int, default=10000)
    ap.add_argument('--nrot', type=int, default=200)
    ap.add_argument('--nsims', type=int, default=400)
    ap.add_argument('--idle_max', type=int, default=90, help='sims: give up after this many 10 s waits with no new file')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--tag', default='')
    ap.add_argument('--cscale', type=float, default=1.0)
    ap.add_argument('--exclusion', action='store_true')
    ap.add_argument('--nonlinear', action='store_true', help='two-halo term on the halofit P(k)')
    ap.add_argument('--cov', default='sims', choices=['sims', 'rotations'], help='fit: covariance source')
    ap.add_argument('--ex_r200', type=float, default=1.5)
    args = ap.parse_args()
    {'prep': prep, 'data': data, 'sims': sims, 'rotations': rotations, 'fit': fit, 'planck': planck, 'grid': grid,
     'transfer': transfer}[args.cmd](args)


if __name__ == '__main__':
    main()
