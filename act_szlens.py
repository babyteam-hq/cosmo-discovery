"""CMB-lensing masses of the ACT DR6 SZ clusters in redshift bins (R-154, cosmo-deep-22).

Does the SZ mass bias keep falling at z > 1? M_lens / M_SZ of ACT DR6 SZ-selected clusters per redshift bin, from the ACT
DR6 lensing convergence stacked with R-122's machinery (`erass1_act.py`: the masks, the annulus stack, the halo model filtered
as the map, the 400 simulations, the rotated-map nulls), at each cluster's catalogue SZ mass.

Catalogue: ACT DR6 SZ clusters (ACT-DES-HSC Collaboration 2025, arXiv:2507.21459; LAMBDA `DR6_cluster-catalog_v1.0.fits`).

usage (python -I; data under data/, gitignored):
  python -I act_szlens.py counts                  # clusters per z bin in the DR6 lensing footprint and in the SZ inpainting holes
  python -I act_szlens.py forecast                # Fisher sigma(M_lens/M_SZ) per bin, calibrated on R-122's realised eRASS1 error
"""
import argparse, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import erass1_act as ea      # noqa: E402
import erass1_cmblens as r121  # noqa: E402

DATA = os.path.join(HERE, 'data')
OUT = os.path.join(DATA, 'act_szlens')
CAT = os.path.join(DATA, 'actcl', 'DR6_cluster-catalog_v1.0.fits')
ZBINS = [(0.0, 0.5), (0.5, 0.8), (0.8, 1.0), (1.0, 2.1)]
# the ACT cluster catalogue's fiducial cosmology (flat, H0 = 70, Om = 0.3), in which its M500c are computed
H0_ACT, OM_ACT, OB_ACT = 70.0, 0.30, 0.049
log = ea.log
ea.OUT = OUT   # the masks, stacking matrices and randoms of erass1_act's helpers live in this run's directory


def load_act(mass='M500c'):
    from astropy.io import fits
    d = fits.open(CAT)[1].data
    c = {k: np.asarray(d[k], dtype=float) for k in ['RADeg', 'decDeg', 'SNR', 'fixed_SNR', 'redshift', 'redshiftErr',
                                                    'M500c', 'M500c_errPlus', 'M500c_errMinus', 'M500cUncorr',
                                                    'A10_M500c', 'A10_M500c_errPlus', 'A10_M500c_errMinus', 'A10_M500cUncorr',
                                                    'theta500Arcmin', 'Q', 'flags']}
    for k in list(c):
        if k.startswith('M500') or k.startswith('A10_M500'):
            c[k] = c[k] * 1e14
    c['M'] = c[mass]
    c['RA'], c['DEC'], c['z'] = c['RADeg'], c['decDeg'], c['redshift']
    c['legacy'] = np.asarray(d['footprint_Legacy'], bool) if 'footprint_Legacy' in d.columns.names else np.zeros(len(d), bool)
    c['spec'] = np.asarray(d['redshiftType']) == 'spec'
    return c


class ACTHaloModel(r121.HaloModel):
    """R-121's NFW + linear two-halo model, in the ACT catalogue's cosmology and to z = 2.2 (sigma8 0.81 for the two-halo)."""

    def __init__(self, lmax=3000, cscale=1.0, s8=0.81):
        import camb
        from astropy.cosmology import FlatLambdaCDM
        from scipy.interpolate import RectBivariateSpline
        self.cscale = cscale
        self.co = FlatLambdaCDM(H0=H0_ACT, Om0=OM_ACT, Ob0=OB_ACT, Tcmb0=2.7255)
        self.lmax = lmax
        pars = camb.CAMBparams()
        h = H0_ACT / 100
        pars.set_cosmology(H0=H0_ACT, ombh2=OB_ACT * h ** 2, omch2=(OM_ACT - OB_ACT) * h ** 2, mnu=0.06, tau=0.054)
        pars.InitPower.set_params(As=2.1e-9, ns=0.965)
        pars.set_matter_power(redshifts=list(np.linspace(0, 2.4, 25))[::-1], kmax=50.0)
        pars.NonLinear = camb.model.NonLinear_none
        res = camb.get_results(pars)
        self.s8_scale = (s8 / res.get_sigma8_0()) ** 2
        kh, zs, pk = res.get_matter_power_spectrum(minkh=1e-4, maxkh=50, npoints=600)
        self.h = h
        self.k = kh * h
        self.zs = np.array(zs)
        self.pk = np.array(pk) / h ** 3 * self.s8_scale
        order = np.argsort(self.zs); self.zs = self.zs[order]; self.pk = self.pk[order]
        self._pk_spl = RectBivariateSpline(self.zs, np.log(self.k), np.log(self.pk))
        self.rho_m0 = (self.co.Om0 * self.co.critical_density0).to('Msun/Mpc3').value
        self.rho_c0 = self.co.critical_density0.to('Msun/Mpc3').value
        self.theta = np.logspace(np.log10(0.05 * ea.ARCMIN), np.log10(6 * np.pi / 180), 900)
        self.ell = np.arange(0, lmax + 1)
        self._cache = {}


def model_grid(hm, tag, nM=20, nz=22):
    """kappa(ell) 1h and 2h on a (ln M500c, z) grid to z = 2.1, cached."""
    os.makedirs(OUT, exist_ok=True)
    fn = os.path.join(OUT, f'modelgrid{tag}.npz')
    if os.path.exists(fn):
        g = np.load(fn)
        return g['lnM'], g['z'], g['k1'], g['k2'], g['bias']
    lnM = np.linspace(np.log(3e13), np.log(3e15), nM)
    zg = np.concatenate([np.linspace(0.02, 0.1, 3), np.linspace(0.15, 2.1, nz - 3)])
    k1 = np.zeros((nM, len(zg), len(hm.ell))); k2 = np.zeros_like(k1); bias = np.zeros((nM, len(zg)))
    t0 = time.time()
    for i, lm in enumerate(lnM):
        for j, z in enumerate(zg):
            a, b, info = hm.cluster_kell(np.exp(lm), z)
            k1[i, j] = a; k2[i, j] = b; bias[i, j] = info['bias']
        log(f'  grid row {i + 1}/{nM} ({time.time() - t0:.0f}s)')
    np.savez(fn, lnM=lnM, z=zg, k1=k1, k2=k2, bias=bias)
    return lnM, zg, k1, k2, bias


def annulus_windows(ell, filt, edges):
    """W[a, l]: the annulus-mean response to a unit kappa(l), times the map filter (annulus_means_from_ell's kernel)."""
    from scipy.special import j1
    e = edges * ea.ARCMIN; l = ell.astype(float)
    W = np.zeros((len(e) - 1, len(l)))
    for i in range(len(e) - 1):
        a, b = e[i], e[i + 1]
        with np.errstate(invalid='ignore', divide='ignore'):
            w = 2 * (b * j1(l * b) - a * j1(l * a)) / (l * (b ** 2 - a ** 2))
        w[l == 0] = 1.0
        W[i] = w * filt
    return W


def single_cov(ell, ctot, W):
    """covariance of one object's annulus means for a Gaussian field of power ctot (flat sky)."""
    dl = np.gradient(ell.astype(float))
    return (W * (ell * dl * ctot / (2 * np.pi))[None, :]) @ W.T


def selection(c, variant='baseline'):
    """the footprint selection of R-122 (clean within 30', >= 80 % clean within 2 deg) and the SZ-hole flag."""
    import healpy as hp
    os.makedirs(OUT, exist_ok=True)
    fn = os.path.join(OUT, f'selection_{variant}.npz')
    if os.path.exists(fn):
        return dict(np.load(fn))
    old = ea.OUT; ea.OUT = OUT
    try:
        b_in, b_out = ea.binary_masks(variant)
    finally:
        ea.OUT = old
    f30, f120, keep = ea.footprint(c['RA'], c['DEC'], b_in, b_out)
    szfn = os.path.join(ea.ACT, 'sz_cluster_mask_radius_6.0_arcmin_nside_8192.fits')
    sz = hp.read_map(szfn, dtype=np.float32)
    insz = sz[hp.ang2pix(8192, c['RA'], c['DEC'], lonlat=True)] < 0.5
    del sz
    np.savez(fn, keep=keep, f30=f30, f120=f120, insz=insz)
    return dict(np.load(fn))


def counts(args):
    c = load_act()
    s = selection(c)
    keep, insz = s['keep'], s['insz']
    log(f'catalogue {len(c["z"])}; in the DR6 lensing footprint (R-122 selection) {keep.sum()}; of these in an SZ inpainting hole '
        f'{(keep & insz).sum()}')
    for lo, hi in ZBINS + [(1.2, 2.1)]:
        b = (c['z'] >= lo) & (c['z'] < hi)
        k = b & keep
        log(f'  z {lo}-{hi}: catalogue {b.sum()}, kept {k.sum()}, kept & in hole {(k & insz).sum()}, kept outside holes {(k & ~insz).sum()}; '
            f'median M500c {np.median(c["M500c"][k]):.2e}, A10 {np.median(c["A10_M500c"][k]):.2e}, median SNR {np.median(c["fixed_SNR"][k]):.1f}; '
            f'in hole by SNR>5: {(k & insz & (c["fixed_SNR"] > 5)).sum()}/{(k & (c["fixed_SNR"] > 5)).sum()}')


def fisher_sigma(K1, K2, Ci, W):
    """per-sample marginal sigma of A_1h (A_2h free) and of A_1h alone, clusters independent: F = sum_i M_i^T C^-1 M_i."""
    lw = np.arange(W.shape[1]) * np.gradient(np.arange(W.shape[1], dtype=float)) / (2 * np.pi)
    m1 = (K1 * lw) @ W.T; m2 = (K2 * lw) @ W.T       # (n, nb): the annulus means, as annulus_means_from_ell
    F = np.array([[np.einsum('ia,ab,ib->', m1, Ci, m1), np.einsum('ia,ab,ib->', m1, Ci, m2)],
                  [np.einsum('ia,ab,ib->', m2, Ci, m1), np.einsum('ia,ab,ib->', m2, Ci, m2)]])
    return np.sqrt(np.linalg.inv(F)[0, 0]), 1 / np.sqrt(F[0, 0])


def forecast(args):
    """Fisher sigma(A_1h) per z bin from the release N_L (plus the lensing signal power) and the annulus windows of R-122,
    clusters independent; calibrated by the same Fisher on R-122's eRASS1 cosmology sample against its realised 400-sim
    error (0.110). sigma(M_lens/M_SZ) ~ sigma(A_1h) / (dlnA/dlnM) with dlnA/dlnM from the model."""
    import healpy as hp
    nl = np.loadtxt(os.path.join(ea.ACT, 'N_L_kk_act_dr6_lensing_v1_baseline.txt'))
    log(f'N_L file shape {nl.shape}; first rows {nl[:3]}')
    hm = ACTHaloModel(lmax=3000)
    ell = hm.ell
    L = nl[:, 0] if nl.ndim == 2 and nl.shape[1] > 1 else np.arange(len(nl))
    N = nl[:, 1] if nl.ndim == 2 and nl.shape[1] > 1 else nl
    NL = np.interp(ell, L, N)
    # the lensing signal power (CAMB, the same cosmology) adds to the per-object variance
    import camb
    pars = camb.set_params(H0=H0_ACT, ombh2=OB_ACT * 0.49, omch2=(OM_ACT - OB_ACT) * 0.49, As=2.1e-9, ns=0.965, lmax=3500, lens_potential_accuracy=1)
    res = camb.get_results(pars)
    cpp = res.get_lens_potential_cls(lmax=3000)[:, 0]   # [L(L+1)]^2 C_L^pp / 2pi
    ckk = np.zeros(len(ell)); Lc = np.arange(len(cpp)); ckk[1:] = (cpp[1:] * 2 * np.pi / 4)[:len(ell) - 1] if len(cpp) >= len(ell) else 0
    ctot = NL + ckk
    f_in, f_out = ea.filters(hm, ea.LMIN, ea.LMAX)
    W = np.vstack([annulus_windows(ell, f_in, ea.EDGES_IN), annulus_windows(ell, f_out, ea.EDGES_OUT)])
    C = single_cov(ell, ctot, W)
    Ci = np.linalg.inv(C)
    grid = model_grid(hm, '')
    lnMg, zg, k1g, k2g, _ = grid

    def K(M, z, bM=1.0):
        lm = np.log(bM * M)
        return r121.interp_grid(lnMg, zg, k1g, lm, z), r121.interp_grid(lnMg, zg, k2g, lm, z)

    # calibration on R-122: the eRASS1 cosmology sample in the footprint (1,907 kept), realised sigma_A1h 0.110
    cal = None
    sel_e = os.path.join(DATA, 'erass1_act', 'selection_baseline.npz')
    try:
        ce = r121.load_cat()
        se = selection({'RA': ce['RA'], 'DEC': ce['DEC']}, variant='baseline_erass1') if not os.path.exists(sel_e) else dict(np.load(sel_e))
        ke = se['keep'] & ce['cosmo']
        # eRASS1 masses are in the eROSITA cosmology; at the 2 % level the model difference is irrelevant to a calibration
        K1e, K2e = K(ce['M500'][ke], ce['BEST_Z'][ke])
        s2, s1 = fisher_sigma(K1e, K2e, Ci, W)
        cal = 0.110 / s2
        log(f'calibration: eRASS1 cosmology kept {ke.sum()}, Fisher sigma_A1h (2h free) {s2:.4f} (2h fixed {s1:.4f}); realised 0.110 '
            f'-> factor {cal:.2f}')
    except Exception as e:  # noqa: BLE001
        log('calibration not computed:', repr(e))
    c = load_act()
    s = selection(c)
    keep, insz = s['keep'], s['insz']
    rows = []
    for mass in ['M500c', 'A10_M500c']:
        for lo, hi in ZBINS + [(1.2, 2.1)]:
          for hole in ('all', 'nosz', 'insz'):
            b = keep & (c['z'] >= lo) & (c['z'] < hi)
            b = b & ~insz if hole == 'nosz' else (b & insz if hole == 'insz' else b)
            M, z = c[mass][b], c['z'][b]
            K1, K2 = K(M, z)
            s2, s1 = fisher_sigma(K1, K2, Ci, W)
            # dlnA/dlnM: the one-halo template's matched amplitude at 1.1 M against M
            K1p, _ = K(M, z, 1.1)
            lw = ell * 1.0 / (2 * np.pi)
            m1 = ((K1 * lw) @ W.T).sum(0); m1p = ((K1p * lw) @ W.T).sum(0)
            slope = np.log((m1 @ Ci @ m1p) / (m1 @ Ci @ m1)) / np.log(1.1)
            sA = s2 * (cal or 1.0)
            rows.append((mass, hole, lo, hi, b.sum(), sA, slope, sA / slope))
            log(f'{mass} [{hole}] z {lo}-{hi}: n {b.sum()} (in holes {(b & insz).sum()}), Fisher sigma_A1h {s2:.3f} (2h fixed {s1:.3f}), '
                f'calibrated {sA:.3f}; dlnA/dlnM {slope:.2f}; sigma(M_lens/M_SZ)/(M_lens/M_SZ) {sA / slope:.3f}')
    np.save(os.path.join(OUT, 'forecast_rows.npy'), np.array(rows, dtype=object), allow_pickle=True)


def prep(args):
    """stacking matrices for the kept clusters and the mask-matched randoms (erass1_act's, on this catalogue)."""
    from scipy.sparse import save_npz
    c = load_act(); s = selection(c); keep = s['keep']
    tag = 'baseline'
    b_in, b_out = ea.binary_masks('baseline')
    log(f'stacking matrices for {keep.sum()} kept clusters')
    S_in, S_out, S30, S120 = ea.stacking_matrices(c['RA'][keep], c['DEC'][keep])
    for nm, M in [('S_in', S_in), ('S_out', S_out), ('S30', S30), ('S120', S120)]:
        save_npz(os.path.join(OUT, f'{nm}_{tag}.npz'), M)
    dec_range = (c['DEC'][keep].min() - 1, c['DEC'][keep].max() + 1)
    W_in, W_out, f120r, pos, tried = ea.random_weights(b_in, b_out, args.nrand, 7, dec_range)
    log(f'randoms: {len(pos)} of {tried} tried; mean f120 randoms {f120r.mean():.3f} vs clusters {s["f120"][keep].mean():.3f}')
    np.savez(os.path.join(OUT, f'randoms_{tag}.npz'), W_in=W_in, W_out=W_out, f120=f120r, pos=pos, dec_range=dec_range)
    log('prep done')


def sims(args):
    """per-cluster annulus profiles on each streamed simulation (data/act_dr6/sims, act_dr6_simstream.sh): a (nsim, n, NB) float32
    memmap, so any subset's covariance (and the joint covariance of the redshift bins) comes from the same 400 maps; each
    simulation file deleted once stacked."""
    import glob
    c = load_act(); s = selection(c); keep = s['keep']
    b_in, b_out = ea.binary_masks('baseline')
    st = ea.Stacker('baseline')
    n = st.n
    fn = os.path.join(OUT, 'sims_prof_baseline.npy'); fr = os.path.join(OUT, 'sims_rand_baseline.npy')
    fd = os.path.join(OUT, 'sims_done_baseline.npy')
    P = np.lib.format.open_memmap(fn, mode='r+' if os.path.exists(fn) else 'w+', dtype=np.float32, shape=(400, n, ea.NB))
    Rr = np.lib.format.open_memmap(fr, mode='r+' if os.path.exists(fr) else 'w+', dtype=np.float64, shape=(400, ea.NB))
    done = np.load(fd) if os.path.exists(fd) else np.zeros(400, bool)
    simdir = os.path.join(ea.ACT, 'sims'); slog = os.path.join(ea.ACT, 'simstream_baseline.log')
    idle = 0
    while done.sum() < 400 and idle < 60:
        files = sorted(glob.glob(os.path.join(simdir, 'kappa_alm_sim_act_dr6_lensing_v1_baseline_*.fits')))
        files = [f for f in files if not done[int(f[-9:-5]) - 1]]
        if not files:
            if os.path.exists(slog) and open(slog).read().count('ALLDONE') >= args.streams:
                break
            time.sleep(10); idle += 1; continue
        idle = 0
        f = files[0]; i = int(f[-9:-5]); t0 = time.time()
        try:
            alm = ea.read_alm(f)
        except Exception as e:  # noqa: BLE001
            log(f'sim {i}: unreadable ({e}); removing'); os.remove(f); continue
        m_in, m_out = ea.maps_from_alm(alm)
        P[i - 1] = st.profiles(m_in, m_out, b_in, b_out).astype(np.float32)
        Rr[i - 1] = st.random_stack(m_in, m_out)
        P.flush(); Rr.flush(); done[i - 1] = True; np.save(fd, done)
        os.remove(f); open(os.path.join(simdir, f'done_baseline_{i:04d}'), 'w').close()
        log(f'sim {i}: {done.sum()} done ({time.time() - t0:.0f}s); all-cluster stack x1e4 {np.round(np.nanmean(P[i - 1], 0) * 1e4, 3)}')
    log(f'sims: {done.sum()} stacked')


def planck_cib(args):
    """the measured CIB amplitude per redshift bin (method check point 3): aperture flux of the Planck PR3 545 and 857 GHz maps
    at the kept clusters (0-8' minus the 15-25' ring), the same at 20,000 random positions in the footprint (null), the A10 tSZ
    increment subtracted at 545 (negligible at 857); per bin the mean flux in Jy, extrapolated to 90 and 150 GHz with a modified
    blackbody (beta, T_dust) grid in the rest frame."""
    import healpy as hp
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = load_act(); s = selection(c); keep = s['keep']
    co = ACTHaloModel.__new__(ACTHaloModel)
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=H0_ACT, Om0=OM_ACT, Ob0=OB_ACT)
    g = SkyCoord(c['RA'] * u.deg, c['DEC'] * u.deg).galactic
    l, b = g.l.deg, g.b.deg
    R = np.load(os.path.join(OUT, 'randoms_baseline.npz'))['pos'][:20000]
    gr = SkyCoord(R[:, 0] * u.deg, R[:, 1] * u.deg).galactic
    out = {}
    for nu in (545, 857):
        m = hp.read_map(os.path.join(DATA, 'planck', f'HFI_SkyMap_{nu}_2048_R3.01_full.fits'), field=0, dtype=np.float64)  # MJy/sr
        nside = hp.get_nside(m); om = hp.nside2pixarea(nside)
        def ap(ll, bb):
            v = hp.ang2vec(ll, bb, lonlat=True); f = np.zeros(len(ll))
            for i in range(len(ll)):
                p0 = hp.query_disc(nside, v[i], 8 * ea.ARCMIN); p1 = hp.query_disc(nside, v[i], 25 * ea.ARCMIN)
                p15 = hp.query_disc(nside, v[i], 15 * ea.ARCMIN); ring = np.setdiff1d(p1, p15)
                f[i] = (m[p0].sum() - len(p0) * np.median(m[ring])) * om * 1e6      # Jy
            return f
        fc = np.full(len(l), np.nan); fc[keep] = ap(l[keep], b[keep])
        fr = ap(gr.l.deg, gr.b.deg)
        out[nu] = (fc, fr)
        log(f'{nu} GHz: randoms mean {np.mean(fr) * 1e3:.2f} ± {np.std(fr) / np.sqrt(len(fr)) * 1e3:.2f} mJy (sd {np.std(fr) * 1e3:.0f})')
    # the A10 tSZ at 545 (cylinder within 5 R500 = 1.796 x the sphere), in Jy
    x = 545 / 56.78; hx = x ** 4 * np.exp(x) / np.expm1(x) ** 2 * (x / np.tanh(x / 2) - 4); I0 = 270.33e6   # Jy/sr per unit y
    Ez = cosmo.efunc(c['z']); DA = cosmo.angular_diameter_distance(c['z']).value
    Y = 10 ** -4.739 * Ez ** (2 / 3) * (c['A10_M500c'] / 3e14) ** 1.79 * 1.796 / DA ** 2       # sr
    tsz545 = I0 * hx * Y
    from scipy.constants import h as hP, k as kB, c as cL
    def mbb(nu_ghz, z, beta, T):
        nr = nu_ghz * 1e9 * (1 + z)
        return nr ** (3 + beta) / np.expm1(hP * nr / (kB * T))
    def dBdT(nu_ghz):   # Jy/sr per K_CMB
        x = nu_ghz / 56.78; return 2 * kB * (nu_ghz * 1e9) ** 2 / cL ** 2 * x ** 2 * np.exp(x) / np.expm1(x) ** 2 * 1e26
    rows = []
    for lo, hi in ZBINS:
        bsel = keep & (c['z'] >= lo) & (c['z'] < hi)
        f857 = out[857][0][bsel] - np.mean(out[857][1]); f545 = out[545][0][bsel] - np.mean(out[545][1]) - tsz545[bsel]
        e857 = np.std(out[857][1]) / np.sqrt(bsel.sum()); e545 = np.std(out[545][1]) / np.sqrt(bsel.sum())
        zm = np.mean(c['z'][bsel])
        line = f'z {lo}-{hi} (n {bsel.sum()}, mean z {zm:.2f}): 857 {f857.mean() * 1e3:.1f} ± {e857 * 1e3:.1f} mJy, 545 (tSZ {tsz545[bsel].mean() * 1e3:.1f} mJy removed) {f545.mean() * 1e3:.1f} ± {e545 * 1e3:.1f} mJy, ratio 545/857 {f545.mean() / f857.mean():.2f}'
        ext = {}
        for beta in (1.5, 2.0):
            for T in (20.0, 25.0, 35.0):
                r545 = mbb(545, zm, beta, T) / mbb(857, zm, beta, T)
                for nu in (150, 90):
                    S = f857.mean() * mbb(nu, zm, beta, T) / mbb(857, zm, beta, T)     # Jy
                    ext[(beta, T, nu)] = (S, S / dBdT(nu) * 1e6, r545)          # Jy, uK sr
        line += '; 150 GHz flux (beta, T: mJy, model 545/857): ' + ', '.join(f'({k[0]},{k[1]:.0f}) {v[0] * 1e3:.2f} [{v[2]:.2f}]' for k, v in ext.items() if k[2] == 150)
        log(line)
        rows.append(dict(lo=lo, hi=hi, z=zm, n=int(bsel.sum()), f857=f857.mean(), e857=e857, f545=f545.mean(), e545=e545,
                         ext={f'{k[0]}_{k[1]:.0f}_{k[2]}': v for k, v in ext.items()}))
    np.save(os.path.join(OUT, 'cib_rows.npy'), np.array(rows, dtype=object), allow_pickle=True)


def data(args):
    """per-cluster annulus profiles of a map variant at the kept clusters, and the random stack (data_<variant>.npz)."""
    b_in, b_out = ea.binary_masks('baseline')
    st = ea.Stacker('baseline')
    for v in args.variants:
        m_in, m_out = ea.maps_from_alm(ea.read_alm(ea.act_file('alm', v) if v != 'curl' else
                                                   os.path.join(ea.ACT, 'kappa_alm_data_act_dr6_lensing_v1_curl.fits')))
        prof = st.profiles(m_in, m_out, b_in, b_out)
        rs = st.random_stack(m_in, m_out)
        np.savez(os.path.join(OUT, f'data_{v}.npz'), prof=prof.astype(np.float32), random=rs)
        log(f'{v}: stacked {st.n} clusters; random stack x1e4 {np.round(rs * 1e4, 2)}')


def rotations(args):
    """rotated-map nulls (map and mask turned together in RA by 25-335 deg, clusters fixed, selection re-applied): per-cluster
    profiles per rotation (NaN where a cluster fails the rotated selection), for the covariance of the variants without
    simulations and the rotation null of the fitted statistics."""
    b_in, b_out = ea.binary_masks('baseline')
    st = ea.Stacker('baseline')
    rng = np.random.default_rng(args.seed)
    for v in args.variants:
        fn_alm = ea.act_file('alm', v) if v != 'curl' else os.path.join(ea.ACT, 'kappa_alm_data_act_dr6_lensing_v1_curl.fits')
        alm = ea.read_alm(fn_alm)
        P = np.full((args.nrot, st.n, ea.NB), np.nan, np.float32)
        for r in range(args.nrot):
            dphi = np.radians(rng.uniform(25, 335))
            m_in, m_out = ea.maps_from_alm(alm, dphi=dphi)
            bi = ea.rotate_binary(b_in, ea.NSIDE_IN, dphi); bo = ea.rotate_binary(b_out, ea.NSIDE_OUT, dphi)
            kr = st.selection(bi, bo)
            pr = st.profiles(m_in, m_out, bi, bo); pr[~kr] = np.nan
            P[r] = pr
            if (r + 1) % 10 == 0:
                log(f'{v}: rotation {r + 1}/{args.nrot}, kept {kr.sum()}')
                np.save(os.path.join(OUT, f'rot_{v}.npy'), P[:r + 1])
        np.save(os.path.join(OUT, f'rot_{v}.npy'), P)


TINKER08 = np.array([[200, 0.186, 1.47, 2.57, 1.19], [300, 0.200, 1.52, 2.25, 1.27], [400, 0.212, 1.56, 2.05, 1.34],
                     [600, 0.218, 1.61, 1.87, 1.45], [800, 0.248, 1.87, 1.59, 1.58], [1200, 0.255, 2.13, 1.51, 1.80],
                     [1600, 0.260, 2.30, 1.46, 1.97], [2400, 0.260, 2.53, 1.44, 2.24], [3200, 0.260, 2.66, 1.41, 2.44]])


def tinker08_dndlnm(hm, lnM, z, s8_scale=1.0):
    """Tinker et al. 2008 dn/dlnM for M500c (Delta_m = 500 / Omega_m(z)), relative units; s8_scale multiplies sigma(M)."""
    Om_z = hm.co.Om(z); Dm = 500 / Om_z
    lD = np.log10(TINKER08[:, 0])
    A0, a0, b0, c0 = [np.interp(np.log10(Dm), lD, TINKER08[:, k]) for k in range(1, 5)]
    alpha = 10 ** (-(0.75 / np.log10(Dm / 75)) ** 1.2)
    A = A0 * (1 + z) ** -0.14; a = a0 * (1 + z) ** -0.06; b = b0 * (1 + z) ** -alpha; c = c0
    M = np.exp(lnM)
    # sigma(M) uses the Lagrangian radius of the mass; M500c -> M_Delta_m is the same halo, so M is used directly
    sig = np.array([hm.sigma_m(m, z) for m in M]) * s8_scale
    f = A * ((sig / b) ** -a + 1) * np.exp(-c / sig ** 2)
    dlns = -np.gradient(np.log(sig), lnM)
    return f * hm.rho_m0 / M * dlns


def selection_mock(args):
    """the SZ selection (Eddington) bias of the fitted R, by mock (method check point 5, and points 4 and 6):
    - halos from Tinker 2008 at the bin's catalogue redshifts; the true A10-scale observable ln y = s ln((1-b) M) + noise-free
      intrinsic scatter sigma_int, with s the catalogue's own y-M slope in the bin (from fixed_y_c against A10_M500cUncorr);
    - each mock halo gets the noise level of a random catalogue cluster in the bin (k_i = q_i / M_unc,i^s, so q = k M_y^s);
      Gaussian noise in y; selected at q > 4;
    - optionally, optical confirmation at z > 1 as a logistic in the true mass (M_50, slope 4);
    - the catalogue's posterior: A10 with b = 0, scatter 0.2, Tinker prior at sigma8 0.80 (the catalogue's), its median as M_SZ;
    - the stack: the mock's data vector is the mean annulus model at the true masses (1h + 2h at the truth), fitted for R with the
      real fit (posterior nodes, two-halo free, the per-cluster covariance over N); B = R_fit / R_true per bin.
    Grid: 1-b in {0.25, 0.4, 0.6, 0.75, 0.9}; sigma_int (truth) 0.1, 0.2, 0.3; sigma8 (truth) 0.77, 0.80, 0.83."""
    from astropy.io import fits
    rng = np.random.default_rng(3)
    c = load_act(); s = selection(c); keep = s['keep']
    d = fits.open(CAT)[1].data
    yc = np.asarray(d['fixed_y_c'], float)
    hm = ACTHaloModel(lmax=3000)
    grid = model_grid(hm, '')
    f_in, f_out = ea.filters(hm, ea.LMIN, ea.LMAX)
    ag = ea.annulus_grid(hm, grid, f_in, f_out)
    nl = np.loadtxt(os.path.join(ea.ACT, 'N_L_kk_act_dr6_lensing_v1_baseline.txt'))
    W = np.vstack([annulus_windows(hm.ell, f_in, ea.EDGES_IN), annulus_windows(hm.ell, f_out, ea.EDGES_OUT)])
    C1 = single_cov(hm.ell, np.interp(hm.ell, nl[:, 0], nl[:, 1]), W)
    lnMg = np.linspace(np.log(3e13), np.log(3e15), 160)
    rows = []
    Rg = np.exp(np.linspace(np.log(0.3), np.log(6.0), 241))
    for lo, hi in ZBINS:
        b = keep & (c['z'] >= lo) & (c['z'] < hi)
        zc = c['z'][b]; Mu = c['A10_M500cUncorr'][b]; q = c['fixed_SNR'][b]
        sl = np.linalg.lstsq(np.vstack([np.ones(b.sum()), np.log(Mu), zc]).T, np.log(yc[b]), rcond=None)[0][1]
        kq = q / Mu ** sl
        Ci = np.linalg.inv(C1 / b.sum())
        zgrid = np.quantile(zc, [0.1, 0.3, 0.5, 0.7, 0.9])
        # the catalogue posterior on a (z, ln y_obs/sigma) grid: prior x likelihood with b = 0, scatter 0.2, Tinker sigma8 0.80
        prior = {zz: tinker08_dndlnm(hm, lnMg, zz, s8_scale=0.80 / 0.81) for zz in zgrid}
        sim_prior = {}
        for case in args.cases:
            onemb, sint, s8t, m50 = case
            Mt_all, Mp_all, zz_all, Sg_all, Q_all = [], [], [], [], []
            for zz in zgrid:
                key = (zz, s8t)
                if key not in sim_prior:
                    sim_prior[key] = tinker08_dndlnm(hm, lnMg, zz, s8_scale=s8t / 0.81)
                pm = sim_prior[key] * np.exp(lnMg); pm = pm / pm.sum()                     # dn/dlnM dlnM, sampled in lnM
                n = args.nmock
                lnMt = rng.choice(lnMg, size=n, p=pm) + rng.uniform(-0.5, 0.5, n) * (lnMg[1] - lnMg[0])
                k = rng.choice(kq, size=n)
                lnMy = lnMt + np.log(onemb) + rng.normal(0, sint, n)          # the A10-scale mass the y corresponds to
                qtrue = k * np.exp(sl * lnMy)
                qobs = qtrue + rng.normal(0, 1, n)
                sel = qobs > 4
                if zz > 1.0 and m50 > 0:
                    sel &= rng.uniform(size=n) < 1 / (1 + (m50 / np.exp(lnMt)) ** 4)
                # the catalogue's posterior median for each selected halo (b = 0, scatter 0.2, its prior)
                lnMt_s, qo, ks = lnMt[sel], qobs[sel], k[sel]
                X = np.linspace(-3 * 0.2, 3 * 0.2, 13); wX = np.exp(-0.5 * (X / 0.2) ** 2); wX /= wX.sum()
                lnMp = np.empty(len(qo)); sgp = np.empty(len(qo))
                for j0 in range(0, len(qo), 2000):
                    sl_ = slice(j0, j0 + 2000)
                    qm = ks[sl_, None, None] * np.exp(sl * (lnMg[None, :, None] + X[None, None, :]))
                    like = (np.exp(-0.5 * (qo[sl_, None, None] - qm) ** 2) * wX).sum(2)
                    cdf = np.cumsum(like * prior[zz][None, :], 1); cdf /= cdf[:, -1:]
                    lnMp[sl_] = [np.interp(0.5, cc, lnMg) for cc in cdf]
                    sgp[sl_] = [0.5 * (np.interp(0.84, cc, lnMg) - np.interp(0.16, cc, lnMg)) for cc in cdf]
                Sg_all.append(sgp)
                Mt_all.append(lnMt_s); Mp_all.append(lnMp); zz_all.append(np.full(sel.sum(), zz)); Q_all.append(qo)
            lnMt_s = np.concatenate(Mt_all); lnMp = np.concatenate(Mp_all); zs = np.concatenate(zz_all); sgp = np.concatenate(Sg_all)
            nodes, wn = ea.posterior_nodes(lnMp, sgp)     # the real fit's posterior averaging (5 Gauss-Hermite nodes)
            # the stack: data at the truth, model from the posterior medians scaled by R; R by the real fit (2h free)
            m1t, m2t = ea.stacked_model_posterior(hm, ag, lnMt_s[:, None], np.ones(1), zs)
            dvec = m1t + m2t
            chi = []
            for R in Rg:
                m1, m2 = ea.stacked_model_posterior(hm, ag, nodes + np.log(R), wn, zs)
                chi.append(ea.profile_chi2(dvec, Ci, m1, m2)[0])
            Rfit = Rg[np.argmin(chi)]
            Rtrue_lens = 1 / onemb
            qall = np.concatenate(Q_all); sub = {}
            if args.qsplit:
                for lab, qs in [('q4-5', qall <= 5), ('q>5', qall > 5)]:
                    m1t, m2t = ea.stacked_model_posterior(hm, ag, lnMt_s[qs][:, None], np.ones(1), zs[qs])
                    ch = [ea.profile_chi2(m1t + m2t, Ci, *ea.stacked_model_posterior(hm, ag, nodes[qs] + np.log(R), wn, zs[qs]))[0] for R in Rg]
                    sub[lab] = Rg[np.argmin(ch)] / Rtrue_lens
                log(f'   q split: B(q 4-5) {sub["q4-5"]:.3f}, B(q > 5) {sub["q>5"]:.3f}')
            ratio_med = np.exp(np.median(lnMp - lnMt_s))
            rows.append(dict(lo=lo, hi=hi, slope=sl, case=case, nsel=len(lnMp), Rfit=Rfit, Rtrue=Rtrue_lens, B=Rfit / Rtrue_lens,
                             post_over_true_med=ratio_med))
            log(f'z {lo}-{hi} (y-M slope {sl:.2f}): 1-b {onemb}, sigma_int {sint}, sigma8 {s8t}, M50 {m50:.1e}: selected {len(lnMp)}, '
                f'median M_post/M_true {ratio_med:.3f}; R_fit {Rfit:.3f} vs R_true {Rtrue_lens:.3f}: B = {Rfit / Rtrue_lens:.3f}')
        np.save(os.path.join(OUT, f'selmock_rows{args.tag}.npy'), np.array(rows, dtype=object), allow_pickle=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('cmd')
    p.add_argument('--nrand', type=int, default=100000)
    p.add_argument('--streams', type=int, default=1)   # sims: how many simstreams write ALLDONE to the log
    p.add_argument('--nmock', type=int, default=40000)
    p.add_argument('--tag', default='')
    p.add_argument('--grid', default='main')
    p.add_argument('--qsplit', action='store_true')
    p.add_argument('--variants', nargs='+', default=['baseline'])
    p.add_argument('--nrot', type=int, default=60)
    p.add_argument('--seed', type=int, default=1)
    a = p.parse_args()
    if a.grid == 'main':
        a.cases = [(ob, 0.2, 0.80, 0) for ob in (0.25, 0.4, 0.6, 0.75, 0.9)] + \
                  [(0.75, si, 0.80, 0) for si in (0.1, 0.3)] + [(0.75, 0.2, s8, 0) for s8 in (0.77, 0.83)] + \
                  [(0.4, 0.3, 0.80, 0), (0.4, 0.2, 0.80, 1.0e14), (0.75, 0.2, 0.80, 1.0e14), (0.75, 0.2, 0.80, 1.5e14), (0.4, 0.2, 0.80, 1.5e14)]
    elif a.grid == 'qsplit':
        a.cases = [(0.75, 0.2, 0.80, 0), (0.4, 0.2, 0.80, 0), (0.75, 0.3, 0.80, 0)]
    else:
        a.cases = [(0.75, 0.2, 0.80, 0)]
    {'counts': counts, 'forecast': forecast, 'prep': prep, 'sims': sims, 'cib': planck_cib, 'selmock': selection_mock, 'data': data, 'rotations': rotations}[a.cmd](a)


if __name__ == '__main__':
    main()
