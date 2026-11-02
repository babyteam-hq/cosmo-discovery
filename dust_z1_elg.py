"""R-152 (cosmo-deep-19b): the ELG limit of R-150 on its own, with the referee's two problems measured for the ELG sample.

R-150 (results-R150.md) found rho = 0.04 +- 0.29 for the reddening of z > 2 quasars behind DESI DR1 ELGs at 0.8 < z < 1.6,
against Menard et al. 2010's z ~ 0.36 profile at the same proper r_p. Its gate asked whether that limit is free of the
stellar (photometric) systematic and of the magnification leak; this script measures both for the ELGs alone and puts the
comparison at matched comoving scales. Nothing in the estimator changes from R-150's M1 (dust_z1_common.py,
dust_z1_stack.py); R-152's additions are labelled in results-R150.md (rows E0-E6).

usage: python -u dust_z1_elg.py prep|prepstars|main|stars|table
  prep   slim ELG LSS + QSO_full (dust_z1_prep.py); prepstars: the DR1 stars (zpix-main-dark + zpix-main-backup, SPECTYPE STAR)
  main   the R-150 ELG stack again (reproduction), plus per-band brightness and the per-quasar colour-magnitude slope;
         saves the pair-weight matrices
  stars  the stellar null around the ELGs with R-150's 58k QSO-target stars and with all DR1 dark+backup PSF stars,
         per colour and through the quasar estimator's own projection
  table  every R-152 number (logs/r152_table.log)
"""
import os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dust_z1_common as C
import dust_z1_stack as S
from dust_z1_stack import jk_cov, fit_powerlaw, HARTLAP, NJK, REF_RP, NB

ELG_BINS = C.ZBINS['ELG']
MEN_ERR = 0.25               # Menard's amplitude error as arXiv:2503.04098 quotes it (R-150 used their 7 %)
Z_MEN = 0.36


# ---------------------------------------------------------------------------------------------------------------- prep
def prep():
    import fitsio
    import dust_z1_prep as P
    P.slim_lss('ELG_LOPnotqso')
    P.slim_qso()


def prep_stars():
    import fitsio
    # Galactic transmission for the zpix stars (zpix has EBV but no MW_TRANSMISSION): per band and imaging system the
    # coefficient R = -2.5 log10(T)/EBV, read off QSO_full itself.
    q = fitsio.read(os.path.join(C.D, 'qso', 'QSO_full.dat.fits'),
                    columns=['EBV', 'PHOTSYS', 'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z',
                             'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2'], rows=np.arange(0, 2_000_000, 7))
    Rcoef = {}
    for ps in ('N', 'S'):
        s = (q['PHOTSYS'] == ps) & (q['EBV'] > 0.01)
        for b in ('G', 'R', 'Z', 'W1', 'W2'):
            Rcoef[(ps, b)] = float(np.median(-2.5 * np.log10(q[f'MW_TRANSMISSION_{b}'][s]) / q['EBV'][s]))
    print('extinction coefficients from QSO_full:', {f'{k[0]}{k[1]}': round(v, 3) for k, v in Rcoef.items()})
    cols = ['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'Z', 'ZWARN', 'SPECTYPE', 'SUBTYPE', 'MORPHTYPE', 'EBV', 'PHOTSYS',
            'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FLUX_IVAR_W1', 'ZCAT_PRIMARY', 'DESI_TARGET']
    parts = []
    for prog in ('dark', 'backup'):
        fn = os.path.join(C.D, 'stars', f'zpix-main-{prog}.fits')
        f = fitsio.FITS(fn)
        n = f[1].get_nrows()
        for i0 in range(0, n, 2_000_000):
            a = f[1].read(columns=cols, rows=np.arange(i0, min(n, i0 + 2_000_000)))
            a = a[(a['SPECTYPE'] == 'STAR') & (a['ZWARN'] == 0) & a['ZCAT_PRIMARY'] & (a['MORPHTYPE'] == 'PSF')]
            parts.append((prog, a))
        print(f'{prog}: {n} spectra, {sum(len(p[1]) for p in parts if p[0] == prog)} PSF stars (ZWARN 0, primary)', flush=True)
    a = np.concatenate([p[1] for p in parts])
    prog = np.concatenate([np.full(len(p[1]), p[0] == 'dark') for p in parts])
    out = dict(targetid=a['TARGETID'], ra=a['TARGET_RA'], dec=a['TARGET_DEC'], ebv=a['EBV'], north=a['PHOTSYS'] == 'N',
               flux_ivar_w1=a['FLUX_IVAR_W1'], dark=prog, desi_target=a['DESI_TARGET'],
               subtype=np.char.strip(a['SUBTYPE'].astype(str)))
    for b in ('G', 'R', 'Z', 'W1', 'W2'):
        Rb = np.where(out['north'], Rcoef[('N', b)], Rcoef[('S', b)])
        out[f'flux_{b.lower()}'] = a[f'FLUX_{b}']
        out[f'mw_transmission_{b.lower()}'] = 10 ** (-0.4 * Rb * a['EBV'])
    _, u = np.unique(out['targetid'], return_index=True)          # a star observed in both programs counts once
    out = {k: v[u] for k, v in out.items()}
    np.savez(os.path.join(C.D, 'stars_dr1.npz'), **out)
    print(f'stars_dr1: {len(u)} unique PSF stars')


# ---------------------------------------------------------------------------------------------------------------- main
def load_elg_bin(f, r, perm, zlo, zhi):
    sel_f = (f['z'] >= zlo) & (f['z'] < zhi)
    ridx = perm[(r['z'][perm] >= zlo) & (r['z'][perm] < zhi)][:S.RAN_FACTOR * sel_f.sum()]
    return S.sub(r, ridx)


def run_main():
    t0 = time.time()
    q, res, qr, reg = S.background()
    nq = len(q['z'])
    mres = S.mag_residuals(q)
    f, r, perm = S.load_fg('ELG')
    for (zlo, zhi) in ELG_BINS:
        zf = 0.5 * (zlo + zhi)
        rr = load_elg_bin(f, r, perm, zlo, zhi)
        Sd, rpm = S.pair_sums(q['ra'], q['dec'], q['z'], f['ra'], f['dec'], f['z'], f['w'], zlo, zhi, nq)
        Sr, _ = S.pair_sums(q['ra'], q['dec'], q['z'], rr['ra'], rr['dec'], rr['z'], rr['w'], zlo, zhi, nq)
        # pair-weighted mean foreground redshift per r_p bin (for the comoving conversion)
        Sz, _ = S.pair_sums(q['ra'], q['dec'], q['z'], f['ra'], f['dec'], f['z'], f['w'] * f['z'], zlo, zhi, nq)
        out = dict(rp=rpm, npair=(Sd > 0).sum(0), zmean=Sz.sum(0) / Sd.sum(0))
        ah, w = C.av_estimator(res, qr, zf, 'SMC')
        full, jk, _ = S.stat(Sd, Sr, ah, w, reg)
        out['A_SMC'] = full / S.RESP(zf); out['jk_SMC'] = jk / S.RESP(zf)
        R = C.response(zf, 'SMC')
        for i, cname in enumerate(C.COLOURS):
            full, jk, _ = S.stat(Sd, Sr, res[:, i] / R[i], np.ones(nq), reg)
            out[f'A_{cname}'] = full / S.RESP(zf); out[f'jk_{cname}'] = jk / S.RESP(zf)
        for i, bname in enumerate(('g', 'r', 'z', 'W1')):
            full, jk, _ = S.stat(Sd, Sr, mres[:, i], np.ones(nq), reg)
            out[f'dm_{bname}'] = full; out[f'jkdm_{bname}'] = jk
        # colour-magnitude slope of A_hat inside each r-quartile (the leak's lever arm), and the pair-and-estimator
        # weighted mean slope over the quasars that enter this bin's pairs
        slopes = np.zeros(4)
        for k in range(4):
            s = (qr == k) & np.isfinite(ah)
            rm = q['mag']['r'][s]
            slopes[k] = np.polyfit(rm - rm.mean(), ah[s], 1)[0]
        wp = np.where(np.isfinite(ah), w, 0.)[:, None] * Sd                      # (nq, NB) weights of the data pairs
        out['cm_slope'] = slopes
        out['slope_pairs'] = np.array([np.sum(wp[:, b] * slopes[qr]) / np.sum(wp[:, b]) for b in range(NB)])
        # the same statistic per jackknife region for the stacked slope (it is a fixed property of the sample)
        out['wsum_q'] = np.array([wp[qr == k].sum() for k in range(4)])
        np.savez(os.path.join(C.D, f'r152_ELG_{zlo}_{zhi}.npz'), **out)
        np.savez_compressed(os.path.join(C.D, f'r152_S_ELG_{zlo}_{zhi}.npz'), Sd=Sd.astype(np.float32),
                            Sr=Sr.astype(np.float32))
        print(f'ELG {zlo}-{zhi} done ({time.time() - t0:.0f} s): pairs {int(out["npair"].sum())}; '
              f'cm slopes {np.round(slopes, 3)}', flush=True)


def run_nowise():
    """The WISE-free estimator (GLS over g-r and r-z only), from the saved pair-weight matrices: the stars whose null fails
    entered the sample through W1-based target cuts and fail in z-W1 (results-R150.md E3), so the limit is also given
    without W1."""
    q, res, qr, reg = S.background()
    for (zlo, zhi) in ELG_BINS:
        zf = 0.5 * (zlo + zhi)
        M = np.load(os.path.join(C.D, f'r152_S_ELG_{zlo}_{zhi}.npz'))
        Sd, Sr = M['Sd'].astype(float), M['Sr'].astype(float)
        R = C.response(zf, 'SMC')[:2]
        ah = np.full(len(res), np.nan); w = np.zeros(len(res))
        ok = np.all(np.isfinite(res[:, :2]), axis=1)
        for k in range(4):
            sel = ok & (qr == k)
            Ci = np.linalg.inv(np.cov(res[sel][:, :2].T))
            ah[sel] = res[sel][:, :2] @ (Ci @ R / (R @ Ci @ R)); w[sel] = R @ Ci @ R
        full, jk, _ = S.stat(Sd, Sr, ah, w, reg)
        fn = os.path.join(C.D, f'r152_ELG_{zlo}_{zhi}.npz')
        out = dict(np.load(fn)); out['A_noW1'] = full / S.RESP(zf); out['jk_noW1'] = jk / S.RESP(zf)
        np.savez(fn, **out)
        print(f'noW1 {zlo}-{zhi} done', flush=True)


# --------------------------------------------------------------------------------------------------------------- stars
def star_sample(kind):
    """kind 'qso_targets': R-150's null-4 sample (QSO targets with a stellar spectrum, from QSO_full);
    'dr1': every DR1 dark+backup PSF star with ZWARN 0 (zpix). Same cuts as the quasars on photometry (16.5 < r < 23,
    W1 > 3 sigma)."""
    if kind == 'qso_targets':
        import fitsio
        cols = ['RA', 'DEC', 'SPECTYPE', 'ZWARN', 'Z_RR', 'DESI_TARGET', 'PHOTSYS', 'FLUX_G', 'FLUX_R', 'FLUX_Z',
                'FLUX_W1', 'FLUX_W2', 'FLUX_IVAR_W1', 'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z',
                'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2']
        a = fitsio.read(os.path.join(C.D, 'qso', 'QSO_full.dat.fits'), columns=cols)
        a = a[(a['SPECTYPE'] == 'STAR') & (a['ZWARN'] == 0) & (np.abs(a['Z_RR']) < 0.002) & ((a['DESI_TARGET'] & 4) > 0)]
        st = {k.lower(): a[k] for k in cols if k not in ('SPECTYPE', 'PHOTSYS')}
        st['north'] = a['PHOTSYS'] == 'N'
        st['subtype'] = np.full(len(a), 'X')
    else:
        st = dict(np.load(os.path.join(C.D, 'stars_dr1.npz')))
        # dr1_qso / dr1_elg / dr1_lrg / dr1_other: the DR1 stars by the dark-program target class that put a fibre on them
        # (DESI_TARGET bits 2, 1, 0): tests whether the null follows the sample's make-up rather than the photometry
        cls = {'dr1_qso': (st['desi_target'] & 4) > 0, 'dr1_elg': ((st['desi_target'] & 2) > 0) & ((st['desi_target'] & 5) == 0),
               'dr1_lrg': ((st['desi_target'] & 1) > 0) & ((st['desi_target'] & 6) == 0),
               'dr1_other': (st['desi_target'] & 7) == 0}
        if kind in cls:
            st = {k: v[cls[kind]] for k, v in st.items()}
    st['mag'] = C.mags(st)
    good = np.ones(len(st['ra']), bool)
    for b in ('g', 'r', 'z', 'W1'):
        good &= np.isfinite(st['mag'][b])
    good &= (st['flux_ivar_w1'] > 0) & (st['flux_w1'] * np.sqrt(np.clip(st['flux_ivar_w1'], 0, None)) > 3)
    good &= (st['mag']['r'] > 16.5) & (st['mag']['r'] < 23)
    st = {k: (v[good] if k != 'mag' else {b: m[good] for b, m in v.items()}) for k, v in st.items()}
    # colour cells: imaging system x r quartile x (Redrock stellar subtype for the DR1 sample: spectral type, which is
    # independent of the photometry, so it takes out the stellar locus without touching any colour being tested)
    sub = st['subtype']
    us, inv = np.unique(sub, return_inverse=True)
    st['z'] = 5.0 + inv.astype(float)               # one 'redshift' cell per subtype (zstep 1 below)
    return st


def run_stars():
    t0 = time.time()
    fn = os.path.join(C.D, 'r152_stars.npz')
    old = np.load(fn) if os.path.exists(fn) else None
    out = {} if old is None else {k.rsplit('___', 1)[0]: tuple(old[f"{k.rsplit('___', 1)[0]}___{i}"] for i in range(3))
                                  for k in old.files}
    kinds = sys.argv[2].split(',') if len(sys.argv) > 2 else ['qso_targets', 'dr1']
    q, res, qr, reg = S.background()
    f, r, perm = S.load_fg('ELG', seed=3)
    for kind in kinds:
        st = star_sample(kind)
        ress, qrs = C.colour_residuals(st, zstep=1.0)
        regs = S.jk_regions(st['ra'], st['dec'])
        ns = len(st['ra'])
        rmag = st['mag']['r']
        print(f'{kind}: {ns} stars, median r {np.median(rmag):.2f}, residual robust sigma '
              f'{np.round(1.4826 * np.nanmedian(np.abs(ress), axis=0), 3)}', flush=True)
        for (zlo, zhi) in [(0.8, 1.6)]:
            zf = 1.2
            sel = (f['z'] >= zlo) & (f['z'] < zhi)
            rr = S.sub(r, perm[(r['z'][perm] >= zlo) & (r['z'][perm] < zhi)][:S.RAN_FACTOR * sel.sum()])
            Sd, rpm = S.pair_sums(st['ra'], st['dec'], st['z'], f['ra'], f['dec'], f['z'], f['w'], zlo, zhi, ns)
            Sr, _ = S.pair_sums(st['ra'], st['dec'], st['z'], rr['ra'], rr['dec'], rr['z'], rr['w'], zlo, zhi, ns)
            tag = f'{kind}_{zlo}_{zhi}'
            out[f'{tag}__n'] = (np.array([ns]), np.array([np.median(rmag)]), np.zeros(1))
            # (a) R-150's form: the stars through their own GLS estimator
            ah, w = C.av_estimator(ress, qrs, zf)
            full, jk, _ = S.stat(Sd, Sr, ah, w, regs)
            out[f'{tag}__own'] = (rpm, full, jk)
            # (b) per colour, equal weights: the colour shift at ELG positions
            for i, cname in enumerate(C.COLOURS):
                full, jk, _ = S.stat(Sd, Sr, ress[:, i], np.ones(ns), regs)
                out[f'{tag}__d{cname}'] = (rpm, full, jk)
            # (c) faint stars alone (r > 20, the quasars' magnitudes), own estimator
            full, jk, _ = S.stat(Sd, Sr, np.where(rmag > 20, ah, np.nan), w, regs)
            out[f'{tag}__faint'] = (rpm, full, jk)
            print(f'  {tag} done ({time.time() - t0:.0f} s)', flush=True)
    # the quasar estimator's own projection vectors (per r-quartile and bin), to carry a photometric colour shift seen in
    # the stars into the quasars' A_hat units: A_sys(r_p) = sum_c f_c dcolour_c(r_p), f the weight-averaged projection
    for (zlo, zhi) in [(0.8, 1.6)] + list(ELG_BINS):
        zf = 1.2 if (zlo, zhi) == (0.8, 1.6) else 0.5 * (zlo + zhi)
        Rv = C.response(zf, 'SMC')
        ok = np.all(np.isfinite(res), axis=1)
        fs, ws = [], []
        for k in range(4):
            s = ok & (qr == k)
            Ci = np.linalg.inv(np.cov(res[s].T))
            fs.append(Ci @ Rv / (Rv @ Ci @ Rv)); ws.append((Rv @ Ci @ Rv) * s.sum())
        out[f'proj_{zlo}_{zhi}'] = (np.average(np.array(fs), axis=0, weights=ws), Rv, np.zeros(1))
    np.savez(os.path.join(C.D, 'r152_stars.npz'), **{f'{k}___{i}': v[i] for k, v in out.items() for i in range(3)})


# --------------------------------------------------------------------------------------------------------------- table
def growth(z):
    """Linear growth D(z)/D(0), flat LCDM (Planck 2018 Om), Heath integral."""
    from scipy.integrate import quad
    om = C.cosmo().Om0
    E = lambda a: np.sqrt(om / a ** 3 + 1 - om)
    g = lambda a: 2.5 * om * E(a) * quad(lambda x: 1 / (x * E(x)) ** 3, 0, a)[0]
    return g(1 / (1 + z)) / g(1.0)


def fixed_amp(rp, a, jk, gamma=0.84):
    cov = jk_cov(jk)
    g = np.isfinite(a) & np.all(np.isfinite(jk), axis=0)
    f1 = fit_powerlaw(rp[g], a[g], cov[np.ix_(g, g)], gamma=gamma)
    t = (rp[g] / REF_RP) ** -gamma
    ci = np.linalg.inv(cov[np.ix_(g, g)]) * HARTLAP
    wv = ci @ t / (t @ ci @ t)                       # the amplitude is linear in the 8 bins: A = wv . a
    return f1, wv, g, jk[:, g] @ wv


def chi2_zero(a, jk):
    from scipy.stats import chi2
    cov = jk_cov(jk)
    g = np.isfinite(a) & np.all(np.isfinite(jk), axis=0)
    ci = np.linalg.inv(cov[np.ix_(g, g)]) * HARTLAP
    x = float(a[g] @ ci @ a[g])
    return x, int(g.sum()), float(chi2.sf(x, g.sum()))


def combine(vals, sigs, jks):
    J = np.array(jks); d = J - J.mean(1, keepdims=True)
    cov = (NJK - 1) / NJK * d @ d.T
    sd = np.sqrt(np.diag(cov)); corr = cov / np.outer(sd, sd)
    cov = corr * np.outer(sigs, sigs)
    ci = np.linalg.inv(cov); one = np.ones(len(vals))
    v = 1 / (one @ ci @ one)
    wts = v * (one @ ci)
    return float(wts @ np.asarray(vals)), float(np.sqrt(v)), wts


def bayes_upper(mu, sig, cl=0.95):
    from scipy.stats import norm
    p0 = norm.cdf(0, mu, sig)
    return float(norm.ppf(p0 + cl * (1 - p0), mu, sig))


def table():
    ref = C.menard_amplitude('SMC', REF_RP)
    print(f'Menard SMC reference at {REF_RP:.0f} proper kpc (z = 0.36): {ref:.3e}; Menard amplitude error now {MEN_ERR:.0%}')
    bins = {}
    for (zlo, zhi) in ELG_BINS:
        bins[(zlo, zhi)] = dict(np.load(os.path.join(C.D, f'r152_ELG_{zlo}_{zhi}.npz')))
    # --- E1 reproduction and the main limit (proper and comoving) ---
    print('\n=== E1: R-150 ELG stack, reproduced (SMC GLS, randoms subtracted, / targeting response) ===')
    rows = {}
    for k, s in bins.items():
        f1, wv, g, ajk = fixed_amp(s['rp'], s['A_SMC'], s['jk_SMC'])
        x0, n0, p0 = chi2_zero(s['A_SMC'], s['jk_SMC'])
        zeff = float(np.sum(wv * s['zmean'][g]) / np.sum(wv)) if np.sum(wv) != 0 else 0.5 * sum(k)
        zeff = float(np.average(s['zmean'][g], weights=np.abs(wv)))
        # comoving: our A at comoving R = r_p (1 + z) against Menard's at the same comoving R = r_p,M (1 + 0.36)
        cfac = ((1 + zeff) / (1 + Z_MEN)) ** 0.84
        # information share of each r_p bin in the fixed-shape amplitude (Fisher weight t_i (C^-1 t)_i)
        t = (s['rp'][g] / REF_RP) ** -0.84
        cov = jk_cov(s['jk_SMC'])[np.ix_(g, g)]
        ci = np.linalg.inv(cov) * HARTLAP
        info = t * (ci @ t); info = info / info.sum()
        rows[k] = dict(A=f1['A'], sA=f1['sA'], ajk=ajk, zeff=zeff, cfac=cfac, info=info, rp=s['rp'][g],
                       err=np.sqrt(np.diag(cov)), npair=s['npair'])
        print(f'ELG {k[0]}-{k[1]}: A_150 = {f1["A"] * 1e3:.2f} +- {f1["sA"] * 1e3:.2f} e-3, rho_proper = {f1["A"] / ref:.2f} '
              f'+- {f1["sA"] / ref:.2f} (stat), chi2 vs 0 {x0:.1f}/{n0} (p {p0:.2f}); z_eff {zeff:.3f}, comoving factor '
              f'{cfac:.3f}; information share r_p<131 / 131-458 / >458 kpc: '
              f'{info[s["rp"][g] < 131].sum():.2f} / {info[(s["rp"][g] >= 131) & (s["rp"][g] < 458)].sum():.2f} / '
              f'{info[s["rp"][g] >= 458].sum():.2f}')
    ks = list(rows)

    def comb(vals_per_bin, sig_per_bin, jk_per_bin):
        return combine(vals_per_bin, sig_per_bin, jk_per_bin)
    rp_ = [rows[k]['A'] / ref for k in ks]; sp_ = [rows[k]['sA'] / ref for k in ks]; jp_ = [rows[k]['ajk'] / ref for k in ks]
    m, sm, wts = comb(rp_, sp_, jp_)
    tot = lambda m, sm: np.hypot(sm, abs(m) * MEN_ERR)
    print(f'combined, proper r_p: rho = {m:.3f} +- {sm:.3f} (stat); with Menard 7 %: +- {np.hypot(sm, m * 0.1 / 1.4):.3f}; '
          f'with 25 %: +- {tot(m, sm):.3f}; 95 % upper (1.645 sigma) {m + 1.645 * tot(m, sm):.2f}')
    # R-150's own combination (cross-bin correlations from an unweighted least-squares amplitude per jackknife sample)
    from dust_z1_table import rho_jk, combine as combine150
    j150 = [rho_jk(bins[k]['rp'], bins[k]['jk_SMC'], ref) for k in ks]
    m150, s150, _ = combine150(rp_, sp_, j150)
    print(f'R-150 combination reproduced: rho = {m150:.3f} +- {np.hypot(s150, m150 * 0.1 / 1.4):.3f}; here the cross-bin '
          f'correlations come from the jackknife of the GLS amplitude itself (the estimator that is quoted)')
    rc_ = [rows[k]['A'] / ref * rows[k]['cfac'] for k in ks]; sc_ = [rows[k]['sA'] / ref * rows[k]['cfac'] for k in ks]
    jc_ = [rows[k]['ajk'] / ref * rows[k]['cfac'] for k in ks]
    mc, smc, wc = comb(rc_, sc_, jc_)
    print(f'combined, matched comoving R: rho_c = {mc:.3f} +- {smc:.3f} (stat) +- 25 % ref -> +- {tot(mc, smc):.3f}; '
          f'95 % upper {mc + 1.645 * tot(mc, smc):.2f}; bin weights {np.round(wc, 2)}')
    print(f'  with a physical prior rho_c >= 0 (flat): 95 % upper {bayes_upper(mc, tot(mc, smc)):.2f}; proper '
          f'{bayes_upper(m, tot(m, sm)):.2f}; sensitivity (1.645 sigma about 0): comoving {1.645 * tot(0, smc):.2f}')
    # inner-halo reach
    inner = [(rows[k]['rp'] < 131, rows[k]) for k in ks]
    n131 = sum(int(r_['npair'][:3].sum()) for _, r_ in inner)
    e_in = [r_['err'][m_] for m_, r_ in inner]
    print(f'pairs inside 131 kpc (3 ELG bins): {n131}; per-bin A_V errors there (x1e-3): '
          + '; '.join(np.array2string(e * 1e3, precision=1) for e in e_in))
    for est in ('noW1', 'g-r', 'r-z', 'z-W1'):
        if f'A_{est}' not in bins[ks[0]]:
            continue
        vv, ss, jj, cc = [], [], [], []
        for k in ks:
            f1, wv, g, ajk = fixed_amp(bins[k]['rp'], bins[k][f'A_{est}'], bins[k][f'jk_{est}'])
            vv.append(f1['A'] / ref); ss.append(f1['sA'] / ref); jj.append(ajk / ref); cc.append(rows[k]['cfac'])
        me, se, _ = combine(vv, ss, jj)
        mec, sec, _ = combine([v * c for v, c in zip(vv, cc)], [v * c for v, c in zip(ss, cc)], [v * c for v, c in zip(jj, cc)])
        print(f'estimator {est}: rho_proper = {me:.2f} +- {tot(me, se):.2f}; rho_c = {mec:.2f} +- {tot(mec, sec):.2f}, '
              f'positive-prior 95 % upper {bayes_upper(mec, tot(mec, sec)):.2f}')
    # what the inner bins (r_p < 131 kpc, inside a 10^12 M_sun halo's R_200c at z = 1.2) reach on their own
    ri, si, ji = [], [], []
    for k in ks:
        s = bins[k]
        inn = s['rp'] < 131
        f1, wv, g, ajk = fixed_amp(s['rp'][inn], s['A_SMC'][inn], s['jk_SMC'][:, inn])
        ri.append(f1['A'] / ref); si.append(f1['sA'] / ref); ji.append(ajk / ref)
    mi, smi, _ = combine(ri, si, ji)
    print(f'inner bins alone (20-131 kpc, fixed 0.84 shape): rho_proper = {mi:.2f} +- {smi:.2f}; a 3 sigma detection of '
          f'rho = 1 there needs ~{n131 * (3 * smi) ** 2:.0f} pairs inside 131 kpc (now {n131}); of rho = 0.3, '
          f'~{n131 * (3 * smi / 0.3) ** 2:.1e}')
    # E4: two-halo scaling of the comparison (dust traced like matter, linear bias, same dust bias):
    # A_V,rest(R_c) ~ Omega_d b_d b_g D(z)^2 (1 + z)^2 w_p,lin(R_c)  [proper dust surface density at fixed comoving R]
    rng = np.random.default_rng(152)
    zc = float(np.average([rows[k]['zeff'] for k in ks], weights=wc))
    bD_elg = rng.normal(0.84, 0.84 * 0.15, 200000)               # DESI fiducial b_ELG(z) D(z) = 0.84 (arXiv:1611.00036)
    b_men = rng.normal(1.1, 0.2, 200000)                          # SDSS i < 21 photometric galaxies at z ~ 0.36 (assumed)
    Sfac = ((1 + zc) / (1 + Z_MEN)) ** 2 * bD_elg * growth(zc) / (b_men * growth(Z_MEN) ** 2)
    rho_draw = rng.normal(mc, tot(mc, smc), 200000)
    ratio = rho_draw / Sfac
    print(f'\n=== E4: two-halo scaling at z_eff {zc:.2f} ===')
    print(f'expected rho_c for an unchanged Omega_d b_d: S = {np.median(Sfac):.2f} [{np.percentile(Sfac, 16):.2f}, '
          f'{np.percentile(Sfac, 84):.2f}]; measured / S = Omega_d b_d (z_eff) / (z = 0.36) = {np.median(ratio):.2f} '
          f'[{np.percentile(ratio, 16):.2f}, {np.percentile(ratio, 84):.2f}], 95 % upper {np.percentile(ratio, 95):.2f}; '
          f'positive-prior 95 % upper {np.percentile(ratio[ratio > 0], 95):.2f}; no-evolution excluded at '
          f'{(np.median(Sfac) - mc) / np.hypot(tot(mc, smc), np.std(Sfac)):.1f} sigma')
    # Menard & Fukugita 2012 (arXiv:1204.1978): Omega_dust(Mg II) = 2.3 +- 0.2 e-6 at z ~ 0.5, ~ (1 + z)^-1; Menard
    # et al. 2010: Omega_dust ~ 5e-6 for all galaxy-correlated (halo + IGM) dust at z ~ 0.3 (2.1e-6 inside R_vir)
    om_mg = 2.3e-6 * (1.5 / (1 + zc))
    for name, om0 in (('all galaxy-correlated dust (5e-6)', 5e-6), ('dust inside R_vir (2.1e-6)', 2.1e-6)):
        pred = om_mg / om0 * Sfac
        rr_ = rho_draw / pred
        print(f'Mg II budget at z_eff ({om_mg:.2e}) over Menard {name}: predicted rho_c {np.median(pred):.2f} '
              f'[{np.percentile(pred, 16):.2f}, {np.percentile(pred, 84):.2f}]; measured / predicted {np.median(rr_):.2f} '
              f'[{np.percentile(rr_, 16):.2f}, {np.percentile(rr_, 84):.2f}], positive-prior 95 % upper '
              f'{np.percentile(rr_[rr_ > 0], 95):.2f}')
    # --- E2 magnification for the ELGs: measured achromatic shift and the leak it implies ---
    print('\n=== E2: magnification leak, ELG sample ===')
    leak_rho, leak_sig, leak_jk, mmu_tab = [], [], [], []
    for k in ks:
        s = bins[k]
        zf = 0.5 * sum(k)
        kb = C.band_k(zf); kb = np.array([kb['g'], kb['r'], kb['z'], kb['W1']])
        t = (s['rp'] / REF_RP) ** -0.84
        bands = ('g', 'r', 'z', 'W1')
        y = np.concatenate([s[f'dm_{b}'] for b in bands])
        J = np.concatenate([s[f'jkdm_{b}'] for b in bands], axis=1)
        cov = jk_cov(J)
        ci = np.linalg.pinv(cov) * (NJK - len(y) - 2) / (NJK - 1)
        X = np.c_[np.concatenate([-t] * 4), np.concatenate([kk * t for kk in kb])]
        F = X.T @ ci @ X
        lin = np.linalg.solve(F, X.T @ ci)           # parameters are linear in y
        p = lin @ y
        pj = J @ lin.T
        emu = np.sqrt(jk_cov(pj)[0, 0])
        sbar = float(np.average(s['slope_pairs'], weights=np.abs(rows[k]['info']) if len(rows[k]['info']) == NB else None))
        # leak in A_hat at 150 kpc: s_bar x m_mu(150), m_mu > 0 = brighter; a lensed quasar is intrinsically fainter at
        # its observed r, hence redder by the within-quartile slope
        L, sL = sbar * p[0] / S.RESP(zf), sbar * emu / S.RESP(zf)
        leak_rho.append(L / ref); leak_sig.append(sL / ref); leak_jk.append(sbar * pj[:, 0] / S.RESP(zf) / ref)
        mmu_tab.append((k, p[0], emu, sbar))
        print(f'ELG {k[0]}-{k[1]}: m_mu(150 kpc) = {p[0]:.4f} +- {emu:.4f} mag (brighter > 0); A_V from brightness '
              f'{p[1] * 1e3:.2f} e-3; cm slopes per quartile {np.round(s["cm_slope"], 3)}, pair-weighted mean '
              f'{sbar:.3f} mag/mag; predicted leak in rho {L / ref:+.3f} +- {sL / ref:.3f}')
    ml, sml, _ = combine(leak_rho, leak_sig, leak_jk)
    print(f'combined predicted leak (proper): {ml:+.3f} +- {sml:.3f} in rho; comoving x{np.average([rows[k]["cfac"] for k in ks], weights=wc):.2f}')
    # leak-corrected limit: subtract per bin, jackknife the difference (same regions)
    rcorr = [rp_[i] - leak_rho[i] for i in range(3)]
    jcorr = [jp_[i] - leak_jk[i] for i in range(3)]
    scorr = [np.sqrt(jk_cov(np.c_[jcorr[i]])[0, 0]) for i in range(3)]
    scorr = [np.hypot(sp_[i], leak_sig[i]) if not np.isfinite(scorr[i]) else max(scorr[i], sp_[i]) for i in range(3)]
    mlc, smlc, _ = combine(rcorr, scorr, jcorr)
    print(f'leak subtracted, proper: rho = {mlc:.3f} +- {smlc:.3f} -> 95 % upper {mlc + 1.645 * tot(mlc, smlc):.2f}; '
          f'comoving: {mlc * np.average([rows[k]["cfac"] for k in ks], weights=wc):.3f}')
    # --- E3 stellar null for the ELGs ---
    fn = os.path.join(C.D, 'r152_stars.npz')
    if os.path.exists(fn):
        n = np.load(fn)
        get = lambda key: tuple(n[f'{key}___{i}'] for i in range(3))
        print('\n=== E3: stars as backgrounds around the ELGs ===')
        for kind in ('qso_targets', 'dr1', 'dr1_qso', 'dr1_elg', 'dr1_lrg', 'dr1_other'):
            for (zlo, zhi) in [(0.8, 1.6)]:
                tag = f'{kind}_{zlo}_{zhi}'
                if f'{tag}__own___0' not in n.files:
                    continue
                rpm, a, jk = get(f'{tag}__own')
                f1, wv, g, ajk = fixed_amp(rpm, a, jk)
                x0, n0, p0 = chi2_zero(a, jk)
                rpf, af, jkf = get(f'{tag}__faint')
                ff, _, _, _ = fixed_amp(rpf, af, jkf)
                proj, Rv, _ = get(f'proj_{zlo}_{zhi}')
                dc = [get(f'{tag}__d{c}') for c in C.COLOURS]
                parts = []
                for i, c in enumerate(C.COLOURS):
                    fc, _, _, _ = fixed_amp(dc[i][0], dc[i][1] / Rv[i], dc[i][2] / Rv[i])
                    parts.append(f'{c} {fc["A"] / ref:+.2f} +- {fc["sA"] / ref:.2f}')
                # the photometric shift carried into the quasar estimator: A_sys = sum_c proj_c dcolour_c
                asys = sum(proj[i] * dc[i][1] for i in range(3)); jsys = sum(proj[i] * dc[i][2] for i in range(3))
                fsys, _, _, sjk = fixed_amp(dc[0][0], asys, jsys)
                prof = ' '.join(f'{r_:.0f}:{v * 1e3:+.1f}+-{e * 1e3:.1f}' for r_, v, e in
                                zip(dc[0][0], asys, np.sqrt(np.diag(jk_cov(jsys)))))
                print(f'{tag}: own estimator chi2 {x0:.1f}/{n0} (p {p0:.2f}), amp/ref {f1["A"] / ref:+.2f} +- {f1["sA"] / ref:.2f}; '
                      f'faint (r > 20) {ff["A"] / ref:+.2f} +- {ff["sA"] / ref:.2f}; per colour as A_V: {"; ".join(parts)}; '
                      f'through the quasar projection {fsys["A"] / ref:+.2f} +- {fsys["sA"] / ref:.2f}; N {int(n[f"{tag}__n___0"][0]) if f"{tag}__n___0" in n.files else -1}'
                      f'\n      projected profile (r_p: A_sys x1e-3): {prof}')
        # propagation: the photometric systematic is the standards' null (stars not selected through the W1-based QSO/LRG
        # cuts), projected through the quasar estimator, scaled to comoving like the measurement; subtracted, error in
        # quadrature. The leak is the E2 prediction (comoving), reported both ways; it can only lower the limit.
        rpm, a, jk = get('dr1_other_0.8_1.6__d' + C.COLOURS[0])
        proj, Rv, _ = get('proj_0.8_1.6')
        dc = [get(f'dr1_other_0.8_1.6__d{c}') for c in C.COLOURS]
        fsys, _, _, _ = fixed_amp(dc[0][0], sum(proj[i] * dc[i][1] for i in range(3)), sum(proj[i] * dc[i][2] for i in range(3)))
        cbar = float(np.average([rows[k]['cfac'] for k in ks], weights=wc))
        sc, ssc = fsys['A'] / ref * cbar, fsys['sA'] / ref * cbar
        lc, slc = ml * cbar, sml * cbar
        e0 = tot(mc, smc)
        v1, e1 = mc - sc, np.hypot(e0, ssc)
        v2, e2 = mc - sc - lc, np.hypot(e1, slc)
        print(f'\npropagated (comoving): measured {mc:.2f} +- {e0:.2f} (95 % upper, rho >= 0: {bayes_upper(mc, e0):.2f}); '
              f'photometric systematic (standards) {sc:+.2f} +- {ssc:.2f} -> {v1:.2f} +- {e1:.2f} (upper {bayes_upper(v1, e1):.2f}); '
              f'leak {lc:+.2f} +- {slc:.2f} also subtracted -> {v2:.2f} +- {e2:.2f} (upper {bayes_upper(v2, e2):.2f}); '
              f'R-150 QSO-target-star null taken literally {-3.81 * cbar:+.2f} -> {mc + 3.81 * cbar:.2f}')
        rd = rng.normal(v1, e1, 200000)
        om_mg = 2.3e-6 * (1.5 / (1 + zc))
        r1, r2 = rd / Sfac, rd / (om_mg / 5e-6 * Sfac)
        print(f'E4 with the photometric systematic included: Omega_d b_d ratio (z_eff / 0.36) {np.median(r1):.2f} '
              f'[{np.percentile(r1, 16):.2f}, {np.percentile(r1, 84):.2f}], positive-prior 95 % upper '
              f'{np.percentile(r1[r1 > 0], 95):.2f}; no-evolution (rho_c = S) excluded at '
              f'{(np.median(Sfac) - v1) / np.hypot(e1, np.std(Sfac)):.1f} sigma; measured / Mg II-budget prediction '
              f'{np.median(r2):.2f} [{np.percentile(r2, 16):.2f}, {np.percentile(r2, 84):.2f}], positive-prior 95 % upper '
              f'{np.percentile(r2[r2 > 0], 95):.2f}')


if __name__ == '__main__':
    what = sys.argv[1]
    {'nowise': run_nowise, 'prep': prep, 'prepstars': prep_stars, 'main': run_main, 'stars': run_stars, 'table': table}[what]()
