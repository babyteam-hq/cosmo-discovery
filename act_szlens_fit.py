"""The fit of R-154 (plan-R154.md and its addendum): per redshift bin R = M_lens/M_SZ, and the evolution exponent zeta, from the
ACT DR6 lensing stack on the ACT DR6 SZ clusters, with the method check's calibrations applied.

The calibrations come from act_szlens_inject.py and act_szlens.py selmock / cib:
- the lensing response r_b, from paired NFW injections through the emulated hardened estimator;
- the tSZ and CIB leakage l_b, in units of the bin's one-halo amplitude at R = 1, for the clusters whose templates were
  subtracted (inside the release's SZ mask) and for the rest, weighted by w_T = N_L^MV / N_L^T;
- the selection bias B(R_true), solved self-consistently;
- the contamination f.
The data vector is corrected as d -> (d - l_b m1_b(R=1)) / r_b, then fitted, then R_fit is mapped to R_true through B, then
divided by (1-f)^(1/0.85).

usage: python -I act_szlens_fit.py fit [--variant baseline] [--sample all] [--cov sims|rotations] [--sys]
"""
import argparse, glob, os, sys, time, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import act_szlens as asl   # noqa: E402
import erass1_act as ea    # noqa: E402
import erass1_cmblens as r121  # noqa: E402

OUT = asl.OUT
log = ea.log
ZB = asl.ZBINS
ZMID = [0.33, 0.64, 0.89, 1.23]
RG = np.exp(np.linspace(np.log(0.2), np.log(8.0), 400))


class Model:
    def __init__(self, cscale=1.0, mis_arcmin=0.0, inner_only=False, fix2h=False):
        self.fix2h = fix2h
        tag = '' if cscale == 1.0 else f'_c{cscale}'
        self.hm = asl.ACTHaloModel(lmax=3000, cscale=cscale)
        lnMg, zg, k1g, k2g, bias = asl.model_grid(self.hm, tag)
        if mis_arcmin > 0:
            s = mis_arcmin * ea.ARCMIN
            k1g = k1g * np.exp(-0.5 * (self.hm.ell * s) ** 2)[None, None, :]
        f_in, f_out = ea.filters(self.hm, ea.LMIN, ea.LMAX)
        self.ag = ea.annulus_grid(self.hm, (lnMg, zg, k1g, k2g, bias), f_in, f_out)
        self.use = np.arange(5) if inner_only else np.arange(ea.NB)

    def m12(self, nodes, wn, z):
        m1, m2 = ea.stacked_model_posterior(self.hm, self.ag, nodes, wn, z)
        return m1[self.use], m2[self.use]


def load_cluster_inputs(mass='A10_M500c'):
    c = asl.load_act(mass)
    s = asl.selection(c)
    keep = s['keep']
    M = c['M'][keep]; z = c['z'][keep]
    ep = c[mass + '_errPlus'] * 1.0 if (mass + '_errPlus') in c else None
    from astropy.io import fits
    d = fits.open(asl.CAT)[1].data
    eP = np.asarray(d[mass + '_errPlus'], float)[keep] * 1e14; eM = np.asarray(d[mass + '_errMinus'], float)[keep] * 1e14
    sig = 0.5 * (np.log(M + eP) - np.log(np.maximum(M - eM, 0.05 * M)))
    sig = np.where(np.isfinite(sig) & (sig > 0), sig, np.nanmedian(sig))
    zerr = c['redshiftErr'][keep]; zerr = np.where(np.isfinite(zerr) & (zerr > 0), zerr, 0.0)
    return dict(c=c, keep=keep, M=M, z=z, sig=sig, zerr=zerr, insz=s['insz'][keep], q=c['fixed_SNR'][keep], spec=c['spec'][keep])


def bin_masks(z, extra=None):
    out = []
    for lo, hi in ZB:
        m = (z >= lo) & (z < hi)
        if extra is not None:
            m &= extra
        out.append(m)
    return out


# --------------------------------------------------------------------------------------------- calibrations (method check)
def calibrations(ci, args):
    """per bin: response r, leakage l (units of m1(R=1)), contamination f; and the selection map B(R_true) per bin.
    Each returned with its variants for the systematic budget."""
    inj = {t: list(np.load(os.path.join(OUT, 'inject', f'inject_rows{t}.npy'), allow_pickle=True)) for t in ('_clhard', '_hard15', '_hard30')
           if os.path.exists(os.path.join(OUT, 'inject', f'inject_rows{t}.npy'))}
    cibi = list(np.load(os.path.join(OUT, 'inject', 'cib_inject_clhard.npy'), allow_pickle=True)) if os.path.exists(os.path.join(OUT, 'inject', 'cib_inject_clhard.npy')) \
        else list(np.load(os.path.join(OUT, 'inject', 'cib_inject.npy'), allow_pickle=True))
    cibrows = list(np.load(os.path.join(OUT, 'cib_rows.npy'), allow_pickle=True))
    sel = list(np.load(os.path.join(OUT, 'selmock_rows.npy'), allow_pickle=True))
    nl = {v: np.loadtxt(os.path.join(ea.ACT, f'N_L_kk_act_dr6_lensing_v1_{v}.txt')) for v in ('baseline', 'tonly')}
    Lw = np.arange(500, 2001)
    wT = float(np.mean(np.interp(Lw, nl['baseline'][:, 0], nl['baseline'][:, 1]) / np.interp(Lw, nl['tonly'][:, 0], nl['tonly'][:, 1])))
    bins = bin_masks(ci['z'])
    fin = [float(ci['insz'][b].mean()) for b in bins]
    out = dict(wT=wT, fin=fin, r={}, l={}, f={}, B={})

    def rows_at(tag, z, msz=None):
        rr = [r for r in inj[tag] if abs(r['z'] - z) < 1e-6 and (msz is None or abs(r['Msz'] - msz) / msz < 0.01)]
        return rr

    MED = [1.96e14, 2.08e14, 1.93e14, 1.76e14]
    base = '_clhard' if '_clhard' in inj else '_hard30'
    for b, (z, msz) in enumerate(zip(ZMID, MED)):
        rr = rows_at(base, z, msz)
        # response: the mean of the R = 1.35 and 3.4 injections (they agree within errors); variants: each alone
        resp = np.array([r['lens'] for r in rr]); resp_e = np.array([r['lens_err'] for r in rr])
        out['r'][b] = dict(central=float(np.average(resp, weights=1 / resp_e ** 2)), err=float(1 / np.sqrt(np.sum(1 / resp_e ** 2))),
                           byR={float(r['R']): float(r['lens']) for r in rr})
        r0 = rr[0]
        # tSZ: outside the mask hardened-only (UPP truth), inside the mask template-subtracted; MV weight w_T
        cb = [x for x in cibi if abs(x['z'] - z) < 1e-6]
        meas = cibrows[b]['f545'] / cibrows[b]['f857']
        best = [x for x in cb if x['shape'] == 'point']
        # the SED nearest the measured 545/857 ratio is listed first by the injection; choose by min |model - measured|
        def cib_val(x, key):
            return x[key]
        cib_central_out = float(np.mean([x['cib_h'] for x in cb]))
        cib_central_in = float(np.mean([x['cib_sub'] for x in cb]))
        cib_spread = float(np.ptp([x['cib_h'] for x in cb]) / 2) if len(cb) > 1 else 0.0
        cib_flux_frac = cibrows[b]['e857'] / max(cibrows[b]['f857'], 1e-9)
        tsz_out = float(r0['upp_h']); tsz_in = float(r0['upp_sub'])
        # per-cluster weighting: the stack's leakage is sum_out kappa_i l(M_i) / sum_all kappa_i with l(M) = l(M_med) (M/M_med)^alpha
        # for the tSZ (alpha 2.0, measured 1.9-2.2 from the two injected masses), kappa_i ~ M_i^0.85; the CIB per unit kappa is
        # taken mass-independent (alpha 0)
        Mb = ci['M'][bins[b]]; ob = ~ci['insz'][bins[b]]; kw = Mb ** 0.85
        def frac(alpha, sel):
            return float(np.sum(kw[sel] * (Mb[sel] / msz) ** alpha) / np.sum(kw))
        A_out = {al: frac(al, ob) for al in (0.0, 1.5, 2.0, 2.5)}; A_in = frac(0.0, ~ob)
        l_central = wT * (A_out[2.0] * tsz_out + A_out[0.0] * cib_central_out + A_in * (tsz_in + cib_central_in))
        # variants for the systematic: the tSZ truth profile (Planck 13), the hardening profile (1.5', 3'), the CIB SED/shape
        # spread and its flux error, the in-mask fraction treated as all-subtracted / none-subtracted is not a variant (known)
        var = {}
        var['tsz_p13'] = wT * (A_out[2.0] * r0['p13_h'] + A_out[0.0] * cib_central_out + A_in * (r0['p13_sub'] + cib_central_in)) - l_central
        var['tsz_alpha'] = wT * max(abs(A_out[1.5] - A_out[2.0]), abs(A_out[2.5] - A_out[2.0])) * tsz_out
        for t in ('_hard15', '_hard30'):
            if t in inj and t != base:
                rt = rows_at(t, z, msz)
                if rt:
                    var['hard' + t] = wT * (A_out[2.0] * rt[0]['upp_h'] + A_out[0.0] * cib_central_out + A_in * (rt[0]['upp_sub'] + cib_central_in)) - l_central
        var['cib_sed'] = wT * A_out[0.0] * cib_spread
        var['cib_flux'] = wT * A_out[0.0] * abs(cib_central_out) * cib_flux_frac
        out['l'][b] = dict(central=float(l_central), A_out=A_out[2.0], A_out0=A_out[0.0], A_in=A_in, tsz_out=tsz_out, tsz_in=tsz_in, cib_out=cib_central_out, cib_in=cib_central_in,
                           var={k: float(v) for k, v in var.items()})
        # contamination: chance association 1-6 % (DR6 paper) on the q < 5.5 candidates whose noise false-positive fraction is
        # large; central 0, the systematic one-sided up to 4 % (6 % at z > 1, where the optical catalogues are denser)
        out['f'][b] = dict(central=0.0, max=0.06 if b == 3 else 0.04)
        # selection: B(R_true) at sigma_int 0.2 (sigma8 0.80, no confirmation term), and its variants
        lo, hi = ZB[b]
        sb = [r for r in sel if r['lo'] == lo]
        base_s = sorted([(r['Rtrue'], r['B']) for r in sb if r['case'][1:] == (0.2, 0.80, 0)])
        out['B'][b] = dict(central=base_s, var={})
        for lab, cond in [('sint0.1', (0.1, 0.80, 0)), ('sint0.3', (0.3, 0.80, 0)), ('s8_0.77', (0.2, 0.77, 0)), ('s8_0.83', (0.2, 0.83, 0)),
                          ('optM50_1.0', (0.2, 0.80, 1.0e14)), ('optM50_1.5', (0.2, 0.80, 1.5e14))]:
            pts = sorted([(r['Rtrue'], r['B'] / dict(base_s)[r['Rtrue']]) for r in sb if r['case'][1:] == cond and r['Rtrue'] in dict(base_s)])
            out['B'][b]['var'][lab] = pts      # the ratio to the central B at the R_true values computed
    return out


def B_of(Bc, Rtrue, ratio_pts=None):
    """B(R_true) by interpolation in ln R over the mock grid (clamped); times a variant's ratio if given."""
    x = np.log([p[0] for p in Bc]); y = [p[1] for p in Bc]
    v = np.interp(np.log(Rtrue), x, y)
    if ratio_pts:
        xr = np.log([p[0] for p in ratio_pts]); yr = [p[1] for p in ratio_pts]
        v = v * np.interp(np.log(Rtrue), xr, yr)
    return v


def Rfit_to_true(Rfit, Bc, ratio_pts=None):
    """invert R_fit = B(R_true) R_true on a fine grid."""
    Rt = np.exp(np.linspace(np.log(0.1), np.log(12), 3000))
    Rf = B_of(Bc, Rt, ratio_pts) * Rt
    return np.interp(Rfit, Rf, Rt)


# --------------------------------------------------------------------------------------------- data, covariance
def data_vectors(ci, variant, bins, use):
    D = np.load(os.path.join(OUT, f'data_{variant}.npz'))
    prof = D['prof'].astype(float); rs = D['random']
    return [np.nanmean(prof[b], 0)[use] - rs[use] for b in bins]


def sims_cov(bins, use):
    P = np.load(os.path.join(OUT, 'sims_prof_baseline.npy'), mmap_mode='r')
    Rr = np.load(os.path.join(OUT, 'sims_rand_baseline.npy'))
    done = np.load(os.path.join(OUT, 'sims_done_baseline.npy'))
    idx = np.where(done)[0]
    V = []
    for i in idx:
        p = np.asarray(P[i], float)
        V.append(np.concatenate([np.nanmean(p[b], 0)[use] - Rr[i][use] for b in bins]))
    V = np.array(V)
    return V


def rot_cov(variant, bins, use, n_eff):
    P = np.load(os.path.join(OUT, f'rot_{variant}.npy'))
    V = []
    for r in range(len(P)):
        v = []
        for b, n in zip(bins, n_eff):
            pr = P[r][b]; ok = np.isfinite(pr).all(1)
            m = np.nanmean(pr[ok], 0)[use]
            v.append(m * np.sqrt(ok.sum() / n))     # noise-dominated: a stack of n_ok rescaled to the bin's n
        V.append(np.concatenate(v))
    V = np.array(V)
    return V


def hartlap_inv(V):
    C = np.cov(V.T); n, p = V.shape
    return np.linalg.inv(C) * (n - p - 2) / (n - 1), C


# --------------------------------------------------------------------------------------------- fits
def fit_bin(d, Ci, mod, nodes, wn, z, Rgrid=RG):
    chi = np.empty(len(Rgrid)); a2 = np.empty(len(Rgrid))
    for i, R in enumerate(Rgrid):
        m1, m2 = mod.m12(nodes + np.log(R), wn, z)
        if mod.fix2h:
            r = d - m1 - m2; chi[i], a2[i] = r @ Ci @ r, 1.0
        else:
            chi[i], a2[i] = ea.profile_chi2(d, Ci, m1, m2)
    return chi, a2


def interval(Rgrid, chi, dchi):
    i = np.argmin(chi); c0 = chi[i]
    ok = chi - c0 < dchi
    lo = Rgrid[np.argmax(ok)]; hi = Rgrid[len(ok) - 1 - np.argmax(ok[::-1])]
    return Rgrid[i], lo, hi, c0


def amp1h(d, Ci, mod, nodes, wn, z):
    m1, m2 = mod.m12(nodes, wn, z)
    X = np.vstack([m1, m2]).T
    F = X.T @ Ci @ X; p = np.linalg.solve(F, X.T @ Ci @ d)
    return p[0], np.sqrt(np.linalg.inv(F)[0, 0]), p[1]


def trend_fit(ds, Ci_joint, mod, ci, bins, corr, eta=0.0, grid=None):
    """ln R(z) = ln R0 + zeta ln((1+z)/1.45) [+ eta ln(M/3e14)], all bins jointly, one free two-halo amplitude per bin."""
    z = ci['z']; M = ci['M']
    nodes_all, wn = ea.posterior_nodes(np.log(ci['M']), ci['sig'])
    lnR0g = np.linspace(np.log(0.3), np.log(5.0), 52) if grid is None else grid[0]
    zg = np.linspace(-4.0, 4.0, 81) if grid is None else grid[1]
    nu = len(mod.use)
    dcat = np.concatenate(ds)
    chi = np.full((len(lnR0g), len(zg)), np.inf)
    for i, l0 in enumerate(lnR0g):
        for j, ze in enumerate(zg):
            m1s, m2s = [], []
            for b, msk in enumerate(bins):
                lnR = l0 + ze * np.log((1 + z[msk]) / 1.45) + eta * np.log(M[msk] / 3e14)
                m1, m2 = mod.m12(nodes_all[msk] + lnR[:, None], wn, z[msk])
                m1s.append(m1); m2s.append(m2)
            m1c = np.concatenate(m1s)
            X = np.zeros((len(dcat), len(bins)))
            for b in range(len(bins)):
                X[b * nu:(b + 1) * nu, b] = m2s[b]
            r = dcat - m1c
            a = np.ones(len(bins)) if mod.fix2h else np.linalg.solve(X.T @ Ci_joint @ X, X.T @ Ci_joint @ r)
            rr = r - X @ a
            chi[i, j] = rr @ Ci_joint @ rr
    prof = chi.min(0)
    # the profile on a fine grid (cubic through the coarse one) for the minimum and the crossings
    from scipy.interpolate import CubicSpline
    zf = np.linspace(zg[0], zg[-1], 1201); pf = CubicSpline(zg, prof)(zf)
    j = np.argmin(pf); c0 = pf[j]
    def cross(d):
        ok = pf - c0 < d
        return zf[np.argmax(ok)], zf[len(ok) - 1 - np.argmax(ok[::-1])]
    lo68, hi68 = cross(1.0); lo95, hi95 = cross(3.84)
    i0 = np.argmin(chi[:, np.argmin(np.abs(zg - zf[j]))])
    prof = pf; zg = zf
    return dict(zeta=zg[j], lo68=lo68, hi68=hi68, lo95=lo95, hi95=hi95, R0=np.exp(lnR0g[i0]), chi2=c0,
                dchi2_zeta0=float(np.interp(0.0, zg, prof) - c0), dchi2_zeta18=float(np.interp(1.8, zg, prof) - c0), prof=(zg, prof))


def corrected_vectors(ds, mod, ci, bins, cal, var=None):
    """d -> (d - l m1(R=1)) / r per bin, with an optional variant {'l': {b: shift}, 'r': {b: value}}."""
    nodes_all, wn = ea.posterior_nodes(np.log(ci['M']), ci['sig'])
    out = []
    for b, msk in enumerate(bins):
        m1, _ = mod.m12(nodes_all[msk], wn, ci['z'][msk])
        l = cal['l'][b]['central'] + ((var or {}).get('l', {}).get(b, 0.0))
        r = (var or {}).get('r', {}).get(b, cal['r'][b]['central'])
        out.append((ds[b] - l * m1) / r)
    return out


def fit(args):
    ci = load_cluster_inputs(args.mass)
    extra = None
    if args.sample == 'nosz':
        extra = ~ci['insz']
    elif args.sample == 'insz':
        extra = ci['insz']
    elif args.sample == 'q>5':
        extra = ci['q'] > 5
    elif args.sample == 'q4-5':
        extra = ci['q'] <= 5
    elif args.sample == 'q>5 nosz':
        extra = (ci['q'] > 5) & ~ci['insz']
    elif args.sample == 'q>5 insz':
        extra = (ci['q'] > 5) & ci['insz']
    elif args.sample == 'q4-5 nosz':
        extra = (ci['q'] <= 5) & ~ci['insz']
    elif args.sample == 'spec':
        extra = ci['spec']
    elif args.sample == 'M<3e14':
        extra = ci['M'] < 3e14
    elif args.sample == 'M>3e14':
        extra = ci['M'] >= 3e14
    bins = bin_masks(ci['z'], extra)
    mod = Model(cscale=args.cscale, mis_arcmin=args.mis, inner_only=args.inner)
    use = mod.use; nu = len(use)
    ds = data_vectors(ci, args.variant, bins, use)
    if args.cov == 'sims':
        V = sims_cov(bins, use)
        nm = V.mean(0)
    else:
        V = rot_cov(args.variant, bins, use, [b.sum() for b in bins])
        nm = np.zeros(V.shape[1])
    Cij, Cj = hartlap_inv(V)
    log(f'{args.variant} {args.sample} mass {args.mass} c x{args.cscale} mis {args.mis}\' inner {args.inner}: n per bin {[int(b.sum()) for b in bins]}, '
        f'{len(V)} {args.cov}; null-mean per bin (central annulus x1e4) {np.round(nm[::nu] * 1e4, 2)}')
    cal = calibrations(ci, args) if not args.nocal else None
    if cal:
        log(f'calibrations: w_T {cal["wT"]:.3f}; in-mask fractions {np.round(cal["fin"], 3)}; response {[round(cal["r"][b]["central"], 3) for b in range(4)]} '
            f'± {[round(cal["r"][b]["err"], 3) for b in range(4)]}; leakage {[round(cal["l"][b]["central"], 3) for b in range(4)]}')
    dsc = corrected_vectors(ds, mod, ci, bins, cal) if cal else ds
    nodes_all, wn = ea.posterior_nodes(np.log(ci['M']), ci['sig'])
    res = dict(bins=[])
    for b, msk in enumerate(bins):
        Ci = np.linalg.inv(Cj[b * nu:(b + 1) * nu, b * nu:(b + 1) * nu]) * (len(V) - nu - 2) / (len(V) - 1)
        a1, sa1, a2 = amp1h(dsc[b], Ci, mod, nodes_all[msk], wn, ci['z'][msk])
        chi, a2g = fit_bin(dsc[b], Ci, mod, nodes_all[msk], wn, ci['z'][msk])
        Rf, lo, hi, c0 = interval(RG, chi, 1.0); _, lo95, hi95, _ = interval(RG, chi, 3.84)
        if cal:
            tr = lambda x: Rfit_to_true(x, cal['B'][b]['central'])
            Rt, lot, hit, lo95t, hi95t = tr(Rf), tr(lo), tr(hi), tr(lo95), tr(hi95)
        else:
            Rt, lot, hit, lo95t, hi95t = Rf, lo, hi, lo95, hi95
        det = np.sqrt(max(chi[np.argmin(np.abs(RG - 1e-9))] - c0, 0)) if False else a1 / sa1
        row = dict(bin=ZB[b], n=int(msk.sum()), zmean=float(ci['z'][msk].mean()), Mmed=float(np.median(ci['M'][msk])), A1h=a1, sA1h=sa1, A2h_at_R1=a2,
                   Rfit=Rf, Rfit68=(lo, hi), R=Rt, R68=(lot, hit), R95=(lo95t, hi95t), chi2=c0, dof=nu - 2)
        res['bins'].append(row)
        log(f'  z {ZB[b][0]}-{ZB[b][1]} (n {msk.sum()}, <z> {row["zmean"]:.2f}): A_1h(R=1) {a1:.2f} ± {sa1:.2f} ({a1 / sa1:.1f} sigma), '
            f'R_fit {Rf:.2f} [{lo:.2f}, {hi:.2f}] -> R = M_lens/M_SZ {Rt:.2f} [{lot:.2f}, {hit:.2f}] (95 % [{lo95t:.2f}, {hi95t:.2f}]); '
            f'1-b {1 / Rt:.2f} [{1 / hit:.2f}, {1 / lot:.2f}]; chi2 {c0:.1f}/{nu - 2}')
    if args.trend:
        dst = dsc
        if cal:
            # the selection correction for the trend: applied per bin as the ratio B at that bin's fitted R (a multiplicative
            # rescaling of the data amplitude by 1/B^0.85 is equivalent to dividing R_fit by B)
            dst = [d / (B_of(cal['B'][b]['central'], res['bins'][b]['R']) ** 0.85) for b, d in enumerate(dsc)]
        tf = trend_fit(dst, Cij, mod, ci, bins, cal, eta=args.eta)
        res['trend'] = {k: v for k, v in tf.items() if k != 'prof'}
        log(f'  trend: zeta {tf["zeta"]:.2f} [{tf["lo68"]:.2f}, {tf["hi68"]:.2f}] (95 % [{tf["lo95"]:.2f}, {tf["hi95"]:.2f}]), R0(z=0.45) {tf["R0"]:.2f}; '
            f'chi2 {tf["chi2"]:.1f}; Delta chi2 at zeta = 0: {tf["dchi2_zeta0"]:.2f}, at zeta = 1.8: {tf["dchi2_zeta18"]:.2f}; eta {args.eta}')
    os.makedirs(os.path.join(OUT, 'fits'), exist_ok=True)
    tag = f'{args.variant}_{args.sample}_{args.mass}_c{args.cscale}_mis{args.mis}_in{int(args.inner)}_eta{args.eta}_{args.cov}{"_nocal" if args.nocal else ""}'
    json.dump(res, open(os.path.join(OUT, 'fits', tag + '.json'), 'w'), default=float, indent=1)
    return res, cal


def run_variant(ci, bins, mod, ds, Cj, nV, cal, var, eta=0.0, trend=True):
    """per-bin R (true, after the corrections of this variant) and the trend; var: {'l': {b: dl}, 'r': {b: r}, 'B': label,
    'f': {b: f}}."""
    nu = len(mod.use)
    dsc = corrected_vectors(ds, mod, ci, bins, cal, var)
    nodes_all, wn = ea.posterior_nodes(np.log(ci['M']), ci['sig'])
    Rs = []
    for b, msk in enumerate(bins):
        Ci = np.linalg.inv(Cj[b * nu:(b + 1) * nu, b * nu:(b + 1) * nu]) * (nV - nu - 2) / (nV - 1)
        chi, _ = fit_bin(dsc[b], Ci, mod, nodes_all[msk], wn, ci['z'][msk])
        Rf = RG[np.argmin(chi)]
        ratio = cal['B'][b]['var'].get(var.get('B')) if var.get('B') else None
        Rt = Rfit_to_true(Rf, cal['B'][b]['central'], ratio)
        f = var.get('f', {}).get(b, 0.0)
        Rs.append(Rt / (1 - f) ** (1 / 0.85))
    out = dict(R=Rs)
    if trend:
        dst = []
        for b, d in enumerate(dsc):
            ratio = cal['B'][b]['var'].get(var.get('B')) if var.get('B') else None
            f = var.get('f', {}).get(b, 0.0)
            # B and the contamination as amplitude rescalings of the bin's data (R_fit = B R_true -> amplitude ~ R^0.85)
            dst.append(d / (B_of(cal['B'][b]['central'], Rs[b], ratio) ** 0.85) / (1 - f))
        Cij = np.linalg.inv(Cj) * (nV - len(Cj) - 2) / (nV - 1)
        out['trend'] = {k: v for k, v in trend_fit(dst, Cij, mod, ci, bins, cal, eta=eta).items() if k != 'prof'}
    return out


def systematics(args):
    """every z-dependent systematic as its own variant of the whole fit (plan addendum point 8); shifts on ln R per bin and on
    zeta, added in quadrature (two-sided variants: half the up/down difference or the single shift)."""
    ci = load_cluster_inputs(args.mass)
    bins = bin_masks(ci['z'])
    base_mod = Model(mis_arcmin=0.5)
    V = sims_cov(bins, base_mod.use); Cj = np.cov(V.T)
    ds = data_vectors(ci, 'baseline', bins, base_mod.use)
    cal = calibrations(ci, args)
    runs = {}
    def go(name, var, mod=None, eta=0.0):
        mod = mod or base_mod
        if mod is not base_mod:
            dsm = data_vectors(ci, 'baseline', bins, mod.use); Vm = sims_cov(bins, mod.use); Cm = np.cov(Vm.T)
        else:
            dsm, Cm, Vm = ds, Cj, V
        t0 = time.time()
        r = run_variant(ci, bins, mod, dsm, Cm, len(Vm), cal, var, eta=eta, trend=not args.notrend)
        runs[name] = r
        log(f'variant {name}: R {np.round(r["R"], 3)}' + (f'; zeta {r["trend"]["zeta"]:.2f} [{r["trend"]["lo68"]:.2f}, {r["trend"]["hi68"]:.2f}]' if 'trend' in r else '') + f' ({time.time() - t0:.0f}s)')
        json.dump(runs, open(os.path.join(OUT, 'fits', f'systematics{args.tag}.json'), 'w'), default=float, indent=1)
    os.makedirs(os.path.join(OUT, 'fits'), exist_ok=True)
    go('central', {})
    for key in sorted({k for b in range(4) for k in cal['l'][b]['var']}):
        go('l_' + key, {'l': {b: cal['l'][b]['var'].get(key, 0.0) for b in range(4)}})
    go('r_hard30', {'r': {b: cal['r'][b]['central'] * 0.956 for b in range(4)}})
    go('r_plus1sig', {'r': {b: cal['r'][b]['central'] + cal['r'][b]['err'] for b in range(4)}})
    for lab in ('sint0.1', 'sint0.3', 's8_0.83', 'optM50_1.5'):
        go('B_' + lab, {'B': lab})
    go('f_max', {'f': {b: cal['f'][b]['max'] for b in range(4)}})
    go('mis0', {}, Model(mis_arcmin=0.0))
    go('c0.7', {}, Model(cscale=0.7, mis_arcmin=0.5))
    go('c1.3', {}, Model(cscale=1.3, mis_arcmin=0.5))
    go('inner30', {}, Model(mis_arcmin=0.5, inner_only=True))
    go('eta-0.18', {}, eta=-0.18)
    go('2h_fixed1', {}, Model(mis_arcmin=0.5, fix2h=True))
    go('2h_fixed1_inner30', {}, Model(mis_arcmin=0.5, fix2h=True, inner_only=True))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('cmd')
    p.add_argument('--variant', default='baseline')
    p.add_argument('--sample', default='all')
    p.add_argument('--mass', default='A10_M500c')
    p.add_argument('--cov', default='sims')
    p.add_argument('--cscale', type=float, default=1.0)
    p.add_argument('--mis', type=float, default=0.5)
    p.add_argument('--inner', action='store_true')
    p.add_argument('--eta', type=float, default=0.0)
    p.add_argument('--trend', action='store_true')
    p.add_argument('--nocal', action='store_true')
    p.add_argument('--notrend', action='store_true')
    p.add_argument('--tag', default='')
    a = p.parse_args()
    {'fit': fit, 'sys': systematics}[a.cmd](a)


if __name__ == '__main__':
    main()
