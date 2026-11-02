"""R-150: the halo-dust reddening profile around DESI DR1 ELGs and LRGs from the colours of DESI DR1 z > 2 quasars.

Everything is as fixed in results-R150.md M1 (and dust_z1_common.py): the statistic is the pair-weighted mean per-quasar
A_V,rest estimate around foregrounds minus the same around random foregrounds, divided by the targeting response, errors by
jackknife over 100 sky regions; then the power-law fit and the ratio to Menard+2010's z ~ 0.36 profile; Mg II flagged
pairs with and without; the achromatic (magnification) term from per-band brightness; the five nulls.

usage: python dust_z1_stack.py [main|nulls|all]   (writes data/dust/stack_*.npz; prints the tables)
"""
import os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dust_z1_common as C

NJK = 100
NB = len(C.RP_EDGES) - 1
RESP = lambda zf: 0.885 if zf >= 0.95 else 0.90          # M1: the F4 targeting response
RAN_FACTOR = 4                                            # randoms per data foreground used (subsampled from file 0)
REF_RP = 150.
HARTLAP = (NJK - NB - 2) / (NJK - 1)


def jk_regions(ra, dec, nside=16):
    import healpy as hp
    pix = hp.ang2pix(nside, ra, dec, lonlat=True, nest=True)
    cnt = np.bincount(pix, minlength=hp.nside2npix(nside))
    cum = np.cumsum(cnt) / cnt.sum()
    grp = np.minimum((cum * NJK).astype(int), NJK - 1)       # nested order keeps neighbouring pixels together
    return grp[pix]


def pair_sums(qra, qdec, qz, fra, fdec, fz, fw, zlo, zhi, nq, behind=False, extra=None, chunk=60000):
    """Per quasar and r_p bin, the summed foreground weight S[q, b]. extra(qi, fi) -> bool mask of pairs to drop."""
    S = np.zeros((nq, NB))
    rpsum = np.zeros(NB)
    wsum = np.zeros(NB)
    for i0 in range(0, nq, chunk):
        sl = slice(i0, min(nq, i0 + chunk))
        qi, fi, rp = C.pairs(qra[sl], qdec[sl], qz[sl], fra, fdec, fz, zlo, zhi, behind=behind)
        qi = qi + i0
        if extra is not None:
            keep = ~extra(qi, fi)
            qi, fi, rp = qi[keep], fi[keep], rp[keep]
        b = np.digitize(rp, C.RP_EDGES) - 1
        np.add.at(S, (qi, b), fw[fi])
        rpsum += np.bincount(b, weights=fw[fi] * rp, minlength=NB)
        wsum += np.bincount(b, weights=fw[fi], minlength=NB)
    return S, rpsum / np.maximum(wsum, 1e-30)


def stat(Sd, Sr, x, wq, reg):
    """Pair-weighted mean of x over data pairs minus over random pairs, per r_p bin; full and leave-one-region-out."""
    ok = np.isfinite(x)
    xx = np.where(ok, x, 0.)
    w = np.where(ok, wq, 0.)
    out = []
    for S in (Sd, Sr):
        num = np.zeros((NJK, NB)); den = np.zeros((NJK, NB))
        np.add.at(num, reg, (w * xx)[:, None] * S)
        np.add.at(den, reg, w[:, None] * S)
        N, Dn = num.sum(0), den.sum(0)
        full = N / Dn
        jk = (N[None] - num) / (Dn[None] - den)
        out.append((full, jk, Dn))
    full = out[0][0] - out[1][0]
    jk = out[0][1] - out[1][1]
    return full, jk, out[0][2]


def jk_cov(jk):
    d = jk - jk.mean(0)
    return (NJK - 1) / NJK * d.T @ d


def fit_powerlaw(rp, a, cov, gamma=None):
    """chi^2 fit of a (r_p/150)^-gamma; gamma fixed -> linear amplitude; free -> grid in gamma."""
    good = np.isfinite(a)
    rp, a, cov = rp[good], a[good], cov[np.ix_(good, good)]
    ci = np.linalg.inv(cov) * HARTLAP

    def amp(g):
        t = (rp / REF_RP) ** -g
        f = t @ ci @ t
        A = (t @ ci @ a) / f
        return A, 1 / np.sqrt(f), float((a - A * t) @ ci @ (a - A * t))
    if gamma is not None:
        A, s, chi2 = amp(gamma)
        return dict(A=A, sA=s, chi2=chi2, dof=len(a) - 1)
    gs = np.linspace(-0.5, 2.5, 301)
    res = [amp(g) for g in gs]
    chi = np.array([r[2] for r in res])
    i = int(np.argmin(chi))
    inside = gs[chi < chi[i] + 1]
    # amplitude interval from the 2-D profile (Delta chi^2 < 1 over the grid in (A, gamma))
    As = np.array([r[0] for r in res]); sA = np.array([r[1] for r in res])
    m = chi < chi[i] + 1
    lo = np.min((As - sA * np.sqrt(np.clip(chi[i] + 1 - chi, 0, None)))[m])
    hi = np.max((As + sA * np.sqrt(np.clip(chi[i] + 1 - chi, 0, None)))[m])
    return dict(A=As[i], gamma=gs[i], g_lo=inside.min(), g_hi=inside.max(), A_lo=lo, A_hi=hi, chi2=chi[i], dof=len(a) - 2)


def mgii_flagger(tid_q, zabs_tab):
    tids, zabs = zabs_tab
    order = np.argsort(tids)
    tids, zabs = tids[order], zabs[order]

    def flag(qi, fi, fz):
        t = tid_q[qi]
        lo = np.searchsorted(tids, t, 'left'); hi = np.searchsorted(tids, t, 'right')
        out = np.zeros(len(qi), bool)
        zf = fz[fi]
        for k in range(int((hi - lo).max()) if len(qi) else 0):
            has = lo + k < hi
            j = np.where(has, lo + k, 0)
            out |= has & (np.abs(zabs[j] - zf) < 1000. / 299792.458 * (1 + zf))
        return out
    return flag


def load_fg(tracer, seed=0):
    f = dict(np.load(os.path.join(C.D, f'fg_{tracer}.npz')))
    r = dict(np.load(os.path.join(C.D, f'ran_{tracer}.npz')))
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(r['z']))
    return f, r, perm


def sub(r, idx):
    return {k: v[idx] for k, v in r.items()}


def background():
    q = C.load_qso(zmin=C.ZQ_MIN)
    res, qr = C.colour_residuals(q)
    reg = jk_regions(q['ra'], q['dec'])
    return q, res, qr, reg


def mag_residuals(q, zstep=0.02):
    """Per-band magnitude residuals against the median in (redshift bin, imaging system), clipped at +-2 mag."""
    m = q['mag']
    M = np.c_[m['g'], m['r'], m['z'], m['W1']]
    out = np.full_like(M, np.nan)
    key = np.floor(q['z'] / zstep).astype(int) * 2 + q['north'].astype(int)
    order = np.argsort(key, kind='stable'); ks = key[order]
    for idx in np.split(order, np.flatnonzero(np.diff(ks)) + 1):
        if len(idx) >= 30:
            out[idx] = M[idx] - np.median(M[idx], axis=0)
    return np.clip(out, -2, 2)


def run_main():
    t0 = time.time()
    q, res, qr, reg = background()
    nq = len(q['z'])
    mres = mag_residuals(q)
    ab = np.load(os.path.join(C.D, 'mgii_tab.npz')) if os.path.exists(os.path.join(C.D, 'mgii_tab.npz')) else None
    if ab is None:
        import fitsio
        a = fitsio.read(os.path.join(C.D, 'mgii', 'MgII-Absorbers-DR1.fits'), columns=['TARGETID', 'Z_MGII'])
        np.savez(os.path.join(C.D, 'mgii_tab.npz'), tid=a['TARGETID'], z=a['Z_MGII'])
        ab = np.load(os.path.join(C.D, 'mgii_tab.npz'))
    print(f'backgrounds {nq}; Mg II absorbers {len(ab["tid"])}; jackknife regions {NJK}')
    results = {}
    for tracer in ('ELG', 'LRG'):
        f, r, perm = load_fg(tracer)
        for (zlo, zhi) in C.ZBINS[tracer]:
            zf = 0.5 * (zlo + zhi)
            sel_f = (f['z'] >= zlo) & (f['z'] < zhi)
            ridx = perm[(r['z'][perm] >= zlo) & (r['z'][perm] < zhi)][:RAN_FACTOR * sel_f.sum()]
            rr = sub(r, ridx)
            flag = mgii_flagger(q['targetid'], (ab['tid'], ab['z']))
            Sd, rpm = pair_sums(q['ra'], q['dec'], q['z'], f['ra'], f['dec'], f['z'], f['w'], zlo, zhi, nq)
            Sr, _ = pair_sums(q['ra'], q['dec'], q['z'], rr['ra'], rr['dec'], rr['z'], rr['w'], zlo, zhi, nq)
            Sd_m, _ = pair_sums(q['ra'], q['dec'], q['z'], f['ra'], f['dec'], f['z'], f['w'], zlo, zhi, nq,
                                extra=lambda qi, fi: flag(qi, fi, f['z']))
            Sr_m, _ = pair_sums(q['ra'], q['dec'], q['z'], rr['ra'], rr['dec'], rr['z'], rr['w'], zlo, zhi, nq,
                                extra=lambda qi, fi: flag(qi, fi, rr['z']))
            out = dict(rp=rpm, npair=(Sd > 0).sum(0), wpair=Sd.sum(0))
            for law in ('SMC', 'MW'):
                ah, w = C.av_estimator(res, qr, zf, law)
                full, jk, _ = stat(Sd, Sr, ah, w, reg)
                out[f'A_{law}'] = full / RESP(zf); out[f'jk_{law}'] = jk / RESP(zf)
                if law == 'SMC':
                    fm, jm, _ = stat(Sd_m, Sr_m, ah, w, reg)
                    out['A_noMgII'] = fm / RESP(zf); out['jk_noMgII'] = jm / RESP(zf)
                    # Galactic-dust split (null 5)
                    med = np.median(q['ebv'])
                    for tag, s in (('lowEBV', q['ebv'] < med), ('highEBV', q['ebv'] >= med)):
                        fs, js, _ = stat(Sd, Sr, np.where(s, ah, np.nan), w, reg)
                        out[f'A_{tag}'] = fs / RESP(zf); out[f'jk_{tag}'] = js / RESP(zf)
            # single colours as A_V (SMC response), equal weights
            R = C.response(zf, 'SMC')
            for i, cname in enumerate(C.COLOURS):
                full, jk, _ = stat(Sd, Sr, res[:, i] / R[i], np.ones(nq), reg)
                out[f'A_{cname}'] = full / RESP(zf); out[f'jk_{cname}'] = jk / RESP(zf)
            # per-band brightness shifts (magnification + extinction)
            for i, bname in enumerate(('g', 'r', 'z', 'W1')):
                full, jk, _ = stat(Sd, Sr, mres[:, i], np.ones(nq), reg)
                out[f'dm_{bname}'] = full; out[f'jkdm_{bname}'] = jk
            results[f'{tracer}_{zlo}_{zhi}'] = out
            np.savez(os.path.join(C.D, f'stack_{tracer}_{zlo}_{zhi}.npz'), **out)
            print(f'{tracer} {zlo}-{zhi} done ({time.time() - t0:.0f} s): pairs {int(out["npair"].sum())}', flush=True)
    return results


def run_nulls():
    t0 = time.time()
    q, res, qr, reg = background()
    nq = len(q['z'])
    out = {}
    # (2) random positions as foregrounds, (3) foregrounds rotated by +2 deg in RA: ELG 1.1-1.3 and the ELG z range
    for tracer, (zlo, zhi) in (('ELG', (0.8, 1.6)), ('LRG', (0.4, 1.1))):
        zf = 1.2 if tracer == 'ELG' else 0.75
        f, r, perm = load_fg(tracer, seed=1)
        ah, w = C.av_estimator(res, qr, zf)
        n = ((f['z'] >= zlo) & (f['z'] < zhi)).sum()
        sel = perm[(r['z'][perm] >= zlo) & (r['z'][perm] < zhi)]
        ra1, ra2 = sub(r, sel[:n]), sub(r, sel[n:(1 + RAN_FACTOR) * n])
        S1, rpm = pair_sums(q['ra'], q['dec'], q['z'], ra1['ra'], ra1['dec'], ra1['z'], ra1['w'], zlo, zhi, nq)
        S2, _ = pair_sums(q['ra'], q['dec'], q['z'], ra2['ra'], ra2['dec'], ra2['z'], ra2['w'], zlo, zhi, nq)
        full, jk, _ = stat(S1, S2, ah, w, reg)
        out[f'random_{tracer}'] = (rpm, full, jk)
        rot = (f['ra'] + 2.0) % 360
        Sd, _ = pair_sums(q['ra'], q['dec'], q['z'], rot, f['dec'], f['z'], f['w'], zlo, zhi, nq)
        full, jk, _ = stat(Sd, S2, ah, w, reg)
        out[f'rotated_{tracer}'] = (rpm, full, jk)
        print(f'nulls 2-3 {tracer} done ({time.time() - t0:.0f} s)', flush=True)
    # (1) foregrounds behind the quasar: ELGs at z_f > z_q + 0.1, quasars 0.6 < z_q < 1.5
    qa = C.load_qso()
    s = (qa['z'] > 0.6) & (qa['z'] < 1.5)
    qa = {k: (v[s] if k != 'mag' else {b: m[s] for b, m in v.items()}) for k, v in qa.items()}
    resa, qra = C.colour_residuals(qa)
    rega = jk_regions(qa['ra'], qa['dec'])
    f, r, perm = load_fg('ELG', seed=2)
    ah, w = C.av_estimator(resa, qra, 1.2)
    rr = sub(r, perm[:RAN_FACTOR * len(f['z'])])
    nqa = len(qa['z'])
    Sd, rpm = pair_sums(qa['ra'], qa['dec'], qa['z'], f['ra'], f['dec'], f['z'], f['w'], 0.8, 1.6, nqa, behind=True)
    Sr, _ = pair_sums(qa['ra'], qa['dec'], qa['z'], rr['ra'], rr['dec'], rr['z'], rr['w'], 0.8, 1.6, nqa, behind=True)
    full, jk, _ = stat(Sd, Sr, ah, w, rega)
    out['behind_ELG'] = (rpm, full, jk)
    print(f'null 1 done ({time.time() - t0:.0f} s): {nqa} quasars at 0.6-1.5, pairs {int((Sd > 0).sum())}', flush=True)
    # (4) stars as backgrounds: QSO targets with a stellar DESI spectrum
    import fitsio
    cols = ['RA', 'DEC', 'SPECTYPE', 'ZWARN', 'Z_RR', 'DESI_TARGET', 'PHOTSYS', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1',
            'FLUX_W2', 'FLUX_IVAR_W1', 'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1',
            'MW_TRANSMISSION_W2']
    a = fitsio.read(os.path.join(C.D, 'qso', 'QSO_full.dat.fits'), columns=cols)
    a = a[(a['SPECTYPE'] == 'STAR') & (a['ZWARN'] == 0) & (np.abs(a['Z_RR']) < 0.002) & ((a['DESI_TARGET'] & 4) > 0)]
    st = {k.lower(): a[k] for k in cols if k not in ('SPECTYPE', 'PHOTSYS')}
    st['north'] = a['PHOTSYS'] == 'N'
    st['mag'] = C.mags(st)
    good = np.ones(len(a), bool)
    for b in ('g', 'r', 'z', 'W1'):
        good &= np.isfinite(st['mag'][b])
    good &= (st['flux_w1'] * np.sqrt(st['flux_ivar_w1']) > 3) & (st['mag']['r'] > 16.5) & (st['mag']['r'] < 23)
    st = {k: (v[good] if k != 'mag' else {b: m[good] for b, m in v.items()}) for k, v in st.items()}
    st['z'] = np.full(good.sum(), 5.0)                  # stars: one colour cell per system/quartile; always 'behind'
    ress, qrs = C.colour_residuals(st, zstep=10.)
    regs = jk_regions(st['ra'], st['dec'])
    ns = len(st['z'])
    for tracer, (zlo, zhi), zf in (('ELG', (0.8, 1.6), 1.2), ('LRG', (0.4, 1.1), 0.75)):
        f, r, perm = load_fg(tracer, seed=3)
        sel = (f['z'] >= zlo) & (f['z'] < zhi)
        rr = sub(r, perm[:RAN_FACTOR * sel.sum()])
        ah, w = C.av_estimator(ress, qrs, zf)
        Sd, rpm = pair_sums(st['ra'], st['dec'], st['z'], f['ra'], f['dec'], f['z'], f['w'], zlo, zhi, ns)
        Sr, _ = pair_sums(st['ra'], st['dec'], st['z'], rr['ra'], rr['dec'], rr['z'], rr['w'], zlo, zhi, ns)
        full, jk, _ = stat(Sd, Sr, ah, w, regs)
        out[f'stars_{tracer}'] = (rpm, full, jk)
    print(f'null 4 done ({time.time() - t0:.0f} s): {ns} stars', flush=True)
    np.savez(os.path.join(C.D, 'stack_nulls.npz'), **{f'{k}__{i}': v[i] for k, v in out.items() for i in range(3)})
    return out


def run_slope():
    """Secondary estimator: the Redrock continuum slope (dust_z1_common.redrock_slope), added to each bin's stack file."""
    q, res, qr, reg = background()
    nq = len(q['z'])
    sres = C.slope_residual(q, C.redrock_slope(q))
    for tracer in ('ELG', 'LRG'):
        f, r, perm = load_fg(tracer)
        for (zlo, zhi) in C.ZBINS[tracer]:
            zf = 0.5 * (zlo + zhi)
            sel_f = (f['z'] >= zlo) & (f['z'] < zhi)
            ridx = perm[(r['z'][perm] >= zlo) & (r['z'][perm] < zhi)][:RAN_FACTOR * sel_f.sum()]
            rr = sub(r, ridx)
            Sd, _ = pair_sums(q['ra'], q['dec'], q['z'], f['ra'], f['dec'], f['z'], f['w'], zlo, zhi, nq)
            Sr, _ = pair_sums(q['ra'], q['dec'], q['z'], rr['ra'], rr['dec'], rr['z'], rr['w'], zlo, zhi, nq)
            full, jk, _ = stat(Sd, Sr, sres / C.slope_response(zf), np.ones(nq), reg)
            fn = os.path.join(C.D, f'stack_{tracer}_{zlo}_{zhi}.npz')
            out = dict(np.load(fn))
            out['A_slope'] = full / RESP(zf); out['jk_slope'] = jk / RESP(zf)
            np.savez(fn, **out)
            print(f'slope {tracer} {zlo}-{zhi} done', flush=True)


def run_quartiles():
    """Post-M1 diagnostic (results-R150.md D2): the SMC estimator stacked per r-magnitude quartile, and each quartile's
    within-quartile colour-magnitude slope dA/dr, to separate dust (the same in every quartile) from the magnification
    leak through the colour-luminosity relation (proportional to dA/dr times the induced magnitude shift)."""
    q, res, qr, reg = background()
    nq = len(q['z'])
    for tracer in ('ELG', 'LRG'):
        f, r, perm = load_fg(tracer)
        for (zlo, zhi) in C.ZBINS[tracer]:
            zf = 0.5 * (zlo + zhi)
            sel_f = (f['z'] >= zlo) & (f['z'] < zhi)
            ridx = perm[(r['z'][perm] >= zlo) & (r['z'][perm] < zhi)][:RAN_FACTOR * sel_f.sum()]
            rr = sub(r, ridx)
            Sd, _ = pair_sums(q['ra'], q['dec'], q['z'], f['ra'], f['dec'], f['z'], f['w'], zlo, zhi, nq)
            Sr, _ = pair_sums(q['ra'], q['dec'], q['z'], rr['ra'], rr['dec'], rr['z'], rr['w'], zlo, zhi, nq)
            np.savez_compressed(os.path.join(C.D, f'S_{tracer}_{zlo}_{zhi}.npz'), Sd=Sd.astype(np.float32), Sr=Sr.astype(np.float32))
            ah, w = C.av_estimator(res, qr, zf)
            fn = os.path.join(C.D, f'stack_{tracer}_{zlo}_{zhi}.npz')
            out = dict(np.load(fn))
            slopes = []
            for k in range(4):
                s = (qr == k) & np.isfinite(ah)
                rm = q['mag']['r'][s]
                slopes.append(np.polyfit(rm - rm.mean(), ah[s], 1)[0])
                full, jk, _ = stat(Sd, Sr, np.where(qr == k, ah, np.nan), w, reg)
                out[f'A_q{k}'] = full / RESP(zf); out[f'jk_q{k}'] = jk / RESP(zf)
            out['cm_slope'] = np.array(slopes)
            np.savez(fn, **out)
            print(f'quartiles {tracer} {zlo}-{zhi}: slopes {np.round(slopes, 3)}', flush=True)


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if what in ('main', 'all'):
        run_main()
    if what in ('nulls', 'all'):
        run_nulls()
    if what in ('slope', 'all'):
        run_slope()
    if what == 'quartiles':
        run_quartiles()
