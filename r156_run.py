"""R-156: the measurement. Run only after the method check's points are closed (results-R156.md F2).

usage: python -u r156_run.py [--mocks 1,2] [--zshuffle TAG]
Writes data/r156/run.json and prints the per-bin table, nulls and bounds.
"""
import sys, json, time
import numpy as np
import healpy as hp
import r156_kappa as K
import r156_measure as M
import r156_theory as T

D = 'data/r156/'
BINS = M.BINS
FIT = slice(1, 14)       # d = 20-150 Mpc/h (bin centres 25 ... 145)
FIT30 = slice(2, 14)     # d >= 30
NSHUF = int(sys.argv[sys.argv.index('--nshuf') + 1]) if '--nshuf' in sys.argv else 300
A_RESP, A_RESP_ERR = 0.97, 0.05
ETA_ZSUCC = 0.02         # redshift-success systematic on eta
C_FRAC = 0.5             # assumed fractional uncertainty on C(z)


def region_sums(sh, order, bounds, w_num, w_den):
    """Per-region (NJK, nd, nmu) sums of w_num Q and w_den R over galaxies sorted by region."""
    num = np.zeros((M.NJK,) + sh.Q.shape[1:])
    den = np.zeros_like(num)
    for g in range(M.NJK):
        ii = order[bounds[g]:bounds[g + 1]]
        if len(ii):
            num[g] = np.tensordot(w_num[ii], sh.Q[ii], axes=(0, 0))
            den[g] = np.tensordot(w_den[ii], sh.R[ii], axes=(0, 0))
    return num, den


def ell_of(sh, num, den, ell):
    return M.Shells.multipoles(sh, num, den, ell)


def mu_factor(sh, ell):
    """Estimator response to a pure L_ell(mu) term after binning in mu (exact for 10 bins)."""
    edges = np.linspace(-1, 1, sh.nmu + 1)
    x = np.linspace(-1, 1, 200001)
    from scipy.special import eval_legendre
    binned = np.array([eval_legendre(ell, x[(x >= edges[i]) & (x <= edges[i + 1])]).mean() for i in range(sh.nmu)])
    return float(binned @ sh.L[ell])


def main():
    t0 = time.time()
    mocks = [int(x) for x in sys.argv[sys.argv.index('--mocks') + 1].split(',')] if '--mocks' in sys.argv else []
    zsh = sys.argv[sys.argv.index('--zshuffle') + 1] if '--zshuffle' in sys.argv else None
    s, sel, k, u, rel = M.prepare(False)
    tag = sys.argv[sys.argv.index('--tag') + 1] if '--tag' in sys.argv else 'data'
    if '--blind' in sys.argv:  # pipeline test: kappa_hat permuted, no signal can survive
        k = k.copy()
        ii = np.where(sel)[0]
        k[ii] = k[np.random.default_rng(99).permutation(ii)]
    sh = M.Shells(tag)
    kid = sh.kidx
    gj = sel[kid]
    zj = s['Z'][kid]
    capj = s['CAP'][kid]
    reg = M.jk_regions(s['RA'][kid], s['DEC'][kid], np.ones(len(kid), bool))
    order = np.argsort(reg, kind='stable')
    bounds = np.searchsorted(reg[order], np.arange(M.NJK + 1))
    f1, f3 = mu_factor(sh, 1), mu_factor(sh, 3)
    kj = np.where(gj, k[kid], 0.0)
    uj = np.where(gj, u[kid], 0.0)
    e = np.load(D + 'eta.npz')
    eta_z = e['eta']
    out = {'mu_factor_dipole': f1, 'mu_factor_oct': f3, 'bins': {}}

    # eta per jackknife region (statistical error of eta, propagated jointly)
    rg_s = M.jk_regions(s['RA'], s['DEC'], sel)
    cuts = np.linspace(19.3, 19.5, 5)
    iz = np.round((s['Z'] - M.ZC[0]) / 0.02).astype(int)
    okz = sel & (np.abs(s['Z'] - M.ZC[np.clip(iz, 0, len(M.ZC) - 1)]) < 0.01) & (iz >= 0) & (iz < len(M.ZC))
    Sk = np.zeros((M.NJK, len(M.ZC), 5)); Sw = np.zeros_like(Sk)
    for ci, c in enumerate(cuts):
        m = okz & (s['RMAG'] < c)
        np.add.at(Sk, (rg_s[m], iz[m], ci), (u * k)[m])
        np.add.at(Sw, (rg_s[m], iz[m], ci), u[m])
    eta_jk = np.zeros((M.NJK, len(M.ZC)))
    for g in range(M.NJK):
        mu = (Sk.sum(0) - Sk[g]) / (Sw.sum(0) - Sw[g])
        eta_jk[g] = [1 + 2.1715 * np.polyfit(cuts, mu[i], 1)[0] - 1 + A_RESP * float(e['clip']) for i in range(len(M.ZC))]
    print('eta jackknife done %.0fs' % (time.time() - t0), flush=True)

    # s(z) and dln n/dr for C(z)
    zg = np.linspace(0, 0.7, 1401)
    rg = np.array([T.background(x)[0] for x in zg])
    nz = np.loadtxt(D + 'BGS_BRIGHT_NGC_nz.txt')
    r_nz = np.interp(nz[:, 0], zg, rg)
    dlnn = np.gradient(np.log(np.clip(nz[:, 3], 1e-8, None)), r_nz)
    ssl = []
    for z0 in M.ZC:
        m = (s['Z'] > z0 - 0.01) & (s['Z'] < z0 + 0.01)
        ne = (m & (s['RMAG'] > 19.4) & (s['RMAG'] < 19.5)).sum() / 0.1
        ssl.append(ne / m.sum() / np.log(10))
    ssl = np.array(ssl)
    Cz = np.array([T.c_rel(z0, sv, np.interp(z0, nz[:, 0], dlnn)) for z0, sv in zip(M.ZC, ssl)])
    out['s_z'] = ssl.tolist()
    out['C_z'] = Cz.tolist()

    # mock-calibrated b_s leak and A_D recovery
    mockres = []
    for mk in mocks:
        msh = M.Shells('mock%d' % mk)
        Mo = np.load(D + 'mock_%d.npz' % mk)
        mz = msh.zk
        mreg = M.jk_regions(Mo['gal_RA'][msh.kidx], Mo['gal_DEC'][msh.kidx], np.ones(len(msh.kidx), bool))
        mo = np.argsort(mreg, kind='stable')
        mb = np.searchsorted(mreg[mo], np.arange(M.NJK + 1))
        KD = Mo['gal_KD'][msh.kidx].astype(float)
        KB = Mo['gal_KB'][msh.kidx].astype(float)
        etam = np.interp(mz, M.ZC, eta_z)
        one = np.ones(len(mz))
        res = {}
        for name, mark in (('doppler', etam * KD), ('bs_unit', KB), ('one', one)):
            for bi, (lo, hi) in enumerate(BINS):
                w = ((mz > lo) & (mz < hi)).astype(float)
                num, den = region_sums(msh, mo, mb, w * mark, w)
                tot, jk = msh.multipoles(num, den, 1 if name != 'one' else 0)
                res[(name, bi)] = (tot, jk)
        small = {}
        for bi, (lo, hi) in enumerate(BINS):
            Bm = (mz > lo) & (mz < hi)
            zm = float(np.mean(mz[Bm]))
            r_, H_, f_, D_ = T.background(zm)
            xim_ = T.xi_matter_mono(zm, msh.dcen)
            bsel_ = (msh.dcen > 20) & (msh.dcen < 40)
            rm = np.mean(res[('one', bi)][0][bsel_] / xim_[bsel_])
            bm = (-2 * f_ / 3 + np.sqrt((2 * f_ / 3) ** 2 - 4 * (f_ ** 2 / 5 - rm))) / 2
            tm = M.bin_template(msh, np.where(Bm)[0], np.ones(Bm.sum()), eta_z, bm, mz[Bm]) * f1
            y, jk = res[('doppler', bi)]
            a_, ae_, _, _ = M.gls(y[FIT], tm[FIT], M.jk_cov(jk)[FIT, FIT], 13)
            yb, jkb = res[('bs_unit', bi)]
            ab_, abe_, _, _ = M.gls(yb[FIT], tm[FIT], M.jk_cov(jkb)[FIT, FIT], 13)
            r30 = np.mean(res[('one', bi)][0][(msh.dcen > 30) & (msh.dcen < 80)] / xim_[(msh.dcen > 30) & (msh.dcen < 80)])
            small[bi] = {'leak_unit': yb, 'rec': {'b_mock': float(bm), 'b_true': float(1.34 / D_),
                         'A_doppler': float(a_), 'err': float(ae_), 'A_velocity_free_unit_bs': float(ab_),
                         'err0': float(abe_)}}
        mockres.append(small)
        del msh, res, Mo
        print('mock %d done %.0fs' % (mk, time.time() - t0), flush=True)

    for bi, (lo, hi) in enumerate(BINS):
        B = gj & (zj > lo) & (zj < hi)
        Ball = (zj > lo) & (zj < hi)
        w = np.where(B, uj, 0.0)
        num, den = region_sums(sh, order, bounds, w * kj, w)
        xi = {l: sh.multipoles(num, den, l) for l in range(4)}
        ng, dg = region_sums(sh, order, bounds, Ball.astype(float), Ball.astype(float))
        gg0, gg0jk = sh.multipoles(ng, dg, 0)
        gg2, _ = sh.multipoles(ng, dg, 2)
        zeff = np.average(zj[B], weights=w[B])
        r, Hc, f, Dz = T.background(zeff)
        xim = T.xi_matter_mono(zeff, sh.dcen)
        bsel = (sh.dcen > 20) & (sh.dcen < 40)  # mocks: the monopole is unbiased only below ~45 Mpc/h (results F3)
        # b from the Kaiser monopole
        rat = np.mean(gg0[bsel] / xim[bsel])
        b = (-2 * f / 3 + np.sqrt((2 * f / 3) ** 2 - 4 * (f ** 2 / 5 - rat))) / 2
        # b_s from the monopole ratio (per separation, and its 30-80 mean)
        kais = (b ** 2 + 2 * b * f / 3 + f ** 2 / 5) / (b + f / 3)
        bs_d = xi[0][0] / gg0 * kais
        bs_jk = np.array([xi[0][1][g][bsel].sum() / gg0jk[g][bsel].sum() * kais for g in range(M.NJK)])
        bs = float(xi[0][0][bsel].sum() / gg0[bsel].sum() * kais)
        bs_err = float(np.sqrt((M.NJK - 1) * np.var(bs_jk)))
        # templates
        tD = M.bin_template(sh, np.where(B)[0], uj[B], eta_z, b, zj[B]) * f1
        tD1 = M.bin_template(sh, np.where(B)[0], uj[B], np.ones_like(eta_z), b, zj[B]) * f1
        tO = M.bin_template(sh, np.where(B)[0], uj[B], eta_z, b, zj[B], kind='oct') * f3
        Cb = float(np.average(np.interp(zj[B], M.ZC, Cz), weights=w[B]))
        tRel = -Cb * bs / ((1 / (r * Hc) - 1) * (b + 3 * f / 5)) * tD1
        # mock leak per unit b_s and recovery
        leak, leak_err, rec = np.zeros(len(sh.dcen)), np.zeros(len(sh.dcen)), []
        if mockres:
            Lm = np.array([mr[bi]['leak_unit'] for mr in mockres])
            leak = bs * Lm.mean(0)
            leak_err = np.abs(bs) * Lm.std(0) + np.abs(bs_err * Lm.mean(0))
            rec = [mr[bi]['rec'] for mr in mockres]
        # covariance: shuffles (noise) + analytic cosmic variance + leak error
        rng = np.random.default_rng(11)
        idxB = np.where(B)[0]
        slices = np.round((zj[idxB] - 0.15) / 0.01).astype(int) * 2 + capj[idxB]
        shuf = []
        denB = den.sum(0)
        QB = sh.Q[idxB]
        for it in range(NSHUF):
            kk = kj[idxB].copy()
            ww = w[idxB].copy()
            for sl in np.unique(slices):
                m = np.where(slices == sl)[0]
                p = rng.permutation(m)
                kk[m], ww[m] = kk[p], ww[p]
            numS = np.tensordot(ww * kk, QB, axes=(0, 0))
            shuf.append((numS / denB) @ sh.L[1])
        shuf = np.array(shuf)
        Cs = np.cov(shuf.T)
        Nk = B.sum()
        nbar = float(np.average(np.interp(zj[B], nz[:, 0], nz[:, 3]), weights=w[B]))
        V = 1.0
        _, Ccv = T.covariance(zeff, b, sh.dcen, Nk, nbar, 1.0, 0.42)
        # Ccv is per unit volume: scale by the bin volume from the n(z) table (both caps)
        Vb = 0
        for cap in ['NGC', 'SGC']:
            a_ = np.loadtxt(D + 'BGS_BRIGHT_%s_nz.txt' % cap)
            mm = (a_[:, 0] > lo) & (a_[:, 0] < hi)
            Vb += a_[mm, 5].sum()
        eta_b = float(np.average(np.interp(zj[B], M.ZC, eta_z), weights=w[B]))
        Ccv = Ccv / Vb * eta_b ** 2
        Cjk = M.jk_cov(xi[1][1])
        rel_err = np.abs(tRel) * np.hypot(C_FRAC, bs_err / max(abs(bs), 1e-9))
        Ctot = Cs + Ccv + np.outer(leak_err, leak_err) + np.outer(rel_err, rel_err)
        y = xi[1][0] - tRel - leak
        # eta error: jackknife of eta folded into the template amplitude
        eta_b_jk = np.array([np.average(np.interp(zj[B], M.ZC, eta_jk[g]), weights=w[B]) for g in range(M.NJK)])
        eta_err_stat = float(np.sqrt((M.NJK - 1) * np.var(eta_b_jk)))
        eta_err_sys = float(np.sqrt((A_RESP_ERR * 0.9966) ** 2 + ETA_ZSUCC ** 2 + (0.5 * np.average(
            np.abs(np.interp(zj[B], M.ZC, e['alt']) - np.interp(zj[B], M.ZC, e['alt2'])), weights=w[B])) ** 2))
        A, Ae, chi2, chi0 = M.gls(y[FIT], tD[FIT], Ctot[FIT, FIT], 13, NSHUF)
        A30, Ae30, _, _ = M.gls(y[FIT30], tD[FIT30], Ctot[FIT30, FIT30], 12, NSHUF)
        Ajk, Aejk, _, _ = M.gls(y[FIT], tD[FIT], Cjk[FIT, FIT], 13)
        # amplitude of the raw dipole in eta=1 units (direct test against zero, independent of eta)
        amp, ampe, _, chi0z = M.gls(y[FIT], tD1[FIT], Ctot[FIT, FIT], 13, NSHUF)
        # b from 20-50 Mpc/h (sensitivity to large-scale power in the counts)
        bsel2 = (sh.dcen > 30) & (sh.dcen < 80)  # the pre-registered range, as a sensitivity
        rat2 = np.mean(gg0[bsel2] / xim[bsel2])
        b2 = (-2 * f / 3 + np.sqrt((2 * f / 3) ** 2 - 4 * (f ** 2 / 5 - rat2))) / 2
        tD_b2 = M.bin_template(sh, np.where(B)[0], uj[B], eta_z, b2, zj[B]) * f1
        A_b2, _, _, _ = M.gls(y[FIT], tD_b2[FIT], Ctot[FIT, FIT], 13, NSHUF)
        # shuffle-based p-value for the zero-signal chi2
        Ci = np.linalg.inv(Ctot[FIT, FIT]) * (NSHUF - 15) / (NSHUF - 1)
        chi_sh = np.array([v[FIT] @ Ci @ v[FIT] for v in (shuf - shuf.mean(0))])
        A_eta_err = abs(A) * np.hypot(eta_err_stat, eta_err_sys) / max(abs(eta_b), 1e-9)
        res = {'z': [lo, hi], 'N': int(Nk), 'zeff': float(zeff), 'b': float(b), 'f': float(f), 'eta': eta_b,
               'eta_err_stat': eta_err_stat, 'eta_err_sys': eta_err_sys,
               'b_s': bs, 'b_s_err': bs_err, 'b_s_by_d': bs_d.tolist(), 'C': Cb,
               'xi1': xi[1][0].tolist(), 'xi1_err_tot': np.sqrt(np.diag(Ctot)).tolist(),
               'xi1_err_jk': np.sqrt(np.diag(Cjk)).tolist(), 'xi1_err_shuffle': np.sqrt(np.diag(Cs)).tolist(),
               'xi0': xi[0][0].tolist(), 'xi2': xi[2][0].tolist(), 'xi3': xi[3][0].tolist(),
               'xi3_err': np.sqrt(np.diag(M.jk_cov(xi[3][1]))).tolist(), 'tO': tO.tolist(),
               'tD': tD.tolist(), 'tRel': tRel.tolist(), 'leak': leak.tolist(), 'gg0': gg0.tolist(),
               'A_D': A, 'A_D_err_stat': Ae, 'A_D_err_eta': A_eta_err, 'chi2_fit': chi2, 'chi2_null': chi0,
               'p_null_shuffle': float(np.mean(chi_sh >= chi0)), 'A_D_d30': A30, 'A_D_d30_err': Ae30,
               'A_D_jkcov': Ajk, 'A_D_jkcov_err': Aejk, 'amp_eta1': amp, 'amp_eta1_err': ampe,
               'mock_recovery': rec, 'b_30_80': float(b2), 'A_D_b_30_80': A_b2, 'dcen': sh.dcen.tolist()}
        # octupole against prediction and zero
        Co = M.jk_cov(xi[3][1])[FIT, FIT]
        ao, aoe, _, co0 = M.gls(xi[3][0][FIT], tO[FIT], Co, 13)
        res['A_oct'] = ao
        res['A_oct_err'] = aoe
        res['chi2_oct_zero'] = co0
        out['bins'][bi] = res
        print('bin %.2f-%.2f: N %d zeff %.3f b %.2f b_s %.4f+-%.4f C %.2f eta %.3f | A_D = %.2f +- %.2f (stat) +- %.2f '
              '(eta) | chi2 %.1f/12, null chi2 %.1f (p_shuffle %.3f) | d>=30: %.2f +- %.2f | jk cov: %.2f +- %.2f | '
              'oct A %.2f +- %.2f | mocks %s  [%.0fs]' % (
                  lo, hi, Nk, zeff, b, bs, bs_err, Cb, eta_b, A, Ae, A_eta_err, chi2, chi0, res['p_null_shuffle'],
                  A30, Ae30, Ajk, Aejk, ao, aoe, rec, time.time() - t0), flush=True)
    json.dump(out, open(D + ('run.json' if tag == 'data' and '--blind' not in sys.argv else 'run_%s.json' % tag), 'w'),
              indent=1, default=float)


if __name__ == '__main__':
    main()
