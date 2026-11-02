"""Injection calibration of the CMB-lensing cluster stack (R-154 method check, points 1, 2 and 8).

The release's reconstruction pipeline is not public as a runnable map-level tool, so this emulates it in flat sky. The
emulated pipeline:
- a 256' box at 0.5' pixels;
- Gaussian unlensed CMB from CAMB in the ACT catalogue cosmology;
- one cluster lensing the CMB by remapping, with NFW kappa from act_szlens's halo model at M_lens = R x M_SZ;
- the cluster's tSZ from the A10 Y500-M relation at M_SZ (UPP shape), at 90 and 150 GHz;
- ACT-like beams (2.1' / 1.4' FWHM) and noise (white plus atmospheric 1/f in temperature, matched so that the emulated
  estimator's N_0 is the release's tonly N_L);
- optionally, the release's treatment of SNR > 5 clusters: a UPP template, at the nearest of 15 Nemo filter scales, with its
  amplitude fitted by a matched filter on each frequency map, subtracted before the coadd;
- the inverse-variance coadd;
- the temperature quadratic estimator (Hu & Okamoto), profile-hardened (Sailer et al. 2020) against a tSZ-shaped source, with
  CMB 600 < l < 3000 as in DR6;
- the reconstruction filtered as the release map is used in the stack (20 <= L <= 2048, taper), with annulus means in the
  stack's annuli inside 30'.

Paired realisations (the same CMB and noise with and without the cluster terms) isolate each term. The outputs per (z, M_SZ, R):
- the lensing response: the stacked statistic's amplitude against the input kappa's, which is 1 for a linear, normalised
  estimator;
- the tSZ leakage: in units of the input kappa amplitude, with and without template subtraction, and for a true pressure
  profile that differs from the template (the Planck 2013 GNFW);
- the CIB leakage: a beam-scale emission term at the cluster with a given 150 GHz flux.

The release's baseline is the minimum-variance combination with polarization, which carries no tSZ; temperature leakage
enters it with the weight w_T(L) = N_L^MV / N_L^T, from the release N_L curves.

usage: python -I act_szlens_inject.py calib | run [--nreal 40] | cib --flux ... (see main)
"""
import argparse, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import act_szlens as asl   # noqa: E402
import erass1_act as ea    # noqa: E402

OUT = os.path.join(asl.OUT, 'inject')
ARCMIN = np.pi / 180 / 60
NPIX, PIX = 512, 0.5 * ARCMIN
TCMB = 2.7255e6  # uK
LCMB_MIN, LCMB_MAX = 600, 3000
FREQS = {'f090': dict(nu=97.0, fwhm=2.1), 'f150': dict(nu=149.0, fwhm=1.4)}
log = ea.log


def tsz_cl_profile():
    """sqrt of the one-halo tSZ power spectrum, C_l = int dz dV/dz/dOmega int dlnM dn/dlnM |y_l|^2 (Tinker 2008, A10 Y-M at
    1-b = 0.8, UPP), as a radial function of l; cached."""
    fn = os.path.join(OUT, 'tsz_cl.npz')
    if os.path.exists(fn):
        g = np.load(fn); return g['l'], g['u']
    import erass1_ysz as r131
    from scipy.special import j0
    hm = asl.ACTHaloModel(lmax=3000)
    up = r131.UPP()
    lr = np.arange(0, 20001, 10.0)
    zg = np.linspace(0.05, 3.0, 30); lnM = np.linspace(np.log(5e13), np.log(2e15), 18)
    Cl = np.zeros(len(lr))
    co = hm.co
    for z in zg:
        dn = asl.tinker08_dndlnm(hm, lnM, z)
        DA = co.angular_diameter_distance(z).value
        dV = (co.differential_comoving_volume(z).value)            # Mpc^3 / sr
        for M, n in zip(np.exp(lnM), dn):
            M_sz = 0.8 * M
            Y500 = 10 ** -4.739 * co.efunc(z) ** (2 / 3) * (M_sz / 3e14) ** 1.79
            r500 = (3 * M / (4 * np.pi * 500 * co.critical_density(z).to('Msun/Mpc3').value)) ** (1 / 3)
            t500 = r500 / DA; y0 = Y500 / (r500 ** 2 * up.S)
            th = np.logspace(np.log10(0.01 * t500), np.log10(5 * t500), 300)
            u = up.profile(th, t500) * y0
            yl = np.array([2 * np.pi * np.trapezoid(u * j0(l * th) * th, th) for l in lr])
            Cl += yl ** 2 * n * dV * (lnM[1] - lnM[0]) * (zg[1] - zg[0])
    u = np.sqrt(Cl); u /= u[0]
    os.makedirs(OUT, exist_ok=True); np.savez(fn, l=lr, u=u)
    return lr, u


def gnu(nu):
    x = nu / 56.78
    return x / np.tanh(x / 2) - 4


class Box:
    def __init__(self):
        lx = 2 * np.pi * np.fft.fftfreq(NPIX, d=PIX)
        self.lx, self.ly = np.meshgrid(lx, lx, indexing='xy')
        self.l = np.hypot(self.lx, self.ly)
        x = (np.arange(NPIX) - NPIX // 2) * PIX
        self.X, self.Y = np.meshgrid(x, x, indexing='xy')
        self.r = np.hypot(self.X, self.Y)
        self.area = (NPIX * PIX) ** 2

    def fft(self, m):
        return np.fft.fft2(m) * PIX ** 2

    def ifft(self, a):
        return np.real(np.fft.ifft2(a)) / PIX ** 2

    def ifftc(self, a):
        return np.fft.ifft2(a) / PIX ** 2

    def radial_to_2d(self, ell, f):
        return np.interp(self.l, ell, f, right=0.0)

    def gauss_field(self, rng, cl2d):
        """real Gaussian field with power cl2d (continuous convention: <|a|^2> = area * C_l)."""
        w = rng.standard_normal((NPIX, NPIX))
        return self.ifft(self.fft(w) / PIX * np.sqrt(cl2d))


def camb_cls(lmax=6000):
    import camb
    h = asl.H0_ACT / 100
    pars = camb.set_params(H0=asl.H0_ACT, ombh2=asl.OB_ACT * h ** 2, omch2=(asl.OM_ACT - asl.OB_ACT) * h ** 2, As=2.1e-9, ns=0.965,
                           lmax=lmax, lens_potential_accuracy=1)
    res = camb.get_results(pars)
    cl = res.get_unlensed_scalar_cls(CMB_unit='muK', raw_cl=True, lmax=lmax)[:, 0]
    cll = res.get_lensed_scalar_cls(CMB_unit='muK', raw_cl=True, lmax=lmax)[:, 0]
    return np.arange(lmax + 1), cl, cll


class Pipeline:
    def __init__(self, noise_scale=1.0, white=(20.0, 15.0), lknee=2500.0, alpha=-3.5, hard_t500=3.0):
        self.b = Box()
        ell, cl, cll = camb_cls()
        self.ell = ell
        self.C = self.b.radial_to_2d(ell, cl); self.C[0, 0] = 0
        self.Cl = self.b.radial_to_2d(ell, cll)
        self.beam = {}; self.N = {}
        for (k, f), w in zip(FREQS.items(), white):
            s = f['fwhm'] * ARCMIN / np.sqrt(8 * np.log(2))
            self.beam[k] = np.exp(-0.5 * self.b.l ** 2 * s ** 2)
            nw = (noise_scale * w * ARCMIN) ** 2
            with np.errstate(divide='ignore'):
                self.N[k] = nw * (1 + (np.maximum(self.b.l, 1) / lknee) ** alpha)   # noise on the sky, before beam deconvolution
        # inverse-variance coadd of beam-deconvolved maps: weights w_f = (B_f^2 / N_f) / sum
        inv = {k: self.beam[k] ** 2 / self.N[k] for k in FREQS}
        tot = sum(inv.values())
        self.w = {k: inv[k] / tot for k in FREQS}
        self.Ncoadd = 1 / tot
        self.Ctot = self.Cl + self.Ncoadd
        sel = (self.b.l >= LCMB_MIN) & (self.b.l <= LCMB_MAX)
        self.F = np.where(sel, 1 / self.Ctot, 0.0)
        # hardening profile: the tSZ shape of a typical cluster in the coadd (a UPP at theta500 = 3', beam-convolved, f150)
        if hard_t500 > 0:
            self.u = self.upp_ell(hard_t500 * ARCMIN) * self.beam['f150']
        else:   # the release's choice: u(l) = sqrt(C_l^tSZ), the one-halo tSZ power (Qu et al. 2024 V.8.2; here a Tinker x A10 halo model)
            self.u = self.b.radial_to_2d(*tsz_cl_profile()) * self.beam['f150']
        self.u /= self.u[0, 0] if self.u[0, 0] else 1
        self._responses()
        self.cmb_rng = None

    # ------------------------------------------------------------------------------------------------- profiles
    def upp_ell(self, theta500, gnfw=None):
        """2d Fourier transform of the projected UPP (unit amplitude at x = 0 in u(x)) for a cluster of size theta500."""
        import erass1_ysz as r131
        if not hasattr(self, '_upp'):
            self._upp = r131.UPP()
        up = self._upp
        if gnfw is not None:
            up = self._gnfw_upp(gnfw)
        th = np.logspace(np.log10(0.01 * ARCMIN), np.log10(5 * theta500 * 1.001), 800)
        u = up.profile(th, theta500)
        from scipy.special import j0
        lr = np.arange(0, 20001, 5.0)
        ft = np.array([2 * np.pi * np.trapezoid(u * j0(l * th) * th, th) for l in lr])
        return self.b.radial_to_2d(lr, ft)

    def _gnfw_upp(self, pars):
        import erass1_ysz as r131
        key = tuple(pars)
        if not hasattr(self, '_gn'):
            self._gn = {}
        if key not in self._gn:
            class G(r131.UPP):
                P0, c500, gam, alp, bet = pars
            self._gn[key] = G()
        return self._gn[key]

    # ------------------------------------------------------------------------------------------------- estimator
    def _responses(self):
        b = self.b; F = self.F; C = self.C; u = self.u
        lx, ly = b.lx, b.ly
        # R_phiphi(L) = L_i L_j FT[a_ij b] + L_i L_j FT[c_i c_j]
        bb = b.ifft(F)
        a = {('x', 'x'): b.ifft(lx * lx * C ** 2 * F), ('x', 'y'): b.ifft(lx * ly * C ** 2 * F), ('y', 'y'): b.ifft(ly * ly * C ** 2 * F)}
        c = {'x': b.ifftc(1j * lx * C * F), 'y': b.ifftc(1j * ly * C * F)}
        L = {'x': lx, 'y': ly}
        R = np.zeros_like(lx)
        for (i, j), aij in a.items():
            m = 1 if i == j else 2
            R = R + m * L[i] * L[j] * b.fft(aij * bb)
        for i in 'xy':
            for j in 'xy':
                # c_i c_j with c = IFT[i l C F]: product gives -(l_i C F)(l'_j C F) convolution
                R = R - L[i] * L[j] * b.fft(c[i] * c[j])
        self.Rpp = np.real(R)
        # R_phis(L) = L_i FT[d_i e], d_i = IFT[l_i C u F], e = IFT[u F]
        e = b.ifft(u * F)
        Rps = np.zeros_like(lx, dtype=complex)
        for i in 'xy':
            d = b.ifftc(1j * L[i] * C * u * F)
            Rps = Rps + (-1j) * L[i] * b.fft(d * e)
        self.Rps = np.real(Rps)
        e2 = b.ifft(u ** 2 * F)
        self.Rss = 0.5 * np.real(b.fft(e2 * e2))

    def qe(self, Tco, harden=True):
        """phi(L) from the coadded temperature map (Fourier, continuous normalisation); bias-hardened if harden."""
        b = self.b
        FT = self.F * Tco
        h = b.ifft(FT)
        gx = b.ifft(1j * b.lx * self.C * FT); gy = b.ifft(1j * b.ly * self.C * FT)
        phib = -1j * (b.lx * b.fft(h * gx) + b.ly * b.fft(h * gy))
        if not harden:
            with np.errstate(invalid='ignore', divide='ignore'):
                return np.where(self.Rpp > 0, phib / self.Rpp, 0)
        sb = 0.5 * b.fft(b.ifft(self.u * FT) ** 2)
        with np.errstate(invalid='ignore', divide='ignore'):
            num = phib - self.Rps / self.Rss * sb
            den = self.Rpp - self.Rps ** 2 / self.Rss
            return np.where(den > 0, num / den, 0)

    # ------------------------------------------------------------------------------------------------- sky
    def lens(self, T, kappa_ell):
        """remap T (real space) by the deflection of kappa (Fourier): T(x + grad phi)."""
        from scipy.ndimage import map_coordinates
        b = self.b
        with np.errstate(invalid='ignore', divide='ignore'):
            phi = np.where(b.l > 0, 2 * kappa_ell / b.l ** 2, 0)
        ax = b.ifft(1j * b.lx * phi); ay = b.ifft(1j * b.ly * phi)
        # pixel coordinates (row = y, col = x)
        jj, ii = np.meshgrid(np.arange(NPIX), np.arange(NPIX), indexing='xy')
        coords = [ii + ay / PIX, jj + ax / PIX]
        return map_coordinates(T, coords, order=5, mode='wrap')

    def sky(self, rng, kappa_ell=None, y_ell=None, cib_ell=None, subtract=None, gnfw_y_ell=None):
        """the coadded temperature (Fourier) for one realisation; rng streams fixed by the caller for pairing.
        y_ell: the cluster's Compton-y (Fourier, centred); cib_ell: dict f -> the CIB term in uK (Fourier, centred);
        subtract: None or a dict(theta500_template, ...) for the template subtraction."""
        b = self.b
        st = rng.bit_generator.state
        T = b.gauss_field(rng, self.C)
        if kappa_ell is not None:
            T = self.lens(T, kappa_ell)
        Tl = b.fft(T)
        out = 0
        for k, f in FREQS.items():
            m = Tl.copy()
            if y_ell is not None:
                m = m + gnu(f['nu']) * TCMB * y_ell
            if cib_ell is not None:
                m = m + cib_ell[k]
            m = m * self.beam[k]
            nrng = np.random.default_rng(rng.integers(1 << 62))
            m = m + b.fft(b.gauss_field(nrng, self.N[k] * self.beam[k] ** 2 / self.beam[k] ** 2))
            if subtract is not None:
                m = m - self.template_fit(m, k, subtract)
            with np.errstate(invalid='ignore', divide='ignore'):
                out = out + self.w[k] * np.where(self.beam[k] > 1e-8, m / self.beam[k], 0)
        return out

    def template_fit(self, m, k, subtract):
        """the release's subtraction, emulated: a beam-convolved UPP at the nearest Nemo filter scale, its amplitude from a
        matched filter on this frequency map (noise = lensed CMB + instrument), subtracted (returned in the map's units)."""
        tpl = subtract['tpl'] * gnu(FREQS[k]['nu']) * TCMB * self.beam[k]          # unit-y template at this frequency
        Nf = (self.Cl * self.beam[k] ** 2 + self.N[k])
        sel = self.b.l > 100
        num = np.sum(np.real(np.conj(tpl) * m / Nf)[sel]); den = np.sum((np.abs(tpl) ** 2 / Nf)[sel])
        return (num / den) * tpl

    # ------------------------------------------------------------------------------------------------- statistic
    def stat_setup(self, ell_band, filt_in):
        """the stack's filter (release band, nside-1024 pixel window) on a 2d grid, and the annulus masks inside 30'."""
        self.filt2 = self.b.radial_to_2d(ell_band, filt_in)
        edges = ea.EDGES_IN * ARCMIN
        rr = self.b.r
        self.ann = [(rr >= edges[i]) & (rr < edges[i + 1]) for i in range(len(edges) - 1)]
        # centre the cluster at the box centre: multiply Fourier by a phase so that the field origin maps to the centre
        c = (NPIX // 2) * PIX
        self.shift = np.exp(-1j * (self.b.lx * c + self.b.ly * c))

    def annuli(self, kappa_ell):
        m = self.b.ifft(kappa_ell * self.filt2)
        return np.array([m[a].mean() for a in self.ann])


def cluster_terms(pl, hm, M_sz, z, R, gnfw=None):
    """the cluster's input kappa (Fourier, centred) at M_lens = R M_sz, and its y (Fourier, centred) from A10 at M_sz."""
    from scipy.special import j0
    k1, k2, info = hm.kappa_theta(R * M_sz, z)
    th = hm.theta
    lr = np.arange(0, 20001, 5.0)
    # 1h only, truncated at the theta grid edge (6 deg); Hankel transform to the 2d grid
    ft = np.array([2 * np.pi * np.trapezoid(k1 * j0(l * th) * th, th) for l in lr])
    kap = pl.b.radial_to_2d(lr, ft) * pl.shift
    # A10: E(z)^(-2/3) Y500 D_A^2 = 10^-4.739 (M/3e14)^1.79 Mpc^2 (h70 = 1)
    co = hm.co
    Ez = co.efunc(z); DA = co.angular_diameter_distance(z).value
    Y500 = 10 ** -4.739 * Ez ** (2 / 3) * (M_sz / 3e14) ** 1.79           # Mpc^2
    r500 = (3 * M_sz / (4 * np.pi * 500 * co.critical_density(z).to('Msun/Mpc3').value)) ** (1 / 3)
    t500 = r500 / DA
    import erass1_ysz as r131
    up = r131.UPP() if gnfw is None else pl._gnfw_upp(gnfw)
    y0 = Y500 / (r500 ** 2 * up.S)          # y(theta) = y0 u(theta/theta500)
    yl = pl.upp_ell(t500, gnfw) * y0 * pl.shift
    return kap, yl, t500, info


def nemo_scale(t500):
    """the nearest of 15 Nemo filter scales (theta500 from 0.5' to 8', log-spaced) to the cluster's own."""
    grid = np.logspace(np.log10(0.5), np.log10(8.0), 15) * ARCMIN
    return grid[np.argmin(np.abs(np.log(grid / t500)))]


def weights_cov(pl, hm):
    """the per-cluster covariance of the 5 inner annulus means (release N_L + C_kk), to weight the amplitude fit as the real fit."""
    nl = np.loadtxt(os.path.join(ea.ACT, 'N_L_kk_act_dr6_lensing_v1_baseline.txt'))
    ell = hm.ell
    f_in, _ = ea.filters(hm, ea.LMIN, ea.LMAX)
    W = asl.annulus_windows(ell, f_in, ea.EDGES_IN)
    C = asl.single_cov(ell, np.interp(ell, nl[:, 0], nl[:, 1]), W)
    return np.linalg.inv(C), f_in


def amp(a_out, a_in, Ci):
    return (a_in @ Ci @ a_out) / (a_in @ Ci @ a_in)


def coadd_y(pl, y_ell):
    """the tSZ of a Compton-y map in the beam-deconvolved coadd (uK, Fourier): sum_f w_f g_f T_CMB y."""
    return sum(pl.w[k] * gnu(f['nu']) * TCMB * y_ell for k, f in FREQS.items())


def run(args):
    """per (z, M_SZ): the deterministic tSZ leakage QE(s, s) (UPP truth and Planck-2013 GNFW truth; hardened and not); the
    template-subtracted leakage by Monte Carlo (the residual depends on the noisy amplitude fit); per R, the lensing response
    by paired Monte Carlo. All in units of the input kappa's stacked amplitude at that R."""
    os.makedirs(OUT, exist_ok=True)
    hm = asl.ACTHaloModel(lmax=3000)
    pl = Pipeline(noise_scale=args.noise_scale, hard_t500=args.hard_t500)
    Ci, f_in = weights_cov(pl, hm)
    pl.stat_setup(hm.ell, f_in)
    gnfw_p13 = (6.41, 1.81, 0.31, 1.33, 4.13)   # Planck 2013 GNFW: a true pressure profile unlike the A10 template
    k_of = lambda T, h=True: pl.b.l ** 2 / 2 * pl.qe(T, harden=h)
    rows = []
    cases = [(0.33, 1.96e14), (0.64, 2.08e14), (0.89, 1.93e14), (1.23, 1.76e14), (0.64, 4.5e14), (1.23, 4.0e14)]
    if args.quick:
        cases = cases[:1]
    for z, Msz in cases:
        t0 = time.time()
        kap1, yl, t500, _ = cluster_terms(pl, hm, Msz, z, 1.0)
        _, yl13, _, _ = cluster_terms(pl, hm, Msz, z, 1.0, gnfw=gnfw_p13)
        a1 = pl.annuli(kap1)
        det = {}
        for nm, y in [('upp', yl), ('p13', yl13)]:
            sco = coadd_y(pl, y)
            det[nm + '_h'] = amp(pl.annuli(k_of(sco)), a1, Ci)
            det[nm + '_noh'] = amp(pl.annuli(k_of(sco, False)), a1, Ci)
        tpl = pl.upp_ell(nemo_scale(t500)) * pl.shift
        sub = dict(tpl=tpl)
        # template-subtracted: <QE(T + r) - QE(T)> with r = s - a_hat tpl, a_hat fitted per realisation (Monte Carlo)
        mc = {'upp_sub': [], 'p13_sub': []}
        for i in range(args.nsub):
            seed = 1000 * i + 11
            T0 = pl.sky(np.random.default_rng(seed))
            k0 = k_of(T0)
            for nm, y in [('upp_sub', yl), ('p13_sub', yl13)]:
                Ts = pl.sky(np.random.default_rng(seed), y_ell=y, subtract=sub)
                Tn = pl.sky(np.random.default_rng(seed), subtract=sub)   # the subtraction's own effect without a cluster
                mc[nm].append(amp(pl.annuli(k_of(Ts) - k_of(Tn)), a1, Ci))
        for nm, v in mc.items():
            det[nm] = np.mean(v); det[nm + '_err'] = np.std(v) / np.sqrt(len(v))
        log(f'z {z} M_SZ {Msz:.2e}: theta500 {t500 / ARCMIN:.2f}\', template {nemo_scale(t500) / ARCMIN:.2f}\'; tSZ leakage / kappa(R=1) amplitude: '
            + ', '.join(f'{k} {v:+.3f}' for k, v in det.items()) + f' ({time.time() - t0:.0f}s)')
        for R in args.R:
            t0 = time.time()
            kap, _, _, _ = cluster_terms(pl, hm, Msz, z, R)
            aR = pl.annuli(kap)
            lens = []; lens_n = []
            for i in range(args.nreal):
                seed = 1000 * i + 7
                T0 = pl.sky(np.random.default_rng(seed)); TL = pl.sky(np.random.default_rng(seed), kappa_ell=kap)
                lens.append(amp(pl.annuli(k_of(TL) - k_of(T0)), aR, Ci))
                lens_n.append(amp(pl.annuli(k_of(TL, False) - k_of(T0, False)), aR, Ci))
            row = dict(z=z, Msz=Msz, R=R, t500=t500 / ARCMIN, aR_over_a1=amp(aR, a1, Ci), lens=np.mean(lens), lens_err=np.std(lens) / np.sqrt(len(lens)),
                       lens_noh=np.mean(lens_n), lens_noh_err=np.std(lens_n) / np.sqrt(len(lens_n)), **det)
            rows.append(row)
            log(f'   R {R}: kappa amplitude vs R=1 {row["aR_over_a1"]:.2f}; lensing response hardened {row["lens"]:.3f}±{row["lens_err"]:.3f}, '
                f'unhardened {row["lens_noh"]:.3f}±{row["lens_noh_err"]:.3f} ({time.time() - t0:.0f}s, {args.nreal} pairs)')
        np.save(os.path.join(OUT, f'inject_rows{args.tag}.npy'), np.array(rows, dtype=object), allow_pickle=True)


def template_fit_det(pl, m_ell, k, tpl):
    """the subtraction's amplitude on a noise-free map (the deterministic part of the matched-filter fit)."""
    t = tpl * gnu(FREQS[k]['nu']) * TCMB * pl.beam[k]
    Nf = pl.Cl * pl.beam[k] ** 2 + pl.N[k]
    sel = pl.b.l > 100
    return np.sum(np.real(np.conj(t) * m_ell * pl.beam[k] / Nf)[sel]) / np.sum((np.abs(t) ** 2 / Nf)[sel]) * t


def cib(args):
    """the CIB leakage per bin from the measured Planck 857 GHz flux (act_szlens.py cib), extrapolated to 90/150 GHz with the
    modified blackbody that matches the bin's measured 545/857 ratio; injected as a point source and as a 1' Gaussian at the
    cluster centre, with the cluster's tSZ; deterministic QE terms: [QE(s + c) - QE(s)] hardened, and with the template
    subtraction (its amplitude fitted on s + c, so the CIB also biases the subtracted amplitude)."""
    hm = asl.ACTHaloModel(lmax=3000)
    pl = Pipeline(noise_scale=args.noise_scale, hard_t500=args.hard_t500)
    Ci, f_in = weights_cov(pl, hm)
    pl.stat_setup(hm.ell, f_in)
    rows = np.load(os.path.join(asl.OUT, 'cib_rows.npy'), allow_pickle=True)
    k_of = lambda T, h=True: pl.b.l ** 2 / 2 * pl.qe(T, harden=h)
    from scipy.constants import k as kB, c as cL
    def dBdT(nu):
        x = nu / 56.78; return 2 * kB * (nu * 1e9) ** 2 / cL ** 2 * x ** 2 * np.exp(x) / np.expm1(x) ** 2 * 1e26
    out = []
    for r, (z, Msz) in zip(rows, [(0.33, 1.96e14), (0.64, 2.08e14), (0.89, 1.93e14), (1.23, 1.76e14)]):
        kap1, yl, t500, _ = cluster_terms(pl, hm, Msz, z, 1.0)
        a1 = pl.annuli(kap1)
        tpl = pl.upp_ell(nemo_scale(t500)) * pl.shift
        meas = r['f545'] / r['f857']
        ext = r['ext']
        # the (beta, T) whose model 545/857 ratio is nearest the measured one; and the brackets beta 1.5 T 25, beta 2 T 20
        keys = sorted({k_.rsplit('_', 1)[0] for k_ in ext})
        best = min(keys, key=lambda k_: abs(ext[k_ + '_150'][2] - meas))
        for sedk in sorted({best, '1.5_25', '2.0_20'}):
            for shape in ('point', 'gauss1'):
                prof = np.ones_like(pl.b.l) if shape == 'point' else np.exp(-0.5 * (pl.b.l * ARCMIN / np.sqrt(8 * np.log(2))) ** 2 * 0 - 0.5 * (pl.b.l * 1.0 * ARCMIN) ** 2 / (8 * np.log(2)))
                cm = {'f150': ext[sedk + '_150'][1], 'f090': ext[sedk + '_90'][1]}     # uK sr per cluster
                c_ell = {k_: cm[k_] * prof * pl.shift for k_ in FREQS}
                s_co = coadd_y(pl, yl); c_co = sum(pl.w[k_] * c_ell[k_] for k_ in FREQS)
                lh = amp(pl.annuli(k_of(s_co + c_co) - k_of(s_co)), a1, Ci)
                lnoh = amp(pl.annuli(k_of(s_co + c_co, False) - k_of(s_co, False)), a1, Ci)
                # with subtraction: per frequency, fit the template on s + c (noise-free part), subtract, coadd
                def resid(with_c):
                    tot = 0
                    for k_, f in FREQS.items():
                        m = gnu(f['nu']) * TCMB * yl + (c_ell[k_] if with_c else 0)
                        m = m - template_fit_det(pl, m, k_, tpl) / pl.beam[k_]
                        tot = tot + pl.w[k_] * m
                    return tot
                lsub = amp(pl.annuli(k_of(resid(True)) - k_of(resid(False))), a1, Ci)
                lsub_tot = amp(pl.annuli(k_of(resid(True))), a1, Ci)
                out.append(dict(z=z, sed=sedk, shape=shape, S150_mJy=ext[sedk + '_150'][0] * 1e3, meas_ratio=meas,
                                cib_h=lh, cib_noh=lnoh, cib_sub=lsub, tsz_cib_sub_total=lsub_tot))
                log(f'z {z}: SED {sedk} (model 545/857 {ext[sedk + "_150"][2]:.2f}, measured {meas:.2f}), {shape}, S150 {ext[sedk + "_150"][0] * 1e3:.2f} mJy: '
                    f'CIB leakage / kappa amplitude hardened {lh:+.3f}, unhardened {lnoh:+.3f}, template-subtracted {lsub:+.3f} '
                    f'(tSZ+CIB residual after subtraction {lsub_tot:+.3f}); 857 S/N {r["f857"] / r["e857"]:.1f}')
    np.save(os.path.join(OUT, f'cib_inject{args.tag}.npy'), np.array(out, dtype=object), allow_pickle=True)


def calib(args):
    """the emulated estimator's response to Gaussian kappa (must be 1) and its N_0 against the release's tonly N_L."""
    hm = asl.ACTHaloModel(lmax=3000)
    pl = Pipeline(noise_scale=args.noise_scale)
    b = pl.b
    ell = np.arange(0, 6001)
    # a Gaussian kappa input, C_kk flat-ish (large) so the paired response is measured precisely
    ckk = np.where((ell > 10) & (ell < 4000), 1e-10, 0.0)   # kappa rms 0.035
    rat = []; n0 = []
    bins = np.array([50, 150, 300, 500, 800, 1200, 1600, 2000])
    for i in range(args.nreal):
        rng = np.random.default_rng(99 + i)
        kin = b.gauss_field(rng, b.radial_to_2d(ell, ckk))
        kl = b.fft(kin)
        s = 1000 * i + 3
        T0 = pl.sky(np.random.default_rng(s)); TL = pl.sky(np.random.default_rng(s), kappa_ell=kl)
        for h in (True, False):
            kr = b.l ** 2 / 2 * (pl.qe(TL, h) - pl.qe(T0, h))
            x = np.real(kr * np.conj(kl)); a = np.abs(kl) ** 2
            idx = np.digitize(b.l, bins)
            rat.append((h, [x[idx == j].sum() / a[idx == j].sum() for j in range(1, len(bins))]))
        k0 = b.l ** 2 / 2 * pl.qe(T0)
        p = np.abs(k0) ** 2 / b.area
        n0.append([p[idx == j].mean() for j in range(1, len(bins))])
    rh = np.mean([r for h, r in rat if h], 0); rn = np.mean([r for h, r in rat if not h], 0)
    nl = {v: np.loadtxt(os.path.join(ea.ACT, f'N_L_kk_act_dr6_lensing_v1_{v}.txt')) for v in ('baseline', 'tonly', 'polonly')}
    lc = 0.5 * (bins[1:] + bins[:-1])
    log(f'bins {bins}; response hardened {np.round(rh, 3)}; unhardened {np.round(rn, 3)}')
    log(f'N0 (hardened, emulated, includes C_kk-free) {np.array2string(np.mean(n0, 0), precision=2)}')
    for v, a in nl.items():
        log(f'release N_L {v} at bin centres {np.array2string(np.interp(lc, a[:, 0], a[:, 1]), precision=2)}')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('cmd')
    p.add_argument('--nreal', type=int, default=40)
    p.add_argument('--noise_scale', type=float, default=1.6)   # emulated N_0 ~ the release tonly N_L
    p.add_argument('--R', type=float, nargs='+', default=[1.35, 3.4])
    p.add_argument('--quick', action='store_true')
    p.add_argument('--nsub', type=int, default=60)
    p.add_argument('--hard_t500', type=float, default=0.0)   # arcmin: a UPP hardening profile of this size; 0 = sqrt(C_l^tSZ), the release's
    p.add_argument('--tag', default='')
    a = p.parse_args()
    {'run': run, 'calib': calib, 'cib': cib}[a.cmd](a)


if __name__ == '__main__':
    main()
