"""DESI DR1 LRG BAO per Galactic cap (sweep-2.md row 1): do the NGC and SGC give the same isotropic BAO scale alpha_iso?
DESI 2024 III (arXiv:2404.03000) says the caps' clustering is consistent and publishes only cap-combined BAO numbers.

Data: DR1 LSScats v1.5 LRG_{NGC,SGC}_clustering.dat.fits and LRG_{NGC,SGC}_{i}_clustering.ran.fits (data.desi.lbl.gov/public),
0.4 < z < 1.1 (or the DESI bins with --zbins), weights WEIGHT x WEIGHT_FKP for data and randoms, DESI fiducial cosmology
(Planck 2018, h = 0.6736, Om = 0.3153), positions in Mpc/h. xi_l(s) by pycorr (Corrfunc, Landy-Szalay, s in 4 Mpc/h bins to
200, 100 mu bins, midpoint LOS), with a jackknife over healpix nside-8 regions for the covariance (Hartlap-corrected).
BAO fit: a linear-theory template (camb P(k) at z = 0.75, de-wiggled against Eisenstein & Hu 1998 no-wiggle, BAO damped with
Sigma_nl = 7 Mpc/h), xi_model(s) = B xi_t(alpha s) + a0 + a1/s + a2/s^2 on the monopole over 50 < s < 150 Mpc/h; alpha scanned
0.8-1.2 with the linear parameters solved analytically; sigma from Delta chi2 = 1; detection from the no-wiggle template's
chi2. Then alpha_NGC - alpha_SGC in sigma, and the caps' mean directions against the CMB apex (the row's second question).
Pair counts are cached in data/desi_lss/xi_<tracer>_<cap>_<zmin>_<zmax>.npy, so the fit can be rerun without them.
usage: <venv>/bin/python desi_bao_caps.py data/desi_lss [--tracer LRG] [--nran 1] [--rmult 4] [--zbins] [--nthreads 4] [--fit-only]
venv: python -m venv v && v/bin/pip install numpy scipy astropy healpy camb git+https://github.com/cosmodesi/Corrfunc@desi
      git+https://github.com/cosmodesi/pycorr (needs libgsl-dev; PyPI Corrfunc lacks the PIP weights pycorr calls)
"""
import sys, os, argparse, time, numpy as np
from astropy.io import fits

H0, OM, OBH2, NS = 67.36, 0.3153, 0.02237, 0.9649
H = H0/100; OMH2 = OM*H*H
CMB_APEX = (264.021, 48.253)  # (l, b)
SIGMA_NL = 7.0

def log(*a): print(time.strftime('%H:%M:%S'), *a, flush=True)

def comoving(z):
    from scipy.integrate import cumulative_trapezoid
    zz = np.linspace(0, 2.0, 4001); Ez = np.sqrt(OM*(1+zz)**3 + 1 - OM)
    chi = cumulative_trapezoid(2997.92458/Ez, zz, initial=0)  # Mpc/h
    return np.interp(z, zz, chi)

def load(path, zmin, zmax, nmax=None):
    d = fits.open(path)[1].data; m = (d['Z'] > zmin) & (d['Z'] < zmax)
    ra, dec, z, w = (np.asarray(d[c][m], float) for c in ('RA', 'DEC', 'Z', 'WEIGHT')); w = w*np.asarray(d['WEIGHT_FKP'][m], float)
    if nmax and len(ra) > nmax:
        i = np.random.default_rng(1).choice(len(ra), nmax, replace=False); ra, dec, z, w = ra[i], dec[i], z[i], w[i]
    return ra, dec, z, w

def xyz(dist, ra, dec):
    ra, dec = np.radians(ra), np.radians(dec); return np.array([dist*np.cos(dec)*np.cos(ra), dist*np.cos(dec)*np.sin(ra), dist*np.sin(dec)])

def regions(ra, dec, nside=8):
    import healpy as hp; return hp.ang2pix(nside, ra, dec, lonlat=True)

def count(ddir, tracer, cap, zmin, zmax, nran, nthreads, rmult):
    from pycorr import TwoPointCorrelationFunction, setup_logging
    setup_logging('warning')
    ra, dec, z, w = load(f'{ddir}/{tracer}_{cap}_clustering.dat.fits', zmin, zmax)
    rr = [load(f'{ddir}/{tracer}_{cap}_{i}_clustering.ran.fits', zmin, zmax) for i in range(nran)]
    rra, rdec, rz, rw = (np.concatenate([r[k] for r in rr]) for k in range(4))
    if len(rra) > rmult*len(ra):  # the pair counts scale with N_randoms^2; rmult x data is enough for the BAO scale
        i = np.random.default_rng(2).choice(len(rra), int(rmult*len(ra)), replace=False); rra, rdec, rz, rw = rra[i], rdec[i], rz[i], rw[i]
    log(f'{cap} {zmin}-{zmax}: {len(ra)} galaxies, {len(rra)} randoms ({nran} files, capped at {rmult} x data); z_eff (w^2-weighted) {np.sum(w*w*z)/np.sum(w*w):.3f}')
    pos = xyz(comoving(z), ra, dec); rpos = xyz(comoving(rz), rra, rdec)
    edges = (np.linspace(0, 200, 51), np.linspace(-1, 1, 201))
    lab, rlab = regions(ra, dec), regions(rra, rdec)
    # jackknife regions need randoms in every data region and vice versa: keep the common set
    common = np.intersect1d(np.unique(lab), np.unique(rlab)); km, krm = np.isin(lab, common), np.isin(rlab, common)
    t0 = time.time()
    res = TwoPointCorrelationFunction('smu', edges, data_positions1=pos[:, km], data_weights1=w[km], randoms_positions1=rpos[:, krm],
                                      randoms_weights1=rw[krm], data_samples1=lab[km], randoms_samples1=rlab[krm], position_type='xyz',
                                      engine='corrfunc', nthreads=nthreads, los='midpoint', estimator='landyszalay')
    log(f'{cap}: pair counts done in {time.time()-t0:.0f} s, {len(common)} jackknife regions')
    return res, np.sum(w*w*z)/np.sum(w*w), len(ra)

def template(zeff=0.75):
    k = np.logspace(-4, 1, 3000)
    try:
        import camb
        pars = camb.set_params(H0=H0, ombh2=OBH2, omch2=OMH2 - OBH2 - 0.00064, ns=NS, As=2.1e-9, mnu=0.06, WantTransfer=True)
        pars.set_matter_power(redshifts=[zeff], kmax=20.0, nonlinear=False)
        r = camb.get_results(pars); kh, _, pk = r.get_matter_power_spectrum(minkh=1e-4, maxkh=10, npoints=3000); plin = np.interp(k, kh, pk[0])
        src = 'camb'
    except Exception as e:
        raise SystemExit(f'camb unavailable ({e}); install it in the venv')
    # Eisenstein & Hu 1998 no-wiggle
    th = 2.7255/2.7; ob = OBH2/H/H; s = 44.5*np.log(9.83/OMH2)/np.sqrt(1 + 10*OBH2**0.75)
    ag = 1 - 0.328*np.log(431*OMH2)*ob/OM + 0.38*np.log(22.3*OMH2)*(ob/OM)**2
    geff = OM*H*(ag + (1 - ag)/(1 + (0.43*k*H*s)**4)); q = k*th*th/geff
    L0 = np.log(2*np.e + 1.8*q); C0 = 14.2 + 731/(1 + 62.5*q); T0 = L0/(L0 + C0*q*q); pnw = k**NS*T0*T0
    m = (k > 0.01) & (k < 0.5); pnw *= np.exp(np.mean(np.log(plin[m]/pnw[m])))
    pdw = (plin - pnw)*np.exp(-0.5*(k*SIGMA_NL)**2) + pnw
    kk = np.linspace(1e-4, 5, 20000); svals = np.linspace(20, 260, 481)
    def xi(p):
        pp = np.interp(kk, k, p); integ = kk*kk*pp*np.exp(-(kk*1.0)**2)/(2*np.pi*np.pi)
        return svals, np.array([np.trapezoid(integ*np.sinc(kk*sv/np.pi), kk) for sv in svals])
    return xi(pdw), xi(pnw), src

def fit(s, xi0, cov, tmpl, smin=50, smax=150, alphas=np.arange(0.8, 1.2001, 0.001)):
    """chi2(alpha) with B, a0, a1, a2 solved linearly; B is held >= 0 (a negative B fits the BAO peak upside down and makes
    spurious minima: the 0.6-0.8 NGC bin had one at alpha 1.21 with B = -4.5); where the free B is negative, B = 0."""
    m = (s > smin) & (s < smax); s, y, C = s[m], xi0[m], cov[np.ix_(m, m)]; Ci = np.linalg.inv(C)
    st, xt = tmpl; chi2 = np.empty(len(alphas)); bb = np.column_stack([np.ones_like(s), 1/s, 1/s**2])
    A0 = bb.T@Ci@bb; p0 = np.linalg.solve(A0, bb.T@Ci@y); r0 = y - bb@p0; chi2_bb = r0@Ci@r0
    for i, a in enumerate(alphas):
        M = np.column_stack([np.interp(a*s, st, xt), bb])
        A = M.T@Ci@M; b = M.T@Ci@y; p = np.linalg.solve(A, b); r = y - M@p; chi2[i] = r@Ci@r if p[0] > 0 else chi2_bb
    i = np.argmin(chi2); a = alphas[i]
    # sigma from Delta chi2 = 1 on each side
    lo = alphas[:i][chi2[:i] <= chi2[i] + 1]; hi = alphas[i:][chi2[i:] <= chi2[i] + 1]
    sig = (hi[-1] - lo[0])/2 if len(lo) and len(hi) else np.nan
    return a, sig, chi2[i], m.sum() - 4, (lo[0] if len(lo) else np.nan, hi[-1] if len(hi) else np.nan)

def cap_direction(ra, dec, w):
    from astropy.coordinates import SkyCoord; import astropy.units as u
    c = SkyCoord(ra*u.deg, dec*u.deg).galactic; l, b = np.radians(c.l.deg), np.radians(c.b.deg)
    v = np.array([np.sum(w*np.cos(b)*np.cos(l)), np.sum(w*np.cos(b)*np.sin(l)), np.sum(w*np.sin(b))]); v /= np.linalg.norm(v)
    la, ba = np.radians(CMB_APEX); apex = np.array([np.cos(ba)*np.cos(la), np.cos(ba)*np.sin(la), np.sin(ba)])
    return np.degrees(np.arctan2(v[1], v[0])) % 360, np.degrees(np.arcsin(v[2])), v@apex

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('ddir'); ap.add_argument('--tracer', default='LRG'); ap.add_argument('--nran', type=int, default=1)
    ap.add_argument('--zbins', action='store_true'); ap.add_argument('--nthreads', type=int, default=4); ap.add_argument('--fit-only', action='store_true')
    ap.add_argument('--rmult', type=float, default=4.0); ap.add_argument('--sigma-nl', type=float, default=7.0); ap.add_argument('--rebin', type=int, default=2)
    ap.add_argument('--smin', type=float, default=50.0); ap.add_argument('--smax', type=float, default=150.0)
    a = ap.parse_args(); global SIGMA_NL; SIGMA_NL = a.sigma_nl
    bins = [(0.4, 0.6), (0.6, 0.8), (0.8, 1.1)] if a.zbins else [(0.4, 1.1)]
    (st, xt), (sn, xn), src = template(); log(f'template from {src}, Sigma_nl {SIGMA_NL}')
    out = {}
    for zmin, zmax in bins:
        for cap in ('NGC', 'SGC'):
            cache = f'{a.ddir}/xi_{a.tracer}_{cap}_{zmin}_{zmax}.npy'
            from pycorr import TwoPointCorrelationFunction
            if os.path.exists(cache):
                res = TwoPointCorrelationFunction.load(cache); log(f'{cap} {zmin}-{zmax}: loaded {cache}')
            elif a.fit_only: continue
            else:
                res, zeff, n = count(a.ddir, a.tracer, cap, zmin, zmax, a.nran, a.nthreads, a.rmult); res.save(cache)
            res = res[::a.rebin]  # 8 Mpc/h bins by default, for a better-conditioned jackknife covariance
            s, xi, cov = res.get_corr(ells=(0, 2), return_sep=True, return_cov=True)
            nb = len(s); xi0, c00 = xi[0], cov[:nb, :nb]
            nreal = len(res.realizations); njk = nreal
            m = (s > a.smin) & (s < a.smax); hart = (njk - m.sum() - 2)/(njk - 1); c00 = c00/hart
            al, sig, chi2, dof, (lo, hi) = fit(s, xi0, c00, (st, xt), a.smin, a.smax); _, _, chi2n, _, _ = fit(s, xi0, c00, (sn, xn), a.smin, a.smax)
            ra, dec, z, w = load(f'{a.ddir}/{a.tracer}_{cap}_clustering.dat.fits', zmin, zmax); l, b, proj = cap_direction(ra, dec, w)
            zeff = np.sum(w*w*z)/np.sum(w*w)
            log(f'{cap} {zmin}-{zmax}: n={len(ra)}, z_eff={zeff:.3f}, jackknife regions {njk}, Hartlap {hart:.2f}; alpha_iso = {al:.4f} +{hi-al:.4f} -{al-lo:.4f}; chi2/dof {chi2:.1f}/{dof}; no-wiggle chi2 {chi2n:.1f} (BAO detection sqrt(dchi2) = {np.sqrt(max(chi2n-chi2,0)):.1f} sigma); mean direction (l,b)=({l:.0f},{b:.0f}), cos to CMB apex {proj:+.2f}')
            out[(zmin, zmax, cap)] = (al, sig, len(ra), zeff, proj, chi2, dof, chi2n, njk)
            np.savetxt(f'{a.ddir}/xi0_{a.tracer}_{cap}_{zmin}_{zmax}.txt', np.column_stack([s, xi0, np.sqrt(np.diag(c00))]), header='s xi0 err (8 Mpc/h bins, jackknife, Hartlap-scaled)')
        if all((zmin, zmax, c) in out for c in ('NGC', 'SGC')):
            aN, sN = out[(zmin, zmax, 'NGC')][:2]; aS, sS = out[(zmin, zmax, 'SGC')][:2]
            d = aN - aS; e = np.hypot(sN, sS); comb = (aN/sN**2 + aS/sS**2)/(1/sN**2 + 1/sS**2); ecomb = 1/np.sqrt(1/sN**2 + 1/sS**2)
            log(f'== {zmin}-{zmax}: alpha_NGC - alpha_SGC = {d:+.4f} ± {e:.4f} ({d/e:+.2f} sigma); inverse-variance combined alpha_iso = {comb:.4f} ± {ecomb:.4f}; the caps\' CMB-apex projections differ by {out[(zmin, zmax, "NGC")][4] - out[(zmin, zmax, "SGC")][4]:+.2f} in cos, so a dipole in alpha of amplitude A would give d = {out[(zmin, zmax, "NGC")][4] - out[(zmin, zmax, "SGC")][4]:+.2f} A: A = {d/(out[(zmin, zmax, "NGC")][4] - out[(zmin, zmax, "SGC")][4]):+.4f} ± {e/abs(out[(zmin, zmax, "NGC")][4] - out[(zmin, zmax, "SGC")][4]):.4f}')

main()
