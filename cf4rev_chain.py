"""R-168 (cosmo-cf4-rev): the CF4-like and 2M++-like mock chain that calibrates the null of the SN-vs-field beta test.

Route 2 of the brief (no public CF4 posterior realisations: the IP2I Cosmicflows page serves only the CF4++ mean and std
grids and the 2023 grids). Baseline realisation (variants below change one thing each):
 1. truth: a LCDM Gaussian density field on a 256^3 grid (L = 1000 Mpc/h, 3.9 Mpc/h cells: power below the reconstruction's
    7.8 Mpc/h cells), supergalactic Cartesian, observer at the centre; Eisenstein-Hu no-wiggle P(k) (Om 0.31, Ob 0.049, h 0.68,
    ns 0.965, sigma8 0.81), linear velocity, f = Om^0.55; plus 250 km/s unresolved velocity per object.
 2. positions in redshift space: every object (CF4 galaxy, SN) sits at its real observed position r_obs; its true position
    solves r_true + (v_r(r_true) + v_nl)/H = r_obs; the truth is read at r_true, the reconstructions at r_obs (as the data).
 3. CF4-like catalogue: the 55,694 real CF4 galaxies with a distance error (VizieR J/ApJ/944/94 table 2), each observed as a
    radial velocity u = v_true,r + v_nl + noise, noise sd = (ln10/5) e_DM Vcmb (the catalogue's own error). Shared inputs: the
    CF4 galaxy matched to each Pantheon+ SN host carries that SN's own mock light-curve noise and velocity.
 4. CF4-like field: the Wiener filter (posterior mean under the prior and the catalogue's stated errors) on the CF4++ grid
    (128^3), preconditioned conjugate gradients in Fourier space; its linear velocity.
 5. 2M++-like field (on the 256^3 grid): galaxies Poisson-sampled in real space from n(r) exp(b delta - b^2 var/2) (log-normal,
    large-scale bias b), b = 1.2, n(r) from a
    2MASS K-band Schechter function (phi* 0.0116 h^3/Mpc^3, M* - 5 log h = -23.39, alpha = -1.09), K < 11.5 (2MRS) or 12.5
    (6dF south, SDSS north cap), within 200 Mpc/h; moved to redshift space; the zone of avoidance (|b| < 5, < 10 within 30 deg of
    the Galactic centre) filled by cloning the mirrored strip above/below it; delta_g smoothed by a 4 Mpc/h Gaussian; three
    iterations of the real-space reconstruction (positions moved back by the predicted velocity), as Carrick+2015; the scale
    and an external dipole V_ext fitted to the mock CF4 velocities within 100 Mpc/h and the Pantheon+ SNe at z < 0.067 (sharing
    their mock noise: Carrick et al. fitted SFI++ and the A1 SNe, which are Pantheon+ objects); field = scale v_pred + V_ext inside
    the Galactic +-200 Mpc/h cube, V_ext outside.
 6. outputs at each sample's SN positions: the truth, the CF4-like and the 2M++-like radial velocities; the Pantheon+ mock
    residual y = corr(v_true + v_nl) + eps_lc, eps_lc ~ N(0, C - D_v) (C without the VPEC systematic, D_v its velocity block).
Variants: s8lo / s8hi (WF prior sigma8 0.70 / 0.90), nslo (WF prior ns 0.90 at fixed sigma8), smooth (WF field smoothed by an
8 Mpc/h Gaussian), grouped (CF4 groups, 1PGC: shared v_nl, one datum per group), lognorm (log-normal distance errors, no Malmquist
correction), err13 (true errors 1.3x the stated ones the filter uses), zpstep (FP-only distances +0.05 mag north of the
equator, -0.05 south: R-119's measured step), noshare (no shared SN inputs), old (128^3 truth, real-space positions: the
first chain), b10 / b14 (2M++-like bias 1.0 / 1.4), vnl150 / vnl400 (unresolved velocity), tr074 / tr086 (the truth's sigma8,
prior matched), clustered (objects moved into the mock overdensities within +-7.8 Mpc/h), wfest (log-normal distance errors
through the real-data Watkins-Feldman estimator and H0 rescale: the null of WF(real CF4)), lpt2 (a 2LPT particle truth: objects
are particles near their real positions, velocities the particles' 2LPT velocities; the 2M++-like galaxies are the particles
thinned to n(r), bias 1). 'cf4' variants skip the 2M++-like part and vice versa.
usage: python cf4rev_chain.py <variant> <first seed> <n realisations> <outdir>
"""
import sys, os, time, numpy as np, pandas as pd
import scipy.fft as sfft
from scipy.integrate import quad
from scipy.ndimage import map_coordinates
from astropy.coordinates import SkyCoord
import astropy.units as u
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cf4rev_common import D, C, dc_h, corr_mag, sg_xyz, gal_xyz, pantheon, des, ztf

L = 1000.0; H = 100.0
OMm, OMb, h, NS, S8 = 0.31, 0.049, 0.68, 0.965, 0.81; F = OMm**0.55
VARIANT = sys.argv[1] if len(sys.argv) > 1 else 'base'
SIG_NL = {'vnl150': 150.0, 'vnl400': 400.0}.get(VARIANT, 250.0)
SIG_GRP = 150.0/250.0*SIG_NL   # the part of v_nl shared by a CF4 group (1PGC); the rest is per object
WORKERS = int(os.environ.get('FFT_WORKERS', '1'))
CF4_ONLY = {'s8lo', 's8hi', 'nslo', 'smooth', 'grouped', 'lognorm', 'err13', 'zpstep', 'noshare', 'vnl150', 'vnl400', 'tr074', 'tr086', 'clustered', 'wfest'}
TWOMPP_ONLY = {'b10', 'b14'}

def eh_nowiggle(k):  # k in h/Mpc
    om = OMm*h*h; fb = OMb/OMm; s = 44.5*np.log(9.83/om)/np.sqrt(1 + 10*(OMb*h*h)**0.75)
    a = 1 - 0.328*np.log(431*om)*fb + 0.38*np.log(22.3*om)*fb**2
    gam = OMm*h*(a + (1 - a)/(1 + (0.43*k*h*s)**4)); q = k/gam
    l0 = np.log(2*np.e + 1.8*q); c0 = 14.2 + 731/(1 + 62.5*q); return l0/(l0 + c0*q*q)
def tophat(x): return 3*(np.sin(x) - x*np.cos(x))/x**3
def make_pk(ns=NS, s8=S8):
    pu = lambda k: k**ns*eh_nowiggle(k)**2
    amp = s8**2/quad(lambda lk: np.exp(lk)**3*pu(np.exp(lk))*tophat(8*np.exp(lk))**2/(2*np.pi**2), np.log(1e-4), np.log(50), limit=400)[0]
    return lambda k: amp*pu(k)

class Grid:
    def __init__(self, N, pk):
        self.N = N; self.dx = L/N; self.vc = self.dx**3
        kx = 2*np.pi*sfft.fftfreq(N, self.dx); kz = 2*np.pi*sfft.rfftfreq(N, self.dx)
        KX, KY, KZ = np.meshgrid(kx, kx, kz, indexing='ij'); K2 = KX**2 + KY**2 + KZ**2; K2[0, 0, 0] = 1; self.K2 = K2
        self.sqp = np.sqrt(pk(np.sqrt(K2))/self.vc); self.sqp[0, 0, 0] = 0
        nyq = np.pi/self.dx
        self.G = [1j*F*H*Kj/K2*(np.abs(Kj) < nyq*0.999) for Kj in (KX, KY, KZ)]
        self.ax = (np.arange(N) + 0.5)*self.dx - L/2
    def rf(self, a): return sfft.rfftn(a, workers=WORKERS)
    def irf(self, a): return sfft.irfftn(a, s=(self.N,)*3, workers=WORKERS)
    def vel(self, dk): return np.stack([self.irf(dk*g) for g in self.G])
    def smooth(self, a, R): return self.irf(self.rf(a)*np.exp(-0.5*self.K2*R*R))
    def at(self, V, xyz):
        """3-vector field V (3,N,N,N) at positions xyz (3,n), trilinear."""
        ci = (xyz + L/2)/self.dx - 0.5
        return np.stack([map_coordinates(V[k], ci, order=1, mode='nearest') for k in range(3)])
    def at_r(self, V, xyz): return np.sum(self.at(V, xyz)*xyz/np.linalg.norm(xyz, axis=0), 0)
    def cell(self, xyz): return np.ravel_multi_index(tuple(np.clip(((xyz + L/2)//self.dx).astype(int), 0, self.N - 1)), (self.N,)*3)

class WF:
    """Wiener filter of delta from radial velocities binned NGP on the grid: per data cell the 3x3 sum of w rhat rhat^T and
    the 3-vector sum of w u rhat. Solves (I + A^T M A) x = A^T b, A = G S^1/2; delta_WF = S^1/2 x."""
    def __init__(self, g, xyz, sig):
        self.g = g; self.cell = g.cell(xyz)
        self.rhat = xyz/np.linalg.norm(xyz, axis=0); self.w = 1/sig**2
        self.uc, self.inv = np.unique(self.cell, return_inverse=True)
        self.M = np.zeros((len(self.uc), 3, 3))
        for a in range(3):
            for b in range(3): self.M[:, a, b] = np.bincount(self.inv, self.w*self.rhat[a]*self.rhat[b], len(self.uc))
    def A(self, x): return self.g.vel(self.g.rf(x)*self.g.sqp)
    def AT(self, v): return self.g.irf(self.g.sqp*sum(np.conj(gg)*self.g.rf(v[j]) for j, gg in enumerate(self.g.G)))
    def Mmul(self, v):
        out = np.zeros_like(v); vv = v.reshape(3, -1)[:, self.uc]; o = out.reshape(3, -1)
        o[:, self.uc] = np.einsum('cab,bc->ac', self.M, vv); return out
    def solve(self, u, tol=1e-5, maxit=500):
        N = self.g.N; bv = np.zeros((3, N**3))
        for a in range(3): bv[a, self.uc] = np.bincount(self.inv, self.w*u*self.rhat[a], len(self.uc))
        rhs = self.AT(bv.reshape(3, N, N, N)); Hop = lambda x: x + self.AT(self.Mmul(self.A(x)))
        x = np.zeros_like(rhs); r = rhs.copy(); p = r.copy(); rr = np.sum(r*r); r0 = np.sqrt(rr)
        for it in range(maxit):
            Hp = Hop(p); al = rr/np.sum(p*Hp); x += al*p; r -= al*Hp; rn = np.sum(r*r)
            if np.sqrt(rn) < tol*r0: break
            p = r + rn/rr*p; rr = rn
        return self.A(x), it + 1

def rot(frame_to):
    """rotation matrix from supergalactic Cartesian to another frame's Cartesian."""
    sc = SkyCoord(sgx=[1, 0, 0], sgy=[0, 1, 0], sgz=[0, 0, 1], unit='Mpc', frame='supergalactic', representation_type='cartesian')
    c = getattr(sc, frame_to).cartesian; return np.stack([c.x.value, c.y.value, c.z.value])  # columns = images of SG axes
R_GAL = rot('galactic'); R_EQ = rot('icrs')

def nbar_fn():
    phi, Ms, al = 0.0116, -23.39, -1.09
    def n_of(K, r): xl = 10**(-0.4*(K - 5*np.log10(r) - 25 - Ms)); return phi*quad(lambda x: x**al*np.exp(-x), xl, 60)[0]
    rr = np.linspace(0.5, 300, 600); n11 = np.array([n_of(11.5, r) for r in rr]); n12 = np.array([n_of(12.5, r) for r in rr])
    return lambda r, deep: np.where(deep, np.interp(r, rr, n12), np.interp(r, rr, n11))

def sky(xyz_sg):
    """galactic l, b and declination of SG positions."""
    gx = R_GAL @ xyz_sg; ex = R_EQ @ xyz_sg; r = np.linalg.norm(xyz_sg, axis=0) + 1e-9
    return np.degrees(np.arctan2(gx[1], gx[0])) % 360, np.degrees(np.arcsin(np.clip(gx[2]/r, -1, 1))), np.degrees(np.arcsin(np.clip(ex[2]/r, -1, 1)))
def masked_deep(xyz_sg):
    l, b, dec = sky(xyz_sg)
    bm = np.where((l < 30) | (l > 330), 10.0, 5.0)
    return np.abs(b) < bm, ((dec < 0) & (np.abs(b) > 10)) | (b > 30), bm, l, b

def load_cf4():
    raw = [l for l in open(os.path.join(D, 'cf4_table2m.tsv')).read().split('\n') if l.strip() and not l.startswith('#')]
    cols = [c.strip() for c in raw[0].split('\t')]
    t = pd.DataFrame([l.split('\t') for l in raw[3:]], columns=cols).apply(pd.to_numeric, errors='coerce')
    return t.dropna(subset=['Vcmb', 'DM', 'e_DM']).query('Vcmb > 300 and e_DM > 0').reset_index(drop=True)

def ztf_hosts(t, Z):
    """ZTF SN -> CF4 galaxy (host coordinates from the DR2 host table), within 30" and 1000 km/s: 113 of 1,122."""
    zt = pd.read_csv(os.path.join(D, 'ztf', 'ztfsniadr2_lite', 'tables', 'snia_data.csv'))
    gh = pd.read_csv(os.path.join(D, 'ztf', 'ztfsniadr2_lite', 'tables', 'globalhost_data.csv')).set_index('ztfname')
    key = pd.Series(zt.ztfname.values, index=list(zip(np.round(zt.ra, 6), np.round(zt.dec, 6))))
    names = key.reindex(list(zip(np.round(Z['ra'], 6), np.round(Z['dec'], 6)))).values
    hra = gh.ra_host.reindex(names).values; hde = gh.dec_host.reindex(names).values
    hra = np.where(np.isnan(hra), Z['ra'], hra); hde = np.where(np.isnan(hde), Z['dec'], hde)
    j, sep, _ = SkyCoord(hra*u.deg, hde*u.deg).match_to_catalog_sky(SkyCoord(t.RAJ2000.values*u.deg, t.DEJ2000.values*u.deg))
    return np.where((sep.arcsec < 30) & (np.abs(t.Vcmb.values[j] - C*Z['zc']) < 1000), j, -1)

def setup(variant):
    t = load_cf4(); S = dict(variant=variant, t=t)
    S['gx'] = sg_xyz(t.RAJ2000.values, t.DEJ2000.values, t.Vcmb.values/H)          # observed positions
    S['gsd'] = np.log(10)/5*t.e_DM.values*t.Vcmb.values                            # stated distance-error in km/s
    S['fp_only'] = t.DMfp.notna().values & t.DMtf.isna().values & t.DMsnIa.isna().values
    S['north'] = t.DEJ2000.values > 0
    samples = [pantheon(), des(), ztf()]; P = samples[0]; S['samples'] = samples
    hra = np.where(P['host_ra'] > -900, P['host_ra'], P['ra']); hde = np.where(P['host_dec'] > -900, P['host_dec'], P['dec'])
    j, sep, _ = SkyCoord(hra*u.deg, hde*u.deg).match_to_catalog_sky(SkyCoord(t.RAJ2000.values*u.deg, t.DEJ2000.values*u.deg))
    P['host'] = np.where((sep.arcsec < 120) & (np.abs(t.Vcmb.values[j] - C*P['zc']) < 1000), j, -1)
    cid = np.asarray(P['cid']).astype(str); ev = P['err_vpec']
    w, U = np.linalg.eigh(P['cov'] - np.outer(ev, ev)*(cid[:, None] == cid[None, :]))
    P['Llc'] = U*np.sqrt(np.clip(w, 0, None)); P['lc_negeig'] = w.min()
    _, P['cid_inv'] = np.unique(cid, return_inverse=True)
    S['gid'] = t['1PGC'].fillna(t.PGC).values
    Z = samples[2]; Z['host'] = ztf_hosts(t, Z)
    for smp in (P, Z): smp['grp'] = np.where(smp['host'] >= 0, S['gid'][np.clip(smp['host'], 0, None)], -1 - np.arange(len(smp['zc'])))
    first = np.zeros(len(cid), bool); _, fi = np.unique(P['cid_inv'], return_index=True); first[fi] = True; P['first'] = first
    for s in samples:
        s['xyz'] = sg_xyz(s['ra'], s['dec'], dc_h(s['zc'])); s['gxyz'] = gal_xyz(s['ra'], s['dec'], dc_h(s['zc']))
    old = variant == 'old'
    pk_true = make_pk(s8={'tr074': 0.74, 'tr086': 0.86}.get(variant, S8))
    S['gT'] = Grid(128 if old else 256, pk_true)
    pk_wf = {'s8lo': make_pk(s8=0.70), 's8hi': make_pk(s8=0.90), 'nslo': make_pk(ns=0.90)}.get(variant, pk_true)  # tr0xx: prior = truth
    S['gW'] = S['gT'] if old else Grid(128, pk_wf)
    if variant == 'grouped':
        gid = t['1PGC'].fillna(t.PGC).values; ug, ginv = np.unique(gid, return_inverse=True)
        wd = 1/S['gsd']**2; sw = np.bincount(ginv, wd)
        S['ginv'] = ginv; S['ng'] = len(ug)
        S['grp_xyz'] = np.stack([np.bincount(ginv, wd*S['gx'][k])/sw for k in range(3)])
        S['grp_sig'] = np.sqrt(1/sw + SIG_NL**2)
        S['wf'] = WF(S['gW'], S['grp_xyz'], S['grp_sig'])
    elif variant not in TWOMPP_ONLY:
        S['wf'] = WF(S['gW'], S['gx'], np.sqrt(S['gsd']**2 + SIG_NL**2))
    S['near'] = t.Vcmb.values < 10000
    if variant not in CF4_ONLY and not old:
        g = S['gT']; S['g2'] = g
        ax = g.ax; X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij'); xyz = np.stack([X.ravel(), Y.ravel(), Z.ravel()]); del X, Y, Z
        r = np.linalg.norm(xyz, axis=0)
        msk, deep, bm, l, b = masked_deep(xyz); del xyz
        S['nbar'] = nbar_fn()
        S['cell_r'] = r; S['cell_deep'] = deep; S['cell_mask'] = msk; S['cell_in'] = r < 200
        mi = np.where(msk & (r < 200))[0]
        bb = np.radians(np.sign(b[mi] + 1e-9)*(2*bm[mi] - np.abs(b[mi]))); ll = np.radians(l[mi])
        gpos = r[mi]*np.stack([np.cos(bb)*np.cos(ll), np.cos(bb)*np.sin(ll), np.sin(bb)])
        S['clone_from'] = g.cell(R_GAL.T @ gpos); S['clone_to'] = mi
        S['incube_s'] = {s['name']: np.all(np.abs(s['gxyz']) < 200, axis=0) for s in samples}
    return S

def true_pos(g, Vt, xyz_obs, vnl, it=3):
    """the truth's radial velocity at r_true, where r_true + (v_r(r_true) + v_nl)/H = r_obs along the observed direction."""
    rh = xyz_obs/np.linalg.norm(xyz_obs, axis=0); r_obs = np.linalg.norm(xyz_obs, axis=0); r = r_obs.copy()
    for _ in range(it):
        v = g.at_r(Vt, rh*r); r = np.clip(r_obs - (v + vnl)/H, 0.5, None)
    return g.at_r(Vt, rh*r)

OFFS = np.stack(np.meshgrid(*[np.arange(-2, 3)]*3, indexing='ij')).reshape(3, -1)
def twompp_like(S, rng, delta, Vt, bias, fit_u, fit_xyz, fit_w, gal=None):
    g = S['g2']; N = g.N; nb = S['nbar']
    if gal is not None:   # lpt2: galaxies are given (real-space positions, true radial velocities, the sampling density)
        pos, vr, nb = gal; rh = pos/np.linalg.norm(pos, axis=0)
        s_r = np.linalg.norm(pos, axis=0) + (vr + rng.normal(0, SIG_NL, pos.shape[1]))/H
        return _recon(S, g, N, nb, pos, rh, s_r, bias, fit_u, fit_xyz, fit_w)
    # log-normal galaxy field: exp(b delta - b^2 var/2) keeps the large-scale bias b without clipping (sigma_delta ~ 1.5 here)
    lam = nb(S['cell_r'], S['cell_deep'])*g.vc*np.exp(bias*delta.ravel() - 0.5*bias**2*delta.var())*(S['cell_in'] & ~S['cell_mask'])
    cnt = rng.poisson(lam); idx = np.repeat(np.arange(N**3), cnt)
    pos = np.stack(np.unravel_index(idx, (N,)*3)).astype(float)*g.dx - L/2 + rng.random((3, len(idx)))*g.dx   # real space
    rh = pos/np.linalg.norm(pos, axis=0)
    s_r = np.linalg.norm(pos, axis=0) + (g.at_r(Vt, pos) + rng.normal(0, SIG_NL, len(idx)))/H                    # redshift space
    return _recon(S, g, N, nb, pos, rh, s_r, bias, fit_u, fit_xyz, fit_w)

def _recon(S, g, N, nb, pos, rh, s_r, bias, fit_u, fit_xyz, fit_w):
    deep = masked_deep(pos)[1]; r = s_r.copy(); cin = S['cell_in'].reshape((N,)*3)
    for it in range(4):
        p = rh*r; keep = r < 200
        dg = np.bincount(g.cell(p[:, keep]), 1/(nb(r[keep], deep[keep])*g.vc), N**3) - 1.0
        dg[~S['cell_in']] = 0.0
        dg[S['clone_to']] = dg[S['clone_from']]
        dg = g.smooth(dg.reshape((N,)*3), 4.0); dg[~cin] = 0.0
        Vp = g.vel(g.rf(dg)/bias)
        if it < 3: r = np.clip(s_r - g.at_r(Vp, p)/H, 0.5, None)
    rhat = fit_xyz/np.linalg.norm(fit_xyz, axis=0)
    A = np.column_stack([g.at_r(Vp, fit_xyz), rhat.T])*fit_w[:, None]
    coef, *_ = np.linalg.lstsq(A, fit_u*fit_w, rcond=None)
    return Vp, coef, pos.shape[1]

def lpt2(g, dk):
    """2LPT particles on the grid's Lagrangian lattice: x = q + psi1 + psi2, v = f1 H psi1 + f2 H psi2 (f2 = 2 Om^(6/11)),
    psi1(k) = i k delta/k^2, psi2(k) = (3/7) i k S2/k^2, S2 = sum_{i<j} (phi_ii phi_jj - phi_ij^2), phi_ij(k) = k_i k_j delta/k^2."""
    N = g.N; kv = [g.G[j]/(1j*F*H) for j in range(3)]           # = k_j/k^2 (zero on the Nyquist planes)
    kk = [g.G[j]/(1j*F*H)*g.K2 for j in range(3)]               # = k_j
    ph = lambda i, j: g.irf(kk[i]*kk[j]/g.K2*dk)
    pxx, pyy, pzz = ph(0, 0), ph(1, 1), ph(2, 2)
    S2 = pxx*pyy + pxx*pzz + pyy*pzz; del pxx, pyy, pzz
    for i, j in ((0, 1), (0, 2), (1, 2)): S2 -= ph(i, j)**2
    s2k = g.rf(S2); del S2
    f2 = 2*OMm**(6/11)
    q = np.stack(np.unravel_index(np.arange(N**3), (N,)*3)).astype(np.float32)*np.float32(g.dx) - np.float32(L/2 - g.dx/2)
    x = np.empty_like(q); v = np.empty_like(q)
    for j in range(3):
        p1 = g.irf(1j*kv[j]*dk).ravel(); p2 = g.irf(1j*kv[j]*s2k*3/7).ravel()
        x[j] = (q[j] + p1 + p2 + L/2) % L - L/2; v[j] = F*H*p1 + f2*H*p2
    del q
    cell = g.cell(x); order = np.argsort(cell, kind='stable'); cnt = np.bincount(cell, minlength=N**3)
    start = np.concatenate([[0], np.cumsum(cnt)[:-1]])
    return x, v, order, start, cnt

def pick_particles(g, P3, xyz, rng):
    """for each object, a random particle from the 5x5x5 cells around its position, cells weighted by their particle count."""
    x, v, order, start, cnt = P3
    c0 = np.clip(((xyz + L/2)//g.dx).astype(int), 2, g.N - 3)
    cand = np.ravel_multi_index(tuple(c0[:, :, None] + OFFS[:, None, :]), (g.N,)*3)
    w = cnt[cand].astype(float); w /= w.sum(1, keepdims=True)
    k = (w.cumsum(1) < rng.random(len(w))[:, None]).sum(1); k = np.minimum(k, w.shape[1] - 1)
    ch = cand[np.arange(len(w)), k]
    return order[start[ch] + (rng.random(len(ch))*cnt[ch]).astype(int)]

OFFS = np.stack(np.meshgrid(*[np.arange(-2, 3)]*3, indexing='ij')).reshape(3, -1)
from scipy.interpolate import interp1d
_zg = np.linspace(0, 0.2, 4001); _zinv = interp1d(dc_h(_zg), _zg)
def cluster_resample(g, delta, xyz, rng, bias=1.2):
    """'clustered' variant: move each object to a cell within +-2 cells (+-7.8 Mpc/h) of its real position, chosen with
    weight exp(bias delta): tracers and hosts then sit in the mock overdensities, the survey selection is kept on larger scales."""
    c0 = np.clip(((xyz + L/2)//g.dx).astype(int), 2, g.N - 3)
    cand = c0[:, :, None] + OFFS[:, None, :]
    w = np.exp(bias*delta[cand[0], cand[1], cand[2]]); w /= w.sum(1, keepdims=True)
    k = (w.cumsum(1) < rng.random(len(w))[:, None]).sum(1)
    ch = cand[:, np.arange(len(w)), k]
    return ch*g.dx - L/2 + rng.random(ch.shape)*g.dx

def realise(S, seed):
    rng = np.random.default_rng(seed); t0 = time.time(); var = S['variant']; gT = S['gT']
    dk = gT.rf(rng.standard_normal((gT.N,)*3))*gT.sqp
    delta = gT.irf(dk); Vt = gT.vel(dk)
    P3 = lpt2(gT, dk) if var == 'lpt2' else None
    del dk
    def place(xyz):   # clustered / lpt2: the object's true position and true radial velocity
        if var == 'lpt2':
            i = pick_particles(gT, P3, xyz, rng); xt = P3[0][:, i].astype(float)
            return xt, np.sum(P3[1][:, i]*xt, 0)/np.linalg.norm(xt, axis=0)
        xt = cluster_resample(gT, delta, xyz, rng); return xt, gT.at_r(Vt, xt)
    P = S['samples'][0]; out = {}; real_space = var == 'old'
    # small-scale velocity: a group part (SIG_GRP) shared by every object in a CF4 group (1PGC), plus a per-object part
    t = S['t']; ng = len(t); ug, ginv = np.unique(S['gid'], return_inverse=True)
    vgrp = rng.normal(0, SIG_GRP, len(ug))[ginv]
    if var == 'grouped': vind_g = np.zeros(ng); vgrp = rng.normal(0, SIG_NL, len(ug))[ginv]
    else: vind_g = rng.normal(0, np.sqrt(SIG_NL**2 - SIG_GRP**2), ng)
    vnl_g = vgrp + vind_g
    Z = S['samples'][2]
    vnl_s = {}; vt = {}; xobs = {s['name']: s['xyz'] for s in S['samples']}; zcm = {s['name']: s['zc'] for s in S['samples']}
    for s in S['samples']:
        nm = s['name']
        if nm == 'Pantheon+':
            vi = (rng.standard_normal(P['cid_inv'].max() + 1))[P['cid_inv']]
            hk = P['host'] >= 0
            vnl = np.where(hk, vgrp[np.clip(P['host'], 0, None)] + vi*np.sqrt(SIG_NL**2 - SIG_GRP**2), vi*SIG_NL)
        elif nm == 'ZTF':
            vi = rng.standard_normal(len(s['zc'])); hk = s['host'] >= 0
            vnl = np.where(hk, vgrp[np.clip(s['host'], 0, None)] + vi*np.sqrt(SIG_NL**2 - SIG_GRP**2), vi*SIG_NL)
        else: vnl = np.zeros(len(s['zc']))   # DES: velocity noise is inside its own covariance, drawn at fit time
        vnl_s[nm] = vnl
        if var in ('clustered', 'lpt2'):
            xt, v = place(s['xyz'])
            xo = xt + xt/np.linalg.norm(xt, axis=0)*(v + vnl)/H; vt[nm] = v; xobs[nm] = xo
            zcm[nm] = _zinv(np.linalg.norm(xo, axis=0))
            sc = SkyCoord(sgx=xo[0], sgy=xo[1], sgz=xo[2], unit='Mpc', frame='supergalactic', representation_type='cartesian').icrs
            out[nm + '_zc'] = zcm[nm]; out[nm + '_ra'] = sc.ra.deg; out[nm + '_dec'] = sc.dec.deg
        else:
            vt[nm] = gT.at_r(Vt, s['xyz']) if real_space else true_pos(gT, Vt, s['xyz'], vnl)
    out['ZTF_vnl'] = vnl_s['ZTF']
    eps = P['Llc'] @ rng.standard_normal(P['Llc'].shape[1])
    yP = corr_mag(zcm['Pantheon+'], vt['Pantheon+'] + vnl_s['Pantheon+']) + eps
    # the matched CF4 galaxy of a ZTF host carries the SN's own small-scale velocity (its distance is independent: TF/FP)
    zh = Z['host'] >= 0; vnl_g[Z['host'][zh]] = vnl_s['ZTF'][zh]
    gx = S['gx']; gsd = S['gsd']
    if var in ('clustered', 'lpt2'):
        gxt, vtr = place(S['gx'])
        gx = gxt + gxt/np.linalg.norm(gxt, axis=0)*(vtr + vnl_g)/H
        gsd = np.log(10)/5*t.e_DM.values*H*np.linalg.norm(gx, axis=0)
        wf = WF(S['gW'], gx, np.sqrt(gsd**2 + SIG_NL**2))
    else:
        vtr = gT.at_r(Vt, S['gx']) if real_space else true_pos(gT, Vt, S['gx'], vnl_g); wf = S.get('wf')
    vtrue_g = vtr + vnl_g
    if var in ('lognorm', 'wfest'):   # log-normal distance error, cz fixed
        cz = t.Vcmb.values; d_obs = np.clip((cz - vtrue_g)/H, 0.5, None)*10**(0.2*t.e_DM.values*rng.standard_normal(ng))
        if var == 'lognorm': u = cz - H*d_obs                               # no Malmquist correction
        else:                                                                # the real-data estimator and H0 rescale
            u = cz*np.log(cz/(H*d_obs)); u = u - np.sum(u*cz)/np.sum(cz*cz)*cz
    else:
        u = vtrue_g + rng.standard_normal(ng)*gsd*(1.3 if var == 'err13' else 1.0)
    if var == 'zpstep':    # FP-only distances: +0.05 mag north, -0.05 south -> u shifts by -(ln10/5) V dDM
        u = u - np.log(10)/5*t.Vcmb.values*np.where(S['north'], 0.05, -0.05)*S['fp_only']
    hh = (P['host'] >= 0) & P['first']
    if var != 'noshare':
        u[P['host'][hh]] = vtr[P['host'][hh]] + vnl_s['Pantheon+'][hh] - np.log(10)/5*C*zcm['Pantheon+'][hh]*eps[hh]
    if var not in TWOMPP_ONLY:
        if var == 'grouped':
            wd = 1/S['gsd']**2; Vwf, nit = S['wf'].solve(np.bincount(S['ginv'], wd*u)/np.bincount(S['ginv'], wd))
        else: Vwf, nit = wf.solve(u)
        if var == 'smooth': Vwf = np.stack([S['gW'].smooth(Vwf[k], 8.0) for k in range(3)])
        for s in S['samples']: out[s['name'] + '_cf4'] = S['gW'].at_r(Vwf, xobs[s['name']])
        out['nit'] = nit; del Vwf
    if var not in CF4_ONLY and var != 'old':
        bias = {'b10': 1.0, 'b14': 1.4}.get(var, 1.2)
        lowz = P['zc'] < 0.067
        fit_u = np.r_[u[S['near']], -np.log(10)/5*C*P['zc'][lowz]*yP[lowz]]
        fit_xyz = np.c_[S['gx'][:, S['near']], P['xyz'][:, lowz]]
        fit_w = np.r_[1/np.sqrt(S['gsd'][S['near']]**2 + SIG_NL**2), 1/(np.log(10)/5*C*P['zc'][lowz]*np.sqrt(np.diag(P['cov'])[lowz]))]
        gal = None
        if var == 'lpt2':   # galaxies: the 2LPT particles thinned to n(r) (bias 1), real-space positions and velocities
            bias = 1.0; x, v = P3[0], P3[1]; r = np.linalg.norm(x, axis=0); sel = np.where(r < 200)[0]
            xs = x[:, sel].astype(float); msk, deep = masked_deep(xs)[:2]; n_p = 1/gT.vc
            nbe = lambda rr, dd: np.minimum(S['nbar'](rr, dd), n_p)
            keep = ~msk & (rng.random(len(sel)) < nbe(r[sel], deep)/n_p)
            xs = xs[:, keep]; vr = np.sum(v[:, sel[keep]]*xs, 0)/np.linalg.norm(xs, axis=0); gal = (xs, vr, nbe)
        Vp, coef, ngal = twompp_like(S, rng, delta, Vt, bias, fit_u, fit_xyz, fit_w, gal)
        for s in S['samples']:
            rh = s['xyz']/np.linalg.norm(s['xyz'], axis=0); vext = coef[1:]@rh
            out[s['name'] + '_2mpp'] = np.where(S['incube_s'][s['name']], coef[0]*S['g2'].at_r(Vp, s['xyz']) + vext, vext)
        out['scale2'] = coef[0]; out['vext2'] = coef[1:]; out['ngal2'] = ngal
    for s in S['samples']: out[s['name'] + '_true'] = vt[s['name']]
    out['Pantheon+_y'] = yP; out['secs'] = time.time() - t0
    return out

def main(variant, seed0, n, outdir):
    os.makedirs(outdir, exist_ok=True); t0 = time.time()
    S = setup(variant); P = S['samples'][0]
    print(f"[{variant}] setup {time.time()-t0:.0f}s: {len(S['t'])} CF4 galaxies; Pantheon+ {len(P['zc'])} SNe, {(P['host'] >= 0).sum()} matched "
          f"to a CF4 galaxy; C - D_v min eigenvalue {P['lc_negeig']:.2e}", flush=True)
    for sd in range(seed0, seed0 + n):
        fn = os.path.join(outdir, f'real_{sd:04d}.npz')
        if os.path.exists(fn): continue
        o = realise(S, sd); np.savez(fn, **o)
        print(f"[{variant}] seed {sd}: {o['secs']:.0f}s" + (f", CG {o['nit']}" if 'nit' in o else '') +
              (f", 2M++-like {o['ngal2']} gal, scale {o['scale2']:.2f}, Vext {np.round(o['vext2']).astype(int)}" if 'scale2' in o else ''), flush=True)

if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
