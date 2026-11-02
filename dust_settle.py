"""Settle the dust candidate of finding.md (SN Hubble residuals track Planck-GNILC minus SFD E(B-V)), R-116.
Four parts, each a labelled block of output:
 (a) emission-free reddening maps that reach past Edenhofer's 1.25 kpc: Bayestar 2019 (Green+2019, Gaia/PS1/2MASS
     stellar 3D map, the far edge of its best-fit profile, DM 18.9 = 60 kpc; Dec > -30) and Lenz+2017 (HI4PI column
     density scaled to E(B-V), defined where N_HI < 4e20, ~40% of the sky), plus Edenhofer 2023 if present. Each map is
     put on SFD's scale by one factor k (least squares MAP ~ k SFD over the sample), then the residual is regressed on
     (k MAP - SFD) alone and jointly with (GNILC - SFD), for Pantheon+ (full STAT+SYS covariance) and ZTF DR2 (diagonal
     errors, alpha and beta free). If GNILC is right and SFD wrong, the emission-free difference carries the GNILC slope
     (~ +1.8) and takes over from GNILC - SFD in the joint fit.
 (c) the physical check: propagate dE = GNILC - MWEBV (the E(B-V) actually applied, 0.86 SFD for Pantheon+) through the
     SALT2 fit analytically: dmB = R(lamB(1+z)) dE and dc = [R(lamB(1+z)) - R(lamV(1+z))] dE with a Fitzpatrick99 R_V=3.1
     law at the observer-frame wavelengths of rest B and V, so dmu = dmB - beta dc, and refit the slopes on the
     GNILC-corrected residuals. If the slope were an extinction error of GNILC's size it would go to zero.
 (d) the Galactic-axis anisotropy terms of frames.py (sin b dipole, P2(sin b) quadrupole, full dipole, CMB-axis
     quadrupole) on the Pantheon+ residuals before and after the GNILC correction.
 (b) trials accounting is in results.md (counted from the scripts, not computed here); the corrected p-values printed at
     the end use the counts given on the command line (--trials N).
usage: python dust_settle.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/ztf/ztfsniadr2_lite/tables/snia_data.csv data/dustmaps
"""
import sys, numpy as np, pandas as pd, healpy as hp, h5py
from scipy.integrate import quad
from scipy.stats import chi2 as chi2d, norm
from astropy.coordinates import SkyCoord
import astropy.units as u
import extinction
C = 299792.458; OM = 0.334; AL, BE = 0.15, 3.1
R_gal = np.array([[-0.0548755604, -0.8734370902, -0.4838350155], [0.4941094279, -0.4448296300, 0.7469822445], [-0.8676661490, -0.1980763734, 0.4559837762]])

def dl(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3+1-OM); return (1+z)*C/70*np.array([quad(e, 0, zi)[0] for zi in z])

def unit(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec); return np.stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)], -1)

_cache = {}
def bayestar_far(c, path):
    """Bayestar 2019 best-fit reddening in the farthest distance bin (cumulative to 60 kpc), multi-resolution nested pixels."""
    if 'bayestar' not in _cache:
        with h5py.File(path) as f:
            pi = f['pixel_info'][:]; far = f['best_fit'][:, -1]
        tabs = {}
        for ns in np.unique(pi['nside']):
            sel = pi['nside'] == ns; order = np.argsort(pi['healpix_index'][sel])
            tabs[int(ns)] = (pi['healpix_index'][sel][order], far[sel][order])
        _cache['bayestar'] = tabs
    g = c.galactic; out = np.full(len(c), np.nan)
    for ns, (keys, vals) in _cache['bayestar'].items():
        idx = hp.ang2pix(ns, g.l.deg, g.b.deg, nest=True, lonlat=True)
        pos = np.searchsorted(keys, idx); pos = np.clip(pos, 0, len(keys)-1); hit = keys[pos] == idx
        out[hit & np.isnan(out)] = vals[pos[hit & np.isnan(out)]]
    return out

def lenz(c, path):
    from astropy.io import fits
    if 'lenz' not in _cache: _cache['lenz'] = fits.getdata(path, 1)['EBV'].astype(float)
    g = c.galactic
    return _cache['lenz'][hp.ang2pix(1024, g.l.deg, g.b.deg, nest=False, lonlat=True)]

def maps(c, data_dir):
    from dustmaps.config import config; config['data_dir'] = data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.planck import PlanckGNILCQuery
    out = {'SFD': SFDQuery()(c), 'GNILC': PlanckGNILCQuery()(c)}
    try: out['Bayestar19'] = bayestar_far(c, f"{data_dir}/bayestar/bayestar2019.h5")
    except Exception as e: print("no Bayestar map:", repr(e)[:150])
    try: out['Lenz17'] = lenz(c, f"{data_dir}/lenz2017/ebv_lhd.hpx.fits")
    except Exception as e: print("no Lenz map:", repr(e)[:150])
    try:
        from dustmaps.edenhofer2023 import Edenhofer2023Query
        q = Edenhofer2023Query(integrated=True); out['Edenhofer'] = q(SkyCoord(c.ra, c.dec, distance=1.24*u.kpc, frame='icrs'))
    except Exception as e: print("no Edenhofer map:", repr(e)[:150])
    return out

def rl(z, lam):
    """F99 R_V=3.1 extinction per unit E(B-V) at rest wavelength lam seen at redshift z (observer frame)."""
    return np.array([extinction.fitzpatrick99(np.array([lam*(1+zi)]), 3.1, 3.1)[0] for zi in z])

# SALT2 response to a unit change of the fixed Milky-Way E(B-V) in the fit (dust_salt_response.py, sncosmo SALT2 T23, F99 dust,
# gri light curves; BVRI and ZTF gri agree to 0.05): d mB, d c, d x1 per unit E(B-V), by redshift. The analytic
# R(lamB) - beta [R(lamB) - R(lamV)] overstates the response because SALT2's colour term absorbs 1.1, not 1.0, per unit E(B-V).
SALT_Z = np.array([0.02, 0.04, 0.07, 0.10, 0.20, 0.40]); SALT_DMB = np.array([4.07, 3.97, 3.83, 3.70, 3.27, 2.57])
SALT_DC = np.array([1.10, 1.07, 1.03, 1.00, 0.88, 0.65]); SALT_DX1 = np.array([0.21, 0.20, 0.17, 0.08, 0.03, -0.01])
def salt_response(z, beta, alpha=0.148):
    """|dmu/dE|: how much a distance modulus moves (brighter) per unit of E(B-V) the fit under-corrected, at fixed alpha, beta."""
    zz = np.clip(z, SALT_Z[0], SALT_Z[-1])
    return np.interp(zz, SALT_Z, SALT_DMB) - beta*np.interp(zz, SALT_Z, SALT_DC) - alpha*np.interp(zz, SALT_Z, SALT_DX1)

class Sample:
    """y: residual; w(X, sel): whitening of a column stack on a subset; base: base columns; extras for the correction."""
    def __init__(self, label, y, w, base, c, z, applied_ebv, beta):
        self.label, self.y, self.w, self.base, self.c, self.z, self.applied, self.beta = label, y, w, base, c, z, applied_ebv, beta
        self.b = c.galactic.b.deg
    def gls(self, y, cols, sel=None):
        sel = np.ones(len(y), bool) if sel is None else sel
        X = np.column_stack([b[sel] for b in self.base] + [cc[sel] for cc in cols]); Xw, yw = self.w(X, sel), self.w(y[sel], sel)
        b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); e = np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))); r = yw-Xw@b
        return b[len(self.base):], e[len(self.base):], r@r
    def slope(self, y, cols, sel=None):
        _, _, c0 = self.gls(y, [], sel); b, e, c1 = self.gls(y, cols, sel); return b, e, c0-c1

def part_a(S, mp, tag=""):
    """For each emission-free map: slope on (k MAP - SFD), GNILC-SFD alone, both; and the slope the emission-free difference
    would show if GNILC were the truth and the map a noisy measure of it: slope_GNILC x cov(x_map, x_GNILC)/var(x_map)."""
    s, g = mp['SFD'+tag], mp['GNILC'+tag]; comb = {}
    for name in [n for n in ('Bayestar19', 'Lenz17', 'Edenhofer') if n+tag in mp]:
        m = mp[name+tag]; ok = np.isfinite(m) & (m > -0.5) & (m < 5)
        for bl, sel in [("all", ok), ("|b|>20", ok & (np.abs(S.b) > 20)), ("|b|>30", ok & (np.abs(S.b) > 30))]:
            if sel.sum() < 50: continue
            k = np.sum(s[sel]*m[sel])/np.sum(m[sel]**2); x = k*m - s; xg = g - s
            b1, e1, d1 = S.slope(S.y, [x], sel); b3, e3, d3 = S.slope(S.y, [xg], sel); b2, e2, _ = S.slope(S.y, [x, xg], sel)
            lam = np.cov(x[sel], xg[sel])[0, 1]/np.var(x[sel]); exp = b3[0]*lam
            print(f"  {S.label} {name+tag:14s} {bl:7s} n={sel.sum():4d} k={k:.3f} rms(kMAP-SFD)={np.std(x[sel]):.4f} rms(GNILC-SFD)={np.std(xg[sel]):.4f}"
                  f" | kMAP-SFD {b1[0]:5.2f}±{e1[0]:4.2f} (p={chi2d.sf(d1,1):.2g}; expected if GNILC true {exp:5.2f}) | GNILC-SFD alone {b3[0]:5.2f}±{e3[0]:4.2f}"
                  f" | joint {b2[0]:5.2f}±{e2[0]:4.2f}, {b2[1]:5.2f}±{e2[1]:4.2f} | corr {np.corrcoef(x[sel], xg[sel])[0,1]:.2f}")
            comb[(name+tag, bl)] = (b1[0], e1[0], exp)
    return comb

def smooth(c, mp, data_dir, R=60.0):
    """Disk averages (radius R arcmin: centre + rings at R/3, 2R/3, R with 12, 20, 28 points) of SFD, GNILC, Bayestar, Lenz."""
    from dustmaps.config import config; config['data_dir'] = data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.planck import PlanckGNILCQuery
    pts = [c]
    for r, npt in ((R/3, 12), (2*R/3, 20), (R, 28)):
        for th in np.linspace(0, 2*np.pi, npt, endpoint=False): pts.append(c.directional_offset_by(th*u.rad, r*u.arcmin))
    qs = {'SFD': SFDQuery(), 'GNILC': PlanckGNILCQuery()}
    if 'Bayestar19' in mp: qs['Bayestar19'] = lambda cc: bayestar_far(cc, f"{data_dir}/bayestar/bayestar2019.h5")
    if 'Lenz17' in mp: qs['Lenz17'] = lambda cc: lenz(cc, f"{data_dir}/lenz2017/ebv_lhd.hpx.fits")
    out = {}
    for name, q in qs.items():
        v = np.array([q(p) for p in pts]); out[name+'_sm'] = np.nanmean(v, axis=0) if name != 'Lenz17' else np.where(np.isfinite(v).all(0), np.nanmean(v, axis=0), np.nan)
    return out

def part_t(S, mp, data_dir):
    """Is GNILC - SFD the dust temperature? Residual on the GNILC temperature at the SN, alone and with GNILC - SFD."""
    import os; from astropy.io import fits
    path = f"{data_dir}/planck/COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits"
    if not os.path.exists(path): print("  no GNILC temperature map"); return
    T = fits.getdata(path, 1)['TEMP'].astype(float); g = S.c.galactic
    t = T[hp.ang2pix(2048, g.l.deg, g.b.deg, nest=False, lonlat=True)]
    s, gn = mp['SFD'], mp['GNILC']; x = gn - s
    print(f"  {S.label} GNILC T at the SNe: median {np.median(t):.2f} K, rms {np.std(t):.2f}; corr(T, GNILC-SFD) {np.corrcoef(t, x)[0,1]:.2f}, corr(T, SFD) {np.corrcoef(t, s)[0,1]:.2f}")
    b, e, d = S.slope(S.y, [t - np.median(t)]); print(f"    residual on T-median: {b[0]*1e3:6.2f}±{e[0]*1e3:4.2f} mmag/K, dchi2 {d:.1f}")
    b, e, d = S.slope(S.y, [x, t - np.median(t)]); print(f"    joint GNILC-SFD {b[0]:5.2f}±{e[0]:4.2f}, T {b[1]*1e3:6.2f}±{e[1]*1e3:4.2f} mmag/K")
    b, e, d = S.slope(S.y, [x*(t - np.median(t))]); print(f"    residual on (GNILC-SFD) x (T-median): {b[0]:6.2f}±{e[0]:4.2f} per E(B-V) K")

def part_c(S, mp):
    s, g = mp['SFD'], mp['GNILC']; dE = g - S.applied
    rB, rV = rl(S.z, 4400.0), rl(S.z, 5500.0)
    # dE > 0 (GNILC says more dust than was applied): the SN is under-corrected, too faint; the correction brightens it,
    # mu' = mu - R_B dE + beta dc, so the corrected residual is y - (R_B - beta (R_B - R_V)) dE
    for beta_lbl, coef in [("analytic F99, beta=3.1", rB - 3.1*(rB - rV)), (f"SALT2 response, sample beta={S.beta:.2f}", salt_response(S.z, S.beta))]:
        dmu = coef*dE; yc = S.y - dmu
        print(f"  {S.label} GNILC correction ({beta_lbl}): dmu per unit dE: median {np.median(coef):.2f} (range {np.min(coef):.2f} to {np.max(coef):.2f}); rms dmu {np.std(dmu):.4f} mag")
        for lbl, x in [("GNILC-SFD", g-s), ("GNILC-applied", dE)]:
            b0, e0, d0 = S.slope(S.y, [x]); b1, e1, d1 = S.slope(yc, [x])
            print(f"    slope on {lbl:14s} before {b0[0]:5.2f}±{e0[0]:4.2f} (dchi2 {d0:5.1f}) -> after {b1[0]:5.2f}±{e1[0]:4.2f} (dchi2 {d1:5.1f})")
        for bl, sel in [("|b|>30", np.abs(S.b) > 30), ("|b|<30", np.abs(S.b) < 30), ("z<0.1", S.z < 0.1), ("z>=0.1", S.z >= 0.1)]:
            if sel.sum() < 50: continue
            b0, e0, _ = S.slope(S.y, [g-s], sel); b1, e1, _ = S.slope(yc, [g-s], sel)
            print(f"    {bl:7s} n={sel.sum():4d} GNILC-SFD slope before {b0[0]:5.2f}±{e0[0]:4.2f} -> after {b1[0]:5.2f}±{e1[0]:4.2f}")
    # how much of the correction the SNe want: the slope of the residual on the predicted correction itself (1 = exactly the swap)
    dmu1 = salt_response(S.z, S.beta)*dE; b, e, d = S.slope(S.y, [dmu1])
    print(f"  {S.label} slope of the residual on the SALT2-response GNILC correction (1 = the SNe want exactly the GNILC-for-SFD swap, 0 = none of it): {b[0]:5.2f}±{e[0]:4.2f}, dchi2 {d:.1f}")
    return dmu1, (b[0], e[0])

def part_d(S, dmu):
    n_eq = unit(S.c.ra.deg, S.c.dec.deg); nn = n_eq@R_gal.T; sb = nn[:, 2]
    cmb = unit(167.942, -6.944); sc = n_eq@cmb
    for lbl, y in [("before", S.y), ("GNILC-corrected", S.y - dmu)]:
        print(f"  {S.label} {lbl}:")
        for tname, cols, dof in [('galactic sin b dipole', [sb], 1), ('galactic P2(sin b) quad', [1.5*sb*sb-0.5], 1), ('full dipole (gal xyz)', [nn[:, 0], nn[:, 1], nn[:, 2]], 3),
                                 ('CMB-axis dipole', [sc], 1), ('CMB-axis quad', [1.5*sc*sc-0.5], 1)]:
            b, e, d = S.slope(y, cols); print(f"    {tname:26s} dchi2={d:6.2f} dof={dof} p={chi2d.sf(d, dof):.3f} coef={np.round(b, 4)} ±{np.round(e, 4)}")

def main(dat, covf, ztfcsv, data_dir):
    # Pantheon+
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f: n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    msk = (d.IS_CALIBRATOR == 0) & (d.zHD > 0.01); idx = np.where(msk)[0]
    z = d.zHD.values[idx]; y = d.m_b_corr.values[idx] - (5*np.log10(dl(z))+25)
    c = SkyCoord(d.RA.values[idx]*u.deg, d.DEC.values[idx]*u.deg, frame='icrs'); mp = maps(c, data_dir)
    print(f"Pantheon+: {len(idx)} SNe; MWEBV/SFD median {np.median(d.MWEBV.values[idx]/mp['SFD']):.3f}; maps: {sorted(mp)}")
    for k_, v in mp.items(): print(f"  {k_:10s} finite {np.isfinite(v).mean()*100:5.1f}%  median {np.nanmedian(v):.4f}  corr with SFD {np.corrcoef(v[np.isfinite(v)], mp['SFD'][np.isfinite(v)])[0,1]:.3f}")
    Lfull = cov[np.ix_(idx, idx)]
    def w(X, sel):
        Lc = np.linalg.cholesky(Lfull[np.ix_(sel, sel)]); return np.linalg.solve(Lc, X)
    P = Sample("P+ ", y, w, [np.ones(len(y))], c, z, d.MWEBV.values[idx], 3.04)
    # ZTF
    t = pd.read_csv(ztfcsv)
    sel = (t.sn_type == 'snia-cosmo') & (t.lccoverage_flag == 1) & (t.fitquality_flag == 1) & (t.source == 'z_gal') \
        & (t.x1.abs() < 3) & (t.c.abs() < 0.3) & (t.redshift > 0.015) & (t.redshift < 0.1) & t.x0.notna()
    t = t[sel].reset_index(drop=True)
    ct = SkyCoord(t.ra.values*u.deg, t.dec.values*u.deg, frame='icrs'); apex = SkyCoord(l=264.021*u.deg, b=48.253*u.deg, frame='galactic')
    zc = (1+t.redshift.values)/(1 - 369.82*np.cos(ct.galactic.separation(apex).rad)/C) - 1
    mb = -2.5*np.log10(t.x0.values); x1 = t.x1.values; cc = t.c.values
    var = (2.5/np.log(10)*t.x0_err.values/t.x0.values)**2 + AL**2*t.x1_err.values**2 + BE**2*t.c_err.values**2 \
        + 2*(-AL)*(-2.5/np.log(10)/t.x0.values)*t.cov_x0_x1.values + 2*BE*(-2.5/np.log(10)/t.x0.values)*t.cov_x0_c.values \
        + 2*(-AL)*BE*t.cov_x1_c.values + (5/np.log(10)*300/(C*zc))**2 + (5/np.log(10)*t.redshift_err.values/zc)**2 + 0.078**2
    yz = mb - (5*np.log10(dl(zc))+25); wz = 1/np.sqrt(var)
    def w2(X, sel): return X*wz[sel][:, None] if X.ndim == 2 else X*wz[sel]
    mpz = maps(ct, data_dir)
    Z = Sample("ZTF", yz, w2, [np.ones(len(yz)), -x1, cc], ct, zc, t.mwebv.values, 3.1)
    bz, ez, _ = Z.gls(yz, []); Z.beta = bz[2] if False else 3.1  # base has alpha, beta as columns; the correction uses 3.1 and the fitted value below
    bfit = np.linalg.lstsq(np.column_stack(Z.base)*wz[:, None], yz*wz, rcond=None)[0]; Z.beta = bfit[2]
    print(f"ZTF DR2: {len(t)} SNe; mwebv/SFD median {np.median(t.mwebv.values/mpz['SFD']):.3f}; fitted alpha {bfit[1]:.3f} beta {bfit[2]:.3f}")
    for k_, v in mpz.items(): print(f"  {k_:10s} finite {np.isfinite(v).mean()*100:5.1f}%  median {np.nanmedian(v):.4f}")

    print("\n(a) emission-free maps: residual on (k MAP - SFD), GNILC-SFD alone, and both; 'expected' = the slope the map difference"
          " would carry if GNILC were right and the map a noisy measure of the same truth (attenuation by cov/var)")
    ca = part_a(P, mp); cz = part_a(Z, mpz)
    print("  the same on 1-degree disk averages of every map (the GNILC-SFD slope lives at >40'; averaging beats down the stellar maps' pixel noise):")
    mp.update(smooth(c, mp, data_dir)); mpz.update(smooth(ct, mpz, data_dir))
    for k_ in ('SFD', 'GNILC', 'Bayestar19', 'Lenz17'):
        if k_+'_sm' in mp: print(f"    P+ {k_:10s} rms point-minus-smoothed {np.nanstd(mp[k_]-mp[k_+'_sm']):.4f}")
    ca.update(part_a(P, mp, '_sm')); cz.update(part_a(Z, mpz, '_sm'))
    for key in ca:
        if key in cz:
            (b1, e1, x1_), (b2, e2, x2_) = ca[key], cz[key]; wt = np.array([1/e1**2, 1/e2**2]); m = (b1*wt[0]+b2*wt[1])/wt.sum(); e = 1/np.sqrt(wt.sum())
            ex = (x1_*wt[0]+x2_*wt[1])/wt.sum()
            print(f"  combined {key[0]:14s} {key[1]:7s}: {m:5.2f}±{e:4.2f} ({m/e:.1f} sigma from 0; expected if GNILC true {ex:5.2f}, {(ex-m)/e:.1f} sigma above the measurement)")
    print("\n(t) GNILC dust temperature at the SN")
    part_t(P, mp, data_dir); part_t(Z, mpz, data_dir)
    print("\n(c) GNILC-corrected distance moduli (dE = GNILC - applied E(B-V), F99 at observer-frame rest-B and rest-V)")
    (dmuP, rP), (dmuZ, rZ) = part_c(P, mp), part_c(Z, mpz)
    wt = np.array([1/rP[1]**2, 1/rZ[1]**2]); m = (rP[0]*wt[0]+rZ[0]*wt[1])/wt.sum(); e = 1/np.sqrt(wt.sum())
    print(f"  combined: the SNe want {m:.2f}±{e:.2f} times the GNILC-for-SFD swap ({(m-1)/e:.1f} sigma above 1, {m/e:.1f} sigma above 0)")
    print("\n(d) Galactic-axis and CMB-axis anisotropy terms, Pantheon+ (zHD), before and after the GNILC correction")
    part_d(P, dmuP)
    print("  ZTF (same terms; northern sample, little leverage on a Galactic quadrupole):"); part_d(Z, dmuZ)

if __name__ == '__main__': main(*sys.argv[1:5])
