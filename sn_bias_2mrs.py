"""Clustering bias of Type Ia supernovae relative to galaxies in the local universe (sweep-4 row B1; the open part of the
trajectory verdict): the projected SN-galaxy cross-correlation w_p,sg(r_p) of the ZTF SN Ia DR2 volume-limited sample
(snia-cosmo, light-curve and fit quality flags, 0.015 < z < 0.06, |b| > 10) against 2MRS galaxies (K < 11.75), divided by the
projected auto-correlation of a 2MRS control sample weighted to the SN redshift distribution ("pseudo-SNe"), on 2-20 Mpc/h;
b_SN / b_2MRS = w_p,sg / w_p,pg. Landy-Szalay with randoms that follow the 2MRS mask (Dec > -28, |b| > 10) and n(z), and for
the SNe the SN sky density smoothed with a 30 deg Gaussian (ZTF's depth is not uniform); jackknife over nside-2 sky regions.
Splits: x1 < -0.5 vs >= -0.5, host log mass < 10 vs >= 10, host-spectrum redshifts only; a K-luminosity-weighted galaxy
control (SNe trace stellar mass); uniform-sky SN randoms as a check.
usage: python -I sn_bias_2mrs.py data/ztf/ztfsniadr2_lite/tables/snia_data.csv data/ztf/ztfsniadr2_lite/tables/globalhost_data.csv data/2mrs/2mrs.tsv [nrand_factor] [vlonly]"""
import sys, numpy as np, pandas as pd, healpy as hp
from scipy.spatial import cKDTree
from scipy.integrate import quad
rng = np.random.default_rng(2)
C = 299792.458; OM = 0.3
RP_EDGES = np.array([0.5, 1, 2, 3.5, 6, 10, 20]); PI_MAX = 30.0; NPI = 15
DEC_MIN, B_MIN = -28.0, 10.0

def comoving(z):
    e = lambda x: 1 / np.sqrt(OM * (1 + x) ** 3 + 1 - OM)
    zz = np.linspace(0, 0.1, 401); dd = np.array([quad(e, 0, zi)[0] for zi in zz]) * C / 100  # Mpc/h
    return np.interp(z, zz, dd)

def unit(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], -1)

def gal_lb(ra, dec):
    r = hp.Rotator(coord=['C', 'G']); th, ph = r(np.radians(90 - dec), np.radians(ra)); return np.degrees(ph), 90 - np.degrees(th)

APEX = unit(167.942, -6.944)
def z_cmb(z_hel, ra, dec):
    return z_hel + 369.82 / C * (unit(ra, dec) @ APEX)

def vizier_tsv(path):
    rows = [l for l in open(path, encoding='utf-8', errors='ignore') if not l.startswith('#') and l.strip()]
    hdr = rows[0].rstrip('\n').split('\t'); body = [r.rstrip('\n').split('\t') for r in rows[3:]]
    return pd.DataFrame(body, columns=hdr).replace('', np.nan)

def in_mask(ra, dec):
    l, b = gal_lb(ra, dec); return (dec > DEC_MIN) & (np.abs(b) > B_MIN)

def random_sky(n):
    out = []
    while sum(len(o) for o in out) < n:
        ra = rng.uniform(0, 360, 2 * n); dec = np.degrees(np.arcsin(rng.uniform(-1, 1, 2 * n)))
        m = in_mask(ra, dec); out.append(np.column_stack([ra[m], dec[m]]))
    return np.vstack(out)[:n]

def pair_counts(xa, wa, ra_, xb, wb, rb_, nreg, same=False, nta=1):
    """Weighted pair counts in (r_p, |pi|) bins, tagged by the jackknife region of each member: H[tag_a, reg_b, rp, pi],
    where tag_a = reg_a * nta + (an optional extra tag of a, e.g. its redshift bin, carried in ra_ as reg*nta + zbin).
    A tree on b, chunks of a. For same=True (auto), self pairs are dropped."""
    tree = cKDTree(xb); R = np.sqrt(RP_EDGES[-1] ** 2 + PI_MAX ** 2)
    nrp = len(RP_EDGES) - 1; H = np.zeros((nreg * nta * nreg, nrp * NPI)); pe = np.linspace(0, PI_MAX, NPI + 1)
    for i0 in range(0, len(xa), 1000):
        A = xa[i0:i0 + 1000]; nb = tree.query_ball_point(A, R)
        ia = np.repeat(np.arange(len(A)), [len(n) for n in nb]); ib = np.concatenate([np.asarray(n, int) for n in nb]) if len(ia) else np.array([], int)
        if same:
            keep = (ib != ia + i0); ia, ib = ia[keep], ib[keep]
        s = xb[ib] - A[ia]; mid = xb[ib] + A[ia]; mid /= np.linalg.norm(mid, axis=1)[:, None]
        pi = np.abs((s * mid).sum(1)); rp = np.sqrt(np.maximum((s * s).sum(1) - pi * pi, 0))
        irp = np.searchsorted(RP_EDGES, rp, side='right') - 1; ipi = np.minimum((pi / (PI_MAX / NPI)).astype(int), NPI)
        ok = (irp >= 0) & (irp < nrp) & (ipi < NPI)
        cell = (ra_[i0:i0 + 1000][ia[ok]] * nreg + rb_[ib[ok]]) * (nrp * NPI) + irp[ok] * NPI + ipi[ok]
        H += np.bincount(cell, weights=wa[i0:i0 + 1000][ia[ok]] * wb[ib[ok]], minlength=H.size).reshape(H.shape)
    return H.reshape(nreg, nta, nreg, nrp, NPI)

class Counts:
    """The four LS pair-count histograms, region-tagged, with the per-region weight sums; w_p for the full sample and for
    each delete-one-region jackknife sample."""
    def __init__(self, d1, w1, r1_, d2, w2, r2_, R1, wR1, rR1_, R2, wR2, rR2_, nreg, same=False, nta=1, t1=None, tR1=None):
        """t1, tR1: the extra tag (0..nta-1) of the first members (data and randoms), e.g. a redshift bin."""
        self.nreg, self.nta = nreg, nta
        t1 = np.zeros(len(d1), int) if t1 is None else t1; tR1 = np.zeros(len(R1), int) if tR1 is None else tR1
        a1, aR1 = r1_ * nta + t1, rR1_ * nta + tR1
        self.DD = pair_counts(d1, w1, a1, d2, w2, r2_, nreg, same, nta); self.DR = pair_counts(d1, w1, a1, R2, wR2, rR2_, nreg, False, nta)
        self.RD = pair_counts(R1, wR1, aR1, d2, w2, r2_, nreg, False, nta); self.RR = pair_counts(R1, wR1, aR1, R2, wR2, rR2_, nreg, same, nta)
        bc = lambda r, w: np.bincount(r, weights=w, minlength=nreg)
        bc2 = lambda r, t, w: np.bincount(r * nta + t, weights=w, minlength=nreg * nta).reshape(nreg, nta)
        self.W = [(bc2(r1_, t1, w1), bc(r2_, w2)), (bc2(r1_, t1, w1), bc(rR2_, wR2)), (bc2(rR1_, tR1, wR1), bc(r2_, w2)), (bc2(rR1_, tR1, wR1), bc(rR2_, wR2))]
    def wp(self, drop=None, wt=None):
        """wt: weights per extra tag applied to the first members (a redshift reweighting), default all ones."""
        wt = np.ones(self.nta) if wt is None else np.asarray(wt); out = []
        for H, (Wa, Wb) in zip([self.DD, self.DR, self.RD, self.RR], self.W):
            Hw = np.einsum('t,atbrp->abrp', wt, H); Waw = Wa @ wt
            if drop is None: h = Hw.sum((0, 1)); n = Waw.sum() * Wb.sum()
            else: h = Hw.sum((0, 1)) - Hw[drop].sum(0) - Hw[:, drop].sum(0) + Hw[drop, drop]; n = (Waw.sum() - Waw[drop]) * (Wb.sum() - Wb[drop])
            out.append(h / n)
        DD, DR, RD, RR = out; xi = (DD - DR - RD + RR) / RR
        return 2 * xi.sum(1) * (PI_MAX / NPI)

def main(snf, hostf, mrsf, nr=4, vlonly=False):
    sn = pd.read_csv(snf).merge(pd.read_csv(hostf), on='ztfname', how='left')
    sn = sn[(sn.sn_type == 'snia-cosmo') & (sn.lccoverage_flag == 1) & (sn.fitquality_flag == 1)].copy()
    sn['zc'] = z_cmb(sn.redshift.values, sn.ra.values, sn.dec.values)
    sn = sn[(sn.zc > 0.015) & (sn.zc < 0.06) & in_mask(sn.ra.values, sn.dec.values)].reset_index(drop=True)
    g = vizier_tsv(mrsf)
    for c in ['RAJ2000', 'DEJ2000', 'Kcmag', 'cz']: g[c] = pd.to_numeric(g[c], errors='coerce')
    g = g[g.cz.notna() & g.Kcmag.notna()].copy(); g['zc'] = z_cmb(g.cz.values / C, g.RAJ2000.values, g.DEJ2000.values)
    g = g[(g.zc > 0.003) & (g.zc < 0.07) & in_mask(g.RAJ2000.values, g.DEJ2000.values)].reset_index(drop=True)
    g['LK'] = 10 ** (-0.4 * (g.Kcmag - 5 * np.log10(comoving(g.zc) * 1e4) - 3.28))
    print(f"SNe: {len(sn)} (host-spectrum z: {(sn.source=='z_gal').sum()}; x1 < -0.5: {(sn.x1 < -0.5).sum()}; host mass known: {sn.mass.notna().sum()}); "
          f"2MRS galaxies in the mask at 0.003 < z < 0.07: {len(g)}")
    xs = unit(sn.ra.values, sn.dec.values) * comoving(sn.zc.values)[:, None]
    xg = unit(g.RAJ2000.values, g.DEJ2000.values) * comoving(g.zc.values)[:, None]
    # randoms: galaxies (mask x shuffled n(z)); SNe (mask x SN-density weighting x shuffled SN n(z))
    rg_sky = random_sky(nr * len(g)); rg_z = rng.choice(g.zc.values, nr * len(g))
    xrg = unit(rg_sky[:, 0], rg_sky[:, 1]) * comoving(rg_z)[:, None]
    NS = 16; dens = np.bincount(hp.ang2pix(NS, sn.ra.values, sn.dec.values, lonlat=True), minlength=hp.nside2npix(NS)).astype(float)
    sky_r = random_sky(200000); pixr = hp.ang2pix(NS, sky_r[:, 0], sky_r[:, 1], lonlat=True); maskmap = np.bincount(pixr, minlength=hp.nside2npix(NS)).astype(float)
    def sn_randoms(smooth_deg, n):
        if smooth_deg is None: wmap = np.ones_like(dens)
        else:
            sm_d = hp.smoothing(dens, fwhm=np.radians(smooth_deg * 2.355), verbose=False) if hasattr(hp, 'smoothing') else dens
            sm_m = hp.smoothing(maskmap, fwhm=np.radians(smooth_deg * 2.355), verbose=False)
            wmap = np.where(sm_m > 0, sm_d / np.maximum(sm_m, 1e-9), 0); wmap = np.maximum(wmap, 0)
        cand = random_sky(6 * n); w = wmap[hp.ang2pix(NS, cand[:, 0], cand[:, 1], lonlat=True)]
        keep = rng.uniform(0, w.max(), len(cand)) < w; cand = cand[keep][:n]
        return cand
    # pseudo-SN weights for galaxies: n_SN(z) / n_g(z)
    ze = np.linspace(0.015, 0.06, 10)
    # jackknife regions: nside-2 pixels, indexed over every pixel any object falls in
    pix = lambda ra, dec: hp.ang2pix(2, ra, dec, lonlat=True)
    allpix = np.unique(np.concatenate([pix(sn.ra.values, sn.dec.values), pix(g.RAJ2000.values, g.DEJ2000.values), pix(rg_sky[:, 0], rg_sky[:, 1])]))
    ridx = {p: i for i, p in enumerate(allpix)}; nreg = len(allpix)
    reg = lambda ra, dec: np.array([ridx[p] for p in pix(ra, dec)])
    reg_s = reg(sn.ra.values, sn.dec.values); reg_g = reg(g.RAJ2000.values, g.DEJ2000.values); reg_rg = reg(rg_sky[:, 0], rg_sky[:, 1])
    regions = np.unique(reg_s); n = len(regions); print(f"jackknife regions with SNe: {n} (of {nreg} pixels)")
    rpc = np.sqrt(RP_EDGES[1:] * RP_EDGES[:-1]); fit = rpc >= 2
    def jk(counts, wt=None):
        full = counts.wp(wt=wt); J = np.array([counts.wp(r, wt=wt) for r in regions]); err = np.sqrt((n - 1) / n * ((J - J.mean(0)) ** 2).sum(0)); return full, J, err
    def ratio(a, b, e): w = 1 / e[fit] ** 2; return (w * a[fit] * b[fit]).sum() / (w * b[fit] ** 2).sum()
    # the galaxy auto-correlation with the first member tagged by its redshift bin, once per control weighting; the
    # pseudo-SN weighting n_SN(z)/n_g(z) of each SN selection is then applied as a reweighting of the tags
    NZ = len(ze) - 1; zb_g = np.clip(np.searchsorted(ze, g.zc.values, side='right') - 1, 0, NZ - 1); zb_rg = np.clip(np.searchsorted(ze, rg_z, side='right') - 1, 0, NZ - 1)
    hg = np.histogram(g.zc, ze)[0]
    def zweights(zsel):
        hs = np.histogram(zsel, ze)[0]; return np.where(hg > 0, hs / np.maximum(hg, 1), 0.0)
    controls = {}
    for name, wg in ([] if vlonly else [('number', np.ones(len(g))), ('LK', g.LK.values / g.LK.values.mean())]):
        wr = np.ones(len(rg_z)) if name == 'number' else np.interp(rg_z, np.sort(g.zc.values), wg[np.argsort(g.zc.values)])
        controls[name] = Counts(xg, wg, reg_g, xg, np.ones(len(g)), reg_g, xrg, wr, reg_rg, xrg, np.ones(len(xrg)), reg_rg, nreg, same=True, nta=NZ, t1=zb_g, tR1=zb_rg)
        pg, Jpg, epg = jk(controls[name], zweights(sn.zc.values))
        print(f"control '{name}' weighted to the full SN n(z): w_p,pg: " + ' '.join(f'{v:6.1f}' for v in pg) + "  ± " + ' '.join(f'{v:4.1f}' for v in epg))
    if not vlonly: pgg, _, epgg = jk(controls['number']); print(f"2MRS auto-correlation unweighted (0.003 < z < 0.07): " + ' '.join(f'{v:6.1f}' for v in pgg) + "  ± " + ' '.join(f'{v:4.1f}' for v in epgg))
    print("r_p (Mpc/h): " + ' '.join(f'{r:6.1f}' for r in rpc))

    def measure(sel, control='number', smooth=30.0, label=''):
        s = sn[sel]; xs_ = xs[sel]
        rs_sky = sn_randoms(smooth, 40 * len(s)); rs_z = rng.choice(s.zc.values, len(rs_sky)); xrs = unit(rs_sky[:, 0], rs_sky[:, 1]) * comoving(rs_z)[:, None]
        cnt = Counts(xs_, np.ones(len(s)), reg_s[sel], xg, np.ones(len(g)), reg_g, xrs, np.ones(len(xrs)), reg(rs_sky[:, 0], rs_sky[:, 1]), xrg, np.ones(len(xrg)), reg_rg, nreg)
        sg, Jsg, esg = jk(cnt); pg, Jpg, epg = jk(controls[control], zweights(s.zc.values))
        A = ratio(sg, pg, esg); Ajk = np.array([ratio(Jsg[i], Jpg[i], esg) for i in range(n)]); eA = np.sqrt((n - 1) / n * ((Ajk - Ajk.mean()) ** 2).sum())
        small = ~fit; As = sg[small] / pg[small]; Asjk = Jsg[:, small] / Jpg[:, small]; eAs = np.sqrt((n - 1) / n * ((Asjk - Asjk.mean(0)) ** 2).sum(0))
        print(f"\n[{label}] N_SN={len(s)}")
        print(f"   w_p,sg: " + ' '.join(f'{v:6.1f}' for v in sg) + "  ± " + ' '.join(f'{v:4.1f}' for v in esg) + f"   (mean z {s.zc.mean():.4f})")
        print(f"   w_p,pg: " + ' '.join(f'{v:6.1f}' for v in pg) + "  ± " + ' '.join(f'{v:4.1f}' for v in epg) + "   (control weighted to this selection's n(z))")
        print(f"   ratio per bin: " + ' '.join(f'{v:6.2f}' for v in sg / pg))
        print(f"   **b_SN / b_control (2-20 Mpc/h) = {A:.3f} ± {eA:.3f}**;  0.5-2 Mpc/h: " + ', '.join(f'{a:.2f} ± {e:.2f}' for a, e in zip(As, eAs)))
        return A, eA

    allsel = np.ones(len(sn), bool)
    res = {}
    if not vlonly: res['all'] = measure(allsel, label='all SNe, 2MRS number-weighted control, 30 deg SN randoms')
    if not vlonly: res['LK'] = measure(allsel, control='LK', label='all SNe, K-luminosity-weighted control')
    if not vlonly: res['zgal'] = measure((sn.source == 'z_gal').values, label='host-spectrum redshifts only')
    if not vlonly: res['x1lo'] = measure((sn.x1 < -0.5).values, label='x1 < -0.5')
    if not vlonly: res['x1hi'] = measure((sn.x1 >= -0.5).values, label='x1 >= -0.5')
    if not vlonly: res['mlo'] = measure((sn.mass < 10).values, label='host log mass < 10')
    if not vlonly: res['mhi'] = measure((sn.mass >= 10).values, label='host log mass >= 10')
    if not vlonly: res['uniform'] = measure(allsel, smooth=None, label='all SNe, uniform-sky SN randoms')
    if not vlonly: res['sm15'] = measure(allsel, smooth=15.0, label='all SNe, 15 deg SN randoms')
    if not vlonly: res['x1lo_LK'] = measure((sn.x1 < -0.5).values, control='LK', label='x1 < -0.5, K-luminosity control')
    if not vlonly: res['x1hi_LK'] = measure((sn.x1 >= -0.5).values, control='LK', label='x1 >= -0.5, K-luminosity control')
    if not vlonly: res['z<0.04'] = measure((sn.zc < 0.04).values, label='z < 0.04')
    if not vlonly: res['z>0.04'] = measure((sn.zc >= 0.04).values, label='z >= 0.04')
    # a volume-limited control: M_K - 5 log h < -23.7 is complete to z = 0.04 at K < 11.75 (distance modulus 36.2); SNe at z < 0.04 against it,
    # so the galaxy population does not change with redshift inside the comparison volume
    MK = g.Kcmag.values - 5 * np.log10(comoving(g.zc.values) * 1e5)  # M_K - 5 log h (D in Mpc/h; K < 11.75 at z = 0.04 is -23.65)
    vl = (MK < -23.7) & (g.zc.values < 0.04); vl_r = np.interp(rg_z, np.sort(g.zc.values), vl[np.argsort(g.zc.values)].astype(float)) > 0.5
    print(f"\nvolume-limited control: {vl.sum()} galaxies with M_K - 5 log h < -23.7 at z < 0.04")
    cvl = Counts(xg[vl], np.ones(vl.sum()), reg_g[vl], xg[vl], np.ones(vl.sum()), reg_g[vl], xrg[vl_r], np.ones(vl_r.sum()), reg_rg[vl_r], xrg[vl_r], np.ones(vl_r.sum()), reg_rg[vl_r], nreg, same=True, nta=NZ, t1=zb_g[vl], tR1=zb_rg[vl_r])
    for label, sel in [('z < 0.04, all', (sn.zc < 0.04).values), ('z < 0.04, x1 < -0.5', ((sn.zc < 0.04) & (sn.x1 < -0.5)).values), ('z < 0.04, x1 >= -0.5', ((sn.zc < 0.04) & (sn.x1 >= -0.5)).values)]:
        s_ = sn[sel]; xs_ = xs[sel]
        rs_sky = sn_randoms(30.0, 40 * len(s_)); rs_z = rng.choice(s_.zc.values, len(rs_sky)); xrs = unit(rs_sky[:, 0], rs_sky[:, 1]) * comoving(rs_z)[:, None]
        # cross with the volume-limited galaxies only
        cnt = Counts(xs_, np.ones(len(s_)), reg_s[sel], xg[vl], np.ones(vl.sum()), reg_g[vl], xrs, np.ones(len(xrs)), reg(rs_sky[:, 0], rs_sky[:, 1]), xrg[vl_r], np.ones(vl_r.sum()), reg_rg[vl_r], nreg)
        hsz = np.histogram(s_.zc, ze)[0]; hgv = np.histogram(g.zc.values[vl], ze)[0]; wtv = np.where(hgv > 0, hsz / np.maximum(hgv, 1), 0.0)
        sg, Jsg, esg = jk(cnt); pg, Jpg, epg = jk(cvl, wtv)
        A = ratio(sg, pg, esg); Ajk = np.array([ratio(Jsg[i], Jpg[i], esg) for i in range(n)]); eA = np.sqrt((n - 1) / n * ((Ajk - Ajk.mean()) ** 2).sum())
        print(f"[{label}, volume-limited control] N_SN={len(s_)}")
        print(f"   w_p,sg: " + ' '.join(f'{v:6.1f}' for v in sg) + "  ± " + ' '.join(f'{v:4.1f}' for v in esg))
        print(f"   w_p,vv: " + ' '.join(f'{v:6.1f}' for v in pg) + "  ± " + ' '.join(f'{v:4.1f}' for v in epg))
        print(f"   **b_SN / b_VL (2-20 Mpc/h) = {A:.3f} ± {eA:.3f}**; per bin " + ' '.join(f'{v:5.2f}' for v in sg / pg))
    if vlonly: return
    # the SN auto-correlation (selection-independent of the galaxies), all SNe, against the number control
    rs_sky = sn_randoms(30.0, 40 * len(sn)); rs_z = rng.choice(sn.zc.values, len(rs_sky)); xrs = unit(rs_sky[:, 0], rs_sky[:, 1]) * comoving(rs_z)[:, None]; rrs = reg(rs_sky[:, 0], rs_sky[:, 1])
    css = Counts(xs, np.ones(len(sn)), reg_s, xs, np.ones(len(sn)), reg_s, xrs, np.ones(len(xrs)), rrs, xrs, np.ones(len(xrs)), rrs, nreg, same=True)
    ss, Jss, ess = jk(css); pg, Jpg, epg = jk(controls['number'], zweights(sn.zc.values))
    print(f"\n[SN auto-correlation, all {len(sn)} SNe] w_p,ss: " + ' '.join(f'{v:6.1f}' for v in ss) + "  ± " + ' '.join(f'{v:4.1f}' for v in ess))
    A = ratio(ss, pg, ess); Ajk = np.array([ratio(Jss[i], Jpg[i], ess) for i in range(n)]); eA = np.sqrt((n - 1) / n * ((Ajk - Ajk.mean()) ** 2).sum())
    print(f"   w_p,ss / w_p,pg (2-20 Mpc/h) = {A:.2f} ± {eA:.2f}, i.e. (b_SN / b_control)^2; sqrt = {np.sqrt(max(A, 0)):.2f} ± {eA / (2 * np.sqrt(max(A, 1e-3))):.2f}")
    print("\nsummary: " + '; '.join(f"{k} {v[0]:.2f} ± {v[1]:.2f}" for k, v in res.items()))
    d = res['x1lo'][0] - res['x1hi'][0]; e = np.hypot(res['x1lo'][1], res['x1hi'][1])
    print(f"x1 split difference {d:+.2f} ± {e:.2f} ({d/e:.1f} sigma); host-mass split {res['mhi'][0]-res['mlo'][0]:+.2f} ± {np.hypot(res['mhi'][1], res['mlo'][1]):.2f}")

if __name__ == '__main__':
    a = sys.argv[1:]; main(a[0], a[1], a[2], int(a[3]) if len(a) > 3 and a[3] != 'vlonly' else 4, vlonly='vlonly' in a)
