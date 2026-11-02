"""CHIME Catalog 2: is the Galactic north-south asymmetry of the extragalactic DM (R-115 row 10) a property of the sky or of
the instrument? At fixed declination the two Galactic hemispheres share CHIME's exposure and beam, so compare DM_exc between
hemispheres within Dec bands, with |b| controlled; then the foreground check on the localized FRBs (Macquart residuals by
hemisphere). usage: python frb_ns_dec.py data/chimefrbcat2.csv data/frb_public_hosts.csv data/frb_json"""
import sys, glob, json, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.stats import mannwhitneyu
rng = np.random.default_rng(1)

def boot_median_diff(a, b, n=2000):
    a, b = np.asarray(a), np.asarray(b)
    d = np.median(a) - np.median(b)
    bs = [np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(n)]
    return d, np.std(bs)

def quantile_reg(y, X, q=0.5, iters=200):
    """Median regression by iteratively reweighted least squares (L1 via IRLS); errors by bootstrap outside."""
    w = np.ones(len(y)); b = np.linalg.lstsq(X, y, rcond=None)[0]
    for _ in range(iters):
        r = y - X @ b; w = 1 / np.maximum(np.abs(r), 1e-3)
        Xw = X * w[:, None]; b = np.linalg.solve(X.T @ Xw, Xw.T @ y)
    return b

def boot_qreg(y, X, n=300):
    bs = np.array([quantile_reg(y[i], X[i]) for i in (rng.integers(0, len(y), len(y)) for _ in range(n))])
    return quantile_reg(y, X), bs.std(0)

def chime(path):
    d = pd.read_csv(path)
    d = d[(d.excluded_flag == 0) & (d.sub_num == 0) & d.dec.notna()].copy()
    # one burst per repeater (the first by MJD); non-repeaters as they are
    rep = d.repeater_name.notna()
    first = d[rep].sort_values('mjd_400').drop_duplicates('repeater_name')
    d = pd.concat([d[~rep], first])
    d = d[(d.dm_exc_ne2001 > 0) & (d.dm_exc_ymw16 > 0)]
    d['absb'] = np.abs(d.gb); d['north'] = (d.gb > 0).astype(int)
    print(f"CHIME Cat 2 after cuts: {len(d)} (repeaters collapsed to {rep.sum()}->{len(first)})")
    for dmcol in ['dm_exc_ne2001', 'dm_exc_ymw16']:
        print(f"\n=== {dmcol} ===")
        for bcut in [10, 20, 30]:
            s = d[d.absb > bcut]
            n, so = s[s.north == 1][dmcol], s[s.north == 0][dmcol]
            diff, err = boot_median_diff(n, so)
            print(f"|b|>{bcut}: N={len(n)} median {np.median(n):.0f}  S={len(so)} median {np.median(so):.0f}  "
                  f"diff {diff:+.0f}±{err:.0f}  p10 N/S {np.percentile(n,10):.0f}/{np.percentile(so,10):.0f}  "
                  f"p90 {np.percentile(n,90):.0f}/{np.percentile(so,90):.0f}  MWU p={mannwhitneyu(n,so).pvalue:.1e}")
        # Dec bands at |b|>10 (south reaches |b|>10 only at Dec < ~53)
        print("Dec bands (|b|>10): band  N(n,median,|b|med)  S(n,median,|b|med)  diff±err  p10 diff")
        s = d[d.absb > 10]; bands = [(-11, 0), (0, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 60), (60, 90)]
        diffs, errs = [], []
        for lo, hi in bands:
            t = s[(s.dec >= lo) & (s.dec < hi)]; n, so = t[t.north == 1], t[t.north == 0]
            if len(n) < 15 or len(so) < 15:
                print(f"  {lo:+3d}..{hi:+3d}: N={len(n)} S={len(so)} (too few)"); continue
            diff, err = boot_median_diff(n[dmcol], so[dmcol]); diffs.append(diff); errs.append(err)
            print(f"  {lo:+3d}..{hi:+3d}: N=({len(n)},{np.median(n[dmcol]):.0f},{np.median(n.absb):.0f})  "
                  f"S=({len(so)},{np.median(so[dmcol]):.0f},{np.median(so.absb):.0f})  {diff:+.0f}±{err:.0f}  "
                  f"p10 {np.percentile(n[dmcol],10)-np.percentile(so[dmcol],10):+.0f}")
        w = 1 / np.array(errs) ** 2
        print(f"  stratified (inverse-variance over bands): {np.sum(w*diffs)/np.sum(w):+.0f} ± {1/np.sqrt(np.sum(w)):.0f}")
        # 2-D strata: Dec band x |b| band, |b| 10-30 and >30
        print("Dec x |b| strata: (dec band, |b| band): diff±err")
        diffs, errs = [], []
        for lo, hi in bands:
            for blo, bhi in [(10, 20), (20, 35), (35, 90)]:
                t = s[(s.dec >= lo) & (s.dec < hi) & (s.absb >= blo) & (s.absb < bhi)]
                n, so = t[t.north == 1][dmcol], t[t.north == 0][dmcol]
                if len(n) < 15 or len(so) < 15: continue
                diff, err = boot_median_diff(n, so); diffs.append(diff); errs.append(err)
                print(f"  ({lo:+3d}..{hi:+3d}, {blo}-{bhi}): N={len(n)} S={len(so)} {diff:+.0f}±{err:.0f}")
        w = 1 / np.array(errs) ** 2
        print(f"  stratified: {np.sum(w*diffs)/np.sum(w):+.0f} ± {1/np.sqrt(np.sum(w)):.0f}")
        # median regression of log DM_exc on hemisphere, csc|b|, Dec-band dummies, log exposure, log S/N
        t = s[s.exp_up.notna() & (s.exp_up > 0) & (s.bonsai_snr > 0)].copy()
        y = np.log10(t[dmcol].values)
        cols = [np.ones(len(t)), t.north.values, 1 / np.sin(np.radians(t.absb.values)),
                np.log10(t.exp_up.values), np.log10(t.bonsai_snr.values)]
        names = ['const', 'north', 'csc|b|', 'logExp', 'logSNR']
        for lo, hi in bands[1:]:
            cols.append(((t.dec >= lo) & (t.dec < hi)).values.astype(float)); names.append(f'dec{lo}')
        X = np.column_stack(cols); b, e = boot_qreg(y, X, 200)
        print("median regression of log10 DM_exc:", "  ".join(f"{nm} {bi:+.3f}±{ei:.3f}" for nm, bi, ei in zip(names, b, e)))
        print(f"  -> north term in DM at the median: x{10**b[1]:.2f} (±{np.log(10)*e[1]*100:.0f} %)")
        # the same without the Dec terms and exposure (the raw asymmetry), for comparison
        b2, e2 = boot_qreg(y, X[:, :3], 200)
        print(f"  raw (hemisphere + csc|b| only): north {b2[1]:+.3f}±{e2[1]:.3f} -> x{10**b2[1]:.2f}")
        # where on the sky: the north split by longitude (anticentre l 90-270 vs centre), and the eROSITA-bubble footprint
        t = s[s.north == 1]
        bub = (np.abs(((t.gl + 180) % 360) - 180) < 30) & (t.gb > 10) & (t.gb < 80)
        anti = (t.gl > 90) & (t.gl < 270)
        print(f"  north: bubble footprint n={bub.sum()} median {np.median(t[bub][dmcol]):.0f} vs rest of north "
              f"{np.median(t[~bub][dmcol]):.0f}; anticentre half (l 90-270) {np.median(t[anti][dmcol]):.0f} vs centre half "
              f"{np.median(t[~anti][dmcol]):.0f}")
        so = s[s.north == 0]; anti = (so.gl > 90) & (so.gl < 270)
        print(f"  south: anticentre half {np.median(so[anti][dmcol]):.0f} vs centre half {np.median(so[~anti][dmcol]):.0f}")
        # scattering and width as a selection proxy by hemisphere
        for col in ['scat_time', 'bc_width', 'bonsai_snr', 'exp_up']:
            v = s[s[col] > 0]; print(f"  {col}: median N {np.median(v[v.north==1][col]):.3g}  S {np.median(v[v.north==0][col]):.3g}")
    return d

def macquart(z, f_igm=0.83, H0=67.7, Om=0.31, Ob=0.049):
    """mean cosmic DM in pc/cm3 (Macquart 2020 form, chi_e = 7/8)."""
    pref = 3 * 299792.458 * H0 * Ob / (8 * np.pi * 6.674e-11 * 1.673e-27) * (7 / 8)  # in SI-ish units, fixed below
    # do it numerically with astropy-free constants: DM = n_e0 c/H0 ∫ (1+z)/E(z) dz, n_e0 = rho_crit Ob chi / m_p
    rho_crit = 1.878e-26 * (H0 / 100) ** 2  # kg/m3
    n_e0 = rho_crit * Ob * (7 / 8) / 1.673e-27  # m^-3
    cH0 = 299792.458e3 / (H0 * 1e3 / 3.0857e22)  # m
    E = lambda x: np.sqrt(Om * (1 + x) ** 3 + 1 - Om)
    I = quad(lambda x: (1 + x) / E(x), 0, z)[0]
    return f_igm * n_e0 * cH0 * I / 3.0857e16 * 1e-6  # m^-3 m -> pc cm^-3

def localized(hosts, jdir):
    h = pd.read_csv(hosts); rows = []
    for _, r in h.iterrows():
        try: j = json.load(open(f"{jdir}/FRB{r.FRB}.json"))
        except FileNotFoundError: continue
        if 'z' not in j or 'DM' not in j: continue
        rows.append(dict(frb=r.FRB, z=float(j['z']), dm=j['DM']['value'], dmism=j['DMISM']['value'], ra=j['ra'], dec=j['dec'],
                         rep=bool(j.get('repeater', False)), proj=str(r.Projects)))
    d = pd.DataFrame(rows)
    from astropy.coordinates import SkyCoord
    c = SkyCoord(d.ra.values, d.dec.values, unit='deg').galactic; d['l'], d['b'] = c.l.deg, c.b.deg
    d['dm_cos'] = [macquart(z) for z in d.z]
    d['res'] = d.dm - d.dmism - 50 - 50 / (1 + d.z) - d.dm_cos  # halo 50, host 50/(1+z)
    d = d[(np.abs(d.b) > 10) & (d.z > 0.005)]
    print(f"\n=== localized FRBs with z: {len(d)} at |b|>10 (CHIME-found: {d.proj.str.contains('CHIME').sum()}) ===")
    for lab, m in [('all', np.ones(len(d), bool)), ('non-CHIME', ~d.proj.str.contains('CHIME').values),
                   ('z<0.3', (d.z < 0.3).values)]:
        s = d[m]; n, so = s[s.b > 0], s[s.b < 0]
        diff, err = boot_median_diff(n.res, so.res)
        print(f"{lab}: N={len(n)} median res {np.median(n.res):+.0f} (mean {n.res.mean():+.0f}±{n.res.std()/np.sqrt(len(n)):.0f})  "
              f"S={len(so)} {np.median(so.res):+.0f} (mean {so.res.mean():+.0f}±{so.res.std()/np.sqrt(len(so)):.0f})  "
              f"median diff {diff:+.0f}±{err:.0f}  MWU p={mannwhitneyu(n.res, so.res).pvalue:.2f}")
    # regression of the residual on hemisphere and csc|b|
    X = np.column_stack([np.ones(len(d)), (d.b > 0).astype(float), 1 / np.sin(np.radians(np.abs(d.b)))])
    b, *_ = np.linalg.lstsq(X, d.res.values, rcond=None); r = d.res.values - X @ b
    cov = np.linalg.inv(X.T @ X) * (r @ r / (len(d) - 3))
    print(f"OLS res ~ north + csc|b|: north {b[1]:+.0f}±{np.sqrt(cov[1,1]):.0f}, csc|b| {b[2]:+.0f}±{np.sqrt(cov[2,2]):.0f} pc/cm3")

if __name__ == '__main__':
    chime(sys.argv[1]); localized(sys.argv[2], sys.argv[3])
