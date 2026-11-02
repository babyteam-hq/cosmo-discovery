"""Are the CHIME/FRB Catalog 2 "DM walls" (Ravi+2026, arXiv:2604.05093: a 150 pc/cm3, 4.3 sigma excess of mean
extragalactic DM over ~30 deg near RA 12h, Dec +55, read as baryons in the Ursa Major supercluster) a feature of the sky or
of the catalogue, and is the DM-excess field correlated with the local galaxy distribution?
Ravi+ already used one burst per source, RA-randomised and DM-scrambled nulls, S/N cuts and declination bins (the
trajectory referee's reading of the paper, results-R120.md). Part 1 rebuilds their smoothed map (Gaussian sigma 7.3 deg at
nside 8, bursts within 3 sigma, S/N >= 9, |b| > 10; NE2001/YMW16 here, NE2025 there) with all bursts and with one burst per
source, adds the look-elsewhere DM-shuffle null, and the one check they did not run: the map after DM_exc is regressed on
log S/N and declination. Part 2, the open part: the per-burst DM_exc and the smoothed map against the 2MRS K-band luminosity
density in distance shells, and DM_exc stacked around Tully 2015 groups.
usage: python -I frb_dm_walls.py data/chime/chimefrbcat2.csv data/2mrs/2mrs.tsv data/tully15/galaxies.tsv [nshuffle]
(data/tully15/galaxies.tsv is VizieR J/AJ/149/171 table3, the group table; the file name is the download's)"""
import sys, numpy as np, pandas as pd, healpy as hp
rng = np.random.default_rng(1)
NSIDE = 8; SIG = np.radians(7.3); RMAX = 3 * SIG
WALL1 = (180.0, 55.0); WALL2 = (30.0, 45.0)  # RA, Dec (deg)

def unit(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], -1)

def vizier_tsv(path):
    rows = [l for l in open(path, encoding='utf-8', errors='ignore') if not l.startswith('#') and l.strip()]
    hdr = rows[0].rstrip('\n').split('\t'); body = [r.rstrip('\n').split('\t') for r in rows[3:]]
    return pd.DataFrame(body, columns=hdr).replace('', np.nan)

def load_chime(path, one_per_source):
    d = pd.read_csv(path)
    d = d[(d.excluded_flag == 0) & (d.sub_num == 0) & d.dec.notna() & (d.bonsai_snr >= 9) & (np.abs(d.gb) > 10)].copy()
    d = d[d.dm_exc_ne2001 > 0]
    nrep = d.repeater_name.notna().sum()
    if one_per_source:
        rep = d.repeater_name.notna()
        first = d[rep].sort_values('mjd_400').drop_duplicates('repeater_name')
        d = pd.concat([d[~rep], first])
    d['src'] = d.repeater_name.fillna(d.tns_name)
    return d.reset_index(drop=True), nrep

class Smoother:
    """Gaussian-weighted mean of a per-burst quantity on nside-8 pixel centres (Ravi+'s kernel), with s.e.m. from the
    weighted scatter and n_eff = (sum w)^2 / sum w^2; pixels with < 50 bursts within 2 sigma are masked."""
    def __init__(self, ra, dec):
        npix = hp.nside2npix(NSIDE); th, ph = hp.pix2ang(NSIDE, np.arange(npix))
        cen = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], -1)
        ang = np.arccos(np.clip(cen @ unit(ra, dec).T, -1, 1))  # (npix, nburst)
        self.w = np.where(ang < RMAX, np.exp(-ang ** 2 / (2 * SIG ** 2)), 0.0)
        self.n2 = (ang < 2 * SIG).sum(1); self.ok = self.n2 >= 50
        self.sw = self.w.sum(1); self.neff = self.sw ** 2 / np.maximum((self.w ** 2).sum(1), 1e-30)
        self.cen = cen
    def mean(self, y):
        m = self.w @ y / np.maximum(self.sw, 1e-30)
        var = (self.w * (y[None, :] - m[:, None]) ** 2).sum(1) / np.maximum(self.sw, 1e-30)
        sem = np.sqrt(var / np.maximum(self.neff, 1))
        return np.where(self.ok, m, np.nan), np.where(self.ok, sem, np.nan)
    def near(self, radec):
        v = unit(*radec); return int(np.argmax(self.cen @ v))

def peak_stats(sm, y, label, nshuf=0):
    m, sem = sm.mean(y); g = np.nanmean(y)
    z = (m - g) / sem
    p1 = sm.near(WALL1); p2 = sm.near(WALL2)
    # the Wall 1 "peak" is the largest z within 15 deg of the nominal centre
    v = unit(*WALL1); near1 = (sm.cen @ v) > np.cos(np.radians(15)); near1 &= sm.ok
    i1 = np.where(near1)[0][np.nanargmax(z[near1])]
    zmax = np.nanmax(z); imax = int(np.nanargmax(z)); th, ph = hp.pix2ang(NSIDE, imax)
    print(f"{label}: N={len(y)} global mean {g:.0f}; Wall 1 pixel z={z[p1]:.2f} (excess {m[p1]-g:+.0f}), peak within 15 deg "
          f"z={z[i1]:.2f} (excess {m[i1]-g:+.0f}, n_eff {sm.neff[i1]:.0f}); Wall 2 z={z[p2]:.2f} ({m[p2]-g:+.0f}); "
          f"map max z={zmax:.2f} at RA {np.degrees(ph):.0f} Dec {90-np.degrees(th):.0f}; pixels {sm.ok.sum()}")
    if nshuf:
        zs1, zsmax = [], []
        for _ in range(nshuf):
            ys = rng.permutation(y); ms, ss = sm.mean(ys); zz = (ms - g) / ss
            zs1.append(np.nanmax(zz[near1])); zsmax.append(np.nanmax(zz))
        zs1, zsmax = np.array(zs1), np.array(zsmax)
        print(f"   DM-shuffle null ({nshuf}): P(peak within 15 deg of Wall 1 >= {z[i1]:.2f}) = {np.mean(zs1 >= z[i1]):.4f} "
              f"(null mean {zs1.mean():.2f}, sd {zs1.std():.2f}); P(map max >= {zmax:.2f}) = {np.mean(zsmax >= zmax):.4f} "
              f"(null mean {zsmax.mean():.2f}, sd {zsmax.std():.2f})")
    return m, sem, z

def detrend(d, y):
    """Residual of log DM_exc about a median regression on log S/N and declination (linear + quadratic), returned as a
    DM on the original scale (10**(resid + median))."""
    X = np.column_stack([np.ones(len(d)), np.log10(d.bonsai_snr.values), d.dec.values / 90, (d.dec.values / 90) ** 2])
    ly = np.log10(y); b = np.linalg.lstsq(X, ly, rcond=None)[0]
    for _ in range(100):
        r = ly - X @ b; w = 1 / np.maximum(np.abs(r), 1e-3); Xw = X * w[:, None]; b = np.linalg.solve(X.T @ Xw, Xw.T @ ly)
    print(f"   median regression of log10 DM_exc: logSNR {b[1]:+.3f}, dec/90 {b[2]:+.3f}, (dec/90)^2 {b[3]:+.3f}")
    return 10 ** (ly - X @ b + np.median(ly))

def main(chime, mrs, groups, nshuf=300):
    print("=== Part 1: the Wall 1 excess vs the catalogue's structure ===")
    d_all, nrep = load_chime(chime, False)
    d_one, _ = load_chime(chime, True)
    print(f"bursts after Ravi+ cuts: {len(d_all)} ({nrep} from repeaters); one per source: {len(d_one)}")
    sm_all = Smoother(d_all.ra.values, d_all.dec.values); sm_one = Smoother(d_one.ra.values, d_one.dec.values)
    for col in ['dm_exc_ne2001', 'dm_exc_ymw16']:
        print(f"-- {col}")
        peak_stats(sm_all, d_all[col].values, "all bursts (Ravi+-like)", nshuf)
        m1, s1, z1 = peak_stats(sm_one, d_one[col].values, "one burst per source", nshuf)
        # which repeaters sit inside the Wall 1 region, and how many bursts they contribute
        v = unit(*WALL1); inw = (unit(d_all.ra.values, d_all.dec.values) @ v) > np.cos(np.radians(15))
        rep_in = d_all[inw & d_all.repeater_name.notna()].groupby('repeater_name').agg(n=('tns_name', 'size'), dm=(col, 'median'), dec=('dec', 'first'))
        print(f"   within 15 deg of Wall 1: {inw.sum()} bursts, {d_all[inw].src.nunique()} sources; repeaters there: "
              + ', '.join(f"{k} (n={r.n}, DM_exc {r.dm:.0f}, Dec {r.dec:.0f})" for k, r in rep_in.sort_values('n', ascending=False).head(8).iterrows()))
        ydt = detrend(d_one, d_one[col].values)
        peak_stats(sm_one, ydt, "one per source, S/N and Dec trends removed", nshuf)
        # the S/N and Dec make-up of the Wall 1 region vs the rest
        w1 = d_one[(unit(d_one.ra.values, d_one.dec.values) @ v) > np.cos(np.radians(15))]
        rest = d_one.drop(w1.index)
        print(f"   Wall 1 region (one per source): n={len(w1)}, median DM_exc {w1[col].median():.0f} vs rest {rest[col].median():.0f}; "
              f"median S/N {w1.bonsai_snr.median():.1f} vs {rest.bonsai_snr.median():.1f}; median Dec {w1.dec.median():.0f} vs {rest.dec.median():.0f}; "
              f"mean DM_exc {w1[col].mean():.0f} vs {rest[col].mean():.0f}; 90th pct {w1[col].quantile(.9):.0f} vs {rest[col].quantile(.9):.0f}; "
              f"10th pct {w1[col].quantile(.1):.0f} vs {rest[col].quantile(.1):.0f}")
        # the same-declination control: bursts at the Wall's declinations (45-65) but other RAs
        ctrl = rest[(rest.dec > 45) & (rest.dec < 65)]
        print(f"   same-Dec control (Dec 45-65, outside): n={len(ctrl)}, median {ctrl[col].median():.0f}, mean {ctrl[col].mean():.0f}")
        # median instead of mean, one per source
        v = unit(*WALL1)
        print(f"   Wall-1 region median-based excess (one per source): {w1[col].median()-rest[col].median():+.0f} "
              f"(bootstrap sd {np.std([np.median(rng.choice(w1[col].values, len(w1))) for _ in range(500)]):.0f}); "
              f"mean-based {w1[col].mean()-rest[col].mean():+.0f} (sd {w1[col].std()/np.sqrt(len(w1)):.0f})")

    print("\n=== Part 2: the DM field vs the 2MRS luminosity density ===")
    g = vizier_tsv(mrs)
    for c in ['RAJ2000', 'DEJ2000', 'Kcmag', 'cz', 'GLAT']: g[c] = pd.to_numeric(g[c], errors='coerce')
    g = g[g.cz.notna() & (g.cz > 300) & g.Kcmag.notna()].copy()
    H0 = 70.0; g['dist'] = g.cz / H0
    g['LK'] = 10 ** (-0.4 * (g.Kcmag - 5 * np.log10(g.dist * 1e6 / 10) - 3.28))  # L_K in solar units (M_K,sun = 3.28)
    d = d_one; sm = sm_one
    y = d.dm_exc_ne2001.values; ydt = detrend(d, y)
    mdm, sdm = sm.mean(y); mdt, _ = sm.mean(ydt)
    shells = [(300, 3000), (3000, 6000), (6000, 9000), (9000, 13000), (13000, 20000), (300, 13000)]
    print("shell (km/s): N_gal, luminosity-density map correlation with the smoothed DM map over the unmasked pixels "
          "(raw DM; detrended DM); DM-shuffle null sd and p for the raw; per-burst slope of DM_exc on the shell's density "
          "contrast at the burst (pc/cm3 per unit delta, bootstrap error, shuffle p)")
    for lo, hi in shells:
        s = g[(g.cz >= lo) & (g.cz < hi)]
        # luminosity-weighted smoothed density over the same kernel, per unit solid angle; contrast about the mean over pixels
        smg = Smoother(s.RAJ2000.values, s.DEJ2000.values)
        dens = smg.w @ s.LK.values  # luminosity within the kernel (not normalised: the kernel is the same everywhere)
        # the 2MRS footprint: |b| > 5 only; CHIME bursts at |b| > 10, so the pixels compared have full 2MRS coverage
        delta = dens / dens[sm.ok].mean() - 1
        r_raw = np.corrcoef(delta[sm.ok], mdm[sm.ok])[0, 1]; r_dt = np.corrcoef(delta[sm.ok], mdt[sm.ok])[0, 1]
        rs = []
        for _ in range(nshuf):
            ms, _ = sm.mean(rng.permutation(y)); rs.append(np.corrcoef(delta[sm.ok], ms[sm.ok])[0, 1])
        rs = np.array(rs)
        # per-burst: density contrast at each burst (kernel-smoothed luminosity at the burst position)
        wb = np.exp(-np.arccos(np.clip(unit(d.ra.values, d.dec.values) @ unit(s.RAJ2000.values, s.DEJ2000.values).T, -1, 1)) ** 2 / (2 * SIG ** 2))
        lb = wb @ s.LK.values; db = lb / lb.mean() - 1
        X = np.column_stack([np.ones(len(d)), db]); b = np.linalg.lstsq(X, y, rcond=None)[0]
        bb = [np.linalg.lstsq(X[i], y[i], rcond=None)[0][1] for i in (rng.integers(0, len(y), len(y)) for _ in range(300))]
        bs = [np.linalg.lstsq(X, rng.permutation(y), rcond=None)[0][1] for _ in range(nshuf)]
        print(f"  {lo:5d}-{hi:5d}: N={len(s):5d}  r = {r_raw:+.3f} ({r_dt:+.3f}); null sd {rs.std():.3f}, p = {np.mean(rs >= r_raw):.3f};  "
              f"per-burst slope {b[1]:+.1f} ± {np.std(bb):.1f}, shuffle p = {np.mean(np.abs(bs) >= abs(b[1])):.3f}; "
              f"delta at Wall 1 pixel {delta[sm.near(WALL1)]:+.2f}, rank {int((delta[sm.ok] > delta[sm.near(WALL1)]).sum())+1} of {sm.ok.sum()}")
    # stacks around Tully 2015 groups
    t = vizier_tsv(groups)
    for c in ['_RA.icrs', '_DE.icrs', '<Vcmba>', 'Mlum', 'Nmb', 'logLK']: t[c] = pd.to_numeric(t[c], errors='coerce')
    t = t[t['<Vcmba>'].between(1000, 13000) & t.Nmb.ge(5) & t.Mlum.notna()]
    print(f"Tully 2015 groups with >= 5 members at 1000-13000 km/s: {len(t)}; stacked mean DM_exc (one per source) of bursts "
          f"within R of any group with M_lum above the cut, minus the rest, with a DM-shuffle null")
    ub = unit(d.ra.values, d.dec.values)
    for mcut in [1e13, 1e14]:
        tt = t[t.Mlum >= mcut]; ug = unit(tt['_RA.icrs'].values, tt['_DE.icrs'].values)
        cosmin = (ub @ ug.T).max(1)
        for R in [2, 5]:
            inn = cosmin > np.cos(np.radians(R))
            if inn.sum() < 20: print(f"  M >= {mcut:.0e}, R < {R} deg: {inn.sum()} bursts (too few)"); continue
            dlt = y[inn].mean() - y[~inn].mean()
            nulls = [(lambda ys: ys[inn].mean() - ys[~inn].mean())(rng.permutation(y)) for _ in range(nshuf)]
            print(f"  M >= {mcut:.0e} ({len(tt)} groups), R < {R} deg: {inn.sum()} bursts, excess {dlt:+.0f} pc/cm3, null sd {np.std(nulls):.0f}, "
                  f"p = {np.mean(np.array(nulls) >= dlt):.3f}; median excess {np.median(y[inn]) - np.median(y[~inn]):+.0f}")

if __name__ == '__main__':
    a = sys.argv[1:]; main(a[0], a[1], a[2], int(a[3]) if len(a) > 3 else 300)
