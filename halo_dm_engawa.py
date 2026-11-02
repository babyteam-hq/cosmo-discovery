"""Does the floor of the CHIME/FRB Catalog 2 extragalactic DM follow the anisotropic Milky Way halo DM map of the ENGAWA
simulations (Lucchini+2026, arXiv:2607.20601; l26halodm maps: nside-128 HEALPix in Galactic coordinates, 16/50/84th
percentiles, median 10 pc/cm3 at the poles rising to 44 at |b| 10-20, with a 13-38 pc/cm3 spread across pixels at |b| 30-60),
or a smooth model (constant, or constant + csc|b|)? Sweep-4 row E (the open part of the trajectory verdict: the map has never
been confronted with FRB data). One burst per source, |b| > 20. (1) Sky cells (nside 4, >= 25 bursts): the 10th and 20th
percentile and the median of DM_exc per cell, bootstrap errors, regressed on the cell-mean ENGAWA median with and without a
csc|b| term; the longitude structure alone is tested with the map rotated in l (which keeps its |b| profile) as the null.
(2) Quantile regression on individual bursts: q-th quantile of DM_exc = a + b ENGAWA + c csc|b| + declination terms, bootstrap
errors; b = 1 if the map's anisotropy is in the data, 0 if the floor is set by something else.
usage: python -I halo_dm_engawa.py data/chime/chimefrbcat2.csv data/engawa/L26_allsky_sim.npy"""
import sys, numpy as np, pandas as pd, healpy as hp
rng = np.random.default_rng(4)

def chime(path):
    d = pd.read_csv(path)
    d = d[(d.excluded_flag == 0) & (d.sub_num == 0) & d.dec.notna() & (d.bonsai_snr >= 9)].copy()
    rep = d.repeater_name.notna(); first = d[rep].sort_values('mjd_400').drop_duplicates('repeater_name')
    d = pd.concat([d[~rep], first]); d = d[(d.dm_exc_ne2001 > 0) & (d.dm_exc_ymw16 > 0)].reset_index(drop=True)
    d['absb'] = np.abs(d.gb); return d

def quantile_reg(y, X, q, iters=300):
    """q-th quantile regression by IRLS on the check function."""
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    for _ in range(iters):
        r = y - X @ b; w = np.where(r > 0, q, 1 - q) / np.maximum(np.abs(r), 1e-3)
        Xw = X * w[:, None]; b = np.linalg.solve(X.T @ Xw + 1e-9 * np.eye(X.shape[1]), Xw.T @ y)
    return b

def wls(y, X, w):
    Xw = X * w[:, None]; cov = np.linalg.inv(X.T @ Xw); b = cov @ (Xw.T @ y); return b, np.sqrt(np.diag(cov))

def main(cat, mapf):
    d = chime(cat); maps = np.load(mapf, allow_pickle=False); p16, p50, p84 = maps
    NS = hp.npix2nside(p50.shape[0])
    at = lambda m, l, b: hp.get_interp_val(m, np.radians(90 - b), np.radians(l))
    d['eng'] = at(p50, d.gl.values, d.gb.values); d['eng16'] = at(p16, d.gl.values, d.gb.values)
    d['csc'] = 1 / np.sin(np.radians(d.absb))
    print(f"CHIME Cat 2 one per source, S/N >= 9: {len(d)}; ENGAWA median at the bursts: mean {d.eng.mean():.1f}, 16-84 pct {np.percentile(d.eng, 16):.1f}-{np.percentile(d.eng, 84):.1f}")
    for bcut in [20, 30]:
        s = d[d.absb > bcut].copy()
        print(f"\n=== |b| > {bcut}: {len(s)} sources; ENGAWA median spread here {np.percentile(s.eng, 16):.1f}-{np.percentile(s.eng, 84):.1f} (sd {s.eng.std():.1f}); corr(ENGAWA, csc|b|) = {np.corrcoef(s.eng, s.csc)[0,1]:.2f} ===")
        for col in ['dm_exc_ne2001', 'dm_exc_ymw16']:
            y = s[col].values
            # (1) cells
            pix = hp.ang2pix(4, s.ra.values, s.dec.values, lonlat=True); cells = [p for p in np.unique(pix) if (pix == p).sum() >= 25]
            rows = []
            for p in cells:
                m = pix == p; yy = y[m]
                qs = [np.percentile(yy, q) for q in (10, 20, 50)]
                bs = np.array([[np.percentile(rng.choice(yy, len(yy)), q) for q in (10, 20, 50)] for _ in range(300)])
                rows.append([s.eng.values[m].mean(), s.csc.values[m].mean(), m.sum()] + qs + list(bs.std(0)))
            R = np.array(rows); print(f"-- {col}: {len(cells)} cells of >= 25 bursts; cell-mean ENGAWA range {R[:,0].min():.1f}-{R[:,0].max():.1f}")
            for k, name in enumerate(['10th pct', '20th pct', 'median']):
                yq, eq = R[:, 3 + k], R[:, 6 + k]; w = 1 / eq ** 2
                b1, e1 = wls(yq, np.column_stack([np.ones(len(R)), R[:, 0]]), w)
                b2, e2 = wls(yq, np.column_stack([np.ones(len(R)), R[:, 0], R[:, 1]]), w)
                b3, e3 = wls(yq, np.column_stack([np.ones(len(R)), R[:, 1]]), w)
                chi1 = (w * (yq - np.column_stack([np.ones(len(R)), R[:, 0]]) @ b1) ** 2).sum(); chi3 = (w * (yq - np.column_stack([np.ones(len(R)), R[:, 1]]) @ b3) ** 2).sum()
                # longitude-rotation null for the ENGAWA slope at fixed csc|b| term
                rot = []
                for _ in range(200):
                    dl = rng.uniform(0, 360); e_rot = at(p50, (s.gl.values + dl) % 360, s.gb.values)
                    er = np.array([e_rot[pix == p].mean() for p in cells])
                    rot.append(wls(yq, np.column_stack([np.ones(len(R)), er, R[:, 1]]), w)[0][1])
                rot = np.array(rot)
                print(f"   {name:8s}: on ENGAWA alone {b1[1]:+.2f} ± {e1[1]:.2f} (chi2 {chi1:.0f}/{len(R)-2}); with csc|b|: ENGAWA {b2[1]:+.2f} ± {e2[1]:.2f}, csc {b2[2]:+.1f} ± {e2[2]:.1f}; "
                      f"csc alone {b3[1]:+.1f} ± {e3[1]:.1f} (chi2 {chi3:.0f}); rotated-map null for the ENGAWA slope: {rot.mean():+.2f} ± {rot.std():.2f}, P(>= data) = {np.mean(rot >= b2[1]):.3f}")
            # (2) per-burst quantile regression with CHIME's declination terms
            X = np.column_stack([np.ones(len(s)), s.eng.values, s.csc.values, s.dec.values / 90, (s.dec.values / 90) ** 2])
            for q in (0.1, 0.2, 0.5):
                b = quantile_reg(y, X, q)
                bb = np.array([quantile_reg(y[i], X[i], q) for i in (rng.integers(0, len(y), len(y)) for _ in range(100))])
                Xr = X.copy(); rot = []
                for _ in range(100):
                    Xr[:, 1] = at(p50, (s.gl.values + rng.uniform(0, 360)) % 360, s.gb.values); rot.append(quantile_reg(y, Xr, q)[1])
                print(f"   per-burst q={q}: ENGAWA {b[1]:+.2f} ± {bb[:,1].std():.2f}, csc|b| {b[2]:+.1f} ± {bb[:,2].std():.1f}; rotated-map null {np.mean(rot):+.2f} ± {np.std(rot):.2f}, P(>= data) = {np.mean(np.array(rot) >= b[1]):.2f}")
        # what the floor looks like against |b| in bands, for the record
        print("   DM_exc (NE2001) 10th percentile by |b| band vs ENGAWA median there:")
        for lo, hi in [(20, 30), (30, 40), (40, 50), (50, 60), (60, 90)]:
            m = (s.absb > lo) & (s.absb <= hi)
            if m.sum() > 30: print(f"     |b| {lo}-{hi}: n={m.sum()}, q10 {np.percentile(s.dm_exc_ne2001[m], 10):.0f}, q20 {np.percentile(s.dm_exc_ne2001[m], 20):.0f}, ENGAWA {s.eng[m].mean():.1f}")

if __name__ == '__main__':
    main(*sys.argv[1:])
