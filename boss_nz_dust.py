"""Sweep-3 row 5: does the BOSS CMASS/LOWZ redshift distribution (SDSS photometry, independent of the Legacy Survey) shift with
the DESI stellar-reddening correction to SFD (Zhou+2024, EBV_GR - EBV_SFD at the galaxy's nside-256 pixel), as DESI's own LRGs
do (R-116 row 4: -0.041 +- 0.011 in mean z per unit E(B-V)), and not with Planck-GNILC - SFD? Weighted regression of z on
[|b|/90, SFD, Zhou - SFD, GNILC - SFD, cap], quintiles, and a Galactic-longitude-rotation null for the sky-correlated error.
usage: python boss_nz_dust.py"""
import numpy as np, pandas as pd, healpy as hp
from astropy.io import fits
from cf4_io import sfd_ebv, gnilc_ebv, zhou_ebv, R_GAL, unit_eq, gal_lb
rng = np.random.default_rng(5)

def load(name):
    parts = []
    for cap in ['North', 'South']:
        t = fits.getdata(f'data/boss/galaxy_DR12v5_{name}_{cap}.fits.gz', 1)
        w = t['WEIGHT_SYSTOT'] * (t['WEIGHT_CP'] + t['WEIGHT_NOZ'] - 1)
        parts.append(pd.DataFrame(dict(ra=t['RA'], dec=t['DEC'], z=t['Z'], w=w, cap=cap == 'South')))
    return pd.concat(parts, ignore_index=True)

def wls(y, X, w):
    Xw = X * w[:, None]; A = X.T @ Xw; b = np.linalg.solve(A, Xw.T @ y); r = y - X @ b
    cov = np.linalg.inv(A) * (np.sum(w * r * r) / (len(y) - X.shape[1])); return b, np.sqrt(np.diag(cov))

def run(name, zlo, zhi):
    d = load(name); d = d[(d.z > zlo) & (d.z < zhi)].copy()
    l, b = gal_lb(unit_eq(d.ra.values, d.dec.values)); d['l'], d['b'] = l, b; d['absb'] = np.abs(b)
    d['sfd'] = sfd_ebv(l, b); d['gnilc'] = gnilc_ebv(l, b); zh, zs = zhou_ebv(l, b); d['zhou'] = zh; d['zsfd'] = zs
    d['g_s'] = d.gnilc - d.sfd; d['z_s'] = d.zhou - d.zsfd
    print(f"\n=== {name} {zlo} < z < {zhi}: {len(d)} galaxies, mean z {np.average(d.z, weights=d.w):.4f}; in the DESI-stars footprint "
          f"{d.zhou.notna().sum()} ({d.zhou.notna().mean()*100:.0f} %) ===")
    t = d[d.zhou.notna() & (d.absb > 20)].copy()
    print(f"|b| > 20 in footprint: {len(t)}; rms of Zhou - SFD {t.z_s.std():.4f}, GNILC - SFD {t.g_s.std():.4f}, corr {np.corrcoef(t.z_s, t.g_s)[0,1]:.2f}")
    X = np.column_stack([np.ones(len(t)), t.absb / 90, t.sfd, t.z_s, t.g_s, t.cap.astype(float)]); names = ['const', '|b|/90', 'SFD', 'Zhou-SFD', 'GNILC-SFD', 'SGC']
    bb, ee = wls(t.z.values, X, t.w.values)
    print("regression of z: " + ", ".join(f"{n} {b:+.4f}±{e:.4f}" for n, b, e in zip(names[1:], bb[1:], ee[1:])))
    # rotation null for the two map-difference coefficients: resample the maps at l rotated (b kept), 200 times
    zs, gs = [], []
    for _ in range(200):
        l2 = (t.l.values + rng.uniform(0, 360)) % 360; b2 = t.b.values * rng.choice([-1, 1])
        zh2, zs2 = zhou_ebv(l2, b2); ok = np.isfinite(zh2)
        if ok.sum() < 5000: continue  # the rotated sky misses the DESI-stars footprint
        X2 = np.column_stack([np.ones(ok.sum()), t.absb.values[ok] / 90, sfd_ebv(l2[ok], b2[ok]), (zh2 - zs2)[ok], (gnilc_ebv(l2[ok], b2[ok]) - sfd_ebv(l2[ok], b2[ok]))])  # no cap term: a rotation can leave one cap only
        b2v, e2v = wls(t.z.values[ok], X2, t.w.values[ok]); zs.append(b2v[3] / e2v[3]); gs.append(b2v[4] / e2v[4])
    print(f"rotation null: z-score spread Zhou-SFD {np.std(zs):.2f}, GNILC-SFD {np.std(gs):.2f}; data z-scores {bb[3]/ee[3]:+.2f}, {bb[4]/ee[4]:+.2f} "
          f"-> {bb[3]/ee[3]/np.std(zs):+.2f}σ, {bb[4]/ee[4]/np.std(gs):+.2f}σ equivalent")
    # quintiles of Zhou - SFD: mean z and the high-z fraction
    q = pd.qcut(t.z_s, 5, labels=False); zq = t.z.quantile(0.75)
    print("quintiles of Zhou-SFD (mean z, top-quartile-z fraction):", [(f"{np.average(t.z[q==k], weights=t.w[q==k]):.4f}",
          f"{np.average((t.z[q==k] > zq), weights=t.w[q==k]):.3f}") for k in range(5)], " (edges", np.round(t.z_s.quantile([0, .2, .4, .6, .8, 1]).values, 3), ")")
    q = pd.qcut(t.g_s, 5, labels=False)
    print("quintiles of GNILC-SFD (mean z):", [f"{np.average(t.z[q==k], weights=t.w[q==k]):.4f}" for k in range(5)])
    # the DESI-LRG expectation: -0.041 per unit E(B-V) of (Zhou - SFD) (R-116 row 4), given here for comparison
    print(f"DESI LRG (R-116): -0.041 ± 0.011 per unit E(B-V); here Zhou-SFD {bb[3]:+.4f} ± {ee[3]:.4f}")

if __name__ == '__main__':
    run('CMASS', 0.43, 0.70); run('LOWZ', 0.15, 0.43)
