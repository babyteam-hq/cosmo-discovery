"""eRASS1 clusters (Bulbul+2024, VizieR J/A+A/685/A106 emain, 12,247 clusters): does the redshift distribution shift with
Galactic latitude or with E(B-V) at fixed X-ray flux (sweep-2.md row 26: the optical confirmation depth, LS DR10, leaking
into the sample), and does the photometric-vs-spectroscopic z fraction? Also the count per unit exposure-weighted area vs
|b| (the simplest selection check) and the dipole of the counts projected on the CMB axis (row 25, expected null).
usage: python erass1_nz.py data/erass1/emain.tsv data/dustmaps
"""
import sys, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from astropy.coordinates import SkyCoord
import astropy.units as u

def main(tsv, data_dir):
    d = pd.read_csv(tsv, sep='\t', comment='#', low_memory=False); d = d[pd.to_numeric(d['zBest'], errors='coerce').notna()].copy()
    for c in ('RAJ2000', 'DEJ2000', 'zBest', 'e_zBest', 'Exp', 'F500', 'pcont', 'ExtLike', 'DetLike0', 'KT'): d[c] = pd.to_numeric(d[c], errors='coerce')
    print(f"eRASS1 clusters: {len(d)}; zType values: {d.zType.value_counts().to_dict()}")
    c = SkyCoord(d.RAJ2000.values*u.deg, d.DEJ2000.values*u.deg, frame='icrs'); b = c.galactic.b.deg
    from dustmaps.config import config; config['data_dir'] = data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.planck import PlanckGNILCQuery
    s = SFDQuery()(c); g = PlanckGNILCQuery()(c)
    z = d.zBest.values; f = d.F500.values; ex = d.Exp.values; spec = d.zType.astype(str).str.contains('spec', case=False).values
    print(f"z median {np.median(z):.3f}; spec fraction {spec.mean():.2f}; F500 median {np.nanmedian(f):.2e}; exposure median {np.nanmedian(ex):.0f} s")
    # by |b| at fixed flux: clusters above the flux where the sample is near-complete (upper half)
    fcut = np.nanmedian(f); hi = f > fcut
    print(f"mean z and spec fraction by |b| (F500 > median, n={hi.sum()}):")
    for lo, up in [(10, 20), (20, 30), (30, 45), (45, 60), (60, 90)]:
        m = hi & (np.abs(b) >= lo) & (np.abs(b) < up)
        if m.sum() < 30: continue
        print(f"  {lo:2d}<=|b|<{up:2d} n={m.sum():5d} mean z {z[m].mean():.3f}±{z[m].std()/np.sqrt(m.sum()):.3f} median {np.median(z[m]):.3f} spec frac {spec[m].mean():.2f} mean E(B-V) {s[m].mean():.3f} mean exposure {np.nanmean(ex[m]):.0f}")
    r, p = spearmanr(np.abs(b[hi]), z[hi]); print(f"  Spearman(|b|, z) at F500>median: {r:+.3f} p={p:.2g}")
    print("mean z by E(B-V) (F500 > median):")
    for lo, up in [(0, 0.02), (0.02, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 2)]:
        m = hi & (s >= lo) & (s < up)
        if m.sum() < 30: continue
        print(f"  {lo:.2f}<=E(B-V)<{up:.2f} n={m.sum():5d} mean z {z[m].mean():.3f}±{z[m].std()/np.sqrt(m.sum()):.3f} spec frac {spec[m].mean():.2f} mean |b| {np.abs(b[m]).mean():.0f}")
    r, p = spearmanr(s[hi], z[hi]); print(f"  Spearman(E(B-V), z): {r:+.3f} p={p:.2g}")
    # at fixed |b| (>30), does z depend on GNILC - SFD (a map-error proxy for the optical depth of confirmation)?
    m = hi & (np.abs(b) > 30); r, p = spearmanr((g-s)[m], z[m]); print(f"  Spearman(GNILC-SFD, z) at |b|>30: {r:+.3f} p={p:.2g} (n={m.sum()})")
    # multiple regression: z on log exposure, E(B-V), |b|, log flux (all clusters)
    ok = np.isfinite(z) & np.isfinite(f) & (f > 0) & np.isfinite(ex) & (ex > 0)
    X = np.column_stack([np.ones(ok.sum()), np.log10(ex[ok]), s[ok], np.abs(b[ok])/90, np.log10(f[ok])]); bb, *_ = np.linalg.lstsq(X, z[ok], rcond=None)
    r = z[ok] - X@bb; cov = np.linalg.inv(X.T@X)*r.var(); e = np.sqrt(np.diag(cov))
    print(f"z regressed on [1, log exp, E(B-V), |b|/90, log F500]: coefficients {np.round(bb, 3)} ± {np.round(e, 3)} (n={ok.sum()})")
    # counts: number per exposure-weighted area in |b| bins is not possible without the exposure map; instead the mean exposure of
    # detected clusters vs |b|, and the count dipole along the CMB axis from the detected positions (western Galactic hemisphere only)
    n_eq = np.stack([np.cos(np.radians(d.DEJ2000.values))*np.cos(np.radians(d.RAJ2000.values)), np.cos(np.radians(d.DEJ2000.values))*np.sin(np.radians(d.RAJ2000.values)), np.sin(np.radians(d.DEJ2000.values))], -1)
    cmb = np.array([np.cos(np.radians(-6.944))*np.cos(np.radians(167.942)), np.cos(np.radians(-6.944))*np.sin(np.radians(167.942)), np.sin(np.radians(-6.944))])
    proj = n_eq@cmb; print(f"mean cos(angle to CMB apex) of the clusters: {proj.mean():+.4f} ± {proj.std()/np.sqrt(len(proj)):.4f} (the footprint's own mean dominates; a kinematic dipole of 0.01 would add ~0.003)")

main(*sys.argv[1:3])
