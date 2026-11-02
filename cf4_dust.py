"""Cosmicflows-4 individual distances (Tully+2023, VizieR J/ApJ/944/94 table2) as an independent check of the dust-map
systematic found in SN residuals (finding.md): do TF, FP, SBF, SNIa distance-modulus residuals track (GNILC - SFD) with the
extinction-vector slope of the band? A map error dE under-corrects a galaxy's magnitude by R_band dE, so its distance
modulus is too large by R_band dE (no colour standardisation in TF/FP): the expected slope if GNILC were right is +R_band
(I/i band ~1.9-2.1 for the optical TF, W1 ~0.19 for the WISE TF, J/H/K ~0.9/0.6/0.35 for 6dFGS FP, r ~2.3 for SDSS FP).
Residual: DM - 5 log10(Vcmb/H0 [Mpc]) - 25 with the intercept free per indicator (so H0 and each zero point drop out),
Vcmb > 3000 km/s to keep peculiar velocities below the indicator scatter, |b| > 10. Diagonal errors (e_DM plus a 250 km/s
velocity term). Also (sweep-2.md rows 6, 9): zero-point offsets between Galactic and CMB-dipole hemispheres per indicator,
and the indicator differences for galaxies with two methods.
usage: python cf4_dust.py data/cf4/cf4_table2_full.tsv data/dustmaps
"""
import sys, numpy as np, pandas as pd
from scipy.stats import chi2 as chi2d
from astropy.coordinates import SkyCoord
import astropy.units as u
H0 = 74.6; C = 299792.458

def main(tsv, data_dir):
    d = pd.read_csv(tsv, sep='\t', comment='#', skiprows=[1, 2] if False else None, low_memory=False)
    d = d[pd.to_numeric(d['Vcmb'], errors='coerce').notna()].copy()
    for c in d.columns:
        if c not in ('PGC', 'CF3', 'Sloan', 'Simbad', 'NED', 'LEDA'): d[c] = pd.to_numeric(d[c], errors='coerce')
    print(f"CF4 table2: {len(d)} galaxies")
    c = SkyCoord(d.RAJ2000.values*u.deg, d.DEJ2000.values*u.deg, frame='icrs'); b = c.galactic.b.deg
    from dustmaps.config import config; config['data_dir'] = data_dir
    from dustmaps.sfd import SFDQuery; from dustmaps.planck import PlanckGNILCQuery
    s = SFDQuery()(c); g = PlanckGNILCQuery()(c); x = g - s; xa = g - 0.86*s
    n_eq = np.stack([np.cos(np.radians(d.DEJ2000.values))*np.cos(np.radians(d.RAJ2000.values)), np.cos(np.radians(d.DEJ2000.values))*np.sin(np.radians(d.RAJ2000.values)), np.sin(np.radians(d.DEJ2000.values))], -1)
    cmb = np.array([np.cos(np.radians(-6.944))*np.cos(np.radians(167.942)), np.cos(np.radians(-6.944))*np.sin(np.radians(167.942)), np.sin(np.radians(-6.944))]); sc = n_eq@cmb
    v = d.Vcmb.values; mu_h = 5*np.log10(np.clip(v, 1, None)/H0) + 25
    print("indicator: slope of the DM residual on (GNILC - SFD) and on (GNILC - 0.86 SFD), expected +R_band if GNILC is right; Vcmb>3000, |b|>10")
    for col, lbl, R in [('DMtf', 'TF (optical/WISE mix)', '1.9 (I) / 0.19 (W1)'), ('DMfp', 'FP (6dFGS JHK, SDSS r)', '0.35-0.9 (K-J) / 2.3 (r)'), ('DMsnIa', 'SN Ia', '~0.65 (after colour)'), ('DMsbf', 'SBF', '~1.9 (I)'), ('DM', 'CF4 combined DM', 'mix')]:
        e = d['e_'+col].values if 'e_'+col in d else d['e_DM'].values
        ok = np.isfinite(d[col].values) & (v > 3000) & (np.abs(b) > 10) & np.isfinite(e) & (e > 0)
        if ok.sum() < 100: print(f"  {lbl}: only {ok.sum()}"); continue
        y = d[col].values[ok] - mu_h[ok]; w = 1/np.sqrt(e[ok]**2 + (5/np.log(10)*250/v[ok])**2)
        def fit(cols):
            X = np.column_stack([np.ones(ok.sum())] + cols)*w[:, None]; yw = y*w; bb, *_ = np.linalg.lstsq(X, yw, rcond=None); r = yw - X@bb
            return bb, np.sqrt(np.diag(np.linalg.inv(X.T@X))), r@r
        _, _, c0 = fit([]); chi_dof = c0/(ok.sum()-1)
        # rescale errors so chi2/dof = 1 (the table errors are nominal)
        w = w/np.sqrt(chi_dof); _, _, c0 = fit([])
        b1, e1, c1 = fit([x[ok]]); b2, e2, c2 = fit([xa[ok]]); b3, e3, c3 = fit([s[ok]])
        print(f"  {lbl:24s} n={ok.sum():5d} rms resid {np.std(y):.3f} (chi2/dof before rescale {chi_dof:.2f}) rms(GNILC-SFD) {np.std(x[ok]):.4f} | GNILC-SFD {b1[1]:6.2f}±{e1[1]:4.2f} (dchi2 {c0-c1:5.1f}) | GNILC-0.86SFD {b2[1]:6.2f}±{e2[1]:4.2f} | SFD itself {b3[1]:6.2f}±{e3[1]:4.2f} | expected {R}")
        for cutlbl, cut in [("|b|>30", np.abs(b[ok]) > 30), ("|b|<30", np.abs(b[ok]) < 30), ("Dec>0", d.DEJ2000.values[ok] > 0), ("Dec<0", d.DEJ2000.values[ok] < 0), ("v<6000", v[ok] < 6000), ("v>6000", v[ok] > 6000)]:
            if cut.sum() < 100: continue
            X = np.column_stack([np.ones(cut.sum()), x[ok][cut]])*w[cut][:, None]; yw = y[cut]*w[cut]; bb, *_ = np.linalg.lstsq(X, yw, rcond=None); ee = np.sqrt(np.diag(np.linalg.inv(X.T@X)))
            print(f"      {cutlbl:7s} n={cut.sum():5d} GNILC-SFD {bb[1]:6.2f}±{ee[1]:4.2f}")
        # hemispheres (rows 6, 9): zero-point difference north-south Galactic, and along the CMB dipole axis (a dipole in the residual)
        bg = b[ok]; bb, ee, _ = fit([(bg > 0).astype(float)]); bc, ec, _ = fit([sc[ok]])
        print(f"      Galactic N minus S zero point {bb[1]*1e3:+.1f}±{ee[1]*1e3:.1f} mmag; residual dipole along the CMB axis {bc[1]*1e3:+.1f}±{ec[1]*1e3:.1f} mmag (a bulk flow of 300 km/s at v=5000 would give ~{5/np.log(10)*300/5000*1e3:.0f})")
    # galaxies with two indicators: TF - FP, TF - SN, SBF - TF by hemisphere
    print("indicator differences for galaxies holding two methods (zero-point consistency across the sky):")
    for a_, b_ in [('DMtf', 'DMfp'), ('DMtf', 'DMsnIa'), ('DMfp', 'DMsnIa'), ('DMsbf', 'DMtf')]:
        ok = np.isfinite(d[a_].values) & np.isfinite(d[b_].values)
        if ok.sum() < 20: print(f"  {a_}-{b_}: {ok.sum()} galaxies"); continue
        diff = d[a_].values[ok] - d[b_].values[ok]; bg = b[ok]
        n_, s_ = diff[bg > 0], diff[bg < 0]
        print(f"  {a_}-{b_}: n={ok.sum()} mean {diff.mean():+.3f}±{diff.std()/np.sqrt(len(diff)):.3f}; Galactic N {n_.mean():+.3f}±{n_.std()/np.sqrt(max(len(n_),1)):.3f} (n={len(n_)}) vs S {s_.mean():+.3f}±{s_.std()/np.sqrt(max(len(s_),1)):.3f} (n={len(s_)}); slope on GNILC-SFD {np.polyfit(x[ok], diff, 1)[0]:+.2f}")

main(*sys.argv[1:3])
