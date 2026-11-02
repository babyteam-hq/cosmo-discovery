"""The NGC 1052 trail of dark-matter-free dwarfs (DF2, DF4, DF9 ...; van Dokkum+2022 bullet-dwarf scenario): tests the
published positions, velocities and distances allow that the papers did not report (sweep-2.md rows 31, 32, 35, 36, 33).
Numbers: positions van Dokkum+2022 Table 1 (Nature 605, 435); velocities and (dx, dy) along the DF2-DF4 line Keim+2025
Table 2 (arXiv:2506.10220); SBF/TRGB distances Tang+2026 Table 1 (arXiv:2606.05144); group mean and dispersion
1447 / 123 km/s (Keim+2025, including LEDA 4014647 and DF7) and 1435 / 115 (van Dokkum+2022, 28 galaxies).
 (1) v(x): straight line with free slope, intercept and intrinsic scatter (profile likelihood) against a constant with
     free scatter, with and without the two objects the papers call interlopers (DF7, LEDA 4014647); F-test, delta-chi2,
     AIC; leave-one-out slopes.
 (2) membership: two-component mixture (group Gaussian vs trail line with its fitted scatter), per-object probabilities,
     with the trail prior = fraction of objects.
 (3) distance vs x with the Tang+2026 SBF distances: slope posterior against the bullet-dwarf prediction (-0.06 Mpc/arcmin
     in RA, i.e. ~ +0.06 per arcmin of x since x runs with RA) and against a constant.
 (4) chance alignment: Poisson expectation of >=k LSB dwarfs in a 69' x 10.4' strip at the Roman+2021 density (1.2/deg^2)
     and at the density the 12 trail objects themselves imply for the group core (a 1 deg radius).
 (5) DESI DR1 dark-program spectra in the trail field (redrock-main-dark-36785/36787.fits): group-velocity objects.
usage: python trail_ngc1052.py [data/desi_trail]
"""
import sys, numpy as np
from scipy.optimize import minimize_scalar, minimize
from scipy.stats import f as fdist, chi2 as chi2d, norm, poisson
# name, dx (arcmin along DF2-DF4 line, from its midpoint; +x toward DF2/RCP32, i.e. increasing RA), dy, v, +err, -err
VEL = [("RCP32", 29.8, 5.8, 1783, 4, 4), ("DF2", 20.8, 0.0, 1805, 1, 1), ("RCP26", 9.7, 1.9, 1670, 10, 10), ("DF9", -6.2, -0.8, 1648, 1, 1),
       ("Ta21-12000", -12.5, -3.6, 1540, 9, 8), ("DF5", -20.8, -1.6, 1440, 30, 20), ("DF4", -20.8, 0.0, 1433.3, 0.3, 0.4),
       ("LEDA4014647", -28.5, 0.9, 1596, 3, 2), ("DF7", -35.9, 5.7, 1690, 10, 10)]
INTERLOPERS = ("LEDA4014647", "DF7")
# Tang+2026 SBF distances (Mpc), with RA (deg) for the x axis: x_arcmin = (RA - RA_mid)*60*cos(dec)
SBF = [("DF2", 40.44531, 17.7, 1.4), ("RCP26", 40.28970, 18.9, 2.2), ("DF9", 40.02927, 18.7, 1.5), ("Ta21-12000", 39.91385, 21.8, 2.6),
       ("DF4", 39.81271, 18.8, 1.7), ("DF5", 39.80284, 20.5, 2.7), ("LEDA4014647", 39.70206, 17.7, 1.4), ("DF7", 39.62401, 17.8, 2.4)]
TRGB = {"DF2 HST": 21.7, "DF2 JWST": 17.6, "DF4 HST": 20.0}
GROUP = (1447.0, 123.0)

def fit_line(x, v, e, const=False):
    """Profile likelihood over intrinsic scatter s for v = a + b x (or v = a), returns (a, b, s, -2lnL, chi2 at s)."""
    def nll(s):
        w = 1/(e**2 + s**2); X = np.column_stack([np.ones_like(x)] + ([] if const else [x]))
        beta = np.linalg.solve((X*w[:, None]).T@X, (X*w[:, None]).T@v); r = v - X@beta
        return 0.5*np.sum(r**2*w) + 0.5*np.sum(np.log(2*np.pi/w)), beta, r, w
    res = minimize_scalar(lambda s: nll(s)[0], bounds=(0, 500), method='bounded'); s = res.x; n2, beta, r, w = nll(s)
    cov = np.linalg.inv((np.column_stack([np.ones_like(x)] + ([] if const else [x]))*w[:, None]).T@np.column_stack([np.ones_like(x)] + ([] if const else [x])))
    return beta, np.sqrt(np.diag(cov)), s, 2*n2, np.sum(r**2*w)

def part1():
    print("(1) velocity vs position along the trail")
    for label, excl in [("all 9 with velocities", ()), ("7 without DF7 and LEDA 4014647 (the papers' interlopers)", INTERLOPERS)]:
        rows = [r for r in VEL if r[0] not in excl]; x = np.array([r[1] for r in rows], float); v = np.array([r[3] for r in rows], float)
        e = np.array([max(r[4], r[5]) for r in rows], float)
        (a, b), (ea, eb), s, m2l, c2 = fit_line(x, v, e); (a0,), (ea0,), s0, m2l0, c20 = fit_line(x, v, e, const=True)
        n = len(x); d = m2l0 - m2l
        # F-test with the scatter fixed at the constant-model value (the classical nested-model test)
        w0 = 1/(e**2 + s0**2); X = np.column_stack([np.ones(n), x]); beta = np.linalg.solve((X*w0[:, None]).T@X, (X*w0[:, None]).T@v)
        rss1 = np.sum((v - X@beta)**2*w0); rss0 = np.sum((v - np.sum(v*w0)/np.sum(w0))**2*w0); F = (rss0 - rss1)/(rss1/(n-2)); pF = fdist.sf(F, 1, n-2)
        print(f"  {label}: n={n}")
        print(f"    line: v = {a:.0f}±{ea:.0f} + ({b:.2f}±{eb:.2f}) x km/s per arcmin, intrinsic scatter {s:.0f} km/s; constant: {a0:.0f}±{ea0:.0f}, scatter {s0:.0f}")
        print(f"    -2 delta lnL = {d:.1f} (line vs constant, 1 extra parameter; p={chi2d.sf(d,1):.2g}); AIC favours the line by {d-2:.1f}; F({1},{n-2})={F:.1f}, p={pF:.2g}; slope {b/eb:.1f} sigma from 0")
        print(f"    bullet-dwarf prediction (van Dokkum+2022): ~ -10 km/s per arcmin of RA, i.e. slope ~ +9 per arcmin of x here (DF2-DF4: {(1805-1433.3)/41.6:.1f})")
        print("    leave-one-out slopes:", ", ".join(f"{rows[i][0]} out: {fit_line(np.delete(x,i), np.delete(v,i), np.delete(e,i))[0][1]:.1f}" for i in range(n)))
    # the sign test that needs no model: are the 7 velocities monotonic in x?
    rows = [r for r in VEL if r[0] not in INTERLOPERS]; x = np.array([r[1] for r in rows]); v = np.array([r[3] for r in rows]); o = np.argsort(x)
    from scipy.stats import spearmanr, kendalltau
    print(f"    rank correlation (7): Spearman {spearmanr(x, v)[0]:.2f} p={spearmanr(x, v)[1]:.3g}; Kendall tau {kendalltau(x, v)[0]:.2f} p={kendalltau(x, v)[1]:.3g}")
    x9 = np.array([r[1] for r in VEL]); v9 = np.array([r[3] for r in VEL])
    print(f"    rank correlation (9): Spearman {spearmanr(x9, v9)[0]:.2f} p={spearmanr(x9, v9)[1]:.3g}; Kendall tau {kendalltau(x9, v9)[0]:.2f} p={kendalltau(x9, v9)[1]:.3g}")

def part2():
    print("(2) membership: group N(1447, 123) vs trail line (fit to the 7, with its scatter), prior = fraction")
    rows7 = [r for r in VEL if r[0] not in INTERLOPERS]; x7 = np.array([r[1] for r in rows7], float); v7 = np.array([r[3] for r in rows7], float); e7 = np.array([max(r[4], r[5]) for r in rows7], float)
    (a, b), _, s, _, _ = fit_line(x7, v7, e7)
    for prior in (0.5, 7/9):
        print(f"  trail prior {prior:.2f}:")
        for name, x, dy, v, ep, em in VEL:
            e = max(ep, em); lt = norm.pdf(v, a + b*x, np.sqrt(s**2 + e**2)); lg = norm.pdf(v, GROUP[0], np.sqrt(GROUP[1]**2 + e**2))
            p = prior*lt/(prior*lt + (1-prior)*lg)
            print(f"    {name:12s} x={x:6.1f} v={v:6.0f}  trail predicts {a+b*x:6.0f}±{s:.0f}  group {GROUP[0]:.0f}±{GROUP[1]:.0f}  P(trail)={p:.3f}  ({(v-GROUP[0])/GROUP[1]:+.1f} sigma from the group mean)")
    # how much of the 7's "fit" is just that the group's velocity range is 2 sigma = 250 km/s wide: a constant-velocity
    # group member set with sigma 123 would give rms about any line of ~123; the observed is s
    print(f"  intrinsic scatter about the line (7): {s:.0f} km/s vs the group's 123: ratio {s/123:.2f}; chance that 7 draws from N(1447,123) have rms <= {s:.0f} about their best line:", end=" ")
    rng = np.random.default_rng(1); cnt = 0; N = 20000
    for _ in range(N):
        vv = rng.normal(GROUP[0], GROUP[1], 7); X = np.column_stack([np.ones(7), x7]); beta = np.linalg.lstsq(X, vv, rcond=None)[0]; r = vv - X@beta
        if np.sqrt(np.mean(r**2)) <= s: cnt += 1
    print(f"{cnt/N:.3f} (positions fixed, any slope allowed)")
    # the selection the papers make: of 9 objects, the 2 farthest from the best line are called interlopers. Null: 9 draws
    # from the group Gaussian at the 9 positions; drop the 2 most deviant from the 9-point line; refit the 7; how often is
    # the 7-point line as significant as observed (-2 dlnL >= observed) and the rms as small?
    x9 = np.array([r[1] for r in VEL], float); e9 = np.array([max(r[4], r[5]) for r in VEL], float); v9 = np.array([r[3] for r in VEL], float)
    def best7(vv):
        X = np.column_stack([np.ones(9), x9]); beta = np.linalg.lstsq(X, vv, rcond=None)[0]; r = np.abs(vv - X@beta); keep = np.argsort(r)[:7]
        return x9[keep], vv[keep], e9[keep]
    xo, vo, eo = best7(v9); (_, bo), _, so, m2o, _ = fit_line(xo, vo, eo); (_,), _, _, m2o0, _ = fit_line(xo, vo, eo, const=True); dobs = m2o0 - m2o
    print(f"  'best 7 of 9' null (drop the 2 most deviant from the 9-point line, refit): observed -2dlnL={dobs:.1f}, slope {bo:.1f}, scatter {so:.0f};", end=" ")
    N = 4000; c1 = c2 = 0
    for _ in range(N):
        xo, vo, eo = best7(rng.normal(GROUP[0], GROUP[1], 9)); (_, bb), _, ss, m2, _ = fit_line(xo, vo, eo); (_,), _, _, m20, _ = fit_line(xo, vo, eo, const=True)
        if m20 - m2 >= dobs: c1 += 1
        if m20 - m2 >= dobs and ss <= so: c2 += 1
    print(f"P(-2dlnL >= obs) = {c1/N:.4f}; and scatter <= {so:.0f}: {c2/N:.4f}  ({N} draws)")
    # the papers' own procedure (Keim+2025): the line is defined by DF2 and DF4; of the other 7, keep the 5 closest to it.
    names = [r[0] for r in VEL]; i2, i4 = names.index("DF2"), names.index("DF4"); others = [i for i in range(9) if i not in (i2, i4)]
    def keim(vv):
        sl = (vv[i2] - vv[i4])/(x9[i2] - x9[i4]); pred = vv[i4] + sl*(x9 - x9[i4]); off = np.abs(vv - pred)
        keep = sorted(others, key=lambda i: off[i])[:5] + [i2, i4]; keep = np.array(keep)
        return x9[keep], vv[keep], e9[keep]
    xo, vo, eo = keim(v9); (_, bo), _, so, m2o, _ = fit_line(xo, vo, eo); (_,), _, _, m2o0, _ = fit_line(xo, vo, eo, const=True); dobs = m2o0 - m2o
    print(f"  Keim-style null (line through DF2 and DF4, keep the 5 of 7 closest; 9 velocities drawn from the group Gaussian): observed -2dlnL={dobs:.1f}, scatter {so:.0f};", end=" ")
    c1 = c2 = 0
    for _ in range(N):
        xo, vo, eo = keim(rng.normal(GROUP[0], GROUP[1], 9)); (_, bb), _, ss, m2, _ = fit_line(xo, vo, eo); (_,), _, _, m20, _ = fit_line(xo, vo, eo, const=True)
        if m20 - m2 >= dobs: c1 += 1
        if ss <= so: c2 += 1
    print(f"P(-2dlnL >= obs) = {c1/N:.4f}; P(scatter <= {so:.0f}) = {c2/N:.4f}  ({N} draws; Keim+2025's uniform-velocity version gave 2%)")

def part3():
    print("(3) distance vs position (Tang+2026 SBF), prediction D ~ 21.5 - 0.06 dRA(arcmin) Mpc")
    ra0 = (40.44531 + 39.81271)/2; cosd = np.cos(np.radians(-8.26))
    x = np.array([(r[1]-ra0)*60*cosd for r in SBF]); D = np.array([r[2] for r in SBF]); e = np.array([r[3] for r in SBF])
    for label, sel in [("8 SBF dwarfs", np.ones(8, bool)), ("6 without LEDA 4014647 and DF7", np.array([r[0] not in INTERLOPERS for r in SBF]))]:
        (a, b), (ea, eb), s, m2l, _ = fit_line(x[sel], D[sel], e[sel]); (a0,), (ea0,), s0, m2l0, _ = fit_line(x[sel], D[sel], e[sel], const=True)
        # the predicted slope: D decreases with RA at 0.06 Mpc/arcmin, x increases with RA, so predicted b = -0.06
        print(f"  {label}: D = {a:.1f}±{ea:.1f} + ({b:+.3f}±{eb:.3f}) x Mpc/arcmin, scatter {s:.1f}; constant {a0:.1f}±{ea0:.1f} (scatter {s0:.1f}); -2dlnL={m2l0-m2l:.1f}; predicted slope -0.06: {(b+0.06)/eb:.1f} sigma away; zero: {b/eb:.1f} sigma")
    print(f"  TRGB: DF2 HST 21.7±1.2 vs JWST 17.6±0.6 ({(21.7-17.6)/np.hypot(1.2,0.6):.1f} sigma apart); DF4 HST 20.0±1.6; the bullet-dwarf DF2-DF4 separation 1.7±0.5 Mpc (Shen+2021) vs SBF DF2-DF4 {17.7-18.8:+.1f}±{np.hypot(1.4,1.7):.1f}")

def part4():
    print("(4) chance alignment of LSB dwarfs on a 69' x 10.4' strip")
    area = 69*10.4/3600
    for dens, lbl in [(1.2, "Roman+2021 mean LSBG density, 1.2/deg^2"), (12/(np.pi*1.0**2), "the 12 trail objects spread over a 1-deg-radius core"), (42/36, "Roman+2021's 42 LSBGs over 6x6 deg")]:
        mu = dens*area; print(f"  {lbl}: expected {mu:.2f} per strip; P(>=7)={poisson.sf(6, mu):.2g}, P(>=11)={poisson.sf(10, mu):.2g}")
    print("  (van Dokkum+2022's radius-preserving randomisation gave 3% for >=11, 0.6% with DF2 and DF4 both on the line; the Poisson numbers above are the no-clustering baseline, not a replacement)")

def part5(desi_dir):
    print("(5) DESI DR1 dark-program spectra in the trail field (RA 39.5-40.8, Dec -8.5 to -7.85)")
    try:
        from astropy.io import fits
    except ImportError: print("  no astropy"); return
    for p in (36785, 36787):
        try: h = fits.open(f"{desi_dir}/redrock-main-dark-{p}.fits")
        except Exception as ex: print("  no file:", ex); return
        d = h['REDSHIFTS'].data; f = h['FIBERMAP'].data; ra, dec, z, zw, sp = f['TARGET_RA'], f['TARGET_DEC'], d['Z'], d['ZWARN'], d['SPECTYPE']
        box = (ra > 39.5) & (ra < 40.8) & (dec > -8.5) & (dec < -7.85); good = box & (zw == 0) & (sp == 'GALAXY') & (z*299792.458 > 900) & (z*299792.458 < 2200)
        print(f"  pixel {p}: {box.sum()} spectra in the box; group-velocity galaxies (900-2200 km/s, ZWARN 0): {good.sum()}")
        for i in np.where(good)[0]:
            # distance from the trail line (van Dokkum+2022: dDEC = 0.45 dRA - 2.8', through DF2/DF4 midpoint) and nearest trail object
            ra0, dec0 = (40.44531 + 39.81271)/2, (-8.40297 - 8.11597)/2; cosd = np.cos(np.radians(dec0))
            dx = (ra[i]-ra0)*60*cosd; dyl = (dec[i]-dec0)*60 - 0.45*dx
            print(f"    RA {ra[i]:.4f} Dec {dec[i]:.4f} cz={z[i]*299792.458:.0f} km/s  morph {f['MORPHTYPE'][i]}  {dyl:+.1f}' off the trail line, x={dx:+.1f}'")
    print("  NGC 1052 itself (40.2700, -8.2558) is the cz~1391 object in pixel 36785 (NED 1488; the DESI fibre is on the LINER nucleus)")

part1(); part2(); part3(); part4(); part5(sys.argv[1] if len(sys.argv) > 1 else "data/desi_trail")
