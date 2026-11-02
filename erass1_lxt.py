"""eRASS1 L_X-T anisotropy (the Migkas+2020/21 "H0 dipole" test) and whether it is the survey's exposure pattern.
Data: eRASS1 primary cluster catalogue (Bulbul+2024, erass1cl_primary_v3.2.fits): L500 (0.2-2.3 keV, erg/s), KT (keV),
BEST_Z, EXP (s), CTS500, RA, DEC. Fit log L500 - log E(z) = a + b log KT (OLS with weights from the T and L errors),
residual r = log L_obs - log L_fit. Then: (1) weighted mean residual in a 75-deg cone about Migkas's axis (l,b)=(280,-15)
and its shuffle-null significance; the full-sky scan's extreme cone (look-elsewhere); (2) residual vs log EXP; (3) the
same cone test with a log EXP term in the relation, and within exposure quartiles. L ∝ H0^-2, so dH0/H0 = -ln10/2 * dlogL.
usage: python erass1_lxt.py data/erass1/erass1cl_primary_v3.2.fits
"""
import sys, numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.stats import spearmanr
rng = np.random.default_rng(7); OM = 0.3
def Ez(z): return np.sqrt(OM*(1+z)**3 + 1-OM)

def main(f):
    t = fits.open(f)[1].data
    kt, ktl, kth = t['KT'], t['KT_L'], t['KT_H']; L, Ll, Lh = t['L500'], t['L500_L'], t['L500_H']; z = t['BEST_Z']; exp = t['EXP']
    ok = np.isfinite(kt) & (kt > 0) & np.isfinite(L) & (L > 0) & (z > 0.05) & (z < 0.6) & np.isfinite(exp) & (exp > 0)
    ekt = np.where(np.isfinite(ktl) & np.isfinite(kth), 0.5*(np.abs(kth - kt) + np.abs(kt - ktl)), np.nan)
    eL = np.where(np.isfinite(Ll) & np.isfinite(Lh), 0.5*(np.abs(Lh - L) + np.abs(L - Ll)), np.nan)
    ok &= np.isfinite(ekt) & np.isfinite(eL) & (ekt/kt < 0.5)
    idx = np.where(ok)[0]
    kt, ekt, L, eL, z, exp = kt[idx], ekt[idx], L[idx], eL[idx], z[idx], exp[idx]
    sc = SkyCoord(ra=t['RA'][idx]*u.deg, dec=t['DEC'][idx]*u.deg).galactic
    l, b = sc.l.deg, sc.b.deg
    n = np.stack([np.cos(np.radians(b))*np.cos(np.radians(l)), np.cos(np.radians(b))*np.sin(np.radians(l)), np.sin(np.radians(b))], -1)
    x = np.log10(kt/4.0); y = np.log10(L/1e44) - np.log10(Ez(z)); ey = eL/L/np.log(10); ex = ekt/kt/np.log(10)
    def fit(cols):
        X = np.column_stack([np.ones_like(x), x] + cols)
        w = 1/(ey**2 + (2.0*ex)**2 + 0.25**2)   # slope ~2 propagates T error; 0.25 dex intrinsic scatter
        beta = np.linalg.solve((X*w[:, None]).T@X, (X*w[:, None]).T@y); r = y - X@beta
        return beta, r, w
    beta, r, w = fit([])
    print(f"N={len(idx)} clusters with T (rel. err < 50%), 0.05<z<0.6. L-T: logL = {beta[0]:.3f} + {beta[1]:.2f} log(T/4keV); residual rms {np.sqrt(np.average(r**2, weights=w)):.3f} dex")
    def cone_mean(r, w, axis, ang=75):
        m = n@axis >= np.cos(np.radians(ang))
        if m.sum() < 60: return np.nan, m.sum()
        return np.average(r[m], weights=w[m]), m.sum()
    axis = np.array([np.cos(np.radians(-15))*np.cos(np.radians(280)), np.cos(np.radians(-15))*np.sin(np.radians(280)), np.sin(np.radians(-15))])
    def scan(r, w, ang=75, nside=600):
        # fibonacci sphere grid of directions; return the min and max cone means and their directions
        i = np.arange(nside) + 0.5; phi = np.arccos(1 - 2*i/nside); th = np.pi*(1 + 5**0.5)*i
        dirs = np.stack([np.cos(th)*np.sin(phi), np.sin(th)*np.sin(phi), np.cos(phi)], -1)
        means = np.array([cone_mean(r, w, d, ang)[0] for d in dirs]); k = np.nanargmin(means); kk = np.nanargmax(means)
        lb = lambda d: (np.degrees(np.arctan2(d[1], d[0])) % 360, np.degrees(np.arcsin(d[2])))
        return means[k], lb(dirs[k]), means[kk], lb(dirs[kk])
    def report(label, r, w, nshuf=2000):
        cm, nc = cone_mean(r, w, axis); allm = np.average(r, weights=w)
        shuf = np.array([cone_mean(rng.permutation(r), w, axis)[0] for _ in range(nshuf)])
        mn, dmn, mx, dmx = scan(r, w)
        # look-elsewhere: extreme of the scan under shuffles (fewer, coarser)
        ext = []
        for _ in range(200):
            rp = rng.permutation(r); a, _, bb, _ = scan(rp, w, nside=200); ext.append(max(abs(a - np.average(rp, weights=w)), abs(bb - np.average(rp, weights=w))))
        ext = np.array(ext); obs_ext = max(abs(mn - allm), abs(mx - allm))
        print(f"  {label}: Migkas cone (75 deg about l=280,b=-15): n={nc}, mean residual {cm-allm:+.3f} dex vs sky ({(cm-allm)/shuf.std():+.1f} sigma, shuffle sd {shuf.std():.3f}); "
              f"as H0 (Migkas convention, lower L = lower H0): {(cm-allm)*np.log(10)/2*100:+.1f}%. Scan: min {mn-allm:+.3f} at (l,b)=({dmn[0]:.0f},{dmn[1]:.0f}), max {mx-allm:+.3f} at ({dmx[0]:.0f},{dmx[1]:.0f}); "
              f"look-elsewhere P(extreme >= {obs_ext:.3f}) = {np.mean(ext >= obs_ext):.3f}")
    report("L-T only", r, w)
    rho, p = spearmanr(r, np.log10(exp)); print(f"  residual vs log EXP: Spearman {rho:+.3f} (p={p:.1e}); exposure in Migkas cone median {np.median(exp[n@axis >= np.cos(np.radians(75))]):.0f} s vs outside {np.median(exp[n@axis < np.cos(np.radians(75))]):.0f} s")
    beta2, r2, w2 = fit([np.log10(exp/np.median(exp))])
    print(f"  with exposure term: logL = ... + {beta2[2]:+.3f} log(EXP/median)")
    report("L-T + log EXP", r2, w2)
    q = np.quantile(exp, [0, 0.25, 0.5, 0.75, 1.0]); q[-1] += 1
    for i in range(4):
        m = (exp >= q[i]) & (exp < q[i+1]); rr = r[m] - np.average(r[m], weights=w[m]); nn = n[m]
        mc = nn@axis >= np.cos(np.radians(75))
        if mc.sum() < 20 or (~mc).sum() < 20: print(f"  EXP quartile {i+1} [{q[i]:.0f},{q[i+1]:.0f}) s: n={m.sum()}, cone n={mc.sum()}: too few"); continue
        cm = np.average(rr[mc], weights=w[m][mc]); sh = np.array([np.average(rng.permutation(rr)[mc], weights=w[m][mc]) for _ in range(1000)])
        print(f"  EXP quartile {i+1} [{q[i]:.0f},{q[i+1]:.0f}) s: n={m.sum()}, cone n={mc.sum()}: cone mean {cm:+.3f} dex ({cm/sh.std():+.1f} sigma)")
    # counts instead of exposure (what actually sets the T measurement)
    cts = t['CTS500'][idx]; okc = np.isfinite(cts) & (cts > 0)
    rho, p = spearmanr(r[okc], np.log10(cts[okc])); print(f"  residual vs log CTS500: Spearman {rho:+.3f} (p={p:.1e})")

main(sys.argv[1])

def extra(f):
    """collinearity check: the exposure-residual relation inside and outside the Migkas cone separately, and the cone
    deficit in a narrow exposure band where both sides are populated."""
    import numpy as np
    from astropy.io import fits
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from scipy.stats import spearmanr
    t = fits.open(f)[1].data
    kt, ktl, kth = t['KT'], t['KT_L'], t['KT_H']; L, Ll, Lh = t['L500'], t['L500_L'], t['L500_H']; z = t['BEST_Z']; exp = t['EXP']
    ekt = 0.5*(np.abs(kth - kt) + np.abs(kt - ktl)); eL = 0.5*(np.abs(Lh - L) + np.abs(L - Ll))
    ok = np.isfinite(kt) & (kt > 0) & np.isfinite(L) & (L > 0) & (z > 0.05) & (z < 0.6) & np.isfinite(exp) & (exp > 0) & np.isfinite(ekt) & np.isfinite(eL) & (ekt/kt < 0.5)
    idx = np.where(ok)[0]; kt, ekt, L, eL, z, exp = kt[idx], ekt[idx], L[idx], eL[idx], z[idx], exp[idx]
    sc = SkyCoord(ra=t['RA'][idx]*u.deg, dec=t['DEC'][idx]*u.deg).galactic; l, b = sc.l.deg, sc.b.deg
    n = np.stack([np.cos(np.radians(b))*np.cos(np.radians(l)), np.cos(np.radians(b))*np.sin(np.radians(l)), np.sin(np.radians(b))], -1)
    axis = np.array([np.cos(np.radians(-15))*np.cos(np.radians(280)), np.cos(np.radians(-15))*np.sin(np.radians(280)), np.sin(np.radians(-15))])
    inc = n@axis >= np.cos(np.radians(75))
    x = np.log10(kt/4.0); y = np.log10(L/1e44) - np.log10(Ez(z)); ey = eL/L/np.log(10); ex = ekt/kt/np.log(10)
    w = 1/(ey**2 + (2.0*ex)**2 + 0.25**2); X = np.column_stack([np.ones_like(x), x])
    beta = np.linalg.solve((X*w[:, None]).T@X, (X*w[:, None]).T@y); r = y - X@beta
    for name, m in [('inside cone', inc), ('outside cone', ~inc)]:
        rho, p = spearmanr(r[m], np.log10(exp[m])); print(f"  residual vs log EXP {name}: n={m.sum()}, Spearman {rho:+.3f} (p={p:.1e})")
    for lo, hi in [(100, 160), (120, 200), (150, 250)]:
        m = (exp >= lo) & (exp < hi); rr = r[m] - np.average(r[m], weights=w[m]); mc = inc[m]
        cm = np.average(rr[mc], weights=w[m][mc]); sh = np.array([np.average(rng.permutation(rr)[mc], weights=w[m][mc]) for _ in range(2000)])
        print(f"  EXP in [{lo},{hi}) s: n={m.sum()} (cone {mc.sum()}, outside {(~mc).sum()}): cone deficit {cm:+.3f} dex ({cm/sh.std():+.1f} sigma)")
    # residual vs sky WITHOUT exposure, restricted to the cone's own exposure range outside it? the converse: the sky deficit predicted by the exposure trend alone
    s1 = np.polyfit(np.log10(exp), r, 1, w=np.sqrt(w)); pred = np.polyval(s1, np.log10(exp))
    print(f"  deficit predicted by the exposure trend alone (cone mean of the fitted trend): {np.average(pred[inc], weights=w[inc]) - np.average(pred, weights=w):+.3f} dex (observed -0.039)")
extra(sys.argv[1])

def robust(f):
    """variants: T error cut, redshift range, cone angle, counts instead of exposure, Lbol, the cosmology subsample cut via EXT_LIKE."""
    import numpy as np
    from astropy.io import fits
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    t = fits.open(f)[1].data
    base = dict(terr=0.5, zlo=0.05, zhi=0.6, ang=75, depth='EXP', lum='L500', extlike=0)
    def run(**kw):
        c = dict(base); c.update(kw)
        kt, ktl, kth = t['KT'], t['KT_L'], t['KT_H']; L, Ll, Lh = t[c['lum']], t[c['lum']+'_L'], t[c['lum']+'_H']; z = t['BEST_Z']; dep = t[c['depth']]
        ekt = 0.5*(np.abs(kth - kt) + np.abs(kt - ktl)); eL = 0.5*(np.abs(Lh - L) + np.abs(L - Ll))
        ok = np.isfinite(kt) & (kt > 0) & np.isfinite(L) & (L > 0) & (z > c['zlo']) & (z < c['zhi']) & np.isfinite(dep) & (dep > 0) & np.isfinite(ekt) & np.isfinite(eL) & (ekt/kt < c['terr']) & (t['EXT_LIKE'] > c['extlike'])
        i = np.where(ok)[0]; kt, ekt, L, eL, z, dep = kt[i], ekt[i], L[i], eL[i], z[i], dep[i]
        sc = SkyCoord(ra=t['RA'][i]*u.deg, dec=t['DEC'][i]*u.deg).galactic; l, b = sc.l.deg, sc.b.deg
        n = np.stack([np.cos(np.radians(b))*np.cos(np.radians(l)), np.cos(np.radians(b))*np.sin(np.radians(l)), np.sin(np.radians(b))], -1)
        axis = np.array([np.cos(np.radians(-15))*np.cos(np.radians(280)), np.cos(np.radians(-15))*np.sin(np.radians(280)), np.sin(np.radians(-15))])
        inc = n@axis >= np.cos(np.radians(c['ang']))
        x = np.log10(kt/4.0); y = np.log10(L/1e44) - np.log10(Ez(z)); ey = eL/L/np.log(10); ex = ekt/kt/np.log(10)
        w = 1/(ey**2 + (2.0*ex)**2 + 0.25**2)
        out = []
        for cols in ([], [np.log10(dep/np.median(dep))]):
            X = np.column_stack([np.ones_like(x), x] + cols); beta = np.linalg.solve((X*w[:, None]).T@X, (X*w[:, None]).T@y); r = y - X@beta
            cm = np.average(r[inc], weights=w[inc]) - np.average(r, weights=w)
            sh = np.array([np.average(rng.permutation(r)[inc], weights=w[inc]) - np.average(r, weights=w) for _ in range(1500)])
            out.append((cm, cm/sh.std()))
        return len(i), inc.sum(), out
    print("\nrobustness (cone deficit, dex and sigma: without / with the depth term):")
    for label, kw in [("baseline", {}), ("T err < 30%", dict(terr=0.3)), ("z < 0.3", dict(zhi=0.3)), ("0.1 < z < 0.6", dict(zlo=0.1)), ("cone 60 deg", dict(ang=60)), ("cone 90 deg", dict(ang=90)),
                      ("depth = CTS500", dict(depth='CTS500')), ("Lbol500", dict(lum='Lbol500')), ("EXT_LIKE > 12", dict(extlike=12))]:
        N, nc, ((c0, s0), (c1, s1)) = run(**kw)
        print(f"  {label:16s} N={N:4d} cone={nc:4d}: {c0:+.3f} ({s0:+.1f} sigma)  ->  {c1:+.3f} ({s1:+.1f} sigma)")
robust(sys.argv[1])

def joint(f):
    """referee's test: residual = a + d cos(theta_axis) + c log EXP, jointly; the d-c correlation; single-term models compared;
    and the exposure trend fitted outside the cone only, predicting the cone mean."""
    import numpy as np
    from astropy.io import fits
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    t = fits.open(f)[1].data
    kt, ktl, kth = t['KT'], t['KT_L'], t['KT_H']; L, Ll, Lh = t['L500'], t['L500_L'], t['L500_H']; z = t['BEST_Z']; exp = t['EXP']
    ekt = 0.5*(np.abs(kth - kt) + np.abs(kt - ktl)); eL = 0.5*(np.abs(Lh - L) + np.abs(L - Ll))
    ok = np.isfinite(kt) & (kt > 0) & np.isfinite(L) & (L > 0) & (z > 0.05) & (z < 0.6) & np.isfinite(exp) & (exp > 0) & np.isfinite(ekt) & np.isfinite(eL) & (ekt/kt < 0.5)
    i = np.where(ok)[0]; kt, ekt, L, eL, z, exp = kt[i], ekt[i], L[i], eL[i], z[i], exp[i]
    sc = SkyCoord(ra=t['RA'][i]*u.deg, dec=t['DEC'][i]*u.deg).galactic; l, b = sc.l.deg, sc.b.deg
    n = np.stack([np.cos(np.radians(b))*np.cos(np.radians(l)), np.cos(np.radians(b))*np.sin(np.radians(l)), np.sin(np.radians(b))], -1)
    axis = np.array([np.cos(np.radians(-15))*np.cos(np.radians(280)), np.cos(np.radians(-15))*np.sin(np.radians(280)), np.sin(np.radians(-15))])
    cth = n@axis; le = np.log10(exp/np.median(exp))
    x = np.log10(kt/4.0); y = np.log10(L/1e44) - np.log10(Ez(z)); ey = eL/L/np.log(10); ex = ekt/kt/np.log(10)
    w = 1/(ey**2 + (1.71*ex)**2 + 0.47**2)   # referee: use the fitted slope and the observed rms
    def gls(cols):
        X = np.column_stack([np.ones_like(x), x] + cols); XtW = (X*w[:, None]).T
        beta = np.linalg.solve(XtW@X, XtW@y); cov = np.linalg.inv(XtW@X); r = y - X@beta
        return beta, cov, np.sum(w*r**2)
    b0, c0, chi0 = gls([]); bd, cd, chid = gls([cth]); be, ce, chie = gls([le]); bj, cj, chij = gls([cth, le])
    print(f"\nreferee's joint test (weights with slope 1.71 and 0.47 dex scatter):")
    print(f"  dipole only:   d = {bd[2]:+.4f} ± {np.sqrt(cd[2,2]):.4f} dex per unit cos(theta); dchi2 = {chi0-chid:.1f}")
    print(f"  exposure only: c = {be[2]:+.4f} ± {np.sqrt(ce[2,2]):.4f} dex per dex; dchi2 = {chi0-chie:.1f}")
    rho = cj[2,3]/np.sqrt(cj[2,2]*cj[3,3])
    print(f"  joint: d = {bj[2]:+.4f} ± {np.sqrt(cj[2,2]):.4f}, c = {bj[3]:+.4f} ± {np.sqrt(cj[3,3]):.4f}; parameter correlation {rho:+.3f}; dchi2 vs exposure-only = {chie-chij:.1f}, vs dipole-only = {chid-chij:.1f}")
    print(f"  raw correlation of cos(theta) and log EXP across the sample: {np.corrcoef(cth, le)[0,1]:+.3f}")
    # exposure trend fitted outside the 75-deg cone only, then predicted inside
    inc = cth >= np.cos(np.radians(75)); out = ~inc
    X = np.column_stack([np.ones(out.sum()), x[out], le[out]]); XtW = (X*w[out][:, None]).T; bo = np.linalg.solve(XtW@X, XtW@y[out])
    pred_in = bo[0] + bo[1]*x[inc] + bo[2]*le[inc]; pred_out = X@bo
    res_in = y[inc] - pred_in; res_out = y[out] - pred_out
    print(f"  exposure slope fitted OUTSIDE the cone only (n={out.sum()}): c = {bo[2]:+.4f}; predicted-vs-observed inside the cone: mean residual {np.average(res_in, weights=w[inc]):+.4f} ± {np.sqrt(1/np.sum(w[inc])):.4f} dex (outside, by construction, {np.average(res_out, weights=w[out]):+.4f})")
    # a dipole model fitted with the same outside-only approach cannot be done (cos theta barely varies outside); say so with its range
    print(f"  cos(theta) range: inside [{cth[inc].min():.2f}, {cth[inc].max():.2f}], outside [{cth[out].min():.2f}, {cth[out].max():.2f}]; log EXP range inside [{le[inc].min():.2f}, {le[inc].max():.2f}], outside [{le[out].min():.2f}, {le[out].max():.2f}]")
joint(sys.argv[1])
