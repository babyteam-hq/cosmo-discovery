"""Cosmicflows-4 at the catalogue level (R-119, sweep-3 rows 1 and 2).
Part 1: do the FP distance moduli (SDSS r-band north, 6dFGSv JHK south) carry a Galactic-extinction-dependent residual?
Part 2: hemispheric zero points between methods (flow-free, from groups and galaxies holding two methods), shell dipoles with
one monopole vs a free monopole per survey, and the bulk flow a measured SDSS/6dF step manufactures on CF4's sky.
usage: python cf4_catalogue.py [part1|part2|all]   (data in data/: cf4_table2.tsv, cf4_groups.tsv, dust maps)"""
import sys, numpy as np, pandas as pd
from cf4_io import *
rng = np.random.default_rng(2)
K = 5 / np.log(10)  # mag per unit v/cz

def wls(y, X, w):
    Xw = X * w[:, None]; A = X.T @ Xw; b = np.linalg.solve(A, Xw.T @ y); r = y - X @ b
    chi2 = np.sum(w * r * r); dof = len(y) - X.shape[1]
    cov = np.linalg.inv(A) * (chi2 / dof)  # errors rescaled to chi2/dof = 1
    return b, np.sqrt(np.diag(cov)), chi2 / dof, cov

def load():
    d = vizier_tsv('data/cf4_table2.tsv'); g = vizier_tsv('data/cf4_groups.tsv')
    d['absb'] = np.abs(d.GLAT); d['n'] = list(unit_eq(d.RAJ2000.values, d.DEJ2000.values))
    d['survey_fp'] = fp_survey(d)
    return d, g

def sample(d, dm, edm, vmin=3000, vmax=30000, bmin=10):
    m = d[dm].notna() & (d.Vcmb > vmin) & (d.Vcmb < vmax) & (d.absb > bmin)
    s = d[m].copy(); s['r'] = s[dm] - hubble_dm(s.Vcmb); s['e'] = s[edm].clip(lower=0.05); return s

def design(s, ebv_cols=(), dipole=True, vbins=True):
    cols = [np.ones(len(s))]; names = ['const']
    if vbins:  # a smooth trend with distance (Malmquist-type), 2nd-order in log Vcmb
        lv = np.log10(s.Vcmb.values) - 4; cols += [lv, lv * lv]; names += ['logV', 'logV2']
    if dipole:
        n = np.stack(s.n.values); f = -K / s.Vcmb.values[:, None] * 1000  # mag per 1000 km/s of flow
        cols += [n[:, 0] * f[:, 0], n[:, 1] * f[:, 0], n[:, 2] * f[:, 0]]; names += ['Bx', 'By', 'Bz']
    for c in ebv_cols: cols.append(s[c].values); names.append(c)
    return np.column_stack(cols), names

def part1(d):
    print("=== Part 1: FP (and TF, SN) residual vs Galactic E(B-V) ===")
    l, b = d.GLON.values, d.GLAT.values
    d['sfd'] = sfd_ebv(l, b); d['gnilc'] = gnilc_ebv(l, b); d['zhou'], d['zsfd'] = zhou_ebv(l, b)
    d['g_s'] = d.gnilc - d.sfd; d['z_s'] = d.zhou - d.zsfd; d['cscb'] = 1 / np.sin(np.radians(d.absb))
    samples = {'FP SDSS': ('DMfp', 'e_DMfp', 'sdss'), 'FP 6dF': ('DMfp', 'e_DMfp', '6df'), 'TF': ('DMtf', 'e_DMtf', None),
               'SN Ia': ('DMsnIa', 'e_DMsnIa', None)}
    out = {}
    for name, (dm, edm, sv) in samples.items():
        s = sample(d, dm, edm)
        if sv: s = s[s.survey_fp == sv]
        w = 1 / s.e.values ** 2
        X, nm = design(s, ['sfd']); bb, ee, c2, _ = wls(s.r.values, X, w)
        i = nm.index('sfd'); out[name] = (bb[i], ee[i])
        print(f"\n{name}: N={len(s)} rms {np.sqrt(np.average(s.r**2, weights=w)):.3f} E(B-V) mean {s.sfd.mean():.3f} rms {s.sfd.std():.3f}")
        print(f"  slope on SFD: {bb[i]:+.3f} ± {ee[i]:.3f} mag per E(B-V) (chi2/dof {c2:.2f}); dipole B=({bb[nm.index('Bx')]*1000:.0f},"
              f"{bb[nm.index('By')]*1000:.0f},{bb[nm.index('Bz')]*1000:.0f}) km/s")
        # rotation null: the map resampled at positions rotated in l (and reflected in b half the time)
        zs = []
        for _ in range(200):
            lr = (s.GLON.values + rng.uniform(0, 360)) % 360; br = s.GLAT.values * rng.choice([-1, 1])
            s2 = s.copy(); s2['sfd'] = sfd_ebv(lr, br); X2, _ = design(s2, ['sfd']); b2, e2, *_ = wls(s.r.values, X2, w)
            zs.append(b2[i] / e2[i])
        zs = np.array(zs); print(f"  rotation null: z spread {zs.std():.2f}, data z {bb[i]/ee[i]:+.2f} -> {bb[i]/ee[i]/zs.std():+.2f} sigma-equivalent")
        # with csc|b| (a latitude trend of any origin) and in |b| and Vcmb bands
        X, nm = design(s, ['sfd', 'cscb']); b3, e3, *_ = wls(s.r.values, X, w)
        print(f"  with csc|b|: SFD {b3[nm.index('sfd')]:+.3f} ± {e3[nm.index('sfd')]:.3f}, csc|b| {b3[nm.index('cscb')]:+.4f} ± {e3[nm.index('cscb')]:.4f}")
        for lab, m in [('|b|>30', s.absb > 30), ('|b|<30', s.absb <= 30), ('V<10000', s.Vcmb < 10000), ('V>10000', s.Vcmb >= 10000),
                       ('E<0.05', s.sfd < 0.05)]:
            t = s[m]; X, nm = design(t, ['sfd']); b4, e4, *_ = wls(t.r.values, X, 1 / t.e.values ** 2)
            print(f"    {lab}: {b4[nm.index('sfd')]:+.3f} ± {e4[nm.index('sfd')]:.3f} (N={len(t)})")
        # map discrimination: SFD plus (GNILC - SFD); SFD plus (Zhou - SFD) where the DESI map has stars
        X, nm = design(s, ['sfd', 'g_s']); b5, e5, *_ = wls(s.r.values, X, w)
        print(f"  SFD + (GNILC-SFD): {b5[nm.index('sfd')]:+.3f} ± {e5[nm.index('sfd')]:.3f}, {b5[nm.index('g_s')]:+.3f} ± {e5[nm.index('g_s')]:.3f}")
        t = s[s.zhou.notna()]
        if len(t) > 500:
            X, nm = design(t, ['sfd', 'z_s']); b6, e6, *_ = wls(t.r.values, X, 1 / t.e.values ** 2)
            X, nm2 = design(t, ['sfd', 'z_s', 'g_s']); b7, e7, *_ = wls(t.r.values, X, 1 / t.e.values ** 2)
            print(f"  in the DESI-stars footprint (N={len(t)}): SFD + (Zhou-SFD): {b6[nm.index('sfd')]:+.3f} ± {e6[nm.index('sfd')]:.3f}, "
                  f"{b6[nm.index('z_s')]:+.3f} ± {e6[nm.index('z_s')]:.3f}; with GNILC-SFD too: Zhou {b7[nm2.index('z_s')]:+.3f} ± {e7[nm2.index('z_s')]:.3f}, "
                  f"GNILC {b7[nm2.index('g_s')]:+.3f} ± {e7[nm2.index('g_s')]:.3f}")
    return out

def hemi_offsets(g, d):
    print("\n=== Part 2a: flow-free zero-point offsets between methods, by celestial hemisphere ===")
    print("(pairs: the same group or galaxy holding both methods; weighted mean of DM_a - DM_b, errors rescaled to chi2/dof=1)")
    res = {}
    for tab, lab in [(g, 'groups'), (d, 'galaxies')]:
        for a, b in [('DMfp', 'DMtf'), ('DMfp', 'DMsnIa'), ('DMtf', 'DMsnIa'), ('DMfp', 'DMsbf')]:
            if a not in tab or b not in tab: continue
            m = tab[a].notna() & tab[b].notna() & (tab.Vcmb > 1000)
            t = tab[m]; diff = (t[a] - t[b]).values; e = np.sqrt(t['e_' + a].fillna(0.4) ** 2 + t['e_' + b].fillna(0.4) ** 2).clip(0.05).values
            # per-group errors in the group table are of the means; for galaxies use the listed errors
            line = f"{lab} {a[2:]}-{b[2:]}: "
            vals = {}
            for h, hm in [('N', t.DEJ2000 > 0), ('S', t.DEJ2000 < 0), ('all', np.ones(len(t), bool))]:
                dd, ee = diff[hm], e[hm]
                if len(dd) < 5: continue
                w = 1 / ee ** 2; mu = np.sum(w * dd) / np.sum(w); c2 = np.sum(w * (dd - mu) ** 2) / max(len(dd) - 1, 1)
                err = np.sqrt(max(c2, 1) / np.sum(w)); vals[h] = (mu, err, len(dd))
                line += f"{h} {mu:+.3f}±{err:.3f} (n={len(dd)}, chi2/dof {c2:.1f})  "
            if 'N' in vals and 'S' in vals:
                dns = vals['N'][0] - vals['S'][0]; ens = np.hypot(vals['N'][1], vals['S'][1])
                line += f"N-S {dns:+.3f}±{ens:.3f} ({dns/ens:+.1f}σ)"; res[(lab, a, b)] = (dns, ens)
            print(line)
    # FP members in SDSS vs 6dF within the same group? (none by construction: different skies). Instead: groups holding FP
    # and TF, split by whether the FP members are SDSS or 6dF
    return res

def shell_dipoles(d, inject=None, verbose=True):
    """Weighted dipole fits of r = DM - DM_H in Vcmb shells. Returns {shell: {fit: (B, eB)}} with B in km/s.
    Fits: 'one' = one monopole; 'method' = a monopole per method (FP tied across hemispheres); 'survey' = a monopole per survey
    (SDSS-FP, 6dF-FP, equatorial FP, TF, SN, other), i.e. the SDSS/6dF tie free."""
    rows = []
    for dm, edm, meth in [('DMfp', 'e_DMfp', 'fp'), ('DMtf', 'e_DMtf', 'tf'), ('DMsnIa', 'e_DMsnIa', 'sn'), ('DMsbf', 'e_DMsbf', 'other'),
                          ('DMsnII', 'e_DMsnII', 'other'), ('DMtrgb', 'e_DMtrgb', 'other')]:
        s = sample(d, dm, edm, vmin=2000, vmax=32000); s['method'] = meth
        s['survey'] = np.where(meth == 'fp', s.survey_fp, meth); rows.append(s[['r', 'e', 'Vcmb', 'n', 'method', 'survey', 'DEJ2000']])
    s = pd.concat(rows); shells = [(2000, 6000), (6000, 10000), (10000, 15000), (15000, 20000), (20000, 32000)]
    out = {}
    for lo, hi in shells:
        t = s[(s.Vcmb >= lo) & (s.Vcmb < hi)].copy(); y = t.r.values.copy(); w = 1 / t.e.values ** 2
        if inject is not None:  # replace the data by isotropic noise plus a zero-point step on the 6dF FP galaxies
            y = rng.normal(0, t.e.values) + np.where(t.survey == '6df', inject, 0)
        n = np.stack(t.n.values); f = -K / t.Vcmb.values * 1000
        dip = np.column_stack([n[:, 0] * f, n[:, 1] * f, n[:, 2] * f])
        res = {}
        for fit, key in [('one', None), ('method', 'method'), ('survey', 'survey')]:
            if key is None: M = np.ones((len(t), 1))
            else:
                cats = sorted(t[key].unique()); M = np.column_stack([(t[key] == c).values.astype(float) for c in cats])
            X = np.column_stack([M, dip]); bb, ee, c2, cov = wls(y, X, w)
            B = bb[-3:] * 1000; eB = ee[-3:] * 1000; res[fit] = (B, eB)
            if verbose:
                ll, bl = gal_lb(B[None, :] / np.linalg.norm(B)); counts = t[key].value_counts().to_dict() if key else len(t)
                print(f"  {lo}-{hi} {fit:6s}: |B|={np.linalg.norm(B):4.0f} toward (l,b)=({ll[0]:.0f},{bl[0]:.0f}); B=({B[0]:+.0f}±{eB[0]:.0f}, "
                      f"{B[1]:+.0f}±{eB[1]:.0f}, {B[2]:+.0f}±{eB[2]:.0f}) km/s; Bz (celestial pole) {B[2]:+.0f}±{eB[2]:.0f}; N={len(t)}"
                      + (f"; zero points {dict(zip(cats, np.round(bb[:len(cats)], 3)))}" if key else ''))
        out[(lo, hi)] = res
    return out

def part2(d, g):
    off = hemi_offsets(g, d)
    print("\n=== Part 2b: shell dipoles, one monopole vs per-method vs per-survey monopoles (B in equatorial km/s) ===")
    data = shell_dipoles(d)
    # the measured step: FP N-S relative to the all-sky methods, from the group pairs
    step = np.nanmean([off[k][0] for k in off if k[0] == 'groups' and k[1] == 'DMfp' and k[2] in ('DMtf', 'DMsnIa')])
    print(f"\n=== Part 2c: injection of an isotropic catalogue with a 6dF-FP zero-point step of {step:+.3f} mag (the measured FP N-S "
          f"offset relative to TF/SN from the group pairs; sign: 6dF lower by that much), 100 mocks, single-monopole fit ===")
    acc = {}
    for i in range(100):
        m = shell_dipoles(d, inject=-step, verbose=False)
        for sh in m:
            acc.setdefault(sh, []).append(m[sh]['one'][0])
    for sh in acc:
        a = np.array(acc[sh]); mean = a.mean(0); sd = a.std(0)
        ll, bl = gal_lb(mean[None, :] / np.linalg.norm(mean))
        print(f"  {sh[0]}-{sh[1]}: manufactured B = ({mean[0]:+.0f}, {mean[1]:+.0f}, {mean[2]:+.0f}) ± ({sd[0]:.0f}, {sd[1]:.0f}, {sd[2]:.0f}) km/s, "
              f"|B| {np.linalg.norm(mean):.0f} toward ({ll[0]:.0f},{bl[0]:.0f}); data one-monopole B = {np.round(data[sh]['one'][0]).astype(int)}, "
              f"per-survey B = {np.round(data[sh]['survey'][0]).astype(int)}")
    acc0 = {}
    for i in range(100):
        m = shell_dipoles(d, inject=0.0, verbose=False)
        for sh in m: acc0.setdefault(sh, []).append(m[sh]['one'][0])
    print("  no-step mocks (noise only): |B| mean per shell", {sh: round(float(np.mean(np.linalg.norm(acc0[sh], axis=1)))) for sh in acc0})

if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    d, g = load()
    if what in ('part1', 'all'): part1(d)
    if what in ('part2', 'all'): part2(d, g)

def part2x(d, g):
    """Diagnostics: is the FP-vs-TF/SN offset a step at the celestial equator (SDSS vs 6dF) or a gradient in Dec; and the
    shell dipoles from each subsample alone (free monopole each) and from everything but the FP."""
    print("\n=== Part 2x: flow-free FP - TF and FP - SN group offsets vs declination ===")
    for b in ['DMtf', 'DMsnIa']:
        t = g[g.DMfp.notna() & g[b].notna() & (g.Vcmb > 1000)]
        diff = (t.DMfp - t[b]).values; e = np.sqrt(t.e_DMfp.fillna(0.4) ** 2 + t['e_' + b].fillna(0.4) ** 2).clip(0.05).values
        line = f"FP-{b[2:]}: "
        for lo, hi in [(-90, -45), (-45, -20), (-20, -5), (-5, 5), (5, 20), (20, 45), (45, 90)]:
            m = (t.DEJ2000 >= lo) & (t.DEJ2000 < hi)
            if m.sum() < 8: line += f"[{lo},{hi}) n={m.sum()}  "; continue
            w = 1 / e[m] ** 2; mu = np.sum(w * diff[m]) / np.sum(w); c2 = np.sum(w * (diff[m] - mu) ** 2) / (m.sum() - 1)
            line += f"[{lo},{hi}) {mu:+.3f}±{np.sqrt(max(c2,1)/np.sum(w)):.3f} (n={m.sum()})  "
        print(line)
        # by Vcmb within each hemisphere (does the step depend on distance?)
        line = f"  by Vcmb: "
        for h, hm in [('N', t.DEJ2000 > 0), ('S', t.DEJ2000 < 0)]:
            for lo, hi in [(1000, 5000), (5000, 10000), (10000, 20000)]:
                m = hm & (t.Vcmb >= lo) & (t.Vcmb < hi)
                if m.sum() < 8: continue
                w = 1 / e[m] ** 2; mu = np.sum(w * diff[m]) / np.sum(w); c2 = np.sum(w * (diff[m] - mu) ** 2) / (m.sum() - 1)
                line += f"{h} {lo}-{hi}: {mu:+.3f}±{np.sqrt(max(c2,1)/np.sum(w)):.3f} (n={m.sum()})  "
        print(line)
    print("\n=== Part 2x: shell dipoles per subsample alone (free monopole), and without any FP ===")
    rows = []
    for dm, edm, meth in [('DMfp', 'e_DMfp', 'fp'), ('DMtf', 'e_DMtf', 'tf'), ('DMsnIa', 'e_DMsnIa', 'sn'), ('DMsbf', 'e_DMsbf', 'other'),
                          ('DMsnII', 'e_DMsnII', 'other'), ('DMtrgb', 'e_DMtrgb', 'other')]:
        s = sample(d, dm, edm, vmin=2000, vmax=32000); s['method'] = meth
        s['survey'] = np.where(meth == 'fp', s.survey_fp, meth); rows.append(s[['r', 'e', 'Vcmb', 'n', 'method', 'survey']])
    s = pd.concat(rows)
    for lo, hi in [(2000, 6000), (6000, 10000), (10000, 15000), (15000, 20000), (20000, 32000)]:
        t0 = s[(s.Vcmb >= lo) & (s.Vcmb < hi)]
        for lab, m in [('SDSS FP', t0.survey == 'sdss'), ('6dF FP', t0.survey == '6df'), ('TF', t0.survey == 'tf'), ('SN', t0.survey == 'sn'),
                       ('no FP', t0.method != 'fp'), ('TF+SN', t0.method.isin(['tf', 'sn']))]:
            t = t0[m]
            if len(t) < 30: print(f"  {lo}-{hi} {lab:8s}: N={len(t)} (too few)"); continue
            n = np.stack(t.n.values); f = -K / t.Vcmb.values * 1000
            X = np.column_stack([np.ones(len(t)), n[:, 0] * f, n[:, 1] * f, n[:, 2] * f]); bb, ee, c2, _ = wls(t.r.values, X, 1 / t.e.values ** 2)
            B, eB = bb[1:] * 1000, ee[1:] * 1000; ll, bl = gal_lb(B[None, :] / np.linalg.norm(B))
            print(f"  {lo}-{hi} {lab:8s}: N={len(t):5d} |B|={np.linalg.norm(B):4.0f} toward ({ll[0]:.0f},{bl[0]:.0f}); B=({B[0]:+.0f}±{eB[0]:.0f}, "
                  f"{B[1]:+.0f}±{eB[1]:.0f}, {B[2]:+.0f}±{eB[2]:.0f})")

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'part2x':
    d, g = load(); part2x(d, g)

def part3(d, g):
    """Tie parts 1 and 2 together: (i) flow-free FP - TF group offsets vs SFD E(B-V) and vs |b|; (ii) the SDSS FP slope on
    E(B-V) in Vcmb bins (a magnitude-limit selection effect grows with distance, an extinction-coefficient error does not);
    (iii) shell dipoles with an E(B-V) term per FP survey; (iv) the N-S step after the E(B-V) term."""
    print("\n=== Part 3(i): FP - TF group pairs vs E(B-V) and |b| (flow-free) ===")
    g = g.copy(); g['sfd'] = sfd_ebv(g.GLON.values, g.GLAT.values); g['absb'] = np.abs(g.GLAT)
    t = g[g.DMfp.notna() & g.DMtf.notna() & (g.Vcmb > 1000)]
    diff = (t.DMfp - t.DMtf).values; e = np.sqrt(t.e_DMfp.fillna(0.4) ** 2 + t.e_DMtf.fillna(0.4) ** 2).clip(0.05).values; w = 1 / e ** 2
    north = (t.DEJ2000 > 0).values.astype(float)
    for lab, X, nm in [('const + north', np.column_stack([np.ones(len(t)), north]), ['c', 'north']),
                       ('const + north + E(B-V)', np.column_stack([np.ones(len(t)), north, t.sfd.values]), ['c', 'north', 'ebv']),
                       ('const + E(B-V)', np.column_stack([np.ones(len(t)), t.sfd.values]), ['c', 'ebv']),
                       ('const + north + csc|b|', np.column_stack([np.ones(len(t)), north, 1 / np.sin(np.radians(t.absb.values))]), ['c', 'north', 'cscb'])]:
        bb, ee, c2, _ = wls(diff, X, w)
        print(f"  {lab}: " + ", ".join(f"{n} {b:+.3f}±{s:.3f}" for n, b, s in zip(nm, bb, ee)) + f"  (E(B-V) mean N {t.sfd[t.DEJ2000>0].mean():.3f} S {t.sfd[t.DEJ2000<0].mean():.3f})")
    print("\n=== Part 3(ii): SDSS FP residual slope on E(B-V) in Vcmb bins ===")
    s = sample(d, 'DMfp', 'e_DMfp'); s = s[s.survey_fp == 'sdss'].copy(); s['sfd'] = sfd_ebv(s.GLON.values, s.GLAT.values)
    for lo, hi in [(3000, 8000), (8000, 12000), (12000, 16000), (16000, 20000), (20000, 25000), (25000, 30000)]:
        t2 = s[(s.Vcmb >= lo) & (s.Vcmb < hi)]
        X, nm = design(t2, ['sfd']); bb, ee, *_ = wls(t2.r.values, X, 1 / t2.e.values ** 2)
        print(f"  {lo}-{hi}: {bb[nm.index('sfd')]:+.3f} ± {ee[nm.index('sfd')]:.3f} (N={len(t2)}, rms {t2.r.std():.3f})")
    s6 = sample(d, 'DMfp', 'e_DMfp'); s6 = s6[s6.survey_fp == '6df'].copy(); s6['sfd'] = sfd_ebv(s6.GLON.values, s6.GLAT.values)
    for lo, hi in [(3000, 8000), (8000, 12000), (12000, 16500)]:
        t2 = s6[(s6.Vcmb >= lo) & (s6.Vcmb < hi)]
        X, nm = design(t2, ['sfd']); bb, ee, *_ = wls(t2.r.values, X, 1 / t2.e.values ** 2)
        print(f"  6dF {lo}-{hi}: {bb[nm.index('sfd')]:+.3f} ± {ee[nm.index('sfd')]:.3f} (N={len(t2)})")
    print("\n=== Part 3(iii): shell dipoles with an E(B-V) term per FP survey (one monopole; per-survey monopoles) ===")
    rows = []
    for dm, edm, meth in [('DMfp', 'e_DMfp', 'fp'), ('DMtf', 'e_DMtf', 'tf'), ('DMsnIa', 'e_DMsnIa', 'sn'), ('DMsbf', 'e_DMsbf', 'other'),
                          ('DMsnII', 'e_DMsnII', 'other'), ('DMtrgb', 'e_DMtrgb', 'other')]:
        t3 = sample(d, dm, edm, vmin=2000, vmax=32000); t3['method'] = meth; t3['survey'] = np.where(meth == 'fp', t3.survey_fp, meth)
        t3['sfd'] = sfd_ebv(t3.GLON.values, t3.GLAT.values); rows.append(t3[['r', 'e', 'Vcmb', 'n', 'method', 'survey', 'sfd']])
    s = pd.concat(rows)
    for lo, hi in [(2000, 6000), (6000, 10000), (10000, 15000), (15000, 20000), (20000, 32000)]:
        t = s[(s.Vcmb >= lo) & (s.Vcmb < hi)]; w = 1 / t.e.values ** 2
        n = np.stack(t.n.values); f = -K / t.Vcmb.values * 1000; dip = np.column_stack([n[:, 0] * f, n[:, 1] * f, n[:, 2] * f])
        ebv = np.column_stack([t.sfd.values * (t.survey == 'sdss').values, t.sfd.values * (t.survey == '6df').values])
        for fit, key in [('one+ebv', None), ('survey+ebv', 'survey')]:
            if key is None: M = np.ones((len(t), 1)); cats = ['all']
            else: cats = sorted(t[key].unique()); M = np.column_stack([(t[key] == c).values.astype(float) for c in cats])
            X = np.column_stack([M, dip, ebv]); bb, ee, c2, _ = wls(t.r.values, X, w)
            B, eB = bb[len(cats):len(cats) + 3] * 1000, ee[len(cats):len(cats) + 3] * 1000; ll, bl = gal_lb(B[None, :] / np.linalg.norm(B))
            zp = {c: round(float(v), 3) for c, v in zip(cats, bb[:len(cats)])} if key else ''
            print(f"  {lo}-{hi} {fit:10s}: |B|={np.linalg.norm(B):4.0f} toward ({ll[0]:.0f},{bl[0]:.0f}); B=({B[0]:+.0f}±{eB[0]:.0f}, {B[1]:+.0f}±{eB[1]:.0f}, "
                  f"{B[2]:+.0f}±{eB[2]:.0f}); E(B-V) slopes sdss {bb[-2]:+.2f}±{ee[-2]:.2f}, 6df {bb[-1]:+.2f}±{ee[-1]:.2f} {zp}")

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'part3':
    d, g = load(); part3(d, g)
