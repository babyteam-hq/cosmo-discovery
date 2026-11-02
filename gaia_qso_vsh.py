"""Gaia DR3 quasar proper motions (the Gaia-CRF3 sources, qso_candidates.gaia_crf_source) fitted with vector spherical
harmonics to degree 2 (glide g and spin w at l = 1, ten l = 2 terms), per-source 2x2 covariance, in bins of SPECTROSCOPIC
redshift (DESI DR1 QSO clustering catalogue and SDSS DR16Q, matched within 1"), and in bins of G, BP-RP, RUWE, astrometric
excess noise, visibility periods and |ecliptic latitude| (sweep-4 row C, the open part of the trajectory verdict: the
redshift-binned glide/spin test with spectroscopic redshifts, and the spin/glide terms against Gaia metadata, with a
permutation null). The glide is the secular aberration (Solar System acceleration, 5.05 uas/yr toward the Galactic centre
in Gaia-CRF3); the spin should be zero by the frame's construction.
usage: python -I gaia_qso_vsh.py data/gaia data/desi_qso data/sdss/dr16q.tsv [nperm]
(data/gaia holds the TAP chunks crf3_*.csv, VOTables despite the name; data/desi_qso the QSO_{NGC,SGC}_clustering.dat.fits)"""
import sys, glob, numpy as np, pandas as pd
from astropy.io.votable import parse_single_table
from astropy.io import fits
from scipy.spatial import cKDTree
rng = np.random.default_rng(3)
GC = (266.405, -28.936)  # RA, Dec of the Galactic centre

def unit(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], -1)

def load_gaia(d):
    parts = []
    for f in sorted(glob.glob(d + '/crf3_*.csv')):
        try: parts.append(parse_single_table(f).to_table().to_pandas())
        except Exception as e: print(f"  skip {f}: {e}")
    g = pd.concat(parts, ignore_index=True)
    g = g[g.pmra.notna() & g.pmdec.notna() & (g.astrometric_params_solved >= 31)].reset_index(drop=True)
    for c in ['pmra', 'pmdec', 'pmra_error', 'pmdec_error']: g[c] = g[c] * 1000.0  # mas/yr -> uas/yr
    return g

def design(ra, dec):
    """Columns: glide (gx, gy, gz), spin (wx, wy, wz), then 5 spheroidal and 5 toroidal l=2 terms; rows: (mu_a*, mu_d)."""
    a, d = np.radians(ra), np.radians(dec); sa, ca, sd, cd = np.sin(a), np.cos(a), np.sin(d), np.cos(d)
    A = [(-sa, -ca * sd), (ca, -sa * sd), (0 * a, cd),               # glide
         (ca * sd, -sa), (sa * sd, ca), (-cd, 0 * a)]                 # spin
    s2 = [(0 * a, 2 * sd * cd), (-sd * sa, np.cos(2 * d) * ca), (sd * ca, np.cos(2 * d) * sa),
          (-2 * cd * np.sin(2 * a), -2 * sd * cd * np.cos(2 * a)), (2 * cd * np.cos(2 * a), -2 * sd * cd * np.sin(2 * a))]
    A += s2 + [(-y, x) for x, y in s2]
    Xa = np.column_stack([c[0] for c in A]); Xd = np.column_stack([c[1] for c in A])
    return Xa, Xd

def fit(g, idx=None, pm=None):
    s = g if idx is None else g.iloc[idx]
    Xa, Xd = design(s.ra.values, s.dec.values)
    ma, md = (s.pmra.values, s.pmdec.values) if pm is None else pm
    ea, ed, r = s.pmra_error.values, s.pmdec_error.values, np.nan_to_num(s.pmra_pmdec_corr.values)
    # whiten each source's 2-vector with the inverse Cholesky of its covariance
    L11 = ea; L21 = r * ed; L22 = ed * np.sqrt(np.maximum(1 - r * r, 1e-6))
    ya = ma / L11; yd = (md - L21 * ya) / L22
    Wa = Xa / L11[:, None]; Wd = (Xd - L21[:, None] * Wa) / L22[:, None]
    X = np.vstack([Wa, Wd]); y = np.concatenate([ya, yd])
    XtX = X.T @ X; b = np.linalg.solve(XtX, X.T @ y); cov = np.linalg.inv(XtX)
    res = y - X @ b; chi2 = res @ res; dof = len(y) - len(b)
    return b, cov, chi2 / dof

NGC = unit(*GC)
def proj(b, cov):
    """Signed projection of the glide on the Galactic-centre direction (unbiased by noise, unlike the amplitude)."""
    return b[:3] @ NGC, np.sqrt(NGC @ cov[:3, :3] @ NGC)

def report(b, cov, red, n, label):
    g, w = b[:3], b[3:6]; eg = np.sqrt(np.diag(cov)[:3]); ew = np.sqrt(np.diag(cov)[3:6])
    G = np.linalg.norm(g); W = np.linalg.norm(w); P, eP = proj(b, cov)
    eG = np.sqrt(g @ cov[:3, :3] @ g) / G; eW = np.sqrt(w @ cov[3:6, 3:6] @ w) / max(W, 1e-9)
    ra = np.degrees(np.arctan2(g[1], g[0])) % 360; dec = np.degrees(np.arcsin(g[2] / G))
    sep = np.degrees(np.arccos(np.clip(unit(ra, dec) @ unit(*GC), -1, 1)))
    print(f"  {label:38s} n={n:7d}  glide {G:5.2f} ± {eG:4.2f} uas/yr toward ({ra:5.1f}, {dec:+5.1f}), {sep:4.1f} deg from the GC, GC projection {P:+5.2f} ± {eP:4.2f}; "
          f"spin |w| {W:5.2f} ± {eW:4.2f} (wx {w[0]:+5.2f}±{ew[0]:.2f}, wy {w[1]:+5.2f}±{ew[1]:.2f}, wz {w[2]:+5.2f}±{ew[2]:.2f}); chi2/dof {red:.2f}")
    return G, eG, W, eW, P, eP

def spec_z(g, desi_dir, dr16q):
    tabs = []
    for f in sorted(glob.glob(desi_dir + '/QSO_*_clustering.dat.fits')):
        t = fits.getdata(f); tabs.append(pd.DataFrame({'ra': t['RA'], 'dec': t['DEC'], 'z': t['Z']}))
    desi = pd.concat(tabs, ignore_index=True)
    rows = [l for l in open(dr16q, encoding='utf-8', errors='ignore') if not l.startswith('#') and l.strip()]
    hdr = rows[0].rstrip('\n').split('\t'); body = [r.rstrip('\n').split('\t') for r in rows[3:]]
    s = pd.DataFrame(body, columns=hdr).replace('', np.nan)
    for c in s.columns: s[c] = pd.to_numeric(s[c], errors='coerce')
    s = s.rename(columns={'_RAJ2000': 'ra', '_DEJ2000': 'dec', 'RAdeg': 'ra', 'DEdeg': 'dec'}); s = s[(s.z > 0) & (s.zwarn == 0)] if 'zwarn' in s else s[s.z > 0]
    z = np.full(len(g), np.nan); src = np.zeros(len(g), int)
    ug = unit(g.ra.values, g.dec.values); tree = cKDTree(ug); tol = 2 * np.sin(np.radians(1 / 3600) / 2)
    for k, cat in [(2, s), (1, desi)]:  # DESI last so it wins where both match
        dd, ii = tree.query(unit(cat.ra.values, cat.dec.values), distance_upper_bound=tol)
        ok = np.isfinite(dd); z[ii[ok]] = cat.z.values[ok]; src[ii[ok]] = k
    print(f"spectroscopic redshifts: DESI DR1 QSO {len(desi)} rows, DR16Q {len(s)} rows; matched {np.isfinite(z).sum()} of {len(g)} "
          f"({(src == 1).sum()} DESI, {(src == 2).sum()} SDSS)")
    return z

def main(gdir, ddir, dr16q, nperm=50):
    g = load_gaia(gdir); print(f"Gaia-CRF3 sources with 5/6-parameter solutions: {len(g)}")
    print("== all sources"); b, cov, red = fit(g); report(b, cov, red, len(g), 'all')
    # permutation null: proper motions (with their errors) shuffled among sources at fixed positions
    Gs, Ws = [], []
    for _ in range(nperm):
        p = rng.permutation(len(g)); gp = g.copy()
        for c in ['pmra', 'pmdec', 'pmra_error', 'pmdec_error', 'pmra_pmdec_corr']: gp[c] = g[c].values[p]
        bp, cp, _ = fit(gp); Gs.append(np.linalg.norm(bp[:3])); Ws.append(np.linalg.norm(bp[3:6]))
    print(f"  permutation null ({nperm}): glide amplitude {np.mean(Gs):.2f} ± {np.std(Gs):.2f}, spin {np.mean(Ws):.2f} ± {np.std(Ws):.2f} (noise bias of the amplitudes)")
    g['zspec'] = spec_z(g, ddir, dr16q)
    def bins(col, edges, label, sel=None, nullperm=0):
        print(f"== by {label}")
        out = []
        base = np.ones(len(g), bool) if sel is None else sel
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = base & (g[col].values >= lo) & (g[col].values < hi)
            if m.sum() < 2000: continue
            b, cov, red = fit(g, np.where(m)[0]); G, eG, W, eW, P, eP = report(b, cov, red, m.sum(), f"{label} [{lo:g}, {hi:g})")
            out.append((0.5 * (lo + hi), P, eP, W, eW, m.sum()))
            if nullperm:
                sub = g[m]; Gn = []
                for _ in range(nullperm):
                    p = rng.permutation(len(sub)); gp = sub.copy()
                    for c in ['pmra', 'pmdec', 'pmra_error', 'pmdec_error', 'pmra_pmdec_corr']: gp[c] = sub[c].values[p]
                    Gn.append(np.linalg.norm(fit(gp)[0][:3]))
                print(f"      noise-only glide amplitude in this bin: {np.mean(Gn):.2f} ± {np.std(Gn):.2f}")
        if len(out) >= 3:
            x = np.array([o[0] for o in out]); y = np.array([o[1] for o in out]); e = np.array([o[2] for o in out])
            w = 1 / e ** 2; xm = (w * x).sum() / w.sum(); slope = (w * (x - xm) * y).sum() / (w * (x - xm) ** 2).sum(); es = 1 / np.sqrt((w * (x - xm) ** 2).sum())
            yw = np.array([o[3] for o in out]); ew = np.array([o[4] for o in out]); ww = 1 / ew ** 2; sw = (ww * (x - xm) * yw).sum() / (ww * (x - xm) ** 2).sum(); esw = 1 / np.sqrt((ww * (x - xm) ** 2).sum())
            chi = (((y - (w * y).sum() / w.sum()) / e) ** 2).sum()
            print(f"  trend of the glide's GC projection with {label}: {slope:+.3f} ± {es:.3f} per unit ({slope/es:+.1f} sigma); constant-projection chi2 {chi:.1f} for {len(out)-1} dof; "
                  f"spin trend {sw:+.3f} ± {esw:.3f} ({sw/esw:+.1f} sigma)")
    zs = g.zspec.values; has = np.isfinite(zs)
    bins('zspec', [0, 0.7, 1.1, 1.5, 1.9, 2.4, 6], 'spectroscopic z', has, nullperm=20)
    # the same sources without using the redshift: footprint control (a partial-sky fit of the matched sources, all z together)
    b, cov, red = fit(g, np.where(has)[0]); report(b, cov, red, has.sum(), 'all matched (footprint control)')
    # redshift at fixed magnitude: G 19-20.5 only, since spin depends on G
    bins('zspec', [0, 0.9, 1.4, 1.9, 6], 'spectroscopic z at G 19.5-20.5', has & (g.phot_g_mean_mag.values >= 19.5) & (g.phot_g_mean_mag.values < 20.5))
    bins('redshift_qsoc', [0, 0.7, 1.1, 1.5, 1.9, 2.4, 6], 'Gaia QSOC photometric z (all sky)')
    print("== QSOC photometric z at fixed G (all sky): the glide's GC projection and the spin per cell, then a joint fit proj = a + b z + c (G - 20)")
    cells = []
    for glo, ghi in [(18.5, 19.25), (19.25, 19.75), (19.75, 20.25), (20.25, 20.75)]:
        for zlo, zhi in [(0, 0.9), (0.9, 1.4), (1.4, 1.9), (1.9, 2.5), (2.5, 6)]:
            m = (g.phot_g_mean_mag.values >= glo) & (g.phot_g_mean_mag.values < ghi) & (g.redshift_qsoc.values >= zlo) & (g.redshift_qsoc.values < zhi)
            if m.sum() < 5000: continue
            b, cov, red = fit(g, np.where(m)[0]); G_, eG_, W, eW, P, eP = report(b, cov, red, m.sum(), f"G [{glo}, {ghi}) z [{zlo}, {zhi})")
            cells.append((g.phot_g_mean_mag.values[m].mean(), g.redshift_qsoc.values[m].mean(), P, eP, b[3], np.sqrt(cov[3, 3]), W, eW))
    C = np.array(cells); X = np.column_stack([np.ones(len(C)), C[:, 1], C[:, 0] - 20])
    for k, (name, yi, ei) in enumerate([('GC projection', 2, 3), ('spin w_x', 4, 5)]):
        y, e = C[:, yi], C[:, ei]; w = 1 / e ** 2; Xw = X * w[:, None]; cv = np.linalg.inv(X.T @ Xw); bb = cv @ (Xw.T @ y); ee = np.sqrt(np.diag(cv))
        chi = (w * (y - X @ bb) ** 2).sum()
        print(f"  joint fit of the {name} over {len(C)} cells: per unit z {bb[1]:+.3f} ± {ee[1]:.3f} ({bb[1]/ee[1]:+.1f} sigma), per mag {bb[2]:+.3f} ± {ee[2]:.3f} ({bb[2]/ee[2]:+.1f} sigma); chi2 {chi:.1f} / {len(C)-3}")
    print("== cells of G x BP-RP x QSOC z (all sky, >= 6000 sources): joint linear fit of the spin components and the GC projection on z, G and colour")
    cells = []
    for glo, ghi in [(18.0, 19.25), (19.25, 19.75), (19.75, 20.25), (20.25, 20.75)]:
        for clo, chi_ in [(-1, 0.55), (0.55, 0.75), (0.75, 1.0), (1.0, 4)]:
            for zlo, zhi in [(0, 0.9), (0.9, 1.4), (1.4, 1.9), (1.9, 2.5), (2.5, 6)]:
                m = (g.phot_g_mean_mag.values >= glo) & (g.phot_g_mean_mag.values < ghi) & (g.bp_rp.values >= clo) & (g.bp_rp.values < chi_) & (g.redshift_qsoc.values >= zlo) & (g.redshift_qsoc.values < zhi)
                if m.sum() < 6000: continue
                b, cov, red = fit(g, np.where(m)[0]); P, eP = proj(b, cov)
                cells.append((g.phot_g_mean_mag.values[m].mean(), g.bp_rp.values[m].mean(), g.redshift_qsoc.values[m].mean(), m.sum(), P, eP, b[3], np.sqrt(cov[3, 3]), b[4], np.sqrt(cov[4, 4]), b[5], np.sqrt(cov[5, 5])))
    C = np.array(cells); print(f"  {len(C)} cells")
    for name, yi, ei in [('GC projection', 4, 5), ('spin w_x', 6, 7), ('spin w_y', 8, 9), ('spin w_z', 10, 11)]:
        y, e = C[:, yi], C[:, ei]; w = 1 / e ** 2
        for label, X in [('z, G, colour', np.column_stack([np.ones(len(C)), C[:, 2], C[:, 0] - 20, C[:, 1] - 0.7])), ('z, G', np.column_stack([np.ones(len(C)), C[:, 2], C[:, 0] - 20])), ('G, colour', np.column_stack([np.ones(len(C)), C[:, 0] - 20, C[:, 1] - 0.7]))]:
            Xw = X * w[:, None]; cv = np.linalg.inv(X.T @ Xw); bb = cv @ (Xw.T @ y); ee = np.sqrt(np.diag(cv)); chi = (w * (y - X @ bb) ** 2).sum()
            terms = ' '.join(f"{t} {bb[k+1]:+.2f} ± {ee[k+1]:.2f} ({bb[k+1]/ee[k+1]:+.1f}σ)" for k, t in enumerate(label.split(', ')))
            print(f"  {name:13s} on [{label:12s}]: {terms}; chi2 {chi:.0f} / {len(C)-X.shape[1]}")
    bins('phot_g_mean_mag', [14, 18, 19, 19.5, 20, 20.5, 21.5], 'G')
    bins('bp_rp', [-1, 0.4, 0.7, 0.9, 1.1, 1.4, 4], 'BP-RP')
    bins('ruwe', [0, 0.95, 1.05, 1.15, 1.3, 5], 'RUWE')
    bins('astrometric_excess_noise', [0, 0.3, 0.6, 1.0, 1.5, 10], 'excess noise (mas)')
    bins('visibility_periods_used', [0, 14, 18, 22, 26, 60], 'visibility periods')
    g['abs_ecl'] = np.abs(g.ecl_lat.values); bins('abs_ecl', [0, 15, 30, 45, 60, 90], '|ecliptic latitude|')

if __name__ == '__main__':
    a = sys.argv[1:]; main(a[0], a[1], a[2], int(a[3]) if len(a) > 3 else 50)
