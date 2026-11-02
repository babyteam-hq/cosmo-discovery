"""Sweep-3 row 25: the kinematic dipole in the mean observed redshift of Quaia (G < 20.0, VizieR J/ApJ/964/69), with the
sample's depth across the sky marginalised. Prediction: dz = (1 + z) (v/c) cos(theta) with v = 369.8 km/s toward
(l, b) = (264.0, 48.25), i.e. a dipole of amplitude (1 + <z>) x 1.233e-3 in the mean z, no free parameters. Singal 2025
(arXiv:2508.20769) fit the raw dipole and found 4.6x that toward the Galactic centre; here the pixel mean z is regressed on a
dipole plus depth templates (the pixel's own source density, SFD E(B-V), csc|b|, |ecliptic latitude|; the signed ecliptic latitude is itself a dipole and is left out), with the (1+z) scaling
test across redshift bins and a template-rotation null. usage: python quaia_zdipole.py"""
import io, numpy as np, pandas as pd, healpy as hp
from cf4_io import sfd_ebv
rng = np.random.default_rng(7)
BETA = 369.82 / 299792.458; APEX = (264.021, 48.253)
NSIDE = 32

def load():
    raw = [l for l in open('data/quaia/quaia_g200.tsv').read().split('\n') if not l.startswith('#')]
    i = [k for k, l in enumerate(raw) if l.startswith('GLON')][0]
    d = pd.read_csv(io.StringIO('\n'.join([raw[i]] + [l for l in raw[i + 3:] if l.strip()])), sep='\t', low_memory=False)
    d.columns = [c.strip() for c in d.columns]
    for c in d.columns: d[c] = pd.to_numeric(d[c], errors='coerce')
    return d.dropna(subset=['GLON', 'GLAT', 'zQuaia'])

def unit_gal(l, b):
    l, b = np.radians(l), np.radians(b); return np.stack([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], -1)

def pixel_table(d, bcut, zlo=0, zhi=10, zerr_max=None):
    s = d[(np.abs(d.GLAT) > bcut) & (d.zQuaia > zlo) & (d.zQuaia < zhi)]
    if zerr_max is not None: s = s[s.e_zQuaia < zerr_max]
    pix = hp.ang2pix(NSIDE, s.GLON.values, s.GLAT.values, lonlat=True)
    g = pd.DataFrame({'pix': pix, 'z': s.zQuaia.values}).groupby('pix')
    t = pd.DataFrame({'n': g.size(), 'zmean': g['z'].mean(), 'zmed': g['z'].median(), 'zstd': g['z'].std()})
    t = t[t.n >= 20].copy(); t['ez'] = t.zstd / np.sqrt(t.n)
    l, b = hp.pix2ang(NSIDE, t.index.values, lonlat=True); t['l'], t['b'] = l, b
    # masks: Magellanic Clouds (LMC (280.5, -32.9) r 8 deg; SMC (302.8, -44.3) r 5 deg)
    n = unit_gal(l, b)
    for (lc, bc, r) in [(280.5, -32.9, 8), (302.8, -44.3, 5)]:
        t = t[np.degrees(np.arccos(np.clip(n @ unit_gal(np.array([lc]), np.array([bc]))[0], -1, 1))) > r]; n = unit_gal(t.l.values, t.b.values)
    t['nx'], t['ny'], t['nz'] = n.T
    t['sfd'] = sfd_ebv(t.l.values, t.b.values); t['cscb'] = 1 / np.sin(np.radians(np.abs(t.b.values)))
    # ecliptic latitude of the pixel
    from astropy.coordinates import SkyCoord
    c = SkyCoord(t.l.values, t.b.values, unit='deg', frame='galactic').barycentricmeanecliptic
    t['ecl'] = np.sin(np.radians(c.lat.deg)); t['absecl'] = np.abs(t.ecl)
    # depth proxy: the pixel's own density relative to the mean (clustering adds noise, the selection function the trend)
    t['dens'] = np.log(t.n / t.n.mean())
    return t

def fit(t, templates, ycol='zmean', weights=True):
    cols = [np.ones(len(t)), t.nx.values, t.ny.values, t.nz.values] + [t[c].values for c in templates]
    X = np.column_stack(cols); w = 1 / t.ez.values ** 2 if weights else np.ones(len(t))
    Xw = X * w[:, None]; A = X.T @ Xw; b = np.linalg.solve(A, Xw.T @ t[ycol].values); r = t[ycol].values - X @ b
    cov = np.linalg.inv(A) * (np.sum(w * r * r) / (len(t) - X.shape[1]))
    D = b[1:4]; eD = np.sqrt(np.diag(cov)[1:4]); amp = np.linalg.norm(D)
    l = np.degrees(np.arctan2(D[1], D[0])) % 360; bb = np.degrees(np.arcsin(D[2] / amp))
    apex = unit_gal(np.array([APEX[0]]), np.array([APEX[1]]))[0]; proj = D @ apex; eproj = np.sqrt(apex @ cov[1:4, 1:4] @ apex)
    return dict(D=D, eD=eD, amp=amp, l=l, b=bb, proj=proj, eproj=eproj, tmpl=dict(zip(templates, zip(b[4:], np.sqrt(np.diag(cov)[4:])))))

def report(lab, t, templates, ycol='zmean'):
    r = fit(t, templates, ycol); zbar = np.average(t[ycol], weights=t.n); pred = (1 + zbar) * BETA
    ang = np.degrees(np.arccos(np.clip(r['D'] @ unit_gal(np.array([APEX[0]]), np.array([APEX[1]]))[0] / r['amp'], -1, 1)))
    print(f"{lab}: N_pix={len(t)} <z>={zbar:.3f} pred {pred:.4f} | dipole {r['amp']:.4f} toward ({r['l']:.0f},{r['b']:.0f}), {ang:.0f} deg from the apex; "
          f"components ±{np.mean(r['eD']):.4f}; projection on the apex {r['proj']:+.4f} ± {r['eproj']:.4f} ({r['proj']/pred:.2f} x predicted)"
          + (("; templates " + ", ".join(f"{k} {v[0]:+.4f}±{v[1]:.4f}" for k, v in r['tmpl'].items())) if templates else ''))
    return r

def main():
    d = load(); print(f"Quaia G<20.0: {len(d)} sources; <z> {d.zQuaia.mean():.3f}; z_err median {d.e_zQuaia.median():.3f}")
    for bcut in [30, 20]:
        t = pixel_table(d, bcut)
        print(f"\n=== |b| > {bcut}: {len(t)} nside-{NSIDE} pixels, {t.n.sum()} sources ===")
        report('raw dipole (no templates)', t, [])
        report('raw dipole, median z', t, [], 'zmed')
        report('+ density', t, ['dens'])
        report('+ density, SFD, csc|b|', t, ['dens', 'sfd', 'cscb'])
        r = report('+ density, SFD, csc|b|, |ecl|', t, ['dens', 'sfd', 'cscb', 'absecl'])
        report('same, median z', t, ['dens', 'sfd', 'cscb', 'absecl'], 'zmed')
        # what the templates themselves look like: the dipole of the density and of mean z explained by density alone
        rd = fit(t, [], 'dens'); print(f"  density dipole: {rd['amp']:.3f} toward ({rd['l']:.0f},{rd['b']:.0f})")
        # (1+z) scaling: the apex projection in redshift bins, full template set
        print("  by redshift bin (apex projection, predicted (1+<z>) beta):")
        for zlo, zhi in [(0, 0.8), (0.8, 1.3), (1.3, 1.8), (1.8, 2.5), (2.5, 5)]:
            tb = pixel_table(d, bcut, zlo, zhi); rb = fit(tb, ['dens', 'sfd', 'cscb', 'absecl']); zb = np.average(tb.zmean, weights=tb.n)
            print(f"    {zlo}-{zhi}: <z>={zb:.2f} proj {rb['proj']:+.4f} ± {rb['eproj']:.4f} (pred {(1+zb)*BETA:.4f}); free dipole {rb['amp']:.4f} toward ({rb['l']:.0f},{rb['b']:.0f})")
        # the Galactic-centre direction: projection on (0, 0) with and without templates
        gc = unit_gal(np.array([0.0]), np.array([0.0]))[0]
        for lab, tm in [('raw', []), ('templates', ['dens', 'sfd', 'cscb', 'absecl'])]:
            rr = fit(t, tm); print(f"  projection on the Galactic centre ({lab}): {rr['D'] @ gc:+.4f}")
        # null: the templates resampled at rotated positions (longitude rotation + reflection), 300 times: spread of the apex projection
        proj = []
        for _ in range(300):
            t2 = t.copy(); ang = rng.uniform(0, 360); refl = rng.choice([-1, 1])
            l2 = (t.l.values + ang) % 360; b2 = t.b.values * refl
            t2['sfd'] = sfd_ebv(l2, b2)
            t2['cscb'] = t.cscb.values  # |b| unchanged by the rotation
            proj.append(fit(t2, ['dens', 'sfd', 'cscb', 'absecl'])['proj'])
        print(f"  rotation null (SFD rotated; density, |b|, ecliptic kept): apex projection spread {np.std(proj):.4f}, mean {np.mean(proj):+.4f}")
        # a stricter depth control: pixels within a narrow band of density (|dens| < 0.1)
        tb = t[np.abs(t.dens) < 0.1]; report('pixels with density within 10 % of the mean, templates', tb, ['dens', 'sfd', 'cscb', 'absecl'])
        # low-uncertainty redshifts only
        tq = pixel_table(d, bcut, zerr_max=0.1); report('z_err < 0.1 only, templates', tq, ['dens', 'sfd', 'cscb', 'absecl'])

if __name__ == '__main__': main()
