"""Is the Migkas+2020 L_X-T anisotropy (eeHIFLUGCS, 313 clusters; A&A 636, A15, VizieR J/A+A/636/A15 via vizier.idia.ac.za)
the ROSAT All-Sky Survey exposure pattern? Follows `erass1_lxt.py`, where the eRASS1 version of the signal toward Migkas's
axis is the eROSITA exposure map. RASS and eRASS share the scan geometry: great circles through the ecliptic poles, so the
exposure runs as ~1/cos(ecliptic latitude) and peaks at the poles; Migkas+2021's axis (280, -15) is 15 deg from the south
ecliptic pole (276, -30). The eeHIFLUGCS L_X are RASS count rates; T are Chandra/XMM.

Fit log L = a + B log(T/4 keV) + log E(z) (Migkas's convention, L/E(z); --gfree lets the E(z) power float) by maximum
likelihood with errors in both and an intrinsic scatter; residual r = log L - model. Then, as in erass1_lxt.py: the weighted
mean residual in the 75 deg cone about (280, -15) [and about Migkas+2020's own L-T axis (303, -27), memory] minus outside,
in dex and as H0_cone/H0_rest = 10^(-dr/2), against a shuffle null; the all-sky scan with its look-elsewhere null; the fit with
an exposure term d log(1/cos beta) (capped at |beta| = 88.5 deg); Spearman(r, exposure proxy); the cone deficit within
exposure (|beta|) bands; and the residual's own dependence on ecliptic latitude.
usage: python migkas_exposure.py data/migkas/J_A_A_636_A15.tsv [--gfree] [--nshuf 4000] [--zmin 0.0]
"""
import sys, argparse, numpy as np, healpy as hp
from scipy.optimize import minimize
from scipy.stats import spearmanr
from astropy.coordinates import SkyCoord, BarycentricMeanEcliptic
import astropy.units as u

OM = 0.3; CONE = 75.0
AXES = {'Migkas+2021 (280,-15)': (280.0, -15.0), 'Migkas+2020 L-T (303,-27) [memory]': (303.0, -27.0), 'south ecliptic pole (276.4,-29.8)': (276.4, -29.8)}

def unit(l, b):
    l, b = np.radians(l), np.radians(b); return np.stack([np.cos(b)*np.cos(l), np.cos(b)*np.sin(l), np.sin(b)], -1)

def read(tsv):
    rows = [l.rstrip('\n').split('\t') for l in open(tsv) if not l.startswith('#') and l.strip()]
    hdr = [h.strip() for h in rows[0]]; data = [r for r in rows[3:] if len(r) == len(hdr)]
    col = lambda n: np.array([r[hdr.index(n)].strip() for r in data])
    f = lambda n: np.array([float(v) if v else np.nan for v in col(n)])
    return dict(name=col('Cluster'), z=f('z'), l=f('GLON'), b=f('GLAT'), T=f('T'), Tp=f('E_T'), Tm=f('e_T'), L=f('LX'), eLpct=f('e_LX'), F=f('FX'), nh=f('NHtot'))

def fit(y, ey, x, ex, zt, gfree, extra=None):
    cols = [np.ones_like(y), x] + ([zt] if gfree else []) + ([extra] if extra is not None else []); X = np.column_stack(cols)
    yy = y if gfree else y - zt   # Migkas: L/E(z)
    def nll(p):
        beta, lns = p[:-1], p[-1]; var = ey**2 + (beta[1]*ex)**2 + np.exp(2*lns); res = yy - X@beta; return 0.5*np.sum(res**2/var + np.log(var))
    b0, *_ = np.linalg.lstsq(X, yy, rcond=None); r = minimize(nll, np.append(b0, np.log(0.2)), method='BFGS'); p = r.x
    err = np.sqrt(np.abs(np.diag(r.hess_inv))); beta, s = p[:-1], np.exp(p[-1]); var = ey**2 + (beta[1]*ex)**2 + s**2
    return beta, err[:-1], s, yy - X@beta, var

def cone_stat(res, w, nvec, axis):
    inside = nvec@axis > np.cos(np.radians(CONE)); mi = np.sum(w[inside]*res[inside])/np.sum(w[inside]); mo = np.sum(w[~inside]*res[~inside])/np.sum(w[~inside]); return mi - mo, inside.sum()

def scan(res, w, nvec, nside=8):
    th, ph = hp.pix2ang(nside, np.arange(hp.nside2npix(nside))); axes = unit(np.degrees(ph), 90 - np.degrees(th))
    inside = (nvec@axes.T) > np.cos(np.radians(CONE)); W = w[:, None]*inside; sw = W.sum(0); swo = w.sum() - sw
    with np.errstate(invalid='ignore', divide='ignore'):
        mi = (W*res[:, None]).sum(0)/sw; mo = ((w*res).sum() - (W*res[:, None]).sum(0))/swo
    d = mi - mo; d[(sw < 20*w.mean()) | (swo < 20*w.mean())] = 0; return d, axes

def analyse(label, res, var, nvec, lexp, rng, nshuf):
    w = 1/var; out = {}
    dirs, axes = scan(res, w, nvec); imin = np.argmin(dirs)
    null = {k: np.empty(nshuf) for k in AXES}; nullmax = np.empty(nshuf)
    for k in range(nshuf):
        p = rng.permutation(len(res)); r_, w_ = res[p], w[p]
        for name, ax in AXES.items(): null[name][k] = cone_stat(r_, w_, nvec, unit(*ax))[0]
        nullmax[k] = np.max(np.abs(scan(r_, w_, nvec)[0]))
    for name, ax in AXES.items():
        d, nin = cone_stat(res, w, nvec, unit(*ax)); sig = d/null[name].std(); pM = np.mean(np.abs(null[name]) >= abs(d)); h0 = 10**(-d/2)
        print(f'[{label}] cone about {name} (n_in={nin}): {d:+.4f} dex ({sig:+.2f} sigma, p={pM:.3g}); apparent H0_cone/H0_rest = {h0:.3f} ({100*(h0-1):+.1f} %)')
        out[name] = (d, sig, h0)
    lmin, bmin = np.degrees(np.arctan2(axes[imin, 1], axes[imin, 0])) % 360, np.degrees(np.arcsin(axes[imin, 2]))
    pLE = np.mean(nullmax >= abs(dirs[imin])); rho, prho = spearmanr(res, lexp)
    print(f'[{label}] lowest-residual direction (l,b)=({lmin:.0f},{bmin:.0f}): {dirs[imin]:+.4f} dex, {np.degrees(np.arccos(np.clip(axes[imin]@unit(280, -15), -1, 1))):.0f} deg from (280,-15); look-elsewhere p = {pLE:.3f}; Spearman(residual, log exposure proxy) = {rho:+.3f} (p={prho:.2g})')
    out['scan'] = (dirs[imin], lmin, bmin, pLE); out['rho'] = (rho, prho); return out

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('tsv'); ap.add_argument('--gfree', action='store_true'); ap.add_argument('--nshuf', type=int, default=4000)
    ap.add_argument('--zmin', type=float, default=0.0); a = ap.parse_args()
    d = read(a.tsv); ok = np.isfinite(d['T']) & np.isfinite(d['L']) & (d['L'] > 0) & np.isfinite(d['z']) & (d['z'] > a.zmin)
    for k in d: d[k] = d[k][ok]
    print(f'eeHIFLUGCS (Migkas+2020): {ok.sum()} clusters with T and L_X; z median {np.median(d["z"]):.3f} (5-95 %: {np.percentile(d["z"], 5):.3f}-{np.percentile(d["z"], 95):.3f}); T median {np.median(d["T"]):.2f} keV')
    c = SkyCoord(l=d['l']*u.deg, b=d['b']*u.deg, frame='galactic'); beta = c.transform_to(BarycentricMeanEcliptic()).lat.deg
    nvec = unit(d['l'], d['b']); expo = 1/np.cos(np.radians(np.clip(np.abs(beta), 0, 88.5))); lexp = np.log10(expo)
    y = np.log10(d['L']); ey = d['eLpct']/100/np.log(10); x = np.log10(d['T']/4); ex = (np.log10(d['T'] + d['Tp']) - np.log10(d['T'] - d['Tm']))/2
    Ez = np.sqrt(OM*(1 + d['z'])**3 + 1 - OM); zt = np.log10(Ez)
    inM = nvec@unit(280, -15) > np.cos(np.radians(CONE))
    print(f'ecliptic |beta|: mean inside the (280,-15) cone {np.abs(beta[inM]).mean():.1f} deg (n={inM.sum()}), outside {np.abs(beta[~inM]).mean():.1f}; mean log exposure proxy {lexp[inM].mean():+.3f} vs {lexp[~inM].mean():+.3f}')
    rng = np.random.default_rng(118)
    beta_, err, s, res, var = fit(y, ey, x, ex, zt, a.gfree)
    print(f'L-T fit: slope {beta_[1]:.2f} ± {err[1]:.2f}' + (f', E(z) power {beta_[2]:.2f} ± {err[2]:.2f}' if a.gfree else ' (L/E(z), Migkas convention)') + f'; intrinsic scatter {s:.3f} dex')
    r0 = analyse('plain', res, var, nvec, lexp, rng, a.nshuf)
    beta_, err, s, res_e, var_e = fit(y, ey, x, ex, zt, a.gfree, extra=lexp)
    print(f'L-T + exposure: slope {beta_[1]:.2f} ± {err[1]:.2f}; exposure term {beta_[-1]:+.3f} ± {err[-1]:.3f} dex per dex of 1/cos(beta); scatter {s:.3f}')
    r1 = analyse('exposure term', res_e, var_e, nvec, lexp, rng, a.nshuf)
    # the residual against ecliptic latitude directly, and the cone deficit within |beta| bands
    w = 1/var
    print('plain residual by ecliptic |beta| band (weighted mean ± error; and the (280,-15) cone deficit within the band, shuffle null within the band):')
    for lo, hi in [(0, 20), (20, 40), (40, 60), (60, 90)]:
        m = (np.abs(beta) >= lo) & (np.abs(beta) < hi)
        if m.sum() < 15: continue
        mean = np.average(res[m], weights=w[m]); e = np.sqrt(1/np.sum(w[m]))
        dq, nin = cone_stat(res[m], w[m], nvec[m], unit(280, -15)) if (nvec[m]@unit(280, -15) > np.cos(np.radians(CONE))).sum() >= 5 and (nvec[m]@unit(280, -15) <= np.cos(np.radians(CONE))).sum() >= 5 else (np.nan, 0)
        null = np.array([cone_stat(res[m][p], w[m][p], nvec[m], unit(280, -15))[0] for p in (rng.permutation(m.sum()) for _ in range(1000))]) if nin else np.array([np.nan])
        print(f'  {lo:2d}<=|beta|<{hi:2d}: n={m.sum():3d}, mean residual {mean:+.4f} ± {e:.4f}; cone deficit {dq:+.4f} ({dq/np.nanstd(null):+.2f} sigma, n_in={nin})')
    # a 2D check: residual on sin(b) and on the exposure proxy together (is the "Galactic" axis the ecliptic one?)
    X = np.column_stack([np.ones_like(res), lexp, nvec@unit(280, -15)]); Wm = np.diag(w); cov = np.linalg.inv(X.T@Wm@X); p = cov@X.T@Wm@res; e = np.sqrt(np.diag(cov))
    print(f'residual on [1, log exposure, cos(angle to (280,-15))] jointly: exposure {p[1]:+.3f} ± {e[1]:.3f}, axis {p[2]:+.4f} ± {e[2]:.4f}')
    print('summary: deficit in the (280,-15) cone, dex / sigma / apparent H0 ratio; Spearman(res, exposure):')
    for k, r in (('plain', r0), ('exposure term', r1)):
        d, sig, h0 = r['Migkas+2021 (280,-15)']; print(f'  {k:14s} {d:+.4f} {sig:+.2f} {h0:.3f}  rho={r["rho"][0]:+.3f} (p={r["rho"][1]:.2g}); lowest direction ({r["scan"][1]:.0f},{r["scan"][2]:.0f}) {r["scan"][0]:+.4f} p_LE={r["scan"][3]:.3f}')

main()
