"""Sweep-3 row 26: the local galaxy density profile in DESI DR1 BGS_BRIGHT, per Galactic cap. A volume-limited subsample
(M_r < -21, k+e corrected, complete to z = 0.2 at r < 19.5) counted in redshift shells against the comoving shell volume of
each cap's footprint (area from the clustering randoms, 2,500 per deg^2 per file); delta(z) relative to the 0.1 < z < 0.2 mean;
the north-south difference; a jackknife over sky patches; a crude linear-theory cosmic-variance expectation. The k-correction
is a colour-dependent linear approximation and the evolution Q is varied (1.0, 1.5, 2.0 mag per unit z), the dominant
systematic for the absolute profile; the cap comparison shares it. usage: python bgs_local_density.py"""
import numpy as np, healpy as hp
from astropy.io import fits
from scipy.integrate import quad
C = 299792.458; H0 = 100.0; OM = 0.31  # Mpc/h
RAN_DENS = 2500.0  # randoms per deg^2 per file

def dc(z):  # comoving distance, Mpc/h
    return C / H0 * np.array([quad(lambda x: 1 / np.sqrt(OM * (1 + x) ** 3 + 1 - OM), 0, zi)[0] for zi in np.atleast_1d(z)])

def load(cap):
    t = fits.getdata(f'data/desi_bgs/BGS_BRIGHT_{cap}_clustering.dat.fits', 1)
    nran = fits.getheader(f'data/desi_bgs/BGS_BRIGHT_{cap}_0_clustering.ran.fits', 1)['NAXIS2']  # the header's row count (the file itself was cut at ~1 GB by the proxy)
    area = nran / RAN_DENS
    ok = (t['flux_r_dered'] > 0) & (t['flux_g_dered'] > 0) & (t['Z'] > 0.005) & (t['Z'] < 0.5)
    t = t[ok]; rmag = 22.5 - 2.5 * np.log10(t['flux_r_dered']); gmag = 22.5 - 2.5 * np.log10(t['flux_g_dered'])
    z = t['Z']; w = t['WEIGHT_COMP'] * t['WEIGHT_ZFAIL'] * t['WEIGHT_SYS']
    return dict(z=z, r=rmag, gr=gmag - rmag, w=w, ra=t['RA'], dec=t['DEC'], area=area)

def absmag(d, Q):
    z = d['z']; dm = 5 * np.log10(dc(z) * (1 + z) / 0.7 * 1e6 / 10)  # h = 0.7 for the magnitude scale (cancels in ratios)
    k = z * np.clip(1.6 * d['gr'] - 0.4, 0, None)  # colour-dependent r-band k-correction, linear in z (crude, ~0.05 mag)
    return d['r'] - dm - k + Q * z  # evolution-corrected to z = 0

def profile(d, Q=1.5, Mcut=-21.0, zmax=0.2, edges=None):
    M = absmag(d, Q); m = (M < Mcut) & (d['z'] < zmax) & (d['z'] > 0.01)
    edges = np.arange(0.01, zmax + 1e-9, 0.01) if edges is None else edges
    rr = dc(edges); vol = d['area'] * (np.pi / 180) ** 2 / 3 * (rr[1:] ** 3 - rr[:-1] ** 3)  # (Mpc/h)^3 per shell
    cnt = np.array([d['w'][m & (d['z'] >= a) & (d['z'] < b)].sum() for a, b in zip(edges[:-1], edges[1:])])
    raw = np.array([(m & (d['z'] >= a) & (d['z'] < b)).sum() for a, b in zip(edges[:-1], edges[1:])])
    n = cnt / vol; ref = (edges[:-1] >= 0.1); nbar = cnt[ref].sum() / vol[ref].sum()
    return edges, n / nbar - 1, raw, vol, m, M

def main():
    caps = {c: load(c) for c in ['NGC', 'SGC']}
    for c in caps: print(f"{c}: {len(caps[c]['z'])} galaxies, area {caps[c]['area']:.0f} deg^2")
    # the completeness of the M_r < -21 cut at z = 0.2: the flux limit in absolute magnitude for a red galaxy
    for Q in [1.0, 1.5, 2.0]:
        z = 0.2; dm = 5 * np.log10(dc(z)[0] * 1.2 / 0.7 * 1e5); print(f"Q={Q}: M_lim(z=0.2, red) = {19.5 - dm - 0.2*1.04 + Q*0.2:.2f}")
    for Q in [1.5, 1.0, 2.0]:
        print(f"\n=== Q = {Q} mag per unit z (k+e), M_r < -21, shells dz = 0.01 ===")
        res = {}
        for c in ['NGC', 'SGC']:
            edges, delta, raw, vol, m, M = profile(caps[c], Q); res[c] = (edges, delta, raw)
        zc = 0.5 * (res['NGC'][0][:-1] + res['NGC'][0][1:]); d_h = dc(zc)
        print("  z     d(Mpc/h)   N_NGC   delta_NGC   N_SGC   delta_SGC   diff")
        for i in range(len(zc)):
            print(f"  {zc[i]:.3f}  {d_h[i]:6.0f}  {res['NGC'][2][i]:6d}  {res['NGC'][1][i]:+.3f}  {res['SGC'][2][i]:6d}  {res['SGC'][1][i]:+.3f}  {res['NGC'][1][i]-res['SGC'][1][i]:+.3f}")
        # integrated contrasts: within 150 Mpc/h (z < 0.05) and within 300 Mpc/h (z < 0.1), vs the 0.1-0.2 reference
        for c in ['NGC', 'SGC']:
            d = caps[c]; M = absmag(d, Q); m = (M < -21) & (d['z'] > 0.01) & (d['z'] < 0.2)
            for zlo, zhi, lab in [(0.01, 0.05, '<150 Mpc/h'), (0.05, 0.1, '150-300'), (0.01, 0.1, '<300')]:
                e = np.array([zlo, zhi, 0.1, 0.2]); rr = dc(e); vol = d['area'] * (np.pi / 180) ** 2 / 3 * (rr[1:] ** 3 - rr[:-1] ** 3)
                n1 = d['w'][m & (d['z'] >= zlo) & (d['z'] < zhi)].sum() / vol[0]; nref = d['w'][m & (d['z'] >= 0.1)].sum() / vol[2]
                N1 = (m & (d['z'] >= zlo) & (d['z'] < zhi)).sum()
                # jackknife over nside-4 patches
                pix = hp.ang2pix(4, d['ra'], d['dec'], lonlat=True); pats = np.unique(pix[m]); jk = []
                for p in pats:
                    k = m & (pix != p); f = 1 - (pix[m] == p).mean()  # area fraction removed ~ galaxy fraction
                    jk.append((d['w'][k & (d['z'] >= zlo) & (d['z'] < zhi)].sum() / (vol[0] * f)) / (d['w'][k & (d['z'] >= 0.1)].sum() / (vol[2] * f)) - 1)
                jk = np.array(jk); ejk = np.sqrt((len(jk) - 1) / len(jk) * np.sum((jk - jk.mean()) ** 2))
                Reff = (3 * vol[0] / (4 * np.pi)) ** (1 / 3); cv = 0.8 * 1.3 * (Reff / 8) ** -1.0  # crude: sigma_8 x bias x (R/8)^-1
                print(f"  {c} {lab}: delta = {n1/nref-1:+.3f} ± {ejk:.3f} (jackknife, {len(jk)} patches; Poisson {1/np.sqrt(N1):.3f}); N = {N1}; "
                      f"R_eff = {Reff:.0f} Mpc/h, crude LCDM sigma(delta) ~ {cv:.3f}")
    # bright cut M_r < -21.5 and the colour split at Q = 1.5
    print("\n=== variants at Q = 1.5: M_r < -21.5; red (g-r > 0.7) and blue, M_r < -21 ===")
    for c in ['NGC', 'SGC']:
        d = caps[c]; M = absmag(d, 1.5)
        for lab, m in [('M<-21.5', M < -21.5), ('red', (M < -21) & (d['gr'] > 0.7)), ('blue', (M < -21) & (d['gr'] <= 0.7))]:
            mm = m & (d['z'] > 0.01) & (d['z'] < 0.2); e = np.array([0.01, 0.05, 0.1, 0.2]); rr = dc(e)
            vol = d['area'] * (np.pi / 180) ** 2 / 3 * (rr[1:] ** 3 - rr[:-1] ** 3)
            n = [d['w'][mm & (d['z'] >= e[i]) & (d['z'] < e[i + 1])].sum() / vol[i] for i in range(3)]
            print(f"  {c} {lab}: delta(<150) = {n[0]/n[2]-1:+.3f}, delta(150-300) = {n[1]/n[2]-1:+.3f}  (N {(mm & (d['z']<0.05)).sum()}, {(mm & (d['z']>=0.05) & (d['z']<0.1)).sum()})")

if __name__ == '__main__': main()
