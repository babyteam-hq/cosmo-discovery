"""Planck PR3 CMB-lensing convergence: hemispherical power asymmetry and dipolar modulation (sweep-2.md row 18).
Data: COM_Lensing_4096_R3.00 from the Planck Legacy Archive (MV/dat_klm.fits, the mean-field-subtracted kappa alm;
MV/nlkk.dat with N_L and the fiducial C_L^kk + N_L; mask.fits.gz, f_sky 0.67). PR3 ships no simulations, so the null is
Gaussian isotropic realizations of the fiducial total spectrum (C_L^kk + N_L) through the same mask and statistics. That null
has isotropic noise; the real reconstruction noise follows Planck's scan depth (lowest at the ecliptic poles), so bands where
N_L > C_L^kk (L > ~100 for PR3 MV) test the noise pattern, not cosmology. The ecliptic-pole axis is run as that diagnostic.

Statistics on the band-filtered kappa map (nside 256): for an axis n, A(n) = (P_N - P_S)/(P_N + P_S) with P the masked
variance in each hemisphere, for the pre-chosen axes (CMB hemispherical-asymmetry axis (227, -27), CMB dipole (264, 48),
Galactic and ecliptic poles) and over a scan of nside-4 directions (look-elsewhere from the sims' max |A|); and a dipole fit
to the local-variance map (nside-8 cells) with its amplitude and direction against the sims.
usage: python pr3_kappa_asym.py data/pr3/COM_Lensing_4096_R3.00 [--nsim 300] [--nside 256]
"""
import sys, argparse, time, numpy as np, healpy as hp
from astropy.io import fits

AXES = {'CMB power asymmetry (227,-27)': (227.0, -27.0), 'CMB dipole (264,48)': (264.021, 48.253), 'Galactic pole (0,90)': (0.0, 90.0),
        'ecliptic pole (96.4,29.8)': (96.38, 29.81)}
BANDS = [(8, 40), (40, 100), (100, 400), (8, 400)]

def unit(l, b):
    l, b = np.radians(l), np.radians(b); return np.array([np.cos(b)*np.cos(l), np.cos(b)*np.sin(l), np.sin(b)])

def read_klm(path, lmax):
    d = fits.open(path)[1].data; idx = np.asarray(d['index'], int); re, im = np.asarray(d['real']), np.asarray(d['imag'])
    # Planck index = l^2 + l + m + 1
    l = np.floor(np.sqrt(idx - 1)).astype(int); m = idx - 1 - l*l - l
    k = (l <= lmax) & (m <= l); alm = np.zeros(hp.Alm.getsize(lmax), complex); alm[hp.Alm.getidx(lmax, l[k], m[k])] = re[k] + 1j*im[k]
    return alm

def band_maps(alm, lmax, nside, bands):
    out = []
    for lo, hi in bands:
        fl = np.zeros(lmax + 1); fl[lo:hi + 1] = 1.0; out.append(hp.alm2map(hp.almxfl(alm, fl), nside, verbose=False) if 'verbose' in hp.alm2map.__code__.co_varnames else hp.alm2map(hp.almxfl(alm, fl), nside))
    return out

def hemi_stats(maps, mask, pix_vec, axes, scan_axes, cell_nside=8):
    """per band: A for each axis, max |A| over the scan, and the local-variance dipole (amplitude, direction)."""
    res = []
    for m in maps:
        v = m*m; good = mask > 0.5
        A = {}
        for name, ax in axes.items():
            north = pix_vec@ax > 0; pn = v[good & north].mean(); ps = v[good & ~north].mean(); A[name] = (pn - ps)/(pn + ps)
        proj = pix_vec@scan_axes.T  # (npix, nscan)
        sc = np.empty(scan_axes.shape[0])
        for j in range(scan_axes.shape[0]):
            north = proj[:, j] > 0; pn = v[good & north].mean(); ps = v[good & ~north].mean(); sc[j] = (pn - ps)/(pn + ps)
        # local variance dipole: nside-8 cells, mean of v over good pixels in the cell, cells with > 50 % coverage
        cell = hp.ud_grade(np.where(good, v, 0.0), cell_nside, power=0)*1.0; cov = hp.ud_grade(good.astype(float), cell_nside)
        ok = cov > 0.5; lv = cell[ok]/cov[ok]; cv = np.array(hp.pix2vec(cell_nside, np.arange(hp.nside2npix(cell_nside)))).T[ok]
        X = np.column_stack([np.ones(ok.sum()), cv]); p, *_ = np.linalg.lstsq(X, lv/lv.mean() - 1, rcond=None)
        dip = p[1:]; amp = np.linalg.norm(dip); ld, bd = np.degrees(np.arctan2(dip[1], dip[0])) % 360, np.degrees(np.arcsin(dip[2]/amp))
        res.append((A, sc, amp, ld, bd))
    return res

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('dir'); ap.add_argument('--nsim', type=int, default=300); ap.add_argument('--nside', type=int, default=256)
    ap.add_argument('--seed', type=int, default=118); a = ap.parse_args()
    nside = a.nside; lmax = 2*nside; rng = np.random.default_rng(a.seed)
    t0 = time.time(); alm = read_klm(f'{a.dir}/MV/dat_klm.fits', lmax); nl = np.loadtxt(f'{a.dir}/MV/nlkk.dat'); NL, CLtot = nl[:lmax + 1, 1], nl[:lmax + 1, 2]
    mask = hp.ud_grade(hp.read_map(f'{a.dir}/mask.fits.gz'), nside); mask = (mask > 0.9).astype(float); fsky = mask.mean()
    print(f'PR3 MV kappa alm to lmax {lmax}; mask f_sky {fsky:.3f} at nside {nside}; read in {time.time()-t0:.0f} s')
    # sanity: the masked pseudo-spectrum / f_sky against the fiducial total
    kmap = hp.alm2map(alm, nside); cl = hp.anafast(kmap*mask, lmax=lmax)/fsky
    for lo, hi in BANDS: print(f'  band {lo}-{hi}: data pseudo-C_L/f_sky over fiducial (C_L^kk+N_L) = {np.mean(cl[lo:hi+1])/np.mean(CLtot[lo:hi+1]):.3f}; N_L/(C_L+N_L) = {np.mean(NL[lo:hi+1]/CLtot[lo:hi+1]):.2f}')
    pix_vec = np.array(hp.pix2vec(nside, np.arange(hp.nside2npix(nside)))).T
    # everything in Galactic coordinates (the Planck maps are Galactic); axes given in (l, b)
    axes = {k: unit(*v) for k, v in AXES.items()}
    sn = 4; th, ph = hp.pix2ang(sn, np.arange(hp.nside2npix(sn))); scan_axes = np.array(hp.ang2vec(th, ph)); scan_axes = scan_axes[scan_axes[:, 2] >= 0]  # half the sphere: A(-n) = -A(n)
    data = hemi_stats(band_maps(alm, lmax, nside, BANDS), mask, pix_vec, axes, scan_axes)
    # the Gaussian isotropic null
    sims = []
    t0 = time.time()
    for i in range(a.nsim):
        salm = hp.synalm(CLtot, lmax=lmax, new=True) if 'new' in hp.synalm.__code__.co_varnames else hp.synalm(CLtot, lmax=lmax)
        sims.append(hemi_stats(band_maps(salm, lmax, nside, BANDS), mask, pix_vec, axes, scan_axes))
        if i == 0: print(f'  one sim takes {time.time()-t0:.1f} s')
    print(f'{a.nsim} sims in {time.time()-t0:.0f} s')
    for ib, (lo, hi) in enumerate(BANDS):
        A, sc, amp, ld, bd = data[ib]
        print(f'band L={lo}-{hi} (N_L share {np.mean(NL[lo:hi+1]/CLtot[lo:hi+1]):.2f}):')
        for name in AXES:
            null = np.array([s[ib][0][name] for s in sims]); z = (A[name] - null.mean())/null.std(); p2 = np.mean(np.abs(null - null.mean()) >= abs(A[name] - null.mean()))
            print(f'  A({name}) = {A[name]:+.4f}; sims {null.mean():+.4f} ± {null.std():.4f}; {z:+.2f} sigma, two-sided p = {p2:.3f}')
        j = np.argmax(np.abs(sc)); l_, b_ = np.degrees(np.arctan2(scan_axes[j, 1], scan_axes[j, 0])) % 360, np.degrees(np.arcsin(scan_axes[j, 2]))
        nullmax = np.array([np.max(np.abs(s[ib][1])) for s in sims]); pLE = np.mean(nullmax >= abs(sc[j]))
        namp = np.array([s[ib][2] for s in sims]); pamp = np.mean(namp >= amp)
        print(f'  scan: max |A| = {abs(sc[j]):.4f} toward (l,b)=({l_:.0f},{b_:.0f}) [sign {np.sign(sc[j]):+.0f}], sims max |A| {nullmax.mean():.4f} ± {nullmax.std():.4f}, look-elsewhere p = {pLE:.3f}')
        print(f'  local-variance dipole: amplitude {amp:.4f} toward ({ld:.0f},{bd:.0f}); sims {namp.mean():.4f} ± {namp.std():.4f}, p = {pamp:.3f}')

main()
