"""Step 0 of R-139: counts, footprint overlaps and the Gaussian forecast of the kappa x quasar cross power at 3.5 < z < 5.

usage: python -I qso_z4_feas.py   (data under data/, gitignored; writes data/qso_z4/feas.npz and prints the table)
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
DATA = os.path.join(HERE, 'data')
OUT = os.path.join(DATA, 'qso_z4'); os.makedirs(OUT, exist_ok=True)
NS = 64


def lowres_masks():
    import healpy as hp
    fn = os.path.join(OUT, 'lens_masks_ns64.npz')
    if os.path.exists(fn):
        g = np.load(fn); return g['act'], g['planck']
    m = hp.read_map(os.path.join(DATA, 'act_dr6', 'mask_act_dr6_lensing_v1_healpix_nside_4096_baseline.fits'), dtype=np.float32)
    act = hp.ud_grade(m, NS); del m
    pm = hp.read_map(os.path.join(DATA, 'planck', 'COM_Lensing_4096_R3.00', 'mask.fits.gz'))
    pm = hp.ud_grade(pm.astype(np.float32), NS)
    th, ph = hp.pix2ang(NS, np.arange(hp.nside2npix(NS)))
    thg, phg = hp.Rotator(coord=['C', 'G'])(th, ph)
    planck = pm[hp.ang2pix(NS, thg, phg)]
    np.savez(fn, act=act, planck=planck)
    return act, planck


def main():
    import healpy as hp, fitsio
    import qso_z4_theory as T
    act, planck = lowres_masks()
    npix = hp.nside2npix(NS); pixarea = hp.nside2pixarea(NS, degrees=True)
    th, ph = hp.pix2ang(NS, np.arange(npix)); _, bgal = hp.Rotator(coord=['C', 'G'])(th, ph)
    gal = np.abs(90 - np.degrees(_)) > 20 if False else None
    thg, _phg = hp.Rotator(coord=['C', 'G'])(th, ph)
    hi_b = np.abs(90 - np.degrees(thg)) > 20
    d16 = fitsio.read(os.path.join(DATA, 'qso', 'DR16Q_v4.fits'), columns=['RA', 'DEC', 'Z'])
    dd = np.load(os.path.join(DATA, 'qso', 'desi_full_slim.npy'))
    zd = dd['Z_not4clus']; good = np.isfinite(zd) & (zd > 0) & (zd < 10)
    samples = {}
    for name, ra, dec, z, allra, alldec in [
            ('DR16Q', d16['RA'], d16['DEC'], d16['Z'], d16['RA'], d16['DEC']),
            ('DESI DR1', dd['RA'][good], dd['DEC'][good], zd[good], dd['RA'], dd['DEC'])]:
        # footprint: nside-64 pixels holding at least 30 % of the mean density of all the survey's targets/quasars
        cnt = np.bincount(hp.ang2pix(NS, allra, alldec, lonlat=True), minlength=npix)
        mean = cnt[cnt > 0].mean()
        foot = (cnt > 0.3 * mean) & hi_b
        s = (z > 3.5) & (z < 5)
        pix = hp.ang2pix(NS, ra[s], dec[s], lonlat=True)
        inf = foot[pix]
        samples[name] = dict(foot=foot, z=z[s][inf], pix=pix[inf])
        print(f'{name}: {s.sum()} at 3.5<z<5, {inf.sum()} in footprint of {foot.sum() * pixarea:.0f} deg2 (|b|>20)')
    # overlap of the two surveys' high-z samples (positional match at 2")
    from scipy.spatial import cKDTree
    def xyz(ra, dec):
        r, d = np.radians(ra), np.radians(dec); return np.c_[np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)]
    s16 = (d16['Z'] > 3.5) & (d16['Z'] < 5); sd = good & (zd > 3.5) & (zd < 5)
    t = cKDTree(xyz(dd['RA'][sd], dd['DEC'][sd]))
    dist, _i = t.query(xyz(d16['RA'][s16], d16['DEC'][s16]))
    print(f'DR16Q high-z with a DESI high-z match within 2": {(dist < np.radians(2 / 3600)).sum()}')
    th_ = T.Theory(lmax=2000)
    nk_act = np.loadtxt(os.path.join(DATA, 'act_dr6', 'N_L_kk_act_dr6_lensing_v1_baseline.txt'))[:, 1]
    nk_pl = np.loadtxt(os.path.join(DATA, 'planck', 'COM_Lensing_4096_R3.00', 'MV', 'nlkk.dat'))
    l_pl = nk_pl[:, 0].astype(int); nkp = np.zeros(4097); nkp[l_pl] = nk_pl[:, 1]   # column 1 = N_L (2 = N_L + C_L)
    rows = []
    for name, S in samples.items():
        zc = np.linspace(3.5, 5.0, 31)
        h, _e = np.histogram(S['z'], bins=np.linspace(3.45, 5.05, 33)); zc = 0.5 * (_e[1:] + _e[:-1])
        sp = th_.spectra(zc, h.astype(float))
        np.savez(os.path.join(OUT, f'spectra_{name.split()[0]}.npz'), zc=zc, nz=h, **sp)
        for lens, lm, nk, lmin, lmax in [('ACT DR6', act > 0.5, nk_act, 40, 1000), ('Planck PR3', planck > 0.5, nkp, 8, 400)]:
            ov = S['foot'] & lm
            fsky = ov.mean(); n_ov = ov[S['pix']].sum()
            nbar = n_ov / (ov.sum() * hp.nside2pixarea(NS))
            for s_mag in (0.4, 0.8):
                sn = {b: T.forecast(sp, b, s_mag, nbar, fsky, nk, lmin, lmax) for b in (3, 6, 13)}
                rows.append((name, lens, n_ov, fsky * 41253, s_mag, sn))
                print(f'{name:9s} x {lens:10s}: N {n_ov:6d} on {fsky * 41253:6.0f} deg2; s={s_mag}: S/N b=3 {sn[3]:.1f}, '
                      f'b=6 {sn[6]:.1f}, b=13 {sn[13]:.1f}')
    # the model difference: chi^2 separation between b=13 and b=4 (Gaussian), per row with s = 0.8
    print('done')


if __name__ == '__main__':
    main()
