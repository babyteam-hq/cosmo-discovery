"""R-156: lognormal light-cone mocks of DESI DR1 BGS_BRIGHT with linear velocities, for injection-recovery.

One Gaussian linear field G (z = 0, Planck 2018 CAMB P(k)) on a 512^3 grid, box 3000 Mpc/h, observer at the centre.
Galaxies: Poisson draws per cell with mean nbar_obs(r) x mask(angle) x V_cell x exp(b D G - (b D sigma_G)^2 / 2),
b(z) = 1.34/D(z); placed uniformly in the cell, then moved to redshift space by s = r + f D u_r (u = grad^-2 grad G,
the linear displacement, so v/c = f H D u). Because the velocity is linear in G and the galaxies are a lognormal
function of G, <delta_g v> is exactly b D <G v> (Stein's lemma): the mock's Doppler dipole is the linear template.
Randoms: the same draws with b = 0, 2.35x the galaxies (matching the data run's random fraction).
Marks written per galaxy: KD = (1/(rH) - 1) f H D u_r (the true Doppler convergence at the galaxy),
KB = D G at its cell (the matter contrast: kappa = b_s KB is the intrinsic size-density term with unit b_s).

The angular mask is the HEALPix (nside 128) density of the DESI randoms, normalised to its median over occupied
pixels and capped at 1; nbar_obs(z) per cap is the data's observed (unweighted) count per redshift slice over the
masked volume.

usage: python -u r156_mock.py SEED
"""
import sys, time
import numpy as np
import healpy as hp
import fitsio
import r156_theory as T

D = 'data/r156/'
NG = 512
L = 3000.0
NSIDE = 128
RANX = 2.35


def main():
    seed = int(sys.argv[1])
    t0 = time.time()
    rng = np.random.default_rng(seed)
    res, kh, pk0, om = T.cosmo()
    cell = L / NG
    # Gaussian field
    kx = np.fft.fftfreq(NG, d=cell) * 2 * np.pi
    kz = np.fft.rfftfreq(NG, d=cell) * 2 * np.pi
    k2 = kx[:, None, None] ** 2 + kx[None, :, None] ** 2 + kz[None, None, :] ** 2
    kk = np.sqrt(k2)
    Pk = np.interp(kk, kh, pk0, left=0, right=0).astype(np.float32)
    Gk = np.fft.rfftn(rng.standard_normal((NG, NG, NG), dtype=np.float32))
    Gk *= np.sqrt(Pk / cell ** 3).astype(np.float32)
    del Pk
    Gk[0, 0, 0] = 0
    G = np.fft.irfftn(Gk, s=(NG, NG, NG)).astype(np.float32)
    sigG = float(G.std())
    print('field: sigma_G %.3f (cell %.2f Mpc/h), %.0fs' % (sigG, cell, time.time() - t0), flush=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        inv = np.where(k2 > 0, 1.0 / k2, 0).astype(np.float32)
    U = []
    for ax in range(3):
        kv = [kx[:, None, None], kx[None, :, None], kz[None, None, :]][ax]
        U.append(np.fft.irfftn(1j * kv * inv * Gk, s=(NG, NG, NG)).astype(np.float32))
    del Gk, inv, k2, kk
    print('displacements done %.0fs' % (time.time() - t0), flush=True)

    # mask and nbar
    zg = np.linspace(0, 0.7, 1401)
    bg = [T.background(z) for z in zg]
    rg = np.array([b[0] for b in bg]); Hg = np.array([b[1] for b in bg])
    fg = np.array([b[2] for b in bg]); Dg = np.array([b[3] for b in bg])
    npix = hp.nside2npix(NSIDE)
    mask = np.zeros(npix)
    capmap = np.full(npix, -1)
    zbins = np.arange(0.08, 0.5601, 0.01)
    nbar = {}
    for cap, cname in enumerate(['NGC', 'SGC']):
        ran = fitsio.read(D + 'BGS_BRIGHT_%s_0_clustering.ran.fits' % cname, columns=['RA', 'DEC'])
        pix = hp.ang2pix(NSIDE, ran['RA'], ran['DEC'], lonlat=True)
        c = np.bincount(pix, minlength=npix).astype(float)
        occ = c > 0
        mask[occ] = c[occ]
        capmap[occ] = cap
    med = np.median(mask[mask > 0])
    mask = np.clip(mask / med, 0, 1)
    pixarea = 4 * np.pi / npix
    for cap, cname in enumerate(['NGC', 'SGC']):
        dat = fitsio.read(D + 'BGS_BRIGHT_%s_clustering.dat.fits' % cname, columns=['Z'])
        h, _ = np.histogram(dat['Z'], zbins)
        omega = (mask * (capmap == cap)).sum() * pixarea
        r_e = np.interp(zbins, zg, rg)
        vol = omega / 3 * (r_e[1:] ** 3 - r_e[:-1] ** 3)
        nbar[cap] = h / vol
    # cells
    ax1 = (np.arange(NG) + 0.5) * cell - L / 2
    rmin, rmax = np.interp(0.08, zg, rg) - 30, np.interp(0.56, zg, rg) + 30
    gal = {k: [] for k in ['RA', 'DEC', 'Z', 'ZTRUE', 'KD', 'KB']}
    ran = {k: [] for k in ['RA', 'DEC', 'Z']}
    for i in range(NG):
        X = ax1[i]
        Y, Zc = np.meshgrid(ax1, ax1, indexing='ij')
        R = np.sqrt(X ** 2 + Y ** 2 + Zc ** 2)
        m = (R > rmin) & (R < rmax)
        if not m.any():
            continue
        y, zc, r = Y[m], Zc[m], R[m]
        pix = hp.vec2pix(NSIDE, np.full_like(y, X), y, zc)
        mk = mask[pix]
        cp = capmap[pix]
        ok = mk > 0
        y, zc, r, pix, mk, cp = y[ok], zc[ok], r[ok], pix[ok], mk[ok], cp[ok]
        zt = np.interp(r, rg, zg)
        izb = np.clip(np.searchsorted(zbins, zt) - 1, 0, len(zbins) - 2)
        nb = np.where(cp == 0, nbar[0][izb], nbar[1][izb])
        nb[(zt < 0.08) | (zt > 0.56)] = 0
        Dz = np.interp(zt, zg, Dg); fz = np.interp(zt, zg, fg); Hz = np.interp(zt, zg, Hg)
        b = 1.34 / Dz
        g = G[i][m][ok]
        rho = np.exp(b * Dz * g - 0.5 * (b * Dz * sigG) ** 2)
        lam = nb * mk * cell ** 3
        for kind, mean in (('gal', lam * rho), ('ran', lam * RANX)):
            n = rng.poisson(mean)
            sel = np.repeat(np.arange(len(n)), n)
            if not len(sel):
                continue
            px = X + (rng.random(len(sel)) - 0.5) * cell
            py = y[sel] + (rng.random(len(sel)) - 0.5) * cell
            pz = zc[sel] + (rng.random(len(sel)) - 0.5) * cell
            rr = np.sqrt(px ** 2 + py ** 2 + pz ** 2)
            if kind == 'gal':
                ur = (U[0][i][m][ok][sel] * px + U[1][i][m][ok][sel] * py + U[2][i][m][ok][sel] * pz) / rr
                ztr = np.interp(rr, rg, zg)
                Dg_ = np.interp(ztr, zg, Dg); fg_ = np.interp(ztr, zg, fg); Hg_ = np.interp(ztr, zg, Hg)
                s = rr + fg_ * Dg_ * ur
                vr = fg_ * Hg_ * Dg_ * ur
                gal['KD'].append(((1 / (rr * Hg_) - 1) * vr).astype(np.float32))
                gal['KB'].append((Dg_ * g[sel]).astype(np.float32))
                gal['ZTRUE'].append(ztr)
                scale = s / rr
                px, py, pz, rr = px * scale, py * scale, pz * scale, s
                tgt = gal
            else:
                tgt = ran
            ra, dec = hp.vec2ang(np.vstack([px, py, pz]).T, lonlat=True)
            tgt['RA'].append(ra); tgt['DEC'].append(dec); tgt['Z'].append(np.interp(rr, rg, zg))
        if i % 64 == 0:
            print('slab %d %.0fs' % (i, time.time() - t0), flush=True)
    gal = {k: np.concatenate(v) for k, v in gal.items()}
    ran = {k: np.concatenate(v) for k, v in ran.items()}
    gal['CAP'] = capmap[hp.ang2pix(NSIDE, gal['RA'], gal['DEC'], lonlat=True)]
    ran['CAP'] = capmap[hp.ang2pix(NSIDE, ran['RA'], ran['DEC'], lonlat=True)]
    np.savez(D + 'mock_%d.npz' % seed, **{'gal_' + k: v for k, v in gal.items()}, **{'ran_' + k: v for k, v in ran.items()})
    print('mock %d: %d galaxies, %d randoms, %.0fs' % (seed, len(gal['Z']), len(ran['Z']), time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
