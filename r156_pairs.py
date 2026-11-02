"""R-156: per-galaxy shell sums for the size-density cross-correlation (numba pair counter).

For every size ("kappa") galaxy j and every density point i (BGS_BRIGHT galaxies, weight WEIGHT; and randoms)
within separation d < DMAX, with s = x_j - x_i (from the density point to the size galaxy) and the line of sight
the midpoint direction l = (x_i + x_j)/2, mu = s.l/|s||l| (mu > 0: the size galaxy is behind the density point):
    Q[j, d, mu] = sum_{i in gal} w_i - alpha sum_{i in ran} w_i      (alpha = sum w_gal / sum w_ran, per cap)
    R[j, d, mu] = alpha sum_{i in ran} w_i
Then for any per-galaxy mark kappa_j with weights u_j:  xi(d, mu) = sum_j u_j kappa_j Q / sum_j u_j R, and the
Legendre multipoles follow by integrating over the mu bins. Storing Q and R per galaxy makes every null (PSF, stars,
shuffles, jackknife) a weighted sum, with no recount.

usage: python -u r156_pairs.py <out-tag> [--zshuffle SEED] [--ranfrac F]
"""
import sys, time
import numpy as np
import numba as nb
import fitsio

D = 'data/r156/'
DEDGES = np.arange(10.0, 151.0, 10.0)
NMU = 10
ZD = (0.08, 0.56)   # density tracers (wider than the size range by > 150 Mpc/h)
ZK = (0.15, 0.45)   # size galaxies
CELL = 15.0


@nb.njit(parallel=True, fastmath=True)
def count(kx, ky, kz, px, py, pz, pw, cstart, cend, nc, xmin, cell, dedges, nmu):
    nk = kx.shape[0]
    nd = dedges.shape[0] - 1
    out = np.zeros((nk, nd, nmu), dtype=np.float32)
    dmax = dedges[-1]
    dmin = dedges[0]
    dd = dedges[1] - dedges[0]
    reach = int(np.ceil(dmax / cell))
    for j in nb.prange(nk):
        xj, yj, zj = kx[j], ky[j], kz[j]
        cx = int((xj - xmin[0]) / cell)
        cy = int((yj - xmin[1]) / cell)
        cz = int((zj - xmin[2]) / cell)
        for ix in range(max(cx - reach, 0), min(cx + reach + 1, nc[0])):
            for iy in range(max(cy - reach, 0), min(cy + reach + 1, nc[1])):
                for iz in range(max(cz - reach, 0), min(cz + reach + 1, nc[2])):
                    c = (ix * nc[1] + iy) * nc[2] + iz
                    for p in range(cstart[c], cend[c]):
                        sx = xj - px[p]
                        sy = yj - py[p]
                        sz = zj - pz[p]
                        s2 = sx * sx + sy * sy + sz * sz
                        if s2 >= dmax * dmax or s2 < dmin * dmin:
                            continue
                        s = np.sqrt(s2)
                        lx = xj + px[p]
                        ly = yj + py[p]
                        lz = zj + pz[p]
                        mu = (sx * lx + sy * ly + sz * lz) / (s * np.sqrt(lx * lx + ly * ly + lz * lz))
                        ib = int((s - dmin) / dd)
                        im = int((mu + 1.0) * 0.5 * nmu)
                        if im >= nmu:
                            im = nmu - 1
                        if im < 0:
                            im = 0
                        out[j, ib, im] += pw[p]
    return out


def cart(ra, dec, r):
    ra, dec = np.radians(ra), np.radians(dec)
    return r * np.cos(dec) * np.cos(ra), r * np.cos(dec) * np.sin(ra), r * np.sin(dec)


def grid(x, y, z, w, xmin, nc):
    ix = ((x - xmin[0]) / CELL).astype(np.int64)
    iy = ((y - xmin[1]) / CELL).astype(np.int64)
    iz = ((z - xmin[2]) / CELL).astype(np.int64)
    c = (ix * nc[1] + iy) * nc[2] + iz
    o = np.argsort(c, kind='stable')
    c = c[o]
    ncell = nc[0] * nc[1] * nc[2]
    cstart = np.searchsorted(c, np.arange(ncell), 'left')
    cend = np.searchsorted(c, np.arange(ncell), 'right')
    return x[o], y[o], z[o], w[o].astype(np.float64), cstart, cend


def main():
    tag = sys.argv[1]
    zshuffle = int(sys.argv[sys.argv.index('--zshuffle') + 1]) if '--zshuffle' in sys.argv else None
    ranfrac = float(sys.argv[sys.argv.index('--ranfrac') + 1]) if '--ranfrac' in sys.argv else 1.0
    mock = int(sys.argv[sys.argv.index('--mock') + 1]) if '--mock' in sys.argv else None
    kfrac = float(sys.argv[sys.argv.index('--kfrac') + 1]) if '--kfrac' in sys.argv else 1.0
    if mock is not None:
        M = np.load(D + 'mock_%d.npz' % mock)
    import r156_theory as T
    zg = np.linspace(0, 0.7, 1401)
    rg = np.array([T.background(x)[0] for x in zg])
    # every clustering galaxy in the size range (the size cuts are applied later, so this does not wait for them)
    S = {'Z': [], 'RA': [], 'DEC': [], 'CAP': []}
    for cap, cname in enumerate(['NGC', 'SGC']):
        a = fitsio.read(D + 'BGS_BRIGHT_%s_clustering.dat.fits' % cname, columns=['RA', 'DEC', 'Z'])
        for k in ['Z', 'RA', 'DEC']:
            S[k].append(a[k])
        S['CAP'].append(np.full(len(a), cap))
    S = {k: np.concatenate(v) for k, v in S.items()}
    if mock is not None:
        S = {k: M['gal_' + k] for k in ['Z', 'RA', 'DEC', 'CAP']}
    ksel = (S['Z'] > ZK[0]) & (S['Z'] < ZK[1])
    if kfrac < 1:
        ksel &= np.random.default_rng(7).random(len(ksel)) < kfrac
    kidx = np.where(ksel)[0]
    zk = S['Z'][kidx].copy()
    if zshuffle is not None:
        # randomised-redshift null: each size galaxy keeps its sky position but takes another size galaxy's z
        rng = np.random.default_rng(zshuffle)
        for cap in (0, 1):
            m = S['CAP'][kidx] == cap
            zk[m] = rng.permutation(zk[m])
    Q = np.zeros((len(kidx), len(DEDGES) - 1, NMU), np.float32)
    R = np.zeros_like(Q)
    t0 = time.time()
    for cap, cname in enumerate(['NGC', 'SGC']):
        dat = fitsio.read(D + 'BGS_BRIGHT_%s_clustering.dat.fits' % cname, columns=['RA', 'DEC', 'Z', 'WEIGHT'])
        ran = fitsio.read(D + 'BGS_BRIGHT_%s_0_clustering.ran.fits' % cname, columns=['RA', 'DEC', 'Z', 'WEIGHT'])
        if mock is not None:
            def rec(p):
                m = M[p + 'CAP'] == cap
                a = np.zeros(m.sum(), dtype=[('RA', 'f8'), ('DEC', 'f8'), ('Z', 'f8'), ('WEIGHT', 'f8')])
                a['RA'], a['DEC'], a['Z'], a['WEIGHT'] = M[p + 'RA'][m], M[p + 'DEC'][m], M[p + 'Z'][m], 1.0
                return a
            dat, ran = rec('gal_'), rec('ran_')
        dat = dat[(dat['Z'] > ZD[0]) & (dat['Z'] < ZD[1])]
        ran = ran[(ran['Z'] > ZD[0]) & (ran['Z'] < ZD[1])]
        if ranfrac < 1:
            ran = ran[np.random.default_rng(1).random(len(ran)) < ranfrac]
        alpha = dat['WEIGHT'].sum() / ran['WEIGHT'].sum()
        m = S['CAP'][kidx] == cap
        kx, ky, kz = cart(S['RA'][kidx][m], S['DEC'][kidx][m], np.interp(zk[m], zg, rg))
        out = []
        for arr, sign in [(dat, 1.0), (ran, alpha)]:
            px, py, pz = cart(arr['RA'], arr['DEC'], np.interp(arr['Z'], zg, rg))
            allx = np.concatenate([px, kx]); ally = np.concatenate([py, ky]); allz = np.concatenate([pz, kz])
            xmin = np.array([allx.min(), ally.min(), allz.min()]) - 1
            nc = (np.array([allx.max(), ally.max(), allz.max()]) - xmin) // CELL + 2
            nc = nc.astype(np.int64)
            gx, gy, gz, gw, cs, ce = grid(px, py, pz, arr['WEIGHT'] * sign, xmin, nc)
            t1 = time.time()
            out.append(count(kx, ky, kz, gx, gy, gz, gw, cs, ce, nc, xmin, CELL, DEDGES, NMU))
            print('%s %s: %d size galaxies x %d points, %.0fs' % (cname, 'gal' if sign == 1.0 else 'ran',
                                                                  m.sum(), len(px), time.time() - t1), flush=True)
        Q[m] = out[0] - out[1]
        R[m] = out[1]
    np.savez(D + 'pairs_%s.npz' % tag, Q=Q, R=R, kidx=kidx, zk=zk, dedges=DEDGES, nmu=NMU)
    print('DONE %.0fs' % (time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
