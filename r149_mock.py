"""R-149: photo-z isolation leakage and its lensing amplitude on the UniverseMachine Bolshoi-Planck mock (cosmo-deep-16c).

One run = one line of sight (box axis), one photo-z seed, one virtual redshift z_v, one magnitude limit. For every lens
candidate it decides photo-z isolation (Bilicki's KiDS-bright Lorentzian, as R-147) and true isolation (true masses, true
redshift space), the R-147 contamination statistic, and the ESD from the mock's own halos (truncated NFW at each halo's
catalogue mass and radius, plus a satellite lens's own subhalo and the lens's baryons) at Brouwer's 15 g_bar radii and at
fixed radii. Everything is fixed in results-R147.md ("R-149 fixed before computing").

Usage: python -I r149_mock.py <data dir> <out dir> --axis 0 --seed 1 --zv 0.20 --mlim 19.5 [--untrunc] [--mhmin 1e11]
Units: physical Mpc at the snapshot, M_sun (h-free), km/s.
"""
import argparse, json, os, sys, time
import numpy as np
from scipy.spatial import cKDTree
from scipy.integrate import cumulative_trapezoid
from scipy.ndimage import map_coordinates
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r147_leak import GBAR, G_SI, MSUN, MPC, NFW, dcom

H, OMM = 0.678, 0.307                       # Bolshoi-Planck
A_SNAP = 0.827654; Z_SNAP = 1 / A_SNAP - 1
LBOX = 250.0 / H * A_SNAP                   # physical Mpc
C = 299792.458
HZ = lambda z: 100 * H * np.sqrt(OMM * (1 + z) ** 3 + 1 - OMM)
RFIX = np.logspace(np.log10(0.03), np.log10(3.0), 11)                 # fixed radii (Mpc) for Brouwer's App. A benchmark
SIGMA_PC = 1e12                                                       # Mpc^2 -> pc^2
DTYPE = np.dtype([('id', 'i8'), ('descid', 'i8'), ('upid', 'i8'), ('flags', 'i4'), ('uparent_dist', 'f4'),
                  ('pos', 'f4', (6)), ('vmp', 'f4'), ('lvmp', 'f4'), ('mp', 'f4'), ('m', 'f4'), ('v', 'f4'), ('r', 'f4'),
                  ('rank1', 'f4'), ('rank2', 'f4'), ('ra', 'f4'), ('rarank', 'f4'), ('A_UV', 'f4'), ('sm', 'f4'),
                  ('icl', 'f4'), ('sfr', 'f4'), ('obs_sm', 'f4'), ('obs_sfr', 'f4'), ('obs_uv', 'f4'), ('empty', 'f4')],
                 align=True)


def cvir(m, z):
    """Dutton & Maccio 2014, BN98 virial; m in M_sun (h-free)."""
    b = -0.097 + 0.024 * z; a = 0.537 + (1.025 - 0.537) * np.exp(-0.718 * z ** 1.08)
    return 10 ** (a + b * np.log10(m * H / 1e12))


class TruncNFW:
    """ESD of an NFW truncated at r_vir, about a point at projected offset d, in units of M / r_vir^2: G(X, Y; c)."""
    XS = np.logspace(-3, np.log10(25), 220); YS = np.r_[0, np.logspace(-2.5, np.log10(20), 160)]
    CS = np.logspace(np.log10(2), np.log10(30), 13)

    def __init__(self, cache=None):
        if cache and os.path.exists(cache):
            self.T = np.load(cache)['T']
        else:
            self.build()
            if cache: np.savez(cache, T=self.T)
        self.lx, self.ly, self.lc = np.log(self.XS), self.YS, np.log(self.CS)

    def build(self):
        th = (np.arange(1024) + 0.5) / 1024 * np.pi
        self.T = np.empty((self.XS.size, self.YS.size, self.CS.size))
        for k, c in enumerate(self.CS):
            xx = np.r_[np.logspace(-5, -0.0001, 1500)]
            # Sigma(X) of the truncated profile, total mass 1, r_vir = 1: rho = rho_s / (c r (1 + c r)^2) for r < 1
            l = np.r_[0, np.logspace(-6, 0, 600)][None, :]
            lmax = np.sqrt(1 - xx ** 2)[:, None]
            r = np.sqrt(xx[:, None] ** 2 + (l * lmax) ** 2)
            rho = 1 / (c * r * (1 + c * r) ** 2)
            sig = 2 * np.trapezoid(rho, l, axis=1) * lmax[:, 0]
            norm = 4 * np.pi / c ** 3 * (np.log(1 + c) - c / (1 + c))      # mass of rho with rho_s = 1, r_s = 1/c
            sig /= norm
            sfun = lambda R: np.interp(R, xx, sig, right=0.0)
            for j, y in enumerate(self.YS):
                rr = np.sqrt(self.XS[:, None] ** 2 + y ** 2 + 2 * self.XS[:, None] * y * np.cos(th)[None, :])
                ring = sfun(rr).mean(axis=1)
                if y == 0:
                    ring = sfun(self.XS)
                xf = np.logspace(-5, np.log10(25), 1500)
                rf = np.sqrt(xf[:, None] ** 2 + y ** 2 + 2 * xf[:, None] * y * np.cos(th[::4])[None, :])
                ringf = sfun(rf).mean(axis=1) if y > 0 else sfun(xf)
                enc = cumulative_trapezoid(2 * np.pi * ringf * xf, xf, initial=0)
                mean_in = np.interp(self.XS, xf, enc) / (np.pi * self.XS ** 2)
                self.T[:, j, k] = mean_in - ring

    def __call__(self, R, m, rv, c, d):
        """ESD (M_sun pc^-2) at radii R (Mpc) from halos of mass m (M_sun), r_vir rv (Mpc), conc c, offset d (Mpc)."""
        R, m, rv, c, d = np.broadcast_arrays(R, m, rv, c, d)
        X = R / rv; Y = d / rv
        pt = m / (np.pi * R ** 2)
        out = np.where(X > Y + 1, pt, 0.0)                                      # exact outside the halo
        inside = (X <= Y + 1) & (Y <= self.YS[-1]) & (X <= self.XS[-1])
        if inside.any():
            Xi = np.clip(X[inside], self.XS[0], self.XS[-1])
            fx = np.interp(np.log(Xi), self.lx, np.arange(self.XS.size))
            fy = np.interp(Y[inside], self.ly, np.arange(self.YS.size))
            fc = np.interp(np.log(np.clip(c[inside], self.CS[0], self.CS[-1])), self.lc, np.arange(self.CS.size))
            g = map_coordinates(self.T, np.vstack([fx, fy, fc]), order=1, mode='nearest')
            out[inside] = (m / rv ** 2)[inside] * g
        far = (X <= Y + 1) & (Y > self.YS[-1])                                  # point mass beyond 20 r_vir
        out[far] = np.where(R[far] > d[far], pt[far], 0.0)
        return out / SIGMA_PC


class Untrunc:
    """R-147's untruncated NFW (Duffy c), as a variant: same call signature."""
    def __init__(self):
        self.n = NFW()

    def __call__(self, R, m, rv, c, d):
        R, m, d = np.broadcast_arrays(R, m, d)
        mm = m * 0.7 / H                                                       # NFW.params takes h70 M_200m; close enough
        return self.n.dsig_off(R, mm, Z_SNAP, np.maximum(d, 1e-4))


def load(datadir):
    h = np.fromfile(os.path.join(datadir, 'um_bp_a0.827654.bin'), dtype=DTYPE)
    return h


def periodic(dx):
    return dx - LBOX * np.round(dx / LBOX)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('data'); ap.add_argument('out')
    ap.add_argument('--axis', type=int, default=0); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--zv', type=float, default=0.20); ap.add_argument('--mlim', type=float, default=19.5)
    ap.add_argument('--untrunc', action='store_true'); ap.add_argument('--mhmin', type=float, default=1e11)
    ap.add_argument('--nesd', type=int, default=60000, help='max lenses per sample and type for the ESD')
    a = ap.parse_args(); t0 = time.time(); os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(1000 * a.seed + 10 * a.axis + int(round(a.zv * 100)) + int(a.mlim * 7))
    h = load(a.data)
    los = a.axis; sky = [k for k in range(3) if k != los]
    pos = h['pos'][:, :3].astype(np.float64) / H * A_SNAP                       # physical Mpc
    vlos = h['pos'][:, 3 + los].astype(np.float64)
    host_id = np.where(h['upid'] == -1, h['id'], h['upid'])
    distinct = h['upid'] == -1
    # ---- galaxies
    sm = h['obs_sm'].astype(np.float64)
    gal = np.where(sm > 10 ** 8.3)[0]
    lm = np.log10(sm[gal]); q = (h['obs_sfr'][gal] / sm[gal]) < 1e-11
    ml = json.load(open(os.path.join(a.data, 'maglim.json')))[f'{a.zv:.2f}']
    zobs = a.zv + (1 + a.zv) * vlos[gal] / C
    dm = 5 * np.log10(dcom(np.maximum(zobs, 1e-3)) * (1 + zobs) * 1e6 / 10)
    coef = np.where(q[:, None], np.array(ml['Q'])[None, :], np.array(ml['SF'])[None, :])
    mag = coef[:, 0] + coef[:, 1] * (lm - 10.5) + coef[:, 2] * rng.standard_normal(gal.size) + dm
    keep = mag < a.mlim
    gal, lm, q, zobs, mag = gal[keep], lm[keep], q[keep], zobs[keep], mag[keep]
    N = gal.size
    X = pos[gal].copy()
    s_true = X[:, los] + vlos[gal] / HZ(Z_SNAP)                                 # true redshift space (physical Mpc)
    # ---- photo-z view
    nu = 2 * 2.613 - 1; sc = 0.0149 * np.sqrt(2 * 2.613 / nu)
    dz = sc * rng.standard_t(nu, N)
    zph = zobs + (1 + a.zv) * dz
    s_ph = s_true + C * (1 + a.zv) * dz / HZ(a.zv) / (1 + a.zv)                 # physical displacement at z_v
    dl = lambda z: dcom(np.clip(z, 1e-3, 2)) * (1 + np.clip(z, 1e-3, 2))
    lmp = lm + 2 * np.log10(dl(zph) / dl(zobs)) + 0.12 * rng.standard_normal(N)
    cand = np.where((zph > 0.1) & (zph < 0.5) & (lmp > 8.5) & (lmp < 11.0))[0]
    print(f'z_v {a.zv} r < {a.mlim} axis {los} seed {a.seed}: {N} galaxies, {cand.size} candidates, '
          f'{q[cand].mean():.2f} quenched ({time.time() - t0:.0f}s)', flush=True)

    def coords(s):
        Y = X.copy(); Y[:, los] = np.mod(s, LBOX); return np.mod(Y, LBOX)

    def isolated(s, m, R):
        P = coords(s); tree = cKDTree(P, boxsize=LBOX)
        bad = np.zeros(cand.size, bool)
        for st in range(0, cand.size, 20000):
            ii = cand[st:st + 20000]
            lists = tree.query_ball_point(P[ii], R, workers=4)
            n = np.fromiter((len(l) for l in lists), int, ii.size)
            i = np.repeat(np.arange(ii.size) + st, n); j = np.fromiter((x for l in lists for x in l), int, n.sum())
            k = (cand[i] != j) & (m[j] > m[cand[i]] - 1)
            bad[np.unique(i[k])] = True
        return ~bad
    iso_p = {R: isolated(s_ph, lmp, R) for R in (3.0, 4.0)}
    iso_t = {R: isolated(s_true, lm, R) for R in (3.0, 4.0)}
    print('photo-z isolated ' + ', '.join(f'{R:g}: {iso_p[R].mean():.3f}' for R in iso_p) + '; truly isolated ' +
          ', '.join(f'{R:g}: {iso_t[R].mean():.3f}' for R in iso_t) + f' ({time.time() - t0:.0f}s)', flush=True)
    # ---- R-147's contamination statistic for photo-z-isolated lenses: real > 10 % neighbour, projected R, 500 km/s
    lens_c = cand[iso_p[3.0] | iso_p[4.0]]
    S = np.c_[X[:, sky[0]], X[:, sky[1]], s_true * HZ(Z_SNAP) / 500.0 * 4.0]   # 500 km/s -> 4 Mpc along the axis
    Sb = np.array([LBOX, LBOX, LBOX * HZ(Z_SNAP) / 500.0 * 4.0])
    S = np.mod(S, Sb); tree = cKDTree(S, boxsize=Sb)
    rpmin = np.full(lens_c.size, np.inf); hidsat = np.zeros(lens_c.size, bool)
    hid_g = host_id[gal]; sat_g = h['upid'][gal] != -1
    for st in range(0, lens_c.size, 20000):
        ii = lens_c[st:st + 20000]
        lists = tree.query_ball_point(S[ii], 4.0 * np.sqrt(2) + 0.01, workers=4)
        n = np.fromiter((len(l) for l in lists), int, ii.size)
        i = np.repeat(np.arange(ii.size) + st, n); j = np.fromiter((x for l in lists for x in l), int, n.sum())
        li = lens_c[i]
        dxy = periodic(S[j, :2] - S[li, :2]); rp = np.hypot(dxy[:, 0], dxy[:, 1])
        dv = np.abs(periodic((s_true[j] - s_true[li])) * HZ(Z_SNAP))
        k = (li != j) & (lm[j] > lm[li] - 1) & (dv < 500) & (rp < 4.0)
        np.minimum.at(rpmin, i[k], rp[k])
        k3 = k & (rp < 3.0) & sat_g[li] & (hid_g[j] == hid_g[li])
        hidsat[np.unique(i[k3])] = True
    print(f'contamination done ({time.time() - t0:.0f}s)', flush=True)
    # ---- ESD from the mock's halos, for photo-z-isolated and truly isolated lenses (subsampled per type)
    want = np.zeros(cand.size, bool)
    for flags in (iso_p[3.0], iso_t[3.0], iso_p[4.0], iso_t[4.0]):
        for t in (q[cand], ~q[cand]):
            idx = np.where(flags & t)[0]
            if idx.size > a.nesd: idx = rng.choice(idx, a.nesd, replace=False)
            want[idx] = True
    el = cand[want]                                                               # galaxy-sample indices
    hx = np.where(distinct & (h['m'] / H > a.mhmin))[0]
    HP = np.c_[pos[hx][:, sky[0]], pos[hx][:, sky[1]], pos[hx][:, los] * (8.0 / 30.0)]
    Hb = np.array([LBOX, LBOX, LBOX * 8.0 / 30.0]); HP = np.mod(HP, Hb)
    htree = cKDTree(HP, boxsize=Hb)
    hm = h['m'][hx].astype(np.float64) / H; hr = h['r'][hx].astype(np.float64) / 1e3 / H * A_SNAP
    hc = cvir(hm, Z_SNAP); hid_h = h['id'][hx]
    prof = Untrunc() if a.untrunc else TruncNFW(os.path.join(a.data, 'tnfw_table.npz'))
    print(f'profile table ready; {hx.size} halos > {a.mhmin:.0e}; ESD for {el.size} lenses ({time.time() - t0:.0f}s)',
          flush=True)
    mgal = 10 ** lmp[el] * (1 + 10 ** (-0.69 * lmp[el] + 6.63))
    Rk = np.sqrt(G_SI * mgal[:, None] * MSUN / GBAR[None, :]) / MPC
    Rall = np.c_[Rk, np.broadcast_to(RFIX, (el.size, RFIX.size))]                  # (n, 15 + 11)
    ds = np.zeros_like(Rall)
    LP = np.c_[X[el][:, sky[0]], X[el][:, sky[1]], X[el][:, los] * (8.0 / 30.0)]; LP = np.mod(LP, Hb)
    own_host = hid_g[el]
    for st in range(0, el.size, 5000):
        sl = slice(st, st + 5000)
        lists = htree.query_ball_point(LP[sl], 8.0 * np.sqrt(2), workers=4)
        n = np.fromiter((len(l) for l in lists), int, len(lists))
        i = np.repeat(np.arange(len(lists)) + st, n); j = np.fromiter((x for l in lists for x in l), int, n.sum())
        dxy = periodic(HP[j, :2] - LP[i, :2]); d = np.hypot(dxy[:, 0], dxy[:, 1])
        dl_ = np.abs(periodic((HP[j, 2] - LP[i, 2]) * 30.0 / 8.0))
        k = (d < 8.0) & (dl_ < 30.0)
        i, j, d = i[k], j[k], d[k]
        # the lens's own host at offset 0 when the lens is its central (galaxies sit at their halo centres)
        d = np.where(hid_h[j] == own_host[i], np.where(sat_g[el[i]], d, 0.0), d)
        v = prof(Rall[i], hm[j][:, None], hr[j][:, None], hc[j][:, None], d[:, None])
        np.add.at(ds, i, v)
    # own host when below the halo-mass floor, own subhalo for satellites, baryons as a point mass
    order = np.argsort(h['id']); ids_sorted = h['id'][order]
    def row(ids):
        return order[np.searchsorted(ids_sorted, ids)]
    hostrow = row(own_host); hm_own = h['m'][hostrow] / H
    low = hm_own <= a.mhmin
    if low.any():
        r_ = h['r'][hostrow[low]] / 1e3 / H * A_SNAP
        dd = np.zeros(low.sum())
        satl = sat_g[el[low]]
        if satl.any():
            dxy = periodic(pos[hostrow[low]][:, sky] - X[el[low]][:, sky]); dd = np.where(satl, np.hypot(*dxy.T), 0.0)
        ds[low] += prof(Rall[low], hm_own[low][:, None], r_[:, None], cvir(hm_own[low], Z_SNAP)[:, None], dd[:, None])
    s_ = sat_g[el]
    if s_.any():
        gr = gal[el[s_]]
        ms = h['m'][gr].astype(np.float64) / H; rs = h['r'][gr] / 1e3 / H * A_SNAP
        ok = (ms > 0) & (rs > 0)
        sub = np.zeros((s_.sum(), Rall.shape[1]))
        sub[ok] = prof(Rall[s_][ok], ms[ok][:, None], rs[ok][:, None], cvir(ms[ok], Z_SNAP)[:, None], np.zeros((ok.sum(), 1)))
        ds[s_] += sub
    ds += (10 ** lm[el] * (1 + 10 ** (-0.69 * lm[el] + 6.63)))[:, None] / (np.pi * Rall ** 2) / SIGMA_PC
    print(f'ESD done ({time.time() - t0:.0f}s)', flush=True)
    jk = (np.floor(X[cand][:, sky[0]] / LBOX * 4) * 4 + np.floor(X[cand][:, sky[1]] / LBOX * 4)).astype(int)
    rp_full = np.full(cand.size, np.inf); rp_full[np.searchsorted(cand, lens_c)] = rpmin
    hs_full = np.zeros(cand.size, bool); hs_full[np.searchsorted(cand, lens_c)] = hidsat
    esd_rows = np.full(cand.size, -1); esd_rows[np.where(want)[0]] = np.arange(el.size)
    tag = f'{"untrunc_" if a.untrunc else ""}{"mh%.1f_" % np.log10(a.mhmin) if a.mhmin != 1e11 else ""}' \
          f'z{a.zv:.2f}_r{a.mlim:g}_ax{los}_s{a.seed}'
    fn = os.path.join(a.out, f'r149_{tag}.npz')
    np.savez_compressed(fn, q=q[cand], lmp=lmp[cand], lm=lm[cand], sat=sat_g[cand], jk=jk, rpmin=rp_full, hidsat=hs_full,
                        isop3=iso_p[3.0], isop4=iso_p[4.0], isot3=iso_t[3.0], isot4=iso_t[4.0], esd_row=esd_rows, ds=ds,
                        Rall=Rall, zv=a.zv, mlim=a.mlim)
    print(f'wrote {fn} ({time.time() - t0:.0f}s)', flush=True)


if __name__ == '__main__':
    main()
