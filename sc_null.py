"""The LCDM null for R-137's supercluster statistic, from the Websky light-cone (R-140, cosmo-deep-13b).

R-137 (sc_lens.py) measured, on ACT DR6 lensing, E - E_control for 365 eRASS1 superclusters: the band-filtered convergence
integrated over each supercluster's aperture, minus its members' NFW + linear two-halo model, in units of that model, against
the same on superpositions of isolated clusters. This script runs the same pipeline on the Websky simulation (Stein et al.
2020, arXiv:2001.08787): its halo light-cone (halos.pksc, sorted by mass; the top 30 M halos are enough) and its CMB-lensing
convergence (kap.fits), so that E, E_control and E - E_control have a LCDM expectation and a distribution.

What is fixed, and why, is in results-R140.md ('Fixed before running'). In short:
  * sky: patch 0 is Websky's own frame read as equatorial (ACT's Websky foreground reconstructions use the same convention);
    patches 1..NPATCH-1 are random rotations of the whole simulated sky (halos and kappa together) into the real footprint
    (the eRASS1 hemisphere, the clusters' optical footprint, the ACT DR6 baseline mask), so each patch is a data-like sample.
  * clusters: halo M200m (Websky's 4/3 pi rho_m R_TH^3) -> M500c with Duffy et al. 2008's c200m; observed redshift = the
    cosmological one + the line-of-sight velocity + a redshift error drawn from the eRASS1 clusters at that redshift; the
    selection a threshold M500c_min(z) abundance-matched to the number density of the supercluster finder's input sample
    (sc_members.load_sample) per 0.02 in z. Each selected halo gets a mass-posterior width drawn from the catalogue's.
  * finder: R-137's friends-of-friends (sc_members.linking_length / the 'min' rule), on the simulated sample.
  * statistic, apertures, member model, control, weights: sc_lens.py's, unchanged (imported).
  * noise: the ACT DR6 noise realisations, recon_i - input_i (the release's simulated reconstructions minus their input kappa),
    added to the noiseless simulated integrals; the weights' variances from the recon_i, as R-137.

usage (python -I; data under data/, gitignored):
  python -I sc_null.py halos                     # read the halo light-cone -> data/websky/halos_sel.npz
  python -I sc_null.py noise --worker k --nworkers n   # stream ACT sims 1..NSIM: band-filtered recon and noise maps
  python -I sc_null.py prep --patch p            # sample, FoF, apertures, controls, models for patch p
  python -I sc_null.py integrate                 # the aperture integrals on the Websky map (and noise maps) per patch
  python -I sc_null.py fit                       # the null table
"""
import argparse, glob, os, subprocess, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sc_lens as L  # noqa: E402

DATA = os.path.join(HERE, 'data')
WS = os.path.join(DATA, 'websky')
OUT = os.path.join(DATA, 'sc_null')
BASE = OUT
NOISE = os.path.join(BASE, 'noise')
# variants (information rows, defined after the main table): SC_NULL_VAR=s0 (a sharp mass threshold, s = 0, the brief's literal
# 'mass-redshift threshold'), SC_NULL_VAR=mnoise (the model and the control matching use a noisy mass, M500 x exp(sig N(0, 1)),
# as the catalogue's masses are noisy estimates, instead of the true one)
VAR = os.environ.get('SC_NULL_VAR', '')
if VAR:
    OUT = os.path.join(BASE, f'var_{VAR}')
NSIM = 100
NPATCH = 25
NCTRL = L.NCTRL
# Websky cosmology (Stein et al. 2020)
WS_H, WS_OM = 0.68, 0.31
RHO_M = 2.775e11 * WS_H ** 2 * WS_OM  # Msun / Mpc^3
NHALO = 30_000_000
CLIGHT = 299792.458
ZCUT = 0.9


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


# --------------------------------------------------------------------------------------------------------------- noise
def noise(args):
    """fetch each ACT DR6 baseline simulation and its input kappa, keep two band-filtered nside-512 maps (recon, recon-input)."""
    import healpy as hp
    os.makedirs(NOISE, exist_ok=True)
    tmp = os.path.join(DATA, 'act_dr6', 'tmp'); os.makedirs(tmp, exist_ok=True)
    for i in range(1, NSIM + 1):
        if i % args.nworkers != args.worker:
            continue
        n = f'{i:04d}'
        fo = os.path.join(NOISE, f'noise_{n}.npy')
        if os.path.exists(fo):
            continue
        fr = os.path.join(tmp, f'recon_{n}.fits'); fi = os.path.join(tmp, f'input_{n}.fits')
        ok = subprocess.call(['sh', os.path.join(HERE, 'act_dr6_fetch.sh'),
                              f'maps/baseline/simulations/kappa_alm_sim_act_dr6_lensing_v1_baseline_{n}.fits', fr]) == 0
        ok = ok and subprocess.call(['sh', os.path.join(HERE, 'act_dr6_fetch.sh'),
                                     f'sim_inputs/kappa_alm/input_kappa_alm_sim_{n}.fits', fi]) == 0
        if not ok:
            log(f'sim {n}: fetch failed'); continue
        mr = L.map_from_alm(L.read_alm(fr)); mi = L.map_from_alm(L.read_alm(fi))
        np.save(os.path.join(NOISE, f'recon_{n}.npy'), mr.astype(np.float32))
        np.save(fo + '.tmp.npy', (mr - mi).astype(np.float32)); os.replace(fo + '.tmp.npy', fo)
        os.remove(fr); os.remove(fi)
        log(f'sim {n} done')


# --------------------------------------------------------------------------------------------------------------- halos
def m200m_to_m500c(M, z):
    """NFW with Duffy et al. 2008's c200m(M200m, z) (full sample): the radius x = r/r200m where the mean enclosed density is
    500 rho_c(z), i.e. m(c x) / (m(c) x^3) = 2.5 / Omega_m(z); M500c = M200m m(c x) / m(c). Websky cosmology."""
    c = 10.14 * (M * WS_H / 2e12) ** -0.081 * (1 + z) ** -1.01
    omz = WS_OM * (1 + z) ** 3 / (WS_OM * (1 + z) ** 3 + 1 - WS_OM)
    target = 2.5 / omz
    mf = lambda y: np.log1p(y) - y / (1 + y)
    lo = np.full(len(M), 0.2); hi = np.full(len(M), 1.0)
    for _ in range(40):
        x = 0.5 * (lo + hi)
        f = mf(c * x) / (mf(c) * x ** 3) - target   # decreasing in x
        lo = np.where(f > 0, x, lo); hi = np.where(f > 0, hi, x)
    x = 0.5 * (lo + hi)
    return M * mf(c * x) / mf(c)


def halos(args):
    """the top NHALO halos of Websky (sorted by mass): sky position (Websky frame), redshift with the line-of-sight velocity,
    M500c; kept if z < ZCUT and M500c > 5e12."""
    from astropy.cosmology import FlatLambdaCDM
    os.makedirs(OUT, exist_ok=True)
    raw = np.memmap(os.path.join(WS, 'halos_top30M.bin'), dtype=np.float32, mode='r', offset=12, shape=(NHALO, 10))
    co = FlatLambdaCDM(H0=100 * WS_H, Om0=WS_OM)
    zg = np.linspace(0, 1.5, 3001); chig = co.comoving_distance(zg).value
    keep = []
    for a in range(0, NHALO, 2_000_000):
        b = np.array(raw[a:a + 2_000_000])
        x = b[:, :3].astype(np.float64); chi = np.linalg.norm(x, axis=1)
        zc = np.interp(chi, chig, zg)
        sel = zc < ZCUT
        x, chi, zc, b = x[sel], chi[sel], zc[sel], b[sel]
        u = x / chi[:, None]
        vr = (b[:, 3:6] * u).sum(1)
        zv = (1 + zc) * (1 + vr / CLIGHT) - 1
        M = 4 / 3 * np.pi * RHO_M * b[:, 6].astype(np.float64) ** 3
        M5 = m200m_to_m500c(M, zc)
        s2 = M5 > 5e12
        keep.append(np.c_[u[s2], zc[s2], zv[s2], M[s2], M5[s2]])
        log(f'halos {a + len(b)}..: kept {s2.sum()} (min M200m in chunk {M.min():.2e}, z < {ZCUT})')
    k = np.concatenate(keep)
    np.savez(os.path.join(OUT, 'halos.npz'), u=k[:, :3].astype(np.float32), zc=k[:, 3], zv=k[:, 4], M200m=k[:, 5], M500=k[:, 6])
    log(f'halos.npz: {len(k)} halos, M500c from {k[:, 6].min():.2e}; the 30M-th halo has M200m '
        f'{4 / 3 * np.pi * RHO_M * float(raw[-1, 6]) ** 3:.2e}')


# ----------------------------------------------------------------------------------------------------------- selection
ZB = np.arange(0.0, ZCUT + 1e-9, 0.02)


def data_sample():
    """the supercluster finder's input sample (sc_members.load_sample): positions, z, M500, posterior widths, its footprint."""
    import healpy as hp
    from sc_members import load_sample
    p, sel = load_sample()
    ra = np.asarray(p['RA_XFIT'], float); dec = np.asarray(p['DEC_XFIT'], float)
    bad = ~np.isfinite(ra) | ~np.isfinite(dec); ra[bad] = np.asarray(p['RA'], float)[bad]; dec[bad] = np.asarray(p['DEC'], float)[bad]
    z = np.asarray(p['BEST_Z'], float); ze = np.asarray(p['BEST_ZERR'], float)
    M = np.asarray(p['M500'], float) * 1e13
    sig = 0.5 * (np.log(np.asarray(p['M500_H'], float)) - np.log(np.asarray(p['M500_L'], float)))
    sig = np.where(np.isfinite(sig) & (sig > 0), sig, np.nanmedian(sig))
    s = sel  # the finder's input: every cluster counts for the density; a mass, where it has one, for the medians
    M = np.where(np.isfinite(M) & (M > 0), M, np.nan)
    foot = np.bincount(hp.ang2pix(16, ra[s], dec[s], lonlat=True), minlength=hp.nside2npix(16)) > 0
    area = foot.sum() * hp.nside2pixarea(16, degrees=True)
    return dict(ra=ra[s], dec=dec[s], z=z[s], ze=np.where(np.isfinite(ze), ze, 0.0)[s], M=M[s], sig=sig[s]), foot, area


def select(args):
    """the simulated cluster sample over the whole Websky sky: observed z = cosmological + velocity + a redshift error drawn
    from the eRASS1 sample at that z; selection probability Phi((ln M500 - ln Mc(z)) / s); s one number, Mc(z) per 0.02 in z,
    so that the expected number density and the median M500 of the selected match the eRASS1 sample's in every bin; s chosen
    to minimise the squared log difference of the medians at 0.05 < z < 0.6."""
    from scipy.special import ndtr
    d, foot, area = data_sample()
    h = np.load(os.path.join(BASE, 'halos.npz'))
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(140)
    zv = h['zv']; M = h['M500']
    # redshift errors and posterior widths: drawn from the eRASS1 sample's clusters in the same 0.02 bin of z
    bi_d = np.clip(np.digitize(d['z'], ZB) - 1, 0, len(ZB) - 2)
    bi_h = np.clip(np.digitize(zv, ZB) - 1, 0, len(ZB) - 2)
    ze = np.zeros(len(zv)); sg = np.zeros(len(zv))
    for b in range(len(ZB) - 1):
        hb = np.where(bi_h == b)[0]; db = np.where(bi_d == b)[0]
        if len(db) == 0:
            db = np.where(np.abs(d['z'] - ZB[b]) < 0.06)[0]
        j = rng.choice(db, len(hb))
        ze[hb] = d['ze'][j]; sg[hb] = d['sig'][j]
    zo = zv + rng.standard_normal(len(zv)) * ze
    uu = rng.random(len(zv))
    bi = np.clip(np.digitize(zo, ZB) - 1, 0, len(ZB) - 2)
    nd = np.array([((d['z'] >= a) & (d['z'] < b)).sum() / area for a, b in zip(ZB[:-1], ZB[1:])])
    medd = np.array([np.nanmedian(d['M'][(d['z'] >= a) & (d['z'] < b)]) if ((d['z'] >= a) & (d['z'] < b)).sum() > 5 else np.nan
                     for a, b in zip(ZB[:-1], ZB[1:])])
    lnM = np.log(M); sky = 41252.96
    order = [np.where(bi == b)[0] for b in range(len(ZB) - 1)]

    def solve(s):
        lnMc = np.zeros(len(ZB) - 1); med = np.zeros(len(ZB) - 1)
        for b, ix in enumerate(order):
            x = lnM[ix]; target = nd[b] * sky
            lo, hi = np.log(1e12), np.log(3e15)
            for _ in range(50):
                c = 0.5 * (lo + hi)
                n = ndtr((x - c) / s).sum() if s > 0 else (x > c).sum()
                lo, hi = (c, hi) if n > target else (lo, c)
            lnMc[b] = 0.5 * (lo + hi)
            w = ndtr((x - lnMc[b]) / s) if s > 0 else (x > lnMc[b]).astype(float)
            o = np.argsort(x); cw = np.cumsum(w[o]); med[b] = np.exp(x[o][np.searchsorted(cw, cw[-1] / 2)])
        return lnMc, med
    zc = 0.5 * (ZB[:-1] + ZB[1:]); fitb = (zc > 0.05) & (zc < 0.6)
    best = None
    for s in ((0.0,) if VAR == 's0' else (0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0)):
        lnMc, med = solve(s)
        chi = np.nansum((np.log(med / medd))[fitb] ** 2)
        log(f'selection scatter s = {s:.1f}: sum (ln median ratio)^2 = {chi:.3f}; median ratio sim/data at z 0.1, 0.3, 0.5: '
            + ', '.join(f'{(med / medd)[np.argmin(np.abs(zc - q))]:.2f}' for q in (0.1, 0.3, 0.5)))
        if best is None or chi < best[0]:
            best = (chi, s, lnMc, med)
    chi, s, lnMc, med = best
    p = ndtr((lnM - np.interp(zo, zc, lnMc)) / s) if s > 0 else (lnM > np.interp(zo, zc, lnMc)).astype(float)
    S = (uu < p) & (zo > 0.005)  # the finder's input has z > 0.005
    log(f'chosen s = {s}; selected {S.sum()} of {len(S)} halos over the full sky (the data: {len(d["z"])} in {area:.0f} deg2, '
        f'i.e. {len(d["z"]) / area * sky:.0f} full-sky equivalent)')
    for q in (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6):
        b = np.argmin(np.abs(zc - q))
        log(f'  z {zc[b]:.2f}: Mc {np.exp(lnMc[b]):.2e}; median M500 sim {med[b]:.2e} data {medd[b]:.2e}; n/deg2 data {nd[b]:.4f} '
            f'sim {(S & (bi == b)).sum() / sky:.4f}')
    np.savez(os.path.join(OUT, 'selected.npz'), u=h['u'][S], z=zo[S], zc=h['zc'][S], M=M[S], sig=sg[S], s=s, lnMc=lnMc, zc_bins=zc,
             medd=medd, meds=med, nd=nd, foot=foot, area=area)


# ------------------------------------------------------------------------------------------------------------ patches
def rotation(patch):
    if patch == 0:
        return np.eye(3)
    from scipy.spatial.transform import Rotation
    return Rotation.random(random_state=1000 + patch).as_matrix()


def vec2radec(v):
    return np.degrees(np.arctan2(v[:, 1], v[:, 0])) % 360, np.degrees(np.arcsin(np.clip(v[:, 2], -1, 1)))


def fof_groups(ra, dec, z):
    """R-137's finder (sc_members.build): FoF on 3D comoving positions (H0 70, Om 0.3), l = (n f)^-1/3 per cluster from the
    clusters in z +- 0.005 (one exposure bin: the simulated selection does not depend on exposure), the 'min' rule."""
    import sc_members as SM
    dc, _ = SM.comoving(z)
    r, d = np.radians(ra), np.radians(dec)
    xyz = dc[:, None] * np.c_[np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)]
    l = SM.linking_length(z, np.zeros(len(z)), 1)
    pairs, dist = SM.fof(xyz, l)
    keep = dist < np.minimum(l[pairs[:, 0]], l[pairs[:, 1]])
    return SM.groups_from(pairs[keep], len(z)), l


def prep(args):
    import healpy as hp
    from scipy import sparse
    patch = args.patch
    po = os.path.join(OUT, f'patch{patch:02d}'); os.makedirs(po, exist_ok=True)
    S = np.load(os.path.join(OUT if VAR == 's0' else BASE, 'selected.npz'))
    R = rotation(patch)
    u = S['u'].astype(np.float64) @ R.T
    ra, dec = vec2radec(u)
    # the eRASS1 footprint of the finder's input sample (nside 16, from the data's own clusters), in the observed frame
    inf = S['foot'][hp.vec2pix(16, u[:, 0], u[:, 1], u[:, 2])]
    idx = np.where(inf)[0]
    cl = dict(ra=ra[idx], dec=dec[idx], z=S['z'][idx], M=S['M'][idx], sig=S['sig'][idx])
    n = len(idx)
    if VAR == 'mnoise':
        cl['M'] = cl['M'] * np.exp(cl['sig'] * np.random.default_rng(2400 + patch).standard_normal(n))
    log(f'patch {patch}: {n} simulated clusters in the eRASS1 footprint ({S["area"]:.0f} deg2)')
    lab, l = fof_groups(cl['ra'], cl['dec'], cl['z'])
    order = np.argsort(lab); ls_ = lab[order]
    st = np.r_[0, np.where(np.diff(ls_))[0] + 1]; en = np.r_[st[1:], len(ls_)]
    groups = [order[a:b] for a, b in zip(st, en) if b - a >= 2]
    member = np.zeros(n, bool)
    for g in groups:
        member[g] = True
    log(f'  FoF: {len(groups)} systems ({sum(len(g) == 2 for g in groups)} pairs, {sum(len(g) >= 3 for g in groups)} with >= 3); '
        f'median linking length {np.nanmedian(l):.1f} Mpc')
    b = L.binary_mask('baseline')
    co = L.r121.cosmo(); chi_of = lambda zz: co.comoving_distance(zz).value
    rng = np.random.default_rng(140 + patch)
    rows = []
    for g in groups:
        v = L.unitvec(cl['ra'][g], cl['dec'][g]).mean(0); v /= np.linalg.norm(v)
        ra0, dec0 = vec2radec(v[None])[0][0], vec2radec(v[None])[1][0]
        zk = cl['z'][g].mean()
        if not (L.ZMIN < zk < L.ZMAX):
            continue
        chi = chi_of(zk)
        x, y = L.tangent(ra0, dec0, L.unitvec(cl['ra'][g], cl['dec'][g]))
        P = np.c_[x, y] * chi
        H = L.Hull(P, L.HULL_MARGIN)
        ap = {}
        for kind, obj in (('hull', H), ('circle', L.CIRCLE_R)):
            pix, _ = L.aperture_pixels(ra0, dec0, chi, (kind, obj))
            if len(pix) == 0 or b[pix].mean() < L.CLEAN_MIN:
                break
            ap[kind] = pix
        if len(ap) < 2:
            continue
        rows.append(dict(z=zk, chi=chi, ra=ra0, dec=dec0, mem=g, P=P, H=H, ap=ap, N=len(g)))
    pixc = hp.ang2pix(L.NSIDE, cl['ra'], cl['dec'], lonlat=True)
    iso_idx = np.where(~member & b[pixc])[0]
    log(f'  sample: {len(rows)} systems at {L.ZMIN} < z < {L.ZMAX} with both apertures >= 80 % clean: '
        f'{sum(r["N"] == 2 for r in rows)} pairs, {sum(r["N"] >= 3 for r in rows)} with >= 3; isolated pool {len(iso_idx)}')
    SETS = ['hull', 'circle', 'ctrl_hull', 'ctrl_circle']
    ap_pix, ap_meta, contrib = [], [], []
    for r_i, r in enumerate(rows):
        for kind in ('hull', 'circle'):
            ap_pix.append(r['ap'][kind]); ap_meta.append((SETS.index(kind), r_i, -1, -1)); contrib.append(r['mem'])
    nfail = 0; pool = iso_idx
    for real in range(NCTRL):
        for r_i, r in enumerate(rows):
            for slot, i in enumerate(r['mem']):
                dz = np.abs(cl['z'][pool] - cl['z'][i]); dm = np.abs(np.log(cl['M'][pool] / cl['M'][i]))
                cand = pool[(dz < L.MATCH_DZ) & (dm < L.MATCH_DLNM)]
                if len(cand) < 5:
                    cand = pool[np.argsort(dz / L.MATCH_DZ + dm / L.MATCH_DLNM)[:20]]
                for attempt in range(20):
                    j = rng.choice(cand)
                    chij = chi_of(cl['z'][j])
                    got = {}
                    for kind, obj in (('hull', r['H']), ('circle', L.CIRCLE_R)):
                        pix, _ = L.aperture_pixels(cl['ra'][j], cl['dec'][j], chij, (kind, obj), offset=r['P'][slot])
                        if len(pix) == 0 or b[pix].mean() < L.CLEAN_MIN:
                            break
                        got[kind] = pix
                    if len(got) == 2:
                        break
                else:
                    nfail += 1; continue
                for kind in ('hull', 'circle'):
                    ap_pix.append(got[kind]); ap_meta.append((SETS.index('ctrl_' + kind), r_i, slot, real)); contrib.append(np.array([j]))
    log(f'  controls: {NCTRL} realisations, {nfail} member draws failed; apertures {len(ap_pix)}')
    grid = L.profile_grid()
    need = np.unique(np.concatenate(contrib).astype(int))
    P1, P2 = L.cluster_profiles(cl['M'][need], cl['sig'][need], cl['z'][need], grid)
    rowof = {c: i for i, c in enumerate(need)}
    omega = 4 * np.pi / hp.nside2npix(L.NSIDE)
    K1 = np.zeros(len(ap_pix)); K2 = np.zeros(len(ap_pix)); nclean = np.zeros(len(ap_pix), int)
    rows_W, cols_W = [], []
    for a, (pix, cs) in enumerate(zip(ap_pix, contrib)):
        pc = pix[b[pix]]; nclean[a] = len(pc)
        rows_W.append(np.full(len(pc), a, np.int32)); cols_W.append(pc.astype(np.int32))
        v = np.stack(hp.pix2vec(L.NSIDE, pc), -1)
        for c in cs:
            th = np.arccos(np.clip(v @ L.unitvec(cl['ra'][c], cl['dec'][c]), -1, 1))
            K1[a] += np.interp(th, L.THETA, P1[rowof[c]]).sum() * omega
            K2[a] += np.interp(th, L.THETA, P2[rowof[c]]).sum() * omega
    W = sparse.csr_matrix((np.full(sum(len(c) for c in cols_W), omega, np.float32),
                           (np.concatenate(rows_W), np.concatenate(cols_W))), shape=(len(ap_pix), hp.nside2npix(L.NSIDE)))
    sparse.save_npz(os.path.join(po, 'W.npz'), W)
    np.savez(os.path.join(po, 'apertures.npz'), meta=np.array(ap_meta), K1=K1, K2=K2, nclean=nclean, sets=np.array(SETS),
             row_N=np.array([r['N'] for r in rows]), row_z=np.array([r['z'] for r in rows]),
             row_ra=np.array([r['ra'] for r in rows]), row_dec=np.array([r['dec'] for r in rows]),
             row_Msum=np.array([cl['M'][r['mem']].sum() for r in rows]), n_clusters=n, n_groups=len(groups),
             n_pairs_all=sum(len(g) == 2 for g in groups), n_rich_all=sum(len(g) >= 3 for g in groups), link_median=np.nanmedian(l),
             iso_n=len(iso_idx))
    s_ = np.array(ap_meta)[:, 0]
    log(f'  wrote {po}: hull model 1h share {K1[s_ == 0].sum() / (K1[s_ == 0] + K2[s_ == 0]).sum():.3f}, circle '
        f'{K1[s_ == 1].sum() / (K1[s_ == 1] + K2[s_ == 1]).sum():.3f}')


# ----------------------------------------------------------------------------------------------------------- integrals
def websky_map():
    """Websky's convergence, R-137's band filter, nside 1024 in the Websky frame (read per patch by bilinear interpolation)."""
    import healpy as hp
    fn = os.path.join(BASE, 'websky_band_1024.npy')
    if os.path.exists(fn):
        return np.load(fn)
    fa = os.path.join(WS, 'kap_alm2048.fits')
    if not os.path.exists(fa):
        k = hp.read_map(os.path.join(WS, 'kap.fits'), dtype=np.float32)
        hp.write_alm(fa, hp.map2alm(k, lmax=L.LMAX, iter=0), overwrite=True); del k
    a = hp.read_alm(fa).astype(np.complex128)
    m = hp.alm2map(hp.almxfl(a, L.band_filter(np.arange(L.LMAX + 1).astype(float))), 1024, lmax=L.LMAX).astype(np.float32)
    np.save(fn, m)
    return m


def integrate(args):
    import healpy as hp
    from scipy import sparse
    m = websky_map()
    npix = hp.nside2npix(L.NSIDE)
    v = np.stack(hp.pix2vec(L.NSIDE, np.arange(npix)), -1)
    Ws = {}
    for p in range(NPATCH):
        po = os.path.join(OUT, f'patch{p:02d}')
        Ws[p] = sparse.load_npz(os.path.join(po, 'W.npz'))
        vs = v @ rotation(p)          # v_obs = R v_sim  ->  v_sim = R^T v_obs (rows: v_obs @ R)
        th, ph = hp.vec2ang(vs)
        mp = hp.get_interp_val(m, th, ph)
        np.save(os.path.join(po, 'K_sky.npy'), Ws[p] @ mp)
        if p == 0:
            fg = L.map_from_alm(L.read_alm(os.path.join(WS, 'kappa_fg_freqcoadd_prh.fits')))
            np.save(os.path.join(po, 'K_fg.npy'), Ws[p] @ fg)
        log(f'patch {p}: noiseless integrals')
    files = sorted(glob.glob(os.path.join(NOISE, 'noise_*.npy')))[:NSIM]
    Kn = {p: [] for p in Ws}; Kr = {p: [] for p in Ws}
    for f in files:
        nm = np.load(f).astype(np.float64); rm = np.load(f.replace('noise_', 'recon_')).astype(np.float64)
        for p, W in Ws.items():
            Kn[p].append(W @ nm); Kr[p].append(W @ rm)
    for p in Ws:
        po = os.path.join(OUT, f'patch{p:02d}')
        np.save(os.path.join(po, 'K_noise.npy'), np.array(Kn[p])); np.save(os.path.join(po, 'K_recon.npy'), np.array(Kr[p]))
    log(f'noise: {len(files)} realisations integrated on {len(Ws)} patches')


# ---------------------------------------------------------------------------------------------------------------- fit
# R-137's measurement (results-R137.md, D4): E, E_control, E - E_control and its simulation error, per sample and aperture
OBS = {('>= 2 members', 'hull'): (0.46, 0.17, 0.29, 0.36), ('>= 2 members', 'circle'): (-0.25, 0.06, -0.31, 0.40),
       ('>= 3 members', 'hull'): (0.30, 0.03, 0.27, 0.49), ('>= 3 members', 'circle'): (-0.18, -0.13, -0.05, 0.56),
       ('pairs', 'hull'): (np.nan, np.nan, 0.32, 0.47), ('pairs', 'circle'): (np.nan, np.nan, -0.55, 0.48)}
DATA_N = {2: 246, 34: 96, 5: 23}
# the 25 patches cover the sky ~3 times over (results-R140.md, D3: union 95 % of the sky, 4,865 deg2 each, mean pairwise overlap
# 12 %): the number of independent patches for the error on the LCDM mean is the sky over the footprint, ~8
N_EFF = 8   # R-137's 365: pairs, 3-4 members, >= 5 members


K2SCALE = float(os.environ.get('SC_NULL_K2SCALE', '1'))  # information: (0.81 / 0.88)^2 puts the model at Websky's sigma8


def stats_many(K, A, rows_sel, w):
    """E and E_control (R-137's definitions, sc_lens.stats) for each row of K (nreal x naperture)."""
    meta = A['meta']; sets = list(A['sets']); Km = A['K1'] + K2SCALE * A['K2']
    K = np.atleast_2d(K); out = {}
    for kind in ('hull', 'circle'):
        s = np.where(meta[:, 0] == sets.index(kind))[0]
        r = meta[s, 1]; use = rows_sel[r]; s = s[use]; ww = w[kind][r[use]]
        out[kind] = (K[:, s] - Km[s][None]) @ ww / (ww * Km[s]).sum()
        sc = meta[:, 0] == sets.index('ctrl_' + kind)
        Ec = []
        for real in np.unique(meta[sc, 3]):
            q = np.where(sc & (meta[:, 3] == real))[0]
            rr = meta[q, 1]; u = rows_sel[rr]; q = q[u]; ww = w[kind][rr[u]]
            Ec.append((K[:, q] - Km[q][None]) @ ww / (ww * Km[q]).sum())
        out['ctrl_' + kind] = np.mean(Ec, 0)
    return out


def fit(args):
    rng = np.random.default_rng(1401)
    res = {}   # (sample, kind) -> dict of arrays over patches
    info = []
    fgline = []
    for p in range(NPATCH):
        po = os.path.join(OUT, f'patch{p:02d}')
        A = dict(np.load(os.path.join(po, 'apertures.npz')))
        meta = A['meta']; sets = list(A['sets']); Km = A['K1'] + K2SCALE * A['K2']
        Ksky = np.load(os.path.join(po, 'K_sky.npy')); Kn = np.load(os.path.join(po, 'K_noise.npy'))
        Kr = np.load(os.path.join(po, 'K_recon.npy'))
        N = A['row_N']; nrow = len(N)
        w = {}
        for kind in ('hull', 'circle'):
            s = meta[:, 0] == sets.index(kind)
            ww = np.zeros(nrow); ww[meta[s, 1]] = Km[s] / Kr[:, s].var(0); w[kind] = ww
        # the multiplicity-matched subsample: R-137's mix of pairs, 3-4 and >= 5 members
        cls = np.where(N == 2, 2, np.where(N >= 5, 5, 34))
        nmax = min((cls == c).sum() / (DATA_N[c] / 365) for c in DATA_N)
        matched = np.zeros(nrow, bool)
        for c in DATA_N:
            ix = np.where(cls == c)[0]; matched[rng.choice(ix, int(round(nmax * DATA_N[c] / 365)), replace=False)] = True
        samples = {'>= 2 members': N >= 2, '>= 3 members': N >= 3, 'pairs': N == 2, 'matched mix': matched}
        info.append((N.size, (N == 2).sum(), (N >= 3).sum(), (N >= 5).sum(), int(A['n_clusters']), int(A['n_groups']), int(A['iso_n'])))
        for name, rs in samples.items():
            d0 = stats_many(Ksky, A, rs, w)
            dn = stats_many(Ksky[None] + Kn, A, rs, w)
            if p == 0:
                df = stats_many(Ksky + np.load(os.path.join(po, 'K_fg.npy')), A, rs, w)
            for kind in ('hull', 'circle'):
                r = res.setdefault((name, kind), dict(E0=[], C0=[], En=[], Cn=[], n=[]))
                r['E0'].append(d0[kind][0]); r['C0'].append(d0['ctrl_' + kind][0])
                r['En'].append(dn[kind]); r['Cn'].append(dn['ctrl_' + kind]); r['n'].append(rs.sum())
                if p == 0:
                    fgline.append(f'patch 0 {name:13s} {kind:6s}: E - E_control noiseless {d0[kind][0] - d0["ctrl_" + kind][0]:+.3f}, '
                                  f'with ACT\'s Websky foreground reconstruction added {df[kind][0] - df["ctrl_" + kind][0]:+.3f} '
                                  f'(foreground shift {df[kind][0] - df["ctrl_" + kind][0] - d0[kind][0] + d0["ctrl_" + kind][0]:+.3f}; '
                                  f'E shift {df[kind][0] - d0[kind][0]:+.3f})')
    lines = []
    info = np.array(info)
    lines.append(f'patches {NPATCH}; per patch (mean, min-max): systems {info[:, 0].mean():.0f} ({info[:, 0].min()}-{info[:, 0].max()}), '
                 f'pairs {info[:, 1].mean():.0f}, >= 3 {info[:, 2].mean():.0f} ({info[:, 2].min()}-{info[:, 2].max()}), '
                 f'>= 5 {info[:, 3].mean():.0f}; clusters in the footprint {info[:, 4].mean():.0f}; FoF systems (all z) '
                 f'{info[:, 5].mean():.0f}; isolated pool {info[:, 6].mean():.0f}  (data: 365, 246, 119, 23; 8,910; 1,338; 1,574)')
    for (name, kind), r in res.items():
        E0, C0 = np.array(r['E0']), np.array(r['C0']); D0 = E0 - C0
        En, Cn = np.array(r['En']), np.array(r['Cn']); Dn = En - Cn
        sp = D0.std(ddof=1); snz = Dn.std(1, ddof=1).mean(); st = Dn.std(ddof=1)
        sm = sp / np.sqrt(N_EFF); sd_ = np.hypot(st, sm)   # the comparison's error: one sample's spread and the mean's
        lo, hi = np.percentile(Dn, [2.5, 97.5])
        o = OBS.get((name, kind))
        line = (f'{name:13s} {kind:6s} n {np.mean(r["n"]):5.0f} | LCDM: E {E0.mean():+.3f} (patches ± {E0.std(ddof=1):.3f}), '
                f'E_control {C0.mean():+.3f} (± {C0.std(ddof=1):.3f}), E - E_control {D0.mean():+.3f} ± {sm:.3f} '
                f'(mean, N_eff {N_EFF}); one data-sized sample: sd {st:.3f} (patches {sp:.3f}, noise {snz:.3f}), 95 % [{lo:+.3f}, {hi:+.3f}]; '
                f'E sd {En.std(ddof=1):.3f}, E_control sd {Cn.std(ddof=1):.3f}')
        if o:
            Eo, Co, Do, so = o
            dl = Do - D0.mean()
            pct = (Dn < Do).mean() * 100
            line += (f' | observed E - E_control {Do:+.2f}: Δ = {dl:+.3f} ± {sd_:.3f} ({dl / sd_:+.1f}σ), percentile {pct:.1f}, '
                     f'95 % [{dl - 1.96 * sd_:+.2f}, {dl + 1.96 * sd_:+.2f}]')
            if np.isfinite(Eo):
                line += (f'; observed E {Eo:+.2f} vs LCDM {E0.mean():+.3f} ± {En.std(ddof=1):.3f} ({(Eo - E0.mean()) / En.std(ddof=1):+.1f}σ), '
                         f'E_control {Co:+.2f} vs {C0.mean():+.3f} ± {Cn.std(ddof=1):.3f} ({(Co - C0.mean()) / Cn.std(ddof=1):+.1f}σ)')
        lines.append(line)
    lines += fgline
    for l_ in lines:
        log(l_)
    tag = '' if K2SCALE == 1 else f'_k2_{K2SCALE:.3f}'
    with open(os.path.join(OUT, f'fit_null{tag}.txt'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    np.savez(os.path.join(OUT, f'fit_null{tag}.npz'), **{f'{n}|{k}|{q}': np.array(v) for (n, k), r in res.items() for q, v in r.items()})


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('cmd'); ap.add_argument('--worker', type=int, default=0)
    ap.add_argument('--nworkers', type=int, default=1); ap.add_argument('--patch', type=int, default=0)
    a = ap.parse_args()
    {'noise': noise, 'halos': halos, 'select': select, 'prep': prep, 'integrate': integrate, 'fit': fit}[a.cmd](a)


if __name__ == '__main__':
    main()
