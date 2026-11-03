"""R-168 (cosmo-cf4-rev): the beta statistic on the mock chain (the calibrated null) and on the data.

modes:
  forecast <mockdir>             Step 0: CF4-like null of the joint beta_rest / beta_dipole on Pantheon+ (published covariance).
  null <variant> <mockdir> <out> Method 1: per realisation, sample and field (truth, CF4-like, 2M++-like) the joint fit, the
                                 paired CF4 - 2M++ difference and the bulk-flow angle. ZTF fit-stage variants: zp01, zp02 (per
                                 ZTF field zero points of 0.01 / 0.02 mag), sel188 (selection m50 = 18.8), mwebv (E(B-V) nuisance).
  data <out.json>                Method 2: observed betas per sample and field, the paired difference (separate fits on the same
                                 SNe; bootstrap over whitened rows; joint fit of both fields), the bulk-flow direction and its
                                 bootstrap error, and data checks (published covariance, E(B-V) and |b| nuisances, ZTF field
                                 jackknife).
Mock residuals: Pantheon+ uses the chain's y (its noise is shared with the CF4-like inputs); DES y = corr(v_true) + L z;
ZTF y = corr(v_true) + sigma z, then the ZTF magnitude selection (logistic in peak m_B, m50 = 19.0, width 0.15) applied by
redrawing the noise of rejected SNe (positions are the real ones), and sigma_int refitted per mock. ZTF fits carry a host-mass
step nuisance (global host mass > 10) in data and mocks.
"""
import sys, os, glob, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cf4rev_common import (D, C, pantheon, des, ztf, beta_fits, corr_mag, whitener, cf4pp_at, twompp_at, split, gls, mu, dc_h,
                           gal_xyz)
from astropy.coordinates import SkyCoord
import astropy.units as u

M50, MW = 19.0, 0.15

def ztf_extras(Z):
    t = pd.read_csv(os.path.join(D, 'ztf', 'ztfsniadr2_lite', 'tables', 'snia_data.csv'))
    gh = pd.read_csv(os.path.join(D, 'ztf', 'ztfsniadr2_lite', 'tables', 'globalhost_data.csv')).set_index('ztfname')
    key = pd.Series(np.arange(len(t)), index=list(zip(np.round(t.ra, 6), np.round(t.dec, 6))))
    rows = key.reindex(list(zip(np.round(Z['ra'], 6), np.round(Z['dec'], 6)))).values.astype(int)
    names = t.ztfname.values[rows]; mass = gh.mass.reindex(names).values
    Z['massstep'] = np.where(np.isnan(mass), 0.5, (mass > 10).astype(float))
    Z['mwebv'] = t.mwebv.values[rows]
    Z['nuis'] = Z['nuis'] + [Z['massstep']]
    fp = os.path.join(D, 'ztf_fields.csv')
    if not os.path.exists(fp):
        o = pd.read_parquet(os.path.join(D, 'ztf', 'ztfsniadr2_lite', 'tables', 'observing_logs.parquet'),
                            columns=['fieldid', 'fieldra', 'fielddec']).drop_duplicates('fieldid'); o.to_csv(fp, index=False)
    o = pd.read_csv(fp); o = o[o.fieldid < 1000]
    j, _, _ = SkyCoord(Z['ra']*u.deg, Z['dec']*u.deg).match_to_catalog_sky(SkyCoord(o.fieldra.values*u.deg, o.fielddec.values*u.deg))
    Z['field'] = o.fieldid.values[j]
    X = np.column_stack(Z['nuis'])*(1/Z['sig'])[:, None]; b, *_ = np.linalg.lstsq(X, Z['y']/Z['sig'], rcond=None)
    Z['r_data'] = Z['y'] - np.column_stack(Z['nuis'])@b
    Z['mB'] = Z['y'] + mu(Z['zc']) + 10.635
    Z['var0'] = Z['sig']**2 - Z['sig_int']**2
    return Z

def sub(s, m, name):
    def cut(v):
        if not isinstance(v, np.ndarray) or v.shape[:1] != m.shape: return v
        return v[np.ix_(m, m)] if v.ndim == 2 and v.shape[1] == len(m) else v[m]
    o = {k: cut(v) for k, v in s.items() if k != 'Lc'}
    if 'cov' in o: o['Lc'] = np.linalg.cholesky(o['cov'])
    o['nuis'] = [c[m] for c in s['nuis']]; o['name'] = name; o['_sel'] = m; return o

SIG_NL = 250.0; SIG_GRP = 150.0
def vmag(zc, v): return 5/np.log(10)*v/(C*zc)

def add_groups(P, Z):
    """CF4 group (1PGC) of each SN's matched host (Pantheon+: 120"/1000 km/s; ZTF: 30"/1000 km/s, 113 SNe), and the shared
    small-scale velocity (150 km/s of the 250) of SNe in the same group added to the covariance; ZTF becomes a full matrix."""
    from cf4rev_chain import load_cf4, ztf_hosts
    t = load_cf4(); gid = t['1PGC'].fillna(t.PGC).values
    hra = np.where(P['host_ra'] > -900, P['host_ra'], P['ra']); hde = np.where(P['host_dec'] > -900, P['host_dec'], P['dec'])
    j, sep, _ = SkyCoord(hra*u.deg, hde*u.deg).match_to_catalog_sky(SkyCoord(t.RAJ2000.values*u.deg, t.DEJ2000.values*u.deg))
    P['host'] = np.where((sep.arcsec < 120) & (np.abs(t.Vcmb.values[j] - C*P['zc']) < 1000), j, -1)
    Z['host'] = ztf_hosts(t, Z)
    for s in (P, Z):
        s['grp'] = np.where(s['host'] >= 0, gid[np.clip(s['host'], 0, None)], -1 - np.arange(len(s['zc'])))
        same = (s['grp'][:, None] == s['grp'][None, :]); np.fill_diagonal(same, False)
        if s is P: same &= np.asarray(P['cid']).astype(str)[:, None] != np.asarray(P['cid']).astype(str)[None, :]
        sg = vmag(s['zc'], SIG_GRP); s['Gcov'] = np.outer(sg, sg)*same
        s['npair'] = int(same.sum()//2)
    P['cov'] = P['cov'] + P['Gcov']; P['Lc'] = np.linalg.cholesky(P['cov'])
    Z['cov'] = np.diag(Z['sig']**2) + Z['Gcov']; Z['Lc'] = np.linalg.cholesky(Z['cov'])

def load_samples(dropm=False):
    P = pantheon(); Z = ztf_extras(ztf()); add_groups(P, Z)
    if dropm:
        Z = sub(Z, Z['host'] < 0, 'ZTF'); Z['cov'] = Z['cov'][np.ix_(Z['_sel'], Z['_sel'])]; Z['Lc'] = np.linalg.cholesky(Z['cov'])
        Z['_dropsel'] = Z.pop('_sel')
    Z6 = sub(Z, Z['zc'] < 0.06, 'ZTF z<0.06'); Z6['cov'] = Z['cov'][np.ix_(Z6['_sel'], Z6['_sel'])]; Z6['Lc'] = np.linalg.cholesky(Z6['cov'])
    return [P, des(), Z, Z6]

def refit_sigint(s, y):
    """sigma_int for chi2/dof = 1 on the base (nuisance) model, as for the data; returns the per-SN sigma."""
    lo, hi = 0.0, 0.5
    for _ in range(30):
        si = (lo+hi)/2; w = 1/np.sqrt(s['var0'] + si**2); X = np.column_stack(s['nuis'])*w[:, None]
        b, *_ = np.linalg.lstsq(X, y*w, rcond=None); r = y*w - X@b
        lo, hi = (si, hi) if r@r > len(y) - X.shape[1] else (lo, si)
    return np.sqrt(s['var0'] + ((lo+hi)/2)**2)

def mock_y(s, rec, rng, fv):
    nm = s['name']
    if nm == 'Pantheon+': return rec['Pantheon+_y'], None
    if nm == 'DES': return corr_mag(s['zc'], rec['DES_true']) + s['Lc'] @ rng.standard_normal(len(s['zc'])), None
    vt = pick(s, rec['ZTF_true']); vnl = pick(s, rec['ZTF_vnl'])
    cv = corr_mag(s['zc'], vt + vnl)               # the chain's small-scale velocity (shared with CF4 hosts and groups)
    sd = np.sqrt(np.clip(s['sig']**2 - vmag(s['zc'], SIG_NL)**2, 1e-6, None))
    y = cv + sd*rng.standard_normal(len(vt))
    if fv in ('zp01', 'zp02'):
        uf, inv = np.unique(s['field'], return_inverse=True); y = y + rng.normal(0, 0.01 if fv == 'zp01' else 0.02, len(uf))[inv]
    m50 = 18.8 if fv == 'sel188' else M50; mw = 0.30 if fv == 'selw30' else MW
    for _ in range(100):   # the magnitude selection: redraw the noise of SNe that would not have been selected
        m = s['mB'] - s['r_data'] + y - np.median(y)
        rej = rng.random(len(y)) > 1/(1 + np.exp((m - m50)/mw))
        if not rej.any(): break
        y[rej] = cv[rej] + sd[rej]*rng.standard_normal(rej.sum())
    sig = refit_sigint(s, y)
    from scipy.linalg import solve_triangular
    Lm = np.linalg.cholesky(np.diag(sig**2) + s['Gcov'])
    return y, (lambda X, Lm=Lm: solve_triangular(Lm, X, lower=True))

def pick(s, v):
    if '_dropsel' in s: v = v[s['_dropsel']]
    return v[s['_sel']] if '_sel' in s else v

def field_of(s, rec, f):
    base = 'ZTF' if s['name'].startswith('ZTF') else s['name']
    return pick(s, rec[f'{base}_{f}'])

def with_positions(s, rec):
    """'clustered' variant: the realisation's own SN positions (redshift, RA, Dec)."""
    base = 'ZTF' if s['name'].startswith('ZTF') else s['name']
    if f'{base}_zc' not in rec: return s
    o = dict(s); o['zc'] = pick(s, rec[f'{base}_zc']); o['ra'] = pick(s, rec[f'{base}_ra']); o['dec'] = pick(s, rec[f'{base}_dec'])
    return o

def bulk_cols(s):
    rad, de = np.radians(s['ra']), np.radians(s['dec'])
    nh = np.stack([np.cos(de)*np.cos(rad), np.cos(de)*np.sin(rad), np.sin(de)], -1)
    return [corr_mag(s['zc'], 100*nh[:, k])/100 for k in range(3)], nh

def angle(a, b): return np.degrees(np.arccos(np.clip(a@b/np.linalg.norm(a)/np.linalg.norm(b), -1, 1)))

def run_null(variant, mockdir, out, fv=None):
    global SIG_NL, SIG_GRP
    SIG_NL = {'vnl150': 150.0, 'vnl400': 400.0}.get(variant, 250.0); SIG_GRP = 0.6*SIG_NL
    S = load_samples(dropm=fv == 'dropm'); rng = np.random.default_rng(5)
    if fv == 'mwebv':
        for s in S[2:]: s['nuis'] = s['nuis'] + [s['mwebv']]
    W = {s['name']: whitener(s) for s in S}; BC = {s['name']: bulk_cols(s) for s in S}
    files = sorted(glob.glob(os.path.join(mockdir, 'real_*.npz'))); res = {}
    for fn in files:
        rec = np.load(fn); fields = [f for f in ('true', 'cf4', '2mpp') if f'Pantheon+_{f}' in rec]
        for s0 in S:
            s = with_positions(s0, rec)
            y, Wm = mock_y(s, rec, rng, fv); Wm = Wm or W[s['name']]
            bb = {}
            for f in fields:
                r = beta_fits(s, field_of(s, rec, f), y, Wm); bb[f] = r
                for k in ('bd', 'br', 'bw', 'corr_rest'): res.setdefault(f"{s['name']}|{f}|{k}", []).append(r[k])
            if 'cf4' in bb and '2mpp' in bb:
                for k in ('bd', 'br'): res.setdefault(f"{s['name']}|diff|{k}", []).append(bb['cf4'][k] - bb['2mpp'][k])
            if 'cf4' in bb:
                b, _, _, _ = gls(s, bulk_cols(s)[0] if s is not s0 else BC[s['name']][0], y, Wm)
                res.setdefault(f"{s['name']}|angle_cf4|deg", []).append(angle(b, bb['cf4']['B']))
                res.setdefault(f"{s['name']}|angle_true|deg", []).append(angle(b, split(s['ra'], s['dec'], s['zc'], field_of(s, rec, 'true'))[2]))
    res = {k: np.array(v) for k, v in res.items()}
    np.savez(out, **res)
    print(f"[{variant}{'/' + fv if fv else ''}] {len(files)} realisations")
    for k, v in res.items(): print(f"  {k:28s} mean {v.mean():7.3f} sd {v.std(ddof=1):6.3f} median {np.median(v):7.3f}")

def wfreal_at(fn, s):
    """the plan's Wiener filter applied to the real CF4 catalogue (cf4rev_realwf.py), read like CF4++."""
    from scipy.ndimage import map_coordinates
    from cf4rev_common import sg_xyz
    V = np.load(fn)['v']; N = V.shape[1]; xyz = sg_xyz(s['ra'], s['dec'], dc_h(s['zc'])); ci = (xyz + 500.0)/(1000.0/N) - 0.5
    return np.sum(np.stack([map_coordinates(V[k], ci, order=1, mode='nearest') for k in range(3)])*xyz/np.linalg.norm(xyz, axis=0), 0)

def run_data(out):
    S = load_samples(); R = {}
    for s in S:
        nm = s['name']; W = whitener(s)
        fld = {'cf4': cf4pp_at(s['ra'], s['dec'], s['zc']), '2mpp': twompp_at(s['ra'], s['dec'], s['zc'], order='ijk')}
        wfr = os.path.join(D, 'cf4rev', 'wf_real.npz')
        if os.path.exists(wfr): fld['wfreal'] = wfreal_at(wfr, s)
        R[nm] = {f: {k: (v.tolist() if isinstance(v, np.ndarray) else float(v)) for k, v in beta_fits(s, fld[f]).items()} for f in fld}
        # paired difference: joint fit of both fields' dipole and rest columns on the same SNe
        cols = []
        for f in ('cf4', '2mpp'):
            d_, r_, _ = split(s['ra'], s['dec'], s['zc'], fld[f]); cols += [corr_mag(s['zc'], d_), corr_mag(s['zc'], r_)]
        b, e, _, cv = gls(s, cols)
        R[nm]['joint4'] = dict(b=b.tolist(), e=e.tolist(), d_rest=float(b[1] - b[3]), e_d_rest=float(np.sqrt(cv[1, 1] + cv[3, 3] - 2*cv[1, 3])),
                               d_dip=float(b[0] - b[2]), e_d_dip=float(np.sqrt(cv[0, 0] + cv[2, 2] - 2*cv[0, 2])))
        # bootstrap over whitened rows: both fields refitted on each draw; and the free bulk flow
        Xn = np.column_stack(s['nuis']); k = Xn.shape[1]
        bc, nh = bulk_cols(s)
        def wcols(c): return W(np.column_stack(c))
        Xc = W(np.column_stack([Xn] + cols[:2])); X2 = W(np.column_stack([Xn] + cols[2:])); XB = W(np.column_stack([Xn] + bc)); yw = W(s['y'])
        Bhat = np.linalg.lstsq(XB, yw, rcond=None)[0][k:]
        cf4B = split(s['ra'], s['dec'], s['zc'], fld['cf4'])[2]
        rng = np.random.default_rng(7); n = len(yw); dr, dd, ang, Bs = [], [], [], []
        for _ in range(1000):
            i = rng.integers(0, n, n)
            bc_ = np.linalg.lstsq(Xc[i], yw[i], rcond=None)[0][k:]; b2_ = np.linalg.lstsq(X2[i], yw[i], rcond=None)[0][k:]
            dd.append(bc_[0] - b2_[0]); dr.append(bc_[1] - b2_[1])
            bb = np.linalg.lstsq(XB[i], yw[i], rcond=None)[0][k:]; Bs.append(bb); ang.append(angle(bb, cf4B))
        Bs = np.array(Bs)
        gl = lambda v: SkyCoord(x=v[0], y=v[1], z=v[2], representation_type='cartesian', frame='icrs').galactic
        g0 = gl(Bhat); offs = [angle(bb, Bhat) for bb in Bs]
        R[nm]['boot'] = dict(d_rest=float(np.mean(dr)), sd_d_rest=float(np.std(dr)), d_dip=float(np.mean(dd)), sd_d_dip=float(np.std(dd)),
                             B=float(np.linalg.norm(Bhat)), sd_B=float(np.std(np.linalg.norm(Bs, axis=1))),
                             l=float(g0.l.deg), b=float(g0.b.deg), dir_err68=float(np.percentile(offs, 68)),
                             angle_to_cf4=float(angle(Bhat, cf4B)), angle_sd=float(np.std(ang)), angle_p16_84=np.percentile(ang, [16, 84]).tolist(),
                             cf4_B=float(np.linalg.norm(cf4B)), cf4_l=float(gl(cf4B).l.deg), cf4_b=float(gl(cf4B).b.deg))
        # data checks
        chk = {}
        if nm == 'Pantheon+':
            P0 = pantheon(drop_vpec_sys=False)
            chk['published_cov'] = {f: {k: float(v) for k, v in beta_fits(P0, fld[f]).items() if k != 'B'} for f in fld}
            chk['VPEC_column_2mpp'] = {k: float(v) for k, v in beta_fits(s, s['vpec']).items() if k != 'B'}
        if nm.startswith('ZTF'):
            s2 = dict(s); s2['nuis'] = s['nuis'] + [s['mwebv']]
            chk['mwebv_nuisance'] = {f: {k: float(v) for k, v in beta_fits(s2, fld[f]).items() if k != 'B'} for f in fld}
            gb = SkyCoord(s['ra']*u.deg, s['dec']*u.deg).galactic.b.deg
            s3 = dict(s); s3['nuis'] = s['nuis'] + [np.abs(gb)/90]
            chk['absb_nuisance'] = {f: {k: float(v) for k, v in beta_fits(s3, fld[f]).items() if k != 'B'} for f in fld}
            # delete-one-group jackknife over ZTF fields (fields grouped into 30 groups by field id order)
            uf = np.unique(s['field']); grp = {f_: i % 30 for i, f_ in enumerate(uf)}; g_ = np.array([grp[f_] for f_ in s['field']])
            jk = {f: [] for f in fld}
            for gi in range(30):
                m = g_ != gi; ss = sub(s, m, nm)
                for f in fld: jk[f].append(beta_fits(ss, fld[f][m])['br'])
            chk['field_jackknife_rest_err'] = {f: float(np.sqrt(29/30*np.sum((np.array(v) - np.mean(v))**2))) for f, v in jk.items()}
        R[nm]['checks'] = chk; R[nm]['n'] = int(len(s['y']))
        print(nm, json.dumps({f: {k: round(R[nm][f][k], 3) for k in ('bd', 'ed', 'br', 'er', 'bw', 'ew')} for f in fld}))
        print('   joint4', {k: (np.round(v, 3).tolist() if isinstance(v, list) else round(v, 3)) for k, v in R[nm]['joint4'].items()})
        print('   boot', {k: (np.round(v, 2).tolist() if isinstance(v, list) else round(v, 2)) for k, v in R[nm]['boot'].items()})
        print('   checks', json.dumps(chk)[:900])
    json.dump(R, open(out, 'w'), indent=1)

def run_forecast(mockdir):
    S = load_samples(); P = S[0]; rng = np.random.default_rng(5); W = whitener(P); br, bd = [], []
    for fn in sorted(glob.glob(os.path.join(mockdir, 'real_*.npz'))):
        rec = np.load(fn); r = beta_fits(P, rec['Pantheon+_cf4'], rec['Pantheon+_y'], W); br.append(r['br']); bd.append(r['bd'])
    print(f"forecast: beta_rest null {np.mean(br):.2f} +- {np.std(br):.2f}; beta_dipole {np.mean(bd):.2f} +- {np.std(bd):.2f}")

if __name__ == '__main__':
    m = sys.argv[1]
    if m == 'forecast': run_forecast(sys.argv[2])
    elif m == 'null': run_null(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5] if len(sys.argv) > 5 else None)
    elif m == 'data': run_data(sys.argv[2])
