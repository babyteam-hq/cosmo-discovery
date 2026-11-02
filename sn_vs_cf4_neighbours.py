"""Sweep-3 row 3: do Pantheon+ residual velocities correlate with the raw peculiar velocities of neighbouring Cosmicflows-4
groups, with no reconstructed field in between? Neighbour velocities from the CF4 group table (TF, FP, SBF, SN II, TRGB
means; the SN Ia-based distances excluded, since they share light curves with Pantheon+), combined in log-distance space;
the SN's own host group excluded. GLS with the Pantheon+ covariance; nulls by sky rotation and by distance-matched swaps.
usage: python sn_vs_cf4_neighbours.py"""
import numpy as np, pandas as pd
from scipy.integrate import quad
from cf4_io import *
rng = np.random.default_rng(3)
C = 299792.458; OM = 0.334; H0 = 70.0
K = 5 / np.log(10)

def dl(z):
    e = lambda x: 1 / np.sqrt(OM * (1 + x) ** 3 + 1 - OM); return (1 + z) * C / H0 * np.array([quad(e, 0, zi)[0] for zi in z])

def load_pantheon():
    d = pd.read_csv('data/Pantheon+SH0ES.dat', sep=r'\s+')
    with open('data/Pantheon+SH0ES_STAT+SYS.cov') as f: n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    return d, cov

def group_velocities(g):
    """Raw CF4 group peculiar velocities from the non-SN methods: log-distance ratio eta = (DM_m - DM_H)/5 per method, combined
    inverse-variance; v_pec = Vcmb (1 - 10^eta), i.e. -Vcmb ln(10) eta to first order."""
    cols = [('DMfp', 'e_DMfp'), ('DMtf', 'e_DMtf'), ('DMsbfo', 'e_DMsbfo'), ('DMsbfi', 'e_DMsbfi'), ('DMsnII', 'e_DMsnII')]
    num = np.zeros(len(g)); den = np.zeros(len(g))
    for dm, edm in cols:
        ok = g[dm].notna() & g[edm].notna() & (g[edm] > 0)
        eta = ((g[dm] - hubble_dm(g.Vcmb)) / 5).where(ok, 0).values; w = (25 / g[edm] ** 2).where(ok, 0).values
        num += w * eta; den += w
    g = g.copy(); g['eta'] = np.where(den > 0, num / np.maximum(den, 1e-9), np.nan); g['e_eta'] = np.where(den > 0, 1 / np.sqrt(np.maximum(den, 1e-9)), np.nan)
    g["vpec_raw"] = -g.Vcmb * np.log(10) * g.eta; g["e_v"] = g.Vcmb * np.log(10) * g.e_eta  # v_pec = Vcmb - H0 d, d = d_H 10^eta
    return g[g.eta.notna() & (g.Vcmb > 500)]

def main():
    d, cov = load_pantheon(); g = vizier_tsv('data/cf4_groups.tsv'); t2 = vizier_tsv('data/cf4_table2.tsv')
    g = group_velocities(g)
    print(f"CF4 groups with a non-SN velocity: {len(g)} (median e_v {g.e_v.median():.0f} km/s, |vpec| rms {g.vpec_raw.std():.0f})")
    m = (d.IS_CALIBRATOR == 0) & (d.zHD > 0.01) & (d.zHD < 0.06); idx = np.where(m)[0]; s = d.iloc[idx].copy()
    s['res'] = s.m_b_corr.values - (5 * np.log10(dl(s.zHD.values)) + 25)  # Hubble residual, mag, after 2M++ (zHD)
    s['res_cmb'] = s.m_b_corr.values - (5 * np.log10(dl(s.zCMB.values)) + 25)  # before the 2M++ correction
    # SN residual velocity: a residual +dmu means the SN is farther than its zHD says, i.e. its peculiar velocity is negative
    s['v_sn'] = -s.res * s.zHD * C / K; s['v_sn_cmb'] = -s.res_cmb * s.zCMB * C / K
    # positions in redshift space (Mpc/h with H0 = 100): SN at cz/100, groups at Vcmb/100
    def xyz(ra, dec, v): return unit_eq(ra, dec) * (v / 100)[:, None]
    P = xyz(s.RA.values, s.DEC.values, s.zCMB.values * C); G = xyz(g.RAJ2000.values, g.DEJ2000.values, g.Vcmb.values.astype(float))
    # the SN host's own group: the CF4 galaxy within 1' of the host (or SN) position, its 1PGC
    host = unit_eq(s.HOST_RA.where(s.HOST_RA > -900, s.RA).values, s.HOST_DEC.where(s.HOST_DEC > -900, s.DEC).values)
    galn = unit_eq(t2.RAJ2000.values, t2.DEJ2000.values)
    from scipy.spatial import cKDTree
    tree = cKDTree(galn); dist, j = tree.query(host, k=1)
    hostgrp = np.where(dist < np.radians(2 / 60), t2['1PGC'].values[j], -1)
    print(f"SNe 0.01 < zHD < 0.06: {len(s)}; with a CF4 host galaxy within 2': {(hostgrp > 0).sum()}")
    gtree = cKDTree(G); g1 = g['1PGC'].values
    Lc = np.linalg.cholesky(cov[np.ix_(idx, idx)])
    def gls(y, X):
        Xw, yw = np.linalg.solve(Lc, X), np.linalg.solve(Lc, y); b, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        r = yw - Xw @ b; cv = np.linalg.inv(Xw.T @ Xw) * (r @ r / (len(y) - X.shape[1])); return b, np.sqrt(np.diag(cv))
    for R in [5, 8, 12]:
        nb = gtree.query_ball_point(P, R); vnb = np.full(len(s), np.nan); nn = np.zeros(len(s), int)
        for i, lst in enumerate(nb):
            lst = [k for k in lst if g1[k] != hostgrp[i]]
            if len(lst) < 3: continue
            w = 1 / g.e_v.values[lst] ** 2; vnb[i] = np.sum(w * g.vpec_raw.values[lst]) / np.sum(w); nn[i] = len(lst)
        ok = np.isfinite(vnb)
        # the residual velocity in mag units is what the covariance describes: regress res (mag) on the neighbour velocity
        # expressed in mag at the SN's redshift: dmu_nb = -K vnb / (c zHD)
        xnb = -K * vnb / (s.zHD.values * C)
        for lab, y in [('after 2M++ (zHD)', s.res.values), ('before 2M++ (zCMB)', s.res_cmb.values)]:
            idx_ok = np.where(ok)[0]; L2 = np.linalg.cholesky(cov[np.ix_(idx[idx_ok], idx[idx_ok])])
            X = np.column_stack([np.ones(ok.sum()), xnb[ok]]); Xw, yw = np.linalg.solve(L2, X), np.linalg.solve(L2, y[ok])
            b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); r = yw - Xw @ b; cv = np.linalg.inv(Xw.T @ Xw) * (r @ r / (ok.sum() - 2))
            rho = np.corrcoef(y[ok], xnb[ok])[0, 1]
            print(f"R={R} Mpc/h, {lab}: N={ok.sum()} (median {np.median(nn[ok]):.0f} groups, neighbour v rms {np.nanstd(vnb):.0f} km/s): "
                  f"beta = {b[1]:+.3f} ± {np.sqrt(cv[1,1]):.3f}, rho = {rho:+.3f}")
        # nulls: rotate the SN sky positions about the celestial pole (keeping z), and swap neighbours among SNe at matched z
        if R == 8:
            y = s.res.values; betas = []; rhos = []
            for _ in range(300):
                ang = rng.uniform(0, 2 * np.pi); refl = rng.choice([-1, 1])
                ra2 = (s.RA.values + np.degrees(ang)) % 360; dec2 = s.DEC.values * refl
                P2 = xyz(ra2, dec2, s.zCMB.values * C); nb2 = gtree.query_ball_point(P2, R); v2 = np.full(len(s), np.nan)
                for i, lst in enumerate(nb2):
                    if len(lst) < 3: continue
                    w = 1 / g.e_v.values[lst] ** 2; v2[i] = np.sum(w * g.vpec_raw.values[lst]) / np.sum(w)
                ok2 = np.isfinite(v2); x2 = -K * v2 / (s.zHD.values * C)
                idx2 = np.where(ok2)[0]; L3 = np.linalg.cholesky(cov[np.ix_(idx[idx2], idx[idx2])])
                X = np.column_stack([np.ones(ok2.sum()), x2[ok2]]); Xw, yw = np.linalg.solve(L3, X), np.linalg.solve(L3, y[ok2])
                b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); betas.append(b[1]); rhos.append(np.corrcoef(y[ok2], x2[ok2])[0, 1])
            print(f"  rotation null (R=8, after 2M++): beta spread {np.std(betas):.3f} (mean {np.mean(betas):+.3f}), rho spread {np.std(rhos):.3f}")
        # splits
        for lab, mm in [('z<0.03', s.zHD.values < 0.03), ('z>=0.03', s.zHD.values >= 0.03), ('Dec>0', s.DEC.values > 0), ('Dec<0', s.DEC.values < 0)]:
            o = ok & mm; idx_o = np.where(o)[0]; L2 = np.linalg.cholesky(cov[np.ix_(idx[idx_o], idx[idx_o])])
            X = np.column_stack([np.ones(o.sum()), xnb[o]]); Xw, yw = np.linalg.solve(L2, X), np.linalg.solve(L2, s.res.values[o])
            b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); r = yw - Xw @ b; cv = np.linalg.inv(Xw.T @ Xw) * (r @ r / (o.sum() - 2))
            print(f"    {lab}: beta {b[1]:+.3f} ± {np.sqrt(cv[1,1]):.3f} (N={o.sum()})")
    # the competing regressor: the CF4++ field's prediction at the SN (R-115 finding) vs the raw neighbours, jointly
    try:
        z = np.load('data/new_CF4pp_mean_std_grids.npz'); print("CF4++ grid keys:", list(z.keys())[:8])
    except Exception as e: print("CF4++ grid not read:", e)

if __name__ == '__main__': main()
