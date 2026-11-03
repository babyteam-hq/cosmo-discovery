"""R-168 (cosmo-cf4-rev): the null chain's Wiener filter applied to the real Cosmicflows-4 catalogue (method check pass 2,
point 1). Same grid (128^3, 1000 Mpc/h, supergalactic), prior (EH no-wiggle, sigma8 0.81) and noise model as the baseline
mocks (sd = (ln10/5) e_DM Vcmb plus 250 km/s), data at r = Vcmb/100 Mpc/h. Velocity estimator: Watkins & Feldman (2015),
u = Vcmb ln(Vcmb / (H0 d)), H0 d from the catalogue's DM at its own H0 = 74.6, which has near-Gaussian errors in log distance;
a global H0 rescale (u -> u - a Vcmb, a fitted) is removed first. Writes data/cf4rev/wf_real.npz (the 3-velocity grid) and
prints: the fitted a; the chi2 per datum of the data about the WF prediction against the same statistic in mock realisations;
and the field-level comparison with CF4++ at the SN positions (rms, correlation) and in shells.
usage: python cf4rev_realwf.py [n_mock_chi2]
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cf4rev_chain as ch
from cf4rev_common import D, cf4pp_at, pantheon, ztf, dc_h, sg_xyz

if __name__ == '__main__':
    nmock = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    t = ch.load_cf4(); gx = sg_xyz(t.RAJ2000.values, t.DEJ2000.values, t.Vcmb.values/ch.H)
    V = t.Vcmb.values; d = 10**((t.DM.values - 25)/5)
    u = V*np.log(V/(74.6*d))
    a = np.sum(u*V)/np.sum(V*V); u = u - a*V
    sig = np.sqrt((np.log(10)/5*t.e_DM.values*V)**2 + ch.SIG_NL**2)
    gW = ch.Grid(128, ch.make_pk()); wf = ch.WF(gW, gx, sig)
    Vwf, nit = wf.solve(u); pred = gW.at_r(Vwf, gx)
    c2 = np.mean(((u - pred)/sig)**2); c2_0 = np.mean((u/sig)**2)
    print(f"real CF4: {len(u)} galaxies; H0 rescale a = {a:+.4f} (u -= a Vcmb); CG {nit}; chi2/datum about the WF {c2:.3f} (about zero {c2_0:.3f})")
    os.makedirs(os.path.join(D, 'cf4rev'), exist_ok=True)
    np.savez(os.path.join(D, 'cf4rev', 'wf_real.npz'), v=Vwf.astype(np.float32), a=a, chi2=c2)
    # the same statistic on mocks of the baseline chain (Gaussian errors at the stated sd) and with true errors 1.3x
    S = ch.setup('noshare'); rng = np.random.default_rng(99)
    for lab, fac in (('stated errors', 1.0), ('true errors 1.3x', 1.3)):
        vals = []
        for i in range(nmock):
            gT = S['gT']; dk = gT.rf(rng.standard_normal((gT.N,)*3))*gT.sqp; Vt = gT.vel(dk); del dk
            vtr = gT.at_r(Vt, gx) + rng.normal(0, ch.SIG_NL, len(u))
            um = vtr + rng.standard_normal(len(u))*S['gsd']*fac
            Vm, _ = wf.solve(um); vals.append(np.mean(((um - gW.at_r(Vm, gx))/sig)**2))
        print(f"  mocks ({lab}): chi2/datum {np.mean(vals):.3f} +- {np.std(vals):.3f} (n={nmock})")
    # field-level comparison with CF4++
    for s in (pantheon(), ztf()):
        xyz = sg_xyz(s['ra'], s['dec'], dc_h(s['zc'])); vr = gW.at_r(Vwf, xyz); vc = cf4pp_at(s['ra'], s['dec'], s['zc'])
        print(f"  {s['name']}: rms radial v at the SNe: WF(real CF4) {vr.std():.0f}, CF4++ {vc.std():.0f} km/s; correlation {np.corrcoef(vr, vc)[0, 1]:.2f}")
    ax = gW.ax; X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij'); R = np.sqrt(X**2 + Y**2 + Z**2)
    C4 = np.load(os.path.join(D, 'new_CF4pp_mean_std_grids.npz'))['v_mean']
    vr_wf = (Vwf[0]*X + Vwf[1]*Y + Vwf[2]*Z)/np.maximum(R, 1e-9); vr_c4 = (C4[0]*X + C4[1]*Y + C4[2]*Z)/np.maximum(R, 1e-9)
    for lo, hi in ((0, 50), (50, 100), (100, 150), (150, 200), (200, 300)):
        m = (R >= lo) & (R < hi)
        print(f"  shell {lo}-{hi} Mpc/h: rms v_r WF(real) {vr_wf[m].std():.0f}, CF4++ {vr_c4[m].std():.0f}; correlation {np.corrcoef(vr_wf[m], vr_c4[m])[0, 1]:.2f}")
