"""Answers to the blind referee on the CF4 result: (1) mocks: inject the CF4++ field at beta=1 into Pantheon+-like residuals
drawn from the full covariance, run the same decomposition and GLS fit, and see what beta_dipole / beta_rest come back;
(2) attenuation: the errors-in-variables correction from the CF4++ per-voxel posterior std; (3) position error: evaluate
the field at positions displaced along the line of sight by the field's own velocity (and by 300 km/s random), and refit.
usage: python cf4_mocks.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/new_CF4pp_mean_std_grids.npz
"""
import sys, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.ndimage import map_coordinates
from astropy.coordinates import SkyCoord
import astropy.units as u
C = 299792.458; OM = 0.334; N = 128; L = 1000.0; DELTA = L/N
rng = np.random.default_rng(11)
def dc_h(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3 + 1-OM)
    return C/100*np.array([quad(e, 0, zi)[0] for zi in np.atleast_1d(z)])
def mu(z): return 5*np.log10((1+z)*dc_h(z)/0.7) + 25
def zhd(zcmb, v): return (1+zcmb)/(1+v/C) - 1

def main(dat, covf, grid):
    g = np.load(grid); V = g['v_mean']; VS = g['v_std']
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f:
        n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    m = (d.IS_CALIBRATOR == 0) & (d.zCMB > 0.01) & (d.zCMB < 0.1); idx = np.where(m)[0]
    zc = d.zCMB.values[idx]; r = dc_h(zc)
    ra, dec = d.RA.values[idx], d.DEC.values[idx]
    def field_at(rr):
        sc = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, distance=rr*u.Mpc, frame='fk5', equinox='J2000').supergalactic
        xyz = np.stack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value]); rhat = xyz/np.linalg.norm(xyz, axis=0)
        ci = (xyz + L/2)/DELTA - 0.5
        v = np.sum(np.stack([map_coordinates(V[k], ci, order=1, mode='nearest') for k in range(3)])*rhat, 0)
        s = np.sqrt(np.sum(np.stack([map_coordinates(VS[k], ci, order=1, mode='nearest') for k in range(3)])**2*rhat**2, 0))
        return v, s, rhat
    vcf, vstd, rhat = field_at(r)
    rad, de = np.radians(ra), np.radians(dec); nh = np.stack([np.cos(de)*np.cos(rad), np.cos(de)*np.sin(rad), np.sin(de)], -1)
    A = np.column_stack([r, nh])
    def split(v):
        coef, *_ = np.linalg.lstsq(A, v, rcond=None); return A@np.r_[0, coef[1:]], v - A@coef
    Cm = cov[np.ix_(idx, idx)]; Lc = np.linalg.cholesky(Cm); w = lambda X: np.linalg.solve(Lc, X)
    def beta(y, v):
        cc = mu(zhd(zc, v)) - mu(zc); X = np.column_stack([np.ones(len(y)), cc]); Xw, yw = w(X), w(y)
        b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); e = np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))); return b[1], e[1]
    dip, rest = split(vcf)
    y = d.m_b_corr.values[idx] - mu(zc)
    print(f"data: beta_dipole {beta(y, dip)[0]:.2f}±{beta(y, dip)[1]:.2f}, beta_rest {beta(y, rest)[0]:.2f}±{beta(y, rest)[1]:.2f}; rms rest {np.sqrt(np.mean(rest**2)):.0f}, mean field std {vstd.mean():.0f} km/s")
    # (1) mocks: y = c_field(beta=1) + N(0, Cov)
    cfull = mu(zhd(zc, vcf)) - mu(zc)
    bd, br = [], []
    for _ in range(200):
        ym = cfull + Lc @ rng.standard_normal(len(idx))
        bd.append(beta(ym, dip)[0]); br.append(beta(ym, rest)[0])
    print(f"mocks (true field injected at beta=1, noise from the Pantheon+ covariance, 200 draws): beta_dipole {np.mean(bd):.2f}±{np.std(bd):.2f}, beta_rest {np.mean(br):.2f}±{np.std(br):.2f}")
    # (2) attenuation from the posterior std: lambda = var(signal)/(var(signal)+var(err)) for the rest component
    lam = np.mean(rest**2)/(np.mean(rest**2) + np.mean(vstd**2))
    b, e = beta(y, rest); print(f"attenuation factor from the CF4++ posterior std (rest): lambda = {lam:.2f}; errors-in-variables beta_rest = {b/lam:.2f}±{e/lam:.2f}")
    # also: mocks with predictor noise at the posterior std: y from (field + noise), fitted with field
    br2 = []
    for _ in range(100):
        vtrue = vcf + vstd*rng.standard_normal(len(idx)); _, rtrue = split(vtrue)
        ym = (mu(zhd(zc, vtrue)) - mu(zc)) + Lc @ rng.standard_normal(len(idx)); br2.append(beta(ym, rest)[0])
    print(f"mocks where the truth = field + N(0, posterior std) and we fit with the field: beta_rest {np.mean(br2):.2f}±{np.std(br2):.2f}")
    # (3) position error: evaluate the field at r shifted by the field's own velocity (real-space position), and by a random 300 km/s
    v_real, _, _ = field_at(np.clip(r - vcf/100, 1, None)); _, rest_real = split(v_real)
    print(f"field read at real-space positions (r - v_cf/H0): corr with redshift-space read {np.corrcoef(vcf, v_real)[0,1]:.3f}; beta_rest {beta(y, rest_real)[0]:.2f}±{beta(y, rest_real)[1]:.2f}")
    br3 = []
    for _ in range(20):
        v_p, _, _ = field_at(np.clip(r + rng.normal(0, 3.0, len(idx)), 1, None)); _, rp = split(v_p); br3.append(beta(y, rp)[0])
    print(f"field read at positions jittered by 3 Mpc/h (300 km/s): beta_rest {np.mean(br3):.2f}±{np.std(br3):.2f}")
    # how much does a 3 Mpc/h jitter change the field itself (its own coherence)?
    v_p, _, _ = field_at(np.clip(r + rng.normal(0, 3.0, len(idx)), 1, None)); _, rp = split(v_p)
    print(f"  (rest component correlation between jittered and unjittered reads: {np.corrcoef(rest, rp)[0,1]:.3f})")
    # mocks with jitter: truth at real positions, fit at jittered: does beta drop to 0.1?
    br4 = []
    for _ in range(50):
        v_t, _, _ = field_at(np.clip(r + rng.normal(0, 3.0, len(idx)), 1, None))
        ym = (mu(zhd(zc, v_t)) - mu(zc)) + Lc @ rng.standard_normal(len(idx)); br4.append(beta(ym, rest)[0])
    print(f"mocks with a 3 Mpc/h position error between truth and read: beta_rest {np.mean(br4):.2f}±{np.std(br4):.2f}")

main(*sys.argv[1:4])
