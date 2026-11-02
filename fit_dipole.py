"""Dipole in Pantheon+ Hubble residuals: bulk flow (1/z) vs anisotropic expansion (const).

usage: python -I fit_dipole.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov
"""
import sys, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.optimize import minimize

C = 299792.458
OM = 0.334

def dl_mpc(z, om=OM, h0=70.0):
    e = lambda x: 1/np.sqrt(om*(1+x)**3 + 1-om)
    return (1+z) * C/h0 * np.array([quad(e, 0, zi)[0] for zi in np.atleast_1d(z)])

def unit(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)], -1)

def radec_of(v):
    v = v/np.linalg.norm(v)
    return np.degrees(np.arctan2(v[1], v[0])) % 360, np.degrees(np.arcsin(v[2]))

def eq_to_gal(ra, dec):
    # J2000 -> galactic
    n = unit(ra, dec)
    R = np.array([[-0.0548755604, -0.8734370902, -0.4838350155],
                  [ 0.4941094279, -0.4448296300,  0.7469822445],
                  [-0.8676661490, -0.1980763734,  0.4559837762]])
    g = R @ n
    return radec_of(g)

def main(dat, covf):
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f:
        n = int(f.readline())
        cov = np.loadtxt(f).reshape(n, n)
    assert n == len(d)
    for zcol in ['zHD', 'zCMB']:
        m = (d.IS_CALIBRATOR == 0) & (d[zcol] > 0.01)
        idx = np.where(m)[0]
        z = d[zcol].values[idx]; mb = d.m_b_corr.values[idx]
        nhat = unit(d.RA.values[idx], d.DEC.values[idx])
        Cm = cov[np.ix_(idx, idx)]
        # add zhel uncertainty? already in cov. invert
        L = np.linalg.cholesky(Cm)
        def w(X): return np.linalg.solve(L, X)
        mu = 5*np.log10(dl_mpc(z)) + 25
        y = mb - mu
        # constant 1/z-scaled design: bulk flow dmu = -(5/ln10) (v.n)/(c z) for low z (use exact-ish: (1+z)^2/(H d_L) factor ~ 1/(cz) at low z)
        f_bulk = (5/np.log(10)) / (C*z)          # per km/s
        f_const = (5/np.log(10)) * np.ones_like(z)  # per unit dH/H
        def fit(cols):
            X = np.column_stack([np.ones_like(z)] + cols)
            Xw, yw = w(X), w(y)
            beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
            r = yw - Xw@beta
            chi2 = r@r
            covb = np.linalg.inv(Xw.T@Xw)
            return beta, np.sqrt(np.diag(covb)), chi2
        print(f"\n=== {zcol}  N={len(z)} ===")
        b0, e0, c0 = fit([])
        print(f"no dipole: chi2={c0:.1f}")
        for name, f in [('bulk 1/z [km/s]', f_bulk), ('const [dH/H]', f_const)]:
            cols = [-f*nhat[:, k] for k in range(3)]
            b, e, c = fit(cols)
            v = b[1:]; ev = e[1:]
            amp = np.linalg.norm(v)
            ra, dec = radec_of(v); l, bg = eq_to_gal(ra, dec)
            print(f"{name}: dchi2={c0-c:.2f} (3 dof)  amp={amp:.4g} ± ~{np.mean(ev):.2g}  "
                  f"dir RA={ra:.0f} Dec={dec:.0f} (l={l:.0f} b={bg:.0f})")
        cols = [-f_bulk*nhat[:, k] for k in range(3)] + [-f_const*nhat[:, k] for k in range(3)]
        b, e, c = fit(cols)
        print(f"both: dchi2={c0-c:.2f} (6 dof) bulk amp={np.linalg.norm(b[1:4]):.4g} const amp={np.linalg.norm(b[4:7]):.4g}")
        # shells: free dipole in mag per shell
        print("shells (free dipole, mag):")
        for lo, hi in [(0.01,0.02),(0.02,0.04),(0.04,0.08),(0.08,0.15),(0.15,0.3),(0.3,0.6),(0.6,2.5)]:
            s = (z>=lo)&(z<hi)
            if s.sum() < 20: continue
            cols = [np.where(s, -f_const*nhat[:, k], 0.0) for k in range(3)]
            b, e, c = fit(cols)
            v = b[1:]; amp = np.linalg.norm(v)
            ra, dec = radec_of(v); l, bg = eq_to_gal(ra, dec)
            vkms = amp*C*np.sqrt(np.mean(z[s]**2)) if hi < 0.2 else float('nan')
            print(f"  z[{lo},{hi}) n={s.sum():4d} dchi2={c0-c:6.2f} amp={amp:.4f} ±{np.mean(e[1:]):.4f} mag-ish  l={l:.0f} b={bg:.0f}  ~v={vkms:.0f} km/s")

main(*sys.argv[1:])
