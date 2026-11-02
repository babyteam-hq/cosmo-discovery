"""Second blind-review response on the CF4 result.
(1) Significance with spatially correlated residuals: delete-one-sky-patch jackknife of beta_rest (24 patches), vs the GLS error.
(2) Reconciling with Stiskalek+2025 (whole-field beta ~ 1 with a free velocity scatter sigma_v and a free external flow V_ext):
    fit y = M + beta*c_field + V_ext.n f(z) with the covariance plus sigma_v^2 f(z)^2 on the diagonal, profiling sigma_v,
    for the CF4++ and the 2023 ungrouped fields, whole field and non-dipole part.
(3) The 'correlation vs amplitude' decomposition the referee asks about: the correlation between the field's non-dipole part
    and the residuals, and the amplitude ratio that would be implied if the correlation were perfect; both reported, neither
    claimed as 'the' amplitude.
usage: python cf4_referee2.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov data/new_CF4pp_mean_std_grids.npz data/CF4_new_64-z008_velocity.fits
"""
import sys, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.ndimage import map_coordinates
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
C = 299792.458; OM = 0.334
def dc_h(z):
    e = lambda x: 1/np.sqrt(OM*(1+x)**3 + 1-OM)
    return C/100*np.array([quad(e, 0, zi)[0] for zi in np.atleast_1d(z)])
def mu(z): return 5*np.log10((1+z)*dc_h(z)/0.7) + 25
def zhd(zcmb, v): return (1+zcmb)/(1+v/C) - 1

def main(dat, covf, grid, fits2023):
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f:
        n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    m = (d.IS_CALIBRATOR == 0) & (d.zCMB > 0.01) & (d.zCMB < 0.1); idx = np.where(m)[0]
    zc = d.zCMB.values[idx]; r = dc_h(zc); ra, dec = d.RA.values[idx], d.DEC.values[idx]
    sc = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, distance=r*u.Mpc, frame='fk5', equinox='J2000').supergalactic
    xyz = np.stack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value]); rhat = xyz/np.linalg.norm(xyz, axis=0)
    g = np.load(grid); V = g['v_mean']; ci = (xyz + 500)/(1000/128) - 0.5
    vcf = np.sum(np.stack([map_coordinates(V[k], ci, order=1, mode='nearest') for k in range(3)])*rhat, 0)
    V23 = fits.getdata(fits2023).astype(float)*52.0; ci2 = (xyz + 500)/(1000/64) - 0.5
    v23 = np.sum(np.stack([map_coordinates(V23[k], ci2, order=1, mode='nearest') for k in range(3)])*rhat, 0)
    rad, de = np.radians(ra), np.radians(dec); nh = np.stack([np.cos(de)*np.cos(rad), np.cos(de)*np.sin(rad), np.sin(de)], -1)
    A = np.column_stack([r, nh])
    def split(v):
        coef, *_ = np.linalg.lstsq(A, v, rcond=None); return A@np.r_[0, coef[1:]], v - A@coef
    y = d.m_b_corr.values[idx] - mu(zc); Cm = cov[np.ix_(idx, idx)]
    f = (5/np.log(10))/(C*zc)   # mag per km/s at low z
    def gls(cols, Cadd=None):
        Cc = Cm if Cadd is None else Cm + Cadd
        Lc = np.linalg.cholesky(Cc); X = np.column_stack([np.ones(len(y))] + cols); Xw, yw = np.linalg.solve(Lc, X), np.linalg.solve(Lc, y)
        b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); e = np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))); rr = yw - Xw@b
        logdet = 2*np.sum(np.log(np.diag(Lc))); return b, e, rr@rr + logdet
    for name, v in [('CF4++', vcf), ('CF4 2023 ungrouped', v23)]:
        dip, rest = split(v); cr = mu(zhd(zc, rest)) - mu(zc); cf = mu(zhd(zc, v)) - mu(zc)
        b, e, _ = gls([cr]); print(f"\n{name}: beta_rest (GLS) = {b[1]:.3f} ± {e[1]:.3f}")
        # (1) sky-patch jackknife
        patches = (np.floor(ra/45).astype(int) % 8) + 8*(dec > 0) + 16*(np.abs(dec) > 30)
        ids = np.unique(patches); bj = []
        for p in ids:
            keep = patches != p
            if keep.sum() < 100: continue
            Ck = Cm[np.ix_(keep, keep)]; Lk = np.linalg.cholesky(Ck); X = np.column_stack([np.ones(keep.sum()), cr[keep]])
            Xw, yw = np.linalg.solve(Lk, X), np.linalg.solve(Lk, y[keep]); bb, *_ = np.linalg.lstsq(Xw, yw, rcond=None); bj.append(bb[1])
        bj = np.array(bj); k = len(bj); ej = np.sqrt((k-1)/k*np.sum((bj - bj.mean())**2))
        print(f"  delete-one-sky-patch jackknife ({k} patches of 8 RA x 2 Dec x 2 |Dec| bins): beta_rest = {bj.mean():.3f} ± {ej:.3f} (jackknife) vs ± {e[1]:.3f} (GLS); 1 - beta is {(1-bj.mean())/ej:.1f} sigma by the jackknife")
        # (2) Stiskalek-style: free sigma_v and V_ext, whole field
        best = None
        for sv in [0, 100, 150, 200, 250, 300, 350, 400, 500]:
            Cadd = np.diag((sv*f)**2)
            bw, ew, nll = gls([cf] + [f*nh[:, k] for k in range(3)], Cadd)
            br_, er_, nllr = gls([cr] + [f*nh[:, k] for k in range(3)], Cadd)
            if best is None or nll < best[0]: best = (nll, sv, bw[1], ew[1], br_[1], er_[1])
            if sv in (0, 250, 400): print(f"  sigma_v={sv:3d} + free V_ext: beta_whole = {bw[1]:.2f} ± {ew[1]:.2f} (|V_ext|={np.linalg.norm(bw[2:5]):.0f} km/s), beta_rest = {br_[1]:.2f} ± {er_[1]:.2f}")
        print(f"  best sigma_v = {best[1]} km/s: beta_whole = {best[2]:.2f} ± {best[3]:.2f}, beta_rest = {best[4]:.2f} ± {best[5]:.2f}")
        # (3) correlation vs amplitude
        Lc = np.linalg.cholesky(Cm); yw = np.linalg.solve(Lc, y - y.mean()); cw = np.linalg.solve(Lc, cr - cr.mean())
        rho = yw@cw/np.sqrt((yw@yw)*(cw@cw))
        print(f"  whitened correlation between residuals and the non-dipole prediction: {rho:+.3f}; beta_rest = rho * sd(y_w)/sd(c_w) = {rho*np.std(yw)/np.std(cw):.3f}. A correlation this size with a correctly scaled field would need the SN residuals to be {1/rho:.0f}x the field's prediction in rms; the data cannot separate 'too large' from 'mostly uncorrelated'.")

main(*sys.argv[1:5])

def joint_fit(dat, covf, grid):
    """referee pass 3: fit the dipole and non-dipole parts jointly (one GLS with both columns), CF4++ and 2M++."""
    import numpy as np, pandas as pd
    d = pd.read_csv(dat, sep=r'\s+')
    with open(covf) as f:
        n = int(f.readline()); cov = np.loadtxt(f).reshape(n, n)
    m = (d.IS_CALIBRATOR == 0) & (d.zCMB > 0.01) & (d.zCMB < 0.1); idx = np.where(m)[0]
    zc = d.zCMB.values[idx]; r = dc_h(zc); ra, dec = d.RA.values[idx], d.DEC.values[idx]
    sc = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, distance=r*u.Mpc, frame='fk5', equinox='J2000').supergalactic
    xyz = np.stack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value]); rhat = xyz/np.linalg.norm(xyz, axis=0)
    g = np.load(grid); V = g['v_mean']; ci = (xyz + 500)/(1000/128) - 0.5
    vcf = np.sum(np.stack([map_coordinates(V[k], ci, order=1, mode='nearest') for k in range(3)])*rhat, 0)
    rad, de = np.radians(ra), np.radians(dec); nh = np.stack([np.cos(de)*np.cos(rad), np.cos(de)*np.sin(rad), np.sin(de)], -1)
    A = np.column_stack([r, nh])
    def split(v):
        coef, *_ = np.linalg.lstsq(A, v, rcond=None); return A@np.r_[0, coef[1:]], v - A@coef
    y = d.m_b_corr.values[idx] - mu(zc); Lc = np.linalg.cholesky(cov[np.ix_(idx, idx)])
    print("\njoint fit of the dipole and non-dipole parts (one GLS, both columns):")
    for name, v in [('CF4++', vcf), ('2M++', d.VPEC.values[idx])]:
        dip, rest = split(v); cd = mu(zhd(zc, dip)) - mu(zc); cr = mu(zhd(zc, rest)) - mu(zc)
        X = np.column_stack([np.ones(len(y)), cd, cr]); Xw, yw = np.linalg.solve(Lc, X), np.linalg.solve(Lc, y)
        b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); cv = np.linalg.inv(Xw.T@Xw); e = np.sqrt(np.diag(cv)); rho = cv[1,2]/np.sqrt(cv[1,1]*cv[2,2])
        print(f"  {name}: beta_dipole = {b[1]:.2f} ± {e[1]:.2f}, beta_rest = {b[2]:.2f} ± {e[2]:.2f}, parameter correlation {rho:+.2f}")
joint_fit(*sys.argv[1:4])
