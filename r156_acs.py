"""R-156 method-check point 1: the Tractor size response, calibrated against HST/ACS sizes in COSMOS.

The response a = d lambda_DR9 / d lambda_true is the OLS slope of the DR9 size residual (the analysis's own kappa_hat,
type-specific relation fitted on the full sample) on the ACS size residual (ln of the ZEST/GIM2D half-light radius
from HST/ACS F814W, Scarlata et al. 2007, IRSA `cosmos_morph_zurich_1`, against the same z and colour design fitted
on the matched galaxies). The ACS PSF (0.1") is far below the sizes here (~1.4"), so ACS sizes stand in for true ones.
Error: bootstrap over galaxies. Check: the ZEST r50 (SExtractor-free 50 % light radius) in place of GIM2D.
Because R ~ L^0.5 along the population at fixed z and colour, the regression direction is the constant-surface-
brightness direction a magnification moves along.
"""
import numpy as np
import r156_kappa as K
import r156_theory as T


def main():
    s = dict(np.load('data/r156/sample.npz'))
    z = s['Z']
    sel = s['GOOD'] & (z > 0.15) & (z < 0.45)
    rel = K.fit_relation(s, sel, False)
    k, u = K.kappa_hat(s, rel, False)
    acs = np.genfromtxt('data/r156/cosmos/zest.csv', delimiter=',', names=True)
    acs = acs[(acs['acs_clean'] == 1) & (acs['acs_mu_class'] == 1) & (acs['r_gim2d'] > 0) & (acs['r50'] > 0)]
    box = s['GOOD'] & (np.abs(s['RA'] - 150.1) < 1.2) & (np.abs(s['DEC'] - 2.2) < 1.2) & (z > 0.05) & (z < 0.5)
    idx = np.where(box)[0]
    from scipy.spatial import cKDTree
    def xyz(ra, dec):
        ra, dec = np.radians(ra), np.radians(dec)
        return np.vstack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]).T
    tr = cKDTree(xyz(acs['ra'], acs['dec']))
    dist, j = tr.query(xyz(s['RA'][idx], s['DEC'][idx]))
    ok = dist < np.radians(1.0 / 3600)
    idx, j = idx[ok], j[ok]
    zz = z[idx]
    zg = np.linspace(0, 0.7, 701)
    DA = np.interp(zz, zg, [T.background(x)[0] for x in zg]) / (1 + zz) * 1000
    out = {}
    for col, scale in (('r_gim2d', 1.0), ('r50', 0.03)):
        lam_acs = np.log(acs[col][j] * scale * np.pi / 180 / 3600 * DA)
        gr = s['GMAG'][idx] - s['RMAG'][idx]
        rz = s['RMAG'][idx] - s['ZMAG'][idx]
        X = K.design(zz, gr, rz, s['RMAG'][idx], False)
        c, *_ = np.linalg.lstsq(X, lam_acs, rcond=None)
        ka = lam_acs - X @ c
        kd = k[idx]
        good = np.abs(ka) < 4 * np.std(ka)
        def slope(ii):
            return np.polyfit(ka[good][ii], kd[good][ii], 1)[0]
        n = good.sum()
        a = slope(np.arange(n))
        rng = np.random.default_rng(3)
        bs = [slope(rng.integers(0, n, n)) for _ in range(500)]
        for zlo, zhi in ((0.05, 0.5), (0.15, 0.45), (0.15, 0.25), (0.25, 0.45)):
            mm = (zz[good] > zlo) & (zz[good] < zhi)
            if mm.sum() > 30:
                aa = np.polyfit(ka[good][mm], kd[good][mm], 1)[0]
                bb = np.std([np.polyfit(ka[good][mm][q], kd[good][mm][q], 1)[0]
                             for q in [rng.integers(0, mm.sum(), mm.sum()) for _ in range(300)]])
                out['%s z%.2f-%.2f' % (col, zlo, zhi)] = (round(float(aa), 3), round(float(bb), 3), int(mm.sum()))
        out['%s all' % col] = (round(float(a), 3), round(float(np.std(bs)), 3), int(n))
        out['%s corr' % col] = round(float(np.corrcoef(ka[good], kd[good])[0, 1]), 3)
    for kk, v in out.items():
        print(kk, v)
    return out


if __name__ == '__main__':
    main()
