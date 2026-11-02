"""R-156: build the size sample (DESI DR1 BGS_BRIGHT clustering galaxies x LS DR9 Tractor shapes).

Writes data/r156/sample.npz with, per galaxy: RA, DEC, Z, WEIGHT, cap, Tractor columns, lambda = ln R_phys
(R_phys = SHAPE_R x D_A(z_obs), kpc/h), the measurement error of lambda, and the quality flag `good`.

Quality cut (fixed before any dipole is computed): Tractor TYPE in REX/EXP/DEV/SER, SHAPE_R > 0, SHAPE_R_IVAR > 0,
FRACMASKED_R < 0.1, MASKBITS free of bits 1 (BRIGHT), 12/13 (GALAXY/CLUSTER: SGA large galaxies), sigma_lambda < 0.5.
"""
import glob, os
import numpy as np
import fitsio

D = 'data/r156/'


def load_tractor():
    parts = [dict(np.load(f)) for f in sorted(glob.glob(D + 'tractor/*.npz'))]
    parts = [p for p in parts if len(p['TARGETID'])]
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


def main():
    import r156_theory as T
    t = load_tractor()
    g = t['ISBGS']
    tt = {k: v[g] for k, v in t.items()}
    _, first = np.unique(tt['TARGETID'], return_index=True)
    tt = {k: v[first] for k, v in tt.items()}
    cats = []
    for ci, c in enumerate(['NGC', 'SGC']):
        a = fitsio.read(D + 'BGS_BRIGHT_%s_clustering.dat.fits' % c,
                        columns=['TARGETID', 'Z', 'RA', 'DEC', 'WEIGHT', 'flux_g_dered', 'flux_r_dered', 'flux_z_dered'])
        cats.append((a, ci))
    out = {}
    for a, ci in cats:
        idx = np.searchsorted(tt['TARGETID'], a['TARGETID'])
        idx[idx >= len(tt['TARGETID'])] = 0
        ok = tt['TARGETID'][idx] == a['TARGETID']
        rec = {'Z': a['Z'], 'RA': a['RA'], 'DEC': a['DEC'], 'WEIGHT': a['WEIGHT'], 'CAP': np.full(len(a), ci, 'i1'),
               'HAS_TRACTOR': ok,
               'GMAG': 22.5 - 2.5 * np.log10(np.clip(a['flux_g_dered'], 1e-3, None)),
               'RMAG': 22.5 - 2.5 * np.log10(np.clip(a['flux_r_dered'], 1e-3, None)),
               'ZMAG': 22.5 - 2.5 * np.log10(np.clip(a['flux_z_dered'], 1e-3, None))}
        for k in ['TYPE', 'SHAPE_R', 'SHAPE_R_IVAR', 'SERSIC', 'PSFSIZE_R', 'PSFSIZE_G', 'PSFSIZE_Z', 'FRACFLUX_R',
                  'FRACMASKED_R', 'FRACIN_R', 'MASKBITS', 'NOBS_R', 'PSFDEPTH_R', 'GALDEPTH_R', 'SHAPE_E1', 'SHAPE_E2',
                  'FIBERFLUX_R', 'FLUX_R', 'MW_TRANSMISSION_R', 'RELEASE', 'NEA_R', 'BLOB_NEA_R', 'RCHISQ_R']:
            v = tt[k][idx]
            if v.dtype.kind != 'S':
                v = np.where(ok, v, 0)
            rec[k] = v
        for k, v in rec.items():
            out.setdefault(k, []).append(v)
    out = {k: np.concatenate(v) for k, v in out.items()}
    z = out['Z']
    zg = np.linspace(0.0, 0.7, 701)
    rg = np.array([T.background(x)[0] for x in zg])
    r = np.interp(z, zg, rg)
    DA = r / (1 + z) * 1000.0  # kpc/h
    R = out['SHAPE_R'] * np.pi / 180 / 3600 * DA
    sr = out['SHAPE_R']
    typ = np.char.strip(out['TYPE'].astype('U3'))
    ext = np.isin(typ, ['REX', 'EXP', 'DEV', 'SER'])
    with np.errstate(divide='ignore', invalid='ignore'):
        lam = np.log(R)
        slam = 1.0 / (sr * np.sqrt(out['SHAPE_R_IVAR']))
    mb = out['MASKBITS'].astype(int)
    good_nosig = (out['HAS_TRACTOR'] & ext & (sr > 0) & (out['SHAPE_R_IVAR'] > 0) & (out['FRACMASKED_R'] < 0.1)
                  & ((mb & ((1 << 1) | (1 << 12) | (1 << 13))) == 0) & np.isfinite(slam))
    good = good_nosig & (slam < 0.5)
    out.update({'LAMBDA': np.where(good, lam, 0.0), 'SIG_LAMBDA': np.where(good, slam, 9.0), 'GOOD': good, 'GOOD_NOSIG': good_nosig,
                'LAMBDA_RAW': np.where(good_nosig, lam, 0.0), 'SIG_LAMBDA_RAW': np.where(good_nosig, slam, 99.0),
                'TYPE': typ, 'RCOM': r})
    np.savez(D + 'sample.npz', **out)
    print('BGS_BRIGHT clustering galaxies %d; with tractor %d; good sizes %d' % (len(z), out['HAS_TRACTOR'].sum(),
                                                                                  good.sum()))
    for zlo in (0.15, 0.25, 0.35):
        m = (z > zlo) & (z < zlo + 0.1)
        print('z %.2f-%.2f: all %d, with tractor %d, good %d; types %s' % (
            zlo, zlo + 0.1, m.sum(), (m & out['HAS_TRACTOR']).sum(), (m & good).sum(),
            {k: int((m & (typ == k)).sum()) for k in ['PSF', 'REX', 'EXP', 'DEV', 'SER']}))


if __name__ == '__main__':
    main()
