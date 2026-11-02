"""R-150 step 0: slim the DESI DR1 catalogues used by the halo-dust reddening stack.

usage: python -I dust_z1_prep.py
reads  data/dust/lss/{ELG_LOPnotqso,LRG}_{NGC,SGC}_clustering.{dat,0 ran}.fits, data/dust/qso/QSO_full.dat.fits
writes data/dust/fg_{ELG,LRG}.npz, data/dust/ran_{ELG,LRG}.npz, data/dust/qso.npz (all gitignored)
"""
import os
import numpy as np
import fitsio

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, 'data', 'dust')


def slim_lss(tracer):
    for kind, tag in (('dat', 'fg'), ('ran', 'ran')):
        parts = []
        for cap in ('NGC', 'SGC'):
            fn = f'{tracer}_{cap}_clustering.dat.fits' if kind == 'dat' else f'{tracer}_{cap}_0_clustering.ran.fits'
            a = fitsio.read(os.path.join(D, 'lss', fn), columns=['RA', 'DEC', 'Z', 'WEIGHT', 'PHOTSYS'])
            parts.append(a)
        a = np.concatenate(parts)
        name = tracer.split('_')[0]
        np.savez(os.path.join(D, f'{tag}_{name}.npz'), ra=a['RA'], dec=a['DEC'], z=a['Z'], w=a['WEIGHT'],
                 north=(a['PHOTSYS'] == 'N'))
        print(f'{tag}_{name}: {len(a)} rows, z {np.percentile(a["Z"], [1, 50, 99]).round(3)}')


def slim_qso():
    cols = ['TARGETID', 'RA', 'DEC', 'Z_not4clus', 'PHOTSYS', 'DESI_TARGET', 'EBV', 'MASKBITS',
            'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1',
            'FLUX_IVAR_W2', 'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1',
            'MW_TRANSMISSION_W2', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'COEFF', 'NCOEFF', 'Z_RR', 'SPECTYPE']
    a = fitsio.read(os.path.join(D, 'qso', 'QSO_full.dat.fits'), columns=cols)
    z = a['Z_not4clus']
    good = np.isfinite(z) & (z > 0) & (z < 7)
    a = a[good]
    out = {k.lower(): a[k] for k in cols if k not in ('PHOTSYS', 'SPECTYPE', 'Z_not4clus')}
    out['z'] = a['Z_not4clus']
    out['north'] = a['PHOTSYS'] == 'N'
    out['spectype_qso'] = a['SPECTYPE'] == 'QSO'
    np.savez(os.path.join(D, 'qso.npz'), **out)
    print(f'qso: {good.sum()} with a good redshift of {len(good)} targets; z>2: {(out["z"] > 2).sum()}')


if __name__ == '__main__':
    slim_lss('ELG_LOPnotqso')
    slim_lss('LRG')
    slim_qso()
