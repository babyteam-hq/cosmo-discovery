"""How much does a SALT2 distance modulus move per unit error in the Milky-Way E(B-V) used in the fit? Checks the analytic
coefficient of dust_settle.py part (c), dmu/dE = R(lamB(1+z)) - beta [R(lamB(1+z)) - R(lamV(1+z))], with real SALT2 fits:
simulate a SALT2 (T23) light curve with F99 Milky-Way dust E(B-V)=E0 in the observer frame, in the band set a low-z survey
uses, then fit it with mwebv = E0 + dE held fixed (as the cosmology fits do) and read off d mB, d c, d x1 per unit dE,
so dmu/dE = dmB - beta dc + alpha dx1 (alpha 0.148, beta 3.04 for Pantheon+; ZTF's fitted beta also shown).
usage: python dust_salt_response.py
"""
import warnings, numpy as np, sncosmo
from astropy.table import Table
import extinction
warnings.filterwarnings("ignore")
AL, BE = 0.148, 3.04
BANDS = {'BVRI (CfA/LOSS/CSP-like)': ['bessellb', 'bessellv', 'bessellr', 'besselli'],
         'gri (Foundation/PS1/SDSS/DES)': ['sdssg', 'sdssr', 'sdssi'],
         'gri (ZTF)': ['ztfg', 'ztfr', 'ztfi']}
E0, DE = 0.05, 0.02

def analytic(z, beta):
    rB = extinction.fitzpatrick99(np.array([4400.0*(1+z)]), 3.1, 3.1)[0]; rV = extinction.fitzpatrick99(np.array([5500.0*(1+z)]), 3.1, 3.1)[0]
    return rB - beta*(rB - rV)

def response(z, bands, c0=0.0, x1_0=0.0):
    src = sncosmo.get_source('salt2')
    def model(ebv):
        m = sncosmo.Model(source=src, effects=[sncosmo.F99Dust()], effect_names=['mw'], effect_frames=['obs'])
        m.set(z=z, t0=0.0, x1=x1_0, c=c0, mwebv=ebv); m.set_source_peakabsmag(-19.3, 'bessellb', 'ab'); return m
    truth = model(E0)
    phases = np.arange(-10, 41, 3.0)*(1+z); rows = []
    for b in bands:
        try: sncosmo.get_bandpass(b)
        except Exception: continue
        f = truth.bandflux(b, phases, zp=25.0, zpsys='ab'); ok = np.isfinite(f) & (f > 0)
        for t, fl in zip(phases[ok], f[ok]): rows.append((t, b, fl, 0.01*np.max(f[ok]), 25.0, 'ab'))
    data = Table(rows=rows, names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'))
    out = {}
    for lbl, ebv in (('E0', E0), ('E0+dE', E0+DE)):
        m = model(ebv); m.set(mwebv=ebv)
        res, fit = sncosmo.fit_lc(data, m, ['t0', 'x0', 'x1', 'c'], bounds={'x1': (-4, 4), 'c': (-0.5, 0.5)})
        out[lbl] = (fit.source_peakmag('bessellb', 'ab'), fit['x1'], fit['c'])
    dmB = (out['E0+dE'][0]-out['E0'][0])/DE; dx1 = (out['E0+dE'][1]-out['E0'][1])/DE; dc = (out['E0+dE'][2]-out['E0'][2])/DE
    return dmB, dx1, dc

print("per unit E(B-V) error in the fixed Milky-Way E(B-V): d mB, d x1, d c from SALT2 refits, and dmu = dmB - beta dc + alpha dx1")
print(f"{'bands':32s} {'z':>5s} {'dmB':>6s} {'dx1':>6s} {'dc':>6s} | {'dmu(beta=3.04)':>14s} {'analytic':>8s} | {'dmu(beta=2.47)':>14s} {'analytic':>8s}")
for name, bands in BANDS.items():
    for z in (0.02, 0.04, 0.07, 0.1, 0.2, 0.4):
        if 'ZTF' in name and z > 0.1: continue
        if 'BVRI' in name and z > 0.1: continue
        try: dmB, dx1, dc = response(z, bands)
        except Exception as e: print(f"{name:32s} {z:5.2f} failed: {repr(e)[:80]}"); continue
        print(f"{name:32s} {z:5.2f} {dmB:6.2f} {dx1:6.2f} {dc:6.2f} | {dmB - BE*dc + AL*dx1:14.2f} {analytic(z, BE):8.2f} | {dmB - 2.47*dc + AL*dx1:14.2f} {analytic(z, 2.47):8.2f}")
