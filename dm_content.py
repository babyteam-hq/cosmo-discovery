"""Dark-matter content of NGC 1052-DF2, DF4, DF9 as a function of distance and of which velocity dispersion one takes.
Wolf+2010: M(<r_1/2) = 930 sigma^2 r_1/2 [Msun, km/s, pc], r_1/2 = 4/3 R_e (3D half-light radius); stellar mass inside
r_1/2 = half the total. Angular sizes and fluxes are what's measured, so M_dyn ∝ D and M_* ∝ D^2: M_dyn/M_* ∝ 1/D.
Published inputs (at the paper's reference distance D0):
  DF2: R_e,c = 2.0 kpc at 20 Mpc (Cohen+2018 via Danieli+2019), M_*(<r_1/2) = 1.0±0.2e8 at 20 Mpc; sigmas: Danieli+2019 stars
       8.5 +2.3 -3.1; vD+2018 GCs 7.8 +5.2 -2.2; Keim+2026 KCWI 6.3 +3.1 -3.7; MUSE (Emsellem+2019 as quoted) 9.2 +3.8 -4.5.
       Normal-halo expectation (Danieli+2019, SMHM+NFW): 35±6 km/s.
  DF4: R_e = 1.6 kpc at 20 Mpc, M_* = 1.5±0.4e8 (vD+2019); sigmas: GCs 4.2 +4.4 -2.2 (vD+2019); KCWI 7.8 +2.3 -2.0 (Keim+2026);
       Shen+2023 continuum spectroscopy also low. Normal-halo expectation ~30 km/s.
  DF9: R_e = 11.1" = 1.1 kpc at 20.6 Mpc, M_*(<R_e) = 0.71±0.16e8 at 20.6 (Keim+2026); sigma 6.5 +3.9 -4.3; normal halo 24±3 km/s.
Distances to scan: 13 (Trujillo+2019), 17.6 (JWST TRGB, DF2), 18.7-18.8 (SBF), 20, 21.7 (HST TRGB, DF2).
usage: python dm_content.py
"""
import numpy as np
rng = np.random.default_rng(3)
def split_normal(mu, up, dn, n):
    s = rng.standard_normal(n); return np.where(s > 0, mu + s*up, mu + s*dn)
def ratio(sig, rhalf_kpc_D0, mstar_half_D0, mstar_err, D0, D, n=200000):
    s = np.clip(split_normal(*sig, n), 0.1, None)
    m = np.clip(rng.normal(mstar_half_D0, mstar_err, n), 1e6, None)
    mdyn = 930*s**2*(rhalf_kpc_D0*1000*D/D0)
    ms = m*(D/D0)**2
    r = mdyn/ms
    return np.percentile(r, [16, 50, 84]), np.percentile(mdyn, [16, 50, 84])/1e8, np.median(ms)/1e8
gal = {
 'DF2': dict(D0=20.0, rhalf=4/3*2.0, mhalf=1.0e8, merr=0.2e8, sigmas={'stars Danieli+19': (8.5, 2.3, 3.1), 'GCs vD+18': (7.8, 5.2, 2.2), 'KCWI Keim+26': (6.3, 3.1, 3.7), 'MUSE': (9.2, 3.8, 4.5)}, halo=35),
 'DF4': dict(D0=20.0, rhalf=4/3*1.6, mhalf=0.75e8, merr=0.2e8, sigmas={'GCs vD+19': (4.2, 4.4, 2.2), 'KCWI Keim+26': (7.8, 2.3, 2.0)}, halo=30),
 'DF9': dict(D0=20.6, rhalf=4/3*1.1, mhalf=0.71e8, merr=0.16e8, sigmas={'KCWI Keim+26': (6.5, 3.9, 4.3)}, halo=24),
}
Ds = [13.0, 17.6, 18.7, 20.0, 21.7]
for name, g in gal.items():
    print(f"\n{name}: M_dyn(<r_1/2)/M_*(<r_1/2); 1 = stars only, {(g['halo']/list(g['sigmas'].values())[0][0])**2:.0f}x = a normal halo at the first sigma")
    print("  sigma source          " + "  ".join(f"D={D:5.1f}" .rjust(22) for D in Ds))
    for sname, sig in g['sigmas'].items():
        row = []
        for D in Ds:
            (lo, med, hi), md, ms = ratio(sig, g['rhalf'], g['mhalf'], g['merr'], g['D0'], D)
            row.append(f"{med:4.1f} (+{hi-med:.1f} -{med-lo:.1f})".rjust(22))
        print(f"  {sname:20s}" + "  ".join(row))
    # how far the normal-halo expectation is: ratio a halo would give, scaled the same way, vs observed
    for D in (17.6, 20.0):
        (lo, med, hi), md, ms = ratio(list(g['sigmas'].values())[0], g['rhalf'], g['mhalf'], g['merr'], g['D0'], D)
        halo_ratio = 930*g['halo']**2*(g['rhalf']*1000*D/g['D0'])/(g['mhalf']*(D/g['D0'])**2)
        print(f"  at D={D}: M_dyn = {md[1]:.2f} (+{md[2]-md[1]:.2f} -{md[1]-md[0]:.2f}) e8, M_* = {ms:.2f} e8; a normal halo (sigma={g['halo']}) would give M_dyn/M_* = {halo_ratio:.0f}")
