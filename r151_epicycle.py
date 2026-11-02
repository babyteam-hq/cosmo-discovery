"""R-151 step 0a: predicted epicyclic density-oscillation scale per stream.

Kuepper et al. 2008/2012 (as used by Ibata et al. 2002.01488, eqs 12-13):
    y_C = (4 pi Omega / kappa) |1 - 4 Omega^2 / kappa^2| x_E,
    x_E = 2.88 (G M / (4 Omega^2 - kappa^2))^(1/3),
with Omega, kappa at the radius where stars are stripped. The orbits are eccentric, so we evaluate at pericentre
(where stripping is concentrated) and at the guiding radius (r_peri + r_apo)/2, and take the pair as the range.
Orbits: Li et al. 2022 (S5, arXiv:2110.06950, Table 2) for the DES streams; GD-1 from Ibata et al. 2020.
Progenitor masses: GC streams 1e4-1e5 Msun (bound mass while stripping; stream stellar mass from Patrick et al.
2206.04241 as the floor); dwarf streams total (DM) mass from Yang et al. 2603.09305 (Indus 3e8, 1e8-1e9 range).
Potential: galpy MWPotential2014 (Banik et al.'s choice), rescaled to v_c(R0)=220 km/s by galpy's natural units.
"""
import json
import numpy as np
from galpy.potential import MWPotential2014, omegac, epifreq
from galpy.util import conversion

ro, vo = 8.0, 220.0
G = 4.30091e-6  # kpc (km/s)^2 / Msun

# name: (type, r_peri, r_apo, helio distance kpc, observed length deg, observed length kpc (Patrick), M range Msun)
STREAMS = {
    "ATLAS-AliqaUma": ("GC", 12.1, 38.2, 23.0, 37.0, 15.6, (1e4, 1e5)),
    "Phoenix": ("GC", 12.7, 19.1, 19.1, 13.6, 7.4, (1e4, 1e5)),
    "Jhelum": ("DG", 9.1, 27.8, 13.2, 29.2, 5.3, (1e7, 1e9)),
    "Indus": ("DG", 12.3, 17.9, 16.6, 20.3, 5.6, (1e8, 1e9)),
    "Turranburra": ("DG", 8.3, 35.2, 27.5, 16.9, 8.2, (1e7, 1e9)),
    "GD-1": ("GC", 14.0, 26.0, 8.5, 80.0, 12.0, (1e4, 1e5)),
}


def freqs(r_kpc):
    R = r_kpc / ro
    om = omegac(MWPotential2014, R) * conversion.freq_in_kmskpc(vo, ro)
    ka = epifreq(MWPotential2014, R) * conversion.freq_in_kmskpc(vo, ro)
    return om, ka  # km/s/kpc


def y_c(r_kpc, M):
    om, ka = freqs(r_kpc)
    xe = 2.88 * (G * M / (4 * om**2 - ka**2)) ** (1 / 3)
    return 4 * np.pi * om / ka * abs(1 - 4 * om**2 / ka**2) * xe


def main():
    out = {}
    print(f"{'stream':16s} type  r(kpc)   Om/ka    yC(Mlo)  yC(Mhi) kpc   deg(lo-hi)   L/yC cycles in window")
    for name, (typ, rp, ra, d, Ldeg, Lkpc, (mlo, mhi)) in STREAMS.items():
        rg = 0.5 * (rp + ra)
        rows = []
        for r in (rp, rg):
            om, ka = freqs(r)
            rows.append((r, om / ka, y_c(r, mlo), y_c(r, mhi)))
        ylo = min(x[2] for x in rows)
        yhi = max(x[3] for x in rows)
        kpc_per_deg = Lkpc / Ldeg
        dlo, dhi = ylo / kpc_per_deg, yhi / kpc_per_deg
        # the scale at the geometric centre of the mass range and the two radii
        ymid = np.sqrt(ylo * yhi)
        out[name] = dict(type=typ, r_peri=rp, r_apo=ra, d=d, L_deg=Ldeg, L_kpc=Lkpc, M=(mlo, mhi),
                         yC_kpc=(ylo, yhi), yC_deg=(dlo, dhi), yC_mid_deg=ymid / kpc_per_deg,
                         cycles=(Lkpc / yhi, Lkpc / ylo))
        for r, q, a, b in rows:
            print(f"{name:16s} {typ}   {r:5.1f}   {q:5.3f}   {a:7.2f}  {b:7.2f}")
        print(f"{'':16s}  -> yC {ylo:.2f}-{yhi:.2f} kpc = {dlo:.1f}-{dhi:.1f} deg along a {Ldeg:.0f} deg window;"
              f" {Lkpc / yhi:.1f}-{Lkpc / ylo:.1f} cycles")
    json.dump(out, open("research/cosmo-discovery/r151_epicycle.json", "w"), indent=1)


if __name__ == "__main__":
    main()
