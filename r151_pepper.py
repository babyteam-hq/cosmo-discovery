"""R-151 step 0b: galpy streampepperdf mocks of the GC streams (subhalo-perturbed density along the stream).

Follows Bovy, Erkal & Sanders 2017 (arXiv:1606.03470) and Banik et al. 1911.02662: streamdf in MWPotential2014,
trailing arm, impacts at `ntimpact` equally spaced times, CDM subhalo number density per mass decade
n(M) = 0.3 (10^6.5 / M) kpc^-3 (their nsubhalo), Hernquist subhalos with r_s = 1.05 kpc (M/1e8)^0.5, impacts out to
5 r_s, subhalo velocity dispersion 120 km/s, masses 1e5-1e9 Msun drawn from dN/dM ~ M^-2 (their default).
The CDM-relative abundance scales the rate.

Each stream is a progenitor on an orbit with the published pericentre/apocentre (Li et al. 2022, Table 2; GD-1 from
Ibata et al. 2020), inclined as its orbital pole implies; sigv is set so the trailing arm's length roughly matches
the observed length in kpc. Output: density along the arm (path length in kpc) for nreal realizations per
abundance, saved as npz for the forecast (r151_forecast.py).

Usage: python -u r151_pepper.py <stream> <nreal> [ntimpact] [nTrackChunks]
"""
import sys
import time
import numpy as np
from scipy.optimize import brentq
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014
from galpy.actionAngle import actionAngleIsochroneApprox
from streamtools.df import streampepperdf  # jobovy/streamtools, galpy.util imports patched to galpy 1.12
from galpy.util import conversion

ro, vo = 8.0, 220.0
TGYR = conversion.time_in_Gyr(vo, ro)
MSOL = conversion.mass_in_msol(vo, ro)

# stream: r_peri, r_apo, inclination (deg), age (Gyr), observed arm length (kpc), sigv (km/s)
CFG = {
    "Phoenix": (12.7, 19.1, 61.0, 5.0, 7.4, 0.25),
    "ATLAS-AliqaUma": (12.1, 38.2, 66.0, 5.0, 15.6, 0.35),
    "GD-1": (14.0, 26.0, 120.0, 5.0, 12.0, 0.25),
}


def make_orbit(rp, ra, inc):
    i = np.radians(inc)

    def peri(v):
        o = Orbit([ra / ro, 0.0, v * np.cos(i) / vo, 0.0, v * np.sin(i) / vo, 0.0], ro=ro, vo=vo)
        return o.rperi(pot=MWPotential2014, analytic=True, type="staeckel") - rp

    v = brentq(peri, 30.0, 230.0, xtol=0.05)
    o = Orbit([ra / ro, 0.0, v * np.cos(i) / vo, 0.0, v * np.sin(i) / vo, 0.0], ro=ro, vo=vo)
    # move the progenitor a quarter period along so it is not sitting at apocentre
    ts = np.linspace(0.0, -0.08 / TGYR * 1.0, 2)
    o.integrate(ts, MWPotential2014)
    return o(ts[-1])


def nsubhalo(m):
    # Bovy et al. 2017's nsubhalo: subhalos per decade of mass per natural volume (8 kpc)^3, passed as is
    return 0.3 * (10.0**6.5 / m)


def rs(m):
    return 1.05 / ro * (m / 1e8) ** 0.5


def sample_mass(n, lo=1e5, hi=1e9):
    # dN/dM ~ M^-2  ->  uniform in 1/M
    u = np.random.uniform(size=n)
    return 1.0 / (1.0 / lo - u * (1.0 / lo - 1.0 / hi))


def main():
    name, nreal = sys.argv[1], int(sys.argv[2])
    ntimpact = int(sys.argv[3]) if len(sys.argv) > 3 else 16
    nchunks = int(sys.argv[4]) if len(sys.argv) > 4 else 26
    rp, ra, inc, age, Larm, sigv = CFG[name]
    t0 = time.time()
    prog = make_orbit(rp, ra, inc)
    aAI = actionAngleIsochroneApprox(pot=MWPotential2014, b=0.8)
    td = age / TGYR
    timpacts = [td * (k + 0.5) / ntimpact for k in range(ntimpact)]
    sdf = streampepperdf(sigv=sigv / vo, progenitor=prog, pot=MWPotential2014, aA=aAI, leading=False,
                         nTrackChunks=nchunks, nTrackIterations=1, sigMeanOffset=4.5, tdisrupt=td,
                         Vnorm=vo, Rnorm=ro, timpact=timpacts, spline_order=1, hernquist=True)
    print(f"{name}: setup {time.time() - t0:.0f}s", flush=True)
    # smooth stream: apar grid, path length per unit apar from the track
    alen = sdf.length(threshold=0.2)
    apars = np.linspace(0.0, alen, 101)
    dens0 = np.array([sdf._density_par(a) for a in apars])
    # path length (kpc) along the smooth track at each apar
    th = sdf._interpolatedThetasTrack
    XY = sdf._interpolatedObsTrackXY[:, :3] * ro
    X = np.array([np.interp(apars, th, XY[:, j]) for j in range(3)]).T
    s = np.concatenate([[0.0], np.cumsum(np.sqrt((np.diff(X, axis=0) ** 2).sum(1)))])
    print(f"{name}: arm length {s[-1]:.1f} kpc (observed {Larm} kpc), apar len {alen:.3f}", flush=True)
    rate1 = np.sum([sdf.subhalo_encounters(sigma=120.0 / vo, nsubhalo=nsubhalo(m), bmax=5.0 * rs(m))
                    for m in 10.0 ** np.arange(5.5, 9.0, 1.0)])
    print(f"{name}: CDM encounter rate (total impacts over the age) {rate1:.1f}", flush=True)
    out = dict(apars=apars, s=s, dens0=dens0, rate1=rate1)
    for ab in (0.2, 1.0):
        D = []
        for k in range(nreal):
            sdf.simulate(rate=rate1 * ab, sample_GM=lambda: sample_mass(1)[0] / MSOL,
                         sample_rs=lambda gm: rs(gm * MSOL), Xrs=5.0, sigma=120.0 / vo)
            D.append([sdf._density_par(a) for a in apars])
            if k % 10 == 0:
                print(f"{name}: abundance {ab} real {k} ({time.time() - t0:.0f}s)", flush=True)
        out[f"dens_{ab}"] = np.array(D)
        np.savez(f"research/cosmo-discovery/data/r151/pepper_{name}.npz", **out)
    print(f"{name}: done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
