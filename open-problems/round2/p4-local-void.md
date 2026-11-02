# Proposed test: Is there a KBC-size local void? Spectroscopic density profile from DESI DR1 BGS

- **Question:** What is the galaxy density contrast delta(<R) of a
  volume-limited sample out to R = 300 Mpc (z~0.07), relative to the same
  sample's density at 0.1<z<0.2, in the DESI footprint? A void of -20% to -33%
  at that scale would bias local H0 high by a few km/s/Mpc.
- **Shape:** disagreement. Galaxy counts say yes: Keenan, Barger & Cowie 2013
  (arXiv:1304.2884, K-band luminosity density ~1.5x higher beyond z~0.07, i.e.
  a ~-33% void to ~300 Mpc); Whitbourn & Shanks 2014/2016 (1307.4405,
  1603.02322) and Wong et al. 2022 (2107.08505: under-density over 90% of the
  sky to ~200 Mpc). Distances and clusters say no void that large:
  Kenworthy, Scolnic & Riess 2019 (1901.08681: SNe exclude |delta|>20% at
  4-5 sigma, delta<27% at 5 sigma beyond 69 Mpc/h); Böhringer et al. 2020
  (1907.12402, CLASSIX clusters: an under-density, but on ~100 Mpc scales).
  The count claims hinge on luminosity-function normalisation and evolution
  in photometric or sparse samples; no-one has measured the profile with a
  large, uniformly selected spectroscopic sample spanning both the local
  volume and z~0.2.
- **Public data:**
  - DESI DR1 BGS LSS catalogues (full catalogues with completeness weights and
    randoms; the "clustering" files start at z=0.1, so use the full ones):
    https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5/
  - DESI DR1 FastSpecFit VAC for k-corrections and stellar masses:
    https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/
  - DESI DR1 mocks (25 AbacusSummit, 1000 EZmocks):
    https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/mocks ; for
    observer-centred cones down to z~0.01, AbacusSummit halo light cones
    (https://abacussummit.readthedocs.io) with a BGS HOD.
  - 2MRS / 2M++ for an all-sky comparison in the same K-band as KBC.
- **Method:**
  1. Build volume-limited BGS samples (several M_r and stellar-mass cuts),
     weighted for fibre completeness, NGC and SGC separately.
  2. Measure n(z) in shells from 30 to 600 Mpc and delta(<R) relative to
     0.1<z<0.2, fitting luminosity evolution Q jointly (it shifts the bright
     and faint bins differently; a real void shifts them equally).
  3. Get the LCDM distribution of delta(<R) for this footprint and selection
     from many observer positions in the mocks.
  4. Report delta(<300 Mpc) with its interval, and the posterior on the KBC
     profile vs the Kenworthy bound.
- **Expected precision:** ~2x10^4 galaxies with M_r < -20.5 inside z<0.07
  in the DR1 footprint (shot noise <1%); cosmic variance for the cone ~5%,
  evolution systematic ~3% (Q +/- 0.3 mag per unit z over dz~0.1): total
  sigma(delta) ~ 6-7% per cap. A -20% to -33% void stands out at 3-5 sigma;
  delta ~ 0 rules out the KBC amplitude in this footprint. It can tell the
  answers apart, with the caveat that DESI sees ~1/5 of the sky.
- **Quick-kill queries (9):** WebSearch extended "local void underdensity
  300 Mpc eROSITA eRASS1 clusters OR DESI Legacy number counts Hubble tension"
  (only theory, BAO-based and CLASSIX papers); WebSearch extended "DESI Bright
  Galaxy Survey local underdensity KBC void luminosity function number density
  z<0.1" (DESI DR2 LF paper 2511.01803 covers 0.002<z<0.6 but measures no
  local-hole profile, checked; DESIVAST voids only); arXiv "local void DESI"
  (2606.02712 and 2604.05466 are BAO/SN reconstructions, not counts); "local
  hole underdensity"; "KBC void" (2506.10518 tests it with CF4 distances, not
  counts); "galaxy number counts local underdensity"; "local underdensity
  spectroscopic redshift survey Hubble tension"; "Local Hole DESI"; "DESI
  luminosity function density evolution local" (0). No DESI BGS measurement of
  the local density profile was found.
