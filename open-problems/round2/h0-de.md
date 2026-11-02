# Round 2: expansion history, H0 and dark energy, quick-kill survivors

Checked 2026-11-01. The arXiv API (export.arxiv.org) returned HTTP 429 on every
call this session, so the arXiv checks used arXiv's own search
(arxiv.org/search, all fields, newest first), which covers the same index.
Query counts below include those searches and WebSearch (extended).

**Killed in the quick kill:**
- *DESI spectroscopic BAO inside the DES Y6 footprint* (to decide whether DES Y6's
  low D_M/r_d at z=0.85 is sample variance or photo-z systematics): the
  DES-independent split exists (2601.14864: overlap region alpha = 0.91 +/- 0.063,
  non-overlap 0.969 +/- 0.030, called a statistical fluctuation), and DESI DR1
  angular BAO exists (2610.10215). A DESI-in-overlap number could not separate
  the two answers: the DES overlap error (6%) dominates.
- *Cosmic chronometers from DESI DR1 LRG spectra*: done (2608.13178).
- *Matter-radiation equality / turnover with DESI*: done (2505.16153, 2511.06962).
- *Tomographic AP with DESI*: 2510.24089. *HII-galaxy Hubble diagram to z~14 with
  JWST*: 2607.14254. *ISW and w0wa*: 2607.22461 (Quaia), and ISW cosmic variance
  can't separate w0wa from LCDM anyway.
- *Quasar (X-ray/UV) Hubble diagram with eRASS1*: 2508.00732 says current eROSITA
  samples give no significant cosmology, so it can't decide Risaliti-Lusso.
- *ZTF DR2 as an independent low-z SN anchor* (Efstathiou's DES5YR low-z offset,
  2502.04212): ZTF's own paper (2509.04073) says DR2 photometry has ~30 mmag
  calibration uncertainty and is not fit for cosmology, so it can't resolve a
  0.04 mag offset; DR2.5 is the collaboration's.
- *JWST SBF distances in the Hubble flow*: unpublished, but it is the TRGB-SBF
  team's own program (GO-3055 and follow-ups, 2603.11160 announces it).
- *Quasar time dilation to z>4.3 with ZTF + DESI* (extends 2306.04053,
  2501.04171): unpublished, but ZTF's ~7 yr baseline is shorter than the
  observed DRW timescales at z>4, which biases the timescales.
- *Cosmic chronometers at z>2 from JWST quiescent galaxies*: unpublished (3
  queries), dropped for precision: progenitor bias and SPS-model ages give H(z)
  no better than ~30%, with no mock-based error budget.

---

## 1. Is there a KBC-size local void? Spectroscopic density profile from DESI DR1 BGS

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

## 2. Omega_m and D_M/D_H from void-galaxy RSD+AP in DESI DR1

- **Question:** What do the void-galaxy cross-correlation multipoles of DESI
  DR1 (BGS, LRG, ELG; 0.1<z<1.6) give for D_M/D_H(z) and fsigma8, and so for
  Omega_m in flat LCDM, independent of r_d and of SN calibration?
- **Shape:** extension (sample: DESI DR1 has never had this measurement;
  the largest so far are SDSS/BOSS/eBOSS) that also bears on a disagreement.
  Prior: Nadathur et al. 2019 (1904.01030, BOSS), Hamaus et al. 2020
  (2007.07895, BOSS, Omega_m = 0.312 +/- 0.020), Nadathur et al. 2020
  (2008.06060, eBOSS), Woodfinden et al. 2022 (2205.06258, SDSS 0.07<z<1.0,
  Omega_m = 0.337 +0.026/-0.029). The disagreement it bears on: Omega_m from
  DESI DR2 BAO (~0.297 +/- 0.009, 2503.14738) vs uncalibrated SNe
  (Pantheon+ 0.334 +/- 0.018, 2202.04077; DES-SN5YR 0.352 +/- 0.017,
  2401.02929), the tension behind the w0wa signal. Void AP measures the same
  D_M/D_H shape as BAO-AP but from a different scale and estimator.
- **Public data:**
  - DESI DR1 LSS clustering catalogues and randoms:
    https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5/
  - DESI DR1 mocks (1000 EZmocks for covariance, 25 AbacusSummit for
    calibrating the model): https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/mocks
  - Void finder with RSD reconstruction: Revolver
    (https://github.com/seshnadathur/Revolver).
- **Method:**
  1. Reconstruct real-space positions per tracer and redshift bin, then find
     voids on the reconstructed field (removes the void-finding RSD bias).
  2. Measure the void-galaxy monopole and quadrupole in redshift space.
  3. Fit the linear-dispersion model for beta and epsilon = D_M/D_H ratio
     relative to fiducial, calibrated on AbacusSummit, with EZmock covariance.
  4. Convert to Omega_m (flat LCDM) and compare with DESI BAO and SN values;
     blind epsilon until the mock tests pass.
- **Expected precision:** DR1's effective volume is ~2.5-3x SDSS's, so
  sigma(Omega_m) ~ 0.015-0.018 from voids alone (from Woodfinden's 0.027).
  That separates 0.297 from 0.352 at ~3 sigma and from 0.334 at ~2 sigma:
  enough to say which side the void AP sits on, not to settle it alone.
- **Quick-kill queries (7):** WebSearch extended "DESI DR1 voids
  redshift-space distortions Alcock-Paczynski void-galaxy cross-correlation
  Omega_m measurement" (found only 2610.18348, DR1 void catalogues: size
  function, shapes, profiles, bias; checked, no AP or fsigma8); arXiv
  "void-galaxy cross-correlation DESI" (only Legacy Survey void lensing, ISW,
  forecasts); "void Alcock-Paczynski DESI DR1" (0); "voids DESI growth rate"
  (0); "void Alcock-Paczynski test 2025" (0); "void redshift-space distortions
  DESI" (2610.18348 again, emulators, finders); "Woodfinden voids" (SDSS
  analysis and Euclid forecasts only). Note: void *lensing* and the void *size
  function* are on our skip list; this is neither.
