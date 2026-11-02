# Round 6: radio, HI and multi-messenger harvest (2026-11-02)

Scope: radio continuum and HI (RACS/SPICE-RACS, EMU, VLASS, MeerKAT, LoTSS DR2/DR3, WALLABY, FAST) and multi-messenger
archives (IceCube, Auger, PTAs), alone or with DESI DR1, ACT DR6, Planck, eRASS1, Gaia DR3. Search and reasoning only:
no code was run and no data were downloaded. Every number below is a desk estimate for the feasibility gate to check.

**Survivors: 3** (ranked). The strongest is R6-1. R6-3 is a marginal pass: it clears 3 sigma only on optimistic counts.

New facts that change earlier kills:

- **LoTSS DR3 is public** (Shimwell et al. 2026, A&A 707, A198, arXiv:2602.15949; https://lofar-surveys.org/dr3.html):
  88 % of the northern sky, 13.7 M sources. The zero-overlap kill of LoTSS x eRASS1 in round 3 no longer holds: they now
  share roughly 3,000-4,000 deg^2 (dec > 0, l > 180). No cosmology survivor came out of this here (see "Killed").
- **SPICE-RACS DR2 is public** (arXiv:2605.16917): 2.5e5 RMs above 8 sigma (3.4e5 above 6 sigma), 799.5-1087.5 MHz,
  median sigma_RM about 2 rad m^-2, from the south pole to dec +49 (87.5 % of the sky), via CSIRO DAP / CASDA. It is the
  largest broadband RM catalogue, and R6-1 and R6-2 are built on it.

---

## R6-1. Faraday-rotation excess of random over physical pairs at 0.9 GHz (SPICE-RACS DR2)

**Question.** Is the extragalactic RM difference between sources that are close on the sky but at different redshifts
(random pairs) and between the two lobes of one source (physical pairs) about 10 rad m^-2, as at 1.4 GHz, or about
1-2 rad m^-2, as at 144 MHz?

**Shape: disagreement.** The two published answers differ by a factor of 6-10:

- **NVSS, 1.4 GHz** (Vernstrom et al. 2019, arXiv:1905.02410): Taylor et al. 2009 two-channel RMs, |b| >= 20 deg,
  polarisation fraction >= 2 %, separations 1'-20', 317 physical and 5,111 random pairs. rms DeltaRM is 14.9 +/- 0.4
  (random) against 4.6 +/- 1.1 (physical), about 10 rad m^-2 apart, "significant at 5 sigma". This gives an IGM field up
  to 40 nG.
- **LoTSS, 144 MHz:**
  - O'Sullivan et al. 2020 (arXiv:2002.06924): no significant difference, excess < 1.9 rad m^-2 (95 %). Trimmed rms is
    1.8 (random) against 1.4 (physical).
  - Pomakov et al. 2022 (arXiv:2208.01336, DR2): median excess 1.65 +/- 0.10 rad m^-2, with no redshift dependence;
    B0 <~ 2 nG.

The review arXiv:2505.18619 says the RM and stacking constraints can't both be fit. Two steps are disputed:

- **NVSS's RM determination.** Two narrow channels leave n-pi ambiguities, and off-axis leakage biases the RMs
  (arXiv:1810.03620).
- **LOFAR's frequency selection.** Depolarisation at lambda^2 = 4.3 m^2 removes sightlines with large Faraday dispersion.

A broadband RM-synthesis catalogue at 0.9 GHz has neither problem: it has no n-pi ambiguity, and depolarisation is about
40x weaker than at 144 MHz. So it can tell which answer is the physical one.

**Data.**
- SPICE-RACS DR2 RM catalogue (arXiv:2605.16917), from CASDA https://research.csiro.au/casda/ and CSIRO DAP
  https://data.csiro.au/. Use the RACS-low/mid component catalogues for morphology and host IDs.
- Host redshifts:
  - Legacy Surveys DR10 photo-z, https://www.legacysurvey.org/dr10/
  - DESI DR1 spec-z, https://data.desi.lbl.gov/public/dr1/
  - Milliquas v8, https://heasarc.gsfc.nasa.gov/W3Browse/all/milliquas.html
- Comparison catalogues:
  - NVSS RM catalogue, VizieR J/ApJ/702/1230
  - LoTSS DR2 RM grid, https://lofar-surveys.org/
- Galactic RM model: Hutschenreuter et al. 2022 map,
  https://wwwmpa.mpa-garching.mpg.de/~ensslin/research/data/faraday2020.html

**Method (four lines).**
1. Build pairs with Vernstrom's cuts (|b| >= 20 deg, polarisation fraction >= 2 %, separation 1'-20').
   - Physical pairs: two components with one host.
   - Random pairs: both hosts with redshifts and |Delta z| > 3 sigma_z.
2. Compute the rms and median |DeltaRM| per class, noise-debiased with the catalogue sigma_RM. Match the separation
   distributions (reweighting) and the morphology (extended-extended only, the O'Sullivan caveat), and take the excess
   in quadrature.
3. Errors come from three sources:
   - a bootstrap over sources (pairs share sources);
   - nulls: label shuffles between classes in separation bins, and same-redshift random pairs;
   - mock catalogues that put the Galactic RM structure function and the catalogue noise into fake pairs.
4. Pre-register the two predictions, about 10-14 rad m^-2 (NVSS-like) against about 1-2 rad m^-2 (LoTSS-like). Then add
   the excess against the redshift of the farther source (the Pomakov null against an evolving field) as a secondary
   result.

**Overlap and feasibility.**
- **Overlap.** SPICE-RACS (dec < +49) with |b| >= 20 and Legacy DR10 photo-z: about 15,000 deg^2. Spec-z only
  (DESI DR1): about 7,000 deg^2.
- **Random pairs.** The density is about 7-8 RMs per deg^2, so a source has about 2.8 neighbours within 20'. That gives
  about 1.7e5 projected pairs. About 15-25 % have two redshifts that differ, so expect **about 25,000-40,000 random pairs
  (5-8x NVSS)**.
- **Physical pairs.** Doubles with both lobes polarised make up about 1-3 % of 2.5e5 RMs, at 15-25" resolution against
  NVSS's 45". Expect **about 1,500-4,000 physical pairs (5-12x NVSS)**.
- **Scaling NVSS's errors with sqrt(N):** sigma(rms_phys) is about 0.4-0.6 and sigma(rms_rand) about 0.2 rad m^-2. The
  catalogue noise (about 2.8 rad m^-2 per pair) is subtracted in quadrature and is known.

**Expected separation and anchor.**
- **Anchor:** the 5-sigma NVSS excess (about 10 rad m^-2, error about 1.2) against LoTSS's 1.65 +/- 0.10.
- **Separation:** with excess errors of about 0.5-1 rad m^-2, the two predictions are **> 8 sigma apart**, and still
  about 5 sigma apart with only 500 physical pairs.
- **What could erase it:**
  - Galactic RM structure on 1'-20' scales at |b| >= 20 (a few to about 10 rad m^-2) enters random pairs more than
    physical pairs if the separations aren't matched. Matching is pre-registered.
  - ASKAP off-axis leakage puts spurious RMs near 0. Cut on the leakage flags in DR2.
  - Photo-z errors mislabel random pairs as physical. Use spec-z only as a check.
  - Lobe-internal Faraday complexity. Use DR2's complexity metrics.
- **Scoop risk:** POSSUM full-survey RM grids will do this, but they aren't public.

**Quick kill (6 queries).**
- arXiv API, abs:"rotation measure" AND abs:pairs AND (physical OR filaments): only 1905.02410 (NVSS), 2002.06924 and
  2208.01336 (LoTSS).
- arXiv API, abs:"rotation measure" AND ("random pairs" OR "physical pairs" OR "non-physical"): the same three.
- arXiv API, abs:SPICE-RACS OR (POSSUM AND "rotation measure"): 26 entries. Two are close but different:
  - 2605.16924 (Mg II CGM, a 4.5 sigma sightline excess, not pairs);
  - 2610.37646 (group RM^2 stacks).
- arXiv API, au:Vernstrom: no 0.9 GHz pair follow-up.
- WebSearch (extended), SPICE-RACS/POSSUM physical vs non-physical pairs: none. POSSUM pilot 2403.15668 has only a
  nearest-neighbour |DeltaRM| quality check, with no redshifts or classes.
- WebSearch (extended), the NVSS-LOFAR pair discrepancy explanation: only 2404.17402 (LoTSS DR2); nothing at an
  intermediate frequency.

**Overlap with our work:** none. We haven't used RM data (the FRB DM x y work is different).

---

## R6-2. Large-scale alignment of radio quasar polarisation angles, re-tested with broadband, Faraday-derotated angles

**Question.** Are the intrinsic (RM-derotated) polarisation angles of radio quasars coherently aligned over ~150-500 Mpc,
as claimed from JVAS/CLASS 8.4 GHz, or randomly oriented?

**Shape: disagreement.** On the same JVAS/CLASS data:
- Tiwari & Jain 2013 (arXiv:1201.5180, 4,400 sources with P > 1 mJy) found alignment above 3 sigma at about 150 Mpc
  (n_v about 10, S_D statistic, sigma_random = 0.005), and 5.4-6 sigma in other procedures.
- Pelgrims & Hutsemekers 2015 (arXiv:1503.03482, 4,155 sources) found significant alignment in quasar-only samples.
- Joshi et al. 2007 found none.
- Pelgrims & Hutsemekers 2016 (arXiv:1604.03937) tie the angles to the axes of large quasar groups.

The disputed step is instrumental: angle biases in a single-band survey (Tiwari & Jain's own von Mises test gives a
non-uniform PA distribution at 4.7 sigma, and NVSS shows PA clustering at multiples of 45 deg). An independent
instrument, frequency and calibration, with derotated angles, removes it. The arXiv:2207.05765 review lists the effect as
an open challenge to the cosmological principle.

**Data.**
- SPICE-RACS DR2: derotated angles and RMs (arXiv:2605.16917, CASDA).
- Quasar redshifts: Milliquas v8 and DESI DR1 QSO (URLs as in R6-1).
- JVAS/CLASS polarisation (Jackson et al. 2007), VizieR J/MNRAS/376/371, for the same-source comparison.
- VLASS Quick Look, https://cirada.ca/, for spectral index (flat-spectrum selection).

**Method (four lines).**
1. Select SPICE-RACS sources matched to spectroscopic quasars, with |b| >= 20 deg and an angle error below 15 deg
   (sigma_RM x lambda0^2). Make two subsets: flat-spectrum (alpha from RACS-VLASS) and all.
2. Compute the S_D / Z_D statistics with parallel transport in 2D and 3D, at the published n_v (2-2500) and scales.
3. Get significance from three nulls:
   - angle shuffles among sources inside each ASKAP tile and beam position, which keeps any instrumental pattern;
   - a beam-position angle regression;
   - mocks with leakage-like angle patterns injected.
4. Same-source check: for JVAS/CLASS sources also in SPICE-RACS, compare the 8.4 GHz angles with the derotated
   0.9 GHz angles, and rerun the JVAS statistic on the derotated angles.

**Overlap and feasibility.**
- **Footprint.** SPICE-RACS (dec < +49) covers the whole of the quasar sky the claim uses below dec +49.
- **Same-source subset.** JVAS/CLASS lies at 0 < dec < 75, so the shared strip is dec 0-49. That holds about 2,500
  JVAS/CLASS sources, about 1,000 of them polarised in SPICE-RACS. This is a 1-2 sigma instrument check only.
- **Main test.** About 5,000-15,000 SPICE-RACS sources match spectroscopic quasars (estimate for the gate;
  arXiv:2605.16924 already matched 2,483 DESI DR1 quasars after strict cuts).
- **Angle errors.** Per-source derotation error is about 2 rad m^-2 x 0.1 m^2, roughly 11 deg, small against the
  uniform 52 deg rms.

**Expected separation and anchor.**
- **Anchor:** the JVAS effect is above 3 sigma at N = 4,400, sigma(S_D) = 0.005, so Delta S_D is about 0.015.
- **Separation:** at N = 10,000, sigma(S_D) is about 0.0033, so the effect at its claimed amplitude would show at about
  **4.5 sigma**. N = 5,000 gives about 3.2 sigma.
- **What could erase it:**
  - Depolarisation selects against compact flat-spectrum cores at 0.9 GHz, so the sample may be lobe-dominated, a
    different population from JVAS. Report the flat-spectrum subset separately.
  - ASKAP leakage makes angle patterns on beam and tile scales (about 1-6 deg, i.e. tens to hundreds of Mpc at
    z ~ 1). The within-tile shuffle null controls this, at the cost of power on those scales.
- **Note:** a physical prior for real alignment is weak. The value is in settling a twenty-year anomaly with an
  independent instrument.

**Quick kill (5 queries).**
- arXiv API, polarization AND alignment AND (radio OR quasars) AND "large scale": 45 hits. The latest radio tests are
  2015-2016 (1503.03482, 1604.03937); 2009.01637 covers LQGs at 2.4 sigma.
- arXiv API, polarization AND alignments AND (quasars OR blazars OR "radio sources"): no test since 2017.
- arXiv API, ("polarization angles" OR "polarization vectors") AND alignment AND radio: only ELAIS-N1 fields and single
  sources.
- WebSearch (extended), SPICE-RACS polarisation-angle alignment: none.
- WebSearch (extended), NVSS/VLASS/LoTSS polarisation alignment 2020-2026: none. The LOFAR 2020 paper (Osinga,
  2008.10947) is about jet position angles, a different question, which round 4 already answered.

**Overlap with our work:** sky anisotropy in counts and dipoles (tried) is a different statistic, and this candidate
uses no counts.

---

## R6-3 (marginal). HI spin-filament alignment from WALLABY pilot kinematic models

**Question.** Are the 3D HI spin axes of gas-rich galaxies preferentially parallel to their nearest filament,
mean |cos psi| about 0.65 as MIGHTEE-HI found, or at most weakly aligned (about 0.5-0.55) as simulations and
stellar-kinematic surveys find?

**Shape: disagreement.**
- MIGHTEE-HI (Tudorache et al. 2022, arXiv:2204.03041): 77 galaxies, <|cos psi|> = 0.66 +/- 0.04 within 5 Mpc of a
  filament, called "tentative". It comes from one COSMOS field.
- A 15 Mpc HI filament (arXiv:2508.13053): 0.64 +/- 0.05, "significantly more strongly aligned than simulations
  predict".
- The other side:
  - MaNGA stellar spins show no alignment (Krolewski et al. 2019, arXiv:1902.09797) or a mass-dependent anti-parallel
    tendency (arXiv:2506.22794);
  - SIMBA predicts only a weak parallel trend for HI-rich galaxies (arXiv:1906.01623).

The disputed step is single-field cosmic variance, together with the tracer (HI against stars). Independent HI
kinematics from other fields remove the first.

**Data.**
- WALLABY pilot kinematic models, about 235 in all:
  - PDR1, 109 models (Deg et al. 2022, PASA). The phase 1 HI data release is arXiv:2211.07094.
  - PDR2, 126 models (arXiv:2409.13130).
  - Access: CASDA project AS102, https://wallaby-survey.org/data/
- Filaments:
  - DisPerSE run on DESI DR1 BGS (https://data.desi.lbl.gov/public/dr1/) for the NGC 4636, NGC 4808 and NGC 5044 fields;
  - Tempel et al. 2014 SDSS filaments, VizieR J/MNRAS/438/3465;
  - 2MRS / 6dFGS for Hydra.

**Method (four lines).**
1. Build the 3D spin vector from the kinematic PA, the inclination and the approaching side. Carry both inclination
   solutions.
2. Find the nearest filament segment and psi, with distance cuts of < 5 Mpc and 5-10 Mpc (MIGHTEE's bins).
3. Get the null from spin-axis randomisation and filament-catalogue rotations. Take errors from a bootstrap over
   galaxies and a jackknife over fields.
4. Compare with MIGHTEE's 0.66 +/- 0.04 and with a mock built by applying the same selection to a simulation (SIMBA or
   TNG HI-rich galaxies).

**Overlap and feasibility.**
- **Usable sample.** Norma and Vela lie in the Galactic plane and have no usable filament catalogue. That leaves about
  150 models in fields with DESI/SDSS/6dF coverage, and about 100 within 5 Mpc of a filament.
- **Precision.** sigma(<|cos psi|>) = 0.289/sqrt(100), about 0.029.
- **Fields.** NGC 4636 and NGC 4808 lie in Virgo's southern filaments (DESI DR1 and SDSS cover them); NGC 5044 is at
  dec -16 (DESI DR1 partly).

**Expected separation and anchor.**
- **Anchor:** MIGHTEE 0.66 +/- 0.04, 4 sigma from random.
- **Separation:** MIGHTEE's value against random (0.50) is about 5.5 sigma apart. Against a simulation-like 0.55 it is
  about **3.8 sigma**, falling to about 3 sigma if only 70 galaxies are usable.
- **What could erase it:**
  - The pilot fields are centred on groups and clusters, where infall and harassment can randomise spins.
  - Kinematic PAs aren't corrected for beam smearing.
  - The 2-fold inclination ambiguity dilutes the signal.
- **Marginal on N.** The full-survey release (about 800 models, team-internal as of the Feb 2026 newsletter) would make
  it decisive. Run it now only if the gate's N count reaches 100 or more.

**Quick kill (4 queries).**
- WebSearch (extended), WALLABY HI spin alignment with filaments: none. The WALLABY survey paper (2002.07311) only
  plans it.
- arXiv API, all:WALLABY AND spin AND filament: 0 results.
- arXiv API, spin AND filament AND HI AND alignment: 2204.03041, 2508.13053, 1912.01062 (CHILES, 10 galaxies),
  2601.16408 (FAST, 16 galaxies, one filament) and 2411.14940. None of them is WALLABY or a multi-field sample.
- WebSearch (extended), Virgo-region HI kinematic spin alignment (WALLABY/ASKAP/MeerKAT): none. 2610.39987 uses MDPL2
  simulations only.

---

## Killed (with reason)

| Idea | Reason | Evidence |
|---|---|---|
| IceCube tracks x LSS (DESI, unWISE, eRASS1) | Done, null | 2510.18119, 2507.14926 (1.9 sigma max), 2405.09633, 2406.00796 |
| IceCube stack on eRASS1 clusters | Same question as the Planck SZ stack, new catalogue, no dispute | 2206.02054, 2107.10080 |
| ACT DR6 lensing x RACS/EMU radio sources | EMU x Planck done at 5.5 sigma, consistent; no disagreement to settle | 2411.05913; LoTSS x eBOSS / EMU x DES n(z), round 4 |
| LoTSS DR3 radio-halo stacks of eRASS1 clusters (the footprint now overlaps) | Cluster CR physics rather than cosmology; a DR3 cluster/diffuse-emission catalogue paper exists | 2605.06400 |
| f_NL from LoTSS DR3 / RACS counts | sigma(f_NL) of about 50 can't compete; large scales overlap our dipole work | round 2 note |
| PTA: MPTA amplitude excess; GWB anisotropy; HD by sky sector | Infeasible (about 1.5 sigma) or done | round 3; 2306.16221; sweep-2 |
| Auger open data, UHECR x LSS | The open set is a fraction of events; underpowered against the collaboration's own dipole | reasoning only, no query spent |
| Radio PA-jet offsets as an isotropic birefringence test against Planck beta = 0.34 +/- 0.09 deg | Needs about 0.1 deg; absolute ASKAP/VLA angle calibration is about 1 deg | reasoning |
| FASHI / WALLABY Tully-Fisher bulk flow | Overlaps our CF4 flow work and the DESI PV survey (2512.03230) | round 4 |
| Radio jet position-angle alignment with DESI 3D | Answered (round 4) | 2008.10947 |
