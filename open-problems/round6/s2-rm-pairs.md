# Proposed test: Faraday-rotation excess of random over physical pairs at 0.9 GHz (SPICE-RACS DR2)

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
