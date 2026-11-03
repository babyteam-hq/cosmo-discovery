# Analysis plan: how much weaker would a real Giant Arc be on new sightlines than on the members that defined it?

## What is already measured (not changed by this plan)

The Giant Arc (GA; Lopez, Clowes & Williger 2022) is a group of 44-46 Mg II absorbers at z ≈ 0.80 found by single-linkage
clustering (95 Mpc, comoving, flat Ωm 0.27, h 0.7) of the Zhu & Ménard (Z&M) SDSS DR7 + DR12 absorber catalogues in the
slice z = 0.802 ± 0.060 (4,752 absorbers). Its published overdensity, δρ/ρ = 1.3 ± 0.3, is measured inside the convex hull
of its member spheres (CHMS; sphere radius half the mean minimum-spanning-tree edge, 30.4 Mpc) on the same absorbers the
algorithm selected.

A rebuild of the GA with that recipe (46 members) was frozen, and DESI DR1 Mg II absorbers were counted on 3,822 DESI quasar
sightlines through its CHMS that SDSS never observed. Expected counts come from DESI dN/dz per stratum of quasar flux,
exposure depth and redshift, measured on all independent sightlines. Result, already in hand: 8 absorbers against 8.86
expected; after referring the GA to its own shift null (the same shape moved rigidly to ~150 places in the SDSS footprint;
its spread, 0.34, is the error and carries Poisson noise, absorber clustering and coverage together) and subtracting a chance
baseline from the slice's other linkage groups (+0.22 ± 0.34), the fitted amplitude of a uniform excess in the CHMS is
A_obs = −0.23 (68 %: −0.71 to +0.26; 95 %: −0.90 to +0.63).

## The question this plan answers

A_obs measures the overdensity a structure in the GA's hull carries on sightlines that played no part in finding it. The
published 1.3 is measured on the members the algorithm picked, which upward fluctuations put there, so even a real structure
would come back lower on new sightlines (selection shrinkage), and the recovered hull need not coincide with where the real
structure is. Comparing A_obs with 1.3 therefore does not say whether the GA is real. The question is: **if the GA were a
real overdensity that, run through the same linkage recipe on the real SDSS sightlines, produces a group like the GA (its
size, extent and member overdensity), what amplitude would that predict on the new sightlines, and does A_obs exclude it?**

## Method: injection into the real SDSS catalogue

- **Background.** The real Z&M slice (4,752 absorbers) on the real SDSS sightlines: real clustering, real coverage and real
  chance groups, with no model of either.
- **SDSS sightlines and path.** The Z&M searched-quasar lists (DR7 ∪ DR12, de-duplicated at 1.5″). A sightline's path in the
  slice runs where Mg II is redward of the quasar's Lyα emission and more than 5000 km/s below z_q. Expected SDSS counts
  E_SDSS use one slice dN/dz per unit path (N_slice / total path); the per-sightline completeness of SDSS is not modelled
  (no S/N column in the searched lists).
- **Member overdensity estimator.** δ_mem = N(CHMS ∩ slice) / E_SDSS(CHMS ∩ slice) − 1, counting every slice absorber in the
  group's own CHMS. Computed on the GA rebuild first (δ_GA), so the injections are conditioned on the value our estimator
  gives the GA, not on a published number our estimator may not reproduce. The published 1.3 ± 0.3 is reported beside it.
- **Injection sites.** The GA's frozen CHMS shape moved rigidly (rotation on the sky, random spin) to random places in the
  SDSS footprint, not overlapping the GA, where E_SDSS(CHMS ∩ slice) is within 0.75-1.33× the GA's.
- **Injected structure.** A real overdensity of amplitude A_true: Poisson(A_true × E_SDSS(region ∩ slice)) extra absorbers
  placed uniformly per unit SDSS path inside the region ∩ slice. Two shapes: profile H, uniform in the GA's CHMS shape (the
  reading under test), and profile S, uniform only in the union of the 30.4 Mpc member spheres (a clumpy arc). A_true grid:
  0 (no injection: chance groups only), 0.25, 0.5, 0.75, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0; 200 sites per value and profile.
- **Recovery.** Linkage at 95 Mpc on the slice plus the injected absorbers. The recovered group is the group of ≥ 20 members
  holding the most injected absorbers. It is "GA-like" if it has 30-70 members, a largest member separation of 0.5-2× the
  GA's 742 Mpc, and δ_mem within ±0.3 of δ_GA.
- **Prediction on new sightlines.** For each GA-like recovery: A_pred = A_true × E_DESI(CHMS_rec ∩ R_inj ∩ slice) /
  E_DESI(CHMS_rec ∩ slice), with the recovered hull and the injected region rotated back onto the GA's position, so the DESI
  weights are those of the GA's actual 3,822 independent sightlines and their strata. This is the amplitude the DESI
  measurement would recover in expectation, the quantity A_obs estimates. (Chance clustering in the recovered hull is already
  removed in A_obs by the baseline.)
- **Inference.** Likelihood of the GA's member overdensity under each A_true: L(A_true) = the fraction of injections giving a
  GA-like recovery. Flat prior on A_true in [0, 3]. The real reading's predictive distribution for A is the L-weighted
  mixture of A_pred over A_true > 0, convolved with the measurement's own error distribution (the GA shift-null spread
  combined in quadrature with the baseline's standard error, as A_obs was). Reported: the predicted A (median, 68 %), the
  one-sided probability of A ≤ A_obs under it, and the same for chance (A_true = 0: GA-like groups among the uninjected
  sites). The two shapes are reported separately.
- **Decision fixed before computing.** If the real reading's predicted A sits ≥ 2σ (one-sided p ≤ 0.023) above A_obs for
  profile H, the claim becomes "a real structure strong enough to produce the GA's member overdensity through this recipe is
  excluded at that level on new sightlines". If not, the claim becomes "this test cannot distinguish a real Giant Arc from
  chance", and the finding is withdrawn.

## Error budget and what the result depends on

- A_obs's error: the shift-null spread (Poisson, clustering and coverage), plus the baseline's standard error; not Poisson
  alone.
- The prediction depends on: the injected shapes (two are run), the SDSS per-sightline completeness (uniform), the GA-like
  window (size, extent, δ_mem ± 0.3; widened to ± 0.5 as a variant), and the flat prior (a prior ∝ 1/A_true is a variant).
- The expected DESI count in the hull, 8.86, is the sum over those 3,822 sightlines of each one's path in the hull times the
  DESI dN/dz of its stratum.

## Method check (blind-review.md 1b; 2026-11-03, fresh context-free Opus, this plan alone), verbatim

METHOD: fixable: replace the fixed A=0 error convolution with count-level Poisson(E_DESI(1+A_pred)) plus clustering scatter per injection; calibrate the chance baseline conditional on GA-like recovery under each hypothesis (or report results with and without it) and check it on A_true=0 GA-like groups; match DESI to the Z&M equivalent-width/completeness selection or measure/vary the DESI/SDSS relative bias; fix minimum GA-like recovery counts and independent-site counts with bootstrap errors on L and p and check L(3)≈0; add an SDSS completeness-weighted injection variant; restrict the exclusion claim to profile H unless S (and a larger diffuse profile) also exclude; exclude known neighbouring groups from injection sites and track merges

How each point was closed before the injections ran (`garev_inject.py`, `garev_infer.py`):

1. **Count level.** Each GA-like injection predicts the hull's DESI count as Poisson(8.86 (1 + bg + A_true f_rec)), times a
   lognormal for clustering beyond Poisson (width 0.058: the GA shift null's 0.341 against Poisson's 0.336 at 8.86). Compared
   with the 8 observed as P(N ≤ 8), not with a Gaussian.
2. **Chance part, per hypothesis, measured.** bg is real DESI N_obs / N_exp in the recovered hulls at the injection sites
   (the injected absorbers are not in DESI), pooled over the GA-like recoveries under each hypothesis; for chance it is also
   checked on the real slice's own GA-like groups. No baseline is subtracted from the GA's count: the comparison is on the
   raw 8 against 8.86, with bg carrying the model's offset and the chance clustering of a GA-like hull.
3. **DESI against SDSS.** The EW cut is already on the Z&M scale (DESI EW / Z&M REW = 1.03 on 21,789 shared absorbers), and
   92 % of DESI's W > 0.3 sample is above 0.6 Å; W > 0.6 is run as the check. The ratio r = δ_DESI / δ_SDSS is scanned
   (0.25-1.25) and the result reported against it.
4. **Monte Carlo.** 200 sites; independent ones counted (greedy, centres > 20° apart); errors on P from a block bootstrap
   over those independent patches. A grid extended to 4.0 to check L at the top. A hypothesis's P is reported only with ≥ 20
   GA-like recoveries.
5. **Completeness variant.** Profile H also run with SDSS placement weighted by the quasar's g magnitude (DR16Q PSF g, seven
   quantile classes plus unmatched, each with its own measured slice dN/dz), conditioned on δ_GA under that model.
6. **Shapes.** Profiles H, S and L (the hull of 60.8 Mpc member spheres: larger, more diffuse) all run; an exclusion is
   claimed only for the shapes that give it.
7. **Neighbours.** Sites are rejected when the moved hull (at twice the sphere radius) holds ≥ 10 members of any real slice
   group of ≥ 40 members; each recovery records how many members it drew from one pre-existing real group.

Disclosure: A_obs and the DESI count were known before this plan was written (R-155). The plan, the GA-like window and the
decision rule were fixed before any injection was run.
