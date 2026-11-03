# Results (R-155, cosmo-deep-23, 2026-11-02)

Does the Giant Arc (GA; Lopez, Clowes & Williger, arXiv:2201.06875) reappear as an excess of Mg II absorbers on
*independent* sightlines, DESI DR1 quasars that SDSS never observed? Proposal: `open-problems/round3/q5-giant-arc.md`
(trajectory check: "VERDICT: novel: no prior work answers this question"). This file is the worker's own.

**Disclosure.** The harvester that proposed this test ran an unblinded rough count at proposal time (delta = -0.16 / -0.28 /
-0.29 at R = 47.5 / 70 / 95 Mpc; partial SDSS exclusion, no EW matching, no TSNR2 weighting). The worker saw it in the
brief. The template, samples, statistic, nulls and primary cell below were fixed before the worker's own first look at the
GA's observed count, and nothing was changed after that look.

Scripts: `r155_common.py` (cosmology, matching), `r155_template.py` (GA rebuild and freeze -> `r155_templates.json`),
`r155_sightlines.py` (DESI sightlines and absorbers, SDSS removal), `r155_core.py` (path, completeness strata, counts),
`r155_run.py` (feasibility / measure, nulls). Data under `data/r155/` (gitignored), logs `logs/r155_*.log`.

## F0. Access (fetch times; `data/` started empty)

| file | source | size | time |
|---|---|---|---|
| MgII-Absorbers-DR1.fits | data.desi.lbl.gov/public/dr1/vac/dr1/mgii-absorber/v1.0 | 66 MB | ~40 s |
| QSO_cat_iron_cumulative_v0.fits (the VAC's parent "afterburner" QSO catalogue) | data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/QSO/iron | 816 MB | ~1 min |
| DR16Q_v4.fits | data.sdss.org/sas/dr16/eboss/qso/DR16Q | 1.7 GB | 5.5 min |
| Z&M Trimmed_SDSS_DR7_107, Trimmed_BOSS_DR12_107, Searched_Quasars_* | s3.us-east-2.amazonaws.com/jhu-sdss-qsoabs | 22 MB | seconds |

All fetches done in 331 s. Python deps (astropy, scipy) installed with pip in ~1 min.

## F1. The GA rebuild, frozen before any DESI file was opened (`r155_template.py`)

Z&M DR7 (35,752) + DR12 (41,895); 8,441 DR7 absorbers duplicate a DR12 one (sightline within 1.5", |dz| < 0.005) and are
dropped: 69,206 merged. Slice z = 0.802 +- 0.060: 4,752 absorbers. Single linkage at 95 Mpc (flat Om 0.27, h 0.7): 1,322
groups. **The GA group: 46 members, RA 205.0-224.6, Dec 36.2-42.2, z 0.745-0.860, W2796 median 0.82 A (10-90 %:
0.25-1.96)**, the same as the proposal's (the paper has 44). It is not the slice's largest group: the six largest have 88-131
members. Members in `r155_templates.json`.

## F2. Sightlines, absorbers, EW scale (`r155_sightlines.py`)

- DESI main-dark quasars, one row per TARGETID: 1,556,507. Every VAC absorber's quasar is main-dark. The VAC's TSNR2_QSO equals
  the sum of the catalogue's per-tile rows for 98.9 % of VAC quasars, so TSNR2 = that sum.
- Within 1.5" of a DR16Q quasar or a Z&M searched quasar: 322,606. **Independent sightlines: 1,233,901.**
- The VAC does not search the Lya forest: dN/dz at z 0.7-0.9 is 0.06-0.09 for z_q < 3.4 and exactly 0 above. Path starts at
  z_abs = 1250 (1 + z_q)/2796 - 1, ends 5000 km/s below z_q.
- **EW scale** (`r155_ewcal`, inline check): 21,789 absorbers sit in both DESI (on re-observed SDSS sightlines) and Z&M (1.5",
  |dz| < 0.003). DESI EW_2796 / Z&M REW_MGII_2796 = 1.03 median (16-84 %: 0.85-1.21), flat across EW (1.02-1.05 in four
  bins) and 1.03 in the slice. The matched cuts W > 0.3 / 0.6 A on the Z&M scale are DESI > 0.309 / 0.618 A.
- All-sky dN/dz in the slice on independent sightlines: 0.0599 (W > 0.3, 7,864 absorbers), 0.0548 (W > 0.6). DR1 detects few
  systems below 0.6 A.

**Completeness is set by the quasar's flux, not the exposure depth.** dN/dz (W > 0.3, z = 0.8) by g-flux quintile:
0.00004 / 0.001 / 0.011 / 0.049 / 0.280, a ~6000x range; by TSNR2 tertile only 0.052 / 0.063 / 0.090. Removing SDSS quasars
removes the brighter ones where SDSS observed, so inside the SDSS footprint the independent sightlines are fainter: through
the GA template (R = 95) 8 % sit in the brightest g quintile, against 20 % sky-wide. A flat dN/dz therefore overpredicts the
GA's expected count by ~2x (path dz 595 x 0.060 = 36 against the stratified 16.4), enough on its own to make delta about -0.5.
**The proposal's rough count (flat dN/dz 0.074) is biased low for this reason**, so its negative delta says nothing.
Expected counts here use dN/dz measured per stratum (TSNR2 tertile x 12 g-flux bins, finest at the bright end x 4 z_q bins);
the all-sky slice check gives N_exp 7,844 against N_obs 7,864.

## F3. Step 0: feasibility (`r155_run.py feasibility`; the GA's observed count is neither computed nor printed)

| template | sightlines | path dz | N_exp (W>0.3) | sigma_shift (n) | shift mean | clust mean (sd, median N_exp) | separation of 1.3 vs baseline |
|---|---|---|---|---|---|---|---|
| R = 47.5 | 3,534 | 152 | 4.13 | 0.499 (378) | -0.045 | +0.092 (0.31, 12.2) | 2.6 sigma |
| R = 70 | 4,807 | 350 | 9.64 | 0.316 (380) | -0.052 | +0.013 (0.18, 37.5) | 4.1 sigma |
| R = 95 | 6,282 | 595 | 16.39 | 0.240 (366) | -0.073 | -0.010 (0.10, 87.7) | 5.4 sigma |
| CHMS hull | 3,823 | 302 | 8.86 | 0.342 (346) | -0.072 | n/a | 3.8 sigma |

W > 0.6 is within 5 % of W > 0.3 everywhere. Against the proposal: its sightline counts (5,711 / 7,795 / 10,104) and N_exp
(17.8 / 39.1 / 62.1) are higher because it removed only Z&M sightlines (not DR16Q), used the LSS NGC clustering catalogue
(z_q > 0.9), and used one flat dN/dz; its null sigma (0.29 at median N_exp ~12) matches ours at equal N_exp. The clustering
baseline over random template positions (clust - shift) is +0.14 / +0.07 / +0.06, as the proposal's +0.07.

**Gate: passed.** "Real" (delta ~ 1.3 above baseline) and "chance" (delta ~ baseline) separate by >= 3 sigma at R = 70, 95 and the
hull; not at R = 47.5 (2.6 sigma), which is reported but cannot decide alone.

**Big Ring and Giant Ring: not run.** Neither paper's abstract (2402.07591, 2604.17534) gives a member list or a reproducible
membership recipe (the Big Ring is a CHMS/visual ring from Anand et al.'s DR16Q catalogue; the Giant Ring uses FilFinder and
ellipse matching). A template frozen from them would be ours, not theirs. Left as a loose end.

## F4. Method check (blind-review.md 1b; 2026-11-02, fresh context-free Opus, plan alone in `r155_plan.md`)

Verdict, verbatim:

> METHOD: fixable: derive the real-structure Δ per cell by injecting the published overdensity (CHMS and profiles, slice versus 0.70–0.90 path) through the real DESI sightlines, and re-pick the primary cell from that; replace the random-46 clustering null with the linkage recipe run on clustered SDSS-selected mocks (or promote null 3, matched in size and extent) as the chance baseline; use the same mocks to quantify the selection-bias shrinkage so an intermediate Δ has a defined reading; match the chance baseline in footprint, SDSS fraction and N_exp to the GA, or refer each to its own shift null, so the −0.07 calibration offset cancels; exclude DESI sightlines within about 200 proper kpc of member sightlines (and vary that radius) in the GA and in all nulls; carry the amplitude over to DESI's effective EW distribution or check it at W > 0.6

The referee's other points: no change was needed for Poisson over-dispersion (the shift null captures it), the fixed cosmology,
or the template freeze. Each point and how it was closed before any GA count was seen (`r155_measure.py`):

1. **Real-structure prediction.** Computed in expectation through the real sightlines and paths (`r155_measure.py predict`,
   `r155_predict.json`). Profile H: A = 1.3 uniform in the CHMS cut to the slice. Profile S: the same excess in member spheres
   of r_sph = 30.4 Mpc cut to the slice (amplitude 9.3). The primary cell was re-picked from this (F5).
2. **Chance baseline.** Null 3 is promoted. The baseline is the 18 other single-linkage groups from the same slice and recipe,
   matched in size (23-92 members) and extent (largest separation 371-1483 Mpc, 0.5-2x the GA's). Lognormal/N-body mocks were
   not built (budget). The random-46 null stays as a secondary row.
3. **Selection-bias shrinkage.** Closed by the reading, not by mocks. The claim is the fitted amplitude A = Delta / f_H of a
   real structure's overdensity in the CHMS, so any value has a defined meaning: 0 for a chance pattern, 1.3 as published, and
   in between "real but weaker than its selected members imply". How much shrinkage selection alone would produce is not
   modelled, so an intermediate A cannot be split into "real" and "inflated by selection". Recorded as a limit.
4. **Calibration offset.** Every template (the GA and each comparison group) is referred to its own shift null (the same shape
   in the SDSS footprint, N_exp 0.75-1.33x its own): delta_cal = delta - mean(shift).
5. **Close pairs.** DESI sightlines within 200 proper kpc of a member's sightline at its redshift are dropped in the GA and in
   every null. 0 and 500 kpc are reported as variants.
6. **EW mix.** DR1 detects few systems below 0.6 A (all-sky dN/dz 0.0599 at W > 0.3 against 0.0548 at W > 0.6, so 92 % of the
   W > 0.3 sample is above 0.6 A). The DESI sample's EW mix is therefore close to the SDSS members' (median 0.82 A). W > 0.6
   is reported as the amplitude check.

## F5. Prediction and primary cell, fixed before the measurement (`r155_measure.py predict`, 200 kpc exclusion)

| cell (W > 0.3) | sightlines | N_exp | f_H | pred. Delta, profile H | pred., profile S | sigma_shift (n) | separation H / S |
|---|---|---|---|---|---|---|---|
| R = 47.5 | 3,533 | 4.13 | 0.878 | 1.14 | 2.58 | 0.495 (159) | 2.3 / 5.2 |
| R = 70 | 4,806 | 9.64 | 0.730 | 0.95 | 1.11 | 0.334 (167) | 2.8 / 3.3 |
| R = 95 | 6,281 | 16.39 | 0.515 | 0.67 | 0.65 | 0.272 (163) | 2.5 / 2.4 |
| **CHMS** | 3,822 | 8.86 | 0.976 | **1.27** | 1.20 | 0.369 (158) | **3.4** / 3.3 |

W > 0.6 gives the same within 0.02 in f and prediction. The referee was right about R = 95: its spheres reach z 0.70-0.90 and
far past the hull, and dilute a real excess to 0.67. **Primary cell: CHMS, W > 0.3** (the largest profile-H separation,
3.4 sigma). V(CHMS & slice) = 3.5e7 Mpc^3.

## F6. The measurement (`r155_measure.py measure`, `logs/r155_measure.log`; 200 kpc exclusion, which drops no sightline)

| cell | N_obs | N_exp | δ (TSNR2-w) | GA shift off (sd, n) | baseline (n groups) | Δ [68 %] [95 %] | A [95 %] | vs 0 | vs published |
|---|---|---|---|---|---|---|---|---|---|
| **CHMS W>0.3** | **8** | **8.86** | −0.097 (−0.190) | −0.100 (0.341, 148) | +0.223 ± 0.340 (4) | **−0.22 [−0.69, +0.25] [−0.88, +0.62]** | **−0.23 [−0.90, +0.63]** | −0.5σ | **3.2σ** |
| CHMS W>0.6 | 8 | 8.40 | −0.048 (−0.147) | −0.068 (0.337) | +0.116 ± 0.360 (4) | −0.10 [−0.59, +0.36] [−0.77, +0.69] | −0.10 [−0.79, +0.70] | −0.2σ | 2.9σ |
| R47.5 W>0.3 | 4 | 4.13 | −0.032 (−0.060) | −0.085 (0.467) | +0.070 ± 0.251 (4) | −0.02 [−0.52, +0.50] [−0.97, +0.98] | −0.02 [−1.10, +1.11] | 0.0σ | 2.3σ |
| R47.5 W>0.6 | 4 | 3.90 | +0.027 | −0.067 | −0.048 ± 0.336 (4) | +0.14 [−0.85, +1.23] (95) | +0.16 | 0.3σ | 1.8σ |
| R70 W>0.3 | 6 | 9.64 | −0.378 (−0.428) | −0.113 (0.336) | +0.065 ± 0.193 (5) | −0.33 [−0.73, +0.05] [−0.91, +0.41] | −0.45 [−1.25, +0.56] | −0.8σ | 3.3σ |
| R70 W>0.6 | 6 | 9.14 | −0.343 | −0.090 | +0.012 ± 0.227 (5) | −0.27 [−0.86, +0.49] (95) | −0.36 | −0.7σ | 3.0σ |
| R95 W>0.3 | 15 | 16.39 | −0.085 (−0.093) | −0.111 (0.265) | −0.001 ± 0.149 (5) | +0.03 [−0.28, +0.31] [−0.43, +0.63] | +0.05 [−0.84, +1.23] | 0.1σ | 2.2σ |
| R95 W>0.6 | 15 | 15.58 | −0.037 | −0.094 | −0.032 ± 0.167 (5) | +0.09 [−0.42, +0.68] (95) | +0.17 | 0.3σ | 2.0σ |

Exclusion radius 0, 200, 500 kpc: identical counts in every cell. GA rank among the calibrated groups: 0.20-0.50.

**Baseline sensitivity (not pre-fixed; reported because the fixed rule kept only 4-5 groups).** The groups whose own shifts
never matched (0 of up to 1,500 tries) are the well-covered ones outside the SDSS-thinned area. Using all covered groups'
raw δ: CHMS baseline +0.128 ± 0.128 (14 groups), Δ = −0.13 ± 0.37, A = −0.13, published excluded at 3.8σ; R95: +0.07 ± 0.28, 2.2σ.
No baseline: CHMS A = 0.00, 3.7σ. Random-46 clustering null (feasibility run, 200 sets): +0.14 / +0.07 / +0.06 above the shift
mean at R = 47.5 / 70 / 95.

## F7. Finding check (2026-11-02, fresh context-free Opus, the `finding.md` section alone)

VERDICT: the findings are confirmed novel and sound if the calculations hold

Strongest objection (referee): the published 1.3 is selection-inflated, so the test falsifies the published amplitude but
cannot separate a modest real structure (A ≲ 0.6) from a chance pattern without the selection shrinkage from clustered mocks.
**Pass on the first finding check.** No revision made.

## F8. Re-check, one revision, withdrawn (R-171, cosmo-ga-rev, 2026-11-03)

**The re-check.** R-170 (2026-11-03) re-ran the finding check on the author email after the page's title was reworded:
VERDICT: novel but not sound: the 3.2σ is a comparison against a member-selected overdensity that even a real structure would not reproduce on new sightlines, and without mock calibration of that shrinkage (and with Poisson-only errors on ~9 expected absorbers) the test cannot discriminate a real Arc from chance, so the title's implied tension is unsupported.

The objection, in one line: the published 1.3 is measured on the members the algorithm picked, so A_obs against 1.3 is not
a test of whether the Arc is real; what a real Arc predicts on new sightlines is unmeasured. It is the first finding check's
strongest objection (F7), and physics/method, not presentation: no rewording makes the comparison answer the question.
Computable on public data: inject real overdensities into the real SDSS catalogue, re-run the recipe, keep the GA-like
recoveries, and read off what DESI would see. Forecast before computing: the test's σ(δ) ≈ 0.34 at 8.86 expected, so it
discriminates only if the shrinkage leaves a real Arc's excess ≳ 0.7.

**Method check on the revision's plan** (blind-review.md 1b; fresh context-free Opus, `ga_rev_plan.md` alone), verbatim:
METHOD: fixable: replace the fixed A=0 error convolution with count-level Poisson(E_DESI(1+A_pred)) plus clustering scatter per injection; calibrate the chance baseline conditional on GA-like recovery under each hypothesis (or report results with and without it) and check it on A_true=0 GA-like groups; match DESI to the Z&M equivalent-width/completeness selection or measure/vary the DESI/SDSS relative bias; fix minimum GA-like recovery counts and independent-site counts with bootstrap errors on L and p and check L(3)≈0; add an SDSS completeness-weighted injection variant; restrict the exclusion claim to profile H unless S (and a larger diffuse profile) also exclude; exclude known neighbouring groups from injection sites and track merges

All seven closed before any injection ran; how each was closed is in `ga_rev_plan.md`.

**The computation** (`garev_inject.py`, `garev_infer.py`; logs `logs/garev_run.log`, `logs/garev_infer.log`; 8,800
injections in 20 min on 3 cores). The member-overdensity estimator gives the GA rebuild 1.32 (48 in its CHMS against 20.7
SDSS-expected), the published 1.3. 200 sites (633 tried), 14 independent patches.

| injected shape | GA-like recoveries | A_true (16/50/84) | predicted δ in hull (16/50/84) | P(N ≤ 8 \| real) | P(N ≤ 8 \| chance) | ratio real/chance |
|---|---|---|---|---|---|---|
| **H, uniform in the CHMS (primary)** | 280 | 0.5 / 1.3 / 2.0 | +0.25 / +0.62 / +1.05 | **0.104 (1.26σ)**, block bootstrap 0.091-0.114 | 0.32 | 0.37 |
| S, member spheres | 170 | 1.6 / 3.0 / 4.0 | +0.15 / +0.26 / +0.42 | 0.22 (0.77σ) | 0.24 | 0.91 |
| L, 60.8 Mpc hull | 200 | 0.5 / 1.0 / 1.3 | +0.22 / +0.56 / +1.05 | 0.12 (1.19σ) | 0.32 | 0.42 |
| H, g-magnitude completeness | 244 | 0.5 / 1.0 / 1.6 | +0.19 / +0.48 / +0.86 | 0.13 (1.11σ) | 0.33 | 0.46 |

Primary-shape variants: window ±0.5, 0.11; prior 1/A, 0.15; W > 0.6, 0.13; r = δ_DESI/δ_SDSS 0.25-1.25, 0.28-0.08. Chance
part bg: +0.04 (real), +0.15 (chance); the slice's own GA-like groups +0.20. f_rec median 0.54 (H). L(4.0) = 0 for H and L.

**Decision, by the rule fixed with the plan:** P = 0.10 > 0.023, so the claim becomes "this test cannot distinguish a real
Giant Arc from chance", and the finding is withdrawn. Section: `finding-ga-rev.md`.

**Finding check on the revised section** (fresh context-free Opus, `finding-ga-rev.md` alone), verbatim:
VERDICT: the findings are confirmed novel and sound if the calculations hold

Strongest objection (referee): the real-Arc prediction rests on uniform or sphere-clumped injections, not the thin filament
the discoverers describe, so "cannot tell" is shown for those shapes; the DR2 forecast inherits that and should be redone
with a filament injection. Closest work: Sawala & Teeriaho 2025 (arXiv:2505.11072), which the note should cite.

The revised note passes, but its claim is that the test cannot decide; by the rule fixed before computing, the page is a
withdrawal notice (`site/research-giant-arc.html`), not a paper. Whether a "cannot decide" note may stand as a paper is a
founder's call, recorded in the work log; the notice is reversible.

## F9. Published as its interval (R-172, 2026-11-03)

The dispatcher overrode R-171's withdrawal default: the plan's rule governs what is claimed (no exclusion of the Arc), not
whether a sound result is published, so the revision that passed (F8) is the paper at `site/research-giant-arc.html`
(`briefs/mkpaper_ga.py`, `briefs/paper_ga_body.html`), its title and abstract the interval.

**Finding check on the revised title + abstract** (fresh context-free Opus, title and abstract alone in a file), verbatim:
VERDICT: the findings are confirmed novel and sound if the calculations hold

Strongest objection (referee): no filament-shaped injection; a thin filament's hull is mostly empty, so its prediction
would be diluted below +0.62, which moves the 1.3σ and the DR2 forecast but strengthens "cannot yet tell". Cautions: the
test is weak (even a perfect count at ~9 expected separates +0.62 from +0.15 by ~1.4σ); the DR2 forecast assumes the count
holds and no systematic in the expected count; the DESI/SDSS bias range is assumed, not measured; a varied hull definition
is missing. It suggested dropping the 1.3σ from the title; kept, since the title leads with "cannot yet be told" and the
1.3σ is the interval's own number. Closest work: 2201.06875, 2402.07591, 2604.17534, 2502.03515, 2505.11072, 2504.14940,
1211.6256, astro-ph/0111135, 1311.0900, 1104.0949.

**Finding check on the author email** (`outreach/authors-giant-arc.md`, rewritten for the revised paper; fresh context-free
Opus, the email alone), verbatim:
VERDICT: the findings are confirmed novel and sound if the calculations hold

Strongest objection (referee): only a uniform-in-hull Arc was injected; for a thin filament a real Arc's prediction falls
toward zero, so neither the 1.3σ nor the DR2 forecast applies to it, while "cannot yet tell" survives as a statement of
low power. Cautions: the prediction's sensitivity to the injected amplitudes and profiles; the DESI-only sightlines'
detection function; the 0.5σ against chance includes clustering scatter beyond Poisson. Closest work: 2201.06875,
2402.07591, 2604.17534, 2502.03515, 2505.11072, 2504.14940, 1306.1700, 2305.20016, 2412.15383. Email not changed after it;
nothing sent.
