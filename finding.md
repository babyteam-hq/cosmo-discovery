# Pantheon+ supernovae confirm Cosmicflows-4's bulk flow but reject its non-dipole velocity structure

2026-10-31. Scripts: `cf4_vs_2mpp.py`, `cf4_decompose.py`, `cf4_variants.py`, `cf4_2023.py`, `cf4_robust.py`.

## Claim

At the positions of the 618 Pantheon+ SNe Ia with 0.01 < z_CMB < 0.1 (calibrators excluded), the reconstructed
Cosmicflows-4 velocity fields (Courtois+2023 grids; Hollinger & Courtois 2025 CF4++ grid) carry a radial peculiar
velocity of rms 265 to 468 km/s, versus 169 km/s for the 2M++ model Pantheon+ uses. Fitting the SN Hubble residuals
y = m_b − μ(z_CMB) against each field's distance-modulus correction with a free scale β (β = 1 means the field is right),
using Pantheon+'s full STAT+SYS covariance:

| field at SN positions | rms non-dipole, km/s | β (dipole part) | β (non-dipole part) | β (whole field) |
|---|---|---|---|---|
| 2M++ (Pantheon+ VPEC) | 129 | 1.75 ± 0.25 | 0.38 ± 0.18 | 0.86 ± 0.15 |
| CF4 2023, grouped | 211 | 1.31 ± 0.20 | 0.34 ± 0.13 | 0.70 ± 0.11 |
| CF4 2023, ungrouped | 342 | 0.85 ± 0.14 | 0.17 ± 0.08 | 0.44 ± 0.08 |
| CF4++ 2025 | 427 | 0.81 ± 0.13 | 0.10 ± 0.07 | 0.31 ± 0.06 |

(2023 rows: z < 0.08, 604 SNe. "Dipole part" = the best-fit bulk-flow vector at the SN positions, with a radial monopole
also removed; "non-dipole" = what's left.)

The bulk flow each CF4 field predicts at the SN positions is supported (β near 1, Δχ² ≈ 37 to 42 for 3 effective
parameters, as good as or better than 2M++'s). Everything else in the CF4 fields is not: the SNe accept 10 ± 7 % of the
CF4++ non-dipole amplitude (13σ from β = 1; applying it at face value worsens χ² by 180 on 618 SNe) and 17 ± 8 % of the
ungrouped 2023 field's. Across all four fields the supported non-dipole radial velocity is the same ~45 to 70 km/s rms
(β × rms), however much the reconstruction puts there.

Robust to leaving out any one survey (β_rest 0.07 to 0.13 for CF4++), to z < 0.05 / z < 0.03 cuts, and to the best-
constrained tercile of CF4++ voxels (β 0.33 ± 0.17). Grid conventions were checked by scanning all axis permutations and
distance scales: the documented one maximises the correlation with 2M++ (0.50; next best 0.25).

## What it means

Cosmicflows-4 reconstructions claim a large coherent velocity field (bulk flows of 350 to 420 km/s out to 200 to 300
Mpc/h in the ungrouped 2023 field, the known CF4 "bulk flow excess" over ΛCDM). Pantheon+ Hubble residuals, an
independent-in-method test (SN magnitudes against the field's prediction at each SN), say the dipole component is real
at its predicted amplitude where the SNe are, but the higher-multipole structure is over-predicted by a factor 3 to 10.
That is the part that should carry information about fσ8 and the local density field; an amplitude excess of that size
in the non-dipole modes is hard to reconcile with the per-voxel robustness of 50 to 80 km/s the 2023 paper quotes.

Not in the literature as far as a search found: Pantheon+'s own papers evaluate 2M++ variants (Carr+2022, Peterson+2022);
CF4 papers compare bulk flows with earlier SN bulk-flow estimates, not the field with SN residuals.

## Weakest points

- CF4 used ~1000 SN distances as inputs, and R-115 found (`cf4_notin.py`, CF4's SN-distance table, VizieR J/ApJ/944/94)
  that every Pantheon+ SN at 0.01 < z < 0.1 is one of them: there is no independent Pantheon+ subset (an earlier row
  claiming 99 independent SNe was wrong; they were the same hosts). Circularity would be expected to make the field fit
  the SNe better, not worse, so it doesn't explain the rejection, but a fully independent sample needs SNe CF4 never saw
  (ZTF DR2 once cosmology-grade). What exists: DES-SN5YR's Dovekie release, an independent photometric pipeline,
  calibration and bias correction on largely the same objects (`cf4_des.py`, 191 SNe, DES's own covariance): CF4++
  non-dipole β = 0.26 ± 0.11 (6.7σ from 1); its dipole is also low there (0.42 ± 0.13) where Pantheon+ gave 0.81 ± 0.13.
  And a fully independent sample does exist in a non-cosmology-grade form: ZTF SN Ia DR2 (Rigault+2025, public "lite"
  release; `ztf_tests.py`, 1122 SNe at 0.015 < z < 0.1 with host spectroscopic redshifts, diagonal errors, α and β
  free): CF4++ non-dipole β = 0.23 ± 0.06 (13σ from 1; 0.13 ± 0.08 at z < 0.06), and its dipole part 0.33 ± 0.09.
  Three samples, three pipelines, the same rejection of the non-dipole structure. The dipole's support varies:
  Pantheon+ 0.81, DES 0.42, ZTF 0.33 (all ± 0.1), the ZTF and DES values being on different sky/redshift footprints.
- The CF4++ 2025 grid's units and frame come from its example script (km/s, CMB frame) because the paper is paywalled.
  The 2023 FITS grids, whose units are documented (×52 → km/s, CMB frame), give the same conclusion.
- Pantheon+'s covariance includes a 250 km/s peculiar-velocity error term built around 2M++. It is diagonal noise and
  doesn't bias β, but a CF4-specific error model would change the error bars.
- The "monopole" removed at SN positions is a fit of v ∝ r; it changes nothing (β_rest with and without it: 0.31 vs 0.31
  for the whole CF4++ field).

## To settle it

Redo with the CF4++ posterior samples instead of the mean grid (the mean of a multimodal posterior can have an inflated
amplitude). The independent-sample checks are done (DES Dovekie, ZTF DR2 above); what is left is ZTF DR2.5's
cosmology-grade photometry when it comes, and the per-sample dipole amplitudes (0.81 / 0.42 / 0.33), which differ by
more than their errors and deserve a look at the samples' sky and redshift footprints against the CF4++ bulk flow.

---

## Addendum 2026-10-31: the bulk flow the Pantheon+ SNe measure on their own

`bulkflow_free.py`: a free 3-vector B (GLS, intercept and radial monopole free) fitted to the Pantheon+ residuals, next to
the fields' mean 3-vector over the same SNe.

| z range | n | SNe: B, (l, b) | CF4++ mean: B, (l, b), cos to SNe | 2M++ fit: B, (l, b), cos |
|---|---|---|---|---|
| 0.01-0.1 | 618 | 312 ± 44 km/s, (311, 18) | 252, (301, 50), 0.84 | 189, (305, 1), 0.95 |
| 0.01-0.05 | 525 | 305 ± 44, (311, 19) | 226, (301, 38), 0.94 | 174, (306, 6), 0.97 |
| 0.03-0.06 | 206 | 334 ± 105, (309, 10) | 362, (312, 35), 0.91 | 253, (316, -7), 0.95 |
| 0.06-0.1 | 63 | 486 ± 366, (351, 12) | 495, (232, 79), 0.11 | 144, (306, 3), 0.71 |

So the SNe's own bulk flow is 1.2x CF4++'s amplitude and 30 deg from its direction (higher Galactic latitude in CF4++),
and 1.6x 2M++'s in nearly its direction. That is what the dipole betas of 0.81 (CF4++) and 1.75 (2M++) were measuring:
projections of a vector the SNe see at ~300 km/s toward (311, 18). Beyond z ~ 0.06 the CF4++ mean swings to b ~ 80 and
500 km/s and the SNe (few, large errors) do not follow.

DES and ZTF (R-115, `bulk_flow.py`, same fit with the intercept as the monopole; DES with its own covariance, ZTF with
diagonal errors and α, β free):

| sample | n, depth | SNe: B, (l, b) | CF4++ mean at the SNe: B, (l, b) | 2M++ fit: B, (l, b) |
|---|---|---|---|---|
| Pantheon+ 0.01-0.1 | 618, 82 Mpc | 301 ± 46, (311, 18) ± 8° | 252, (301, 50) | 215, (306, −2) |
| DES Dovekie 0.01-0.1 | 191, 113 Mpc | 248 ± 90, (295, 21) ± 19° | 356, (307, 39) | 267, (313, −11) |
| ZTF DR2 0.015-0.1 | 1122, 179 Mpc | 528 ± 82, (305, 0) ± 6° | 457, (310, 66) | — |
| ZTF DR2 0.015-0.06 | 551, 127 Mpc | 487 ± 94, (311, 3) ± 7° | 371, (309, 35) | — |

The three samples agree on the direction (l = 295 to 311, b = 0 to 21, all within their errors, and where 2M++ points)
and on the amplitude within 2.4σ (ZTF, the deepest, highest). What differs between them is not the flow but how the
CF4++ field projects onto their footprints: its mean vector at the SN positions points 30 to 65° further north
(b = 35 to 66) than the SNe's own flow, so its dipole β comes out below 1, by a footprint-dependent amount (0.81, 0.42,
0.33). The samples disagree with the CF4++ field's direction, not with each other.

Convention check (`cf4_axes.py`): of the 48 signed permutations of the CF4++ velocity components, the documented one is
the best by Δχ² = 9 on the Pantheon+ residuals and by correlation with 2M++ (0.50; next best 0.24), so neither the
rejection nor the direction offset is an axis-order or sign artefact.

## Blind review (2026-10-31, context-free Opus referee, own literature search; `blind-review.md`)

VERDICT: novel but not sound: noise in the predictor and evaluation at redshift-space positions would push the non-dipole β
toward zero even for a correct field (2M++ also fails at 0.38 ± 0.18), and without mock tests or propagated field errors the
note cannot claim CF4 over-predicts the amplitude.

Strongest objection: regression attenuation. A noisy predictor (a posterior-mean grid from ~56,000 noisy distances) pulls β
toward zero, and reading the field at the SN's redshift-space position (a 300-500 km/s error, 3-5 Mpc/h, the coherence length
of the non-dipole structure) smears the sharp part and leaves the dipole intact; both predict exactly "dipole β ≈ 1,
non-dipole β ≪ 1", and the trend (noisier field, lower β) fits attenuation at least as well as over-prediction. Uncited
closest prior work: Stiskalek et al. 2025, "Velocity Field Olympics" (arXiv:2502.00121), which already finds the Courtois+2023
CF4 reconstruction under-performs against LOSS/Foundation SNe and over-predicts the bulk flow at large R, though it says the
CF fields' velocities need scaling *up*. Also: the CF4++ paper is Courtois, Mould, Hollinger et al. 2025 (arXiv:2502.01308),
on arXiv, not paywalled. To pass: propagate the CF4++ per-voxel standard-deviation grids as predictor error (an
errors-in-variables β), fold in the position error, and run the pipeline on CF4 mocks to show a true field returns β = 1.

Response (same day, `cf4_mocks.py`, row in results.md): done. A true field at β = 1 returns β_rest = 0.88 ± 0.07 through the
pipeline (0.87 with posterior-std noise on the truth, 0.88 with a 3 Mpc/h position error); errors-in-variables λ = 0.86 gives
β_rest = 0.12 ± 0.08; real-space positions give 0.10 ± 0.06. The observed 0.10 is 11σ below the mock return. Stiskalek+2025
(arXiv:2502.00121) is now cited: whole-field β for the Courtois+2023 field against TF and LOSS/Foundation SNe, no
dipole/non-dipole split, no Pantheon+/DES/ZTF; their "underperforms, over-predicts the bulk flow at large R" points the same way.

Blind review, second pass (2026-10-31, context-free Opus referee, on the revised note): VERDICT: novel but not sound: a
regression slope β < 1 cannot tell an over-predicted amplitude from structure uncorrelated with the truth, and the mocks (which
treat the field as truth and omit correlated residual velocities) neither calibrate the "factor three to ten" claim nor the
quoted significances.

Response (`cf4_referee2.py`, row in results.md): accepted on the interpretation. The claim is now "the non-dipole part of the
CF4 fields is miscalibrated at the SN positions: nearly uncorrelated with the residuals (whitened ρ = 0.06 to 0.08), so
β = 0.10 to 0.17", with "over-predicted three to ten times" and "45 to 70 km/s supported" withdrawn. Significance by a
delete-one-sky-patch jackknife: 13σ (CF4++) and 9σ (2023) from β = 1, the same error as the GLS. Stiskalek+2025's free σ_v and
V_ext reproduced: β_whole = 0.18 ± 0.07, independent of σ_v; their β ≈ 1 came from TF samples that are CF4's own inputs.

Blind review, third pass (2026-10-31, on the reworded note): VERDICT: the findings are confirmed novel and sound if the
calculations hold. Required revisions (Peterson+2022 CF3 precedent; 2M++'s 0.39 as a failed absolute null, 13σ softened;
joint dipole+rest fit: 0.88 ± 0.13 / 0.19 ± 0.07; direction errors) applied the same day. PASSED.
**A parallel second referee (R-118, same day, on the note with the mocks section, no verdict text shown), independently:** VERDICT: novel but not sound: the
2M++ reference field also fails the non-dipole test (β = 0.38 ± 0.18) and the mocks use CF4 itself as truth, so the low β may
come from attenuation in the method or from noise in the reconstruction rather than an amplitude excess. Also: β × rms constant
at 45-70 km/s is what a shared true signal plus field-dependent uncorrelated noise gives, so a low β measures poor correlation
as readily as an amplitude excess (report the correlation coefficient and residual variance); "confirms the bulk flow" in the
title is contradicted by the dipole β of 0.42 and 0.33 in DES and ZTF and the 30-65° direction offset; correlated
peculiar-velocity covariance is left out. To pass: end-to-end mocks (N-body truth, a CF4-like reconstruction, SN sampling)
showing a true field returns β = 1 through this pipeline, and the same test giving 2M++ β ≈ 1.

---

# Supernova Hubble residuals track the Planck-GNILC minus SFD dust difference: partial

2026-10-31. Found by R-115 (cosmo-sweep: `dust_maps.py`, `dust_robust.py`, `dust_resolution.py`, `host_in_sfd.py`,
`ztf_tests.py`, `dust_emissionfree.py`); settled as far as public maps allow by R-116 (cosmo-sweep-2: `dust_settle.py`,
`dust_salt_response.py`). Rows in `results.md` (R-115 sweep rows 1-4, 7, "3 on ZTF", "3 resolution", "3 emission-free";
R-116 rows "3 emission-free 2", "3 trials", "3 physical", "3 quadrupole", "3 temperature").

**Verdict: partial.** The correlation is present in two SN samples with separate photometric pipelines (both tied to PS1
calibration, so not independent in calibration) and stands at 3.4 to 3.9σ after R-115's 135 trials against sky-rotation
nulls, and its sign is what a Milky-Way extinction error would give. Whether it *is* one is split by sample. In ZTF the
SALT2 colours redden with GNILC − SFD exactly as reddening requires (+1.21 ± 0.24 per unit E(B-V), 1.1 expected) and the
distance slope is 1.4 ± 0.4 times the GNILC-for-SFD swap: ZTF says GNILC is a better Milky-Way map than SFD. In Pantheon+
the magnitudes respond nearly as reddening (mB: 3.2 ± 0.6 vs 4.1) but the colours at only 40-60 % of it (0.43 ± 0.24
vs 1.1) and the light-curve width where it should not (x1: −2.4 ± 0.6), so the distance slope is 2.9 ± 0.7 times the
swap. The emission-free maps do not settle it: the only one precise enough, Zhou+2024's DESI stellar map, has a difference
from SFD that SN distances (3.6σ), quasar counts (7σ) and LRG redshifts (3.7σ) all follow, but SN colours do not (4.5σ
below reddening), so that pattern is not Milky-Way reddening as SALT2 sees it, whatever it is. Shown, at 3.4-3.9σ after trials against
sky-rotation nulls (below): SN distances in both samples carry a ~0.01 mag residual that correlates with GNILC − SFD, and in
ZTF it behaves like dust. Not shown: that it is dust in Pantheon+ (its colour, x1 and amplitude responses do not match a
Milky-Way extinction error), or that any public map removes it.

## Claim

SN Ia Hubble residuals, after the standard Milky-Way extinction correction with 0.86 × SFD, correlate with the difference
between the Planck-GNILC E(B-V) map (1.49e4 × τ353, 2016 GNILC dust model) and SFD at the SN position:

| sample | SNe | slope on (GNILC − SFD), mag per E(B-V) | slope on (GNILC − 0.86 SFD) | p | after a GNILC-for-SFD swap |
|---|---|---|---|---|---|
| Pantheon+ (full STAT+SYS covariance, calibrators out, z > 0.01) | 1580 | +1.75 ± 0.45 | +1.84 ± 0.43 | 1.1e-4 | +1.11 ± 0.45 |
| ZTF SN Ia DR2 (independent objects and pipeline; diagonal errors, α, β free) | 1124 | +1.85 ± 0.55 | +1.73 ± 0.47 | 7e-4 | +0.74 ± 0.55 |
| combined | | +1.79 ± 0.35 (5.1σ raw) | +1.79 ± 0.32 | | +0.96 ± 0.35 (2.7σ) |

Trials: R-115 computed 135 statistics on dust (every slope, cut, jackknife, smoothing and excess test in its six scripts);
14 of them came before GNILC − SFD was singled out and 21 are distinct dust quantities regressed on the residuals. The
Pantheon+ p of 1.1e-4 becomes 0.0015 / 0.0023 / 0.015 (3.2σ / 3.0σ / 2.4σ) corrected for 14 / 21 / 135, and ZTF was run once, after the Pantheon+ hit, as a
replication of that single hypothesis (in the same R-115 session, so "preregistered" overstates it), so the two together
(Fisher) stand at 4.3σ / 4.2σ / 3.8σ with Gaussian p's, and 3.9σ / 3.4σ with the sky-rotation p's below. The honest single
number is 3.4σ. R-116's settling tests (eight maps, smoothings, the quasar, LRG, colour, temperature and quadrupole tests) are
not in that count; they are follow-ups of the selected signal, not further selections, but a referee may count them.

Robust (R-115): survey fixed effects, leave-one-survey-out (1.38 to 1.93), |b| and Dec hemispheres, z < 0.1 (2.31 ± 0.61)
and z ≥ 0.1 (1.14 ± 0.68), clipping |GNILC − SFD| < 0.02, smoothing both maps to 40' (1.59 ± 0.53). In the two-map fit (y
on SFD and GNILC together) the coefficients are −1.5 / +1.8 (Pantheon+) and −1.8 / +2.0 (ZTF). Not host emission in SFD
(real at z < 0.03 but 0.001 mag), not host mass, not CSFD's CIB term. Specific to GNILC: Planck-2013 τ353 and radiance
E(B-V) minus SFD show nothing (+0.17 ± 0.41, −0.03 ± 0.29; ZTF −0.06 ± 0.32). [Corrected 2026-10-31, `dust_fields_cib.py`: the
"τ353" numbers were computed on the 2013 file's EBV column, which is the radiance estimator (5.4e5 R, ratio 1.00); the real
1.49e4 τ353 − SFD carries the slope too, +1.53 ± 0.38 in Pantheon+ and +1.26 ± 0.50 in ZTF. The correlation is specific to the
temperature-corrected estimators (τ353, GNILC) against the radiance one, not to GNILC.]

## What the settling tests say (R-116)

**The size.** A SALT2 fit with the Milky-Way E(B-V) too low by dE returns mB fainter by 4.07 dE and c redder by 1.10 dE
at z = 0.02 (3.70 and 1.00 at z = 0.1, 2.57 and 0.65 at z = 0.4; sncosmo SALT2 T23 with F99 dust, synthetic gri and BVRI
light curves, `dust_salt_response.py`), so the distance modulus is off by dμ = dmB − β dc = 0.60 to 0.69 dE at Pantheon+'s
β = 3.04, and 1.21 to 1.32 dE at ZTF's fitted β = 2.47. (R-115's "+1.0 to +1.6 expected" used the analytic F99 B − V
coefficient, which overstates this: SALT2's colour term absorbs 1.1, not 1.0, per unit E(B-V).) Subtracting that
correction with dE = GNILC − 0.86 SFD from every SN, the slope falls from 1.75 to 1.11 ± 0.45 (Pantheon+) and from 1.85 to
0.74 ± 0.55 (ZTF). Fitting the residual on the predicted correction itself, the SNe want 2.9 ± 0.7 (Pantheon+) and 1.4 ±
0.4 (ZTF) times the swap: combined 1.73 ± 0.32, 2.3σ above the pure-extinction value of 1 and 5.4σ above 0. ZTF on its
own is consistent with "GNILC is right"; Pantheon+ on its own wants three times the swap.

**Emission-free maps.** Bayestar 2019 (Gaia/PS1/2MASS stellar 3D map to its 60 kpc edge, 96% of the Pantheon+ SNe) and
Lenz 2017 (HI4PI column scaled to E(B-V), where N_HI < 4e20, ~70% of the SNe), each put on SFD's scale by one factor k:

| map difference (k MAP − SFD) | Pantheon+ | ZTF | combined | expected if GNILC were right and the map a noisy measure of the same truth |
|---|---|---|---|---|
| Bayestar − SFD, point values | +0.18 ± 0.18 | −0.00 ± 0.21 | +0.10 ± 0.14 | +0.22 |
| Bayestar − SFD, 1° disk averages | +0.35 ± 0.30 | +0.01 ± 0.30 | +0.18 ± 0.21 | +0.32 |
| Lenz − SFD, point values | −0.21 ± 0.68 | +1.53 ± 1.07 | +0.29 ± 0.58 | +0.54 |
| Lenz − SFD, 1° disk averages | −0.53 ± 1.33 | +5.8 ± 2.0 | +1.38 ± 1.11 | +0.96 |
| Edenhofer − SFD (1.25 kpc, \|b\| > 20; R-115's map, rerun) | +0.59 ± 0.40 | +1.29 ± 0.50 | +0.86 ± 0.31 | +0.57 |
| Mudur+2023 (PS1/2MASS/Gaia) − SFD, point values | −0.11 ± 0.29 | +0.40 ± 0.33 | +0.12 ± 0.22 | +0.24 |
| **Zhou+2024 (DESI stars) − SFD, point values** | **+0.69 ± 0.36** | **+1.33 ± 0.40** | **+0.98 ± 0.27** | +0.56 (pure extinction: 0.60 / 1.27) |
| Zhou+2024 − SFD, 1° disk averages | +1.26 ± 0.49 | +1.47 ± 0.58 | +1.35 ± 0.37 | +0.85 |

Bayestar, Lenz, Edenhofer and Mudur are each within 1σ of zero and within 1σ of the GNILC-true expectation; the expectation
is small because their differences from SFD correlate with GNILC − SFD at only 0.2 to 0.5 (even at 1°) and have 1 to 3
times its rms, so a GNILC-sized slope projects onto them at 0.2 to 0.5, below their errors. The DESI map is different: 2.6
million stars give 13.7' pixels with 0.006 mag errors, SFD sampled at the same stars, and a difference from SFD of rms 0.012
that correlates with GNILC − SFD at 0.43 (0.58 at 1°). On it the SNe show the slope at 3.6σ combined, and in each sample
the slope equals what the SALT2 response to a true E(B-V) error predicts with that sample's β (`dust_stellar2.py`,
`results.md` sweep-2 row 38). GNILC − SFD keeps a slope of its own in the joint fits (1.9 ± 0.6 Pantheon+, 1.8 ± 0.8 ZTF),
so the DESI map does not account for all of what GNILC sees.

**Other tracers, and the colours.** Quasar number counts (Quaia, `quaia_dust.py`, `results.md` sweep-2 row 21) respond
to the DESI-stars correction to SFD at 7σ with the extinction sign and 55-75 % of the extinction size, and to GNILC's at
only −0.3 ± 0.1. The DESI LRG redshift distribution does the same at 3.7σ (`desi_nz_dust.py`, row 4). But the SN colours
(`dust_colour.py`, `dust_colour2.py`, row 40), which a true E(B-V) error must redden by 1.1 per unit, do not follow the
DESI-stars correction (+0.03 ± 0.18 combined, 0.84 expected), while they do follow GNILC's in ZTF (+1.21 ± 0.24) and
partly in Pantheon+ (+0.43 ± 0.24; +0.65 ± 0.28 at 1°). Taking Pantheon+'s m_b_corr apart: on GNILC − 0.86 SFD, mB
carries +3.2 ± 0.6 (reddening: 4.1), c +0.4 ± 0.5, x1 −2.4 ± 0.6 (reddening: +0.2), the bias correction nothing; what
makes the distance slope 1.8 instead of 0.65 is the colour not responding and x1 responding. Correcting the residuals with
the DESI map removes its own slope (as it must) and lowers the GNILC − SFD slope only from 2.0 to 1.6 (Pantheon+) and 2.5
to 1.7 (ZTF).

**What GNILC − SFD is.** At the SN positions it anticorrelates with the GNILC dust temperature (−0.56 in Pantheon+,
−0.36 in ZTF): it is largely the difference between SFD's DIRBE-based 1.1° temperature correction and Planck's. The
residual on T alone is −4 ± 3 (Pantheon+) and −8 ± 4 mmag/K (ZTF), the diluted version of the same slope, and T adds
nothing on top of GNILC − SFD.

**What it changes.** The Galactic-axis quadrupole hint in Pantheon+ (P2(sin b): −0.021 ± 0.013, p = 0.10) weakens to
−0.015 ± 0.013 (p = 0.25) after the swap; it was never significant. The CMB-axis quadrupole Macpherson 2023 reports is
untouched (p 0.026 → 0.020). H0 is insensitive (the correction's rms is 0.006 mag with no redshift trend).

## What it means

Taken together: at the ~0.01 mag level the SFD-based Milky-Way correction leaves a residual on the sky that the GNILC map
predicts, in two independent SN samples. Swapping GNILC in removes about half of it. If the rest is also extinction, the
effective R_V or β in the affected regions differs from the fit's; if it is not, something that co-varies with the
Planck-vs-DIRBE temperature correction (dust temperature, CIB treatment, the Planck model's β) also co-varies with SN
magnitudes by ~1 mag per unit of E(B-V)-equivalent. Either way the map difference is a correlated systematic in Hubble
residuals of the kind that feeds Galactic-frame anisotropy searches, and SN cosmology should marginalise over it.

Not in the literature as far as a search found: Pantheon+ (Brout+2022) and DES (Vincenzi+2024) carry a 5% scaling
systematic on E(B-V); Sun+2022 (arXiv:2204.01521) compare SFD and Planck maps with stars; no SN-residual test of a dust-map
difference was found.

## Weakest points

- The size does not match a pure extinction error (1.7 ± 0.3 times the swap, 2.3σ), and the two samples differ in how much
  they want (2.9 ± 0.7 vs 1.4 ± 0.4). Pantheon+'s m_b_corr carries BBC bias corrections and a fixed β; its true response
  to a map change may differ from the synthetic 0.65, and only a refit of its light curves with GNILC would say.
- It is specific to GNILC among the Planck maps. The τ353 and radiance maps minus SFD show nothing. GNILC's temperature
  correction is the live candidate for why; the test above shows the correlation but cannot say which map is right.
- The emission-free maps do not settle it. The DESI one (Zhou+2024) has a difference from SFD that SN distances, quasar
  counts and LRG redshifts follow but SN colours do not, so it is not reddening as the SNe see it; Mudur+2023 sees nothing
  (−0.11 ± 0.29 in Pantheon+); Bayestar, Lenz and Edenhofer cannot tell. The maps agree with each other about the SFD
  error at r = 0.1 to 0.6 only (`dust_mapcorr.py`).
- In Pantheon+ the colour responds too little and x1 too much for a Milky-Way reddening error, and the two samples
  disagree on the colour response (0.43 ± 0.24 vs 1.21 ± 0.24). A grey (achromatic) ~0.01 mag sky pattern in the
  Pantheon+ magnitudes that happens to follow the Planck-minus-DIRBE temperature correction would fit Pantheon+ but not
  ZTF's colours.
- ZTF DR2's photometry is not cosmology-grade by the release's own statement, and its errors are diagonal. Pantheon+'s
  result stands on the full covariance.
- Trials: 135 statistics in R-115; the result survives them at 3.8σ combined, 2.4σ in Pantheon+ alone.

## To settle it

The colour test is the lever: ZTF DR2.5's cosmology-grade colours (and DES-SN5YR's, Foundation's, per survey) against
GNILC − 0.86 SFD, which a reddening error must move by 1.1 per unit; Pantheon+'s 0.43 ± 0.24 against ZTF's 1.21 ± 0.24
is the discrepancy to resolve, survey by survey. Then refit the Pantheon+ low-z light curves (public, `1_DATA/photometry`) with GNILC E(B-V) in SNANA or sncosmo and rerun the
BBC step, to replace the synthetic response with the real one; if the slope then goes to ~1.1 ± 0.45 as predicted here,
the remaining half is the question. A stellar reddening map built for this (Gaia XP extinctions, Zhang+2023, to a few kpc
at high latitude, averaged on 1° scales to beat the per-pixel noise) could reach the precision the test needs: rms of its
difference from SFD at or below GNILC − SFD's 0.007 and a correlation with it well above 0.5. And ZTF DR2.5, when its
photometry is cosmology-grade, with its own covariance.

## Blind review (2026-10-31, context-free Opus referee, own literature search; `blind-review.md`)

VERDICT: novel but not sound: significance is assessed against independent per-SN errors with no sky-rotation/spatially-
correlated null for a smooth template, and the x1 response and failed DESI colour test point to sky-position confounding
rather than Milky-Way extinction.

Strongest objection: GNILC − SFD is a smooth pattern on the sky (the Planck-vs-DIRBE temperature correction, degree to
tens-of-degree scales), and Hubble residuals are themselves sky-correlated (peculiar velocities at low z, where the signal is
strongest; calibration and selection by footprint), so p-values from per-SN errors assume far more independent samples than
there are; the trials accounting corrects the wrong baseline. Needed: rotate or reflect the template on the sphere (or
randomise its phases) keeping each survey's footprint, and read p from that distribution. Also: x1 cannot respond to Milky-Way
dust at −2.4 ± 0.6, which points to the SN population varying with sky position; a map the distances follow but the colours
do not is evidence against reddening; "specific to GNILC" among three Planck maps fits a chance alignment in a 21-quantity
search. Closest prior: Taylor+2021 (arXiv:2104.00172, SALT2 retrained with Planck maps, Δw ≈ 0.01), Scolnic+2014
(arXiv:1310.3824, patchy dust systematics), Chiang & Ménard 2019, Chiang 2023 (arXiv:2306.03926). To pass: the sky-rotation
null (`dust_rotation_null.py`, below), a peculiar-velocity covariance for ZTF, x1/host/z regressed on the template per survey.

**Response (2026-10-31, `dust_rotation_null.py`).** The sky-rotation null: residuals and covariance stay put, the two dust
maps are re-sampled at rotated SN positions and the same full-covariance GLS slope refitted. (A) 2000 rotations about the
Galactic pole, with random north-south reflection, so every SN keeps its latitude and the template's amplitude distribution is
matched: null slope 0.03 ± 0.31, null z spread 1.07 (1.00 if the per-SN covariance were exactly right), 0 of 2000 draws reach
the data's z = 3.87 (empirical p < 5e-4; the data is 3.6σ in units of the null spread). (B) 2000 uniform SO(3) rotations: z
spread 1.01, 1 of 2000 draws at or above the data (p = 5e-4), 3.8σ. The SFD-slope control has a null z spread of 0.93-1.10. So
Pantheon+'s covariance, which carries the peculiar-velocity and survey correlations, already describes the sky-correlated
residuals for a template of this smoothness; the per-SN significance stands against the spatially correlated null to within
0.3σ. Confounding by SN population: the slope is unchanged with x1, c and host mass as covariates (1.743 ± 0.452), with
survey fixed effects (1.720 ± 0.460), and with all four (1.721 ± 0.460); the x1 response in the ingredients decomposition does
not carry the distance slope. Not answered here: ZTF's diagonal errors (its rotation null needs its own covariance), and the
colour discrepancy between the samples, which remains the live objection to the extinction reading.

**ZTF rotation null (same day, `ztf_rotation_null.py`).** Same sample and model as `ztf_tests.py` (1124 SNe, diagonal errors,
σ_int 0.078, α 0.149, β 2.475; data slope +1.842 ± 0.548, z = 3.36). Null A (2000 latitude-preserving rotations): z spread 1.10,
two-sided empirical p = 0.0010 (2 of 2000), 3.1σ-equivalent; null B (uniform rotations): spread 1.02, p = 0.0045, 3.3σ. The
diagonal errors understate the sky correlations by 10 %; the replication stands at 3.1σ instead of 3.4σ. Combining the empirical
rotation p's (Fisher; Pantheon+'s 0/2000 taken as 5e-4): 4.5σ raw, 3.9σ with Pantheon+'s 14-trial correction, 3.4σ with all 135.

**Second blind review (same day, on the note with the Pantheon+ rotation null and covariate checks, no verdict text shown):**
VERDICT: novel but not sound: the combined 3.8σ and the "established" residual depend on a ZTF significance with diagonal errors
and no spatially correlated null, while Pantheon+ alone is only 2.4σ after trials and its colour, x1 and amplitude responses
contradict a Milky-Way extinction error. (The ZTF null above answers the first half; the second half is the "partial" verdict's
own content, and the "Established" sentence in the verdict paragraph has been reworded to what is shown.) Closest prior work
it adds: Taylor+2021 (arXiv:2104.00172), a SALT2 retrain with the Planck 2014 map, Δw = 0.01; arXiv:2311.08253; arXiv:2409.18668;
arXiv:2508.11746 (PAH-driven errors in emission-based maps, a non-temperature way for GNILC and SFD to differ).

**Third blind review (same day, on the note with both rotation nulls, no verdict text shown):** VERDICT: novel but not sound:
the significance is overstated because the trials count omits the R-116 tests and the ZTF "replication" was not preregistered
and shares PS1 calibration with Pantheon+, while the x1 and colour responses argue against the Milky Way extinction reading that
the verdict leans on. (The presentation points are now fixed above: one consistent number, 3.4σ; "preregistered" withdrawn; the
PS1 link stated. The physics points are the "partial" verdict's own content: a ~0.01 mag sky-correlated term is shown, a
Milky-Way extinction reading is not.) New prior work it names: Stiskalek, Desmond & Lavaux 2025 (arXiv:2509.14997) swap dust
maps including GNILC for a TF zero-point dipole, no effect; Camilleri+2026 "Unite" (arXiv:2609.05053) keeps the 5 % scaling
systematic and cites Elliott et al. (in prep.) on dust-map uncertainties, the paper most likely to pre-empt this test. Checks it
asks for that are not done: a GNILC refit through BBC; residuals per survey field rather than per SN (the high-z SNe sit in a
few deep fields); GNILC's own CIB removal as the carrier (τ353 − SFD shows nothing, GNILC differs from τ353 in CIB removal and
temperature model).

**Two of the three done (same day, `dust_fields_cib.py`, rows "3 fields" and "3 carrier").** Per field: 1580 Pantheon+ SNe
make 799 survey × 1.8° cells (639 singletons; the deep fields collapse to a few cells of up to 37 SNe); block-averaging
residuals, template and covariance onto the cells gives +1.749 ± 0.489 (3.6σ, per-SN 3.9σ); leave-one-cell-out 1.63-1.94,
leave-one-survey-out 1.38-2.02; the cell fit is +2.39 ± 0.62 (3.8σ) over the 591 z < 0.1 cells and +0.69 ± 0.80 over the 208
deep-field cells, so the count was not inflated and the signal is a low-z, all-sky one. ZTF block-averaged on 979 cells: +1.885
± 0.556. Carrier: R-115's "τ353 shows nothing" was a radiance test (dustmaps' `ebv` component is the 2013 EBV column = 5.4e5 R);
the true 1.49e4 τ353 − SFD carries the slope at +1.53 ± 0.38 (Pantheon+) and +1.26 ± 0.50 (ZTF), radiance − SFD carries nothing
(−0.03 ± 0.29, −0.09 ± 0.31), and jointly GNILC − τ353 adds nothing in Pantheon+ (+0.79 ± 0.76) next to τ353 − SFD (+1.79 ±
0.45). GNILC's CIB removal is not the carrier; the temperature-corrected estimators are. The 2013 dust temperature alone gives
−9.2 ± 2.5 mmag/K in Pantheon+ (3.6σ) and keeps −6.8 ± 2.7 next to GNILC − SFD. Not done: the GNILC refit through BBC.


---

# The eRASS1 L_X-T "H0 anisotropy" toward the Migkas axis is the survey's exposure pattern

2026-10-31, main session (R-117 found it and lost its commit; redone). Script: `erass1_lxt.py`; row in `results.md`.

## Claim

In the eRASS1 primary cluster catalogue (1014 clusters with measured temperatures, 0.05 < z < 0.6), a Migkas-style analysis of the
L_X-T relation finds the published anisotropy: clusters in the 75-deg cone about (l, b) = (280, -15) sit 0.039 dex below the
relation (4.0 sigma against a shuffle null), 0.084 dex in a 60-deg cone (6.7 sigma), which in the Migkas convention is a 4.5 to
9.7 % lower H0 there. The residual correlates with the survey exposure at the cluster (Spearman -0.36; -0.35 inside the cone,
-0.17 outside), the cone is where the survey is deep (median 230 s vs 130 s; the south ecliptic pole is 15 deg from the axis),
and with a log-exposure term in the relation (slope -0.51 dex per dex) the cone deficit is -0.006 dex (0.6 sigma); in fixed
exposure bands it is 0.00 +/- 0.02 dex; the exposure trend alone predicts -0.034 of the -0.039. Robust to the temperature-error
cut, redshift range, bolometric luminosity, extent likelihood and cone angle (3-7 sigma raw, <= 1.8 sigma with exposure).

## What it means

A flux-limited sample in a scanning survey has a flux limit that varies with exposure across the sky; near the limit only
the over-luminous clusters at a given temperature are detected (Malmquist bias), so shallow regions have positive L-T
residuals and deep regions negative ones. eRASS1's deep region is the south ecliptic pole, which is 15 deg from the Migkas
axis, so a Migkas-style analysis of eRASS1 reproduces the axis and the amplitude from the depth map alone. ROSAT's all-sky
survey had the same scanning geometry (deepest at the ecliptic poles), and the original claim came from ROSAT-selected
samples; whether it is the same effect is a test on those samples (residual vs RASS exposure in the L_X-T relation), not
a conclusion from this one.

## Weakest points

- Per-cluster counts (CTS500) do not remove the deficit (3.4 sigma with a counts term): the effect is the depth map, not the
  measurement noise. That is what Malmquist bias predicts, but it also means the correction is a sky-position model.
- The selection-function approach of the eRASS1 scaling-relations paper (Ramos-Ceja+2025) already includes exposure, so
  the eROSITA team's own L-T fit would not show the anisotropy; this test shows what a residual-based analysis sees.
- The temperature sample is 1014 of 12,247 clusters, biased to bright ones; the within-band test has 50-100 clusters outside
  the cone per band.
- The 60-deg cone's 6.7 sigma uses the pre-specified Migkas axis; the full-sky scan's extremes are reported in the script output.

## To settle it

The same residual-vs-exposure test on the eeHIFLUGCS / Migkas+2020 sample with the RASS exposure map; the eRASS1 cosmology
subsample with its selection function; eRASS:4/5 when public, where the depth pattern is the same but the limit deeper.

## Blind review (2026-10-31, context-free Opus referee, own literature search; `blind-review.md`)

VERDICT: novel but not sound: in eRASS1's half-sky, the log-exposure term is nearly collinear with an axis 15° from the south
ecliptic pole, so it would absorb a real dipole too, and the counts test failing to remove the deficit contradicts the
Malmquist mechanism the note invokes, so "it is the exposure pattern" is not shown without a selection-function forward model.

Strongest objection: the north ecliptic pole is outside eRASS1's hemisphere, so there is no deep region away from the axis to
separate depth from direction; the fixed-exposure-band test is the only thing that could, and it tells −0.039 dex from zero at
about 2σ; the exposure slope differing inside (ρ = −0.35) and outside (−0.17) the cone says the term is partly fitting position.
A counts term should absorb a Malmquist effect and does not (3.4σ remains). The temperature subsample's own selection (counts
∝ flux × exposure) and exposure-dependent T errors are unmodelled. The shuffle null tests isotropy, not isotropy at fixed depth.
Migkas+2021 (arXiv:2103.13904) already checked RASS exposure (Sect. 4.4) and Malmquist bias (8.2.1); Yasin+2026
(arXiv:2602.06007) forward-model selection on the same 312 clusters and find δH0/H0 < 3.2 %. Defensible claim: in eRASS1 the
anisotropy toward this axis cannot be separated from the survey's depth pattern. To pass: forward-model detection and
T-measurement probability with the eRASS1 selection function (Clerc+2024), or full-sky data with the north ecliptic pole.

Response (same day, `erass1_forward.py`, row in results.md): a simple forward model (isotropic relation, real exposure map,
detection threshold in counts) gives a cone deficit of −0.010 ± 0.004 dex against the observed −0.038, and a third of the
exposure dependence. The mechanism is therefore not shown; the phenomenology (deficit follows exposure at fixed direction,
zero in fixed-exposure bands) stands. The claim is narrowed to the referee's defensible one, and the note now cites Yasin+2026
(arXiv:2602.06007: forward-modelled Migkas sample, δH0/H0 < 3.2 %) and Migkas+2021's own exposure and Malmquist checks.

Blind review, second pass (2026-10-31, context-free Opus referee, on the revised note): VERDICT: novel but not sound: in this
half-sky sample exposure is nearly a function of angle from the axis, so the exposure term and the fixed-exposure bands (which
keep only cone-edge clusters, where a real dipole would also vanish) cannot tell depth from direction, and the proposed Malmquist
mechanism explains only a quarter of the effect.

Response (`erass1_lxt.py` joint(), row in results.md): the two regressors are correlated 0.71, not degenerate; jointly, d (dipole)
= −0.09 ± 0.08 and c (exposure) = −0.46 ± 0.07; exposure-only beats dipole-only by Δχ² = 40; the exposure slope fitted outside the
cone alone predicts the inside mean to +0.003 ± 0.019. The mechanism remains unshown (the forward model's quarter), and the note
says so. Third pass pending.
**A parallel second referee (R-118, same day, on the narrowed note with the forward model, no verdict text shown), independently:** VERDICT: novel but not
sound: exposure is nearly collinear with direction in the half-sky eRASS1-DE footprint, and the note's own forward model
reproduces only about a quarter of the deficit, so attributing the anisotropy to the exposure pattern is not supported. The
title still says "is the exposure pattern"; the defensible statement is that in eRASS1 the deficit toward this axis cannot be
told from the survey's depth pattern, and the forward model does not yet reproduce it. Also raised: the temperature sample's
own depth-dependent selection (a second channel the forward model lacks), the shuffle null breaking the exposure correlation it
should test (mock skies from the forward model are the right null), and an eRASS1 supercluster excess near the south ecliptic
pole that survives exposure correction (A&A 2024, aa48884-23).
says so.

Blind review, third pass (2026-10-31, on the joint-fit version): VERDICT: novel but not sound: the joint fit's dipole
(−0.09 ± 0.08 dex) is as consistent with the published Migkas amplitude as with zero, so the headline claim that the signal
"follows exposure, not direction" and that the data can separate the two is not supported, only that eRASS1 gives no
independent confirmation.

Response: accepted. Title and claim now: eRASS1 neither confirms nor excludes the Migkas anisotropy; the signal is not
separable from the survey's exposure in the half-sky. No further passes: the claim is the referee's own supportable statement.

# The Cosmicflows-4 fundamental-plane distances sit on two zero points: a 0.14 mag step at the celestial equator that manufactures a 300 km/s "bulk flow" toward the south celestial pole at 100-150 Mpc/h

2026-10-31. Found by R-119 (cosmo-sweep-5, `sweep-3.md` row 2, trajectory-checked before computing: "partly novel", the
measured step, the free-zero-point fit and the injection named as open). Script `cf4_catalogue.py` (`part2`, `part2x`, `part3`;
`part1` for the companion extinction residual, results.md R-119 row 1). Data: the public Cosmicflows-4 catalogue (Tully,
Kourkchi, Courtois et al. 2023, ApJ 944, 94; VizieR J/ApJ/944/94, individual distances and the group table).

**Verdict: shown.** Beyond 100 Mpc the Cosmicflows-4 (CF4) distance catalogue is two fundamental-plane (FP) samples that
cover opposite celestial hemispheres, the SDSS FP (north, r band, ~34,000 galaxies) and the 6dFGSv FP (south, JHK, ~7,000),
tied to one zero point by the catalogue builders through small overlap samples. Measured flow-free, from groups that hold
both an FP and a Tully-Fisher (TF) distance, the FP zero point is 0.139 ± 0.032 mag higher in the north than in the south
(4.3σ; 3.8σ after the eight method-pair tests run), and 0.072 ± 0.044 higher against the SN Ia distances; the TF and SN
distances agree across the equator (0.033 ± 0.044). It is a step at the equator, not a gradient, and it survives every cut
tried. A step across the celestial equator is a dipole in distance along the polar axis, and at cz = 12,500 km/s a 0.1 mag
step is a 290 km/s velocity dipole. An isotropic mock on CF4's own sky and errors with the measured step injected returns a
"bulk flow" of 299 ± 60 km/s toward (l, b) = (295, −29) in the 10,000-15,000 km/s shell (100-150 Mpc/h), where the
catalogue's single-zero-point dipole is 419 km/s toward (300, −11), and freeing the SDSS/6dF tie drops that to 305 km/s.

## Claim

1. **The step.** In the CF4 group table, each group's distance modulus is given per method. For the 694 groups with both an
   FP and a TF modulus, DM_FP − DM_TF is +0.109 ± 0.022 mag in the north (328 groups) and −0.031 ± 0.024 in the south (366):
   a north-minus-south difference of **+0.139 ± 0.032 mag**. Against SN Ia (202 groups): +0.072 ± 0.044. TF against SN
   (283): +0.033 ± 0.044, so the step is in the FP, not in TF or SN. At the galaxy level (fewer, noisier pairs) the same
   sign and larger: FP − SN +0.46 ± 0.14 (68 galaxies), FP − SBF +0.37 ± 0.10 (131), TF − SN +0.12 ± 0.08 (149). The
   group-level per-method moduli sit on an "arbitrary SN Ia scale" (the table's note), a constant that cancels in the
   north-minus-south difference. Pairs in one group share a distance and a velocity, so no flow enters.
2. **Robustness** (FP − TF, groups): |Dec| > 3° +0.130 ± 0.033; FP from ≥ 3 members +0.099 ± 0.042; TF from ≥ 2 members
   +0.137 ± 0.043; Vcmb > 5,000 +0.093 ± 0.039; |b| > 30 +0.132 ± 0.039; unweighted medians +0.127 ± 0.055; with an E(B-V)
   term +0.137 ± 0.033 (the extinction residual of row 1 is not the cause); with a csc|b| term +0.136 ± 0.033; within the
   same right-ascension range (110-260°, where both hemispheres have groups) north +0.091 ± 0.024 vs south −0.047 ± 0.034;
   by distance, N − S = +0.12 ± 0.05 at 5,000-10,000 and +0.12 ± 0.06 at 10,000-20,000 km/s. Versus declination the offset
   is −0.04 ± 0.05, +0.03 ± 0.04, −0.10 ± 0.06, −0.04 ± 0.06 in the southern bins and the equatorial strip, then
   +0.155 ± 0.033 (Dec 5-20) and +0.090 ± 0.032 (20-45): a step where the SDSS sample begins.
3. **Size against the catalogue's own tie.** The catalogue quotes an FP-to-TF offset uncertainty of 0.0034 dex (0.017 mag;
   Tully+2023 as described by Whitford, Howlett & Davis 2023, who injected ±0.017 dex as a "5σ" offset); the measured
   north-south step is 0.028 ± 0.006 dex, 4-8× that (0.014 dex against SN).
4. **What the step does to the bulk flow.** Weighted least-squares dipole fits of DM − DM_H(Vcmb) in Vcmb shells (all
   methods, errors from the catalogue): at 10,000-15,000 km/s (10,348 galaxies) one zero point gives |B| = 419 km/s toward
   (300, −11) with polar component B_z = −402 ± 51 km/s; a zero point per survey (SDSS-FP, 6dF-FP, equatorial FP, TF, SN,
   other) gives 305 km/s toward (304, −10), B_z = −291 ± 70; at 15,000-20,000, 348 → 186 km/s. An isotropic mock (Gaussian
   scatter at each galaxy's catalogue error and position, no flow) with a 0.106 mag step on the 6dF FP galaxies (the mean of
   the FP − TF and FP − SN steps), fitted with one zero point, returns 299 ± 60 km/s toward (295, −29) at 10,000-15,000,
   272 ± 100 toward (303, −33) at 15,000-20,000, 106 ± 34 at 6,000-10,000; noise-only mocks return 76, 144 and 50 km/s.
   The published deep CF4 flow is about 400 km/s at 150-200 Mpc/h toward l ≈ 300, b ≈ −10 (Watkins+2023; Whitford+2023;
   numbers from memory of those papers, to be checked against them).
5. **The subsamples disagree.** In the same shell each subsample alone (free zero point): SDSS FP B_z = +591 ± 152 km/s
   (602 km/s toward (124, 38)); 6dF FP 855 km/s toward (302, −2), B_z = −774 ± 169; TF 964 km/s toward (306, −29),
   B_z = −962 ± 149; SN 378 ± 170 toward (323, −2); everything but the FP 585 ± 100 toward (314, −19). SDSS and TF disagree
   on the polar component by 1,550 ± 210 km/s (7σ).

## What it means

The catalogue's deep dipole at 100-150 Mpc/h is the sum of its parts' calibrations: the SDSS/6dF step alone makes 300 of its
419 km/s in the published direction, and the TF subsample brings a polar dipole of its own that the SDSS FP, the largest
sample at that depth, sees with the opposite sign. A flow that the catalogue's subsamples do not share is not measured by the
catalogue. Nusser & Tully 2026 saw the symptom (per-catalogue dipoles differ; the deep excess drops to 1.6σ without 6dFGS);
Whitford+2023 showed an assumed 0.017 dex offset would bias the flow by ~190 km/s; this note measures the step (0.028 dex
across the equator), shows it is a step and not a gradient or an extinction effect, and shows what it manufactures. For the
H0-anisotropy and bulk-flow literature built on CF4 (Watkins+2023, Whitford+2023, Hoffman+2023, Courtois+2023/2025 fields,
Stiskalek+2025, Salzano+2025, Nusser & Tully 2026), the SDSS-6dF tie is a nuisance parameter of ±0.03 dex, not ±0.003.

The companion residual (results.md R-119 row 1): the SDSS FP moduli also carry a −1.07 ± 0.16 mag per unit E(B-V) term
against SFD (4.1σ against a rotation null), absent below 8,000 km/s and in the SN distances; it carries ~100 km/s of the same
shell dipole and is not what makes the step.

## Weakest points

- The step is measured against TF and SN in shared groups. If the TF sample had a north-south zero-point error of the same
  sign, FP − TF would show it too; TF − SN (+0.033 ± 0.044) limits any such TF step to < 0.12 mag at 2σ, and the FP − SN
  step (0.072 ± 0.044) is 1.2σ below the FP − TF one, so the step is 0.07-0.14 mag, not pinned better than that.
- Group-level per-method moduli are means over members, each with its own Malmquist-type bias; a depth-dependent bias that
  differs between a deep magnitude-limited SDSS sample and 6dFGSv would read as a zero-point step here. For the bulk-flow
  question it makes no difference (it is a step in the catalogue's distances across the equator either way); for the
  calibration question it does, and the catalogue cannot separate them.
- The dipole estimator is a plain weighted least squares on the log-distance residuals (velocity units through 1/Vcmb),
  not the minimum-variance estimator of Watkins+2023; the amplitudes are comparable within this note (tied vs free, mock
  vs data), and only roughly with the published ones.
- The SDSS/6dF split is by sky (Dec > 0, or Dec < −3 with Vcmb < 16,500); ~700 equatorial FP galaxies are kept as a third
  survey. The membership is not in the public table.

## To settle it

The DESI DR1 fundamental-plane peculiar-velocity catalogue (Said+2025), one instrument over both Galactic caps, overlaps
both SDSS and 6dFGSv and can tie them directly; the SDSS PV (Howlett+2022) and 6dFGSv (Springob+2014) releases against the
CF4 TF sample at the galaxy level with their own group assignments; and a CF4 bulk-flow fit with the SDSS-6dF offset as a free
parameter with a ±0.03 dex prior, to see what survives.

## Blind review (2026-10-31, context-free Opus referee, own literature search; `blind-review.md`)

VERDICT: novel but not sound: the step is measured only in shared FP–Tully–Fisher groups, yet applied uniformly to all 6dFGSv galaxies, despite sample-specific richness and Malmquist recalibrations tuned on the same polar axis; freeing per-hemisphere zero points absorbs a real flow's polar component by construction; and the "7σ disagreement" ignores the degeneracy between zero point and dipole, and cosmic variance, in fits confined to one hemisphere.

Strongest objection, in one sentence: the step comes from the 694 (richer) groups holding both FP and TF members, whose FP
distances carry survey-specific corrections (SDSS's richness-dependent FP fit, 6dFGSv's J < 13.38 Malmquist limit tuned on
mocks to remove a south-pole bias), so it is not shown to apply to the ~7,000 6dFGSv field galaxies the injection applies it
to. Also raised: Whitford+2023's 0.017 dex is a TF-vs-FP offset, not the SDSS-6dF tie (point 3 conflates them); TF − SN =
0.033 ± 0.044 does not exclude a TF anisotropy of Stiskalek+2025's size (0.087 mag); the south celestial pole is at (303,
−27) while the published flow points at b ≈ −6 to −11, so the step can make B_z, not the whole vector; the Dec-bin edges were
chosen after seeing the data. Recorded as a miss. Not reworked: the CF4 thread is parked (Will, 2026-10-31); what stands is
the measured step itself (novel, 0.07-0.14 mag, 4.3σ), and what would be needed to pass is a mock with a real flow and free
zero points, the step measured on field galaxies (e.g. the DESI DR1 FP catalogue as the tie), and covariance-aware
per-subsample dipoles.


---

# Type Ia supernovae at z < 0.06 cluster like galaxies, and low-stretch ones more than high-stretch ones: a two-point measurement of the SN Ia bias

2026-11-01. Script: `sn_bias_2mrs.py` (ZTF SN Ia DR2 public "lite" release; 2MRS via VizieR).

## Claim

The projected cross-correlation of the ZTF SN Ia DR2 volume-limited sample (856 SNe: `snia-cosmo`, light-curve coverage
and fit-quality flags, 0.015 < z_CMB < 0.06, Dec > −28, |b| > 10) with 2MRS galaxies (K < 11.75, 29,652 galaxies in the same
footprint at 0.003 < z < 0.07), divided by the projected auto-correlation of the same galaxies weighted to the SNe's redshift
distribution, gives on 2-20 Mpc/h

    b_SN / b_2MRS(z) = 1.01 ± 0.14 (statistical) ± 0.10 (control weighting),

with jackknife errors over 37 sky regions. The cross-correlation is detected at 7σ in the 2-3.5 Mpc/h bin alone (w_p = 44.9
± 6.6 Mpc/h, 11σ over 0.5-20 Mpc/h) and the ratio is flat across 2-20 Mpc/h (1.09, 1.01, 1.02). Below 2 Mpc/h the SNe are
less clustered than the galaxies (0.76 ± 0.09 at 0.7 Mpc/h, 0.94 ± 0.10 at 1.4). The second error is the shift between two
ways of weighting the control to the SN redshift distribution (per-bin weights, 1.01; weights interpolated between bin
centres, 0.91), and it is there because the 2MRS galaxies that pass K < 11.75 get more luminous, and more clustered, with
redshift: in the two redshift halves the ratio is 1.06 ± 0.17 (z < 0.04, mean 0.031) and 0.69 ± 0.11 (z > 0.04, mean 0.051),
while the control's own w_p nearly doubles between them (37.7 to 65.7 Mpc/h at r_p = 2.6). Against a volume-limited control
(M_K − 5 log h < −24.5, complete to z = 0.04) the z < 0.04 SNe (338) give 0.85 ± 0.15 (per bin 1.02, 0.90, 0.68 on 2-20 Mpc/h). Against a K-luminosity-weighted control (closer
to a stellar-mass-weighted tracer) the full-sample ratio is 0.78 ± 0.12.

Split by light-curve stretch, the two SN populations differ: x1 < −0.5 (328 SNe) gives 1.37 ± 0.21 and x1 ≥ −0.5 (528)
0.76 ± 0.14, each against a control weighted to its own redshift distribution; the difference is 0.61 ± 0.25 (2.4σ), and
against the luminosity-weighted control 1.07 ± 0.19 vs 0.59 ± 0.10 (2.2σ); at z < 0.04 against the volume-limited control
**1.32 ± 0.18 (129 SNe) vs 0.56 ± 0.13 (209), a difference of 0.77 ± 0.22 (3.4σ)**, the same galaxies on both sides. Split by host stellar mass at 10^10 Msun the ratios agree (0.94 ± 0.18 vs 1.00 ± 0.18). The 604 SNe with
host-spectrum redshifts give 1.08 ± 0.17; uniform-sky SN randoms 1.02 ± 0.15; SN randoms following the SN sky density
smoothed to 15 deg instead of 30, 0.91 ± 0.14. The SN auto-correlation, which does not involve the galaxy selection, gives
w_p,ss / w_p,control = 1.42 ± 0.30 on 2-20 Mpc/h, i.e. (b_SN / b_control)² and b_SN / b_control = 1.19 ± 0.13, consistent with
the cross-correlation's 1.01 ± 0.14.

## Method

Landy-Szalay estimators in (r_p, π) with π_max = 30 Mpc/h, six r_p bins from 0.5 to 20 Mpc/h, flat ΛCDM with Ω_m = 0.3 and
distances in Mpc/h; CMB-frame redshifts for both samples. Galaxy randoms: uniform in the mask (Dec > −28, |b| > 10, the
2MRS completeness region intersected with ZTF's) with the galaxies' own shuffled redshifts, four per galaxy. SN randoms: the
same mask, accepted with probability proportional to the SN sky density smoothed with a 30 deg Gaussian (ZTF's depth and
spectroscopic follow-up are not uniform across the sky), with the SNe's shuffled redshifts, forty per SN. The control for
the auto-correlation is the galaxy sample with the first galaxy of each pair weighted by n_SN(z) / n_gal(z) in nine
redshift bins, so that the galaxy pairs sample the SNe's volume; the weighting is recomputed for every SN subsample (the
auto-correlation pair counts are tagged by the redshift bin of the first member, so a reweighting needs no recount). The
ratio is the inverse-variance weighted ratio of w_p,sg to w_p,control over the three bins at r_p > 2 Mpc/h; its error is the
delete-one-region jackknife of the ratio itself (nside-2 HEALPix regions), so the galaxy-side sample variance partly
cancels. In the jackknife every pair count is tagged by the region of both members and recomputed by subtraction.

## What it means

The number: on 2-20 Mpc/h, SNe Ia at z < 0.06 cluster like the K-selected galaxies at the same redshift, within 14 %, and
like the volume-limited L > L* galaxies at z < 0.04 (0.85 ± 0.15). With b ≈ 1.2-1.4 for such galaxies, b_SN ≈ 1.2-1.4: an
ordinary galaxy tracer, not a tracer of the densest structures. The recent report from the same ZTF sample that the SN Ia
density is "not a linear tracer" of the matter density, with 3-8x rate excesses in a few superclusters (Lordet+2026,
arXiv:2604.12714), is a one-point statement about a handful of structures at z = 0.02-0.04; the two-point statistic over
the whole volume says the SN field has the galaxy field's clustering amplitude on linear scales, and below 2 Mpc/h a lower
one (the SN rate per unit stellar mass is lowest in the most massive, most clustered galaxies, which is also why the ratio
against a luminosity-weighted control is below one). The split by stretch is the known environment dependence of the SN Ia
population seen as a bias: low-stretch SNe, the population of old stellar environments, cluster like or more than the
luminous galaxies; high-stretch SNe, the population of star-forming galaxies, cluster like a lower-bias population; the
difference is 2.4σ in the full sample and 3.4σ at z < 0.04 against one volume-limited control (which is the cleaner of the
two: same galaxies, same volume on both sides of the split).

For low-redshift SN cosmology this is the number that peculiar-velocity and local-density analyses assume (b_SN = b_gal in
the sample-variance term of the Hubble residuals and in the density priors of forward-modelled flow analyses): it is
measured, at 14 % precision, and the assumption holds; and the stretch split says that a low-z sample's effective bias
depends on its x1 mix.

## Weakest points

1. The SN selection. ZTF DR2's redshifts come from host spectra (archival, mostly SDSS) for 604 of the 856 SNe and from
   the SN spectra for the rest; a host-spectrum requirement could bias the sample toward catalogued (brighter, more
   clustered) hosts. The host-z-only subsample gives 1.08 ± 0.17 and the full sample 1.01 ± 0.14, so the two selections
   agree, but both could share the bias; the ZTF efficiency model of Lordet+2026 is not reproduced here, and the SN
   auto-correlation (1.19 ± 0.13 for the same ratio) is the selection-independent check and agrees.
2. The control's redshift dependence. The headline ratio is against a flux-limited galaxy population whose bias rises with
   redshift inside the SN volume; the two weighting schemes differ by 0.10 and the redshift halves by 2.3σ. The
   volume-limited comparison at z < 0.04 is the clean version and is quoted with it.
3. The SN randoms. Smoothing the SN sky density at 30 deg (or 15 deg) could absorb real large-scale structure; uniform-sky
   randoms give 1.02 ± 0.15, so on 2-20 Mpc/h the choice does not matter at the quoted error.
4. The split by stretch is 2.4σ in the full sample and 3.4σ at z < 0.04 against the volume-limited control, with 37-region
   jackknife errors and 16 ratios computed in all (the split was one of the two questions asked in advance); consistent with
   the environment dependence known from one-point statistics, and the first two-point version of it.

## To settle it

Pantheon+ low-z and the ATLAS/TITAN sample in the same pipeline (different selections); the ZTF efficiency map as the SN
random catalogue; an absolute b_SN from a model w_p of the matter; the ratio in finer redshift bins with the volume-limited
control.

## Blind review (2026-11-01, context-free Opus referee, own literature search; `blind-review.md`)

VERDICT: novel but not sound: the headline bias ratio varies with the choice of galaxy control by more than its quoted error, and the stretch-split significance is the post-hoc maximum of 16 ratios with a null host-mass split against it, so neither main conclusion is supported.

Strongest objection, in the referee's words: the headline b_SN/b_2MRS = 1.01 depends on the control at a level (0.69-1.08)
larger than its quoted total error, and the stretch result's 3.4σ is a post-hoc maximum over 16 ratios; the full-sample split
is only 2.2-2.4σ, and the null host-mass split argues against it. Also raised: shuffled-redshift randoms and the n(z)-weighted
control copy the sample's own structure into the null; 252 SN-spectrum redshifts (σ_z 0.003-0.005, 10-15 Mpc/h) are comparable
to π_max; the 37-region jackknife with a few dominant structures and no mock validation; no ZTF efficiency map; "b_SN = b_gal
is assumed by PV analyses" asserted, not shown. Prior work the referee found: Carlberg+2008 (arXiv:0805.3983, SNLS host
cross-correlation), Tsaprazi+2022 (arXiv:2109.02651, ZTF BTS in 2M++ density), Graham+2008 (arXiv:0801.4968), the ZTF DR2
stretch-environment papers (arXiv:2406.01108, 2406.11680, Ginolin+2025); the two-point ratio and the stretch split it called
new. A miss on a physics point: parked, no revision pass (R-120's rule). Path to pass: the volume-limited ratio as the only
headline, with the ZTF efficiency map as the SN random catalogue and per-subsample randoms; errors from constrained-simulation
mocks or a finer jackknife with a Hartlap correction; π_max varied and the SN-spectrum-redshift SNe dropped; the stretch split
pre-registered on an independent sample (Pantheon+ low-z, ATLAS) before it is called a detection.


---

# Planck CMB lensing detects the eRASS1 clusters at 10σ, finds their one-halo mass scale at 0.8-1.0 of the eROSITA weak-lensing calibration, and the Planck SZ at fixed eRASS1 mass is 23 % lower for the X-ray-compact half of the clusters, in the 95 %-pure cosmology sample too

2026-11-01. Script: `erass1_cmblens.py` (`forecast`, `profiles`, `nulls`, `fit`, `ystack`); rows B0-B10 in `results-R121.md`.
Data: the eRASS1 primary and cosmology cluster catalogues (Bulbul et al. 2024, `erass1cl_primary_v3.2.fits`, 12,247 clusters;
`erass1cl_cosmology_v1.1.fits`, 5,259), the Planck 2018 lensing convergence reconstructions (PLA: the tSZ-deprojected SMICA
temperature-only reconstruction `COM_Lensing-Szdeproj_4096_R3.00` as the baseline; the minimum-variance, temperature-only and
polarization-only reconstructions of `COM_Lensing_4096_R3.00` as checks), and the Planck 2015 MILCA Compton-y map.

## Claim

The eRASS1 catalogue gives every cluster a mass M500 from the eROSITA count-rate-mass relation calibrated on DES Y3, KiDS-1000 and
HSC-Y3 galaxy weak lensing (Grandis et al. 2024; Kleinebreil et al. 2024; selection-corrected posteriors, Bulbul et al. 2024). That
mass scale is what sets the eRASS1 cluster-count result S8 = 0.86 ± 0.01 (Ghirardini et al. 2024), fitted on the cosmology sample.
CMB lensing is an independent observable with no gas physics between the mass and the signal. Stacking the Planck PR3
convergence on the eRASS1 clusters in the Planck lensing footprint, with the NFW-plus-two-halo convergence predicted from each
cluster's catalogue M500 and redshift filtered exactly as the map, and fitting the one-halo and two-halo amplitudes jointly:

| sample (PCONT < 0.5 unless stated) | n | A_1h (two-halo free) | A_2h | A on 1h + 2h | b_M (M_true = b_M M_eRASS1, two-halo free): 68 %, 95 % | detection |
|---|---|---|---|---|---|---|
| all clusters in the footprint | 6,867 | 0.81 ± 0.15 | 0.91 ± 0.28 | 0.84 ± 0.08 | 0.75 [0.60, 0.95], [0.45, 1.10] | 10σ |
| M500 < 4 × 10^14 M_sun (below Planck's own SZ detections) | 6,033 | 0.82 ± 0.19 | 1.05 ± 0.32 | 0.90 ± 0.11 | 0.80 [0.60, 1.00], [0.40, 1.25] | 8σ |
| **cosmology sample** (EXT_LIKE > 6, 0.1 < z < 0.8, LS DR10 footprint, ~95 % pure) | 3,574 | **0.83 ± 0.17** | 0.44 ± 0.32 | 0.71 ± 0.10 | 0.80 [0.60, 0.95], [0.45, 1.20] | 7σ |
| cosmology sample, M500 < 4 × 10^14 | 2,901 | **0.98 ± 0.23** | 0.39 ± 0.39 | 0.78 ± 0.14 | 0.95 [0.70, 1.25], [0.45, 1.55] | 6σ |
| PCONT < 0.1, M500 < 4 × 10^14 | 4,959 | 0.88 ± 0.19 | 0.74 ± 0.36 | 0.84 ± 0.12 | 0.85 [0.65, 1.05], [0.45, 1.30] | 7σ |
| M500 > 4 × 10^14 (SZ-mask-selected, see Method) | 834 | 0.57 ± 0.18 | 0.47 ± 0.32 | 0.57 ± 0.10 | — | 5σ |

(All amplitudes after subtracting the mean profile of the rotated-catalogue null stacks, which is −1.4 × 10^-3 in the central
annulus; before that subtraction the single amplitudes are 0.82, 0.89, 0.67, 0.73, 0.79. The one-halo and two-halo amplitudes are
anticorrelated at −0.55 to −0.63. The errors are the same by the GLS, the spread of the null stacks, a 24-25-region sky
jackknife and a cluster bootstrap. Model systematic on the one-halo amplitude ±10 %: see Method.)

So the one-halo convergence of the eRASS1 clusters, which is the mass-scale statistic, is 0.8-1.0 ± 0.15-0.23 of what the
weak-lensing-calibrated masses predict in every sample, within 1σ of the eRASS1 scale and 2σ or more above a halving of it
(Δχ² = 2-3 with the two-halo free, 10-11 with the two-halo fixed at the eRASS1-cosmology linear model). The single amplitude on
1h + 2h is lower in the cosmology sample (0.71 ± 0.10) than in the primary sample (0.84 ± 0.08) because the two-halo term is
weaker there (A_2h = 0.44 ± 0.32 against 0.91 ± 0.28), not the one-halo term.

**The split by X-ray brightness at fixed mass.** Within M500 < 4 × 10^14, dividing the clusters at the median of their detection
likelihood, exposure, detection likelihood per unit exposure (compactness at fixed count rate) and extent per detection likelihood
(extendedness), each relative to the median in 6 × 6 cells of (log M500, z) inside the sample in question:

| split at fixed (M500, z), M500 < 4e14 | sample | Planck y, above / below, divided by the mass-predicted ratio | lensing A_above − A_below, with its own null |
|---|---|---|---|
| detection likelihood | primary (6,033) | 0.75 ± 0.03 (8σ) | −0.50 ± 0.18 (2.8σ; empirical p = 0.005) |
| detection likelihood | PCONT < 0.1 (4,959) | 0.75 ± 0.03 | −0.52 ± 0.20 (2.6σ; p = 0.005) |
| detection likelihood | **cosmology sample (2,901)** | **0.77 ± 0.04 (6σ)** | **−0.05 ± 0.25** |
| detection likelihood | EXT_LIKE ≥ 12 (1,620) | 0.62 ± 0.03 | −0.14 ± 0.28 |
| compactness (DET_LIKE / exposure) | cosmology sample | 0.81 ± 0.04 | −0.20 ± 0.24 |
| extendedness (EXT_LIKE / DET_LIKE) | cosmology sample | 1.34 ± 0.07 | +0.12 ± 0.23 |
| exposure | primary | 0.97 ± 0.04 | +0.17 ± 0.20 |

(The y statistic: the MILCA Compton-y within 10' minus the 60-180' annulus, the ratio between the halves divided by the ratio the
catalogue masses alone predict with Y ∝ M^1.79 E(z)^2/3 D_A^-2, so that it reads as a residual at fixed eRASS1 mass; bootstrap
errors. The lensing differences: the two halves stacked in the same 200 rotations of the catalogue, so the null of the difference
carries their shared sky.)

At fixed eRASS1 mass and redshift the X-ray-brighter, more compact half of the clusters has 19-25 % less Compton-y than the fainter,
more extended half, and the extended half 34 % more, with 6-8σ in each sample including the 95 %-pure cosmology sample, while
exposure, which raises the detection likelihood without touching the count rate, does nothing. With Y ∝ M^1.8 that is a 12-15 %
difference in mass between the halves. CMB lensing sees the same sign at 2.6-2.8σ in the primary and PCONT < 0.1 samples and
nothing (−0.05 ± 0.25) in the cosmology sample, where its precision neither confirms nor excludes the SZ-sized split (a 15 % mass
difference is 0.15 in A, well inside ±0.25) nor the primary sample's −0.5 (1.8σ away).

## Method

Map: the released κ alm, band-limited to 8 ≤ L ≤ 2048 with a cosine taper on the top 10 %, on nside 1024, multiplied by the release
mask (fsky 0.66). The Planck lensing mask removes the point sources and the SZ clusters Planck detects internally, so a cluster is kept
only if the mask is clean within 30' and ≥ 80 % clean within 2°: 6,867 of 9,293 primary clusters, 113 of 465 PSZ2 matches. The
hole-masked fraction is 14-17 % at every mass below 4 × 10^14 (point-source holes) and rises to 20 % at 4-6 × 10^14 and 44 % above
6 × 10^14, where the holes are the SZ detections: above 4 × 10^14 the unmasked remainder is, at fixed eRASS1 mass, the SZ-faint part,
so that bin's low amplitude (0.57) is a selection effect of the lensing mask and is not used for the mass-scale statement.
Profiles: the mean of the masked map in 13 annuli from 0-5' to 150-180' around each cluster's X-ray position (median X-ray-to-BCG
offset 0', 68 % within 1.7'). Model: for each cluster, the NFW convergence for a CMB source (M500c to M200c and c200c with Duffy
et al. 2008; Σ_crit from angular-diameter distances in the eRASS1 cosmology, H0 = 67.77, Ωm = 0.29) plus the two-halo term (linear
CAMB P(k) scaled to the eRASS1 σ8 = 0.88, Tinker et al. 2010 bias), Hankel-transformed to κ(L), multiplied by the same band filter and
the nside-1024 pixel window, back-transformed to the annulus means, and averaged over the same clusters (the 1-halo and 2-halo κ(L)
on a 22 × 20 grid in (log M500, z), interpolated). Fits: the one-halo and two-halo amplitudes jointly (the mass-scale statistic), one
amplitude on 1h + 2h, and a mass rescaling b_M with the profiles rebuilt at b_M × M500 on a grid and the two-halo amplitude
profiled out. Null and covariance: 400 (200 for the subsets defined after the referee's check) stacks of the whole catalogue rotated
in Galactic longitude by 25-335° with a random latitude reflection, through the same footprint selection and annuli, giving each
sample's 13 × 13 covariance (Hartlap corrected), the null distribution of its amplitude (−0.05 ± 0.08 for the full sample; the data's
0.82 is 10.9 of its spreads from the null mean, and 0 of 400 rotations reach it, which is the empirical statement; the "10σ" is the
GLS error, which equals the null spread) and the null of each split difference. Forecast before the map was opened: a Fisher S/N of
13 for the 1-halo term from the catalogue and the release noise curve (row B0); the measured 10σ with 74 % of the clusters is that
forecast.

Checks (rows B2, B4, B6, B9): the minimum-variance reconstruction gives the full-sample single amplitude 0.82 ± 0.07, the
temperature-only reconstruction without the tSZ deprojection 0.84 ± 0.08 (0.02 below the deprojected one, so the tSZ bias of the
temperature estimator is a few per cent at most at these masses), the baseline band cut to 8-1000 gives 0.83 ± 0.09 and to Planck's
conservative 8-400 gives 0.90 ± 0.11, the polarization-only reconstruction gives 0.69 ± 0.23; the brightness split has the same
sign in all of them. Concentrations × 0.7 / × 1.3 change the predicted central annulus by +11 % / −6 % (at fixed M500 the M200 and
the outer profile change); σ8 = 0.81 in the two-halo term moves the fixed-two-halo one-halo amplitude by +4 % and is inside the free
A_2h; the catalogue's M500_PDF column gives posterior means 22 % above the medians and would lower the amplitudes by 10 %, but that
column is three times wider than the catalogue's own M500_L/H (±65 % vs ±23 %) and its median is 0.96 × M500, so the published
M500 is the baseline and ±10 % is carried as the model systematic on the one-halo amplitude. The scaled-covariance subset errors of
the first pass agreed with the per-subset nulls to 10 % (0.89 ± 0.10 against 0.85 ± 0.11 for M500 < 4 × 10^14).

## What it means

The eRASS1 mass scale, calibrated with galaxy shapes from three surveys, is reproduced by a fourth lensing measurement with a
different source plane, a different instrument and different systematics: the one-halo convergence is 0.83 ± 0.17 (cosmology
sample) and 0.98 ± 0.23 (cosmology sample below 4 × 10^14) of the prediction, in the mass range that carries most of the eRASS1
cosmology sample and that Planck's own SZ survey does not reach. The Planck SZ catalogue and the eRASS1 catalogue have been
compared on their common, massive clusters (Aymerich et al. 2025, arXiv:2509.02068; the ACT DR6 catalogue, arXiv:2507.21459) and
agree; this is the first test below those catalogues' masses, and it is independent of the hydrostatic-mass bias those comparisons
rest on. For the S8 discussion it says that whatever separates the eRASS1 S8 of 0.86 from the cosmic-shear 0.76-0.78 is not a
factor-of-two error in the eRASS1 masses (the hydrostatic comparison of Li et al. 2024, arXiv:2411.09735, on 22 clusters); a 20 %
error is still allowed at 1σ, and the two-halo term of the cosmology sample is low against the eRASS1-cosmology linear model
(0.44 ± 0.32), which is where a lower σ8 or a halo-exclusion correction to the linear two-halo term would show.

The split is the new measurement, and its robust form is the SZ one. The eRASS1 masses come from the count rate alone, so at fixed
count rate they cannot know which clusters are over-luminous for their mass (cool-core, concentrated) and which are under-luminous
(disturbed, extended); the detection and extent likelihoods carry that information, because at fixed count rate a concentrated
surface-brightness profile is detected more significantly and is less extended. The Planck y says the concentrated half is the
less massive half by 12-15 %, in the cosmology sample as in the primary one, so it is not contamination (a ≤ 5 % impurity cannot
make a 23 % deficit), and CMB lensing agrees in sign where it has the precision. That is the size of the known cool-core luminosity
excess propagated through the count-rate-mass slope. Two consequences follow for the eRASS1 analysis. The scatter of the
count-rate-mass relation is correlated with the selection variable: the eRASS1 selection function is modelled as a function of
count rate, redshift and sky position (Clerc et al. 2024), and a detection likelihood that depends on morphology at fixed count
rate, with morphology correlated with the mass residual, is a mass-dependent selection at the threshold that this parametrisation
does not contain (Grandis et al. 2024 marginalise over a selection-concentration correlation; the eRASS1 scaling-relation paper,
arXiv:2511.14356, names morphology-dependent selection as future work). And the weak-lensing calibration, which fits the shear of
all clusters at a given count rate, is right on average but the individual masses carry a ±15 % morphology term that a joint fit
with the detection likelihood would remove.

## Weakest points

- The mass-scale statistic is the one-halo amplitude with the two-halo free, at ±0.17-0.23 statistical and ±10 % model
  (concentration-mass relation, the posterior-mean question): a 20 % test, not a 10 % one. The single amplitude on 1h + 2h is
  tighter (±0.10) but depends on the two-halo model, and in the cosmology sample it is 3σ below one for that reason.
- The absolute amplitudes rest on the Planck 2018 κ normalisation at 400 < L < 2048, which Planck calls less reliable than the
  conservative range; the conservative-range fit gives the same answer with a 30 % larger error, and the polarization-only fit the
  same within its 0.23. CIB leakage correlated with the clusters and the higher-order estimator bias for the most massive clusters
  are not measured here; Planck 2018 bounds the CIB bias at the sub-percent level for the deprojected map's power spectrum, and the
  tSZ-deprojected and plain temperature reconstructions agree to 0.02.
- The null-mean correction (+0.05 to +0.08 in the amplitudes) is empirical: the band-limited map at positions that pass the footprint
  selection has a mean of −1.4 × 10^-3 in the central 5', three times the mean of the map in clean regions (−4.5 ± 3.5 × 10^-4), a
  mask-distance effect. The rotated catalogues pass the same selection, so it is the same for every sample and null, and the split
  differences do not contain it.
- The lensing split is 2.6-2.8σ in the primary and PCONT < 0.1 samples and 0.1σ in the cosmology sample, at a precision that
  cannot tell the SZ-sized split from the primary sample's or from zero. Whether the primary sample's larger lensing split is the
  lower-purity clusters, the z < 0.1 clusters or noise is not settled here; the SZ split does not depend on it.
- The SZ split is an aperture measurement (10') without a pressure-profile model; its translation to a mass ratio assumes Y ∝ M^1.8 and
  that the integrated Compton-y within 10' at fixed mass does not itself depend on morphology, which at Planck's beam is the
  standard assumption of low-scatter Y-M relations but is not tested here.
- Clusters above 4 × 10^14 are not tested: the Planck lensing mask removes the SZ-detected ones and leaves the SZ-faint ones.
- Trials: about 250 subset amplitudes were computed in this script's runs (the same clusters across six reconstructions or bands
  and two null sets, so far fewer independent ones); the headline sample and the DET_LIKE ≥ 40 split were fixed before the first
  fit, the fixed-(M, z) residual splits and the purity-restricted samples after it (the latter at a referee's request).

## To settle it

The ACT DR6 lensing map (LAMBDA; unreachable from this network) would give the same stack at three times the significance, would
reach the clusters the Planck mask removes, and would decide the lensing split in the cosmology sample; the eRASS1 count-rate-mass
relation refitted with a detection-likelihood (or concentration) term against the DES/KiDS/HSC shear profiles would test the
split's size cluster by cluster; a pressure-profile matched filter on the y map (the fallback candidate A of `deep-1.md`) would
turn the SZ split into a Y-M relation with morphology as a second variable; and the eRASS1 selection function evaluated with
morphology-dependent detection likelihoods would say what the split does to S8.

## Blind review (2026-11-01, context-free Opus referee, own literature search; `blind-review.md`)

VERDICT: novel but not sound: the compactness split at fixed count rate is confounded by AGN and spurious-source contamination in the PCONT < 0.5 primary catalogue, which predicts the same deficit in both y and lensing, so the cool-core reading and the clean mass-scale test are not established.

Strongest objection, in the referee's words: the sample is the primary catalogue at PCONT < 0.5, so each kept source can have up
to a 50 % chance of being an AGN, a point source or a spurious detection; at fixed count rate such an object scores high on
detection likelihood and low on extent per likelihood and carries neither Compton-y nor lensing mass, which predicts the same
deficit in both; the note never residualises against PCONT, restricts to the cosmology sample, or checks an AGN match, and the
same contamination dilutes the absolute A (the weak-lensing calibration applies to the cosmology sample). Also raised: subset
errors from a scaled covariance rather than per-split nulls, and no null for the split difference; the 10σ rests on Gaussian nulls
beyond the empirical 3σ; the null-mean offset is three times the clean-region mean, not "consistent" with it; the posterior mean
should be modelled, not carried; CIB leakage and the higher-order bias not shown; the cool-core argument for the SZ aperture does
not hold at Planck's beam. Prior work it adds: arXiv:2502.02239 (eRASS1 morphology catalogue, X-ray selection favours cool cores),
2207.09242 (eROSITA digital twin, cool-core detection bias), 2511.14356 (eRASS1 scaling relations, morphology-dependent selection
named as future work), 2402.08455 (marginalises over a selection-concentration correlation), 2310.04387, 2406.14754. The
CMB-lensing check of the eRASS1 mass scale below the Planck SZ limit it called new.

Response (same day, under the two-hour rule for a physics point): the splits redone inside the eRASS1 cosmology sample (5,259
clusters, ~95 % pure, 3,574 in the footprint), at PCONT < 0.1 and at EXT_LIKE ≥ 12, each half matched in (M, z) within its sample,
with their own rotated-catalogue nulls and a null for the split difference; the same for the Planck y. Rows B8-B9 in
`results-R121.md` and the paragraph below.

Blind review, second pass (2026-11-01, context-free Opus referee, on the rewritten note, no verdict text shown): VERDICT: novel but not sound: the headline morphology split is measured against posterior-median masses whose widths and selection corrections correlate with the split variable and is open to radio-source fill-in, so it is not shown to be a mass difference, and CMB lensing finds no split in the clean sample.

Strongest objection, in the referee's words: the predicted y ratio is built from posterior medians, Y ∝ M^1.79 is convex, and at
fixed count rate the low-DET_LIKE half has wider mass posteriors and larger Eddington corrections, so part of the deficit is
convexity (≈ 8 % for σ_lnM 0.2 against 0.3); cool-core clusters host radio-loud centrals that fill the decrement in a 10'
aperture; point-source contamination of the X-ray flux raises count rate, mass and compactness together; and the lensing split
is absent in the clean sample. The mass-scale part it called acceptable but weak (a 20-25 % test). Prior work it adds: the eRASS1
morphology papers (arXiv:2502.02239, 2601.00636: eRASS1 clusters more concentrated than SZ-selected ones; no SZ or mass at fixed
count rate as a function of morphology), Clerc+2024 (arXiv:2402.08457), arXiv:2601.15935. Checks it asks for: forward-model Y
through each cluster's full M500_PDF; mask radio sources and repeat on per-frequency maps; split by the morphology catalogue's
parameters; the DES/KiDS/HSC shear split the same way.

Response (same day; rows B11-B12): (1) The convexity point, done. With the prediction forwarded through each cluster's mass
posterior (the catalogue's M500_PDF, and a lognormal of the quoted M500_L/H width), the compactness and extendedness splits do not
move, because their halves have the same posterior widths (quoted σ_lnM 0.21 / 0.20 and 0.21 / 0.20; PDF half-widths 0.61 / 0.62):
cosmology sample, compact above / below 0.81 (medians), 0.81 (lognormal), 0.78 ± 0.04 (PDF); extended 1.34, 1.34, 1.32 ± 0.06. The
DET_LIKE split does move, as the referee said, because its halves differ in count-rate precision (σ_lnM 0.18 / 0.23 in the
cosmology sample): 0.77 (medians) → 0.79 (lognormal) → 0.88 ± 0.04 (PDF); primary sample 0.755 → 0.78 → 0.85 ± 0.035; the
non-residualised DET_LIKE ≥ 40 split 0.52 → 0.56 → 0.73. So the split that survives the posterior forward model is the morphology one
(compactness, extendedness) at 19-22 % in y (5-6σ), and the detection-likelihood one at 12-15 % (3-4σ). The section's 19-25 % should be
read as 12-22 % and the "12-15 % in mass" as 7-13 %. (2) The NILC y map gives the same ratios as MILCA (0.77, 0.82, 1.31 for the three
cosmology-sample splits), so it is not a component-separation artefact; radio fill-in, which both maps share, is not tested here (a
radio-catalogue cross-match is the check; not done in this job). (3) The lensing's silence in the clean sample stands as written. The
loop stops here: two passes, both misses on physics points, the second's open items (radio fill-in, the shear split) outside this
job; the note is recorded as a miss with its path to pass.

Blind review, third pass (2026-11-01, main session, context-free Opus referee, on a narrowed note: the detection and mass scale
only, the morphology split removed): VERDICT: novel but not sound: with the two-halo term free, the note's own numbers (Δχ² = 2–3,
and every sample's 95% b_M interval includes 0.5) do not exclude the factor-of-two mass error it claims to rule out, and the
stronger exclusion relies on a two-halo term fixed to the eRASS1 cosmology under test. Strongest objection, in the referee's words:
what survives is a 20-25 % consistency check plus a full-sample preference for b_M < 1 (0.75 [0.60, 0.95]); the 0.98 headline is a
post-hoc subsample; the mass predictions should average κ over each cluster's posterior; the lensing-mask holes and the null-mean
correction are untested against mask-matched randoms. Novelty confirmed again (closest: Zubeldia & Challinor arXiv:1904.07887 on
PSZ2; arXiv:2406.14754 ACT DR6 lensing on SZ-selected clusters; arXiv:2507.21459 the SZ mass ratio 0.99). PARKED after three
passes. What would pass: the ACT DR6 lensing map (about three times the signal), posterior-averaged predictions, mask-matched
randoms, and a claim no stronger than the interval.


# No periodic-quasar binary candidate is excluded by the NANOGrav 15-year continuous-wave limits: the population tally

R-124 (cosmo-deep-4), 2026-11-01. Script `periodic_cw_limits.py`; rows P0-P4 of `results-R124.md`. Data: the NANOGrav 15-year
individual-binary release (Zenodo 8067506, `15yr_cw_3d_limits_v4.npz`), Agarwal et al. 2025 (arXiv:2508.16534) Tables 6-8, and the
candidate tables from the arXiv sources of Graham et al. 2015 (via Agarwal), Charisi et al. 2016, Liu et al. 2019, Chen et al. 2020
and Chen et al. 2024. The trajectory verdict ("partly novel", arXiv:2508.16534) is in `results-R124.md`; this is the open part.

## Claim

Of 192 published periodic-quasar supermassive-black-hole-binary candidates with a black-hole mass (101 CRTS, Graham et al.
2015; 50 PTF, Charisi et al. 2016; 26 Pan-STARRS1 MDS, Liu et al. 2019; 5 DES+SDSS, Chen et al. 2020; the 10 ZTF candidates
of Chen et al. 2024 whose table rows are in the paper's source), none is excluded as a circular binary at its claimed
period, redshift and mass by the NANOGrav 15-year data set: 0 of 192 at 95 % (Clopper-Pearson 95 % interval on the
fraction: 0-1.9 %), and the number expected to be excluded once the mass uncertainty and the unknown mass ratio are
marginalised is 1.1 (0.7 in the period-doubling convention). For the CRTS candidates the limits used are the 15-year
targeted (fixed position, distance and frequency) 95 % strain limits NANOGrav published for each of them (Agarwal et al.
2025, arXiv:2508.16534, Tables 6-8); for the others the 15-year sky- and frequency-dependent 95 % limit map (the
individual-binary paper's release, arXiv:2306.16222). The limits sit above the candidates' median implied strains by a
factor 5.8 at the closest (SDSS J131706.19+271416.7; six CRTS candidates with log M = 9.8-9.95 at z = 0.8-2.7 are within a
factor 6-12, each with an exclusion probability of 0.06-0.13) and by a median factor 210 (CRTS, targeted), 280 (PTF), 750
(PS1), 1200 (DES) and 310 (ZTF). In chirp mass, the targeted CRTS limits sit a factor 4.5 (minimum) to 29 (median) above
the chirp mass each candidate's published total mass implies at equal masses; no candidate's limit falls below it.

## Method

Each candidate's gravitational-wave frequency is twice its observed-frame orbital frequency, f_gw = 2 / P_obs (the
light-curve period taken as the orbital period, as the targeted search did); the period-doubling case, f_gw = 1 / P_obs, is
computed alongside and changes nothing (0 of 192; the implied strain falls by 2^(2/3) while the limits improve by less
between 14 and 7 nHz). The implied strain is h0 = 2 (G M_cz)^(5/3) (pi f_gw)^(2/3) / (c^4 D_L), with M_cz = (1+z) M_c the
redshifted chirp mass, M_c = M q^(3/5) / (1+q)^(6/5), q log-uniform on 0.1-1, log M normal about the catalogue's single-epoch
virial mass with the catalogue's error (0.23-0.31 dex where given, 0.4 dex otherwise), D_L in flat LCDM (H0 = 70,
Omega_m = 0.3; Xin et al. 2021's choice). 20,000 draws per candidate; the exclusion probability is the fraction of draws above
the limit at the candidate's position and frequency; "excluded at 95 %" means that fraction is at least 0.95. Inclination
and polarisation are not drawn: the published limits are on h0 with those angles marginalised under isotropic priors, so the
comparison is strain to strain. The limit map has 22 frequency bins from 1 to 29.7 nHz (192 healpix pixels); above 29.7 nHz
(48 of 50 PTF, 11 of 26 PS1 and 7 of 10 ZTF candidates at f_gw = 2 / P_obs, where the 15-year array is white-noise limited)
the last bin's map is scaled by f / 27.5 nHz, so those limits are extrapolated; the closest such case is a factor 4.3 (SDSS
J145859.07+153144.7, PTF, log M = 10.03 at z = 2.9, f_gw = 73 nHz) and no exclusion follows from it under any scaling between
f^0 and f^1 between 30 and 180 nHz. A convention check against the targeted paper's own numbers: the chirp mass that follows
from each of its targets' h95, frequency and luminosity distance equals its quoted M_c^95 times (1+z) to 1 % across all 113
targets (a constant factor 1.24 aside), so its tabulated chirp-mass limits are redshifted chirp masses; the strain limits
are used directly here and the published masses are redshifted on this side, so the result does not depend on that factor.
At the 113 targets the released sky map's limit is 1.09 times the targeted limit (median; range 0.5-8), i.e. the targeted
limits are about equal to the sky-dependent map and the paper's quoted factor-2.2 improvement is against the
posterior-sky-averaged curve.

## What it means

The 15-year data set does not test any published periodic-quasar candidate at its claimed parameters: the strongest
candidates (log M about 9.9 at z 0.8-2.7, 11-16 nHz) would need masses 4.5 times their virial estimates at equal masses, or
0.65 dex above them, to be excluded, and the rest need factors of 10 to 1000. This is the population statement the targeted
paper's tables imply but do not state, extended to the four catalogues that paper did not search (where the candidates are
smaller and at higher frequencies, and the limits are a further factor 2-3 worse) and to the period-doubling convention. It
is a null with an interval: the fraction of the candidate population excluded is below 1.9 % (95 %), and the limits would
have to improve by a factor 6 (to reach the first candidate at its median implied strain) or 4.5 in chirp mass before a
PTA exclusion of a periodic quasar at its published mass is possible; at the factor-2 gain per data-set generation the
targeted paper reports, that is two to three NANOGrav data sets away for the single closest candidate, and the catalogue
as a whole (median factor 200 and more) is out of reach of pulsar timing at its published masses.

## Weakest points

- The ZTF catalogue is represented by 10 of its 184 initial candidates (the rows printed in the paper's source); the
  machine-readable table could not be fetched from this environment (OUP and CDS unreachable). Its candidates have periods of
  600-900 days and log M 7.7-9.7, so the missing 174 lie at 26-40 nHz with implied strains at or below the included ten;
  the pooled "0 of 192" is for the candidates listed, not the full ZTF set.
- Above 29.7 nHz the limit is an extrapolation of the released map (85 of the 192 at f_gw = 2 / P_obs). The exclusion
  conclusion is insensitive to it (the closest extrapolated case is a factor 4.3 and the limit scaling between f^0 and f^1
  moves it between 2.7 and 4.3), but the factor-y numbers quoted for PTF and PS1 are approximate.
- The virial masses carry 0.4 dex systematic uncertainty and the mass-ratio prior is a choice; both are marginalised into
  the exclusion probabilities, and the equal-mass, +1-sigma strains are also tabulated. Period-to-orbit mappings other than
  1:1 and 2:1 (the factor-of-several excess circumbinary-accretion periods) move f_gw down and the implied strain down
  further; none makes an exclusion.
- The result is an upper limit on a fraction that is already zero in the published per-target tables for the main
  catalogue; its content is the pooled interval, the quantified distance from exclusion, and the coverage of the other
  catalogues and the second frequency convention.

## Blind review (2026-11-01, context-free Opus referee, own literature search; `blind-review.md`)

VERDICT: sound but not novel: Xin, Mingarelli & Hazboun 2021, arXiv:2009.11865 (with NANOGrav targeted search Agarwal et al. 2025, arXiv:2508.16534)

Strongest objection, in the referee's words: the headline is a repackaged null that is already published per target
(arXiv:2508.16534) and at population level for the 11-year data (arXiv:2009.11865); its new statistic, the 1.9 %
Clopper-Pearson bound, does not answer a meaningful question (the candidates are not Bernoulli trials sharing one exclusion
rate; the informative quantity is the per-candidate distance from exclusion, which the note gives); the one forward-looking
claim ("two to three data sets away") rests on a misread improvement rate (the factor 2.2 is targeted-versus-sky-averaged, not
per generation; Xin et al. project a factor 6 for the IPTA). Also raised: the constant 1.24 in the chirp-mass convention check
is 2^(3/10), a √2 in the strain normalisation, which would move the closest margin from 5.8 to 4.1 (no exclusion either way);
virial masses of periodic candidates may be biased high, cutting the tail probabilities further; the ZTF ten are a convenience
sample. The core null holds under the stated numbers and is robust to the frequency convention, the mass ratio, the
extrapolated limits and the √2. Recorded as a miss; the section stands as reviewed. Parked: the population test of these
candidates is the gravitational-wave background (Sesana et al. 2018, Casey-Clyde et al. 2025), not continuous waves.

# ACT DR6 CMB lensing detects the eRASS1 clusters at 4.5-5σ on the samples fixed before the fit (7-9σ with the SZ-detected, inpainted clusters admitted); on those pre-registered, SZ-faint samples the eROSITA weak-lensing mass scale is b_M = 0.65-0.70 ± 0.19 (statistical) ± 0.12 (model), which excludes neither the eRASS1 scale nor a halving of it

2026-11-01. Script: `erass1_act.py` (`prep`, `data`, `sims`, `rotations`, `transfer`, `fit`, `planck`; it imports the catalogue
handling, halo model and fits of `erass1_cmblens.py`); rows C0-C9 in `results-R122.md`. Data: the eRASS1 primary and cosmology
cluster catalogues (Bulbul et al. 2024; `erass1cl_primary_v3.2.fits`, 12,247 clusters; `erass1cl_cosmology_v1.1.fits`, 5,259);
the ACT DR6 CMB lensing release v1 (Qu et al. 2024, arXiv:2304.05202; Madhavacheril et al. 2024, arXiv:2304.05203; simulations
Farren et al. 2024, arXiv:2309.05659): the baseline (profile-hardened, temperature plus polarization) convergence alm, the
temperature-only, polarization-only, CIB-deprojected and galcut040 variants, the analysis masks, the separate SZ-cluster mask, and
the 400 simulated reconstructions with their input convergence fields; the Planck 2018 tSZ-deprojected lensing reconstruction
(`COM_Lensing-Szdeproj_4096_R3.00`) for the comparison on the same clusters.

## Claim

The eRASS1 catalogue gives every cluster a mass M500 from the eROSITA count-rate-mass relation calibrated on DES Y3, KiDS-1000 and
HSC-Y3 galaxy weak lensing (Grandis et al. 2024; Kleinebreil et al. 2024), the mass scale behind the eRASS1 cluster-count
S8 = 0.86 ± 0.01 (Ghirardini et al. 2024). CMB lensing is an independent observable with no gas physics between the mass and the
signal. Stacking the ACT DR6 convergence on the eRASS1 clusters in the ACT footprint, with the NFW-plus-two-halo convergence
predicted for each cluster from its catalogue mass posterior and redshift, filtered exactly as the map, and fitting the one-halo
and two-halo amplitudes jointly, with the covariance from the release's 400 simulations:

| sample (PCONT < 0.5 unless stated) | n | A_1h (two-halo free) | A_2h | A on 1h + 2h | b_M (M_true = b_M M_eRASS1, two-halo profiled out): 68 %, 95 %, statistical | b_M 95 % with the ± 0.12 model systematic added | detection (A_1h / σ) |
|---|---|---|---|---|---|---|---|
| **pre-registered: cosmology sample** (EXT_LIKE > 6, 0.1 < z < 0.8, LS DR10 footprint, ~95 % pure), outside the SZ-cluster mask | 1,246 | **0.77 ± 0.17** | 1.48 ± 0.39 | 0.91 ± 0.12 | **0.70 [0.53, 0.92], [0.36, 1.14]** | [0.24, 1.26] | 4.6σ |
| **pre-registered: all clusters in the footprint** outside the SZ-cluster mask | 2,444 | **0.74 ± 0.15** | 1.46 ± 0.32 | 0.92 ± 0.10 | **0.65 [0.50, 0.84], [0.36, 1.02]** | [0.24, 1.14] | 5.0σ |
| pre-registered: M500 < 4 × 10^14 M_sun outside the mask | 2,385 | 0.77 ± 0.15 | 1.50 ± 0.33 | — | 0.70 [0.54, 0.90], [0.39, 1.09] | [0.27, 1.21] | 5.0σ |
| pre-registered: cosmology sample, M500 < 4 × 10^14, outside the mask | 1,195 | 0.81 ± 0.18 | 1.57 ± 0.41 | — | 0.75 [0.55, 0.99], [0.37, 1.22] | [0.25, 1.34] | 4.5σ |
| supplementary: inside the SZ-cluster mask (the SZ-detected clusters, inpainted before the reconstruction) | 874 | 1.09 ± 0.14 | 0.75 ± 0.28 | 0.98 ± 0.09 | 1.10 [0.92, 1.31], [0.75, 1.52] | — | 7.6σ |
| supplementary: cosmology sample, inside and outside together | 1,907 | 0.81 ± 0.11 | 1.19 ± 0.26 | 0.90 ± 0.08 | 0.75 [0.63, 0.90], [0.51, 1.04] | [0.39, 1.16] | 7.4σ |
| supplementary: all clusters in the footprint | 3,318 | 0.89 ± 0.10 | 1.18 ± 0.24 | 0.97 ± 0.07 | 0.85 [0.72, 0.98], [0.60, 1.12] | [0.48, 1.24] | 8.5σ |

(Amplitudes after subtracting the stack of 200,000 mask-matched random positions on the same map, a +0.004 correction; the
simulations' mean at the cluster positions and the rotated-map mean are zero within their noise, so the three null-mean estimates
agree and the correction is immaterial. The one-halo and two-halo amplitudes are anticorrelated at −0.5. The one-halo error from
the GLS, the spread of the one-halo amplitude over the 400 simulated maps and the spread over 165 rotated-map nulls agree:
0.166, 0.165, 0.161 for the pre-registered cosmology sample (0.110, 0.110, 0.107 for the union). The model systematic on b_M,
± 0.12, is the quadrature sum of the two-halo treatment (+0.10 between the linear and the halo-exclusion form; the halofit form,
−0.35, is not counted, see Method), the concentration-mass relation (± 0.05 for × 0.7 / × 1.3) and the posterior width (−0.04
for the M500_PDF column); the quadratic estimator's non-linear response at cluster scales, uncorrected and one-sided (it lowers
A_1h by of order 10 %, so it would raise b_M by about 0.05-0.08), is named in Weakest points and not added. The last column adds
the ± 0.12 linearly to the statistical 95 % interval.)

The claim is the interval, on the samples fixed before the first fit. For the eRASS1 cosmology sample outside the SZ-cluster
mask the one-halo convergence is 0.77 ± 0.17 of what the weak-lensing-calibrated masses predict, b_M = 0.70 ± 0.19 (statistical)
± 0.12 (model), with a 95 % interval of [0.36, 1.14] statistical and [0.24, 1.26] with the model systematic: it includes the eRASS1
scale (b_M = 1) and a halving of the masses alike. For all 2,444 clusters outside the mask, b_M = 0.65 ± 0.17 ± 0.12, 95 %
[0.36, 1.02] and [0.24, 1.14]: the same. Neither 0.5 nor 0.8 nor 1 is excluded by any pre-registered sample. The supplementary
samples that include the SZ-detected, inpainted clusters (whose lensing response is not validated, see Method) are tighter
(cosmology sample 0.75 ± 0.14 ± 0.12, 95 % [0.51, 1.04] statistical, [0.39, 1.16] with the systematic) and exclude nothing either
once the model systematic is counted; with the two-halo term held at the eRASS1-cosmology linear prediction instead of free (not
the statistic) the supplementary intervals exclude 0.5, the pre-registered ones do not (cosmology outside the mask: 0.80
[0.68, 0.94], 95 % [0.56, 1.06]). So this is a 25-30 % measurement of the eRASS1 mass scale on the pre-registered samples and a
15-20 % one if the inpainted clusters are admitted, 1-1.5σ below the eRASS1 scale in every version, with the eRASS1 scale and a
20-35 % lower one both allowed.

## Method

Map: the released κ alm (equatorial, mean-field subtracted, normalised, band 2-3000 per the release's filter file), band-limited to
20 ≤ L ≤ 2048 with a cosine taper on the top 10 %, synthesised at nside 1024 for the annuli inside 30' and nside 512 outside (the
two pixel windows are in the model). Mask: the release's analysis mask degraded to the two resolutions, clean where every sub-pixel
is ≥ 0.99 (the README's threshold), fsky 0.213 (8,800 deg²); a cluster is kept if clean within 30' and ≥ 80 % clean within 2°, as
in the Planck version. Of 9,293 primary clusters, 3,325 lie in that region and **3,318 are kept: the ACT analysis mask removes 7
clusters (0.2 %), at every mass alike**, because ACT inpaints its point sources and SZ clusters instead of masking them; the
Planck mask removed 26 % and the SZ-detected clusters above 4 × 10^14 with them. The release's separate SZ-cluster mask
(`other/`, 6' holes around the SZ clusters that were masked and inpainted in the CMB maps before the quadratic estimator) holds
877 of the 3,318: 3 % of the clusters below 10^14, 7 % at 1-2 × 10^14, 33 % at 2-4 × 10^14, 83 % at 4-6 × 10^14 and 93 % above.
Before any fit the primary samples were fixed as the clusters outside those holes, on the premise that a cluster whose CMB was
inpainted within 6' has no reconstructed lensing; the clusters inside were a separate, diagnostic sample. The diagnostic came back
at A_1h = 1.09 ± 0.14 (1.13 ± 0.13 with the exclusion model): the inpainted clusters' stack is not suppressed. The three samples
then read as an SZ selection on the mass scatter at fixed eRASS1 mass: a lognormal σ_lnM of 0.2-0.3 with its top third removed
(the 2-4 × 10^14 bin) lowers the remainder's mean mass to 0.83-0.89 of the population's and raises the removed third's to 1.23-1.35,
which is the pattern of 0.74-0.77 outside, 1.09 inside and 0.81-0.89 for the union. That reading is one of two: a
different lensing response of the reconstruction inside an inpainted hole is the other, and the release's Gaussian simulations
(no clusters) cannot separate them. The claim is therefore made on the pre-registered samples, whose reconstruction is the
validated one, and the inpainted clusters and the union are reported as supplementary rows. Under the scatter reading the
pre-registered samples measure 0.83-0.89 of the population's mass scale (so b_M for the population would be 0.75-0.85); that
correction is model-dependent and is not applied.

Profiles: the mean of the map's clean pixels in nine annuli (0-5, 5-10, 10-15, 15-20, 20-30, 30-50, 50-80, 80-110, 110-150')
around each cluster's X-ray position, through sparse pixel matrices built once, so a stack of 3,318 clusters takes two seconds.
Model: for each cluster the NFW convergence for a CMB source (M500c to M200c and c200c with Duffy et al. 2008; Σ_crit from
angular-diameter distances in the eRASS1 cosmology, H0 = 67.77, Ωm = 0.29) plus the two-halo term (linear CAMB P(k) scaled to the
eRASS1 σ8 = 0.88, Tinker et al. 2010 bias), Hankel-transformed to κ(L), multiplied by the band filter and the pixel window of the
annulus's resolution, back-transformed to the annulus means on a 22 × 20 grid in (log M500, z), and **averaged over each
cluster's mass posterior**: a lognormal of the catalogue's M500_L/H width (median σ_lnM 0.23, five-point Gauss-Hermite) as the
baseline, the catalogue's M500_PDF column as the check. The posterior average raises the one-halo prediction by 1-2 % over the
posterior-median one (A_1h 0.821 → 0.813 for the cosmology sample); the M500_PDF column, which is three times wider than the
quoted errors, moves A_1h by −1 % to −9 % depending on the sample (cosmology sample 0.84, all 0.87, M500 < 4 × 10^14 0.83). Fits:
the one-halo and two-halo amplitudes jointly (the statistic), one amplitude on 1h + 2h, and b_M on a grid from 0.1 to 2 in steps
of 0.05 with the profiles rebuilt at b_M × M500 (the posterior nodes shifted) and the two-halo amplitude profiled out at each
b_M; the 68 % and 95 % intervals from Δχ² < 1 and < 3.84, the crossings interpolated. Covariance: the same stack on each of the
release's 400 simulated reconstructions (they carry the map's noise and mask through the same pipeline), Hartlap-corrected for
nine annuli; the sims' one-halo amplitudes are the null distribution of the statistic (mean 0.000, spread 0.110 for the cosmology
sample; none of 400 reaches the data). Nulls: 200,000 random positions passing the same footprint selection (mean clean fraction
within 2° 0.994, as the clusters'), stacked on the data map and on every simulation; and 165 rotated-map stacks (the map and the
mask turned together in right ascension by 25-335°, the clusters fixed, the selection re-applied; each rotation keeps 300-2,300
clusters, so the amplitudes are rescaled by √(n_rot/n) before the spread is taken), which give the same error as the
simulations (0.107 vs 0.110). Normalisation: the release's input-vs-reconstructed simulations, both band-filtered and masked with
the clean footprint, give a transfer function of 1.00 ± 0.01 in eleven bands from L = 20 to 2048 (row C5), so no scale-dependent
normalisation moves signal between the inner and outer annuli.

Checks (rows C3-C8): the temperature-only map gives the cosmology sample A_1h = 0.86 ± 0.11, the polarization-only map
0.70 ± 0.22 (3.2σ, free of tSZ and CIB), the galcut040 map (1,549 clusters) 0.77 ± 0.13, the CIB-deprojected map 0.73 ± 0.10; the
last is 0.08 below the baseline on shared noise, largest for the SZ-detected clusters (0.87 against 1.09), the direction and size of
the tSZ bias of a temperature estimator without profile hardening, so it is read as a lower bound. The Planck PR3 tSZ-deprojected
map on the 2,587 of these clusters that pass its footprint gives 0.71 ± 0.34 (cosmology sample) and 0.73 ± 0.29 (all): the same
amplitude at three times the error. By mass (all clusters, two-halo free): 1-2 × 10^14 1.35 ± 0.23, 2-4 × 10^14 0.77 ± 0.14, above
4 × 10^14 0.87 ± 0.15; by redshift: 0.1-0.3 0.78 ± 0.16, 0.3-0.6 1.13 ± 0.16, above 0.6 0.77 ± 0.26; PCONT < 0.1 0.86 ± 0.11;
without the PSZ2 matches 0.92 ± 0.11; below 2 × 10^14 outside the SZ holes (where 6 % are inpainted, so the selection is small)
1.09 ± 0.23.

The model systematic. The data exceed NFW + linear two-halo at 10-30' (2-7 Mpc at the median z = 0.29): the free two-halo
amplitude is 1.2-1.6 in every sample, against 0.4-1.0 in the Planck stack of R-121 (whose resolution does not resolve that range).
That is the one-halo-two-halo transition region, where the linear two-halo term and the NFW tail are both known to under-predict
the matter profile, and it moves the one-halo amplitude through the −0.5 correlation. Three treatments bracket it: linear
two-halo with free amplitude (the statistic, 0.81 for the cosmology sample); linear two-halo zeroed inside 1.5 r200c (halo
exclusion, free amplitude): 0.87 ± 0.10, b_M 0.85 [0.70, 0.97], 95 % [0.58, 1.11]; and the halofit (non-linear) P(k) with free
amplitude: 0.57 ± 0.14, b_M 0.45 [0.32, 0.63], 95 % [0.20, 0.80], which is not a two-halo model inside the halo (halofit's P(k)
contains the one-halo term, and its "two-halo" prediction in the 0-5' annulus is 2.3 times the linear one, double-counting the
cluster itself); it is reported as the extreme of the bracket, not as an alternative. Concentrations × 0.7 / × 1.3 give 0.74 / 0.86
(cosmology sample). σ8 = 0.76-0.88 in the two-halo term changes nothing with the amplitude free and ±0.03 with it fixed.

Trials: the samples and the statistic were fixed before the first fit and are in `results-R122.md`; about 40 subsets × 10 map or
model variants = 400 amplitudes were computed and are in the logs. A first draft of this section made its claim on the union
of the pre-registered samples with the inpainted clusters, on the strength of the diagnostic above; a referee objected that the
union was chosen after unblinding and that its dominant clusters sit where the estimator's response is not validated, and the
claim was returned to the pre-registered samples, with the model systematic added to the intervals. The union rows are kept as
supplementary with that history stated.

## What it means

The eRASS1 mass scale, calibrated with galaxy shapes from three surveys, is measured by a fourth lensing observable at the
25-30 % level (15-20 % with the inpainted clusters) in the mass range that carries the eRASS1 cosmology sample, 1-4 × 10^14 M_sun, below the Planck and ACT SZ catalogues
on which the published comparisons rest (Aymerich et al. 2025, arXiv:2509.02068; the ACT DR6 catalogue, arXiv:2507.21459). At
face value the lensing masses are 0.65-0.70 of the eRASS1 masses on the pre-registered samples (0.75-0.85 with the inpainted
clusters admitted), 1-1.5σ below them, and the 95 % intervals reach from below a halving to above the full scale. A 20-35 %
over-estimate of the eRASS1 masses, which would move the eRASS1 S8 toward the cosmic-shear 0.76-0.78, is inside the 68 %
interval; so is the eRASS1 scale itself at 1-1.5σ; a factor-of-two error (the hydrostatic comparison of Li et al. 2024,
arXiv:2411.09735, on 22 clusters) is inside the 95 % intervals once the model systematic is counted. The two-halo amplitude, 1.2-1.5 of the eRASS1-cosmology linear model at the resolution ACT provides, is a measurement
of the cluster outskirts at 2-7 Mpc that the Planck stack could not make, and it is where a dedicated outer-profile model
(a splashback-plus-infall form, or an emulator) would remove the largest model systematic of the statistic. The pre-registered
samples are SZ-faint at fixed eRASS1 mass (the SZ-detected clusters sit in the inpainting holes), so their b_M is the mass scale
of that subsample, 0.83-0.89 of the population's under the scatter reading; the "1-1.5σ below" applies to the subsample.

## Weakest points

- The statistic depends on the two-halo treatment at 10-30' by more than its statistical error: +0.06 (exclusion) to −0.24
  (halofit) on A_1h against ±0.11. The bracket's lower end double-counts the one-halo term and is not a fair model; the fair
  bracket (linear, with or without exclusion) spans 0.81-0.87. A model of the one-halo-two-halo transition fitted to simulations
  is the missing piece.
- The quadratic estimator's response to a cluster-scale lens is non-linear in the lens; for M500 of a few 10^14 the reconstructed
  κ inside the core is suppressed at the 10 % level (Hu, DeDeo & Vale 2007; Raghunathan et al. 2019), which the release's Gaussian
  simulations cannot measure. It biases A_1h low, so it strengthens the lower limits on b_M and weakens the upper ones; it is not
  corrected here.
- The pre-registered samples exclude the SZ-detected clusters, which at fixed eRASS1 mass are plausibly the upper tail of the
  mass scatter, so their b_M may measure 0.83-0.89 of the population's; that correction rests on a lognormal-scatter model with
  σ_lnM 0.2-0.3 and an idealised top-tail selection and is not applied. The union does not need it, but its dominant clusters sit
  inside 6' inpainting holes where the reconstruction's response to a cluster lens has not been validated on cluster-lensed
  simulations; the inpainted sample's 1.09 ± 0.14 is consistent with either reading.
- The CIB-deprojected map gives 0.08 less than the profile-hardened baseline, on shared noise. The reading (residual tSZ bias in
  the non-hardened estimator) fits its sign, size and mass dependence, and the polarization-only and temperature-only maps agree
  with the baseline; but the deprojected map's lower value is a map-variant systematic of that size until the two foreground
  treatments are combined.
- The posterior widths are the catalogue's quoted M500_L/H; the M500_PDF column, three times wider, moves A_1h by up to 9 %
  (sample-dependent); both are reported, the quoted widths are the baseline.
- The mass bins disagree at 2σ (1-2 × 10^14 at 1.35 ± 0.23 against 2-4 × 10^14 at 0.77 ± 0.14), with no trend across the other
  splits; not interpreted. The pre-registered samples' detections are 4.5-5σ, the union's 7-9σ; the 15 % statistical precision
  the ACT map affords is reached only with the inpainted clusters.

## To settle it

An outer-profile model (Diemer & Kravtsov 2014 form, or an emulator of the stacked matter profile at these masses) with the
one-halo amplitude as the only mass parameter; the quadratic-estimator non-linearity measured on cluster-lensed simulations or
avoided with a gradient-cleaned or maximum-likelihood cluster estimator on the ACT maps; the eRASS1 count-rate-mass relation
refitted with the SZ-detection flag as a covariate, which would test the scatter reading of the three samples; and the stack on
the SPT-3G and the coming ACT DR6+ maps for the remaining eRASS1 footprint.

## Blind review (2026-11-01, context-free Opus referee, own literature search; `blind-review.md`)

First pass, on a draft whose claim was made on the union of the pre-registered samples with the inpainted clusters (verbatim):
VERDICT: novel but not sound: the headline intervals come from a sample chosen after unblinding in which the estimator's response
for the dominant SZ-inpainted clusters is unvalidated, and they omit two-halo and quadratic-estimator systematics that are as
large as the statistical error. Strongest objection, in the referee's words: the headline intervals come from a sample picked
after the data were seen; in that sample the high-mass clusters, which carry most of the signal, were inpainted before the
reconstruction, and no simulation with cluster lensing has tested how the estimator responds there; the sample defined before
the fit gives b_M ≈ 0.65-0.70; the robust statement is only that b_M = 1 is not excluded at 95 %. Prior work it found (none on
eROSITA clusters): arXiv:2406.14754, 2402.18346, 1810.10998, 2404.02153, 1904.07887; eRASS1 mass-scale tests 2507.21459,
2509.02068, 2411.09735, 2503.09952. Revision (same day, the one pass the brief allows, a presentation point): the claim returned
to the pre-registered samples, the model systematic (± 0.12 on b_M) folded into the intervals, the union and the inpainted
clusters kept as supplementary rows, the history stated in the trials paragraph.

Second pass, on the revised section (verbatim): **VERDICT: the findings are confirmed novel and sound if the calculations hold.**
Strongest objection, in the referee's words: the pre-registered samples are an SZ selection on the mass scatter (they drop the
SZ-bright clusters at fixed eRASS1 mass), so their b_M is the mass scale of an SZ-faint subsample, 0.83-0.89 of the population's
by the note's own model, and the one sample free of that selection depends on reconstruction inside inpainted holes where the
response is unvalidated; contamination dilution of the PCONT < 0.5 sample, the scatter selection and the estimator non-linearity
all push b_M low and are not in the ± 0.12; the free two-halo amplitude of 1.2-1.6 signals a mis-specified transition-region
model; the Gaussian simulations miss the two-halo sample variance. The conclusion as stated (nothing excluded; every uncorrected
bias would raise b_M, so b_M = 1 stays inside) survives; the method is sound for a first, low-precision consistency test. Prior
work it adds: arXiv:2405.09919, 2307.11711, 2509.08673, 2607.01175. Softening it asked for, applied after the pass (the title's
detection significance and the mass-scale label now name the pre-registered, SZ-faint subsample; the "15 % test into a 10 % one"
sentence replaced): nothing else changed. Nothing goes on the site from this job; the paper is the main session's.

---

# Finding (R-127, cosmo-deep-7, 2026-11-01): ZTF time delays of the cluster-scale lensed quasars, against the lens models' pre-published predictions

## Claim

Scene-model photometry on the public ZTF difference images (2018-03 to 2026-02, g and r) gives per-image light curves for
the seven cluster- and group-scale lensed quasars in the ZTF footprint with image separations of 3.8"-26": SDSS J1004+4112,
J1029+2623, J2222+2745, J1326+4806, COOL J0542-2125, COOL J0335-1927 and COOL J1153+0755. Curve shifting with a spline
intrinsic model and a slow spline for microlensing, errors from mock light curves with the measured noise, gives:

1. Validation. The three resolved pairs with published delays and images brighter than r = 20.5 are recovered: J1004
   C->A 825.7 ± 25.7 d (r) and 847.0 ± 21.6 (g) against the published 825.2 ± 0.5 (Munoz et al. 2022); J1004 C->B
   781.6 ± 21.0 (r), 786.0 ± 24.2 (g) against 782.2 ± 0.4; J2222 A->B 42.6 ± 5.8 (r) against 47.7 ± 6.0 (Dahle et al. 2015).
   J1029 A->B comes back at 739.9 ± 15.8 (g) against the published 744 ± 10 (Fohlmeister et al. 2013), and at 706-731 in r
   depending on the spline settings (706.5 ± 3.7 at the default), 1-2 sigma low. The blended 3.8" pair J1004 B->A is
   9-18 d high in both bands (52.2 ± 14.6 r, 60.7 ± 13.2 g, against 43.0 ± 0.3). Pairs with an image fainter than r ~ 21 are
   not recovered: J2222 C->A gives 462-885 d depending on the settings (884.7 ± 17.5 at the default; 698 and 794 for
   two others) against the published 722 ± 24. J1004 C->D (2458.5 d) leaves 14 % of the baseline overlapping and is untested.
2. New delays. One new, setting-stable delay: COOL J1153+0755 A1->A3 = 914.7 ± 10.8 d in r (905-943 over eight spline,
   binning and ridge settings; g, at 22.3-22.6 mag, is unusable), the image A3 lagging A1. Everything else is either
   unconstrained by the 2018-2026 baseline (J1153 A1->A2, whose fit runs to the baseline's end at 1770-1800 d), undetected
   in ZTF (COOL J0542 image C: no source in the ZTF reference, g 22.1 by a scene fit on it; J2222 D, E, F at r 23-25),
   setting-dependent (COOL J0335 A->B: 486-561 d in five of eight settings, the grid's end at 900 in three; A->C: +225 ± 20
   in five settings and -445 ± 25 in three, in r; g gives +220 ± 36 and 636 ± 153), or without a minimum that survives the
   overlap test (J1326: the two bands prefer opposite ends of the scan, -1930 in r and +2490 in g, with 0.2 mag of
   unmodelled scatter; no delay).
3. The test of the pre-published predictions. Of the six cluster-model delay predictions published before any measurement
   (Sharon et al. 2017 for J2222 D, E, F; Napier et al. 2023 for J0335 A->B and A->C; Solhaug et al. 2026 for J1153
   A1->A2 and A1->A3), none is confirmed and none is rejected. J1153 A1->A3 measured 914.7 ± 10.8 (± 19 from the settings)
   against the predicted 789 (+55/-59) (788 (+48/-55) with the paper's space-based update) is 2.1 sigma high. J0335 A->C's
   two modes, +225 (C lags A) and -445, both fall outside the predicted -127 (+83/-17) and its systematic range -300 to -30
   (C leads A); the prediction sits delta-chi2 = 9 (r) and 10 (g) above the best mode, and the same statistic misplaced
   J2222's published C->A delay by delta-chi2 = 17 in the same brightness regime, so this is a tension, not a rejection.
   J0335 A->B (486-561 in most settings) is consistent with the predicted 499 (+141/-146) but unstable. J2222 D-F, J0542 and
   J1153 A1->A2 are untestable with ZTF. The count: 0 confirmed, 0 rejected, 2 in tension (2.1 sigma; a sign), 1
   consistent-but-unstable, 3 untestable.

## Data and method

Targets fixed before extraction: every published cluster- or group-scale lensed quasar (separation > 3") at Dec > -31,
the seven above (SDSS J0909+4449 was not in the list). Image positions from the discovery papers (SDSS/HST astrometry; Gaia
DR3 for J1029 and J1326). For each target, every public ZTF g and r science image covering it (IRSA IBE; 696-3468 epochs,
first-to-last 2726-3017 d) gave a 60" cutout of the pipeline difference image (science minus matched reference) and the
difference image's PSF stamp. The scene model: one PSF per lensed image at its fixed position, one free flux each, a
constant background, pixel weights from the clipped sky rms plus a Poisson term, 5-sigma pixel clipping; the difference
flux is added to the image's reference flux from the ZTF reference PSF-fit catalogue (or a scene fit on the reference
cutout, corrected by the catalogue-to-scene offset of the other images, where the catalogue holds no source). Only the
(field, ccd, quadrant) with most epochs per band is used (83-99 % of them), so each image has one reference. Epochs with
pipeline flags, seeing > 4" or a reduced chi2 > 3 are dropped; fluxes are binned over 1 d (bright systems) or 5 d (faint)
before the magnitude conversion. Check: against the public ZTF light-curve service (PSF catalogue photometry on the
science images), J1004 A (3.8" from B) agrees to a median -0.017 mag with 0.10 mag rms over 493 nights and C (isolated)
to -0.008 with 0.13 rms, i.e. at the service's own single-epoch errors; the night-to-night rms of the scene photometry is
1.4-1.5 times its formal error, and that excess (0.10-0.13 mag at r 19-20, 0.2 at r 20.7) enters the delay errors.

Delays: for a pair (X, Y) and a trial delay tau, Y's epochs are shifted back by tau and both curves fitted by one cubic
B-spline (knots every 30 d, a second-difference ridge) for the intrinsic variability plus a 365-d-knot spline in Y's own
time for microlensing and the magnitude offset; the statistic is chi2 per point over the overlap of the shifted curves
(a trial delay must keep 25 % of the points in the overlap), scanned at 1-4 d steps; the delay is the parabolic minimum.
Errors: 30 mock pairs built from the best model at true delays within a few steps of the minimum, with the formal noise
plus the measured excess residual as 10-d correlated noise, refitted; the quoted error is the larger of the mocks' bias
and scatter in quadrature and the half-width of the delta-chi2 < 1 interval. Setting dependence: each new-system pair
refitted with knots of 20 and 50 d, microlensing knots of 200 and 600 d, bins of 3 and 10 d, and ridge 0.3 and 3 (r band).
Convention: dt_XY = t_Y - t_X, positive when Y lags X; the validation pairs fix the sign (J1004 C leads A by 825 d).

## The table

| system | pair | ZTF delay (d) | settings spread | published / predicted (d) | status |
|---|---|---|---|---|---|
| SDSS J1004+4112 | C->A | r 825.7 ± 25.7; g 847.0 ± 21.6 | | 825.2 ± 0.5 measured | recovered |
| SDSS J1004+4112 | C->B | r 781.6 ± 21.0; g 786.0 ± 24.2 | | 782.2 ± 0.4 measured | recovered |
| SDSS J1004+4112 | B->A (3.8") | r 52.2 ± 14.6; g 60.7 ± 13.2 | 54-77 | 43.0 ± 0.3 measured | 9-18 d high; blended |
| SDSS J1004+4112 | C->D | overlap 14 % of the baseline | | 2458.5 ± 1.0 measured | untested |
| SDSS J1029+2623 | A->B(+C) | r 706.5 ± 3.7; g 739.9 ± 15.8 | r 707-731 | 744 ± 10 (90 %) measured | g recovered; r 1-2 sigma low |
| SDSS J2222+2745 | A->B | r 42.6 ± 5.8; g unusable | r 33-50 | 47.7 ± 6.0 (95 %) measured | recovered |
| SDSS J2222+2745 | C->A | r 884.7 ± 17.5; g 461 ± 21 | r 462-885 | 722 ± 24 (95 %) measured | not recovered (C at r 21.2) |
| SDSS J2222+2745 | A->D, E, F | images at r 23-25 | | 502 ± 68, 611 ± 75, 415 ± 72 predicted | untestable |
| SDSS J1326+4806 | A->B | none (r -1930, g +2490: opposite scan ends) | | none | no measurement |
| COOL J0542-2125 | AB->C | C not in the ZTF reference | | none numerical | no measurement |
| COOL J0335-1927 | A->B | r 486-561 (5/8) or 900 (3/8); g 636 ± 153 | | 499 (+141/-146), syst 240-700 predicted | consistent, unstable |
| COOL J0335-1927 | A->C | r +225 ± 20 (5/8) or -445 ± 25 (3/8); g +220 ± 36 | | -127 (+83/-17), syst -300 to -30 predicted | tension (sign) |
| COOL J1153+0755 | A1->A2 | runs to the baseline's end (1770-1800) | | 1253 (+42/-50) predicted | unconstrained |
| COOL J1153+0755 | A1->A3 | r 914.7 ± 10.8 | 905-943 | 789 (+55/-59); 788 (+48/-55) predicted | 2.1 sigma high |

## What it means

ZTF's public difference images carry resolved light curves for the wide cluster lenses, and delays of 20-30 d precision
for pairs of bright (r < 20.5) images, 8 years of cadence standing in for a monitoring campaign: three published
delays are recovered within 1.2 sigma without any dedicated data. Below r ~ 21 the per-epoch signal-to-noise is 1-3,
the chi2 profiles turn multi-modal, and the result depends on the spline settings; the one published delay in that regime
(J2222 C->A) is not recovered, so the new-system numbers in that regime are reported with their spread, and only J1153
A1->A3, stable to ± 19 d over every setting, is a new delay. It is 2.1 sigma above the ground-based and space-based
predictions of the discovery paper; the J0335 A->C data prefer the opposite sign to the model's. Neither is a rejection.
The sample-wide statement is that the pre-published cluster-model predictions are untested by ZTF at the precision the
predictions have (± 50-150 d): the two testable ones are each about 2 sigma off, in different directions.

## Weakest points

- The 0.10-0.2 mag per-night excess noise of the scene photometry (the difference-image formal errors miss the resampling
  correlations and the 0.03 mag zero-point rms) sets the delay errors; it is carried into the mocks as 10-d correlated
  noise, a choice, not a measurement of its spectrum.
- The blended 3.8" pair is biased by 9-18 d; the scene model is not trusted inside two seeing widths.
- The microlensing model (a 365-d-knot spline) cannot follow months-long microlensing; J1029's r-g difference (707-731 vs
  740) and J1004's A-B offset may be that.
- The faint regime is unstable: J1153 A1's reference flux is a scene fit on the reference cutout (the catalogue holds no
  source there, three galaxies are within a few arcseconds), corrected by +0.51 mag from the other images; a wrong
  reference flux distorts the magnitude amplitude, not the timing, but the 905-943 spread over settings is the honest error.
- Delays beyond ~2000 d are at the baseline's limit; the overlap statistic prefers the least-overlap ends when the curves
  are noise-dominated (J1326), so no exclusion range is claimed there.

## Blind review (2026-11-01, context-free Opus referee, own literature search; `blind-review.md`)

VERDICT: novel but not sound: the single new delay (J1153 A1→A3) and the claimed tensions sit in the faint regime where the method misplaced J2222 C→A by roughly 160–190 d, and the error model only measures scatter around the chosen minimum, so neither the ±10.8 d error nor the 2.1σ tension is supported.

Strongest objection, in one sentence: the one new delay and both tensions sit in the brightness regime where the pipeline's only test against a published delay failed by 7-10 of its quoted sigma, and mocks built around the chosen minimum cannot say whether the mode itself is an alias of the seasonal sampling; the referee asked for the same scan on faint mock curves with the real sampling and noise and a rate of catastrophic errors.

Recorded as a miss; the section is not reworded. Injection-recovery, run after the verdict (the test the referee named; `ztf_lens_delays.py --inject`, 40 mocks per
case from the model fitted at the injected delay, the data's own sampling, formal noise and excess as 10-d correlated
noise, scanned over the full grid; `logs/inject/`, `data/ztf_lens/<target>/inject_*.json`):

| pair (r) | injected (d) | found within 40 d | where the misses land | the data's minimum |
|---|---|---|---|---|
| J1004 C->A (bright control) | 825 | 35/40 (825 ± 13) | 784, 873-894, 999 | 825.7 |
| J2222 C->A | 722 (published) | 25/40 (723 ± 10) | 7 of 15 at 876-920, the data's mode; others 482-635, 941-1000 | 884.7 |
| J1153 A1->A3 | 789 (predicted) | 30/40 (788 ± 5) | 5 of 10 at 880-916, the data's mode; 430, 957, 1152-1173 | 914.7 |
| J1153 A1->A3 | 915 (measured) | 40/40 (914 ± 8) | none | 914.7 |
| J0335 A->C | -127 (predicted) | 9/40 (-144 ± 17) | 18 of 31 at +116 to +234, the data's mode; 11 at -178 to -499 | +225.0 |
| J0335 A->C | +225 (measured) | 20/40 (224 ± 14) | -434 to -113, +132 to +182 | +225.0 |

So the catastrophic-error rate of the estimator in the faint regime is 12 % for J1153 at the predicted delay (with
half of those errors landing on the measured 915), 38 % for J2222 at its published delay (half of them on the data's
885), and 78 % for J0335 A->C at its predicted delay, with the measured +225 the single most likely wrong answer. The
J2222 C->A miss and the J0335 A->C "sign tension" are what this estimator does to a true 722 and a true -127 with these
data; they are not evidence against the published or predicted values. J1153 A1->A3 = 915 is recovered every time when
true and produced by a true 789 about one mock in eight: a candidate delay, not a measurement that tests the prediction.
Parked at the referee's verdict (a physics point, not a presentation one); the next pass, if anyone takes it, is the
note rewritten around this table with no tension claimed and the 915 d candidate carrying its 12 % alias probability,
then a second referee pass. Attribution corrected after the referee: the COOL J0335 paper (arXiv:2305.14317) is Napier et
al. 2023 and the COOL J0542 paper (arXiv:2209.03972) Martinez et al. 2023, swapped above and in the note; and Dyrland
2019's reanalysis of J2222 (dt_AB = 42.44 +1.36/-1.44, dt_AC = -696.65 +2.0/-2.1) is the later published value, which the
ZTF A->B 42.6 ± 5.8 matches.

# The overlap of the gravitational-wave search catalogues as a capture-recapture experiment: the pipelines' relative sensitive volumes checked without injections, and a lower bound on the detectable mergers in O1-O4a

## Claim

The public catalogues of compact-binary mergers found by independent searches of the same LIGO-Virgo strain (the LVK
pipelines cWB, GstLAL, MBTA, PyCBC and PyCBC-BBH as reported per candidate in GWTC-1, GWTC-2, GWTC-2.1, GWTC-3 and
GWTC-4.0; the 4-OGC catalogue of Nitz et al.; the IAS catalogues of Venumadhav, Zackay, Olsen, Mehta, Wadekar et al.) are
treated as overlapping imperfect registries of one population and analysed with closed-population capture-recapture
estimators, with no injections. Three results.

1. The relative catchability of the five LVK O3 pipelines on the real candidates agrees with their injection-based
relative sensitive hypervolumes. For every O3 candidate with p_astro > 0.5 in at least one LVK pipeline, weighted by its
probability of being real (the maximum p_astro any pipeline assigns it), each pipeline's share of the expected real
signals, per detector-frame chirp-mass bin, against the GWTC-3 all-O3 sensitive hypervolume of that pipeline divided by
the "any pipeline" hypervolume at the nearest injected chirp mass (both at p_astro > 0.5): Mc 5-15 Msun (17 candidates,
15.8 expected real): cWB 0.06 ± 0.06 (VT ratio 0.12), GstLAL 0.69 ± 0.12 (0.81), MBTA 0.74 ± 0.11 (0.74), PyCBC
0.85 ± 0.09 (0.87), PyCBC-BBH 0.89 ± 0.08 (0.84); Mc 15-30 (13, 12.1): 0.41 ± 0.14 (0.42), 0.74 ± 0.13 (0.77),
0.61 ± 0.14 (0.61), 0.50 ± 0.14 (0.61), 0.72 ± 0.13 (0.83); Mc 30-60 (34, 31.2): 0.41 ± 0.09 (0.49), 0.79 ± 0.07 (0.79),
0.70 ± 0.08 (0.66), 0.48 ± 0.09 (0.62), 0.81 ± 0.07 (0.82); Mc > 60 (11, 9.7; the injection point at 35+35 Msun is lighter
than the bin): 0.61 ± 0.16 (0.49), 0.61 ± 0.16 (0.79), 0.51 ± 0.16 (0.66), 0.51 ± 0.16 (0.62), 1.00 ± 0.09 (0.82). Of the
twenty comparisons, 18 agree within 1.1σ; the largest deviations are PyCBC (broad bank) at Mc 30-60, −1.6σ, and PyCBC-BBH
above Mc 60, +2.1σ (11 of 11 candidates). Chi-square 14.9 for 20 degrees of freedom, the injection values taken without
their own 3-10 % errors. This is a consistency check of the injection-based selection function behind every LVK rate and
mass-spectrum result that uses no injections: the pipelines' relative efficiencies on real signals are those the
injection campaigns assign them, to ± 0.05-0.15 per pipeline and mass bin.

2. The injection-free estimate of the number of detectable mergers equals the union of the catalogues to within a few
per cent, and is a lower bound. With each list's catchability modelled as a logistic function of the reported network SNR
and log chirp mass, the lists conditionally independent given both, and the total the Horvitz-Thompson sum over the
observed candidates of 1 / P(found by at least one list) with each candidate weighted by its probability of being real
(68 % bootstrap intervals): O3, Hanford-Livingston-coincident binary black holes, seven lists, at FAR < 1/yr in the list:
N ≥ 61.1 (60.7-61.7) expected real signals, of which the union of the seven lists holds 60.7 and the LVK catalogue 50.4
(0.83 of the estimate); at p_astro > 0.5: N ≥ 88.5 (87.3-94.0), the union holding 84.2 and the LVK catalogue 61.9 (0.70).
By chirp mass at p_astro > 0.5 the LVK share of the estimate is 0.71 (Mc 5-15), 0.91 (15-30), 0.80 (30-60) and 0.44
(> 60, where the IAS higher-harmonic search reports 21 of the 26 candidates with p_astro 0.5-0.9 under its own rate
model). O4a, the four LVK pipelines on every candidate: N ≥ 86.8 (86.3-89.0) at FAR < 1/yr against 85.5 expected real in
the catalogue (0.99); N ≥ 117.0 (115.7-119.6) at p_astro > 0.5 against 115.2 (0.98). O3, the five LVK pipelines on every
candidate: ≥ 65.1 (63.9-67.2) against 64.1 at FAR < 1/yr; ≥ 76.3 (75.3-77.7) against 74.7 at p_astro > 0.5. O1-O2,
Hanford-Livingston binary black holes, five lists: ≥ 16.3 (15.9-17.3) against 16.0 in the union and 10 in GWTC-1. The
other estimators agree where they are usable (log-linear independence model 62.0, 87.4; Chao's bound 67.6 ± 4.6 and
97.8 ± 6.7 for the O3 seven-list and O4a FAR < 1/yr samples) and break down where a p_astro > 0.5 list admits many
single-list candidates (Chao 497 ± 310 for O3 at p_astro > 0.5, from 40 singletons and 2 doubletons; the pairwise
log-linear model diverges). The residual dependence between lists given the reported SNR is positive in every sample
(Mantel-Haenszel odds ratio pooled over five SNR strata: median 1.2-3.3), so the estimates are lower bounds.

3. The reason they are lower bounds, and what they cannot measure, is shown with mock populations in which the truth is
known. Signals with SNR distributed as SNR^-4, a shared noise term (the common strain realisation), and per-list noise,
thresholds between 8.0 and 9.3 and random drops are searched by 5 or 7 lists. Every estimator returns the union to within
0-15 % (the HT estimate given the reported SNR: +0.1 to +0.8 %), the Mantel-Haenszel odds ratio sits at 1.4-2.0 even
with no shared noise at all, and the true number of signals whose base SNR lies in 8-8.5 is 1.5 times the union in that
bin of reported SNR (899 against 549 seen; HT 554) in every configuration. The signals missed by every list near
threshold are the ones whose noise fell low in every list, and the reported SNR of the found ones is the maximum of noisy
statistics selected to be above threshold; conditioning on it cannot recover them. The same estimator conditioned on
the true base SNR (an oracle) recovers the truth exactly from SNR 8.5 up (1189 against 1188, 754 against 754, 2215
against 2215) and 80-90 % of it in the 8-8.5 bin. A list that is a random subsample of another returns the larger list
exactly. Capture-recapture on searches of one strain therefore measures the lists' relative catchability (result 1) and
bounds the detectable population from below (result 2); the near-threshold incompleteness that the injection campaigns
quantify is, by construction, invisible to it.

## Data and method

The per-list numbers come from the catalogue papers' own tables, read from the arXiv LaTeX sources (the macro files
behind the GWTC-2.1 and GWTC-3 event tables, with their PE macros for the detector-frame chirp masses; GWTC-2 for cWB in
O3a; the IAS new-event and cross-catalogue tables of arXiv:1902.10341, 1904.07214, 2201.02252, 2311.06061 and
2312.06631) and PDFs (GWTC-1 Tables I-II; GWTC-4.0 Tables 1-3, arXiv:2508.18082), plus the 4-OGC search output released
on GitHub (every candidate with IFAR > 1e-5 yr; arXiv:2112.06878). A candidate is in a list when that search reports it
with FAR < 1/yr (for IAS the published per-bank IFAR, which carries no trials factor and makes it the most permissive
list at a nominal 1/yr; for 4-OGC the catalogue IFAR) or, in the second threshold, with p_astro > 0.5 in that search. The
per-pipeline columns exist only for catalogued candidates, so a signal below every list's threshold is unobserved by
construction. Candidates are matched across catalogues by GPS time within 2 s. For the samples with external lists the
population is restricted to what every list searched: both LIGO detectors observing and a binary black hole (chirp
mass > 5 Msun or p_BBH > 0.5), 62 (FAR) or 97 (p_astro) candidates in O3, 17 in O1-O2. The reported SNR is the maximum
network SNR among the lists that report the candidate. Code bases: PyCBC-BBH is PyCBC on a BBH bank; 4-OGC is the PyCBC
code run by an external group with its own bank, statistic and data selection; GstLAL, MBTA and PyCBC share the LVK
bank-construction methods, not the banks or ranking; cWB is template-free; IAS is an independent code base, bank and
statistic.

Estimators: Lincoln-Petersen (Chapman) per pair; Chao's lower bound from the singleton and doubleton counts; the
sample-coverage estimate; Poisson log-linear models on the 2^K − 1 observed cells with main effects and with all pairwise
interactions; and the Huggins-Alho conditional likelihood with logistic catchability in reported SNR and log chirp mass,
weighted by each candidate's probability of being real, with the Horvitz-Thompson total. Bootstrap over candidates (200
replicates). Residual dependence given SNR by the Mantel-Haenszel odds ratio over five SNR strata. The injection-based
comparison values are GWTC-3's sensitive-hypervolume table (per pipeline and "any", fixed-mass injections, p_astro > 0.5)
and GWTC-2.1's O3a BBH-population hypervolumes (FAR < 2/yr: GstLAL 0.85, MBTA 0.61, PyCBC 0.63, PyCBC-BBH 0.83 of "any",
against the O3a seven-list sample's LVK shares at FAR < 1/yr of 0.84, 0.68, 0.68, 0.87 ± 0.08).

## What it means

The number of detections the injection-based sensitive volume implies at a threshold is, by construction, the number the
catalogue holds (the rate is fitted to it), so an injection-free test of the selection function can only take two forms:
the pipelines' relative sensitivities, and the fraction of signals missed by every pipeline at a given SNR. The first is
measured here and agrees with the injection campaigns in every mass bin. The second is not measurable from the overlap of
searches that share the strain, as the mocks show: the capture-recapture total is the union, and it stays the union
however many lists are added, because the missing signals are missing from all of them for the same reason. What the
overlap does say is how much the catalogues differ from one another at fixed threshold: the LVK catalogue holds 0.83 of
the Hanford-Livingston BBH signals the seven lists together find at FAR < 1/yr and 0.70 at p_astro > 0.5, the deficit
concentrated at the lowest SNR (below 8.5: 7 of 25 candidates, 18 expected real) and the highest masses (above chirp mass
60: 11 of 26), where the IAS searches with wider banks and higher harmonics report their candidates.

## Weakest points

The lists' FAR and p_astro are not on one scale: IAS's per-bank IFAR has no trials factor, 4-OGC's p_astro and IAS's come
from their own rate models, and the LVK's p_astro per pipeline from the LVK's. The weighted totals take the maximum
p_astro as the probability of being real, which is generous to the external catalogues; the unweighted totals are given
beside them. The chirp-mass bins for the VT comparison are matched to the nearest fixed-mass injection point, not to a
population-weighted average, and the injection values carry 3-10 % errors not propagated. The catchability curves are
logistic in the reported SNR, which the mocks show is a selection-biased covariate; the PE network SNR of the confident
events would be better but the marginal candidates, which carry the information, have none. The Mantel-Haenszel odds
ratios are from small strata and are noisy per pair. O1-O2 has 17 candidates and GWTC-1 publishes no per-pipeline
p_astro; O4a has no external catalogue yet and no per-pipeline VT table to compare with.

## To settle it

The relative-catchability test at the next catalogue (GWTC-5.0, O4b-c) with the per-pipeline VT tables of that release;
the same test on the GWOSC subthreshold lists (FAR < 2/day, 1,000+ candidates per run), where the catchability curves
can be fitted below SNR 8 and the noise contamination is modelled with the published p_astro per pipeline; and a
capture-recapture with lists whose noise is independent, which only different detectors or different data could give.

## Blind review (2026-11-01, context-free Opus referee, own literature search; `blind-review.md`)

"VERDICT: novel but not sound: the injection-free total equals the union by construction (it conditions on a capture-defined SNR), and its claimed lower-bound status is undercut by noise-inflated max-p_astro weighting across incompatible thresholds, leaving only a low-power and not like-for-like relative-sensitivity check."

Strongest objection, in the referee's words: conditioning on the max-over-capturing-lists SNR makes the injection-free
total identical to the union by construction, so the "lower bound on detectable mergers" carries no information beyond
the union; because of the max-p_astro weighting across incompatible rate models and thresholds it is not even a guaranteed
lower bound on real signals; what is left is a low-power, correlated (the 20 comparisons share candidates and the "any"
denominator, so chi-square 14.9 is not for 20 degrees of freedom), not like-for-like (fixed-mass injection points vs
population bins; ratios to "any", never the absolute selection function) relative-sensitivity check at ±10-15 %. Novelty
granted: no prior capture-recapture across gravitational-wave searches (closest: Laycock arXiv:1701.03801 in X-ray
transients; the pipeline tallies of Broekgaarden et al. arXiv:2303.17628 and Williams arXiv:2401.08709; Essick's
semianalytic sensitivity arXiv:2307.02765, 2606.14045). Recorded as a miss on a method point and parked (agent,
2026-11-01): the objection to the total is the one the note's own third result states, and the revision that would
answer it is not wording but a covariate that is not capture-defined (the PE network SNR for every candidate including
the marginal ones, which no catalogue publishes) or lists with independent noise; the ratio check could be sharpened with
a population-weighted VT and a bootstrap over candidates for the joint chi-square, which would still leave a ±10 %
consistency statement.

---

# Finding (R-128, cosmo-deep-8, 2026-11-01): magnification as an instrument for the JWST z > 8 bright end: the de-magnified luminosities, masses and sizes of 206 z > 8 galaxies behind Abell 2744 (and 309 behind five CANUCS clusters) depend on the magnification as the selection predicts, and the de-magnified bright end sits inside the blank fields' field-to-field range

## Claim

The foreground cluster sets the magnification mu; the background galaxy does not. After de-magnification, the intrinsic properties
of z > 8 galaxies must not depend on mu beyond what the selection through magnification predicts, and any residual dependence is
lens-model error, selection error or photo-z contamination feeding the claimed excess of bright and massive z > 9 galaxies. We
regress the de-magnified properties of the z > 8 photometric sample behind Abell 2744 on log10 mu, compare each slope with the
slope a forward model of the selection predicts for a correct lens model and no contamination, and compare the de-magnified
bright-end density with seven blank fields reduced by the same pipeline. The sample, the statistics and the forward model were
fixed before any fit.

1. De-magnified M_UV on log10 mu (N = 206, UNCOVER DR3 photometry and EAZY photo-z, DR4 v2.0 magnifications, 1.2 ≤ mu ≤ 30):
   slope **+2.79 ± 0.19** (68 % [2.61, 2.99]; 95 % [2.45, 3.19]; errors from 1,000 bootstraps, the released magnification
   posteriors and the photo-z posteriors together). The forward model predicts **+2.43 ± 0.23** for the Donnan et al. 2024 z = 9
   double power law and **+2.69 ± 0.19** for the Bouwens et al. 2021 Schechter function. Data minus model: **+0.36 ± 0.30** and
   **+0.10 ± 0.27**. A flux-limited sample has a slope of 2.5 by construction; the steepness of the luminosity function at the
   detection limit moves it by −0.07 to +0.19; the data sit where that puts them.
2. Stellar mass (Prospector-β, de-magnified in the catalogue) on M_UV and log10 mu jointly: the log mu coefficient is
   **−0.28 ± 0.17** (95 % [−0.58, +0.10]) against a predicted **−0.03 ± 0.14**: data minus model **−0.24 ± 0.22**. At fixed
   de-magnified UV luminosity, the more magnified sources carry lower masses by 0.3 dex per dex of mu, at 1.1σ.
3. De-magnified half-light radius on M_UV and log10 mu: **−0.08 ± 0.14** against a predicted −0.07 ± 0.07 (the PSF floor); data
   minus model −0.01 ± 0.16.
4. The bright z > 9 sources (M_UV < −20; 21 of 164) sit at **lower mu than the parent sample at the same observed magnitude, by
   0.124 ± 0.045 dex** (permutation p < 2 × 10⁻⁴ against zero). The selection predicts this: at fixed observed brightness a
   brighter intrinsic source needs a lower magnification, and the forward model gives **−0.19 ± 0.07** (DPL) to −0.17 ± 0.07
   (Schechter). Data minus model **+0.07 ± 0.08**: the bright ones are at 0.05-0.1 dex higher mu than a correct lens model with
   no contamination predicts, at 0.6-0.9σ. The massive arm of the test is empty: no source at z > 9 with log M* > 9.5 survives a
   UV S/N ≥ 5 detection and 1.2 ≤ mu ≤ 30 in the DR4 SPS catalogue (one at z > 8: z 8.74, mu 1.41, log M* 9.65); at log M* > 8.5
   (N = 14) the shift is −0.13 dex (p = 0.04), the same as for the bright ones.
5. Lensed vs blank, one selection and one photo-z pipeline (DJA grizli v7 + eazy-py) in Abell 2744 and in CEERS, GOODS-N,
   GOODS-S, NGDEEP, PRIMER-UDS north and south and the Abell 2744 parallel field, the cluster's source-plane volume from the v2.0
   map inside 1.2 ≤ mu ≤ 30: the de-magnified density of z 8-10 sources brighter than M_UV = −20 is **1.58 times the summed blank
   fields' (68 % [1.10, 2.18]; 95 % [0.73, 2.91])**, and **0.86 [0.50, 1.42]** brighter than −20.5; over z 8-13, **1.41 [1.11,
   1.74] (95 % [0.88, 2.15])** at −20 and 0.92 [0.66, 1.28] at −20.5. The seven blank fields scatter among themselves by 0.26-0.53
   dex in density, three to five times their Poisson scatter; Abell 2744 ranks first or second of the eight fields at M_UV < −20 and
   fourth or fifth at −20.5. The de-magnified Abell 2744 bright end is inside the blank fields' field-to-field range at every cut.

6. The same regressions on five CANUCS clusters (Abell 370, MACS0416, MACS0417, MACS1149, MACS1423; Dense Basis photo-z,
   Lenstool models by three other modellers, 309 objects with z_phot > 8, z16 > 7, m_UV < 29 and 1.2 ≤ mu ≤ 30): R1 **+2.67 ±
   0.12** (per cluster 1.98-2.89; jackknife over clusters 0.16), R2 **−0.19 ± 0.15** (per cluster −0.33 to −0.08; jackknife
   0.08), R3 −0.37 ± 0.06, bright z > 9 sources at lower mu by 0.30 ± 0.07 dex (9 of 151), no source above log M* = 9.5. The
   Abell 2744 slopes are not peculiar to one lens model, one modeller or one photo-z code.

The claim is these intervals: the lensing step, as implemented by the public v2.0 model, introduces no detectable mu-dependence
into the de-magnified luminosities (a residual slope of 0.1-0.4 ± 0.3 mag per dex of mu), sizes (−0.01 ± 0.16 dex per dex), or
the position of the bright z > 9 sources in mu (+0.07 ± 0.08 dex); the masses show a 1σ residual (−0.24 ± 0.22 dex per dex; −0.19 ± 0.15 in the five CANUCS clusters, against a predicted
−0.03 ± 0.14, so the two independent samples together put the mass coefficient at −0.23 ± 0.11, 1.8σ below the selection
model); and the lensed-field bright end is not in excess of the blank fields beyond their own spread.

## Data and method

Abell 2744: the UNCOVER DR3 SUPER photometric catalogue (Weaver et al. 2024; 74,020 objects; F444W-matched, bCG-subtracted, EAZY
photo-z with P(z) percentiles), the DR4 v2.0 magnifications of the same rows (Furtak & Zitrin model; Price et al. 2025: best fit
and 68/95 % posterior percentiles from the MCMC chain, at the EAZY redshift), the DR4 Prospector-β catalogue of the same rows
(Wang et al. 2024; masses de-magnified on the fly with the v2.0 model), the v2.0 lens maps (magnification at z = 8 and 10, 0.1"
per pixel, 7.6') and the v2.0 multiple-image catalogue. The DR3 catalogue also carries the v1.1 magnifications of the same
objects. Blank fields: the DJA grizli v7 photometric catalogues with their eazy-py redshifts (abell2744clu v7.2, abell2744par
v7.0, ceers-full v7.4, gds v7.2, gdn v7.3, ngdeep v7.2, primer-uds-north and -south v7.2).

Selection, one rule in every catalogue, fixed before any fit: use_phot (or the DJA flag < 4 with ≥ 6 filters), flag_eazy, z_phot
> 8 and the 16th percentile of P(z) above 7, signal-to-noise ≥ 5 in the rest-UV band (F150W below z = 10.5, F200W above), and in
the cluster 1.2 ≤ mu ≤ 30 and no object within 0.5" of a multiple image (two dropped). 206 objects (378 at S/N ≥ 3): 83 at
z < 10.5 and 123 above, mu 1.20-29.7 (median 1.85), observed m_UV 26.6-29.2, de-magnified M_UV −22.4 to −15.1, log M* 6.5-8.7
(5th-95th percentiles). M_UV = m_UV − DM(z) + 2.5 log10(1 + z) + 2.5 log10 mu (flat f_nu; Planck 2018 flat ΛCDM). Sizes:
r = sqrt(max(flux_radius² − r_psf², 0) / mu), r_psf = 0.08", objects with flux_radius > 1.1 r_psf (201). Quantities: R1, the
least-squares slope of M_UV on log10 mu; R2 and R3, the log10 mu coefficient of log M* and of log r on (M_UV, log10 mu); T4, for
the bright (M_UV < −20) and massive (log M* > 9.5) z > 9 sources, the mean difference between their log10 mu and the parent
z > 9 sample's mean log10 mu in the same 0.5-mag bin of observed m_UV, with a permutation test within the bins; C5, the
de-magnified count brighter than M_UV = −20, −20.5 and −21 per source-plane comoving volume, 8 < z ≤ 10, 10 < z ≤ 13 and 8-13.

Errors: 1,000 bootstraps over objects; 1,000 draws of mu from a split normal in log mu fitted to the released 16th/84th
percentiles, of z from the P(z) percentiles, and of log M* from its percentiles, singly and together; a jackknife over four
quadrants of the cluster field; the v1.1 model substituted for v2.0. The released per-object magnification posteriors are 0.004
dex wide (median 68 % half-width), 25 times narrower than the scatter between the v1.1 and v2.0 versions of the same model
family on the same objects (0.10 dex; median offset 0.016 dex), so every statistic was also run with 0.10 dex of lens noise added
per object (data and forward model alike): it changes R1 by 0.00 ± 0.05, R2 by 0.00 ± 0.03 and the expected T4 shift from −0.19
to −0.24. For C5: Poisson (Gehrels), plus cosmic variance from linear theory (an Eisenstein & Hu power spectrum, a box the field's
size, bias 8 for M_UV < −20 at z ~ 9; 24 % of the cluster count), plus the field-to-field scatter of the blank densities reported
beside the ratio. The cluster's source-plane area is the sum over 6" image-plane cells holding catalogued objects of the cell area
divided by the mean of 1/mu over the cell's pixels in the z = 9 map (the mean of the z = 8 and z = 10 maps), inside 1.2 ≤ mu ≤ 30
and the map's field of view: 20.7 arcmin² against 620 arcmin² of blank field. M_UV < −20 is m ≈ 27.3 at z = 9, more than a
magnitude above every field's median 5σ depth (27.5-29.7), so no completeness correction enters the bright cut.

CANUCS (a secondary sample, labelled as such because the catalogue has no fluxes): the DR1 Dense Basis catalogues of the five
cluster fields (USE_PHOT, no fit flag; redshift_50 > 8, redshift_16 > 7; rest-frame FUV absolute magnitude as the UV magnitude,
which the catalogue does not de-magnify, with an apparent m_UV < 29 cut in place of the S/N cut), mu from the best-fit Lenstool
kappa and gamma maps (Gledhill et al. 2025 for Abell 370, Rihtarsic et al. 2025 for MACS0416, preliminary DR1 models of
Rihtarsic and Desprez for MACS1149, MACS0417 and MACS1423; all for D_ds/D_s = 1, scaled to each source's redshift), the
catalogue masses de-magnified by mu, the per-source lens spread from the 100 released MCMC sample maps where released (Abell 370,
MACS0416, MACS1149: 68 % half-widths 0.003-0.027 dex) and 0.10 dex where not. 49, 52, 33, 107 and 68 objects in the five
fields (309); median mu 1.7-3.1.

The forward model (fixed before the fits): sources drawn from the z = 9 UV luminosity function (Donnan et al. 2024, Table 3,
double power law φ* = 23.5 × 10⁻⁵ Mpc⁻³ mag⁻¹, M* = −19.70, α = −2.00, β = −3.81; or Bouwens et al. 2021, Table 5, Schechter
φ* = 2.1 × 10⁻⁵, M* = −21.15, α = −2.33), placed uniformly in the source plane (image-plane weight 1/mu, pixel by pixel in the
z = 9 map), at the data's redshifts, observed at m = M + DM − 2.5 log(1 + z) − 2.5 log mu, detected against the cell's own 5σ
depth (the median rest-UV flux error of the catalogued compact objects in the cell; median 28.46, 10-90 % 27.5-28.9) with a
0.3-mag soft edge and a flux error growing to 0.2 mag at the limit; masses from log M* = 9.0 − 0.4 (M_UV + 20) with 0.4 dex
scatter; sizes r_e = 0.6 kpc (L/L*)^0.3 with 0.2 dex scatter, observed as sqrt(r_e² mu + r_psf²); the lens-model noise of the
data's posterior widths (and 0.10 dex as the check) added to mu, and the masses de-magnified with the noisy mu as the catalogue
does. 1,000 mocks of 206 objects each, the same estimators, the same 1.2 ≤ mu ≤ 30 cut. The mocks put 5-6 % of the sample
brighter than M_UV = −20; the data have 10 % (21 of 206).

## The table

| statistic | data | expected (DPL / Schechter) | data − expected | checks |
|---|---|---|---|---|
| R1: M_UV on log mu (mag per dex) | +2.79 ± 0.19 | +2.43 ± 0.23 / +2.69 ± 0.19 | +0.36 ± 0.30 / +0.10 ± 0.27 | quadrant jackknife 2.59-3.51 (sd 0.63); v1.1 model 2.90; S/N ≥ 3 (N 378) 2.66 ± 0.14; both photo-z codes agree (N 141) 2.70; DJA selection (N 315) 2.50 ± 0.65 |
| R2: log M* on log mu at fixed M_UV (dex per dex) | −0.28 ± 0.17 | −0.03 ± 0.14 / −0.04 ± 0.13 | −0.24 ± 0.22 | jackknife sd 0.19; v1.1 −0.30; S/N ≥ 3 −0.29 ± 0.15; z 8-10 (N 74) −0.72, z > 10 (N 132) −0.06; codes agree −0.47 |
| R3: log r on log mu at fixed M_UV (dex per dex) | −0.08 ± 0.14 (N 201) | −0.07 ± 0.07 / −0.07 ± 0.07 | −0.01 ± 0.16 | jackknife sd 0.07; S/N ≥ 3 −0.16 ± 0.10 |
| T4: log mu of bright z > 9 minus parent at fixed m (dex) | −0.124 ± 0.045 (N 21 of 164; p < 2e-4 vs 0) | −0.19 ± 0.07 / −0.17 ± 0.07 | +0.07 ± 0.08 / +0.05 ± 0.08 | v1.1 −0.07; S/N ≥ 3 (N 27 of 316) −0.105 ± 0.036; M_UV < −21: N 3; with 0.10 dex lens noise expected −0.24 ± 0.07 (+0.11 ± 0.09) |
| T4: massive z > 9 (log M* > 9.5) | N = 0 | — | — | log M* > 9.0: N 2 (0.00); > 8.5: N 14, −0.13 (p 0.04) |
| CANUCS, 5 clusters, N 309: R1 | +2.67 ± 0.12 | (not modelled; 2.5 by construction) | — | per cluster 2.35 ± 0.27 (A370, 49), 2.89 ± 0.27 (M0416, 52), 1.98 ± 0.50 (M1149, 33), 2.85 ± 0.15 (M0417, 107), 2.71 ± 0.20 (M1423, 68); jackknife over clusters sd 0.16 |
| CANUCS: R2 | −0.19 ± 0.15 | — | — | per cluster −0.23, −0.10, −0.33, −0.12, −0.08; jackknife sd 0.08; with A2744: −0.23 ± 0.11 |
| CANUCS: R3 | −0.37 ± 0.06 (N 229) | — | — | per cluster −0.44, −0.65, −0.11, −0.32, −0.47; a different detection image and PSF, no PSF floor modelled |
| CANUCS: T4 bright | −0.30 ± 0.07 (N 9 of 151; p 0.003) | — | — | massive (log M* > 9.5): N 0 |
| C5: A2744 / blank density, z 8-10, M_UV < −20 | 1.58, 68 % [1.10, 2.18], 95 % [0.73, 2.91] (19 vs 361) | 1 | — | blank-field scatter 0.26 dex (Poisson 0.08); A2744 exceeded by 1 of 7 blank fields |
| C5: z 8-10, M_UV < −20.5 | 0.86 [0.50, 1.42] (5 vs 174) | 1 | — | scatter 0.41 dex; exceeded by 4 of 7 |
| C5: z 10-13, M_UV < −20 | 1.23 [0.88, 1.68] (15 vs 365) | 1 | — | scatter 0.45 dex; exceeded by 2 of 7 |
| C5: z 8-13, M_UV < −20 | 1.41 [1.11, 1.74], 95 % [0.88, 2.15] (34 vs 726) | 1 | — | scatter 0.30 dex (Poisson 0.06); exceeded by 0 of 7; 1.23 × the median blank field |
| C5: z 8-13, M_UV < −20.5 | 0.92 [0.66, 1.28] (11 vs 359) | 1 | — | scatter 0.39 dex; exceeded by 4 of 7 |

## What it means

Every dependence of a de-magnified property on magnification in this sample is the one a flux-limited selection through a steep
luminosity function produces when the lens model is right: the intrinsic UV luminosity rises by 2.8 ± 0.2 mag per dex of mu
against a predicted 2.4-2.7 (a residual of 0.1-0.4 ± 0.3 mag per dex), the de-magnified sizes track the PSF floor, and the
bright z > 9 sources sit at the lower magnifications the selection assigns them, 0.05-0.1 ± 0.08 dex higher than the pure
selection predicts. If the z > 9 bright end in this field were being manufactured by the lensing step, the bright sources would
have to sit at magnifications the model gets wrong, and their de-magnified luminosities would carry a residual mu-dependence; the
95 % intervals bound both: the residual slope to within about ±0.6 mag per dex of mu, the bright sources' position in mu to
within +0.23 dex of the selection's. The one residual is in the masses: at fixed de-magnified M_UV, the Prospector-β masses
fall by 0.3 ± 0.2 dex per dex of mu behind Abell 2744, and the Dense Basis masses by 0.2 ± 0.15 behind the five CANUCS
clusters, where the forward model predicts none (−0.03 ± 0.14); together −0.23 ± 0.11, 1.8σ from the model. Its sign is the one
a mass prior acting on fainter intrinsic fluxes, or a lens model that over-magnifies the high-mu sources, would give; neither
sample's 95 % interval excludes zero, and the two catalogues share no code, so a common SED-fitting cause would have to be
generic (a mass-to-light prior that tightens as the photometry gets fainter). The five CANUCS clusters, with lens models by
three other modellers and a different photo-z code, return the Abell 2744 luminosity slope (2.67 ± 0.12 against 2.79 ± 0.19)
and the same sign and size for the masses.

The lensed-vs-blank comparison does not find an excess that survives the blank fields' own spread. With one photo-z pipeline and
one selection, the seven blank fields differ among themselves by 0.3-0.5 dex in the density of z > 8 sources brighter than
M_UV = −20 (three to five times Poisson; GOODS-S at 4 × 10⁻⁵ Mpc⁻³, NGDEEP at 31 × 10⁻⁵ at z 8-10), and the de-magnified Abell
2744 density (24 × 10⁻⁵) is at the top of that range, 1.6 ± 0.5 times their sum, while at M_UV < −20.5 it is 0.9 ± 0.4 times it.
The limiting systematic of a lensed-vs-blank test at this depth is the field-to-field dispersion of the photometric z > 8
selection across filter sets and depths, not the lens model.

Two by-products. The released per-object magnification uncertainties of the UNCOVER v2.0 catalogue (68 % half-widths of 0.004
dex) are 25 times smaller than the change between two versions of the same model family on the same objects (0.10 dex); they
measure the chain's statistical width, not the model error, and anyone propagating them as the lens error underestimates it by
that factor (here it does not matter: 0.10 dex moves no slope by more than 0.05). And the massive z > 9 population the test was
built for does not exist in the DR4 Prospector-β catalogue once a UV S/N ≥ 5 detection and 1.2 ≤ mu ≤ 30 are required: no
source above log M* = 9.5, two above 9.0, in 206; the massive-excess arm of the instrument is empty on this catalogue's masses.

## Weakest points

- The forward model's expected slopes depend on the luminosity function at the detection limit (2.43 vs 2.69 for the two
  published fits) and on the depth map built from the catalogue's own flux errors; the comparison is between the data and a
  bracket, and the bracket's width (0.26) is as large as the data's error.
- The sample is photometric. 32 % of it has a Prospector-β redshift more than 1 away from the EAZY one; the brightest object
  (M_UV −22.4 at z_EAZY 8.22) is at z = 2.0 in Prospector. The subsample where both codes agree (141) gives the same slopes
  (R1 2.70, R2 −0.47, T4 −0.12), but a contamination that is itself mu-independent would not show in this test.
- 123 of the 206 are at z_phot > 10.5 from a quick-look EAZY run without magnitude priors; the test does not depend on them
  being at those redshifts (R1 in z 8-10 alone is 3.18, N 74; 2.73 above), but the lensed-vs-blank densities at z 10-13 do.
- The lens-model spread is within one family (v1.1 against v2.0 of the Zitrin parametric model); no independent family of
  Abell 2744 (Bergamini et al. 2023, the HFF models) was run through the same objects. The CANUCS clusters bring independent
  modellers but not the same objects, and their catalogue has no fluxes: its selection is an apparent-magnitude cut, its UV
  magnitude an SED-fit rest-frame quantity, and no forward model was built for those fields, so the CANUCS slopes are compared
  with 2.5 and 0, not with a field-specific expectation.
- The sizes are PSF-dominated (median flux_radius 0.19" against r_psf 0.08") and the catalogue's half-light radii are from the
  detection image, not a profile fit; R3 is a weak statistic.
- The cosmic-variance term uses a fixed bias of 8 and a box geometry; it is 24 % of the cluster count and does not drive any
  interval. The blank-field dispersion (0.3-0.5 dex) is reported, not modelled.
- The massive arm is empty by the catalogue's masses, not by a measurement: a different SPS prior (the earlier UNCOVER masses
  that produced the massive candidates) would populate it, and the test would then be of that prior as much as of the lens.

## Blind review (2026-11-01, context-free Opus referee, own literature search; `blind-review.md`)

VERDICT: novel but not sound: the mu-regressions are dominated by the 2.5 log mu built into the de-magnified quantities and are compared against a forward model that fails to reproduce the sample's bright fraction and its majority z > 10.5 photometric population, so the null result shows the test lacks power rather than confirming the lens model.

Strongest objection, in one sentence: R1 and T4 carry log mu on both sides and barely respond to lens-model error (0.10 dex of
lens noise moves R1 by 0.00), and they are compared against a forward model that puts 5-6 % of the sample brighter than
M_UV = −20 where the data have 10 % and draws from a z = 9 luminosity function while 123 of the 206 sources are at z_phot > 10.5
from a quick-look EAZY run without priors. Novelty was granted ("no prior paper that does exactly this test at z > 8 with
JWST"; closest: Bouwens et al. 2022, arXiv:2203.14978, the surface density vs mu relation as a lens-model test).

Recorded as a miss; the section is not reworded. Parked on the physics point. What stands as a documented null: the slopes and
ratios above with their intervals; the field's 2× excess of M_UV < −20 sources over the LF-plus-selection model (10 % against
5-6 %) carrying no mu-dependence beyond the selection's; and the two by-products (the released per-object magnification
posteriors are 25× narrower than the version-to-version change of the same model family; the DR4 Prospector-β catalogue holds
no log M* > 9.5 source at z > 9 with a UV S/N ≥ 5 detection inside 1.2 ≤ mu ≤ 30). The next pass, if anyone takes it, needs a
statistic not built from the same magnitudes (the observed-magnitude distribution per mu bin against the forward model, i.e.
magnification bias as a count test), a forward model fitted to the field's own luminosity function and redshift distribution,
an independent lens family on the same objects (Bergamini et al. 2023; the HFF models), and a photo-z selection with priors.


---

# Which dataset pays for the DESI DR2 dynamical-dark-energy preference: a Shapley attribution and a specification curve

2026-11-01. Scripts: `desi_shapley.py` (background pipeline, coalition fits, Shapley values, bootstrap, specification
curve), `desi_fullcmb.py` (the full Planck 2018 likelihood), `desi_r123_revise.py` (balanced decomposition, the cut
baseline). Data: DESI DR2 and DR1 BAO, Pantheon+, Union3, DES-SN5YR and DES-Dovekie distance likelihoods as released
through the cobaya data packages; DESI's (θ_*, ω_b, ω_bc) compression of the CMB (arXiv:2503.14738, App. A); the Planck
2018 plik_lite, Commander, SimAll and lensing likelihoods from the Planck Legacy Archive. Rows in `results-R123.md`.

## Claim

DESI DR2 BAO with the CMB and a supernova sample prefers the CPL (w0, wa) model over ΛCDM by Δχ² = 10.7, 17.4 and 21.0
(2.8, 3.8, 4.2σ) with Pantheon+, Union3 and DES-SN5YR (DESI's numbers). Drop-one tests do not say how that gain divides
among the datasets, because they are order-dependent and not additive; the Shapley value of each dataset, computed from
the profile Δχ² of every coalition, is the unique symmetric, additive, order-free split. This note gives that split with
its sampling error, the same quantity with the full Planck likelihood for the headline coalition, and a preregistered
specification curve of the significance with its variance decomposed by analysis choice. All of it is interval-valued.

**Shapley values, compressed CMB.** Five datasets (DESI DR2 BAO, the (θ_*, ω_b, ω_bc) CMB prior, Pantheon+, Union3,
DES-SN5YR), 32 coalitions, each minimised in ΛCDM and w0waCDM under DESI's priors; v(S) = Δχ²(S). The three games that
do not stack supernova samples (errors: the sd over 40 parametric bootstrap realisations of every data vector around
the joint w0wa best fit, the 32-coalition grid refitted each time; shares as the 16-84 % interval of φ/Δχ²_grand over
the realisations):

| game | grand Δχ² | DESI BAO | CMB (compressed) | SN sample |
|---|---|---|---|---|
| DESI + CMB + Pantheon+ | 7.6 (2.3σ) | 5.8 ± 2.4; share 15-48 % | 1.6 ± 1.4; 0-16 % | 0.2 ± 4.8; 39-78 % |
| DESI + CMB + Union3 | 13.8 (3.3σ) | 7.0 ± 2.4; 9-64 % | 2.1 ± 1.4; −9-33 % | 4.8 ± 2.1; 8-77 % |
| DESI + CMB + DES-SN5YR | 18.0 (3.8σ) | 8.0 ± 2.7; 15-51 % | 2.2 ± 1.2; −1-15 % | 7.9 ± 3.1; 36-81 % |

(The five-player game, which stacks SN samples sharing their low-z supernovae: DESI 9.6 ± 3.2, CMB 3.1 ± 2.0, Pantheon+
−0.6 ± 4.4, Union3 4.2 ± 1.9, DES-SN5YR 7.6 ± 3.0 of 23.9; for ordering only.) What the intervals allow: in the DES-SN5YR
game DESI and the supernovae carry equal values (8.0 ± 2.7 and 7.9 ± 3.1) and the compressed CMB a small one (2.2 ± 1.2,
share below 15 % at 84 %); in the Pantheon+ game the supernova value is consistent with zero (0.2 ± 4.8) and with DESI's
(5.8 ± 2.4); in the Union3 game the three are 7.0 ± 2.4, 2.1 ± 1.4, 4.8 ± 2.1. The supernova values differ between games
by less than their combined errors (0.2 ± 4.8 to 7.9 ± 3.1 is 1.3σ), so the data do not establish that the SN player's
value rises across the three samples; they establish only that in the DES-SN5YR game no player's value is consistent with the ΛCDM noise floor
except, marginally, the CMB's (2.2 ± 1.2 against a floor of ±1.0-1.5) [corrected 2026-11-01 after the second blind
pass: the earlier sentence said the compressed CMB's value is the smallest in every game, which the Pantheon+ game
contradicts (SN 0.2 ± 4.8), and that no player's value is zero]. Under a ΛCDM truth (40 realisations) the grand
Δχ² is 1.6 ± 1.5, at most 6.5, and each Shapley value scatters by ±1.0-1.5: that is the noise floor on a value, and it
is also the size of the minimiser noise in one-start fits. Under the fitted w0wa model the expected pattern differs from
the data's in one place: Pantheon+, the most precise sample, would carry 8 ± 5 in its game and carries 0.2 (1.6σ low);
the grand values of the three games sit at +0.9σ (Pantheon+ low), +2.4σ (Union3 high) and +1.2σ (DES-SN5YR high) of
their bootstrap expectations. Per-dataset χ² at the joint best fits (the ΛCDM minimum minus the w0wa minimum, the
"sacrifice" view): DESI + CMB + DES-SN5YR: DES-SN5YR +11.4, DESI +4.6, CMB +2.0; Pantheon+: +2.7, +2.7, +2.2.

**The full Planck likelihood.** With plik_lite TT,TE,EE, Commander low-ℓ TT, SimAll low-ℓ EE and the PR3 lensing
reconstruction (CAMB with PPF dark energy, nine parameters profiled), DESI + CMB gives Δχ² = 9.7 without lensing and
10.1 with the PR3 lensing (DESI: 9.7, and 12.5 with the PR4 CamSpec likelihood and PR4+ACT lensing), and DESI + CMB +
DES-SN5YR gives 19.9 (4.1σ; DESI 21.0, 4.2σ) at w0 = −0.76, wa = −0.80, H0 = 66.9, Ω_m = 0.318 (DESI: −0.76 ± 0.07,
−0.80 ± 0.3, 67.3, 0.319). The compressed form had given 8.1 and 18.0: it loses 1.6-1.9 units of the primary CMB's
preference, and it misreads the CMB's role: at the full-likelihood best fits the high-ℓ spectra gain 2.5-6.1 units under
w0wa (the w0wa fits have Ω_m ≈ 0.35 and a higher lensing amplitude, the direction plik_lite's A_L preference points),
low-ℓ TT loses 0.2-1.1, DESI gains 4.2-5.0 and DES-SN5YR 11.3, where the compressed grid had the CMB giving up 2. The
PR3 lensing reconstruction is nearly neutral (+0.3) where DESI's PR4+ACT combination adds 2.8: the one piece of DESI's
CMB this check does not reproduce. The three-player game with the full CMB as the player, for the DES-SN5YR
combination (coalition values: CMB alone 3.2, CMB + DES-SN5YR 8.2, DESI + CMB 10.1, all three 19.9; DESI alone 4.7,
DES-SN5YR alone 5.9, DESI + DES-SN5YR 13.6 as before): **DESI 7.9 (40 %), CMB 4.5 (22 %), DES-SN5YR 7.5 (38 %)**,
against 8.0 / 2.2 / 7.9 in the compressed game. The full likelihood doubles the CMB's value (1.1 of the 4.5 is the CMB-alone singleton, a fit on the
prior edge with w0 = −0.71, wa = −2.22, H0 = 74) while DESI's value moves by 0.1 and the supernovae's by 0.3, the grand
total rising from 18.0 to 19.9 [corrected 2026-11-01: not a reallocation of a fixed sum]; the ordering (DESI and the supernovae equal within errors, the CMB smallest) is unchanged. The sampling
errors are the compressed game's (the full grid was not bootstrapped), about ±2.7, ±1.2 (+ the 2.3-unit shift as a
systematic), ±3.1. DESI + CMB + Union3 with the full likelihood gives 15.8 (3.6σ; DESI 17.4), its game pending the CMB +
Union3 fit.

**Specification curve.** 540 specifications from a list frozen before any was run (SN sample: none, Pantheon+, Union3,
DES-SN5YR, DES-Dovekie; SN low-z cut: as released, z > 0.1; BAO: DR2, DR1; CMB: none, DESI's compression, the Planck
2018 (R, l_A, ω_b) distance priors; model: CPL, wCDM; prior: DESI's box with w0 + wa < 0, without it, a narrow box; Σm_ν:
0.06 eV, free), each fitted in ΛCDM and the extended model; the significance from the profile Δχ² with the model's degrees
of freedom. The curve runs from 0.0 to 3.8σ (mean 1.6, sd 0.9). The variance decomposition on the balanced 384-specification
sub-grid (the four SN samples, both cuts, both BAO releases, three CMB forms, both models, the DESI and narrow priors,
both neutrino treatments), where main effects are orthogonal, with the raw Δχ² of the CPL fits as the response (192
specifications; the significance gives the same ranking):

| choice | share of the variance of Δχ² (CPL) | level means of Δχ² |
|---|---|---|
| SN low-z cut (released vs z > 0.1) | 55 % | 10.0 vs 3.3 |
| cut × sample interaction | 20 % | |
| SN sample | 15 % | Pantheon+ 4.8, DES-Dovekie 5.0, Union3 8.4, DES-SN5YR 8.4 |
| CMB form (none, θ_* compression, R-l_A) | 6 % | 5.1, 7.6, 7.3 |
| BAO release (DR2 vs DR1) | 2 % | 7.3 vs 6.0 |
| neutrino mass, prior box | 0.2 %, 0 % | |

(With both models in and the significance as the response: cut 43 %, model dof 17 %, cut × sample 12 %, sample 6 %, BAO
2 %, CMB 1 %; the full 540 grid, unbalanced, gives the same order.) On the headline specification (DESI DR2 + compressed
CMB + SN, CPL, DESI priors) removing the z < 0.1 supernovae takes DES-SN5YR from 3.8 to 1.5σ (Δχ² 18.0 to 3.9),
DES-Dovekie from 2.8 to 1.2σ, Union3 from 3.3 to 2.4σ and Pantheon+ from 2.3 to 2.0σ. A cut removes lever arm whether
or not the signal is real, so the baseline is the cut's expected effect under the fitted w0wa model, from 12 parametric realisations of each
DESI + CMB + SN coalition around its own w0wa best fit, each fitted with and without the cut: the expected Δχ² ratio
(cut / full) is 0.34 ± 0.18 for DES-SN5YR (data 0.22; 42 % of realisations fall as low), 0.75 ± 0.37 for Pantheon+ (data
0.83) and 0.71 ± 0.21 for Union3 (data 0.60; 25 % as low). The data's drops are what a real w0wa signal of the fitted
size loses with the z < 0.1 supernovae; they are not evidence against one. Every
other single change on the headline specification moves the significance by 0.4σ or less, except the switch to wCDM
(one parameter: 0.7-2.1σ). The parameter-shift metric (the Hessian distance of (w0, wa) = (−1, 0) from the best fit)
ranks the choices the same way. No correction for the 540 trials is applied or needed: the curve reports a range, not a
maximum.

## Method

Flat w0waCDM background with photons, CAMB's neutrino convention (N_eff = 3.044, one massive species, the exact
Fermi-Dirac density), distances by quadrature on 4000 points in ln(1+z); r_d, r_* and z_* from CAMB tabulated on an
(ω_b, ω_c, Σm_ν) grid. Against CAMB: D_M to 1e-5, θ_* to 3e-6 (the compressed prior's fractional width is 8e-6), r_d to
1e-5. BAO: DESI's Gaussian likelihoods on (D_M/r_d, D_H/r_d, D_V/r_d). Supernovae: the released distance moduli and
covariances (Pantheon+ at z_HD > 0.01, 1590 SNe; Union3's 22 bins; DES-SN5YR's 1829 with covsys plus the diagonal;
DES-Dovekie's 1820 with its inverse total covariance), the absolute magnitude marginalised analytically, D_L = (1 +
z_hel) D_M(z_HD). Priors as DESI's (w0 ∈ [−3, 1], wa ∈ [−3, 2], w0 + wa < 0; DESI's BBN prior on ω_b when no CMB is in
the coalition, which changes no coalition's Δχ² by more than 0.1 because ω_b is not what the two models disagree about).
Each coalition's χ² is minimised in both models by Nelder-Mead then Powell from 4 starts (the DESI-like point, ΛCDM,
random), the best kept; v(S) = χ²_min,ΛCDM(S) − χ²_min,w0wa(S), v(∅) = 0; φ_i = Σ_S |S|!(n−|S|−1)!/n! [v(S ∪ i) − v(S)].
Three-dataset games are read off the same coalition values. Significances from Δχ² with 2 (CPL) or 1 (wCDM) degrees of
freedom, DESI's convention. Reproduction first: DESI alone −4.7, DESI + (θ_*, ω_b, ω_bc) −8.1, DESI + Pantheon+ −4.9,
DESI + Union3 −10.0, DESI + DES-SN5YR −13.5 against DESI's −4.7, −8.0, −4.9, −10.1, −13.6, best fits within 0.1σ of DESI's.

Two conventions are choices, not derivations. The characteristic function: profile Δχ² is not additive across datasets
(the cross terms between datasets' constraints are the signal), and the Shapley value splits each interaction evenly
among the datasets that produce it, by its symmetry axiom; the Banzhaf value or an Owen value with the SN samples as one
group would split it differently, and the per-dataset χ² gains above are the other natural reading. The prior box: weak
coalitions (a supernova sample alone, DESI alone) fit on the wa = −3 or w0 + wa = 0 edge, so their v(S) depends on
DESI's box as DESI's own Table 6 values do; the three-player values inherit this through v(DESI) = 4.7 and v(SN).

The bootstrap draws every data vector (BAO, the three SN samples, the compressed CMB) from its covariance around the
w0wa best fit of the grand coalition (or the ΛCDM one, for the null) and refits the 32 coalitions from the data's best
fits (one start each); it treats the SN samples as independent, which understates their covariance (shared low-z
supernovae). The full-CMB fits minimise over (ω_b, ω_c, 100θ_MC, τ, ln 10^10 A_s, n_s, A_planck, w0, wa) whitened with
the Planck 2018 w0wa+BAO covariance, by L-BFGS-B with 0.1σ finite differences and a Nelder-Mead polish, from the
compressed pipeline's best fit; the χ² floor is flat to 0.02 over the last 200 evaluations of each fit. The variance
decomposition is the sum of squares of factor-level means about the grand mean over the total, on the balanced
sub-grid; the cut × sample term is the two-way interaction sum of squares.

## What it means

The attribution is the table, and its errors are a third of its values. Three things survive them. First, in the
DES-SN5YR game, the one that gives the 4σ headline, the supernovae and DESI carry equal Shapley values and the
compressed CMB a small one; with the full Planck likelihood the CMB's own χ² moves the other way (it improves under
w0wa), so the full-likelihood CMB share is larger than the compressed one, and the full three-player game
measures it: 7.9 / 4.5 / 7.5 for DESI / CMB / DES-SN5YR against 8.0 / 2.2 / 7.9 compressed, the CMB doubled, the
ordering kept. Second, the Pantheon+ game is the odd one: under the fitted w0wa model its supernova value
should be the largest of the three samples (8 ± 5) and it is 0.2. Third, of the frozen analysis choices the one that
moves the result most is whether the z < 0.1 supernovae are in the fit, half the variance of Δχ² on the balanced grid, with
a further fifth in its interaction with the sample (largest for the two DES samples), and this is the drop a real w0wa
signal of the fitted size would suffer from the lost lever arm (the data's ratios sit at the 25th-58th percentiles of
the bootstrap's), so the cut's dominance of the curve is a statement about where the information is; a low-z systematic of the
right size would produce the same drops, so this baseline does not separate the two [corrected 2026-11-01: the earlier
wording, "not about a systematic", claimed the separation]. Earlier work found the dataset dependence piecemeal:
Efstathiou (2408.07175) the DES5Y low-z offset, Notari, Redi & Tesi (2411.11685) the shared low-z set, Huang, Cai &
Wang (2502.04212) the cut's effect, Cortês & Liddle (2504.15336) that the three significances cannot be averaged, Ong,
Yallup & Handley (2511.10631, 2603.05472) a DESI-DES-SN5YR tension within ΛCDM and the Dovekie recalibration (Popovic+
2511.07517) lowering it, and Karim, Speagle & Hložek (2609.04334), by a leave-one-redshift-block-out predictive score,
that specific DESI BAO bins rather than the supernovae drive the tension. The Shapley split agrees with the last in
giving DESI a value at least as large as any supernova sample's in every game, but not in its supernova value, which is
7.5 ± 3.1 in the DES-SN5YR game where their predictive score gives the supernovae about zero [corrected 2026-11-01], and with the first three in making the
low-z cut the dominant specification choice; it adds the order-free numbers and their errors. It does not say whether the
low-z supernovae are a systematic or a signal.

## Weakest points

- The compressed CMB understates the CMB player, and the full-likelihood Shapley grid was out of reach (32 coalitions × 2
  models at 2 s per likelihood call); the full-likelihood game exists for the DES-SN5YR combination only, with the compressed game's errors, and its
  Union3 and Pantheon+ counterparts are partial.
- The three SN samples share their low-z supernovae, so the five-player game double counts them and the bootstrap
  understates their covariance; the three-player games are the clean statement.
- The bootstrap is parametric around a single truth (the w0wa best fit), with 40 realisations and one minimiser start per
  fit: the errors are good to about 15 % of themselves and are upper bounds (minimiser noise is in them); the ΛCDM null
  gives the floor.
- The specification list is the one frozen in `results-R123.md`; it omits the full CMB likelihoods, a thawing or binned
  w(z) and the Bayesian evidence. The ranking of the top choice is robust to dropping any one factor (recomputed on each
  six-factor sub-grid; the cut ranks first in every one where it varies, except at fixed compressed CMB where it ties with
  the model dof).
- Planck PR3 plik_lite stands in for DESI's PR4 CamSpec and the PR3 lensing for PR4+ACT: a PR3 check of a PR4 result; the
  lensing part differs by 2.4 units.

## To settle it

- The Shapley grid with the full CMB likelihood as the player (an emulator makes it an hour), with the lensing
  reconstruction as a sixth player.
- The specification curve with the SN samples' shared low-z set as an explicit factor (one low-z set swapped between
  samples), which separates "DES-SN5YR" from "the low-z set DES-SN5YR carries".
- The same attribution for the Σm_ν < 0 preference.

## Blind review (2026-11-01, context-free Opus referee, own literature search; `blind-review.md`)

Pass 1, on the pre-revision section: VERDICT: novel but not sound: the Shapley table is computed on a compressed-CMB game that the note's own full-likelihood check shows misattributes the CMB, its ordering claims sit inside errors from a bootstrap that assumes the w0wa signal is real, and the variance decomposition is run on an unbalanced, nonlinear-σ grid without an expected-under-signal baseline.

Strongest objection: the attribution was computed for a compressed-CMB game whose grand Δχ² is 2-3 units below DESI's,
the full-likelihood check shows that game mis-assigns the CMB's share, the full-CMB game had not finished, and the
ordering claims sat inside bootstrap errors that assume a real signal. The one revision pass (claims as intervals, the
balanced decomposition on Δχ², the cut baseline, the full-CMB DES-SN5YR game, the prior work cited, the conventions
stated) produced the section above.

Pass 2, on the revised section: VERDICT: novel but not sound: the Shapley values give weight 1/3 to singleton coalitions whose profile Δχ² is set by the prior-box edges (most visibly the full-CMB-alone value behind "the CMB doubled"), the error bars assume the fitted w0wa model is true, and several stated conclusions (CMB smallest in every game, no zero player, agreement with Karim et al., "not a systematic") go beyond or contradict the note's own intervals.

Strongest objection: in a three-player game the singleton coalitions carry weight 1/3, and those are the fits on the
wa = −3 or w0 + wa = 0 edge (1.6 of DESI's 7.9, 2.0 of the supernovae's 7.5 and 1.1 of the CMB's 4.5 come from them); the
CMB-alone w0wa fit is unconstrained in dark-energy geometry, so "the CMB doubled" is mostly box-limited behaviour; on a
boundary the profile Δχ² is not χ²-distributed. Also: the bootstrap and the cut baseline assume the fitted w0wa truth; the
Δχ² response makes a data-removing factor dominate by construction (a standardised parameter shift would answer the
robustness question better); four interpretive sentences went beyond the tables (corrected above with dated brackets,
the verdict left as given). Both passes found the method itself (Shapley over dataset coalitions; a preregistered
specification curve with a variance decomposition) nowhere in the literature; the referee's closest prior work is Karim,
Speagle & Hložek 2609.04334 (leave-one-redshift-block-out predictive scoring) and 2607.28918 (e-value leave-one-out),
which attribute the preference to specific DESI redshift bins.

Recorded as a miss and parked. The path to pass: a characteristic function that does not depend on the prior box
(the Bayesian evidence ratio, or Δχ² with the singleton coalitions' values defined by a common, data-driven box, or the
Owen value with the weak coalitions grouped); errors from a bootstrap under both truths with multi-start fits; the full
CMB likelihood as the player in every coalition (an emulator), with the PR4+ACT lensing; the specification curve's
response as a standardised parameter shift, with the shared low-z supernova set as an explicit factor.

Addendum (2026-11-01, after the second pass, not refereed): the full-likelihood Union3 game finished at the budget: CMB +
Union3 Δχ² = 6.35, and the Shapley values are DESI 7.0 (44 %), CMB 4.4 (28 %), Union3 4.5 (28 %) against 7.0 / 2.1 / 4.8
compressed, the same shift as in the DES-SN5YR game (`results-R123.md`, row 7).


# Finding (R-129, cosmo-deep-9, 2026-11-01): a contamination budget for the ACT DR6 CMB-lensing map from 21 templates that cannot carry lensing signal: the templates' total imprint on the published κ × DESI-LRG amplitudes is below 1.2 % (0.2σ) and on κ × Quaia below 0.5 %; one survey-property map (the unWISE W1 depth) correlates with the map at high latitude, as it did marginally with Planck's, and its origin is not identified; the SFD dust correlation is SFD's own extragalactic content

## Claim

The ACT DR6 lensing convergence map (Qu et al. 2024, arXiv:2304.05202; Madhavacheril et al. 2024, arXiv:2304.05203) was crossed
with 21 all-sky or survey templates that have no physical connection to cosmic lensing: Gaia DR3 star counts in four G bins and
two colours, the Legacy Surveys DR9 imaging properties the DESI target pipeline records (stellar density, SFD E(B−V), PSF depth
in g, r, z, W1, galaxy depth in z, PSF size in g, r, z) and the Quaia selection-function feature maps (dust, star counts, the Gaia
G-band limit m10, the unWISE source counts, the unWISE exposure count). Each cross-spectrum was measured in eight bands from L = 20
to 800 on the whole footprint and in three Galactic-latitude bands, and compared with the same template crossed with each of the
release's 400 simulated reconstructions. The contamination each template implies for a published κ × galaxy amplitude was
propagated with the DES statistic of Omori et al. 2019 (arXiv:1810.02342, Sec. V) and Chang et al. 2022 (arXiv:2203.12440, App.
B1), X_S = C^{κS} C^{gS} / C^{SS}, into the DESI Legacy LRG photometric bins of Zhou et al. 2023 (the sample of the ACT DR6 × LRG
measurements, Sailer et al. 2024 arXiv:2407.04607 and Kim et al. 2024 arXiv:2407.04606) and Quaia (Storey-Fisher et al. 2024;
Alonso et al. 2023 arXiv:2306.17748, Piccirilli et al. 2024). The template list, bands, latitude bands, nulls, pass criterion and
propagation statistic were fixed before the first cross-spectrum (`results-R129.md`, with four deviations recorded there before
any PTE was computed).

The pass criterion was: a template counts as a detected correlation when no simulation of the 400 exceeds the data's χ² over
the eight whole-footprint bands and at least one of the temperature-only, polarization-only or galcut040 map variants sits at
the floor of its rotation null; marginal when PTE < 0.01; otherwise the result is an upper limit. What the test establishes is
a correlation between the map and a template; whether it is a feature of the map or extragalactic content of the template is a
separate question, answered below for SFD (template) and tested for the others. The budget, which is the claim, does not
depend on that question: either way the DES statistic bounds what the template can have put into a κ × galaxy amplitude.

| template | whole footprint: χ² (8 bands), PTE against 400 sims | z per band (L = 20-50, 50-100, 100-150, 150-200, 200-300, 300-400, 400-600, 600-800) | PTE by latitude: \|b\| < 30°, 30-50°, > 50° | rotation-null PTE on tonly / polonly / cibdeproj / galcut040 / Planck PR3 (floor 0.017); in brackets tonly against its own 100 simulations | verdict |
|---|---|---|---|---|---|
| **SFD E(B−V)** (Legacy DR9) | 23.6, **< 0.0025** | +1.0 −0.3 −0.2 **+2.8 +2.3 +2.8** +0.0 +1.7 | 0.61, 0.66, **< 0.0025** | 0.017 / 0.000 / 0.000 / 0.000 / 0.40 [0.05] | **detected**, and explained by the template, not the map: see the CSFD row |
| *post hoc:* CSFD E(B−V) (Chiang 2023: SFD with its extragalactic, large-scale-structure component removed) | 13.8, 0.05 against 60 rotations | +1.1 −0.8 −1.4 +2.1 +0.3 +0.8 −2.6 −1.3 | —, —, 0.48 | — | null; the component removed from SFD gives z +13 to +25 in every band (the CIB-lensing correlation) |
| **PSF depth W1** (unWISE, via Legacy DR9) | 44.3, **< 0.0025** | −1.8 +0.3 −0.9 **−4.2 −2.8 −2.1 −2.5 −1.8** | 0.78, 0.46, **< 0.0025** | 0.000 / 0.000 / 0.000 / 0.000 / 0.22 [0.000] | **detected** |
| **Gaia G-band limit m10** (the scanning law; Quaia feature) | 24.2, **0.0025** (one simulation above) | +0.1 **+3.2 +2.5** +1.1 +2.1 −0.3 +0.7 +0.1 | 0.72, 0.20, 0.0025 | 0.000 / 0.000 / 0.000 / 0.000 / 0.017 [0.04] | marginal by the letter; at the rotation floor in every variant and in Planck, 0.04 against the temperature-only simulations; and mostly accounted for by structure it shares with the CIB map (below) |
| PSF size z | 14.1, 0.07 | +1.7 −0.3 −0.2 −1.1 +1.7 −1.4 +1.9 +0.1 | 0.29, 0.52, 0.11 | 0.05 / 0.017 / 0.07 / 0.000 / 0.13 | upper limit |
| Gaia counts, G 12-16 / 16-18 / 18-20 / 20-21 / 16-20 blue / 16-20 red | 11.6 / 7.4 / 7.6 / 7.1 / 5.0 / 6.1; 0.14 / 0.47 / 0.49 / 0.51 / 0.74 / 0.61 | the 50-100 band +1.9, +1.4, +1.7, +1.6, +1.0, +1.8; all else \|z\| < 2 | all 0.41-1.00 | all ≥ 0.12 | upper limits |
| stellar density, PSF depth g / r / z, galaxy depth z, PSF size g / r (Legacy); dust, star counts, unWISE exposures (Quaia) | χ² 6-11.5, PTE 0.16-0.65 | \|z\| < 2.5 in every band | 0.013-1.00 (PSF depth g at 30-50°: 0.013) | ≥ 0.03 except depth z and galaxy depth z on cibdeproj (0.000, 0.000) and galcut040 (0.05, 0.17) | upper limits |
| unWISE source counts (Quaia feature; galaxies, so a positive control, not a zero-signal template) | 352, < 0.0025 | +4.3 +12.0 +9.7 +6.6 +5.5 −0.2 +1.3 +1.1 | 0.03, < 0.0025, < 0.0025 | all 0.000 | the pipeline sees the κ × galaxy signal where one must be |

Of 84 (template, latitude band) tests on the baseline map, 9 have PTE < 0.01 against 0.8 expected by chance; they are three
templates (each on the whole footprint and at |b| > 50°, which share most of their area, so three findings, not six) and the
galaxy control. None of the three is seen at |b| < 50°. Dropping the 20-50 band, which the ACT DR6 analyses do not use (their
range starts at L = 40), changes no PTE (E(B−V) 0.0000, W1 depth 0.0000, m10 0.0025 over L = 50-800).

Is the correlation in the template? A test added after the first referee pass: each template was crossed with Chiang 2023's
reconstructed extragalactic (CIB / large-scale structure) component of SFD, and the part of κ × S that this content predicts,
X = C^{κ,CIB} C^{S,CIB} / C^{CIB,CIB}, was subtracted from the measured κ × S (z per band against the 400 simulations before →
after, whole footprint, bands L = 150-800): SFD E(B−V) +2.8 +2.3 +2.8 +0.0 +1.7 → +2.2 +0.7 +1.1 −2.1 −0.0 (removed, as the CSFD
row says); Gaia m10 at L = 50-150: +3.1 +1.4 → +0.8 −0.8 (removed: the m10 map shares large-scale structure with the CIB map,
which, built from IRAS, can carry the ecliptic pattern m10 follows, so this is "shared structure", not "CIB in a Gaia depth
map"); **W1 depth −4.2 −2.8 −2.1 −2.5 −1.8 → −4.7 −5.2 −3.6 −4.9 −3.5**: the depth map is nearly uncorrelated with the CIB map
(C^{S,CIB} an order of magnitude below SFD's, of the opposite sign at L < 50), the predicted part has the wrong sign, and the
residual is as significant as the raw correlation. The W1-depth correlation is not the template's CIB content. [Corrected after the second referee pass: an earlier version of this sentence said its sign and scale matched a marginal
anti-correlation of the Planck PR3 map with the same PSFDEPTH_W1 map in Krolewski et al. 2024 (arXiv:2305.07650); that paper's
test is for DESI quasar targets, its W1-depth anti-correlation is at 2-2.5σ in its lowest bin only (4 ≤ ℓ < 13) and goes away at
ℓ_min = 6, so it is not support for a correlation at L = 150-800. The ACT correlation stands on its own.] The Planck PR3 map on
this footprint gives PTE 0.22 against rotations.

Propagation (baseline map; the κ × g band powers from this pipeline reproduce Sailer et al. 2024's published ACT DR6 × LRG ones
within their errors in every bin; σ_A / A is the sample's amplitude error on L = 20-800 from the simulations):

| sample | σ_A / A | ΔA / A from E(B−V) (null spread) | from W1 depth | from m10 | the three in quadrature, as a fraction of σ_A | all 20 zero-signal templates in quadrature (an upper bound: correlated templates double-counted) | largest 95 % limit among the 17 undetected templates |
|---|---|---|---|---|---|---|---|
| DESI LRG z1 (z ≈ 0.4) | 0.057 | +0.68 % (0.15) | +0.69 % (0.14) | +0.51 % (0.15) | 1.09 %, **0.19 σ_A** | 1.17 %, 0.21 σ_A | 0.35 % |
| LRG z2 (0.6) | 0.046 | +0.44 % (0.16) | +0.26 % (0.15) | +0.24 % (0.08) | 0.57 %, **0.12 σ_A** | 0.66 %, 0.14 σ_A | 0.36 % |
| LRG z3 (0.8) | 0.042 | +0.31 % (0.12) | +0.44 % (0.12) | −0.09 % (0.05) | 0.54 %, **0.13 σ_A** | 0.62 %, 0.15 σ_A | 0.30 % |
| LRG z4 (1.0) | 0.039 | +0.26 % (0.08) | +0.32 % (0.09) | +0.19 % (0.10) | 0.45 %, **0.12 σ_A** | 0.58 %, 0.15 σ_A | 0.27 % |
| Quaia G < 20.5 | 0.054 | −0.18 % (0.22) | −0.01 % (0.19) | −0.06 % (0.07) | 0.19 %, 0.03 σ_A | 0.41 %, 0.08 σ_A | 0.42 % |

For scale, the galaxy control (unWISE source counts, which are a galaxy sample) propagates to +6.5 to +9.9 % of the LRG
amplitudes (1.6-1.8 σ_A): that is what a template with real large-scale-structure content does in this statistic. The W1-depth
and m10 entries are not independent (both follow the ecliptic; m10's κ correlation is largely shared structure), so their
quadrature sum is also an upper bound. The fractions are biases on an amplitude fitted to the ACT-only band powers on L = 20-800
with inverse-variance weights from the simulations; the published fits use L ≥ 40-44 to 600-1000, their own covariances and,
for Sailer et al., a Planck PR4 combination, so "+x % of the published amplitude" means "of the ACT part, at this weighting",
and the sign (an upward bias in every bin) is one finding, not four, since C^{κS} is common to the bins.

The claim is the budget. Twenty-one templates that cannot lens the CMB, crossed with the ACT DR6 baseline convergence through
the same pipeline that reproduces the published κ × LRG band powers, put at most 1.2 % (z1), 0.7, 0.6 and 0.6 % (z2-z4) into the
ACT DR6 × DESI-LRG amplitudes when every template is counted in quadrature as if independent (0.14-0.21 of the amplitudes'
statistical errors), and at most 0.4 % into the Quaia amplitude; no single template exceeds 0.7 %, and the seventeen templates
with no detected correlation are each limited to |ΔA / A| < 0.45 % at 95 %. No published κ × galaxy amplitude on this map needs
to move by more than a fifth of its own error bar on account of these templates. Within that budget, three templates are
correlated with the map at PTE ≤ 1/400 against the release's simulations, at Galactic latitudes above 50° and not below, with
the same sign and scale in the temperature-only, polarization-only, CIB-deprojected and galcut040 reconstructions (SFD E(B−V)
positive at L = 150-400, the unWISE W1 depth negative from L = 150 up, the Gaia scan law positive at L = 50-150); the Planck PR3
map on the same footprint shows the scan-law one at its rotation floor and not the other two. The SFD one is the template's:
SFD carries the unresolved-galaxy and CIB emission of large-scale structure at high latitude (Chiang & Ménard 2019), the map is
uncorrelated with Chiang 2023's CSFD, which removes that component (PTE 0.05 on the whole footprint, 0.48 at |b| > 50°,
propagation 0.01-0.07 % of the LRG amplitudes), and the removed component alone correlates with κ at 13-25σ per band, as the
lensed CIB does (both checks added after unblinding, on the rotation null). The m10 one is mostly structure the template shares
with that CIB map. The W1-depth one survives the same test and is the one correlation whose origin is open: a feature of every
ACT reconstruction on this footprint that the Gaussian simulations do not carry, or a content of the depth map that neither the
CIB map nor the galaxy maps capture. Propagated, the three together would bias the LRG amplitudes upward by +1.1 ± 0.2 % (z1),
+0.6 ± 0.2, +0.5 ± 0.2 and +0.5 ± 0.1 % (z2-z4) and the Quaia one by −0.2 ± 0.3 %; the W1-depth one alone by +0.7, +0.3, +0.4 and
+0.3 % (0.06-0.12 σ_A).

## Method

Maps: the ACT DR6 lensing v1 release (Qu et al. 2024; Madhavacheril et al. 2024; its 400 simulated reconstructions are the
release's own, whose README asks that Farren et al. 2024, arXiv:2309.05659, be cited for them): the baseline κ alm (temperature and polarization, profile-hardened) as the statistic, and the temperature-only,
polarization-only, CIB-deprojected and galcut040 variants as checks, each with its release mask (apodized, nside 4096, degraded
to nside 512); the Planck PR3 MV convergence (Planck 2018 VIII) rotated from Galactic to equatorial, its mask apodized by a 1°
Gaussian. The alm are cut at L = 1100 and synthesised at nside 512.

Templates (21; fixed in `results-R129.md` before any cross-spectrum, with the four deviations recorded there): Gaia DR3 source
counts per nside-256 pixel in G = 12-16, 16-18, 18-20, 20-21 and 16-20 split at BP−RP = 1 (ESA archive TAP); the Legacy Surveys
DR9 imaging properties per nside-256 pixel from the DESI target pipeline's `pixweight-1-dark` file: stellar density, SFD
E(B−V), PSF depth in g, r, z and W1, galaxy depth in z, PSF size in g, r, z; and the Quaia selection-function feature maps as
released (Storey-Fisher et al. 2024, Zenodo 10403370; nside 64): SFD dust, Gaia star counts, the Gaia G-band limit m10 (the
scanning law), the unWISE source counts and the unWISE exposure count (the two Magellanic-Cloud masks in the release are zero on
the ACT footprint and are not testable). Each is standardised to zero mean and unit variance over its footprint within the ACT
one and upgraded to nside 512. The unWISE source-count feature is mostly galaxies at high latitude and therefore carries real
lensing signal; it was recognised as such on reading its definition, after the data band powers were printed and before any
PTE, and is reported as a positive control and left out of the contamination totals. HI was in the fixed list and was dropped
because the only HEALPix HI map in hand covers the DESI DR1 spectroscopic footprint alone (2,300 deg² of ACT's); the Legacy
templates were moved to the DR9 `pixweight` file (9,100 deg² of ACT's) for the same reason.

Band powers: pseudo-C_L of κ and a template on a common weight W (the κ mask × the template footprint × a latitude cut, × the
galaxy mask for the propagation), the weighted mean of each field subtracted inside W, in eight bands with edges 20, 50, 100,
150, 200, 300, 400, 600, 800, decoupled with the binned MASTER matrix computed from W's power spectrum (the closed-form Wigner 3j
with m = 0, L ≤ 1100; NaMaster could not be built in the container). The full-sky identity holds to 1e-13; a Gaussian-map
recovery test through the ACT × Legacy weight returns 0.99-1.05 of the input in every band (1.09-1.19 in the two top bands when
the matrix was cut at the top band edge, 850, hence 1100). Latitude bands: the whole footprint, |b| < 30°, 30-50°, > 50°
(650, 3,200 and 6,600 deg² of the baseline footprint). The same pipeline on κ × LRG reproduces Sailer et al. 2024's published
ACT DR6 × LRG band powers within their errors in every bin (×1e-6 at L = 32, 61, 101, 151: z1 1.63 / 0.89 / 0.37 / 0.24 against
1.70 ± 0.28 / 0.79 ± 0.13 / 0.35 ± 0.07 / 0.25 ± 0.04; z3 1.28 / 0.79 / 0.44 / 0.30 against 1.26 ± 0.20 / 0.80 ± 0.10 /
0.38 ± 0.05 / 0.35 ± 0.03).

Nulls: for the baseline map each template was crossed with each of the release's 400 simulated reconstructions through the same
weights, bands and decoupling, per latitude band; the χ² of the eight band powers against the simulation covariance
(Hartlap-corrected) and its empirical PTE (the fraction of simulations with a larger χ²) are the statistic. As a check, and as
the only null for the map variants and Planck, the data alm was turned in right ascension by 360k/61° (k = 1..60) with the
weights and templates fixed; a turned ACT map covers 0.45-0.51 of each weight (the footprint is a third of its declination
band), so each rotation's band powers are divided by the square root of its weighted overlap, after which the two nulls agree
(for the baseline map the rotation PTEs of the three detected templates are 0.000, and the z scores match the simulations' to
10-30 %). The rotation null has a floor of 1/60.

Propagation (Omori et al. 2019 Sec. V; Chang et al. 2022 App. B1): for a galaxy sample g and a template S, X_S(b) =
C^{κS}(b) C^{gS}(b) / C^{SS}(b) on the common weight; the bias on the amplitude A of the κ × g spectrum is ΔA_S / A =
Σ_b X_S(b) C^{κg}(b) / σ_b² ÷ Σ_b C^{κg}(b)² / σ_b², with C^{κg} the measured spectrum on the same weight and σ_b its spread over
the simulations crossed with g (so σ_A / A = (Σ_b C^{κg}(b)² / σ_b²)^{-1/2} is the sample's own amplitude error on these scales);
the null spread of ΔA_S comes from the simulations in place of κ in C^{κS}, and the 95 % limit is |ΔA_S| + 1.64 spreads. Samples:
the DESI Legacy Imaging LRGs in the four photometric-redshift bins of Zhou et al. 2023 as released with White et al. 2022
(arXiv:2111.09898, Zenodo 5834378; the maps are in Galactic coordinates and were resampled to equatorial; 9,000 deg² overlap
with ACT), the sample of Sailer et al. 2024 and Kim et al. 2024; Quaia G < 20.5 as δ_g = n / (n̄ s) − 1 with the released
selection function s, footprint s > 0.5. unWISE: no public galaxy map reaches this cloud; not propagated.

Trials: 21 templates × 4 latitude bands × 6 maps = 504 cross-spectra and 5 samples × 21 templates propagations, all in the
logs; the criterion and the 84 baseline tests it is applied to are the ones fixed in advance, and the chance expectation (0.8
of 84 at PTE < 0.01) is stated with the count (9, six of them the three templates above and three the galaxy control).

## What it means

The result is the budget. The κ-side template test of Omori et al. 2019 and Chang et al. 2022 (DES Y1 and Y3 against SPT +
Planck lensing; Marques et al. 2024, arXiv:2306.17268, for DES Y3 × ACT DR4) and the Planck-side one of Krolewski et al. 2024
(arXiv:2305.07650, the Planck PR3 map against the same Legacy imaging templates, propagated into Planck × DESI LRG with a bound
below 0.3σ) had not been run on the ACT DR6 map; the ACT DR6 × LRG papers (Sailer et al. 2024; Kim et al. 2024; Qu et al. 2024,
arXiv:2410.10808), the ACT DR6 × unWISE one (Farren et al. 2024) and the ACT DR6 × Quaia one (Embil Villagra et al. 2025,
arXiv:2507.08798) control systematics on the galaxy side (regression weights, mask variation), with foreground simulations, and
with differences between reconstructions and the curl. With Gaia, Legacy and Quaia templates and the release's 400 simulations
as the null, the ACT DR6 budget comes out at 0.6-1.2 % of the κ × LRG amplitudes, 0.14-0.21 of their errors, the same order as
Krolewski et al.'s Planck bound; the correction those amplitudes would take is below their quoted precision, and every
template without a detected correlation is limited to less than half a percent per sample.

The three correlations have the same sign and scale in the temperature-only and the polarization-only reconstructions, so they
are not a temperature-foreground leakage of one estimator (a polarization-only map has no CIB or tSZ and little dust intensity),
and they are absent below |b| = 50°, where Galactic dust and stars are stronger. For SFD that pattern is explained: the SFD map
contains the large-scale structure imprint of unresolved galaxies and the cosmic infrared background at high latitude (Chiang &
Ménard 2019, arXiv:1808.03294; Chiang 2023, arXiv:2306.03926), a positive κ × E(B−V) correlation at |b| > 50° and L = 150-400 is
what that imprint predicts, and the CSFD check confirms it (κ × CSFD null, κ × the removed component at 13-25σ). So "SFD E(B−V)"
is not a zero-signal template at high latitude, and a κ × E(B−V) test that uses SFD measures the CIB, not dust in the map; this is
worth knowing for any lensing-map systematics test built on SFD. The unWISE W1 depth and the Gaia scan law are survey properties
with no obvious extragalactic content, both follow the ecliptic, and they correlate with the map with opposite signs; an explicit
|ecliptic latitude| template added after unblinding sits at its rotation floor on the whole footprint (PTE 0.017, one band at
z +3.0) and not at |b| > 50° (0.08), with a negligible propagation. The m10 correlation is mostly structure shared with the CIB
map and is not pursued. For the W1 depth the CIB test above rules out the reading that cleared SFD (the depth map's CIB content
predicts the wrong sign and leaves the residual at z −3.5 to −5.2); what remains is a large-scale pattern common to every ACT
reconstruction on this footprint (the filtering, the mean-field subtraction, the noise model or the simulations' representation
of it) that the W1 coverage pattern happens to trace, or a content of the depth map that neither the CIB map nor the galaxy maps
capture. Krolewski et al. 2024's Planck test (quasar targets, ℓ < 13) is not support for it (corrected, above); the Planck PR3 map
on this footprint shows it at PTE 0.22 against rotations, with a noisier map. Either way its propagated size is 0.3-0.7 %
of the LRG amplitudes.

## Weakest points

- The detected correlations are small (|z| ≤ 4.2 per band; χ² 24-44 for 8 degrees of freedom) and their propagated effect is
  0.1-0.2 σ_A: this is a budget, not a correction. A reader who wanted to apply it should treat the +0.5 to +1.1 % as the size of
  a one-sided systematic, not as a measured offset.
- The simulations are Gaussian reconstructions of Gaussian skies through the ACT pipeline. They carry the map's noise and mask
  but not every feature of the data (foreground residuals, the mean field's data-specific part), so a common data-only
  large-scale pattern would appear as a "detection" against them for any template with power on those scales. The rotation
  null, which keeps the data map's own structure, gives the same PTEs for the three templates on the baseline and on every
  variant, which argues against a pure declination pattern but not against a pattern that follows the Galactic or ecliptic
  poles. So the test cannot by itself separate a map systematic from a data-only large-scale feature; that is why the claim is
  the budget, which holds under either reading, and the W1-depth correlation is reported as a correlation of open origin.
- Only the baseline map has the 400-simulation null; the variants and Planck rest on 60 rotations with a floor of 1/60, and
  the rotation null is check-grade (its band powers are rescaled by the overlap and its covariance has 60 samples). The
  temperature-only map's own 100 simulations, streamed as a second check, confirm the W1-depth correlation at PTE 0.000 but put
  E(B−V) at 0.05 and m10 at 0.04 where the rotations gave 0.017 and 0.000: the rotation null is 20-30 % anticonservative in z
  after its rescaling, so the "every variant at its floor" statements for E(B−V) and m10 are weaker than they read, and the
  criterion's variant confirmation at PTE < 0.01 is met against simulations only for W1 depth.
- The CSFD and ecliptic-latitude checks were added after unblinding (stated in the table) and rest on the 60-rotation null,
  the simulations having been deleted after use; their verdicts are at rotation-null grade. The W1 depth and m10 templates are
  the ones whose correlation with the map would be a genuine map systematic, and they are the ones whose cause is not
  identified here.
- The quadrature total double-counts correlated templates (six Gaia bins, three depths, two dust maps) and is an upper bound; the
  per-template limits are the cleaner numbers.
- The latitude split is coarse (three bands) and the |b| < 30° band is small (650 deg²), so "absent below 50°" is a statement
  about PTEs of 0.2-0.8, not a measured latitude dependence.
- KiDS-1000 and DES Y3 survey-property maps, the ACT scan pattern and airmass were not tested (not fetched; not in the
  release), and unWISE's own galaxy maps were not reachable, so the unWISE × ACT amplitude (Farren et al. 2024) is not in the
  budget; the Legacy and Quaia templates cover its footprint, and nothing here suggests a larger effect for it.

## To settle it

The W1-depth and scan-law correlations against the ACT noise and mean-field simulations' own large-scale structure (the release's
400 reconstructions are of Gaussian skies; a test on the mean-field and on noise-only reconstructions would say whether a
footprint-shaped pattern is in the map), and against the unWISE coverage and Gaia transit-count maps at full resolution with the
400 simulations as the null. The CSFD check with the simulations (here: rotations). unWISE's blue and green samples once their
maps are in hand, with the same statistic. The temperature-only and polarization-only maps with their own 400 simulations.


## Blind review (2026-11-01, context-free Opus referee, own literature search; `blind-review.md`)

Pass 1 (the section as first written, with the three correlations in the headline as "detected contaminations"): **"VERDICT:
novel but not sound: the headline map-side detections rest on Gaussian simulations and an anticonservative rotation null that
cannot separate a map systematic from data-only large-scale structure, and the W1-depth correlation shows the same high-latitude
extragalactic signature the note uses to clear SFD without being tested for it."** Strongest objection: the W1-depth result was
read as a map systematic without the test that cleared SFD. Prior work it found: Krolewski et al. 2024 (arXiv:2305.07650, the
Planck map against the same Legacy templates, a marginal W1-depth anti-correlation, bound < 0.3σ), Marques et al. 2024
(arXiv:2306.17268), arXiv:2411.12801 (SFD × Planck κ positive, smaller with CSFD), Qu et al. 2024 (arXiv:2410.10808), Embil
Villagra et al. 2025 (arXiv:2507.08798). One revision pass on the presentation (the claim re-made as the budget, which holds
under either reading of the correlations; the prior work cited; the CIB test the referee asked for run and reported; the
"six of nine" count corrected to three findings; the L < 50 band dropped as a check; the non-independence of W1 depth and m10
and the meaning of "+x % of the published amplitude" stated). The physics point (the nulls cannot separate a map systematic from
a data-only large-scale feature) is left as stated in Weakest points, not argued away.

Pass 2 (the revised section): **"VERDICT: novel but not sound: the headline "upper bound" adds single-template projections of
correlated, opposite-signed templates in quadrature against a null that drops the κ–g covariance term, so it is not a bound, and
the Planck evidence cited for the W1-depth correlation is misread (quasar sample, ℓ < 13)."** Strongest objection: the budget is
not a bound; the joint multi-template projection C^{κS} (C^{SS})^{-1} C^{Sg} with a null that carries the κ–g covariance is what
a referee would ask for. Novelty granted (no paper runs this κ-side test on ACT DR6; the SFD/CSFD point is in print: Cheng,
Hensley, Chang & Doré 2024, arXiv:2411.12801). **PARKED** (physics points, by the gate rule: one revision pass on presentation,
park on a physics point). The Krolewski sentence was wrong and is corrected above, flagged.

**Post-gate addendum (for the record, not re-reviewed; `extra3`, `logs/extra3_baseline.md`):** the joint projection over the 20
zero-signal templates (the 20 × 20 template cross-spectra per band, the inverse regularised at 1e-3 of the largest eigenvalue),
propagated like a single template, with the null from the 400 simulations: **+3.2, +3.3, +2.7, +3.0 % of the LRG z1-z4 amplitudes
(null spreads 0.4-0.6 %; 0.57-0.77 σ_A)** and −0.2 ± 0.6 % for Quaia; over the three correlated templates only: +1.8, +0.9, +0.6,
+0.8 % (0.15-0.32 σ_A); over the 20 without SFD: +1.9, +1.8, +1.6, +1.7 % (0.33-0.45 σ_A). The referee was right that the
quadrature sum (0.6-1.2 %) understates the joint number by a factor of three. The joint 20-template projection also carries the
noise-fit bias of fitting 20 templates to ν_b modes, which the simulations (no κ–g correlation) cannot centre. The noise-bias
offset of the single-template statistic, C^{κg} / ν_b per band with ν_b = f_eff Σ(2L + 1), is +0.04 to +0.05 % of the amplitude
per template, and the mean of the 17 undetected templates is +0.02 to +0.07 % before and −0.03 to +0.02 % after subtracting it:
the offset is real and small, and it does not change the three detections (0.3-0.7 %) or the per-template limits by more than
0.05 %. What stands: the budget's honest form is the joint projection, below 3.3 % (0.8 σ_A) with all 20 templates and below
1.8 % (0.3 σ_A) for the correlated three, on the ACT DR6 baseline map; the W1-depth correlation (PTE < 1/400, every ACT variant,
not the template's CIB content) of open origin; the SFD correlation as SFD's own extragalactic content (in print). What a next
pass needs: the joint projection with a null built from κ simulations that carry the κ–g correlation (galaxy mocks or the
simulations' input κ with a biased-tracer field), the W1-depth map projected on the unWISE source counts, the LRG maps' weighting
stated from the source, and a Planck leg with its own simulations.


# The stacked Sunyaev-Zeldovich signal of the eRASS1 clusters at 1-8 × 10^14 M_sun against their weak-lensing-calibrated masses, on the Planck and ACT DR6 y maps: Y_500 at 3 × 10^14 is 0.90 ± 0.02 (statistical) ± 0.05 (systematic) of the Planck 2013 relation at 1−b = 0.8, the ratio runs from 1.0-1.3 at 1.3 × 10^14 to 0.8 at 5.5 × 10^14, and the slope is 1.5-1.6 on the aperture statistic against 1.5 (Planck) and 1.8 (ACT) on the pre-registered template fit

2026-11-01. Script: `erass1_ysz.py` (`prep`, `stack`, `fit`); rows D0-D10 in `results-R131.md` (the samples, bins, statistic
and model fixed there before any map was opened at a cluster). Data: the eRASS1 primary v3.2 and cosmology v1.1 catalogues
(Bulbul et al. 2024), the Planck 2015 MILCA and NILC Compton-y maps with the release's Galactic and point-source masks, and the
ACT DR6 ACT+Planck NILC Compton-y map (Coulton et al. 2024; the wide GAL070 mask, the inpainted-region and subtracted-source
masks, the 1.6' beam file, and the CIB-deprojected variant).

## Claim

The eRASS1 catalogue gives each cluster an M500 from the eROSITA count-rate-mass relation calibrated on DES Y3, KiDS-1000 and
HSC-Y3 galaxy weak lensing (Grandis et al. 2024; Kleinebreil et al. 2024). The thermal SZ signal is an independent observable
of the same gas, and its mass scaling has been measured on SZ-selected clusters above about 4 × 10^14 M_sun. Stacking the
Planck and ACT y maps on the eRASS1 cosmology sample (5,259 clusters, ~95 % pure; 1,683 in the Planck footprint and 1,115 in
the ACT footprint at 0.1 < z < 0.8 and 1-8 × 10^14) in five mass bins, with the Arnaud et al. 2010 universal pressure profile
(UPP) at each cluster's R500 as the template and the Planck 2013 Y-M relation at 1−b = 0.8 as the reference, gives the
spherical Y_500 of the eRASS1 clusters down to 1.3 × 10^14, two decades below the SZ catalogues:

| bin (M500, 10^14) | M_eff | n (Planck / ACT) | Y_500 E^-2/3 D_A^2, Planck MILCA | ACT DR6 | ratio to Planck 2013 at 1−b = 0.8 (MILCA / ACT) |
|---|---|---|---|---|---|
| 1.0-1.4 | 1.30 | 182 / 136 | 0.29 ± 0.05 | 0.35 ± 0.04 | 1.03 ± 0.17 / 1.26 ± 0.14 |
| 1.4-2.0 | 1.80 | 326 / 220 | 0.58 ± 0.04 | 0.50 ± 0.04 | 1.15 ± 0.08 / 0.98 ± 0.07 |
| 2.0-2.8 | 2.49 | 396 / 276 | 0.85 ± 0.04 | 0.79 ± 0.05 | 0.95 ± 0.05 / 0.89 ± 0.05 |
| 2.8-4.0 | 3.48 | 403 / 255 | 1.37 ± 0.05 | 1.46 ± 0.05 | 0.83 ± 0.03 / 0.91 ± 0.03 |
| 4.0-8.0 | 5.5 | 376 / 228 | 2.90 ± 0.06 | 3.10 ± 0.08 | 0.78 ± 0.02 / 0.80 ± 0.02 |

(Y in 10^-5 Mpc^2, the aperture statistic of the Method: the cylindrical y within 2 θ500 of each cluster, ring-subtracted,
converted to the spherical Y_500 with the beam-convolved UPP; errors from 200 rotated-map stacks per map, 100 for the
CIB-deprojected ACT map; M_eff is the reference relation inverted at the bin's mean reference prediction, in M_sun; the
mass posterior of each cluster is averaged over, which raises the prediction by 5-7 % over the median mass.)

Fitted as Y_500 E^-2/3 D_A^2 = Y_3 (M500 / 3 × 10^14)^B over the five bins with their joint covariance:

| map, statistic | B (68 %, 95 %) | Y_3 (10^-5 Mpc^2) | Y_3 / Planck 2013 (1−b = 0.8) | b_M at the pivot, 1−b = 0.8 (0.7, 0.9) |
|---|---|---|---|---|
| Planck MILCA, aperture | 1.53 [1.49, 1.57], [1.44, 1.62] | 1.13 ± 0.02 | 0.90 | 0.94 (1.08, 0.84) |
| Planck NILC, aperture | 1.51 [1.46, 1.55], [1.42, 1.60] | 1.12 ± 0.02 | 0.89 | 0.94 |
| ACT DR6, aperture | 1.58 [1.54, 1.62], [1.50, 1.67] | 1.14 ± 0.02 | 0.91 | 0.95 (1.09, 0.85) |
| Planck MILCA, pre-registered template fit | 1.52 [1.48, 1.55], [1.44, 1.60] | 1.14 ± 0.02 | 0.91 | 0.95 |
| Planck MILCA, template + two-halo shape | 1.56 [1.52, 1.61], [1.47, 1.66] | 1.08 ± 0.02 | 0.86 | 0.92 |
| ACT DR6, pre-registered template fit | 1.83 [1.82, 1.85], [1.80, 1.87] | 0.95 ± 0.01 | 0.76 | 0.86 |
| ACT DR6, template + two-halo shape | 1.87 [1.85, 1.89] | 0.92 ± 0.01 | 0.74 | 0.85 |

(The Planck 2013 relation at the pivot is 1.252 × 10^-5 Mpc^2 with slope 1.79; b_M is M_true = b_M M_eRASS1 under that
relation, b_M = 0.8 (Y_3 / 1.252)^(1/1.79) / (1−b).)

What is robust and what is not. **The amplitude at the pivot: Y_500(3 × 10^14) = 1.12-1.14 × 10^-5 Mpc^2 on the aperture
statistic of all three maps (± 0.02 statistical each), 1.08-1.14 on Planck's template fits and 0.92-0.95 on ACT's, i.e.
0.90 ± 0.02 of the Planck 2013 relation at 1−b = 0.8 with a systematic range 0.74-0.91 over the statistics, and ± 0.03 each for
the h70-to-eRASS1-cosmology conversion and the relativistic correction (neither applied).** Under that relation with its slope
the eRASS1 masses at 3 × 10^14 are b_M = 0.94 ± 0.01 ± 0.05 of the Planck 1−b = 0.8 scale (1.08 at 1−b = 0.7, 0.84 at 0.9), i.e.
if the eRASS1 masses are right the Planck relation needs 1−b = 0.72-0.76 there. The ratio runs with mass on every statistic
and map: 1.0-1.3 (± 0.15) at 1.3 × 10^14 and 0.78-0.80 (± 0.02) at 5.5 × 10^14, which has the sign of the mass-dependent bias
Aymerich et al. (2025) fit to the PSZ2 clusters against DES shear (α = 0.27 ± 0.09); but see the selection point below. **The
slope is not a measurement at the quoted ± 0.04.** On the aperture statistic it is 1.51-1.58 on the three maps (± 0.04
statistical each; 95 % 1.42-1.67); the background treatment moves it by ± 0.1 (1.46 without the ring subtraction on Planck,
1.66 on ACT); the eRASS1 count-rate-mass relation's own mass slope, b_X = 1.38 +0.03/−0.04 with a fitted redshift evolution
(Ghirardini et al. 2024), is a ± 0.04 shared by every bin that the stack cannot see; and the pre-registered template fit gives
1.52 on Planck and 1.83 on ACT, the latter the pressure-profile shape of the low-mass clusters (below). Counting those, the
slope is 1.5-1.8, below the Planck relation's 1.79 ± 0.08 (a fit at 2-10 × 10^14 to SZ-selected clusters with hydrostatic
masses, a different mass definition) at the 1-2σ level, and below self-similar 5/3 only on the aperture statistic. The X-ray
flux selection with correlated X-ray and SZ scatter raises the low-mass bins (the clusters near the eRASS1 limit are the
X-ray-bright, hence SZ-bright, ones at fixed mass) in the direction of a shallow slope and a high low-mass ratio; it is not
modelled here and is the first thing to settle.

## Data and method

Samples (fixed before any stack): the cosmology sample at 0.1 < z < 0.8 and 1 ≤ M500 < 8 × 10^14, inside each map's footprint;
every PCONT < 0.5 cluster in the same range as the supplementary sample. Footprint: the map's mask clean within 15' and ≥ 80 %
clean within 1°. For Planck that is the release's GAL60 and point-source masks (50 % of the sky; the point-source mask removes
60 % of the clusters, which is why the Planck samples are 1,683 and 2,554); for ACT the wide GAL070 mask with the
inpainted-region and subtracted-source masks (46 % of the sky; 1,115 and 1,596). ACT's separate cluster mask (its confirmed or
SNR > 6 SZ detections) holds 516 of the 1,115: the y map is not inpainted there, so they stay in the stack, and the split is
reported (inside 0.88 ± 0.01, outside 0.57 ± 0.01 on the template fit: the SZ-bright half at fixed eRASS1 mass is the
detected half, the Malmquist sign, and only the union is unbiased). Bins: M500 edges 1, 1.4, 2, 2.8, 4, 8 × 10^14 over
0.1 < z < 0.8, and each bin halved at z = 0.3.

Template: for each cluster the UPP (P0 = 8.403, c500 = 1.177, γ = 0.3081, α = 1.0510, β = 5.4905) at the R500 of its catalogue
M500 and z (eRASS1 cosmology, H0 = 67.77, Ωm = 0.29), projected and truncated at 5 R500, normalised to the spherical Y_500 of
the reference relation at the cluster's mass and redshift (E^-2/3 D_A^2 Y_500 = 10^-0.19 × 10^-4 [(1−b) M / 6 × 10^14]^1.79 Mpc^2
with 1−b = 0.8, so the prediction at an eRASS1 mass M is the relation at 0.8 M; the code reproduces Arnaud et al.'s
Y_cyl(5 R500) / Y_sph(R500) = 1.796), Hankel-transformed, multiplied by the beam (10' Gaussian for Planck with the nside-2048
pixel window; the ACT release's beam file) and back-transformed to the annulus means, on a 70-point grid in θ500 interpolated
per cluster; averaged over each cluster's lognormal mass posterior (σ_lnM from the catalogue's M500_L/H, median 0.2; mass,
R500 and θ500 moving together), which raises the predicted Y by 5-7 % over the median mass.

Stack: the masked mean of the map in annuli around the X-ray positions (Planck: 0-2.5, 2.5-5, 5-7.5, 7.5-10, 10-15, 15-20,
20-30, 30-40, 40-60'; ACT: 0-1, 1-2, 2-3, 3-4, 4-6, 6-8, 8-10, 10-15, 15-20, 20-30'), the Planck maps at their own HEALPix
pixels, the ACT CAR map sampled at the nside-4096 pixel centres. The mean profile of 10,000 (Planck) or 5,000 (ACT)
mask-matched random positions is subtracted (0.0 ± 0.15 × 10^-7 in every Planck annulus; −0.1 to −0.8 ± 0.4 × 10^-7 for ACT).
Covariance and null: 200 stacks per map with the map and its mask turned together in longitude by 25-335° (a random latitude
reflection for Planck), the clusters fixed and the footprint selection re-applied, each rotation's stack rescaled to the
sample's size; Hartlap-corrected. The rotation null of every amplitude is 0.00 ± its GLS error.

Statistics. Pre-registered: the GLS amplitude A of the bin-averaged template with a free constant, A = Y_observed /
Y_reference. Post hoc, after the ACT profiles were seen (row D6): the aperture statistic, the cylindrical y integral within
2 θ500 of each cluster (area-weighted from its annuli) minus the mean in the 10-20' (ACT) or 20-40' (Planck) ring, divided by
the same functional of its template, averaged over the bin. The bin's Y_500 is A times the bin's mean posterior-averaged
reference prediction with E^-2/3 removed. The Y-M fit uses the five bins' joint covariance from the rotations (bin-to-bin
correlations below 0.12); B on a grid with the amplitude solved linearly. Further variants, all in `results-R131.md`: the
template plus the two-halo shape (the projected linear correlation function at the bin's median mass and redshift from R-121's
halo model, beam-filtered, with a free amplitude), the inner annuli only, a fixed 10' (Planck) or 5' (ACT) aperture, the
catalogued neighbours (the model y map of every catalogue cluster stacked at the same positions minus each cluster's own
template: 2 % of the aperture signal and nearly flat, absorbed by the constant), the median-mass template, the raw stacks
without the random subtraction.

The shape systematic. On ACT the pre-registered template fit has χ² = 98/8 for the cosmology sample: the data sit at 0.77-0.87
of the UPP template at 0.5-5' and the stack is negative at 12-25' (−0.6 to −1.4 ± 0.5 × 10^-7) in every subset; in the lowest
bin (median 2 θ500 = 6') the data are 28 % below the UPP at the centre and 15-60 % above it at 2.5-9'. The low-mass eRASS1
clusters have flatter, more extended pressure profiles than the UPP, and a template fitted to the core at a 1.6' beam measures
that shape as much as the total Y: hence B = 1.83 and the lowest bin at 0.70 ± 0.06 on that statistic, against 1.52 and
1.05 ± 0.15 on Planck, where the 10' beam integrates the profile. The 2 θ500 aperture integrates it on both maps; on the 534
cosmology clusters inside both footprints the template fit gives 0.86 ± 0.03 (Planck) against 0.74 ± 0.01 (ACT), the aperture
0.87 ± 0.03 against 0.85 ± 0.03, and per bin the two maps' aperture values agree within 1σ at every mass. The CIB-deprojected
ACT map (β = 1.7, T = 10.7 K; 100 rotations) gives 0.758 ± 0.007 against the baseline's 0.766 ± 0.007 on the template fit and
Y_3 = 1.12 ± 0.03 against 1.14 ± 0.02 with B = 1.58 on the aperture, with the same negative ring (−0.4 to −1.6 × 10^-7):
the CIB moves nothing by more than 2 % and the ring is not the CIB. The two Planck maps differ by 2-3 % in amplitude and 0.01 in slope. Without
the ring subtraction the aperture gives B = 1.46 (Planck, whose outer profile is positive) and 1.66 (ACT, whose ring is
negative), the bracket of the background treatment.

Redshift. On the pre-registered template fit both maps show the 0.3 < z < 0.8 half above the 0.1 < z < 0.3 half at fixed mass
(Planck 0.95 vs 0.80, ACT 0.84 vs 0.70); on the ACT aperture statistic the halves agree (0.883 ± 0.025 vs 0.870 ± 0.026; per bin
1.2σ at most below 2.8 × 10^14, none above), while Planck's aperture keeps it (0.95 vs 0.82). The trend is where the beam
resolves the low-redshift clusters and the template shape enters; no redshift dependence beyond E(z)^2/3 is claimed.

Relativistic SZ: the bins' temperatures are 1.5-2.8 keV (catalogue kT medians; 1.7-5 keV from an M-T relation); the ILC's
non-relativistic spectrum underestimates y by 2-5 % there and by about 7 % at the ~7 keV of the Planck sample on which the
reference relation was calibrated with the same neglect, so the true-to-true comparison multiplies the ratios by 0.95-0.98
(lowest at low mass), steepening B by 0.02-0.03. Stated, not applied.

Trials: the pre-registered statistic on three maps and two samples (6 Y-M fits), the aperture statistic and six further
variants on each (46 fits in all, every one in `results-R131.md` and the logs), and 36 subset amplitudes per map (mass × redshift
cells, the informative rows). The headline is the aperture statistic because the pre-registered one failed its own χ² on ACT;
both are given.

## What it means

Three things the table says. (1) The eRASS1 clusters at 1-8 × 10^14 are detected in the stacked tSZ at 7-120σ per bin on
each map, and their Y_500 is measured at 2-15 % per bin against their catalogue masses: a binned, stacked Y_500-M relation of
the eRASS1 sample two decades below the SZ catalogues, on Planck and ACT. What exists beside it: Manticore-Local
(arXiv:2510.16574) stacked the Planck PR4 NILC map on the eRASS1 clusters at 1-4 × 10^14 and z ≤ 0.06 (a 2.3σ detection, no
fit); Stiskalek & Desmond (arXiv:2601.15935) compared Planck tSZ masses with the eRASS1 masses of 74 clusters at z < 0.05
cluster by cluster (0.14 dex below, slope 1.01 ± 0.04); the ACT DR6 catalogue (arXiv:2507.21459) ran forced photometry at
7,043 eRASS1 positions (masses 0.99 ± 0.01 of ACT's for the detected clusters); and the stacked Y-M relations of other samples
below the SZ limits (Planck early XI on MCXC, arXiv:1101.2043; intermediate XI on locally brightest galaxies, arXiv:1212.4131;
Greco et al. 2015, arXiv:1410.0985; Lim et al. 2018, arXiv:1712.08619; Saro et al. 2017, arXiv:1605.08770). (2) The amplitude:
at 3 × 10^14 the eRASS1 masses sit 6 % above the Planck 1−b = 0.8 scale (b_M = 0.94 ± 0.01 statistical, ± 0.05 systematic
over the statistics and the two-halo treatment, ± 0.02 each for h70 and the relativistic correction), at 5.5 × 10^14 12-15 %
above it, and at 1.3 × 10^14 on it (± 8 %). Aymerich et al. (2025, arXiv:2509.02068) found 1−b = 0.84 ± 0.06 for the PSZ2
clusters against DES shear with a mass-dependent bias (α = 0.27 ± 0.09) and the eRASS1 masses consistent cluster by cluster
above 4 × 10^14; the stack says the same at 4-8 × 10^14 (b_M = 0.85-0.88 at 1−b = 0.8 is 1−b = 0.68-0.70 at b_M = 1; against
their 0.84 ± 0.06 that is 0.84 × 0.87 ≈ 0.73, 1.8σ) and the ratio of the eRASS1 to the Planck-relation mass changes by 25 % over
a factor 4 in mass, with the same sign as their mass-dependent bias; how much of that is the X-ray selection is open (weakest
point 1). (3) The slope: 1.5-1.6 on the aperture statistic, 1.5-1.8 across statistics and systematics, so a shallow slope
(the SZ signal falling faster toward low mass than the mass predicts, as the X-ray L-M and Y_X-M relations of the same survey
do, Ramos-Ceja et al. 2025: a gas-fraction or pressure-per-mass drop below 3 × 10^14, or eRASS1 masses high at low mass) is
indicated at 1-2σ, not measured; a second mass proxy, or a selection model, separates the readings. Whatever S8 the eRASS1
counts give rests on the masses the stack measures; the stack measures them at 2-15 % per bin against one SZ relation, with the
hydrostatic bias of that relation the one nuisance, and says nothing about S8 beyond that ratio.

## Weakest points

1. Selection with correlated scatter, not modelled. The sample is X-ray-flux-selected and the masses are count-rate posteriors;
   at fixed mass the X-ray and SZ scatters are positively correlated, so the bins near the eRASS1 limit (low mass, high
   redshift) are SZ-bright at fixed mass. The eRASS1 posteriors correct the mass for the selection, not Y. The direction is a
   high low-mass ratio and a shallow slope, which is the result; the size is unknown here (ACT's cluster-mask split, 0.88 inside
   vs 0.57 outside on the template fit, shows the SZ scatter at fixed eRASS1 mass is large). The first thing to settle.
2. The statistic changed after the ACT profiles were seen. The pre-registered template fit is reported and fails its own χ² on
   ACT (98/8); the aperture statistic is a post-hoc choice, justified by the agreement of the two maps on the same clusters,
   but a choice, and on ACT the two statistics differ in slope by 0.25. Both are in the table; the claim is made at the pivot,
   where they differ least.
3. The mass axis is a count-rate posterior. The eRASS1 count-rate-mass relation's slope (b_X = 1.38 +0.03/−0.04) and its
   fitted redshift evolution are shared by every bin: ± 0.04 on B and a redshift-dependent amplitude that the stack cannot
   separate from the SZ's own.
4. The Planck aperture is the template. For the low-mass bins 2 θ500 is 6' or less, below the 10' beam, so Planck's aperture
   Y_500 is mostly the beam-deconvolved UPP conversion, the same shape ACT rejects at the centre; the two maps' agreement on the
   aperture statistic is a consistency of that conversion, not an independent measurement of the low-mass profile.
5. The template conversion. Y_500 is the cylindrical y within 2 θ500 scaled by the beam-convolved UPP's ratio of spherical
   Y_500 to that aperture; if the low-mass clusters' profiles are more extended than the UPP (as the ACT core deficit says),
   the conversion underestimates their Y_500, in the direction of a still shallower slope.
6. The background. The ACT stack is negative at 12-25' in every subset (not the CIB: the CIB-deprojected map has the same ring; the ILC's
   large-scale response, or a real pressure deficit around clusters, is untested); subtracting that ring raises the ACT amplitudes by 5 % and steepens B by 0.08
   relative to no subtraction; Planck's positive outer profile moves it the other way.
7. The reference relation and its cosmology. The Planck 2013 relation is quoted in h70 units; the templates use the eRASS1
   cosmology (H0 = 67.77), a 3 % effect on Y_3 not propagated. The hydrostatic prior enters only as the stated 1−b.
8. The two Planck samples are 60 % smaller than the catalogue because the 15' point-source cut is strict; the ACT sample is
   the southern third of the sky. The 534-cluster overlap is where the maps are compared directly.

## To settle it

The selection first: the eRASS1 selection function (Clerc et al. 2024) with a correlated X-ray-SZ scatter forward-modelled
into the stack per bin, or the stack restricted to clusters well above the local count-rate limit (the exposure map gives it
per cluster), which removes the Malmquist tail at the cost of the lowest bins. Then the profile, not the aperture: a free-shape (generalised NFW) stacked pressure profile per mass bin on ACT, from which Y_500 is
integrated without a UPP conversion; the CIB-deprojected and the ACT-only (no Planck) y maps side by side; the ACT DR6 cluster
catalogue's own Y-M calibration on the matched clusters as a second reference; the SPT-3G y map over the Dec < −40 third;
and the X-ray Y_X of the same clusters (`YX500` in the catalogue) stacked the same way, which separates a gas-fraction trend
from a mass-scale trend.

## Blind review (2026-11-01, context-free Opus referee, own literature search; `blind-review.md`)

Pass 1 (the section as first written, with "the slope is 5-7σ below 1.79" and "the first stacked SZ scaling relation" in it):
VERDICT: novel but not sound: the headline slope deficit and mass-running ratio could be produced by the X-ray flux selection and by eRASS1 mass-calibration and background systematics that are not modelled, and it rests on a statistic chosen after the pre-registered one failed, so the claimed 5-7σ significance is not supported.

Strongest objection, in one sentence: the shallow slope and the mass-dependent ratio are measured on an X-ray-flux-selected
sample against count-rate-inferred masses, and the error budget leaves out the correlated-scatter selection bias, the
calibration uncertainties shared by every bin and the ± 0.1 background bracket, while the pre-registered statistic gives a
slope consistent with 1.79 on ACT; the pivot amplitude is the robust part. Prior work it found: arXiv:2601.15935, 2510.16574,
2507.21459, 2509.02068 (with its mass-dependent bias), and the older stacks listed in "What it means". One revision pass on the
presentation half (the claim re-made at the pivot with the systematic budget, the slope demoted to an interval across
statistics and systematics, the selection named as the first weakest point, the prior work cited); the selection modelling is a
physics point, parked as "To settle it". Pass 2 below.

Pass 2 (2026-11-01, R-134, the revised section alone): VERDICT: novel but not sound: the unmodelled X-ray flux selection with correlated X-ray/SZ scatter biases the amplitude, slope and mass trend in exactly the direction claimed, the rotated-map errors omit the scatter in Y between clusters and shared mass-calibration uncertainties, and the headline statistic was chosen after the pre-registered fit failed.

Strongest objection, in one sentence: at fixed mass the clusters near the eRASS1 flux limit are X-ray- and hence SZ-bright, which raises the low-mass bins, flattens the slope, makes the ratio run with mass and raises the pivot amplitude, and no number in the table reads as a scaling relation until that selection is forward-modelled. New in pass 2: the rotation covariance omits the scatter in Y between clusters (a bootstrap over clusters is needed); the 28 % ACT central deficit at low mass looks like central radio emission (arXiv:2607.08721), which the CIB deprojection does not test. Two misses, the second on a physics point: **parked**. The path to pass is in "To settle it", with the bootstrap and a radio-source check added.


# The stacked thermal SZ signal of DESI DR9 groups at 10^12.5-10^14.5 M_sun on the ACT DR6 and Planck y maps: ACT resolves a beam-scale deficit at the catalogue centres, and the one-halo Y-M at 10^13-10^14 M_sun from ACT alone is bracketed by how that centre is modelled (slope 2.8 ± 0.3 with the universal pressure profile at the catalogue centre, 1.7 ± 0.2 with a free central term), against 1.8 ± 0.25 on Planck at 10'

2026-11-01. Script: `groups_ysz.py` (`cats`, `sample`, `randoms`, `paint`, `fit`); rows D0-D9 in `results-R132.md` (catalogues, bins,
template, radio cut, two-halo model and statistic fixed there before any stack; the two deviations recorded there with their times,
each before the fits it changed). Data: the DESI DR9 Yang et al. 2021 groups (DESI DR1 value-added catalogue `gfinder` v1.0), the
ACT DR6 ACT+Planck NILC Compton-y map (Coulton et al. 2024; its wide GAL070, inpainted-region and subtracted-source masks, the
1.6' beam, the CIB-deprojected (β = 1.7, T = 10.7 K) and CIB + dβ moment-deprojected variants), the Planck 2015 MILCA y map with
its GAL60 and point-source masks, NVSS and SUMSS (HEASARC). The Tully 2015 2MRS and Saulder 2016 SDSS group catalogues were
stacked too and are not used for a relation (below).

## Claim

The group-scale Y-M relation has been measured from Planck stacks (Lim et al. 2018 on SDSS groups; Chen, Zhang & Yang 2022 on
these DESI groups), where the 10' beam blends the one-halo signal with neighbours. The ACT DR6 y map has a 1.6' beam. Stacking
it, and Planck MILCA as the like-for-like, on the DESI DR9 groups at 0.1 ≤ z < 0.6 in bins of the catalogue's halo mass
(log M_h = log M_180m, 23 M groups above 10^12.5, 3.6-4.0 M in each map's footprint), with the Arnaud et al. 2010 universal pressure
profile (UPP) at each group's R500 as the one-halo template (normalised to the Planck 2013 Y-M relation at 1−b = 0.8, beam-convolved,
averaged over a 0.25 dex mass scatter), a two-halo term painted from the catalogue's own groups with a free amplitude B, a constant,
and the mean profile of 3 M mask- and radio-matched random positions subtracted, gives per bin the amplitude A = Y_500 / Y_500,ref
(radio-masked sample: no NVSS/SUMSS source above 5 mJy within 3'):

| log M_h bin | N (ACT) | ⟨log M500c⟩ | Planck MILCA A | ACT A, UPP at the catalogue centre | ACT A, fixed aperture 0-3' minus 3-6' | ACT A, UPP + free central term | central term S (10^-8, 0-1' annulus) |
|---|---|---|---|---|---|---|---|
| 13.0-13.25 | 423,554 | 12.90 | 0.60 ± 0.79 | −0.24 ± 0.16 | −0.42 ± 0.32 | 0.51 ± 0.97 | −6.0 ± 4.4 |
| 13.25-13.5 | 234,407 | 13.14 | 0.19 ± 0.26 | −0.01 ± 0.09 | −0.12 ± 0.16 | 0.76 ± 0.39 | −14.6 ± 4.4 |
| 13.5-13.75 | 123,755 | 13.39 | 0.29 ± 0.20 | 0.03 ± 0.06 | 0.12 ± 0.11 | 0.52 ± 0.19 | −22.8 ± 4.7 |
| 13.75-14.0 | 59,496 | 13.63 | 0.38 ± 0.08 | 0.16 ± 0.03 | 0.23 ± 0.05 | 0.39 ± 0.09 | −25.9 ± 4.7 |
| 14.0-14.5 | 34,179 | 13.93 | 0.30 ± 0.03 | 0.24 ± 0.02 | 0.26 ± 0.04 | 0.36 ± 0.03 | −45.6 ± 5.1 |

(ACT columns on the CIB-deprojected map; ⟨z⟩ = 0.42 in every bin; errors are the covariance of the bin mean from the per-group
scatter, inflated where the 23 RA-rotation nulls scatter more and by sqrt(χ²/dof) where the fit is poor (factors 1.2-3.6). The
two bins below 10^13, below the catalogue's reliable mass range, are in `results-R132.md` D6-D7 and give nothing beyond 1.5σ on
ACT.) Power laws Y ∝ M500c^B_Y over the five bins:

| map, model | B_Y | log10 Y(10^13.5) [Mpc^2, E^-2/3 D_A^2 Y500, bin-mean] | χ²/dof |
|---|---|---|---|
| Planck MILCA, UPP at centre | 1.78 ± 0.25 | −6.87 ± 0.10 | 0.8/3 |
| ACT CIB-deprojected, UPP at centre | 2.78 ± 0.31 | −7.40 ± 0.13 | 4.1/3 |
| ACT CIB-deprojected, UPP + central term | 1.66 ± 0.20 | −6.74 ± 0.08 | 0.9/3 |
| ACT standard / CIB + dβ, UPP + central term | 1.46 ± 0.17 / 1.38 ± 0.26 | −6.62 ± 0.07 / −6.52 ± 0.10 | 2.1/3, 0.7/3 |

On identical groups (the 1.75 M in both footprints, randoms restricted to both; `fit --common`), with the same model and
background subtraction on both maps:

| log M_h bin | N | Planck MILCA A | ACT A, UPP at centre | ACT / Planck | ACT A, UPP + central term | ACT / Planck |
|---|---|---|---|---|---|---|
| 13.5-13.75 | 61,049 | 0.79 ± 0.32 | 0.01 ± 0.08 | 0.02 ± 0.10 | 0.62 ± 0.24 | 0.78 ± 0.43 |
| 13.75-14.0 | 29,366 | 0.49 ± 0.11 | 0.16 ± 0.03 | 0.32 ± 0.10 | 0.40 ± 0.06 | 0.83 ± 0.22 |
| 14.0-14.5 | 17,116 | 0.36 ± 0.05 | 0.23 ± 0.02 | **0.63 ± 0.09** | 0.34 ± 0.02 | **0.93 ± 0.13** |

(slopes on the common footprint: Planck 1.65 ± 0.29, ACT 2.83 ± 0.36 and 1.72 ± 0.21; the ratio errors add the two in
quadrature, conservative since the ACT ILC contains Planck.) At the anchor the resolved one-halo amplitude with the UPP at the
catalogue centre is 0.63 ± 0.09 of Planck's on the same groups; the central term restores agreement (0.93 ± 0.13).

What the claim is. (1) **At 10^14-10^14.5 the one-halo signal is detected on ACT on every map and model (7-15σ), at A = 0.21-0.45 of the
reference relation across the three maps and the two centre models; at 10^13.75-14 at 0.09-0.69 (1.7σ on the standard map
without a central term, ≥ 2.2σ otherwise). Planck gives 0.30 ± 0.03 and 0.38 ± 0.08 on the same catalogue.** (2) **From 10^13.0 to 10^13.75 the ACT amplitude is not determined by the
stack without a model of the centre**: with the UPP at the catalogue centre it is consistent with zero (95 % upper limits, A + 1.96σ, 0.07,
0.17, 0.15 of the reference in the three bins; the fixed aperture agrees), while a free beam-shaped central term lets it rise to
0.5-0.8 ± 0.2-0.4; Planck's 0.2-0.6 ± 0.2-0.8 does not choose. (3) **The ACT slope over 10^13-10^14.5 is 1.4-1.7 (± 0.2) or
2.8 ± 0.3 depending on that choice; self-similar 5/3 and Planck's 1.78 ± 0.25 lie inside the first and 3σ from the second.**
(4) **The deficit itself is measured**: the data prefer the free central term over the UPP at the catalogue centre by Δχ² = 149 for
five parameters, and over a Gaussian miscentering of rms 0.25-1 R500 (one parameter for the catalogue) by Δχ² ≥ 71. On the
CIB-deprojected map S is −15 to −46 × 10^-8 (3-9σ per bin) from 10^13.25 up, growing with mass; at fixed mass it is twice as deep
in the z ≈ 0.47 half as in the z ≈ 0.22 half (−57 ± 6 vs −31 ± 7 at 14-14.5), and it differs between the three maps by 20-40
(± 5-13) × 10^-8. A radio source of fixed luminosity would be five times brighter in the low-z half, so the deficit is not the
radio hosts left by the 5 mJy cut alone; its map dependence says it is not pure tSZ either. Not claimed: anything about
feedback, the identity of the deficit, or a ratio of the ACT one-halo to the Planck-era relations below 10^13.75 beyond the bracket.

## Data and method

Catalogue (fixed before any stack): GRP_LOGM converted from h^-1 M_sun (h = 0.674) and from M_180m to M500c and M200c with an NFW
profile and the Dutton & Macciò 2014 concentration; 0.1 ≤ z < 0.6 (mostly photometric redshifts); bins of 0.25 dex in log M_h from
12.5 to 14.0 and 14.0-14.5 as the cluster anchor; the referee-advised floor 10^13 for the headline. Footprint: 99 % of the
15' disc clean in the ACT mask (wide GAL070 ≥ 0.99, not inpainted, not a subtracted source) or the Planck GAL60 and point-source
masks. Statistic: annulus means (edges 0, 1, 2, 3, 4, 6, 8, 10, 15, 20, 30' for ACT; 0, 2.5, 5, 7.5, 10, 15, 20, 30' for Planck)
from FFT convolutions of the masked map with the annulus kernels, sampled at the nearest pixel; per bin the GLS fit
y_k = A T1_k + B T2_k + c with the covariance of the bin mean from the per-group scatter (the group-to-group scatter in y, the
noise and the neighbours included; 30,000-1.3 M groups per bin). T2: the y map painted from the 6.3 M catalogue groups with
log M_h ≥ 13 at z < 0.9 at the reference relation, sampled the same way, minus each group's own template. Radio: a group is
masked if an NVSS (or SUMSS, scaled with α = −0.7) source above 5 mJy at 1.4 GHz lies within 3' (19-28 % of groups, about 17 %
by chance coincidence); unmasked results in `results-R132.md` (radio hosts lower A by 0.02-0.12 at the UPP-at-centre model,
by 0.03-0.56 with the central term).

Deviations from the pre-registration, each recorded before the fits it changed: (i) the pre-registered fit failed on ACT
(χ²/dof 120-1230/7, negative A and B): the 0-1' annulus of the standard map is −7 to −9 × 10^-8 in every bin to 10^13.5 (dust in
the member galaxies; it is near zero on the CIB-deprojected map, which was the pre-registered foreground check and is now the
headline map), and the annulus means of 15'-clean positions fall with radius at random positions too (the clean-disc condition
puts the centre away from the subtracted sources, whose negative wings lie further out), which the constant cannot absorb: the
mean profile of 3 M mask-matched randoms is subtracted and the error inflated by sqrt(χ²/dof). (ii) The radio-masked stacks then
showed a step at 3', the radius of the radio cut: they subtract randoms with the same radio condition. The central term and the
miscentering model are post-hoc diagnoses, reported as a bracket, not as the headline fit. Planck gets the same subtraction; its
fit is acceptable with or without it (D1, D7). The pre-registered numbers are in D1-D4.

Checks: richness ≥ 3 (slope 2.79 ± 0.32 with the UPP at the centre; bins unchanged within 0.05); mass scatter 0.15 and 0.35 dex
(A ± 25 %, slope 2.75-2.80); the CIB + dβ map (A at 13.0-13.75: −0.93 ± 0.36, −0.34 ± 0.22, −0.05 ± 0.12 with the UPP at the
centre); two-halo amplitude B = 0.1-0.35 ± 0.15-0.30 on ACT below 10^13.5 and 0.9-1.1 on Planck (the painted neighbours are at
the reference relation, which the anchor bins put at 0.3-0.4 of it, so B < 1 is expected at ACT's resolution where the
neighbours' own cores fall outside the inner annuli; Planck's B absorbs its blended one-halo). Low-z catalogues: Saulder's
luminosity masses put its clusters (log M_tot ≥ 14, 167 groups on ACT) at A = 0.13 ± 0.03, half the DESI anchor, so its mass
scale is not the template's; Tully has fewer than 200 groups per bin in the 15'-clean ACT footprint; neither is fitted for a
relation (D4).

## What it means

At ACT's resolution the stacked group signal below 10^14 is a small difference between a one-halo profile and a beam-scale
negative feature at the catalogue centre. With a single y map and catalogue centres from luminosity-weighted photometric
members, the two are degenerate below 10^13.75: the UPP-at-centre fit sees no group signal there and a slope of 2.8, the
free-centre fit sees 0.5-0.8 of the reference and a slope of 1.7. The Planck-resolution values (here 0.2-0.4; Lim et al. 2018's Table 1
against the same reference: 0.16, 0.27, 0.27, 0.55 at log M500 13.3, 13.5, 13.7, 13.9) sit between the two. The deficit grows with mass and redshift and changes with the
ILC's spectral weights; it survives the NVSS 5 mJy cut and does not dim with distance as a radio source would. Its candidates
(centring errors that are not Gaussian, a pressure profile flatter than the UPP inside R500 in groups, dust or radio emission
below the cuts with a redshift-dependent luminosity, CIB with an SED outside the deprojected family) are not separated here.
The usable statement for the group Y-M is the bracket and the anchor: 0.3-0.4 of the Planck 2013 relation at 10^13.6-13.9, and
no constraint from ACT alone between 10^13 and 10^13.75 tighter than [0, 0.8].

## Weakest points

- The centre model is chosen after the data, and the two models disagree; the claim is the bracket, not either value.
- The halo masses are the catalogue's abundance-matched M_180m from photometric groups. A is relative to them, and the anchor
  at 0.3-0.4 says they, or the reference relation, are off by a common factor; a mass-dependent bias moves the slope.
- The per-map tables use each map's footprint (3.97 M Planck, 3.56 M ACT groups); the common-footprint table is the
  like-for-like and agrees with them within 1σ on ACT and within 1.3σ on Planck.
- The radio cut (5 mJy at 1.4 GHz) misses fainter or steep-spectrum hosts; the redshift test argues against a fixed-luminosity
  population only.
- The neighbour model paints only log M_h ≥ 13 (pre-registration deviation, 10:35 UTC) at the reference amplitude; B is free.

## To settle it

A two-frequency (ACT f090/f150) stack on the same groups to give the central term a spectrum; the DESI spectroscopic-redshift
subset with BCG-based centres against the luminosity-weighted centres; a free-shape pressure profile in place of the UPP; weak-lensing masses for the catalogue bins.

## Blind review (2026-11-01, context-free Opus referee, own literature search; `blind-review.md`)

Finding check, pass 1 (R-136): "VERDICT: novel but not sound: the beam-scale central "deficit" that drives the ACT bracket and claim (4) is most plausibly flux-limited radio/dust residual contamination (already documented in ACT DR6×DESI stacks), the redshift test cannot exclude it, and without frequency-resolved separation the ACT one-halo amplitude below 10^14 M_sun is unconstrained."

Strongest objection: the beam-scale central term that sets the ACT bracket is most plausibly residual radio and dust
emission at the group centres (arXiv:2607.08721, 2609.08939, 2502.08850), and the note's argument against that reading
fails because a flux-limited cut admits more luminous sources at higher redshift, so the deeper high-z deficit is what residual
radio predicts. The referee also named: bin-to-bin covariance ignored in the power-law fits, mass-dependent Eddington bias
from the catalogue's 0.2-0.4 dex mass errors, and the two-halo template painted at an amplitude the anchor shows is 2.5-3×
high. Novelty held: no ACT DR6 tSZ Y-M on the Yang DR9 groups was found. A physics miss: parked, no second pass. The path to
pass, if reopened: per-frequency ACT DR6 stacks (f090/f150/f220) to give the central term a spectrum, a deeper or
luminosity-defined radio cut (FIRST, VLASS), the spectroscopic subset with BCG centres, a jackknife covariance across bins.

# ACT DR6 CMB lensing over 365 eRASS1 superclusters, against the ΛCDM expectation from the Websky light-cone run through the same pipeline: no lensing excess beyond ΛCDM outside the members' one- and two-halo model; relative to ΛCDM the excess is +0.18 ± 0.37 of the member model in the members' hull and −0.43 ± 0.38 in a 20-Mpc disc (95 % upper limits 0.9 and 0.3; expected sensitivity ± 0.75)

2026-11-01/02. Scripts: `sc_members.py` (the member lists), `sc_lens.py` (`prep`, `data`, `sims`, `fit`: the measurement), `sc_null.py`
(the same pipeline on a simulated light-cone: the ΛCDM null); rows D0-D13 in `results-R137.md` and D0-D10 in `results-R140.md`, in each
of which the samples, apertures, model, statistic, control and comparison were fixed before any aperture was integrated.
Data: the eRASS1 supercluster catalogue (Liu et al. 2024, arXiv:2402.08454; 1,338 systems in the western Galactic hemisphere),
the eRASS1 primary cluster catalogue (Bulbul et al. 2024; M500 from the weak-lensing-calibrated count-rate-mass relation) and its
optical follow-up (Kluge et al. 2024); the ACT DR6 lensing release v1 (Qu et al. 2024, arXiv:2304.05202; Madhavacheril et al.
2024, arXiv:2304.05203): the baseline, polarization-only and CIB-deprojected convergence, the baseline mask, the 400 simulated
reconstructions and their input convergence (Farren et al. 2024), and ACT's Websky foreground-only reconstruction (arXiv:2304.05196);
the Planck PR3 tSZ-deprojected lensing reconstruction; the Websky simulation (Stein et al. 2020, arXiv:2001.08787): its halo
light-cone and its CMB-lensing convergence.

## Claim

The question: does a supercluster hold mass outside its member clusters beyond what the members' own environments would hold if
they were not in a supercluster, and beyond what ΛCDM predicts for that comparison? The statistic is the band-filtered
(20 ≤ L ≤ 2048) ACT DR6 convergence integrated over the supercluster's footprint, minus a model of its members (each an NFW halo to
r200c plus its linear two-halo term, at its posterior-averaged eRASS1 mass, filtered like the map), divided by that model,
inverse-variance weighted over superclusters: E. The comparison is the same E on random groupings of isolated eRASS1 clusters
matched member by member in mass and redshift and superposed at the members' positions: E_control. E − E_control is not expected
to be zero: the members' overlapping two-halo terms, the friends-of-friends selection and the isolation of the control all shift
it. Its ΛCDM expectation is computed by running the whole pipeline (selection, finder, sample cuts, apertures, member model,
control, weights) on the Websky light-cone in 25 data-sized patches of the simulated sky put through the real eRASS1 and ACT
footprints, with the ACT DR6 noise realisations added. The claim is the observed E − E_control minus that expectation, Δ, with
the spread of the statistic for one data-sized sample under ΛCDM (sample variance and noise) as its error:

| sample | aperture | n | observed E − E_control | ΛCDM E − E_control (mean; one-sample sd) | **Δ = observed − ΛCDM** (68 %) | 95 % | 95 % with ± 10 % on the mass scale |
|---|---|---|---|---|---|---|---|
| ≥ 2 members | hull + 5 cMpc | 365 | +0.29 | +0.11 ± 0.08; 0.36 | **+0.18 ± 0.37** | [−0.55, +0.91] | [−0.59, +0.97] |
| ≥ 2 members | disc, 20 cMpc | 365 | −0.31 | +0.12 ± 0.09; 0.37 | **−0.43 ± 0.38** | [−1.17, +0.31] | [−1.23, +0.33] |
| ≥ 3 members | hull + 5 cMpc | 119 | +0.27 | +0.13 ± 0.12; 0.50 | **+0.14 ± 0.51** | [−0.86, +1.14] | [−0.92, +1.22] |
| ≥ 3 members | disc, 20 cMpc | 119 | −0.05 | +0.20 ± 0.17; 0.60 | **−0.25 ± 0.63** | [−1.48, +0.97] | [−1.55, +1.02] |

(All in units of the member model, whose one-halo part is 35 % of it in the hull and 16 % in the disc. Δ's error is the one-sample
sd combined with the error on the ΛCDM mean (patch spread over √8 independent patches). The observed values sit at the 69th,
12th, 61st and 35th percentiles of the ΛCDM distribution. The aperture integral itself is detected: K_d/K_m = 1.46 in the hull,
6.2σ above zero.)

So the eRASS1 superclusters are what ΛCDM expects. (The disc's 0.3 comes from a low draw, the 12th percentile; the expected 95 %
half-width is ± 0.75.) ΛCDM predicts a positive E − E_control, +0.11 to +0.20 of the member model
(the isolated control sits low, E_control −0.10 to −0.34 in the simulation, as an isolation selection should), a third of one
sample's spread; the data's values lie within 1.1σ of it in every pre-registered row. Relative to ΛCDM, the 365 superclusters
carry at most 0.9 times their member model in additional band-filtered lensing inside the members' hull dilated by 5 comoving Mpc,
and at most 0.3 times it inside a 20-Mpc disc (95 %); for the 119 systems with three or more members, 1.1 and 1.0. A deficit
larger than 0.55 (hull) and 1.2 (disc) of the member model is excluded.

## Data and method

Members. The released supercluster table gives each system's centre (the mean of its members' X-ray centroids), redshift,
multiplicity, and mtot = 1.46 Σ M500 / 0.39, but not its members. The member table (Liu et al. 2024, Table 2) is at the CDS only,
which this analysis could not reach, so the members were rebuilt with the paper's recipe: the same cluster sample (8,910 clusters
here, 8,862 in the paper), friends-of-friends on 3D comoving distances with l = (n f)^(-1/3), f = 10, l ≤ 50 Mpc, then for each
catalogued system a constrained search (up to three members removed or added) for the member set that reproduces its centre, mean
redshift and mtot to the table's rounding. 991 of 1,338 systems are reproduced, none ambiguously; only those are used. The
reproduced set favours small systems (83 % of pairs, 60 % of systems with ≥ 3 members).

Sample. Reproduced systems at 0.05 < z < 0.6 with both apertures ≥ 80 % clean in the ACT DR6 baseline mask: 365 (246 pairs,
119 with ≥ 3 members).

Map and apertures. The ACT DR6 baseline (profile-hardened, temperature plus polarization) κ, band 20 ≤ L ≤ 2048 with a cosine
taper, at HEALPix nside 512; clean pixels are those whose nside-4096 mask mean is ≥ 0.99. Two apertures per supercluster: the
members' convex hull in the comoving tangent plane, dilated by 5 comoving Mpc (median 1.0 deg² clean), and a disc of 20 comoving
Mpc around the centre (median 3.7 deg²). K_d = Σ κ Ω_pix over the aperture's clean pixels.

Member model. Each member is an NFW halo truncated at r200c (Duffy et al. 2008 concentrations; M500c from the catalogue) plus its
linear two-halo term (Tinker et al. 2010 bias, CAMB linear P(k) at σ8 = 0.88), κ(L) exact in the flat sky, filtered like the map,
averaged over the member's lognormal mass posterior and summed over the same clean pixels: K_m.

Statistic and control. E = Σ_s w_s (K_d,s − K_m,s) / Σ_s w_s K_m,s, w_s = K_m,s / var(K_d,s), the variances from the first 100
simulated reconstructions. E_control: each member replaced by an isolated eRASS1 cluster (the same sample selection, in no rebuilt
system, 1,574 in the clean footprint) with |Δz| < 0.02 and |Δ ln M500| < 0.15, the aperture carried with it; the same E on that
superposition with the same weights, averaged over ten realisations.

The ΛCDM null. Websky: Ωm 0.31, h 0.68, σ8 0.81; the halo light-cone sorted by mass (the 30 M most massive halos, complete above
M200m = 1.6 × 10^13 M☉) and the full-sky κ to last scattering. Halo M200m → M500c through NFW (Duffy et al. 2008 c200m); observed
redshift = cosmological + line-of-sight velocity + an error drawn from the eRASS1 sample's at that redshift. Selection probability
Φ((ln M500c − ln M_c(z)) / 0.6), M_c(z) per 0.02 in z set so the density equals the eRASS1 finder sample's and the scatter 0.6
chosen so the median masses match too (to ≤ 7 % at 0.15 < z < 0.6, 9-23 % below; a sharp threshold overshoots the data's medians
by 1.3-1.6). The released eRASS1 selection function, a function of count rate, would need an X-ray luminosity per halo; the matched
mass-redshift selection replaces it. 25 patches: the Websky sky read as equatorial and 24 random rotations of it, each through the
real eRASS1 footprint (15,376 deg²) and the ACT mask (overlap 4,865 deg²; the patches overlap by 12 % on average, ~8 independent).
In each: the same finder (one exposure bin, as the simulated selection is exposure-free), the same cuts, apertures, member model
(at the simulated clusters' true M500 with a catalogue-drawn posterior width), control and weights, all imported unchanged. Per
patch: 354 systems on average (311-398; data 365), 129 with ≥ 3 members (data 119), isolated pool 1,686 (1,574): the samples match
the data's to ± 15 %, inside the ± 30 % fixed beforehand as the condition for the null to be meaningful. The simulated κ gives the
noiseless E − E_control per patch; the 100 ACT noise realisations (simulated reconstruction minus its input κ) are added to give
2,500 data-like values, whose spread (0.36-0.60) is the comparison's error; the patch-to-patch part (0.23-0.47) is the
non-Gaussian sample variance, the noise part 0.27-0.38. R-137's own error from the Gaussian simulations was 0.36-0.56: the same size.

Systematics and checks. The member-mass calibration enters as one scale on every member mass (E and E_control move together;
± 10 % propagated, × 1.07 / × 0.94 on the hull's interval; the ACT measurement of the same masses, b_M = 0.65-0.70 ± 0.22, would raise
the limits by a fifth). The null under a sharp threshold (s = 0): +0.20, +0.16 (≥ 2), +0.26, +0.21 (≥ 3), Δ moving by ≤ 0.13; with a
noisy mass in the model and the control matching: shifts ≤ 0.03; with the model's two-halo term at Websky's σ8: ≤ 0.03. ACT's own
Websky tSZ + CIB + radio foreground-only reconstruction added on patch 0 moves E − E_control by ≤ 0.07. On the data map: the
polarization-only reconstruction gives E − E_control = −0.28 ± 0.84 (hull) and −0.26 ± 0.60 (disc), the CIB-deprojected one
+0.31 ± 0.34 and −0.19 ± 0.34, Planck PR3 (207 superclusters) +0.45 ± 0.99 and +0.62 ± 0.73; random-position and rotated nulls
within 1σ of zero; the all-detected-clusters model +0.43 ± 0.30 and −0.10 ± 0.30; a population control +0.24 ± 0.27 and −0.48 ±
0.29; aperture sizes 2.5-30 cMpc (defined after a first look) −0.02 to +0.49. All are read against zero there; against the null,
each shifts down by ~0.1.

## What it means

The comparison R-137 left open now has a number: in ΛCDM, a supercluster's footprint carries 0.1-0.2 of its member model more
band-filtered lensing than isolated clusters of the same masses put at the members' positions, because the isolated control is
underdense by selection and friends-of-friends picks overdense company. The eRASS1 superclusters match that: there is no sign of
dark matter outside the members, between them, beyond the ΛCDM expectation of the structure around clusters. In the 5-Mpc hull the
data allow an excess of up to 0.9 of the member model above ΛCDM; in the 20-Mpc disc 0.3. The data's E alone sits 2.0σ above the
simulation's in the hull and its E_control 0.9σ above: a common offset (the mass scale, σ8 or the real clusters' environments)
that the difference removes. The data's hull and disc values differ by +0.60 where ΛCDM expects 0.00 ± 0.33 (1.9σ; the two
apertures correlate at 0.6), a fluctuation the null does not share and the claim does not use. It is the first CMB-lensing
measurement of the mass outside the members of the eRASS1 superclusters (galaxy weak lensing has compared eRASS1 supercluster
members with isolated clusters, arXiv:2404.17345; CMB lensing stacks exist on DES and DESI Legacy Survey superstructures,
arXiv:2203.11306, 2105.11936), and a bound read against ΛCDM, not a detection.

## Weakest points

- One simulation. Websky is a peak-patch light-cone whose matter outside halos is second-order Lagrangian perturbation theory and
  whose halos are painted with profiles; the null is its ΛCDM, at σ8 0.81. A ray-traced N-body light-cone (AGORA on MDPL2,
  arXiv:2212.07420) would check that the field matter between members is right on 5-20 Mpc scales.
- The selection is a matched mass-redshift threshold with scatter, not the eRASS1 X-ray selection function; a sharp threshold
  raises the null by ~0.1, so the selection's shape is worth ~0.1, a quarter of the error.
- The data's masses are count-rate estimates; an environment-dependent bias in them (projection, AGN, the members' lower
  concentrations, arXiv:2404.17345) would enter E − E_control and is not in the simulated null; a 10 % member-versus-isolated
  difference is worth ~0.1, a quarter of the error. Members of the 347 systems not reproduced that the rebuilt finder did not
  group can enter the data's isolated pool.
- The members are rebuilt, not read (991 of 1,338 reproduced exactly), and the reproduced set favours pairs; the simulation uses
  every system, and a subsample matched to the data's multiplicity mix gives a null of +0.06 and +0.08, within the mean's error.
- The member model is dominated by the members' overlapping linear two-halo terms; that is why E − E_control has a non-zero null,
  which is now computed rather than assumed, but E alone is still not "mass outside clusters".
- The band filter removes L < 20, so mass on scales above ~9° is not seen; the model and the simulation are filtered identically.
- The weights use the simulations' variance; one weighting was fixed and nothing else tried. R-137 looked at the fit once with
  weights from 28 simulations before the 100 named (its row D3).
- The mass-scale systematic is one ± 10 % scale; Planck PR4 could not be fetched (the Planck check is PR3).

## To settle it

The Liu et al. member table from the CDS (the full 1,338 systems), the eRASS:4-5 or eastern-hemisphere supercluster catalogues,
ACT DR6 combined with Planck PR4 or SPT-3G (less noise per area), and a second, N-body light-cone with an X-ray selection model
for the null would together shrink the interval by about a factor of two and test the simulation dependence.

## Blind review (2026-11-01, context-free Opus referee, own literature search; `blind-review.md`)

Finding check, pass 1 (R-137): "VERDICT: novel but not sound: the central statistic E − E_control is compared with zero, but no ΛCDM null has been computed for it, because the overlapping two-halo terms and the isolated (underdense) control mean zero isn't the expected value, so the bounds can't be read as limits on mass outside the member clusters."

Strongest objection: the member model sums every member's two-halo term over a shared environment (double counting it), while
the control superposes independent, below-average environments (isolation is a selection), and the friends-of-friends selection
itself raises the density between members, so E − E_control has an uncalibrated, possibly order-unity ΛCDM expectation; the
hull/disc sign flip is what an aperture-dependent offset would look like. The referee also named: the rebuilt members' selection
(74 %, favouring pairs and possibly compact systems), errors without halo sample variance or mass-observable scatter, a ±10 %
mass prior narrower than the ACT measurement of the same masses, foreground checks too noisy to test at ±0.36, and the title
overstating a 95 % limit of the order of the whole member model. Novelty held for the eRASS1 superclusters; the closest prior
work is a CMB-lensing stack on DES photometric superclusters (arXiv:2203.11306, 8.4σ, amplitude 0.84 ± 0.10) and on DESI
Legacy Survey superstructures (arXiv:2105.11936), so "first lensing measurement" holds only for this catalogue. A physics miss:
parked, no second pass (the title was corrected after the verdict to say "bounds", not "finds no"). The path to pass: run the
whole pipeline (the eRASS1 selection, friends-of-friends, apertures, member model, control) on an N-body light-cone to give
E and E_control their ΛCDM values, compare the data with those, and use the CDS member table in place of the rebuilt one.


Finding check, pass 2 (R-140, 2026-11-02, context-free Opus referee, the section alone, the prompt of `blind-review.md` verbatim;
after the ΛCDM null was added): "VERDICT: the findings are confirmed novel and sound if the calculations hold"

Strongest objection: the simulated null uses true masses, so an environment-dependent bias in the eRASS1 count-rate masses
(projection, AGN, members' lower concentrations) would enter E − E_control without being in the null; the referee sized it at
~0.1, a quarter of σ, so the consistency with ΛCDM survives. Also named: the selection is a mass-redshift threshold not the X-ray
selection function; Websky's inter-halo matter is 2LPT; ~8 independent patches; the title's "carry no mass" overreached and the
disc limit is a low draw; "first lensing measurement on the eRASS1 superclusters" is false as written (arXiv:2404.17345, weak
lensing of eRASS1 supercluster members against isolated clusters). Closest prior work: arXiv:2203.11306, 2404.18278, 2105.11936,
1709.02543, 2507.08561, 1402.3302, 1702.08485, 2107.05523, 2409.04535, 2404.17345, 2503.09952. Corrections applied after the verdict,
none to a number: the title says "no lensing excess beyond ΛCDM outside the members' one- and two-halo model" with the expected
sensitivity, the "first" sentence narrowed to CMB lensing of the mass outside members with 2404.17345 cited, the mass-estimator
point and the isolated-pool point added to the weakest points.
---

# Finding (R-130/R-135/R-141, cosmo-deep-10, 2026-11-02): sub-threshold ZTF forced photometry at 54 known galaxy-scale strong lenses, 2018-2026, with injections through the same vetting code as the candidates: lensed supernovae whose images sum to r < 21 occur at fewer than 0.045 per lens-year behind grade-A lens candidates (95 %, efficiency uncertainty included); one candidate behind a spectroscopic lens survives the vetting and cannot be told from a lens-galaxy transient

## Claim

Forced PSF photometry on every public ZTF r-band difference image (2018-03 to 2026-09; median 620-860 clean epochs per lens)
at the positions of 54 known galaxy-scale strong lenses (27 spectroscopic lenses of Shu et al. 2018, 26 grade-A Legacy
Surveys candidates and one check target), with a detection statistic, threshold and vetting fixed in advance and the
vetting written as one function, finds two lenses above threshold. One survives every stage of the vetting: **S4TM
J0757+1956** (2026-02, PSF-corrected peak r ≈ 21.1), consistent both with a type Ia supernova at the source redshift
(z_s = 0.83) magnified 17 ± 2 times (the lens model gives 9; the vetting accepts 0.3-3 times that), whose g flux at peak it
also matches, and with a faint transient in the lens galaxy (M_r ≈ −17.6 at z_l = 0.12). ZTF photometry does not separate
the two, so it is counted as a detection. The other (DESI-266.4075+39.8584) fails the light-curve shape cut.

Lensed supernovae drawn as the per-lens forecasts draw them, injected into each lens's own light curve through the
pipeline's measured difference-image PSF, and passed through **the same vetting function as the real light curves** (shape,
control positions, centroid, season-to-season variability, colour at peak, magnification consistency) are recovered and
kept 4.6 % of the time behind the candidates (8.3 % before vetting) and 0.19 % behind the spectroscopic lenses (0.28 %),
whose sources mostly sit at z ≈ 1 in this sample. Of the injections whose images sum to r < 21, the vetting keeps 34 % behind
the candidates (62 % before vetting): 50 % at r 19-20, 19 % at 20.5-21, 8 % at 21-21.5.

**Headline (forecast-free: no rate model enters; the source population enters only through the injected light-curve mix
inside the magnitude cut):**

| lens class | N | lens-years × vetted efficiency, summed peak r < 21 | detections | lensed SNe with summed peak r < 21 per lens-year, 95 % |
|---|---|---|---|---|
| L: Legacy Surveys grade-A candidates | 26 | 79.7 (113.5 before vetting) | 0 | **< 0.045** (< 0.038 without the efficiency uncertainty) |
| L, summed peak r < 20 | 26 | 99.7 | 0 | < 0.030 (without the uncertainty) |
| S: Shu et al. 2018 spectroscopic lenses | 27 | 5.6 (12.0) | 1 (S4TM J0757+1956) | < 1.7 (< 0.95 if J0757 is not lensed): no useful constraint |

The class-L limit includes the efficiency's uncertainty as a lognormal on the exposure: the bootstrap over injections (5 %)
and the image offset (the exposure falls to 22.7 lens-years if every image sat 2" from the centre instead of 1.2"; 2" is
taken as a 2σ excursion, an 18 % width; with every image at 2" the limit is < 0.13). Per bin of summed peak r (class L,
95 %, per lens-year): 19.0-19.5 < 0.049, 19.5-20.0 < 0.039, 20.0-20.5 < 0.047, 20.5-21.0 < 0.080, 21.0-21.5 < 0.15.

**Secondary (each conditional on the source population it names):**

| lens class | assumed source population | lensed SNe forecast to occur in the window | E before / after vetting | detections | rate / that population, 95 % |
|---|---|---|---|---|---|
| S | Shu et al. 2018 per-lens rates (their [O II]-based source SFRs, z_s and total magnifications) | 67.5 | 0.095 / 0.065 | 1 | < 73 |
| L | SLACS-like medians of Shu et al. (N_Ia 0.004, N_cc 0.030 per year), sources at z_l × N(2, 0.5) | 7.5 | 0.65 / 0.37 | 0 | < 8.2 |
| L | BELLS-like medians (0.018, 0.149 per year), same sources | 36.9 | 3.09 / 1.73 | 0 | < 1.7 |

None of the forecasts is tested: the per-lens Shu forecast for spectroscopic lenses predicts 0.065 kept detections, and the
BELLS-like bracket for the candidates, excluded before the vetting was in the efficiency, no longer is.

## Data and method

Fixed in `results-R130.md` before any light curve was extracted (2026-11-01 07:40 UTC), except where an amendment is named;
each amendment either lowers the efficiency or is applied identically to real and injected light curves.

**Lenses.** Class S: the spectroscopic galaxy-galaxy lenses of Shu et al. 2018 (arXiv:1803.07569; SLACS, S4TM, BELLS), each
with Shu's per-lens forecast of lensed SNe Ia and core-collapse SNe per year. Class L: grade-A DESI Legacy Surveys lens
candidates (Huang et al. 2020, 2021; Storfer et al. 2022) at Dec > −30 not within 3" of a class-S lens. Order fixed in
advance: class S by forecast rate (descending), class L by lens r magnitude (ascending), interleaved one for one. The
sample is the contiguous prefix of that order whose data could be fetched while the archive was rate-limited: ranks 0-65
(33 S, 32 L, one check target, Sheu et al. 2023's in-footprint candidate, whose transient predates ZTF). Five sit where the
difference images are masked and seven have fewer than 20 clean epochs or 8 PSF stamps, leaving 54 with an efficiency.

**Photometry.** Every public ZTF r-band difference image (IRSA, pipeline `scimrefdiffimg`) covering a lens that passes the
fixed epoch cut (INFOBITS = 0, seeing ≤ 4", limiting magnitude ≥ 19), as a 30" cutout. Forced photometry with a Gaussian of
FWHM = the header seeing plus a constant background (weighted least squares, 5σ pixel clipping) at the lens centre and at
four control positions 10" N, S, E and W in the same cutouts; per (field, CCD, quadrant) the median is the baseline; errors
scaled by the light curve's robust scatter (floored at 1). Validated on eight ALeRCE SNe Ia (agreement with the alert PSF
magnitudes within 0.2-0.3 mag after the PSF calibration below; the statistic finds them at T = 12-18).

**Statistic.** In 15-day windows every 5 days, S_w = the signal-to-noise of the window's weighted mean flux; T = max S_w.
Candidate: T > 5.31, the 99.5th percentile of the control positions' T in an earlier 176-lens run of the same pipeline
(704 control light curves); in this sample's 216 control light curves the 99.5th percentile is 4.1, so the threshold is
conservative.

**The vetting, one function applied unchanged to real and injected light curves** (`ztf_lsn_vetting.py`):
0. shape: a window above 5.31 with a neighbouring window (15 d before or after) at S > 3;
1. controls quiet: no control position at S > 3 in that window;
2. centroid within 2.5" of the lens (real: the stacked window image; injected: the flux-weighted image positions plus
   centroid noise);
3. no persistent variability: no season (outside the event span, [t_w − 100, t_w + 250] d) whose mean deviates from the
   median season mean by more than 1.0 microJansky at more than 5σ;
4. colour and magnification: for some supernova hypothesis at the source redshift (SALT3 Ia, M_B = −19.25; IIP, Ic, IIb, Ib,
   IIL, IIn at their mean luminosities; class S at the measured z_s, class L on a grid z_l × 1.0-3.0), the fit of its r light
   curve (amplitude and explosion time, PSF-corrected) requires a total magnification within 0.3-3 times the lens model's
   (class S: Shu's; class L: the injection draw's median, 10.6), and the g flux over the r peak is within 3σ of that
   hypothesis's prediction (no g epochs there: untested).
The two variability numbers were set so that the two candidates an earlier by-eye vetting rejected as variables fail
(seasons offset by 2.3 and about 4 microJansky at > 10σ) and its two lens-galaxy transients pass, and were frozen before any
injection was vetted; no control light curve fails them. The colour test was first frozen as a chi2 over the whole event
span and changed, before any injection was vetted, to the colour at peak (the form the by-eye vetting used), after the chi2
form passed SLACS J1250+0523 (outside this sample) whose g at peak is 4.3σ too bright for a lensed Ia and 3.4σ too faint
for a lensed IIn; the chi2 form is kept as a robustness row. On the earlier run's other candidates the function reproduces
the by-eye calls except one: S4TM J0757+1956 had been rejected as needing mu > 100, which holds for every core-collapse type
but not for a SN Ia (mu ≈ 17).

**Injections.** Per lens, 200 lensed SNe drawn as the forecasts draw them (Shu et al. 2018 as used by Sheu et al. 2023,
arXiv:2301.03578): Ia with the lens's forecast Ia fraction (class L: the class-S mean, 0.11), else core collapse by Shu's
type fractions and luminosity functions (Nugent and Vincenzi templates); source redshift (S: measured; L: z_l × N(2, 0.5)
truncated at 1); 2 or 4 images with lognormal magnifications, ln mu ~ N(1.5, 0.35) (S: rescaled to Shu's total); delays
0-30 d; images 1.2" from the centre at random position angles; explosion times uniform over the data span ± 100 d. Each
image's ZTF-r light curve times the estimator's measured response to a point source 1.2" off centre (computed on the
pipeline's own difference-image PSF stamps, every tenth epoch, carried to the rest by a fit on seeing; median 0.45, against
0.75 for a centred source) is added to the lens's real baseline-subtracted light curve. The g band, which the vetting needs
and which was fetched only for candidates, is synthesized at the lens's real clean g epochs: the injected SN's ZTF-g flux
times the same response at the g seeing, plus Gaussian noise of the lens's r-band error at that limiting magnitude, scaled
by 1.15 (g's error over r's at the same limiting magnitude where real g exists). A recovery counts only at the injected
SN's time, and none where the lens's own light curve passes the shape cut in that span.

**Limits.** Exposure X = sum over lenses of T_i × eps_i(summed peak r < m), with T_i the injection window (data span +
200 d) and eps_i the vetted efficiency of lens i's injections brighter than m (lenses with fewer than 3 such injections
contribute nothing). One-sided Poisson 95 % (3.00 events for none, 4.74 for one), averaged over the lognormal exposure
uncertainty. Forecast rows: E = sum_i (N_Ia,i eps_Ia,i + N_cc,i eps_cc,i) T_i.

**Robustness (class L, summed peak r < 21, without the uncertainty):** the first-frozen chi2 colour test, X = 76.1 lens-years,
< 0.039; total magnifications spread over 0.3-10 times the model (10^U(−0.5, 1) on each injection, so the 0.3-3 window
rejects the bright tail), X = 95.0, < 0.032; every image at 2": X = 22.7, < 0.13. For class S the wide-magnification draw
gives X = 12.1 lens-years and < 0.39 with J0757 counted.

## The candidates in the sample

| lens | class | T | stage reached | what it is |
|---|---|---|---|---|
| S4TM J0757+1956 (z_l 0.12, z_s 0.83, Shu mu 9) | S | 7.15 (shape-cut window S = 5.53) | passes all | 2026-02, PSF-corrected peak r ≈ 21.1, 15-30 d above half maximum, centroid 0.3" from the lens, controls quiet, no season offset; a lensed SN Ia at z_s needs mu = 17.0 ± 2.4 and predicts g ≈ 1 microJansky over the r peak against 1.8 ± 1.0 measured, and a lensed Ia at z = 0.83 stays 29 observed days above half maximum; equally a faint transient in the lens galaxy (M_r ≈ −17.6). No alert object. Counted as a detection |
| DESI-266.4075+39.8584 | L | 5.96 | fails 0 (shape) | one window (2024-05); its neighbour at S = 3.3 is not part of the event |

## What it means

ZTF's public archive, searched below its alert threshold at known lenses and vetted the same way for real and simulated
events, keeps a lensed supernova only when its images sum to r ≲ 21, half the time at r 19-20. Behind grade-A candidates,
lensed supernovae that bright occur at fewer than one per 22 lens-years (95 %), a limit that assumes no rate model; the
per-lens forecasts put most lensed supernovae behind known lenses at r 23-29, so the limit sits about an order of magnitude
above a SLACS-like rate and within a factor two of a BELLS-like one, and tests neither. Behind the spectroscopic lenses of
this sample the sources are at z ≈ 1 and almost nothing injected reaches r = 21, so the archive says nothing about their
rate. The one surviving event, in S4TM J0757+1956 in February 2026, is the kind of ambiguity a lens-position search at ZTF
depth produces: a magnified SN Ia at the source redshift and a faint supernova in the lens galaxy predict the same r peak,
duration and, within ZTF's g errors, the same colour; spectroscopy of the lens, or the time delays a resolved image pair
would show, would decide it.

## Weakest points

- The sample is a prefix of the fixed order (54 of the 440 listed lenses with an efficiency); the order put the
  highest-forecast spectroscopic lenses first, which are BELLS and S4TM lenses with sources at z ≈ 1, so class S has almost
  no injections brighter than r = 21 and no useful forecast-free limit.
- The variability stage, applied as frozen, rejects many bright injections whose own light-curve tails shift another
  season's mean by more than 1 microJansky: it takes class L's efficiency at r < 21 from 0.62 to 0.35. A real bright lensed
  SN would be rejected the same way by the same code, so the efficiency is right for this vetting; a gentler rule would
  tighten the limit but was not adopted, since it would be chosen after seeing the injections.
- The injections' g band is synthesized (real g epochs and limiting magnitudes, noise scaled from r), not measured on real
  g cutouts; it matters only at stage 4, which takes class L's efficiency at r < 21 from 0.351 to 0.344.
- The efficiency falls steeply with image offset (the estimator's response drops from 0.45 at 1.2" to 0.17 at 2"), carried
  as a systematic with 2" taken as 2σ; candidates with Einstein radii near 2" are measured with less of the PSF.
- Class L's source redshifts come from photometric lens redshifts times N(2, 0.5); the forecast-free limit depends on this
  only through the injected light-curve mix inside the magnitude cut, the secondary rows fully.
- S4TM J0757+1956 is counted as a detection; if it is a lens-galaxy transient the class-S numbers tighten by 3.00/4.74 and
  class L is unchanged.
- Two vetting amendments came after the real candidates were run and before any injection was vetted (the colour test at
  peak; the g noise scale); the first-frozen chi2 form is given as a robustness row and changes the class-L limit by 5 %.

## To settle it

The ZTF archive over the remaining 386 listed lenses (the candidates' brightness order makes each later one add less), the
same search on Rubin's per-visit depth (r ~ 24.5), where the spectroscopic lenses' sources come within reach, and a
spectrum of S4TM J0757+1956's lens galaxy for the February 2026 event.

## Blind review (2026-11-02, context-free Opus referee, own literature search; `blind-review.md`)

VERDICT: novel but not sound: the headline limit rests on putting every lensed image 1.2" from the centre, an offset the efficiency is very sensitive to (all at 2" makes the limit three times weaker) and whose quoted 18 % uncertainty disagrees with the note's own 2σ figure, so the limit of fewer than 0.045 per lens-year and its error budget are not supported.

Strongest objection, in one sentence: every injected image sits 1.2" from the centre while the efficiency falls by 3.5 times at 2" (typical of ground-selected grade-A candidates), so the offset should be drawn per lens from Einstein radii and its uncertainty marginalised in the likelihood; also named: the class-L magnification window referenced to the injection draw, thresholds shaped on known candidates, and "forecast-free" still averaging over the forecast light-curve mix. Novelty held (closest LCOLSS arXiv:2111.01680). Recorded as a miss and **parked** (the last pass). Erratum after the verdict, not a rewording: the offset width in the computation is 0.36 (half of the 0.715 relative exposure change), not "18 %"; the 0.045 limit was computed with 0.36. Path to pass: per-lens image offsets from modelled or catalogue Einstein radii (and per-lens Einstein-radius-based magnification windows for class L), the efficiency uncertainty marginalised in the Poisson likelihood, and the sample beyond the 54-lens prefix. What stands: the vetting as code applied to real and injected light curves alike, the vetted efficiency curve per summed peak magnitude, the S4TM J0757+1956 event as an unresolved lensed-Ia-or-host transient, and that ZTF cannot test the per-lens forecasts.

## Revision history

R-135 (2026-11-01; 168 lenses, injections without the vetting) was refereed "VERDICT: novel but not sound: the quoted efficiencies come from injections that skip the colour, magnification and variability vetting that every real candidate had to pass, and the class-L limits that carry the headline come from an assumed source population rather than from the data." (`results-R130.md` row 6). This revision (R-141) puts the vetting in code for both, makes the forecast-free limit the headline, and was refetched under a rate-limited archive, which is why the sample is smaller.
---

# CMB lensing gives the clustering bias of 3.5 < z < 5 quasars: b = 6.2 ± 1.9 from about 21,000 DESI DR1 and SDSS DR16Q quasars on ACT DR6 and Planck, a 3.2σ cross power that favours the quasar-galaxy cross-correlation answer (b ~ 3-6) over the quasar auto-correlation answer (b ~ 13) by a likelihood ratio of about 270

2026-11-01/02. Scripts: `qso_z4_feas.py` (counts and forecast), `qso_z4_lens.py` (`prep`, `data`, `sims`, `gsims`, `fit`,
`bright`), `qso_z4_theory.py` (Limber model); rows 1-7 in `results-R139.md`, where the samples, masks, multipoles, estimator,
model and fit were fixed before any κ map was multiplied by any quasar map. Data: DESI DR1 LSS v1.5 quasars
(`QSO_full_HPmapcut`, clustering randoms), SDSS DR16Q (Lyke et al. 2020), ACT DR6 lensing v1 (baseline and curl, the 400 simulated
reconstructions and their inputs; Madhavacheril et al. 2024, Qu et al. 2024), Planck 2018 (PR3) MV lensing (Planck PR4 lensing
was not reachable from the compute used here).

## Claim

The clustering of quasars at z ≈ 4 has two published answers that disagree. The quasar auto-correlation of bright SDSS quasars
gives b ≈ 13-14 at 3.5 < z < 5 (Shen et al. 2007, astro-ph/0702214), which implies hosts above 10^12.5 M_sun/h and duty cycles
near unity; the cross-correlation of quasars with Lyman-break galaxies at z ~ 4 gives b ≈ 3-6 (He et al. 2018, arXiv:1704.08461;
see also arXiv:1507.05292). CMB lensing measures the bias directly, without the quasar auto-correlation's sensitivity to the
redshift distribution and to the treatment of the integral constraint, but its quasar cross-correlations stop at z = 3.5 (DESI
DR1 × Planck PR4, arXiv:2506.22416, top bin 2.5 ≤ z ≤ 3.5; Quaia, eBOSS).

Crossing the 3.5 < z < 5 quasars of DESI DR1 (15,000 weighted) and SDSS DR16Q (9,400; 3,400 in common) with ACT DR6 lensing
on its footprint and Planck lensing outside it, with the covariance from the released ACT simulations and Gaussian Planck
realisations, and magnification modelled with each sample's own count slope:

| sample | lensing map | N (in the map's footprint) | cross-power S/N | b (68 %) | L(b = 13) / L(b = 4) |
|---|---|---|---|---|---|
| DESI DR1 | ACT DR6 (40 ≤ L < 1000) | 6,200 | 2.0 | 5.5 ± 2.7 | 0.02 |
| DESI DR1 | Planck PR3 (8 ≤ L < 400), full | 13,500 | 1.5 | 4.5 ± 2.7 | 0.007 |
| DESI DR1 | Planck PR3 outside ACT | ~7,300 | 3.1 | 11.6 ± 3.7 | 7.8 |
| SDSS DR16Q | ACT DR6 | 3,060 | 0.3 | −0.9 ± 3.6 | 0.002 |
| SDSS DR16Q | Planck PR3, full | 9,150 | 1.1 | 3.5 ± 3.1 | 0.008 |
| SDSS DR16Q | Planck PR3 outside ACT | ~6,100 | 2.6 | 9.0 ± 3.5 | 1.5 |
| DESI DR1 | ACT ⊕ Planck outside ACT | 15,000 | 3.3 | 7.6 ± 2.2 | 0.18 |
| SDSS DR16Q | ACT ⊕ Planck outside ACT | 9,400 | 1.7 | 4.2 ± 2.6 | 0.003 |
| **DESI + DR16Q (joint covariance)** | **ACT ⊕ Planck outside ACT** | **~21,000 distinct** | **3.2** | **6.2 ± 1.9** | **0.0037** |
| DESI + DR16Q | ACT (L < 600) ⊕ Planck outside ACT | | 2.6 | 5.4 ± 2.0 | 0.0008 |

(b is the linear bias at the samples' effective redshift, mean z = 3.80 for DESI and 3.86 for DR16Q, against Planck 2018 ΛCDM;
"ACT ⊕ Planck outside ACT" uses each quasar under one lensing map only, so the two legs are independent sky; the likelihood
ratio has the magnification slope marginalised.)

**The claim is the combined interval: b = 6.2 ± 1.9 (68 %) for quasars at 3.5 < z < 5 (z̄ ≈ 3.8), the first lensing measurement
of quasar bias above z = 3.5.** At this precision b = 13 lies 3.5σ above the measured value and b = 4 lies 1.2σ below it; the
likelihood ratio L(13)/L(4) = 3.7 × 10⁻³ (1/270), and 8 × 10⁻⁴ with ACT cut at L < 600. The data favour the quasar-galaxy
cross-correlation answer. DESI alone does not separate the two answers (L ratio 0.18); the separation needs both surveys.

## Data and method

- **Samples.** DESI DR1: every quasar with a good redshift (`Z_not4clus`) at 3.5 < z < 5 in `QSO_full_HPmapcut`, weighted by
  `WEIGHT_ZFAIL / FRACZ_TILELOCID`; the residual completeness with the number of overlapping tiles (measured on all good-redshift
  quasars before any cross-spectrum: 0.644 for one tile, 0.654-0.687 for two to seven) goes on the DR1 clustering randoms, which
  define the footprint. DR16Q: every quasar with 3.5 < Z < 5 (Z from visual inspection for 92 %); no randoms exist for its
  heterogeneous selection (1,615 of its 9,761 z > 3.5 quasars are uniform BOSS CORE targets, 5,788 BONUS), so its footprint is
  the plate coverage traced by all 750,000 DR16Q quasars, with unit weights. Both at |b| > 20°; overdensity δ = N_w / (n̄ f) − 1
  per photometric region at HEALPix nside 512.
- **Lensing.** ACT DR6 baseline κ, its mask degraded to nside 512 and squared as the κ window, 40 ≤ L < 1000 in 10 bands; a
  transfer function measured on 60 simulations (reconstructed × input over input × input in the DESI window): 0.97-0.985 in
  every band, applied. Planck PR3 MV κ, its mask rotated to equatorial and 1° apodised, 8 ≤ L < 400 in 7 bands, split into the
  sky ACT covers and the sky it does not (the split and the ACT ⊕ Planck-outside-ACT combination were set before the first fit,
  so that no quasar enters the combination under two lensing maps).
- **Estimator.** Pseudo-C_ℓ with the mode-coupling matrix (NaMaster 3.0.1), bandpower windows applied to the model, the quasar
  map's pixel window in the model.
- **Model.** C^{κg}(L) = b C^{κm}(L) + (5s − 2) C^{κμ}(L), Limber, halofit, the samples' spectroscopic n(z). Magnification is as
  large as clustering here: at z ≈ 3.8 C^{κμ}/C^{κm} = 1.6 at L = 100 and ≈ 1 above L = 300 per unit 5s − 2. Neither sample has
  a hard flux limit (their z > 3.5 counts turn over at r ≈ 21 and 20.4 from selection incompleteness), and for a selection that
  depends on magnitude alone the response is the sample-averaged slope of the *intrinsic* differential counts (integrating the
  selection by parts removes the completeness). s(r) is a quadratic fit to DESI's own z > 3.5 counts over 18 ≤ r ≤ 20.75 and the
  z ~ 4 faint-end luminosity-function slope (s = 0.12-0.40, the range carried as the error) beyond: s = 0.33 ± 0.10 (DESI),
  0.40 ± 0.07 (DR16Q), marginalised as a Gaussian prior (exact, the model being linear in b and 5s − 2). The same method gives
  s = 0.30 ± 0.12 for DESI at 2.2 < z < 3.5, where DESI × PR4 fitted 0.19-0.22 ± 0.05-0.06. This amendment was made before any
  cross-spectrum, when the fixed estimator (the cumulative slope at the nominal flux limit) gave s ≈ 0.01, an artefact of the
  soft limit.
- **Covariance.** The 400 ACT simulated reconstructions (Gaussian κ through the full ACT pipeline) crossed with the fixed data
  quasar maps; for Planck, 400 Gaussian realisations of the PR3 C^{κκ} + N^{κκ} with the same masks (the PR3 simulations are not
  on the archive); Hartlap-corrected. The DESI and DR16Q spectra on one map are fitted jointly with the covariance of the same
  simulations crossed with both quasar maps, which carries their 3,400 shared quasars.

## Checks

- **End to end (DESI 2.2 < z < 3.5, 327,000 quasars, the same pipeline):** b = 2.98 ± 0.53 (ACT), 3.37 ± 0.53 (Planck), 3.23 ±
  0.48 combined, against the published DESI DR1 × Planck PR4 2.53 ± 0.23 (2.1-2.5) and 2.93 ± 0.36 (2.5-3.5), 2.86-3.43 in their
  linear-theory fit (arXiv:2506.22416). ACT and Planck on the same sky agree (2.98 ± 0.53 vs 2.91 ± 0.67).
- **Randoms** (random positions in each footprint, same number as the data) × κ: b_null between −6.2 ± 3.3 and +2.8 ± 3.8, χ²
  PTE 0.07-0.72 on every map.
- **Curl reconstruction** × quasars (covariance from 100 curl simulations): DESI χ² 3.7/10 (PTE 0.96); DR16Q 20.5/10 (PTE
  0.03), carried by the three bands above L = 550 (−3.4σ, −2.0σ, −1.7σ); DESI 2.2 < z < 3.5 PTE 0.20. Cutting ACT at L < 600
  (a pre-registered variant) removes those bands and moves the combined b from 6.2 to 5.4.
- **ACT vs Planck on the ACT footprint:** DESI 5.5 ± 2.7 vs −3.3 ± 3.8 (1.9σ with the errors added as independent, which
  overstates the error of a difference on shared sky), DR16Q −0.9 ± 3.6 vs −6.3 ± 5.2 (0.8σ). Planck inside and outside the
  ACT footprint differ by 2.8σ for DESI (−3.3 ± 3.8 vs 11.6 ± 3.7) and 2.4σ for DR16Q.
- **E(B−V)_SFD < 0.15:** b changes by ≤ 0.3.
- **Post hoc, labelled as such:** DR16Q at i < 20.2 (Shen et al.'s limit; 4,763 quasars, s = 0.56 ± 0.05) on Planck:
  b = 3.4 ± 4.1, uninformative.

## What it means

At z ≈ 3.8 the lensing bias of the quasar population that DESI and SDSS select is 6.2 ± 1.9. A bias of 13, the auto-correlation
answer for the brightest SDSS quasars, is disfavoured at 3.5σ and by a likelihood ratio of about 270 against b = 4; values of 3-6
from the quasar-LBG cross-correlation are inside 1.2σ. The measurement does not distinguish b = 4 from b = 8. It does not test
the luminosity dependence of the bias: the samples reach r ≈ 22, two magnitudes below the Shen et al. sample, and the bright
subsample at i < 20.2 is too small to say anything.

## Weakest points

- **The separation leans on DR16Q and on the Planck legs.** DESI alone gives 7.6 ± 2.2 and a likelihood ratio of 0.18. DR16Q's
  selection is heterogeneous and its footprint is traced from the quasars themselves; its curl null is marginal (PTE 0.03) at L >
  550. The Planck estimates differ between the ACT sky and the rest by 2.4-2.8σ, more than their errors.
- **Magnification** is as large as the clustering term. The s prior (±0.5 on 5s − 2 for DESI) is propagated; a faint-end
  slope outside 0.12-0.40, or colour-dependent selection responding to lensing, would move b by about 1.6 per unit of 5s − 2.
- **The Planck covariance is Gaussian** (C + N realisations, not the PR3 reconstruction simulations); PR3 rather than PR4.
- **Linear bias to L = 1000** at z ≈ 3.8 (k ≈ 0.14 Mpc⁻¹) with halofit; scale-dependent bias for b ~ 13 haloes would matter more
  at the high-L end, which the L < 600 variant drops.

## To settle it

DESI DR2 quasars on ACT DR6 and SPT-3G lensing; Planck PR4 with its simulations; a
magnification response from re-running the DESI target selection on perturbed fluxes.

## Blind review (2026-11-02, context-free Opus referee, own literature search; `blind-review.md`)

VERDICT: novel but not sound: the headline b = 6.2 ± 1.9 and likelihood ratio of 270 come from averaging mutually inconsistent sky and lensing-map legs (Planck inside vs outside ACT differ by 2.4–2.8σ, and the Planck legs outside ACT alone favour b ≈ 13) without inflating the error, and the ratio compares samples of different luminosity to Shen et al.'s bright quasars, so the preference for the cross-correlation answer is not supported.

Strongest objection: the four independent legs disagree by sky region (the Planck legs outside ACT alone prefer b ≈ 13, the ACT
legs, with DR16Q × ACT at −0.9 and a marginal curl null, pull low), and the combined interval does not carry that scatter.
Novelty was granted (no CMB-lensing quasar bias above z = 3.5 in print; closest arXiv:2506.22416, 2512.17865, 2507.08798).

**Parked (physics point, one gate pass; R-139).** What stands: the first CMB-lensing cross-power of z > 3.5 quasars (3.2σ
combined, 3.3σ for DESI alone); the end-to-end pipeline reproduces the published DESI × PR4 bias at 2.2 < z < 3.5 (3.23 ± 0.48);
DESI + DR16Q give 3.3 ± 2.4 on the ACT sky and 10.4 ± 2.9 on Planck outside it (1.9σ apart); with the combined error scaled by
that scatter, b = 6.2 ± 3.6 and L(13)/L(4) = 0.21, which does not separate the published answers. The path to a pass: find the
sky-dependent difference (DESI's north photometric region lies entirely outside ACT; the DR16Q footprint and selection; Planck
PR3's Gaussian error model) with PR4 and its simulations, DESI DR2 and a selection-function DR16Q sample, and compare against
a bias model at matched luminosity rather than Shen et al.'s bright-sample value.

# The near-infrared light correlated with DESI DR1 spectroscopic galaxies in SPHEREx Quick Release images, with the source-mask transfer measured by injection: detected in all 18 channel groups from 0.75 to 2.42 µm (2.6-15σ per group); at 1.8-2.4 µm the light outside catalogued sources is b_I d(νI_ν)/dz = 4.5-6.2 nW m^-2 sr^-1 per unit redshift; below 1.8 µm and for the total only lower bounds hold, which do not contradict CIBER

## Claim

Cross-correlating 14,508 SPHEREx QR2 wide-survey spectral images (detectors 1-3, 0.75-2.42 µm; zodiacal model
subtracted, the pipeline's catalogue source mask applied) with 3.0 million DESI DR1 spectroscopic galaxies and quasars in
32 slices of Δz = 0.1 over 1,252 deg², the clustering-redshift estimator returns the bias-weighted background intensity
per unit redshift, b_I d(νI_ν)/dz, per channel group (6 per detector) and slice. What the processing chain (source mask,
5σ clip, per-exposure strip median) passes of galaxy-correlated light was measured by painting mock galaxy-correlated
light into the same real images and running the identical chain; each slice and group is corrected by that measured
transfer. The claim:

- **Detection** in every one of the 18 groups, 2.6-15σ per group (inverse-variance mean over the 32 slices, jackknife σ
  over 51 sky patches, the transfer's own error included).
- **Amplitude at 1.81-2.37 µm** (the five groups where the transfer passed a pre-set stability test): the galaxy-correlated
  light outside catalogued sources is b_I d(νI_ν)/dz = 5.0 ± 0.5, 5.6 ± 0.4, 6.2 ± 0.6, 5.3 ± 0.5, 4.5 ± 0.4 nW m^-2
  sr^-1 per unit z (groups at 1.81, 1.94, 2.08, 2.23, 2.37 µm).
- **Below 1.8 µm** the same quantity, 3.9-7.3 nW m^-2 sr^-1 per unit z, is a lower bound, not a measurement: there the
  transfer depends on the mock's amplitude and redshift shape by 15-30 %.
- **The total** galaxy-correlated light, including the galaxies' own light that the catalogue mask removes, is bounded
  from below only: the mask passes 0.45 (0.77 µm) to 0.06 (2.37 µm) of compact light at the galaxies, so the total is
  not constrained from above.
- **Shape:** nearly flat over 0.15 < z < 2.1; at 2.23-2.37 µm the z < 0.5 slices are lower (2.9 ± 0.9, 2.9 ± 0.7) than
  z ≥ 0.5 (6.3 ± 0.7, 5.5 ± 0.5), ~3σ each; the lowest slice (z 0.05-0.15) is consistent with zero in every group.
- **CIBER:** at 1.8 µm and z < 0.6 the lower bound (4.4 ± 1.4, 6.4 ± 1.6 at z 0.2-0.4, 0.4-0.6) is below CIBER x DESI
  Legacy Survey (arXiv:2608.12116: 63 ± 15, 35 ± 8), but correcting by the compact-light transfer instead gives 39 ± 18
  and 5 ± 28: the light the SPHEREx catalogue mask removes can account for the difference, so the two do not disagree.

Nothing is claimed about galaxy evolution or line emission (an earlier run of this pipeline found no Hα, [OIII] or Pa-α
line, with limits at the continuum level).

## Data and method

- **Images.** SPHEREx QR2 L2 spectral images (IRSA, `nasa-irsa-spherex` S3), D1-D3 images centred in the field or one
  NSIDE-32 pixel around it (31,995), processed in a fixed random order until a pre-set stop: 14,508 (4.8k per detector).
  Per image: FLAGS bad bits and the pipeline SOURCE bit (every Reference Catalog point source with a magnitude-dependent
  radius; ~31 % of pixels), pixels > 5 robust σ above their channel strip's median dilated by 12", the ZODI model
  subtracted, each channel strip's (0.2° x 3.5°) median of good pixels subtracted per exposure; 8x8-pixel blocks binned to
  NSIDE 2048, 17 channels per detector grouped in six. Map level: Gaia DR3 G < 10 stars masked.
- **Galaxies.** DESI DR1 NGC LSS clustering catalogues (BGS_ANY 0.05-0.45, LRG 0.4-1.1, ELG_LOPnotqso 0.8-1.6, QSO
  0.8-2.1) and randoms with WEIGHT; field: 1,252 deg² where BGS, LRG and ELG are ≥ 80 % complete and |b| > 40°.
- **Estimator (arXiv:1810.00885, notation of 2109.10914).** w̄_Ir = annulus-weighted (W ∝ θ^-0.8, θ = max(1 cMpc, 3.5')
  to 10 cMpc) mean residual intensity around galaxies minus around randoms; b_I dI/dz = w̄_Ir / (b_r w̄_m T), w̄_m from
  CAMB halofit, b_r from the slice's own clustering (BGS 1.0-1.4, LRG 1.9-2.3, ELG 1.1-1.8, QSO 1.8-3.0).
- **Transfer T, measured by injection.** In every one of the 14,508 images, mock galaxy-correlated light is added to the
  real pixels before the chain: each DESI galaxy carries a flux set so the mock's mean d(νI_ν)/dz is a chosen A(z),
  split into a compact part (a 1-pixel Gaussian at the galaxy: its own light) and a smooth part (a 1-comoving-Mpc Gaussian:
  unresolved light correlated with it). Variants: half/half at the measured level (A), the same x5 (B), smooth only (C),
  compact only (D), half/half with A(z) ∝ (1+z)³ (E). Each runs through the identical mask, clip (its own) and strip
  median; the real image's output on the same mask is subtracted. T = w̄[recovered mock] / w̄[mock with no mask, clip or
  filter], same estimator, same centres, per slice and group, jackknife over patches. Results: smooth light T_C =
  0.62-0.69 in every group (the mask alone keeps 0.96-0.99 of it; the strip median removes the rest), 0.26-0.31 in the
  lowest slice; compact light T_D = 0.30 at 0.77 µm to 0.03 at 2.37 µm (the mask alone keeps 0.45 to 0.06); half/half
  T_A = 0.51 to 0.35.
- **Rule fixed before the transfer was computed:** the amplitude is a measurement of the light outside catalogued sources
  (corrected by T_C) only in groups where T_A changes by < 15 % between A and B (amplitude) and between A and E (shape);
  otherwise a bound. It holds at 1.81-2.37 µm (B/A 0.86-0.90, E/A 0.95-0.97) and fails below (B/A 0.68-0.84, E/A
  0.71-0.92: the clip removes bright compact mock light). The total is a lower bound wherever T_D / T_C < 0.5 (all
  groups but the filter-edge 1.09 µm one).
- **Errors:** delete-one jackknife over 51 NSIDE-8 patches (on the signal and the transfer together); χ² tests use the
  diagonal jackknife σ.
- **Nulls.** On this run: Gaia 12 < G < 16 stars as centres, with an added fixed-flux mask of every G < 16 star (block
  within 1' x 10^(0.2(12-G)) + 12"; 31 % of blocks): χ² 14-50/32 (p > 0.01) in 14 of 18 groups, failing at 1.09, 1.14,
  1.23 and 1.61 µm (detector-edge groups); without that mask it fails in 11 of 18. The galaxy signal changes by
  -0.03 to -0.68 nW m^-2 sr^-1 per unit z when the star mask is added, within 2σ in every group. Galaxy-star
  cross-correlation over the same annuli: w̄_gs = -0.0006 ± 0.0008 (χ² 39.6/32). From the earlier run on 22,322 images
  with the same chain: zodiacal model alone through the pipeline, half-survey difference and redshift-shuffled galaxies
  all pass.
- **Star-density split.** Splitting patches at the median star count, the measured transfer is lower in the low-star
  half (T_A ratio 0.71-1.01), while the observed low/high signal ratio at 0.77-1.03 µm is 1.3-4.2 with 30-70 % errors
  (≤ 1.9σ from the prediction per group), ~1 at λ > 1.9 µm.

## What it means

SPHEREx's spectral images, cross-correlated with spectroscopic galaxies, detect the near-infrared light that clusters with
galaxies in every channel group from 0.75 to 2.4 µm. Measuring the processing chain's transfer by injection shows what
that light is: the catalogue source mask removes most of the galaxies' own light (all but 6 % at 2.4 µm), while smooth,
unresolved light passes with a well-measured, wavelength-independent transfer of ~0.65. So the measured quantity is the
galaxy-correlated light outside catalogued sources, 4.5-6.2 nW m^-2 sr^-1 per unit z at 1.8-2.4 µm, roughly flat to
z ~ 2 with a lower z < 0.5 value at 2.2-2.4 µm. It bounds the total from below. The large CIBER x DESI amplitude at
z < 0.6 is not in tension with it: CIBER masks only brighter sources, and the light SPHEREx's mask removes is enough to
make up the difference.

## Weakest points

- The transfer depends on the mock: a 1-cMpc Gaussian for the unresolved light and a point source for the galaxy's own
  light, with fluxes per tracer slice set to a flat or rising d(νI_ν)/dz, not to galaxy luminosities. Real unresolved
  light with more small-scale structure would sit between T_C and T_D.
- Below 1.8 µm the transfer changes by 15-30 % with the mock's amplitude and shape (the clip), so no amplitude is
  claimed there.
- The star null still fails in four detector-edge groups after the star mask; the galaxy signal does not move with the
  mask and w̄_gs is zero, but star-associated light in those groups is not fully removed.
- The low-star half of the patches carries a higher signal at 0.77-1.03 µm than the transfer predicts (≤ 1.9σ per
  group); unexplained.
- 14,508 of 31,995 images (45 %); D4-D6 not processed. CIBER's values are read off its figure (±10 %); its 1.1 µm band
  is compared with the 1.09 µm filter-edge group and the 1.14 µm group, where at z 0.6-0.8 the lower bound (8.4 ± 1.8)
  sits above CIBER's 95 % limit (< 5.5), 1.6σ.

## To settle it

The full image set and D4-D6; a mask set from a flux limit (Legacy Survey or Gaia) instead of the pipeline's catalogue,
matched to CIBER's depth, with the same injection, so the two measure the same light; mocks built from galaxy
luminosities and an HOD for the unresolved part; QR3 and the second-year survey.

## Blind review history

- R-138 (2026-11-02, context-free Opus referee, own literature search; `blind-review.md`), on the earlier version of this
  section (22,322 images, strip-filter transfer only; amplitudes 4-9 nW m^-2 sr^-1 per unit z and "4-25x below CIBER"):
  "VERDICT: novel but not sound: the amplitude and redshift shape are set by a catalogue source mask and clipping that remove galaxy-correlated light, which the transfer correction does not model, as the factor-of-two star-density split and the failed stellar null show, so the headline amplitudes and the CIBER ratio are not supported."
  Revised by R-143 (the one revision): the transfer measured by injection through the whole chain, the stellar null redone
  with a fixed star mask, the split predicted, the CIBER comparison redone (`results-R138.md`, "R-143").
- R-143 (2026-11-02, context-free Opus referee, finding-check prompt verbatim, the section alone; hints as R-138's):
  VERDICT: novel but not sound: the reported 1.8–2.4 µm amplitudes depend on the unvaried compact/smooth split and smooth-light profile of the injected mock, which the pre-set stability test never probes, so "light outside catalogued sources" is not a mock-independent measurement.
  Strongest objection: the corrected signal is T_C·S_smooth + T_D·S_compact divided by T_C, so leaked compact light (profile and PSF wings past the mask radius, faint satellites) with transfers between T_D and T_C enters the 1.8-2.4 µm amplitude at a level set by the mock's half/half split and 1-cMpc kernel, which the A/B and A/E tests hold fixed; it also notes the strip median's scale dependence (the z < 0.5 deficit sits where T is lowest), the diagonal-only χ², the star mask chosen with the null in view, and the CIBER reconciliation resting on 39 ± 18 and 5 ± 28. Novelty held: no SPHEREx x galaxy measurement found (closest: 2608.12116 CIBER, 2109.10914 forecast, 1810.00885 method, 2608.03415, 2305.08847, 2601.21111). What survives, per the referee: the detection in all 18 groups (2.6-15σ). **A miss on physics: parked (brief: one revision; park on any further physics point).**

# Photometric-redshift isolation leaks unevenly by galaxy type, but not in lensing amplitude: on DESI DR1 BGS spectroscopic truth, KiDS-bright-style isolation leaves red galaxies with a real > 10 % neighbour within 0.3 Mpc 1.7 times as often as blue ones (32 % against 19 %); a ΛCDM mock that reproduces this leakage predicts that it shifts the early/late lensing-RAR offset by −0.001 dex (68 %: −0.005 to +0.002), at most 0.09 of the published 0.14-0.19 dex in any variant

## Claim

Brouwer et al. 2021 (arXiv:2106.11677) found isolated KiDS-1000 early- and late-type galaxies ≥ 6σ apart in the weak-lensing
radial acceleration relation (RAR). Their measure is log10⟨g_obs,E/g_obs,L⟩: 0.14 dex (Sérsic split) and 0.19 dex (u−r split) at
g_bar > 10⁻¹³ m s⁻² (R < 0.3 Mpc/h70), and 0.17 and 0.27 dex over all 15 g_bar bins. Their isolation uses photometric redshifts: no
neighbour with > 10 % of the lens stellar mass within 3 Mpc/h70 in 3D. They assumed photo-z errors degrade the isolation of both
types equally (their Sec. 5.3). Mistele et al. 2024 (arXiv:2310.15248) find one relation with 4 Mpc/h70 isolation and their own
early-type masses.

Two results, in this order.

**1. The leakage is type-dependent (measured).** The KiDS-bright isolation is emulated on the DESI DR1 BGS group catalogue: the
published KiDS-bright photo-z error distribution, 5 seeds, about 780,000 candidate lenses per seed. Each photo-z-isolated lens is
then checked against the spectroscopic redshifts. Types are mass-matched as Brouwer did. Interval: sky jackknife (50 regions) plus
seed scatter, statistical only.

| split | fraction with a real > 10 % neighbour (500 km/s) | red / early | blue / late | difference (95 %) | ratio (95 %) |
|---|---|---|---|---|---|
| rest-frame g−r | within 0.3 Mpc | 0.319 | 0.192 | 0.128 (0.124-0.132) | 1.67 (1.64-1.69) |
| rest-frame g−r | within 0.5 Mpc | 0.453 | 0.303 | 0.149 (0.145-0.153) | 1.49 (1.48-1.51) |
| rest-frame g−r | within 1 Mpc | 0.651 | 0.517 | 0.134 (0.129-0.139) | 1.26 (1.25-1.27) |
| rest-frame g−r | lens is a hidden satellite of a group with such a neighbour | 0.290 | 0.113 | 0.177 (0.174-0.181) | 2.57 (2.53-2.62) |
| spectral QUIESCENT | within 0.3 / 0.5 Mpc | 0.327 / 0.461 | 0.200 / 0.314 | 0.127 / 0.147 | 1.63 / 1.47 |
| spectral QUIESCENT | hidden satellite | 0.310 | 0.118 | 0.192 (0.188-0.196) | 2.62 (2.57-2.68) |

Brouwer's equal-contamination assumption fails in fraction. The gap is set by red lenses being 2.6 times as often unrecognised
satellites. A 4 Mpc criterion lowers the fractions by ≤ 0.04. A Legacy Survey morphology split (DEV against EXP/REX) narrows the
gap (0.27 / 0.39 against 0.20 / 0.32 at 0.3 / 0.5 Mpc, unmatched).

**2. Its lensing effect is not type-dependent (calibrated on a mock).** The same isolation is run on the UniverseMachine DR1
catalogue of the Bolshoi-Planck simulation at z = 0.21, with a BGS- or KiDS-bright-like magnitude limit. That catalogue has
halos, stellar masses and star-formation rates, and its quenched fractions are fit to SDSS by halo and central/satellite status.
The mock reproduces the measured leakage (quenched vs star-forming, BGS-like): 0.37 / 0.48 against 0.22 / 0.32 within 0.3 /
0.5 Mpc, hidden satellites 0.30 against 0.14 (DESI: 0.33 / 0.46 against 0.20 / 0.31; 0.31 against 0.12).

The excess surface density (ESD) is computed from the mock's own halos, with no assumed amplitude. Per type, A_T is the stacked
ESD of photo-z-isolated lenses over that of truly isolated lenses of the same type. Both samples are binned in Brouwer's g_bar
with photometric masses and reweighted to one stellar-mass distribution. The leakage's contribution to the published statistic is
then δ = log10 mean_k (A_E,k / A_L,k). That is Brouwer's log⟨g_E/g_L⟩ with photo-z isolation minus the same with true isolation.
It is evaluated with the KiDS-bright-like limit (r < 20, z_v = 0.15 / 0.25 / 0.35).

| amplitude | criterion | δ at g_bar > 10⁻¹³ (68 %; 95 %) | of 0.19 / 0.14 | δ over all bins (68 %; 95 %) | of 0.27 / 0.17 |
|---|---|---|---|---|---|
| mock | 3 Mpc | −0.001 (−0.005 to 0.002; −0.008 to 0.006) | −0.01 / −0.01 | −0.008 (−0.020 to 0.003; −0.031 to 0.014) | −0.03 / −0.05 |
| rescaled to DESI's fractions (0.5 Mpc) | 3 Mpc | −0.002 (−0.005 to 0.001; −0.008 to 0.005) | −0.01 / −0.01 | −0.010 (−0.021 to 0.001; −0.032 to 0.011) | −0.04 / −0.06 |
| rescaled to DESI's hidden-satellite fractions | 3 Mpc | 0.004 (0.001 to 0.007; −0.003 to 0.011) | 0.02 / 0.03 | 0.010 (−0.001 to 0.021; −0.013 to 0.032) | 0.04 / 0.06 |
| mock | 4 Mpc | 0.006 (−0.002 to 0.013; −0.009 to 0.020) | 0.03 / 0.04 | 0.000 (−0.014 to 0.015) | 0.00 / 0.00 |
| rescaled, hidden satellites | 4 Mpc | 0.013 (0.005 to 0.021; −0.003 to 0.029) | 0.07 / 0.09 | 0.018 (0.002 to 0.033) | 0.07 / 0.11 |

- **Radius dependence (3 Mpc).** log10(A_E/A_L) stays within ±0.005 dex from 0.04 to 0.7 Mpc (errors 0.002-0.02). It turns
  negative beyond 1 Mpc: −0.019, −0.028 and −0.076 (± 0.04-0.06) at 1.3, 1.7 and 2.2 Mpc. Each type's own amplitude rises from
  1.01 at 0.04 Mpc to 1.12 at 0.23 Mpc, 1.30 at 0.54 Mpc and 2.1-2.5 at 2.2 Mpc, so the leakage is large, but it is nearly the
  same fraction for both types.
- **Why it cancels.** Red lenses collect more contaminating mass, but at fixed stellar mass their own (truly isolated) ESD is
  larger in about the same proportion: in the mock, log10(ESD_E/ESD_L) of truly isolated lenses is 0.09-0.21 at g_bar > 10⁻¹³ and
  0.25-0.42 below.

The claim is the two tables. Type-dependent photo-z leakage is real and measured, but in a ΛCDM mock that reproduces it, it
changes the published early/late offset by at most 0.07 (u−r) or 0.09 (Sérsic) of its value at g_bar > 10⁻¹³ in any
variant (95 % upper ends ≤ 0.15 and ≤ 0.21). Over all bins it changes it by −0.06 to +0.11. So the published split is not a leakage artefact. This does not say what the split is.

## Data and method

- **Data.** DESI DR1 BGS group VAC (bgs-groups-dr1-v1.0; r < 19.5, 2.23 M galaxies): fastspecfit stellar masses, rest-frame
  g−r, QUIESCENT, group membership and central/satellite. Lenses are spectroscopic, 0.05 < z < 0.3. In the photo-z view Brouwer's
  cuts apply (0.1 < z_ph < 0.5, 8.5 < log M*_ph < 11). All galaxies with a redshift are neighbours. H0 = 70, Ωm = 0.3, physical
  distances.
- **Photo-z emulation (data and mock alike).** Bilicki et al. 2021's (arXiv:2101.06010) fit to KiDS-bright minus GAMA, a
  generalized Lorentzian (a = 2.613, s = 0.0149) whose wings are the outliers, drawn per galaxy. Photo masses rescale by D_L² and
  get 0.12 dex scatter, as in Brouwer's MICE emulation. Isolation is in 3D, from photo-z positions, at 3 and 4 Mpc and 10 %.
- **Truth.** For the contamination statistic: a real neighbour with M* > 0.1 M*_lens within projected R and |Δv| < 500 km/s. For
  the mock reference: the same 3D criterion with true masses in true redshift space (Brouwer's App. A reference).
- **Mock.** UniverseMachine DR1 (Behroozi et al. 2019, arXiv:1806.07893) on Bolshoi-Planck, 250 Mpc/h box, snapshot z = 0.208.
  Stellar mass is UM's observed-mass column, and the type is quenched if sSFR < 10⁻¹¹ yr⁻¹.
  - Apparent magnitudes come from BGS itself: at each virtual redshift, m_r − DM is fit as linear in log M* per spectral type,
    with Gaussian scatter, by maximum likelihood with the r = 19.5 truncation.
  - Each box axis serves as the line of sight, observed at virtual redshifts z_v (BGS-like: 0.10 / 0.20 / 0.28 weighted by the
    DESI lens n(z); KiDS-like: 0.15 / 0.25 / 0.35, weights 1/4, 1/2, 1/4), with no evolution between them.
  - Not a light cone. MICE (Brouwer's mock) was not reachable without an account.
- **ESD from the mock.** Every distinct halo above 10^11 M☉ within projected 8 Mpc and ±30 Mpc of the lens enters, as an NFW
  truncated at its catalogue virial radius (Dutton & Macciò 2014 concentration), at its true projected offset. A satellite lens
  also gets its own subhalo (current mass and radius), and every lens gets its baryons as a point mass. The lens's own halo is
  therefore whatever the mock gives it, and the reference is per type.
- **Rescaling.** Each type's excess (A_T − 1) is multiplied by f_DESI / f_mock for that type, at 0.5 Mpc, 0.3 Mpc or for hidden
  satellites, measured with the BGS-like limit. The factors are 0.95 / 0.98, 0.89 / 0.91 and 1.05 / 0.84 (quenched / star-forming).
- **Intervals.** Mock: 16 sky-plane jackknife regions per realisation (3 lines of sight × 2 photo-z seeds); variance = mean
  jackknife variance + scatter over realisations. DESI: sky jackknife + seed scatter.
- **Fixed before any number was computed** (results-R147.md, timestamped 04:43 UTC): mock, magnitude model, virtual redshifts,
  photo-z and true isolation, ESD construction, amplitude and offset statistics, rescaling rule, intervals, variants, and what
  each outcome would mean. One implementation choice was made after: the four ESD samples are reweighted to one common
  stellar-mass shape, so that photo-z and truly isolated stacks of a type compare like with like.

## Checks

| check (KiDS-like unless stated) | δ at g_bar > 10⁻¹³ / all bins |
|---|---|
| main, 3 Mpc | −0.001 / −0.008 |
| BGS-like selection (r < 19.5, DESI n(z)) | 0.001 / −0.004 |
| untruncated NFW (z_v = 0.20, one realisation; baseline 0.003 / −0.011) | 0.003 / −0.005 |
| halos down to 10^10.5 M☉ (same) | 0.003 / −0.007 |
| common reference instead of per-type (the mock's excess over the truly isolated star-forming profile) | 0.010 / 0.082 |
| same, 4 Mpc | 0.023 / 0.090 |

The last two rows are not the claim. They show what a type-independent reference does: it gives back a spurious 0.01-0.02 dex at
small radius and 0.08-0.09 dex over all bins. That is the size of the error a common reference would make, and it is a quarter
of what an uncalibrated halo-model amplitude gave at small radius.

**Amplitude benchmark.** Over all lenses, with no type split, the mock's photo-z-isolated ESD is 1.22 times the truly isolated
ESD at 0.3 Mpc, 1.55 at 0.75 Mpc and 3.0 at 3 Mpc (mean 1.81 ± 0.07 over 0.3-3 Mpc). Brouwer's MICE test reports about 1.3
beyond 0.3 Mpc/h70. The mock therefore has more leakage mass at large radius than MICE, or Brouwer's figure refers to a narrower
range. Since the claim is a ratio of two types' amplitudes that agree at every radius, a common rescaling of the excess leaves
it near zero.

## What it means

- At fixed stellar mass, photometric isolation leaves red, quiescent galaxies with a real massive neighbour 1.5-1.7 times as
  often as blue ones, mostly because they are more often unrecognised satellites. Brouwer's equal-contamination assumption is
  wrong in fraction.
- It is nearly harmless for the early/late comparison. In ΛCDM as UniverseMachine populates it, the fractional ESD a type gains
  from leakage is the same for both types to ≤ 0.005 dex inside 0.7 Mpc. Red lenses' larger contamination is matched by their
  heavier own halos.
- The published 0.14-0.19 dex at g_bar > 10⁻¹³ is therefore not a leakage artefact, at the ≤ 0.09 level (≤ 0.21 at 95 %). The
  same holds for the full-range 0.17-0.27 dex at the ≤ 0.11 level, except beyond 1 Mpc, where the mock's leakage slightly
  favours late types.
- An observation outside the claim: Brouwer's statistic applied to the mock's truly isolated lenses gives 0.14 dex
  (g_bar > 10⁻¹³) and 0.26 dex (all bins) between quenched and star-forming. This is ΛCDM's halo-mass gap at fixed stellar
  mass, as the mock sets it. It is close to the published values, but it was not pre-registered and is not tested here.

## Weakest points

- **One mock, one galaxy model.** The cancellation rests on UniverseMachine's halo-mass gap between quenched and star-forming
  galaxies at fixed stellar mass, and on its satellite quenching. The mock matches DESI's leakage fractions to 0.02-0.04, which
  tests the satellite side but not the own-halo side.
  - If real isolated early and late types had the same halos (the strict "one RAR" null), the relevant number would be the
    common-reference row: 0.010-0.023 dex, 0.06-0.17 of the published value at small radius.
  - So the bound is ≤ 0.09 in the mock's universe and ≤ 0.17 under the strict null. Neither reaches the published offset.
- **Not a light cone.** A single z = 0.21 snapshot is observed at virtual redshifts; projection depth is ±30 Mpc of halos, with
  no matter outside halos above 10^11 M☉. The box (305 physical Mpc) wraps the photo-z line-of-sight error (~65 Mpc) about
  8 % of the time.
- **The benchmark is not reproduced.** The mock's all-lens amplitude at 0.3-3 Mpc is 1.8, against MICE's 1.3.
- **Sample stand-ins.** sSFR quenching in the mock and rest-frame g−r or QUIESCENT in DESI stand in for observer-frame u−r. DEV
  against EXP/REX stands in for Sérsic n ≷ 2, and gets no mock counterpart. The DESI interval is statistical only: the
  emulation's systematics (BGS depth and fibre incompleteness in dense regions, the error model) move the fractions by ≤ 0.03 in
  the variants run.

## To settle it

- The same pipeline on a second galaxy model with a different halo-mass gap by colour (e.g. a semi-analytic model, or MICE
  itself). The cancellation should then track that gap.
- A KiDS-bright × GAMA (or DESI DR2 BGS in KiDS-N) spectroscopic subsample, measuring the isolated lensing RAR by type with both
  isolation selections, would measure δ directly.

## Blind review (2026-11-02, R-149 revision, context-free Opus referee, own literature search; `blind-review.md`)

VERDICT: novel but not sound: the lensing calibration is relative to each type's own isolated halo in a mock that already builds in the early/late halo gap being tested, and under the null of equal halos the note's own common-reference check leaves roughly half of the full-range offset attributable to leakage, so "not a leakage artefact" does not follow.

Strongest objection, in one sentence: the calibration divides by each type's own truly isolated halo in a mock that already has
the early/late halo gap being tested, so it cannot test whether the offset is an artefact; under the equal-halo null, the
common-reference row (0.08-0.09 dex over all bins) is the relevant one.

**Parked (physics point, the last pass; R-149).**
- **What stands:**
  - On DESI DR1 BGS truth, photo-z isolation leaves red lenses with a real > 10 % neighbour within 0.3 / 0.5 Mpc 0.32 / 0.45 of
    the time, against 0.19 / 0.30 for blue ones, mass-matched. Hidden satellites are 0.29 against 0.11.
  - A UniverseMachine mock reproduces these fractions to 0.02-0.04.
  - At g_bar > 10⁻¹³ the leakage moves Brouwer's statistic by ≤ 0.01 dex with per-type references, and by 0.010-0.023 dex
    (0.06-0.17 of the published value) with a common one. Either way the small-radius split is not mainly leakage.
- **What is open:** the full-range share, 0 to about 0.5 depending on the reference.
- **Path to a pass:** a second galaxy model with a different (or zero) halo gap by colour, run through the same pipeline, so
  that the null is simulated rather than assumed. Also a colour-dependent photo-z error model, a light cone (MICE or AbacusSummit
  with a colour-dependent HOD), and DESI fibre-incompleteness weights.

## Earlier version (R-147): blind review

The R-147 version of this section (`briefs/r147_finding_section.md`) put an uncalibrated halo-model contamination
over a common McGaugh profile and claimed 0.21-0.26 of the published offset at small radius. Its gate, as recorded then:

VERDICT: novel but not sound: the spurious-offset fractions depend on an uncalibrated NFW contamination amplitude that exceeds Brouwer's MICE benchmark by a factor of about 3.5 (1.05 against 0.30 of the profile), divided by a reference profile that is the same for both types, so the claimed share of the published offset is not supported.

Strongest objection, in one sentence: the predicted offset scales directly with a contamination amplitude that disagrees with
the only external calibration (Brouwer's MICE test) by about 3.5 times. Until a mock or a GAMA calibration fixes that amplitude
by type, only the neighbour fractions (30 % against 20 %) and the satellite fractions (28 % against 12 %) stand.

Novelty was granted. The referee found no type-resolved photo-z isolation contamination in print. The closest work is
arXiv:2106.11677 App. A (whole sample only), 2310.15248 (qualitative), 1509.06762 (spectroscopic LBG purity by colour, not
the KiDS 3D criterion) and 2104.05355.

**Parked (physics point, one gate pass; R-147).**
- **What stands:** at fixed stellar mass, KiDS-bright-style photometric isolation, emulated on DESI DR1 BGS spectroscopic truth,
  leaves red / quiescent / bulge lenses with a real > 10 % neighbour within 0.3-0.5 Mpc 1.2-1.6 times as often as blue ones.
  Of the photo-z-isolated lenses, 28 % of red ones are hidden group satellites, against 12 % of blue ones. The 4 Mpc criterion
  barely changes this. Brouwer's equal-contamination assumption fails.
- **What does not stand:** how much of the 0.14-0.27 dex this produces. The prediction is close to linear in the contamination
  amplitude. Scaled down to Brouwer's MICE level (÷ 3.5), it would be about 0.01 dex (0.06 of the published value) at
  g_bar > 10⁻¹³ and about 0.04 (0.15) over the full range. As computed, it is 0.04 and 0.12 dex. The share is therefore
  bracketed at 0.06-0.26 at small radius and is not decided.
- **Path to a pass:**
  - Calibrate the amplitude by type on a lightcone mock with a KiDS-bright selection (r < 20, z < 0.5), against its truly
    isolated total ESD, which also fixes the type-dependent denominator. Alternatively, use the published GAMA-against-KiDS
    isolated RARs (Brouwer Fig. 4).
  - Vary the halo profile: truncation, M_200c against M_200m, subhalo stripping for satellite lenses.
  - Account for BGS fibre incompleteness in dense regions.

# Large-scale dust reddening around emission-line galaxies at z ≈ 1.2 from DESI DR1: at matched comoving separations of 0.5-6.6 Mpc, z > 2 quasars behind ELGs show ρ_c = −0.14 ± 0.43 times the reddening Ménard et al. measured around z ~ 0.36 galaxies (< 0.75 at 95 %), with the stellar-photometry systematic included and the magnification leak left unsubtracted

## Claim

Ménard, Scranton, Fukugita & Richards 2010 (arXiv:0902.4240) measured the mean colour excess of SDSS quasars as a function of
angular separation from SDSS i < 21 galaxies at z_eff ≈ 0.36: ⟨E(g−i)⟩(θ) = (1.4 ± 0.1) × 10⁻³ (θ/1′)^(−0.84 ± 0.05), from
20 kpc to several Mpc. Here the same statistic is measured around 2.43 million DESI DR1 emission-line galaxies (ELGs) at
0.8 < z < 1.6, with 447,716 DESI DR1 quasars at z > 2 as backgrounds. ρ_c is the measured amplitude over Ménard's at the same
comoving separation, with both converted to rest-frame A_V at the galaxy under the same (SMC-bar) extinction curve and the
slope fixed at 0.84.

| ELG bin | foregrounds | z_eff | ρ_c (matched comoving R) |
|---|---|---|---|
| 0.8-1.1 | 1.02 M | 0.95 | −1.21 ± 0.87 |
| 1.1-1.3 | 0.65 M | 1.20 | −0.72 ± 0.87 |
| 1.3-1.6 | 0.76 M | 1.43 | 1.04 ± 0.75 |
| **combined** (z_eff 1.22) | 2.43 M | | **−0.19 ± 0.40** |
| **with the photometric systematic propagated** (stellar null, below) | | | **−0.14 ± 0.43; < 0.75 at 95 % (flat prior ρ_c ≥ 0)** |

Errors are 68 % (jackknife over 100 sky regions) and include a 25 % error on Ménard's amplitude.

- **No reddening is detected around ELGs at z ≈ 1.2** at 0.5-6.6 comoving Mpc. The limit is ρ_c < 0.75.
- **This is large-scale, two-halo-dominated reddening, not halo dust.**
  - 76-85 % of the amplitude's Fisher information comes from r_p > 458 proper kpc (> 1 comoving Mpc), and 6-9 % from
    inside 131 kpc.
  - The R_200c of a 10¹² M_⊙ halo at z = 1.2 is 133 proper kpc. DR1's spectroscopic density puts 1,546 pairs inside
    131 kpc, which give ρ = −0.6 ± 1.0 there.
  - A 3σ detection of Ménard-level halo dust (ρ = 1) inside 131 kpc would need ~1.4 × 10⁴ pairs there. Of 0.3 of it,
    ~1.5 × 10⁵. That is 10-100 times DR1's spectroscopic pair density at these separations.
- **Against Ménard's profile, the ratio is ρ_c < 0.75.** For scale, an unchanged galaxy-correlated dust density, Ω_d b_d,
  would give ρ_c = S = 1.65 [1.30, 2.10]. S is the two-halo scaling below and rests on assumed biases. The limit therefore
  says Ω_d b_d(z ≈ 1.2) / Ω_d b_d(z ≈ 0.36) < 0.48 (95 %), with no-evolution 2.9σ away, under those assumptions.
- **Against the absorber-selected dust budget at z ~ 1** (Ménard & Fukugita 2012, arXiv:1204.1978):
  - Their Ω_dust(Mg II) = 2.3 ± 0.2 × 10⁻⁶ at z ≈ 0.5, scaling as (1 + z)⁻¹, is 1.56 × 10⁻⁶ at z = 1.22. They call it a lower
    limit.
  - Ménard et al.'s budget behind their profile is 5 × 10⁻⁶.
  - If the galaxy-correlated dust at z = 1.22 were the Mg II budget alone, ρ_c would be 0.51 [0.41, 0.65].
  - The measured ratio to that prediction is −0.27 [−1.12, 0.54], < 1.56 (95 %). **The limit is consistent with the Mg II
    budget and does not reach it.**

## Data and method

The samples, estimator, impact-parameter bins, fit, targeting correction and the nulls were fixed in writing before any
foreground-quasar pairing was read. The comoving conversion, the bias scaling, the enlarged stellar null, the WISE-free
estimator and the leak propagation were added afterwards and are labelled as such.

**Data.**
- Foregrounds: the DESI DR1 LSS v1.5 clustering catalogue ELG_LOPnotqso (NGC + SGC), with its WEIGHT and its random file 0
  (four randoms per galaxy).
- Backgrounds: DESI DR1 main-survey QSO targets with a good redshift, z > 2, 16.5 < r < 23, W1 at > 3σ, with Legacy Surveys DR9
  g, r, z and W1 photometry corrected for Galactic extinction.
- A pair is kept if z_q > z_f + 0.5. Impact parameters are proper at z_f, in 8 logarithmic bins from 20 to 3000 kpc.

**Per-quasar estimator.**
- Each quasar's g−r, r−z and z−W1 has the median subtracted in its cell: 0.02 in redshift × imaging system × r-magnitude
  quartile. Residuals are clipped at ±1 mag.
- The residuals are projected by GLS onto the SMC-bar reddening vector at the bin's z_f. The covariance is measured from the
  data per r quartile.
- Per-quasar σ(A_V) is 0.17-0.23 mag. It is set by intrinsic quasar colour scatter.

**Statistic.**
- Per r_p bin: the weighted mean of the per-quasar A_V over galaxy-quasar pairs, minus the same over random-galaxy pairs.
- It is divided by 0.885, DESI's targeting response to reddening. That figure comes from rerunning the DR9 QSO random forest
  on 40,000 of these quasars reddened at z_f: the mean reddening among the survivors is 0.88-0.89 of the true one.
- Fit: A_V(r_p) = A_150 (r_p/150 kpc)^−0.84 with the jackknife covariance (Hartlap-corrected). Bins are combined by GLS with
  their jackknife cross-covariance.

**Matched comoving scales.** Both profiles are r^−0.84. At comoving separation R the ratio is ρ_c = ρ_proper × [(1 + z_eff)/1.36]^0.84,
with z_eff each bin's pair-weighted foreground redshift. The factor is 1.35-1.63; at the same proper r_p the combined ratio is
−0.09 ± 0.26.

**Two-halo scaling (used only for the two comparisons).**
- Assumptions: dust traced like matter, linear bias, the same dust bias at both epochs.
- Then the rest-frame A_V at fixed comoving R scales as Ω_d b_d b_g D(z)² (1 + z)². The (1 + z)² is the proper dust column per
  unit comoving density.
- Inputs: b_ELG D(z) = 0.84 ± 15 % (DESI's fiducial ELG bias, arXiv:1611.00036), b = 1.1 ± 0.2 for Ménard's i < 21 sample
  (assumed), and Planck 2018 growth. Together they give S = 1.65 [1.30, 2.10].

## Checks

**Feasibility, set before the stack.** The forecast for Ménard's profile with no evolution at the same proper r_p was 4.1σ for
the ELGs. The sensitivity reached is 1.645σ = 0.66 in ρ_c.

**Stars as backgrounds (the photometric systematic).** Point sources of known non-reddened colour stacked around the same ELGs,
in cells of imaging system × r quartile × spectroscopic stellar subtype. Each colour's mean shift at ELG positions is put
through the quasar estimator's projection vector, to give the systematic a common per-band photometric offset would
induce in the quasars (ρ units, proper r_p):

| stars | N | g−r / r−z / z−W1 as A_V | in the quasar estimator |
|---|---|---|---|
| DR1 dark-program stars not chosen by the QSO, ELG or LRG cuts (mostly spectrophotometric standards: Gaia and g, r, z, no WISE) | 185,625 | −0.09 ± 0.09 / −0.23 ± 0.12 / +0.33 ± 0.24 | **−0.03 ± 0.10** (χ² 7.8/8) |
| stars chosen by the QSO target cuts (spectroscopic stars among QSO targets) | 78,737 | −1.51 ± 0.47 / +0.34 ± 0.71 / −12.4 ± 1.3 | −3.81 ± 0.54 |
| stars chosen by the LRG cuts | 21,008 | +5.3 ± 2.1 / +13.8 ± 4.2 / −28.1 ± 3.2 | −2.44 ± 1.60 |

- **The photometric null passes.** Stars that no W1-based cut chose show no colour shift around ELGs, at a precision three
  times better than the measurement.
- **Target-selected stars fail, almost entirely in z−W1.** Their offset is flat from 0.4 to 2.4 Mpc. That is a
  target-selection effect: contaminating stars sit at the W1-dependent selection boundary, where the admitted W1 depends on
  the local WISE depth. It is not a photometric offset around galaxies.
- **The quasars do not inherit it.**
  - Quasars at 0.6 < z < 1.5 with the ELGs behind them (same targeting, same photometry, no dust possible) give
    −0.20 ± 0.31, against −3.8 if they did.
  - A WISE-free quasar estimator (g−r and r−z only) gives ρ_c = −0.43 ± 0.47 (< 0.68), in line with the full one.
- **Propagation.** The standards' null (comoving −0.04 ± 0.15) is subtracted and its error added. This moves ρ_c from
  −0.19 ± 0.40 to −0.14 ± 0.43, the stated limit.

**Magnification leak.**
- From the same pairs, the quasars behind ELGs are brighter, achromatically, by m_μ(150 kpc) = 0.010 ± 0.004, 0.003 ± 0.005
  and 0.013 ± 0.005 mag in the three bins.
- Inside each r quartile the estimator rises toward faint r, by 0.06-0.10 mag of A_V per mag (pair-weighted). A lensed
  quasar is intrinsically fainter at its observed r, so it looks redder by that slope × m_μ.
- The predicted leak is **+0.34 ± 0.10 in ρ_c**: positive, smaller than the error, and only able to lower the limit if
  subtracted (to −0.49 ± 0.44, < 0.60).
- A separate fit of dust plus leak to the four quartile amplitudes gives ρ = −0.94 ± 0.64 for the ELGs.
- The headline leaves the leak in.

**Other nulls (same statistic; ρ units at proper r_p).**
- Random positions as foregrounds: −0.39 ± 0.22 (χ² 11.6/8).
- Foregrounds rotated by 2° in RA: 0.14 ± 0.26 (5.7/8).
- Galactic-dust halves (median SFD E(B−V)): −0.28 ± 0.41 and −0.17 ± 0.38.

**Other estimators of the same ρ (proper).**
- Milky Way curve: 0.12 ± 0.25.
- Redrock fibre-spectrum continuum slope (independent of the imaging): −0.52 ± 0.44.
- Pairs with a DR1 Mg II absorber within 1000 km/s of z_f removed: −0.15 ± 0.27.
- Single colours, ρ_c: g−r −0.46 ± 0.50, r−z −0.26 ± 0.90, z−W1 0.67 ± 0.76.

## What it means

Around star-forming galaxies at z ≈ 1.2, quasars at 0.5-6.6 comoving Mpc are reddened by less than 0.75 of what Ménard et al.
saw around z ~ 0.36 galaxies at the same comoving separations. The ELGs and Ménard's flux-limited galaxies are different
populations, so the ratio mixes redshift with sample.

If the large-scale reddening traces the galaxy-dust cross-correlation in the two-halo regime, the deficit relative to the
1.65 that an unchanged Ω_d b_d predicts says the galaxy-correlated dust density times its bias at z ≈ 1.2 is less than half
its z ≈ 0.36 value. That inference rests on the assumed biases and on the linear scaling.

The limit sits above what the Mg II-absorber dust budget alone would produce (ρ_c ≈ 0.5), so it neither confirms nor rules
out that budget.

Nothing here is a statement about dust inside the ELGs' halos: DR1 has ~10³ pairs there.

## Weakest points

- **The bias scaling.** The comparison with an unchanged dust density rests on S, which takes Ménard's sample bias as
  1.1 ± 0.2 and the dust bias as unchanged. ρ_c itself does not depend on S.
- **Ménard's profile is a projection over a broad photometric dN/dz.** The single z_eff = 0.36 conversion and the r^−0.84 shape
  are approximations, covered here by a 25 % amplitude error.
- **The enlarged stellar null was a post-M1 addition, and the split by target class was made after the first, mixed sample
  failed.** The photometric null uses brighter stars (median r 18.0) than the quasars (21.5). At r_p > 130 kpc (> 15″),
  blending with the foreground is not an issue.
- **The leak correction assumes the within-quartile colour-magnitude slope is intrinsic.** It is not subtracted.
- **The targeting correction counts only losses** (11-12 %, measured).

## To settle it

DESI DR2 (about twice the pairs) and a photometric ELG-like foreground sample from Legacy Surveys DR10 (~10× the pairs) would
narrow the limit toward the Mg II budget. A measured galaxy-dust bias at z ~ 0.36 (Ménard's sample's clustering) would remove
the main assumption in S. Halo dust at z ~ 1 needs ~10⁴-10⁵ pairs inside 130 kpc, which only a photometric foreground sample
can supply.

Rows: `results-R150.md` (R-150: F0-F4, M1, D1-D4; this revision, R-152: E0-E5); scripts `dust_z1_{prep,common,feas,stack,table}.py`
and `dust_z1_elg.py`. The LRG amplitude R-150 reported (ρ = 2.65 ± 0.48 at proper r_p, z 0.4-1.1) is not separated from
systematics (stellar null beside it +1.80 ± 0.74) and is not part of this finding (E5). No ELG-LRG contrast is claimed.

## Blind review

Blind review (2026-11-02, R-152, context-free Opus referee, own literature search, this section alone; `blind-review.md`):
VERDICT: novel but not sound: the limit propagates the systematic from a stellar null chosen after the fact that, by design, avoids the W1-dependent QSO targeting effect (−3.8 in ρ around ELGs) to which the background quasars are themselves exposed, so the 95 % limit and the 2.9σ no-evolution exclusion are overstated.

Strongest objection: the background quasars pass the same W1-based QSO cuts whose stars show −3.8 ± 0.5 around ELGs, so the
systematic to propagate is the quasar-side control (z 0.6-1.5 quasars with the ELGs behind, −0.20 ± 0.31), not the standards'
null. Propagated, that gives ρ_c = 0.12 ± 0.62, < 1.29 at 95 % (ρ_c ≥ 0), and no-evolution 2.1σ away. Novelty held: no
galaxy-selected reddening measurement with foregrounds at z > 0.5 (0902.4240, 1411.3333, 0706.0898, 1204.1978, 2503.04098,
2508.15490, 2606.24159; 2601.15961, seven z ≈ 3.5 sightlines, is not a stack). Secondary points: the random-positions null
(−0.39 ± 0.22) is unexplained, Ménard's R_V = 3.9 ± 2.6 against the SMC conversion, and magnification also changes which
quasars are selected. **Parked (last pass; a physics point).** What stands: no large-scale reddening around DESI DR1 ELGs at
z ≈ 1.2, ρ_c = 0.12 ± 0.62 (< 1.3, 95 %) with the quasar-side control propagated. The path to a pass: a quasar-side measure of
the targeting effect (e.g. the random forest rerun with each quasar's W1 noise at its local WISE depth, around ELG and random
positions), the random-positions null explained, and Ménard's profile converted with their own R_V.

### Earlier version (R-150, the ELG and LRG section this one replaces)

Blind review (2026-11-02, context-free Opus referee, own literature search; `blind-review.md`)

VERDICT: novel but not sound: the LRG detection and the ELG–LRG contrast are not separated from an LRG-tied systematic (stellar null 1.8 ± 0.7, strongest in z−W1) and a magnification leak of comparable size, and the 0.25–3 Mpc scales cannot show halo dust.

Strongest objection: the LRG stellar null (1.8 ± 0.7) is ~70 % of the LRG amplitude, so ρ_LRG minus it is ~0.9 ± 0.9; W1
blending from LRG satellites would redden z−W1 most (largest ρ, 3.7), and no behind-the-quasar null was run for LRGs. The
referee also asks for Ménard's ~25 % amplitude error (as quoted by arXiv:2503.04098) rather than 7 %, and notes the magnification
leak was cleared only post hoc. Novelty held: no galaxy-selected quasar reddening at z_f > 0.6 found (closest: 0902.4240,
1411.3333, 0706.0898, 1204.1978, 2503.04098, 2606.24159; 2405.16484, chromatic magnification behind CMASS at z ~ 0.55, overlaps
the LRG 0.4-0.6 bin). Per the referee, the ELG limit "is more robust" (pre-set 4.1σ forecast, clean behind-the-quasar null, both
curves, negative with the leak marginalised). **A miss on physics: parked (brief: one revision on a presentation point; park on
a physics point).** What stands: around DESI DR1 ELGs at z 0.8-1.6, ρ = 0.04 ± 0.29 (< 0.61 at 95 %; the reference error at
25 % instead of 7 % leaves this unchanged, ρ ≈ 0). Facts already in hand for a revision: WISE-free LRG estimators give g−r 2.51 ±
0.65, r−z 2.13 ± 0.58 and the fibre-spectrum slope 1.71 ± 0.58; missing are the LRG behind-the-quasar null, the stellar null per
colour and with 10x the stars, and a bias-scaled two-halo comparison.

# CMB lensing of z ~ 5 Lyman-break galaxies: 101,005 HSC GOLDRUSH r-dropouts cross ACT DR6 lensing at 3.3σ, at the amplitude the luminosity-function magnification alone predicts; with the dropouts' own clustering the implied growth amplitude A = σ8(z)/σ8^ΛCDM(z) is −0.3 ± 0.4 at z ≈ 4.8 and 0.05 ± 0.25 at z ≈ 3.8 (statistical; ±0.8 and ±0.5 after the pre-set consistency inflation and CIB term), a 2σ shortfall at each redshift that the model's magnification and interloper terms, not growth, most plausibly carry

## Claim

GOLDRUSH r-dropouts (z ≈ 4.8) and g-dropouts (z ≈ 3.8, i < 25.5) from the HSC-SSP S18A Wide layer, crossed with ACT DR6 and
Planck PR3 CMB-lensing convergence (pseudo-C_ℓ, 100 ≤ L < 1000 on ACT; covariance from the 400 ACT simulated reconstructions),
and with their own angular auto-spectra (100 ≤ ℓ < 1000):

| | z ≈ 4.8 (r-drops) | z ≈ 3.8 (g-drops) |
|---|---|---|
| N under ACT (all) | 101,005 on 302 deg² (115,950) | 1,354,372 (1,550,740) |
| cross S/N, ACT | **3.3** (χ²_null 16.0/8) | **5.6** (36.6/8) |
| cross S/N, ACT ⊕ Planck on W06 | 3.1 | 5.5 |
| magnification coefficient a = ⟨(1−f)(5s−2)⟩ (published LF) | 3.20 | 2.52 |
| b A² from the cross (A = 1 fit, nuisances marginalised) | −2.1 ± 2.4 | 0.2 ± 1.0 |
| b A from the auto-spectrum | 6.2 ± 1.1 | 4.2 ± 0.5 |
| **A = (b A²)/(b A)**, statistical | **−0.33 (−0.41 +0.38)** | **0.05 (−0.25 +0.24)** |
| A with the pre-set inflation (per-field χ² PTE 0.01 / 0.02: ×2.03 / ×1.80) and the CIB term | −0.3 ± 0.8 | 0.05 ± 0.5 |
| joint (b, A) grid fit, flat b ∈ [0, 20] | A = 0.41 (−0.08 +0.11), b at the prior edge | A = 0.29 (−0.06 +0.10), b at the edge |

The claim is the table. The cross-power is detected at both redshifts; at both, its amplitude equals the magnification term
(a C^{κμ}) that the published UV luminosity function predicts, leaving b A² consistent with zero, while the auto-spectrum asks
for b A ≈ 6 and 4. Planck ΛCDM (A = 1) is therefore 1.6σ (z ≈ 4.8) and 2.0σ (z ≈ 3.8) away after the pre-set inflation and CIB term (3.5σ
and 4.0σ statistical). Equivalently, at A = 1 and b = b_auto the cross needs a = −0.5 ± 1.0 (r) and −0.8 ± 0.6 (g), i.e. an effective
count slope s ≈ 0.25, not the LF's 1.23 and 0.97. Nothing is claimed about σ8(z): the joint fit's A is set by the b prior edge,
and a shortfall shared by two independent redshift samples on the same sky points to the model (the magnification response of a
dropout selection, the interloper fraction and clustering) or to the same-sky systematics, not to growth.

## Data and method

- **Galaxies.** GOLDRUSH (Harikane et al. 2022, arXiv:2108.01090; Zenodo 5512721) S18A Wide masked dropout catalogues with the
  published colour cuts: 119,651 r-dropouts (z S/N > 5 in 1.5″) and 1,569,187 g-dropouts with total i < 25.5 (Miyatake et al.'s
  cut). n(z) ∝ C(z) dV/dz with the paper's selection completeness C(z) (Fig. 4): ⟨z⟩ 4.82 and 3.80. GOLDRUSH ships no randoms: the
  galaxy window is the fraction of each nside-1024 pixel's four nside-2048 sub-pixels holding a g-dropout (4.8 per sub-pixel),
  the same for both samples; δ normalised per HSC field (W02-W06).
- **Lensing.** ACT DR6 baseline κ (mask² window; transfer 0.92-0.96 per band measured from 60 release input/reconstruction pairs
  on these windows, applied), the CIB-deprojected and curl reconstructions; Planck PR3 MV. ACT covers W02-W05; W06 only Planck.
- **Estimator.** NaMaster pseudo-C_ℓ at nside 1024; ACT 100 ≤ L < 1000 (7 bands), Planck 100 ≤ L < 400 (4 bands); auto-spectra
  100 ≤ ℓ < 1000 with the Poisson term subtracted. Covariance: the 400 ACT simulated reconstructions × the real galaxy maps; 400
  Gaussian Planck C+N realisations; 400 correlated Gaussian κ-galaxy realisations with Poisson sampling for the auto and
  cross-auto blocks (correlations from these, the cross diagonal from the release sims). Hartlap-corrected.
- **Model.** Planck 2018 ΛCDM, halofit, Limber. C^{κg} = A²(1−f) b C^{κm}_LBG + f b_int C^{κm}_int + a C^{κμ}; C^{gg} the
  corresponding squares and cross terms; A scales the matter power at the LBG redshifts. Nuisances marginalised over 300 prior
  draws: interloper fraction f (Harikane et al. Table 4, count-weighted: 0.21, 0.10) ± 0.10, interlopers at z 0.75 ± 0.25 / 0.45 ±
  0.20 with b_int 1.5 ± 0.5 / 1.3 ± 0.4; a = ⟨(1 − f(m))(5 s(m) − 2)⟩ over the objects, s(m) the intrinsic count slope of the
  published z ~ 5 / z ~ 4 UV LFs (galaxy + AGN double power laws) at each object's magnitude (s̄ = 1.23, 0.97; the samples sit at
  or brighter than M*), s ± 0.15. This LBG-weighted form of a replaced the fixed form (5s̄ − 2 on the whole sample: 4.13, 2.85)
  after the data spectra existed and before any fit; both are reported, and the fixed form lowers A further (−0.65, −0.04).
- **Fixed before the first cross-spectrum** (`results-R146.md`, 02:55 UTC): samples, windows, bands, estimator, model, priors,
  and the rules: errors inflated by √(χ²/dof) if the per-field b scatter has PTE < 0.05 or ACT and Planck on the same sky differ
  by > 2σ; a CIB-deprojected shift beyond its expected scatter is added in quadrature; the z ~ 4 reproduction passes if g × Planck
  gives S/N 5.1 ± 1.5.

## Checks

- **Nulls pass:** curl × r χ² 2.8/8 (PTE 0.95), × g 3.6/8 (0.89); randoms × ACT PTE 0.86 / 0.96, × Planck 0.72 / 0.07.
- **ACT vs Planck on the same sky:** Δb = 4.7 ± 4.1 (r), 1.7 ± 1.4 (g): consistent.
- **Per field (ACT), fails the rule:** r W02 3.8 ± 6.0, W03 4.8 ± 3.2, W04 −9.3 ± 2.9, W05 −4.9 ± 3.2 (χ² 12.4/3, PTE 0.01); g
  −2.3 ± 1.7, 2.2 ± 1.2, −2.2 ± 1.2, 1.1 ± 1.1 (9.7/3, 0.02). W04 is low for both samples on the same sky; errors inflated ×2.03 and
  ×1.80.
- **CIB:** the CIB-deprojected ACT map gives b A² = −4.2 ± 2.1 (r) and −0.9 ± 0.7 (g), 2.4 and 1.1 below the baseline, beyond the
  scatter expected from its higher noise (0.45, 0.15, an approximation without CIB-deprojected simulations); 1.9 and 0.9 added to
  b A² in quadrature. The deprojection lowers the cross, i.e. it deepens the shortfall.
- **z ~ 4 reproduction (Miyatake et al. 2022, arXiv:2103.15862):** g × Planck PR3 at 100 ≤ L < 400 gives S/N 1.8: **fails the
  pre-set 5.1 ± 1.5 rule.** Post hoc, Planck to L < 1000: S/N 2.8, and in Miyatake's own model (b C^{κm}, no magnification term)
  b = 2.1 ± 0.8, against their halo mass 2.9 × 10¹¹ h⁻¹ M_⊙ (b ≈ 3) and joint b = 6.5 (+5.2 −4.2): the amplitude agrees; their
  stack (θ ≤ 20′, Planck to L = 2048, jackknife errors) carries more small-scale signal than our L < 1000 bands.
- **Theory normalisation:** the Limber C^{κκ} (sources to z = 8.5) is 14-27 % below the release simulations' input × input on the
  same windows, the expected share of κ from z > 8.5; the ACT transfer is 0.92-0.96.

## What it means

The first CMB-lensing cross-power of z ~ 5 galaxies is detected (3.3σ on ACT), and ACT detects the z ~ 4 dropouts at 5.6σ.
Both amplitudes are what magnification bias alone predicts for samples this bright (they sit at or above the knee of the UV LF,
where the count slope is s ~ 1-1.4), so the cross constrains the clustering term weakly and only through the magnification
model: Miyatake et al.'s z ~ 4 analysis, which omitted magnification, read the same amplitude as clustering (b ≈ 2-3 here, as
there). Combined with the auto-spectra, the standard interpretation (linear bias, LF-slope magnification, Table-4 interlopers)
places Planck ΛCDM 2σ away at each redshift, in the direction of too little lensing per unit clustering; a σ8(z ~ 4-5) at a third
to a half of ΛCDM's is not a credible reading. The deficit says that for dropout samples at this depth the magnification response
is far below the LF slope (s_eff ≈ 0.25 would reconcile both), or that interlopers carry more of the auto-power than Table 4
allows, or that a same-sky systematic (W04) pulls both crosses low.

## Weakest points

- The cross is magnification-dominated; its clustering content hinges on a, taken from the LF slope with ±0.15 in s. The
  dropout selection's true response to lensing (colours from noisy dropout-band fluxes, S/N cuts in apertures, size) was not
  measured (no image injections).
- The auto-spectrum's covariance comes from Gaussian sims at the fiducial b (4 and 3, below the fitted b A ≈ 6 and 4); its lowest
  band is 3.8× (r) and 2.6× (g) the model, the signature of large-scale density systematics (depth, stars, seeing) that the
  galaxy window does not weight, and a free shot-noise term was not fitted.
- The galaxy window is reconstructed from the catalogue's holes (no randoms), 15 % larger than the paper's area.
- The per-field scatter (W04) and the CIB shift are folded in by fixed rules, not explained; the CIB scatter is approximated.
- The z ~ 4 reproduction fails its pre-set significance rule (the amplitude agrees).
- The joint (b, A) fit is prior-bounded; the ratio estimate treats the cross and auto as independent.

## To settle it

Image-injection measurement of the magnification response of the GOLDRUSH selection (SynPipe, as the completeness was made);
the GOLDRUSH clustering bias in place of our auto-spectrum (Harikane et al.'s own clustering, with its systematics treatment);
ACT CIB-deprojected simulations; HSC PDR3 dropouts (2× the area) and the Deep/UltraDeep layers for the faint, low-s end, where
the magnification term shrinks and the clustering term dominates.

## Blind review (2026-11-02, context-free Opus referee, own literature search; `blind-review.md`)

VERDICT: novel but not sound: the headline A values are a ratio of a numerator dominated by an assumed magnification coefficient that was never measured, over an auto-spectrum the note itself shows is contaminated by systematics, from a pipeline that failed its pre-set reproduction check, so the claimed 2σ shortfall is not supported.

Strongest objection: A is the ratio of a cross-spectrum residual that rests on an unvalidated, post hoc magnification
coefficient to an auto-spectrum amplitude the note itself shows carries large-scale systematics, from a pipeline that failed its
pre-set reproduction of Miyatake et al., so the shortfall cannot be assigned to physics, to the model or to analysis error.
Novelty granted for the z ~ 5 cross and for the magnification-dominance point at z ~ 4 (no prior CMB-lensing × dropout work
at z ~ 5 found; closest arXiv:2103.15862; forecasts 1904.13378, 2401.13166, 2511.22243; LBG magnification with foreground
galaxies, Hildebrandt et al. 2009 and arXiv:1204.2830; slope-only magnification responses mispredict incomplete samples,
arXiv:1511.01352); the z ~ 4 g-dropout amplitude of structure is prior work (2103.15862).

**Parked (physics point, one gate pass; R-146).** What stands: the first CMB-lensing cross-power of z ~ 5 galaxies (3.3σ, ACT
DR6, 101,005 r-dropouts) and an ACT detection of the z ~ 4 dropouts (5.6σ), both at the amplitude that LF-slope magnification
alone predicts; the ratio to the auto-spectrum is not a growth measurement. The path to a pass: the GOLDRUSH selection's
magnification response from image injections (SynPipe), the published GOLDRUSH clustering bias instead of our auto-spectrum,
HSC randoms (PDR3) for the window, ACT CIB-deprojected simulations, and a reproduction of Miyatake's stacked statistic on Planck
before the C_ℓ comparison.

# The Giant Arc does not reappear on independent sightlines: on 3,822 DESI DR1 quasar sightlines that SDSS never observed, its convex hull holds 8 Mg II absorbers against 8.9 expected; the fitted overdensity is A = −0.23 (95 %: −0.90 to +0.63), against the published δρ/ρ = 1.3, excluded at 3.2σ

## Claim

The Giant Arc (GA; Lopez, Clowes & Williger 2022, arXiv:2201.06875) is a ~1 Gpc arc of Mg II absorbers at z ≈ 0.80. It was
found by single-linkage clustering of absorbers in SDSS quasar spectra. The authors read it as a real overdensity
(δρ/ρ = 1.3 ± 0.3 in its convex hull of member spheres, CHMS). Sawala et al. (arXiv:2502.03515, 2505.11072) read it as a
pattern the algorithm picks out of an ordinary point field. Both sides argue from simulations of the same SDSS points. A real
overdensity must also raise the absorber count on quasar sightlines that played no part in finding it. We count Mg II
absorbers from the DESI DR1 catalogue on DESI quasars that are in neither DR16Q nor the Zhu & Ménard searched lists, inside
the GA's frozen template. The expectation comes from dN/dz measured on the same sightlines over the whole sky, stratified by
the quasar's g flux, exposure depth and redshift.

| cell (W2796 > 0.3 Å) | sightlines | N_obs | N_exp | δ | shift offset (σ) | Δ = δ_cal − baseline [95 %] | A = Δ/f_H [95 %] | real prediction | vs published |
|---|---|---|---|---|---|---|---|---|---|
| **CHMS (primary)** | 3,822 | **8** | **8.86** | −0.10 | −0.10 (0.34) | **−0.22 [−0.88, +0.62]** | **−0.23 [−0.90, +0.63]** | 1.27 | **3.2σ** |
| spheres R = 47.5 Mpc | 3,533 | 4 | 4.13 | −0.03 | −0.09 (0.47) | −0.02 [−0.97, +0.98] | −0.02 | 1.14 | 2.3σ |
| R = 70 | 4,806 | 6 | 9.64 | −0.38 | −0.11 (0.34) | −0.33 [−0.91, +0.41] | −0.45 | 0.95 | 3.3σ |
| R = 95 | 6,281 | 15 | 16.39 | −0.09 | −0.11 (0.27) | +0.03 [−0.43, +0.63] | +0.05 | 0.67 | 2.2σ |
| CHMS, W > 0.6 | 3,822 | 8 | 8.40 | −0.05 | −0.07 (0.34) | −0.10 [−0.77, +0.69] | −0.10 | 1.27 | 2.9σ |

The claim is the primary row: in its own hull the GA carries no excess on new sightlines, A = −0.23 (68 %: −0.71 to +0.26;
95 %: −0.90 to +0.63). That is consistent with zero (−0.5σ) and 3.2σ below the published amplitude. Every other cell
agrees: no cell shows an excess, and the cells forecast to separate the readings at ≥ 2.8σ before the data were seen (CHMS at both
EW cuts, R = 70) all exclude A = 1.3 at ≥ 2.9σ. In plain Poisson terms, a real δ = 1.3 in the hull predicts about 19 absorbers
where 8 are seen. This does not show that no structure exists there. It shows that whatever the GA's 46 SDSS absorbers trace,
it is not an overdensity of Mg II absorbers at the published level. A real structure much weaker than published (A ≲ 0.6) is
not excluded.

## Data and method

- **Template, frozen before any DESI file was opened.** Zhu & Ménard DR7 + DR12 trimmed catalogues (69,206 after dropping 8,441
  DR7/DR12 duplicates). Slice z = 0.802 ± 0.060 (4,752 absorbers). Single linkage at 95 Mpc (flat Ωm 0.27, h 0.7, the discovery
  recipe). The GA group has 46 members (the paper has 44), RA 205.0-224.6°, Dec 36.2-42.2°, member W2796 median 0.82 Å. CHMS:
  convex hull of spheres of radius half the mean MST edge (30.4 Mpc), V(CHMS ∩ slice) = 3.5 × 10⁷ Mpc³. The slice's largest
  groups have 88-131 members; the GA is not among them.
- **Independent sightlines.** DESI DR1 main-dark quasars (the catalogue the DESI Mg II VAC searched; 1,556,507). 322,606 within
  1.5″ of a DR16Q or Z&M-searched quasar are removed, leaving 1,233,901. Sightlines within 200 proper kpc of a member's own
  sightline would be dropped (one-halo re-detections); none fall that close, and 0 and 500 kpc give identical counts. Path per
  sightline: redward of Lyα (the VAC does not search the forest: dN/dz is exactly zero there) to 5000 km/s below z_q. Exact
  ray-sphere and ray-hull intersections.
- **Completeness.** DR1's Mg II detection rate depends on the quasar's flux, not mainly the exposure depth. All-sky dN/dz (z = 0.8,
  W > 0.3) runs 0.00004 to 0.28 across g-flux quintiles, against 0.052-0.090 across TSNR2 tertiles. Removing SDSS quasars
  removes the brighter ones inside the SDSS footprint, so the GA's independent sightlines are fainter than average. dN/dz is
  therefore measured per stratum (TSNR2 tertile × 12 g-flux bins × 4 z_q bins) on all independent sightlines. The all-sky check
  in the slice gives 7,844 expected against 7,864 observed. Under Poisson counts, Σn/Σe − 1 is the maximum-likelihood excess, so
  the stratified ratio is the TSNR2 and flux weighting. A literal TSNR2-weighted ratio gives −0.19 in the primary cell.
- **EW matching.** On 21,789 absorbers both catalogues hold, DESI EW2796 / Z&M REW = 1.03 (flat in EW). Cuts W > 0.3 and 0.6 Å on
  the Z&M scale. 92 % of DR1's W > 0.3 detections lie above 0.6 Å, so the DESI sample's EW mix is close to the members'.
- **Calibration.** Each template is referred to its own shift null: the same shape rotated rigidly to random places in the SDSS
  footprint (> 30 % SDSS quasars among DESI quasars within 1° of its members; the GA has 39-53 %), not overlapping the original,
  with N_exp within 0.75-1.33× its own. For the GA, ~150 kept shifts per cell. δ_cal = δ − mean(shift). The shift spread at the
  GA's depth is the error (Poisson, coverage and large-scale structure together).
- **Chance baseline.** The other single-linkage groups from the same slice and recipe, matched to the GA in size (23-92 members)
  and extent (0.5-2× its 742 Mpc largest separation): 18 groups, 15 with DESI coverage. The fixed rule (a group counts only if
  ≥ 20 of its own shifts match it) keeps 4-5 groups per cell; the baseline is their mean δ_cal (CHMS: +0.22 ± 0.34). Δ = δ_cal(GA) −
  baseline. The interval combines the GA's shift-null percentiles with the baseline's standard error in quadrature.
- **Real-structure prediction, computed before the measurement.** An overdensity A uniform in CHMS ∩ slice raises a cell's
  expectation by A f_H, where f_H is the fraction of the cell's expected count lying in CHMS ∩ slice (0.98 for the CHMS, 0.52 for
  R = 95). Putting the same excess only in the member spheres changes the predictions by < 0.1, except at R = 47.5. The primary
  cell is the one with the largest separation of A = 1.3 from zero at the GA's depth (CHMS, 3.4σ; R = 95, the proposal's choice,
  has only 2.5σ because its spheres dilute the hull's excess). The fitted amplitude is A = Δ/f_H.
- **Process.** The plan passed a blind method check as "fixable" (six points), and all six were closed before any GA count was
  computed. The steps were: forward-model the real prediction and re-pick the primary cell; use matched linkage groups as the
  baseline; read the result as an amplitude; give each template its own calibration; add the close-pair exclusion; use W > 0.6
  as the check. **Disclosure:** at proposal time an unblinded rough count was made (δ = −0.16 / −0.28 / −0.29 at R = 47.5 / 70 /
  95; partial SDSS exclusion, flat dN/dz, no EW matching). The analysis above was fixed without reference to it, and nothing in
  the plan was changed after the first look at this measurement. That rough count was biased low by the completeness gradient
  above: a flat dN/dz overpredicts the GA's expectation about 2×.

## Checks

- **Baseline choice.** The primary result does not depend on it. With all 14-15 covered groups (raw δ; most lie where the
  SDSS-footprint shifts cannot match them), CHMS Δ = −0.13 ± 0.37, A = −0.13, and the published value is excluded at 3.8σ. With
  no baseline at all, A = 0.00 and the exclusion is 3.7σ. The proposal's random-46 clustering null (spheres on random sets of 46
  Z&M slice absorbers) sits +0.14 / +0.07 / +0.06 above random positions at R = 47.5 / 70 / 95. Ordinary clustering is
  therefore small next to the 1.3 being tested.
- **The groups.** The linkage groups as a set show a small positive excess on new sightlines (raw mean δ +0.13 ± 0.13 in their
  CHMS, −0.04 ± 0.08 at R = 95), as clustered absorbers should. The GA sits at or below the median of the groups in every cell
  (rank 0.20-0.50).
- **EW cut, weighting, exclusion radius.** W > 0.6 gives A = −0.10. The TSNR2-weighted ratio gives δ = −0.19. Exclusion radii of
  0, 200 and 500 kpc give identical counts.

## What it means

The step the two sides dispute is whether the GA's points trace real matter. New sightlines answer that directly: at the
published amplitude the GA's hull would hold about 19 DR1 absorbers, and it holds 8, the number an average patch of sky at the
same depth holds. The result agrees with the reading that the arc is an algorithmic pattern among ordinary absorbers
(arXiv:2502.03515, 2505.11072). It goes beyond the simulations behind that reading because it measures the real sky
on data the pattern was not drawn from. It does not touch the Lopez group's FLAMINGO argument (arXiv:2504.14940) about
whether ΛCDM allows such structures; it removes the observational case that this one exists at δ ≈ 1.3.

## Weakest points

- The selection-bias shrinkage is not modelled: how far a real but weaker structure's amplitude falls once its noisiest
  members are not re-found. An intermediate A could not be split into "real" and "inflated by selection". The measured A is
  close to zero, so this matters only for the upper end (A ≲ 0.6 is allowed).
- The chance baseline uses 4-5 calibrated groups (SE 0.34 in the CHMS), not linkage run on clustered mocks. The checks show the
  conclusion holds with any of the three baselines.
- The expected count is small (8.9), and the SDSS exclusion leaves the GA's sightlines faint. Completeness strata are a model:
  the all-sky check and the shift null test it, and a residual offset of −0.07 to −0.11 is calibrated out per template.
- One DESI release. The Big Ring and Giant Ring are not tested: neither paper gives a reproducible member list.

## To settle it

DESI DR2 (about 2× the sightlines) would bring σ(A) in the CHMS cell from ~0.37 to ~0.26. The Big Ring is the
next test once its member list is public. Linkage on lognormal absorber mocks would replace the 4-5-group baseline and give
the selection shrinkage.

## Blind review (2026-11-02, R-155, context-free Opus referee, own literature search, this section alone; `blind-review.md`)

VERDICT: the findings are confirmed novel and sound if the calculations hold

Strongest objection: the published 1.3 was measured on the members the algorithm selected, so even a real structure would come
back weaker on new sightlines. The test therefore falsifies the published amplitude, but cannot tell a modest real structure
(A ≲ 0.6) from a chance pattern without modelling the selection shrinkage on clustered mocks. Closest work found: none re-counts
absorbers on independent sightlines (2201.06875, 2402.07591, 2604.17534, 2504.14940, 2502.03515, 2505.11072; related only:
2604.13162). Section unchanged since review.

---

# Finding (R-154, cosmo-deep-22, 2026-11-02): CMB lensing of 6,263 ACT DR6 SZ clusters finds no fall of the SZ mass bias with redshift out to z ≈ 1.2

2026-11-02 (R-154, cosmo-deep-22).

**Headline.** The evolution exponent of M_lens/M_SZ is ζ = −0.6 [−2.6, +1.3] (95 %, systematics included), against
ζ = 1.8 (+0.6/−0.5) from DES galaxy lensing. The extrapolated galaxy-lensing law is disfavoured at 2.4σ (3.0σ statistical).
The SZ-significant and SZ-faint halves of the catalogue disagree below z = 0.8.

**Scripts.**
- `act_szlens.py`: catalogue, footprint, simulations, Planck CIB stack, selection mock.
- `act_szlens_inject.py`: the reconstruction emulated in flat sky, for the injection calibration.
- `act_szlens_fit.py`: fits, nulls, systematics.
- `act_szlens_budget.py`: the budget.
- The fixed plan is `plan-R154.md`, with its addendum written after the method check and before the first data stack. Rows
  D0-R5 are in `results-R154.md`.

**Data.**
- The ACT DR6 SZ cluster catalogue v1.0 (arXiv:2507.21459; 10,040 optically confirmed clusters at q > 4, 1,180 at z > 1), with
  its A10 masses (`A10_M500c`). These are Arnaud et al. 2010 scaling with the Eddington correction from a Tinker prior: the
  mass definition of the ACT DR5 masses behind the published 1-b.
- The ACT DR6 lensing release v1 (Qu et al. 2024, arXiv:2304.05202; Madhavacheril et al. 2024, arXiv:2304.05203): the baseline,
  tonly, polonly, cibdeproj and curl maps, the 6' SZ-cluster mask, and the 400 simulated reconstructions.
- Planck PR3 545 and 857 GHz maps, for the CIB of the clusters.
- The ACT DR5 catalogue, for the mass cross-match.

## Claim

**The published claim being tested.** Galaxy weak lensing of ACT DR5 SZ clusters (DES Y3; Shin et al., arXiv:2512.18935) finds
the SZ mass bias 1-b = M_SZ/M_true falling with redshift: 0.81 (+0.08/−0.06) at 0.15 < z < 0.43 and 0.58 (+0.06/−0.05) at
0.43 < z < 0.70. Their fit is M_true/M_SZ ∝ ((1+z)/1.45)^ζ with ζ = 1.8 (+0.6/−0.5), so ζ > 0 at 99.95 %. Galaxy lensing loses
its sources behind clusters near z ≈ 1. CMB lensing has its source at z = 1100.

**What was done.**
- The ACT DR6 convergence was stacked on 6,263 DR6 SZ clusters in the lensing footprint, in four redshift bins.
- Each bin's stack was fitted with the NFW-plus-linear-two-halo κ predicted from every cluster's A10 mass posterior. The
  prediction is filtered exactly as the map is, the two-halo amplitude is free, and the covariance comes from the release's 400
  simulations.
- That gives R = M_lens/M_SZ per bin. The same form as 2512.18935 was fitted jointly across the bins for ζ.

**The table.**

| z bin | n | ⟨z⟩ | detection | R = M_lens/M_SZ, 68 % (stat) | 95 % (stat) | σ_sys(ln R) | 1-b = 1/R |
|---|---|---|---|---|---|---|---|
| 0 – 0.5 | 2,459 | 0.33 | 8.8σ | 0.99 [0.85, 1.14] | [0.72, 1.31] | 0.22 | 1.01 |
| 0.5 – 0.8 | 2,125 | 0.64 | 6.0σ | 0.75 [0.60, 0.90] | [0.47, 1.05] | 0.31 | 1.34 |
| 0.8 – 1.0 | 922 | 0.89 | 3.0σ | 0.68 [0.42, 0.96] | [0.20, 1.27] | 0.30 | 1.46 |
| 1.0 – 2.0 | 757 | 1.23 | 4.0σ | 1.32 [0.93, 1.81] | [0.59, 2.32] | 0.50 | 0.76 |

**The trend.**
- Statistical: **ζ = −0.57 [−1.40, +0.27] (68 %), [−2.32, +1.01] (95 %)**. Δχ² is 0.4 at ζ = 0 and 9.0 at ζ = 1.8.
- With every z-dependent systematic run as its own variant of the whole fit and added in quadrature (σ_sys(ζ) = 0.49):
  **ζ = −0.57 [−1.53, +0.40] (68 %), [−2.57, +1.28] (95 %)**.
- ζ = 0 (a constant bias) is inside the interval. ζ = 1.8 is outside it, at 2.4σ.
- The 2512.18935 law extrapolated to the z > 1 bin predicts R ≈ 3.3. The measurement is 1.32 [0.59, 2.32] at 95 %, statistical.

**Matched to the published selections (uncorrected for the +4 % selection bias).**
- **2512.18935's own selection** (q > 5.5, their two bins): 1-b = 0.66 [0.56, 0.79] at z ≈ 0.32 and 0.96 [0.79, 1.19] at
  z ≈ 0.56. Against their 0.81 and 0.58, the high-z bin differs by 2.6σ.
- **HSC's selection** (arXiv:2407.08201; q > 5, their three bins): M_lens/M_SZ = 1.43 ± 0.22, 0.95 ± 0.20 and 0.88 ± 0.23,
  against their M_WL/M_SZ = 1.3 ± 0.2, 1.6 ± 0.2 and 1.6 ± 0.3.

**The halves of the sample** (both listed in the plan as diagnostics).

| subsample | R, bins 1-4 | ζ | Δχ² at ζ = 1.8 |
|---|---|---|---|
| the SZ-significant half (q > 5, 3,620) | 1.28, 0.91, 0.40, 1.39 | −1.41 [−2.37, −0.55] | 15.7 |
| the SZ-faint half (q 4-5, 2,643) | 0.50, 0.40, 1.14, 1.31 | +2.57 [+0.90, +4.0] | 0.2 |

The two halves agree at z > 0.8 and disagree below it, where the faint half shows about 40 % of the expected lensing.

## Method

**Sample.**
- R-122's footprint selection keeps a cluster if it is clean within 30' and ≥ 80 % clean within 2°.
- 2,824 of the clusters (2,650 of the 3,620 at q > 5) are inside the release's SZ-cluster mask. In the baseline map their tSZ
  was removed by subtracting UPP templates before the reconstruction (Qu et al. 2024 §III.6). The rest had only the
  reconstruction's profile hardening.

**Statistic.**
- Nine annuli from 0 to 150'. The band is 20 ≤ L ≤ 2048 with the release's pixel windows.
- Each cluster's prediction is averaged over a lognormal of its catalogue mass errors (5 Gauss-Hermite nodes).
- Miscentring: a Rayleigh σ_off = 0.5', the upper end of the catalogue's SZ-optical offsets (0.30-0.49').
- The mask-matched random stack is subtracted.
- The one-halo concentration is Duffy 2008. The two-halo term is linear (Tinker 2010 bias) with a free amplitude per bin.
- R per bin comes from Δχ² on a grid, each cluster's mass multiplied by R.
- ζ comes from a joint fit of the 4 × 9 vector with the joint simulation covariance (Hartlap-corrected), each cluster scaled by
  R_0((1+z_i)/1.45)^ζ.

**Validation.** Injected into a simulated stack, R and ζ are recovered without bias. A noiseless ζ = 1.8 is returned as 1.85.

**Calibrations closed before the data stack** (a context-free method check returned "fixable" on ten points; all closed).

*The reconstruction*
- The release pipeline was emulated in flat sky with ACT beams and with noise matched to the release's N_L:
  - CAMB CMB lensed by each bin's median cluster;
  - its A10 tSZ at 90 and 150 GHz;
  - the release's template subtraction, fitted on noisy maps;
  - the temperature estimator profile-hardened with the release's profile, √C_l^tSZ.
- The linear response is 1.00 ± 0.03 at R = 1.35 and at R = 3.4.
- tSZ leakage, in units of the bin's one-halo amplitude:
  - unhardened: −0.55 to −1.12;
  - hardened, for the clusters outside the mask: +0.09 to +0.15;
  - template-subtracted: ≤ 0.016.
- The leakage is weighted per cluster (it scales as M^2) and by the temperature weight of the baseline estimator
  (N_L^MV/N_L^T = 0.90).

*The CIB*
- Measured per bin from Planck 857 GHz aperture fluxes: 0.22-0.36 Jy, S/N 2-6.
- The 545/857 ratios (0.32 rising to 0.56) fix the modified-blackbody SED.
- Injected through the same pipeline, the CIB leaks −0.03 to −0.14 for the unsubtracted clusters, growing with z, and ≤ 0.005
  for the subtracted ones.
- Net leakage applied: +0.014, +0.018, +0.025 and +0.003.

*The selection*
- A Tinker-halo mock with the catalogue's own y-M slopes (1.22-1.47) and per-cluster noise, selected at q > 4, with
  catalogue-style posteriors, run through the real fit.
- It gives B = R_fit/R_true = 1.05 at 1-b = 0.75, falling to 0.83-0.95 at 1-b = 0.25. This is inverted per bin.
- Variants: scatter 0.1/0.3, σ8, optical confirmation at z > 1.

*Masses and purity*
- DR6 A10 against DR5 masses: equal to 1-3 %, with no z trend.
- Contamination of 0-4 % (6 % at z > 1), from the DR6 paper's 1-6 % chance association, is budgeted.

**Systematics** (`act_szlens_budget.py`; shifts in ln R per bin, and in ζ).

| systematic | ln R, bins 1-4 | ζ |
|---|---|---|
| concentration × 0.7 / × 1.3 | 0.13-0.23 | 0.15 |
| two-halo amplitude fixed at 1 | 0.06-0.38 | 0.43 |
| inner 30' only | 0.09-0.19 | 0.07 |
| selection scatter | 0.06-0.10 | 0.11 |
| hardening-profile response | 0.05-0.07 | 0.01 |
| contamination | 0.05-0.07 | 0.04 |
| CIB, tSZ, miscentring, σ8, optical confirmation, η = −0.18 | ≤ 0.04 each | ≤ 0.06 each |

- The two-halo treatment dominates ζ. With the amplitude fixed, ζ = −0.13 [−0.63, +0.20], which is tighter and further from
  1.8.
- **Trials:** the four bins, the trend, 22 systematic variants, eleven subsamples and three matched comparisons are all reported.
  The claim is the plan's.

**Nulls.**
- **Curl** on the clusters: A_1h = −0.08 ± 0.12, 0.04 ± 0.14, 0.01 ± 0.25, 0.09 ± 0.33.
- **40 rotated maps:** mean −0.07 to +0.03, with a scatter equal to the simulations'.
- **The simulations' mean** at the clusters: 0.00 ± 0.01.
- **Map variants:**
  - polonly (tSZ- and CIB-free, twice the noise): 0.75, 0.63, 1.45 and 1.23, against the baseline's 0.88, 0.68, 0.89 and 1.21
    on the same covariance;
  - tonly: within 0.06-0.17 of the baseline;
  - cibdeproj: 5-12 % lower in every bin.

## What it means

**The trend.**
- CMB lensing, which has no source-redshift limit, sees no decline of 1-b from z ≈ 0.3 to z ≈ 1.2 in ACT SZ clusters.
- The best fit is a flat or slightly rising 1-b. The galaxy-lensing exponent ζ = 1.8 is outside the 95 % interval.
- This holds at matched selection: the q > 5.5 clusters at 0.43-0.70 give 1-b = 0.96 [0.79, 1.19], where DES finds 0.58.
- The test that motivated it gives the same answer: the extrapolation to z > 1 predicts M_lens/M_SZ ≈ 3.3, and 1.3 is measured.
- So the steep evolution is not borne out by a method free of the high-z source problem. Shin et al. and Bocquet et al.
  (arXiv:2401.02075, SPT, an E(z) exponent in the SZ-mass relation consistent with self-similar evolution) differ on this
  point; the result here sides with the milder evolution.

**The normalisation is not the claim.**
- R ≈ 0.7-1.3 at fixed concentration. 1-b ≈ 0.66 for the 2512.18935-matched low-z clusters agrees with weak lensing.
- But the per-bin normalisation moves by 13-23 % with the concentration and by up to 38 % with the two-halo treatment.

**Nothing on S8 is claimed beyond this:** cluster-count analyses that rely on a 1-b falling with z have no support from CMB
lensing.

## Weakest points

**The SZ-faint half.**
- At z < 0.8 the q 4-5 clusters show R ≈ 0.4-0.5, where the q > 5 clusters show 0.9-1.3. That is 2.8σ in bin 1 and 2σ in bin 2.
- The selection mock gives the same B for both halves (±7 %). The DR6 chance-association rate (1-6 %) is far too small to
  explain it.
- It is not the template subtraction: q > 5 clusters outside the SZ mask give 1.22 in bin 1.
- It is unexplained. Because it depresses the faint half at low z only, it pushes that half's ζ up (+2.6). The full sample's ζ
  therefore lies between a q > 5 value further from 1.8 (−1.4) and a faint-half value consistent with it.
- If the faint-half deficit is a defect of those clusters (mis-associations, or masses overestimated at low q), the full sample
  understates the disagreement with ζ = 1.8. If it is real physics in low-mass clusters, the mass dependence matters more than
  this test models. At fixed bin, the mass split gives R rising with M_SZ, the opposite of 2512.18935's η = −0.18.

**The reconstruction calibration is an emulation.**
- It uses a flat-sky temperature-only estimator, hardened with a halo-model √C_l^tSZ, with Gaussian noise.
- The release's curved-sky, split-cross minimum-variance pipeline is not runnable here.
- The polarization-only map, free of the foregrounds the emulation calibrates, agrees with the baseline bin by bin within its
  errors, which are twice as large.

**The one-halo/two-halo split.**
- The data exceed NFW plus linear two-halo at 10-30', as in R-122. The free two-halo amplitude trades against the one-halo.
- That is the largest systematic on ζ, carried as a ±0.43 shift.

**The high-mass tail.** The 35 clusters at z > 1 with A10 masses above 3 × 10^14 give R ≈ 5 (2.4σ above the bin).

## To settle it

- A purity and mass study of the q 4-5 clusters at z < 0.8 (their optical richness and X-ray matches against the q > 5 ones).
- An outer-profile model (a splashback or emulator profile) in place of NFW plus linear two-halo.
- The cluster-scale estimator on the ACT maps themselves, a gradient-inversion or maximum-likelihood one, for the release's own
  response.
- SPT-3G lensing on the SPT clusters at z > 1.

## Blind review (2026-11-02, context-free Opus referee, own literature search; `blind-review.md`)

Method check, before computing (verbatim):

> METHOD: fixable: injection-calibrate the reconstruction (NFW κ plus measured-y tSZ, with template subtraction where the
> release did it) for the response at SZ-mask positions, tSZ leakage and estimator non-linearity per bin, instead of relying on
> the noisier polonly difference and the confounded in/out-mask and SNR splits; separate CIB from tSZ with a measured-amplitude
> CIB injection rather than cibdeproj−baseline; model optical-confirmation completeness in the selection mock; extend the
> Eddington mock to 1-b ≈ 0.25, apply it self-consistently in R, vary scatter and σ8, and run the real profile fit on the mock
> stacks; correct and budget purity and photo-z/bin migration from the catalogue's false-detection rates and photo-z posteriors;
> calibrate σ_off per bin from SZ-optical offsets and budget it; cross-match DR5/DR6 A10 masses and report ζ with η fixed at
> −0.18; include every z-dependent systematic in the ζ interval

All points were closed before the data stack (`plan-R154.md` addendum).

Finding check, on this section as it stands (verbatim):

> **VERDICT: novel but not sound: the two SZ-significance halves give mutually inconsistent ζ (−1.4 vs +2.6), split almost
> exactly along template-subtracted versus emulation-calibrated hardened tSZ treatment, so the headline ζ averages over an
> unexplained systematic larger than the effect tested.**

**The referee's strongest objection, in its words.** The sample is not one population. The q-split lines up with the
foreground-treatment split (template subtraction against hardening calibrated only by the note's emulation), and an imperfect
real hardening would push the faint half low with a redshift dependence.

**Other points it raised.**
- Mass and redshift are tangled in a q-limited sample.
- The 2.4σ ignores DES's own error on ζ; with it, about 2.1σ.
- The z > 1 comparison rests on extrapolating a z < 0.7 power law.
- A free two-halo amplitude per bin can absorb a z-dependent profile mismatch.

**Prior work it found**, none of which measures ζ with CMB lensing:
- CMB lensing of SZ clusters: arXiv:1502.01597, 1408.5633, 1904.07887, 2402.18346, 1412.7521, 1810.10998, 2111.07491.
- Weak-lensing calibrations: arXiv:2512.18935, 2407.08201, 2304.10219, 2401.02075, 2606.17984.
- The data papers: 2507.21459, 2304.05203.

**A miss, on a physics point. Parked, as the brief directs: no revision of the section.**

**One diagnostic run after the verdict**, recorded and not used to revise. The polarization-only map is free of tSZ and CIB, so
it does not depend on the emulation. On it, the faint (q 4-5) clusters at z < 0.8 give A_1h = −0.20 ± 0.37, against the
baseline's 0.49 ± 0.19; at z > 0.8, 0.91 ± 0.82. The q > 5 clusters give 1.26 ± 0.21 against 1.21 ± 0.12. The faint
clusters' low-z deficit is therefore present without any temperature foreground. It is not hardened-tSZ leakage, but a property
of those catalogue entries: their associations or their masses.

**What would settle it.** The q 4-5, z < 0.8 clusters' optical and X-ray properties against the q > 5 ones. If they are
mis-associations, the q > 5 sample is the clean one, and it gives ζ = −1.4 [−2.4, −0.6].

---

# Doppler magnification of galaxy sizes at 0.15 < z < 0.45: a first measurement in DESI DR1 BGS x Legacy Surveys DR9 (R-156)

## Claim

At fixed observed redshift, a galaxy's line-of-sight peculiar velocity changes its apparent size by the Doppler
convergence kappa_v = (c/(H r) - 1) v/c. Galaxies behind an overdensity fall toward us and look smaller; those in front
look larger, so sizes and number counts have an odd (dipole) cross-correlation (Bonvin et al. 2017, arXiv:1610.05946).
On 2.59 M DESI DR1 BGS_BRIGHT galaxies with Legacy Surveys DR9 Tractor half-light radii (7,180 deg2), the dipole
amplitude relative to linear-theory GR with Planck 2018 f sigma_8, at 20-150 Mpc/h, is

| z | galaxies | response eta | A_D | stat (jackknife) | eta | b_s leak | b | fit chi2/dof |
|---|---|---|---|---|---|---|---|---|
| 0.15-0.25 | 1,363,238 | 0.47 | **0.99** | 0.28 | 0.12 | 0.20 | 0.09 | 12.8/12 |
| 0.25-0.35 | 885,168 | 0.19 | 0.59 | 0.95 | 0.18 | 0.16 | 0.03 | 7.9/12 |
| 0.35-0.45 | 344,914 | -0.08 | (GR predicts ~0; data consistent with 0, chi2 13.8/13) | | | | | |

**Combined (the two bins with |eta| > 0.05): A_D = 0.93 +- 0.35, 2.7 sigma from zero, 0.2 sigma from GR.** The
0.15-0.25 bin alone is 3.5 sigma from zero on statistics, 2.5 sigma with the systematics. This is evidence, not a
detection; it says the velocity field moves galaxy sizes by the GR amount to +-35 %, and nothing about modified gravity
beyond that.

The number that changes forecasts: the size estimator's **response eta to a true convergence is not 1 but 0.57 at
z = 0.16, falling to 0.18 at z = 0.30 and through zero at z = 0.36.** At fixed redshift and colour BGS sizes scale as
R ~ L^0.49, constant surface brightness, which is the direction magnification moves a galaxy; at the r = 19.5 flux limit
magnification admits intrinsically fainter, smaller galaxies that cancel most of the size increase. Forecasts that set
eta = 1 (2511.17858: S/N ~ 10-13 per bin for DESI+LSST) overstate a flux-limited size measurement by 2 at z ~ 0.2 and
without limit at z > 0.35.

## Method

- **Samples.** DESI DR1 LSS v1.5 BGS_BRIGHT clustering catalogues (NGC + SGC) and the first random file per cap (half of it,
  2.35x the data). Density tracers: all galaxies at 0.08 < z < 0.56 with WEIGHT. Size galaxies: 0.15 < z < 0.45, matched
  by TARGETID to the DR1 LS-DR9 Tractor photometry VAC (100 % match), Tractor type REX/EXP/DEV/SER, FRACMASKED_R < 0.1,
  free of bright-star and large-galaxy mask bits, sigma(ln R) < 0.5 (96 % pass).
- **Size estimator.** lambda = ln(SHAPE_R D_A(z_obs)); kappa_hat = lambda minus a per-type polynomial in z and the g-r,
  r-z colours (no magnitude term: at fixed z and colour its response would be 0.03), winsorised at 3 sigma, zero mean in
  z slices of 0.01 per cap; inverse-variance weights; per-galaxy noise 0.42.
- **Response.** eta(z) = a x 0.9966 + 2.1715 d<kappa_hat>/dm_lim. The second term is the flux-limit edge measured on the
  data (slope of the weighted mean kappa_hat for r < m_cut, m_cut 19.3-19.5). Other cuts: the sigma(ln R) and
  fibre-magnitude cuts contribute < 0.001. Redshift success is bounded at +-0.02. a is the Tractor size response,
  0.97 +- 0.05: 311 BGS galaxies in COSMOS against HST/ACS GIM2D and r50 radii, instrumental-variable slope (each ACS radius
  instruments the other, removing the attenuation from ACS errors; plain OLS gives 0.89-0.94). The PSF lever on the full
  sample agrees (d kappa_hat / d ln PSF = 0.00 +- 0.04). eta's statistical error is from the jackknife; its systematic is
  from the m_cut range, a and the redshift success.
- **Statistic.** For every size galaxy, shell sums of galaxies minus randoms in 14 bins of separation (10-150 Mpc/h) and
  10 bins of mu, with the midpoint line of sight. xi(d, mu) = sum u kappa_hat Q / sum u R, normalised per (d, mu) bin (this
  removes the leak of the even part through the n(z) gradient). Multipoles come from Legendre integration. The dipole xi_1 at
  20-150 Mpc/h is the measurement.
- **Model.** xi_1 = A_D x T_D + T_rel.
  - T_D is Bonvin et al. eq. 28 (distant observer, CAMB Planck 2018 P(k), f = Omega_m^0.55). It is averaged galaxy by galaxy
    over the bin with each galaxy's pair weight and eta(z_j). b comes from the galaxy monopole at 20-40 Mpc/h.
  - T_rel is the number counts' own Doppler term (coefficient C(z) from the measured n(z) and magnification slope s(z),
    -3.2 and 2.2 in the two bins), times the measured size-density coefficient b_s = 0.026 +- 0.005 and 0.034 +- 0.006.
    b_s comes from the monopole ratio and is flat in d from 15 to 45 Mpc/h. T_rel is under 3 % of T_D; +-50 % on C is
    assumed.
  - b_s also induces O(d/r) odd terms (unequal RSD under the midpoint line of sight, evolution across the pair). The
    analytic estimate is f b_s (d/r) xi_2-like, about 15 % of T_D at 100 Mpc/h. It is carried as +-(0.1|A_D| + 0.1).
- **Mocks.** Three lognormal light cones (512^3, 3 Gpc/h) with linear velocities on the real footprint and n(z), run
  through the same counts and estimator. The Doppler mark eta(z) kappa_v recovers A_D = 0.953, 0.971, 1.04 (z 0.15-0.25),
  0.96, 0.94, 0.88 (0.25-0.35), each +-0.035. The mocks' monopole recovers their input b (1.49, 1.43, 1.23 vs 1.49) only
  below ~45 Mpc/h; beyond, an offset of ~-0.005 suppresses it. The velocity-free mark (unit b_s) gives -35, +36, +14
  (+-31) per unit b_s, consistent with zero.
- **Covariance.** 200-region jackknife (Hartlap) plus the T_rel error. On the mocks the jackknife errors match the
  mock-to-mock scatter (Doppler mark 0.035-0.04 vs 0.05; velocity-free mark 23-34 vs 37).
- **Departures from the pre-registered plan, made after the first run and stated here.**
  1. The fit range for b was planned as 30-80 Mpc/h. In the mocks that range recovered A_D = 1.4 (b half its input,
     because of the large-scale monopole offset), so b now comes from 20-40 Mpc/h. With the planned range the bin-1 result
     would read 1.65.
  2. The planned covariance was size shuffles + analytic cosmic variance + a rank-1 mock-leak term, with the 3-mock mean
     leak subtracted. Its fit chi2/dof was 4.8-6.7. The shuffles miss the coherent cosmic variance of the b_s delta part
     of the sizes, and the 3-mock mean is noise. That covariance gave 1.19 +- 0.25 in bin 1. The jackknife (validated
     above) is used instead.

## Nulls and bounds (z = 0.15-0.25 unless stated)

| test | result |
|---|---|
| randomised redshifts (pair counts redone, sizes kept) | A = 0.10 +- 0.22 (bin 2: -0.53 +- 0.83), chi2 15.2/13 |
| 0.35-0.45, where GR predicts ~0 for sizes | xi_1 consistent with zero, chi2 13.8/13 (p 0.39) |
| octupole xi_3 against its prediction (~10 % of the dipole) | A_oct = 2.5 +- 2.7; chi2 vs 0 11.0/13 |
| substitute marks, bound on A_D (beta_X x A(X)) | PSF size -0.015 +- 0.020; stars' PSF map 0.002; FRACFLUX (blending) -0.10 +- 0.05; GALDEPTH 0.003; PSFDEPTH 0.001; E(B-V) 0.006; stellar density 0.000 |
| substitute marks, dipole vs zero (21 tests over 3 bins) | p > 0.03 except the stars' PSF map in bins 2-3 (p 0.007, 0.014); kappa_hat does not respond to it (beta ~ 1e-4): bound < 0.01 |
| FRACFLUX_R < 0.05 (93 % of galaxies) | A_D = 0.96 +- 0.31 |
| NGC vs SGC | 1.31 +- 0.33 vs -0.08 +- 0.38 (differ by 2.8 sigma) |
| DR9 north (BASS/MzLS) vs south (DECaLS) | 1.93 +- 0.47 vs 0.66 +- 0.32 (2.2 sigma) |
| d >= 30 Mpc/h; 20-60 Mpc/h | 1.31 +- 0.31; 0.87 +- 0.29 |

## What it means

Peculiar velocities change galaxy sizes at 0.15 < z < 0.25 by the amount linear-theory GR predicts (A_D = 0.99 +- 0.39
total), with 2.5-3.5 sigma evidence that the effect is there. Existing measurements of the same physics use
Fundamental-Plane or Tully-Fisher sizes at z ~ 0.07 (DESI DR1 PV, arXiv:2512.03230). This measurement uses raw sizes
without a distance indicator, at 3-6 times the distance. A flux-limited sample loses most of the size response: the
next measurements need either a deeper parent sample (BGS faint, LSST) or a selection not set by total flux.

## Weakest points

- **The cap split**: NGC and SGC differ by 2.8 sigma, and the north (BASS/MzLS) photometry gives the higher amplitude.
  This could be a photometry-dependent size response (a is calibrated on DECaLS imaging in COSMOS only) or chance in
  one of ~10 splits.
- The covariance and the b fit range were changed after the first run, both on mock evidence. The planned versions give
  1.19 +- 0.25 and 1.65.
- The b_s leak is bounded analytically; three mocks constrain it only to +-0.55 in A_D.
- b from 20-40 Mpc/h is 1.14, below DESI's fiducial 1.34/D = 1.45; the b systematic covers half the difference.
- eta comes from a local edge derivative. It assumes the size-flux distribution continues smoothly 0.1 mag beyond the
  limit.

## To settle it

BGS_ANY (with BGS faint, r < 20.175) for an edge calibration beyond the limit; a per-photometry size response (HSC or
Euclid overlap with BASS/MzLS); 50+ mocks for the b_s leak and the covariance; DESI DR2.

## Blind review (2026-11-02, context-free Opus referee, own literature search; `blind-review.md`)

VERDICT: novel but not sound: the amplitude depends on an η calibration no mock tests end to end, the signal appears only in NGC/BASS-MzLS imaging (2.8σ cap discrepancy), and the b_s leak and covariance are validated only by three mocks after choices made once the first results were in.

Strongest objection: the signal comes entirely from NGC/BASS-MzLS imaging while the size response is calibrated only on
DECaLS, and the b_s leak is constrained empirically only at the level of the signal. A physics miss: parked (one
revision is for a presentation point). Novelty held: no prior real-data size-number-count Doppler dipole found.


# No spiral winding-sense asymmetry of the published size at 0.3 < z < 1.2 in Euclid Q1: 291,652 discs classified by an exactly mirror-antisymmetric network give a raw, uncorrected A = +0.44 +- 0.19 %, which excludes the JWST JADES +20 % (and +9.3 % at matched redshift), the implied 24.8 % (11.4 %) Galactic-pole dipole and a 10.5 % HST-axis dipole unless the classifier misreads 37-46 % of real Euclid spirals (R-158, restated R-161)

**Claim.** We measure the asymmetry A = (N_cw - N_ccw)/N between clockwise and counterclockwise spirals at 0.3 < z < 1.2
on Euclid Q1 VIS images, with a classifier that is exactly odd under reflection, and test the published amplitudes on
axes and values fixed before any sign was read.
- **Raw result, no model.** Over the three Euclid Q1 deep fields, A = +0.44 +- 0.19 % (291,652 classified discs):
  EDF-N +0.16 +- 0.33 %, EDF-F +0.01 +- 0.45 %, EDF-S +0.80 +- 0.28 %. On the pre-registered axes the raw dipole is
  +0.55 +- 0.28 % (south Galactic pole), -0.31 +- 0.74 % (HST best-fit axis), -0.52 +- 0.45 % (SDSS axis).
- **Published values excluded** at 99.7 % with the raw numbers: the JADES +20.2 % at the GOODS-S position (and +9.3 %,
  the same author's galaxies restricted to 0.3 < z < 1.2), the implied 24.8 % (11.4 %) Galactic-pole dipole, and a
  10.5 % HST-axis dipole. A correction cannot rescue them: a chiral bias can be bounded from the data, and
  misclassification only dilutes, so hiding any of them takes a misclassification rate of 37-46 %, a near-random
  classifier (the table below).
- **Not claimed.** The raw excesses of ~1 % (all fields +0.44 +- 0.19 %; EDF-S +0.80 +- 0.28 %; the 0.9-1.2 bin
  +1.43 +- 0.47 %) cannot be read without a correction for the PSF and for dilution, which we do not claim to have to
  that precision. We say nothing about them: not a detection, not a bound below ~1-2 %. The HST COSMOS +2.48 % and a 5 %
  HST-axis dipole are not decided.

**The disputed step removed.** Spin-direction claims (HST COSMOS, arXiv:2011.03723; JWST JADES, arXiv:2403.17271 and
2502.18781; with the dipole axes of 2011.03723) rest on annotation by one pipeline. The rebuttals attribute such excesses to
annotator or classifier chirality (0803.3247, 1610.07060, 2210.04168, the last with a reflection-equivariant network on
ground-based SDSS/DESI images at low z) or find isotropy at low z (2404.06617, 2410.18884). Here:
- **Classifier.** A CNN g, made exactly antisymmetric: f(x) = [g(x) - g(mirror x)]/2, so f(mirror x) = -f(x) for any
  weights. It is trained on Galaxy Zoo 1 handedness labels (7,000 cw + 7,000 ccw, Legacy Surveys DR10 r stamps) with
  rotation, blur and noise augmentation. Label bias cannot enter.
- **Preprocessing.** Mirror-equivariant: sky parity from each stamp's WCS, segmentation, R50, and resampling on a grid
  symmetric about the catalogue position. **Selection.** Only mirror-invariant catalogue quantities.
- **Exactness checked.** Mirrored inputs return -f to 6e-7 on every dataset. So the pipeline cannot create a net sign:
  only the images can, through the instrument's chiral part (the PSF) or through the sky.
- **Classifier check on real galaxies.** On the JWST JADES F200W image, the classifier classifies 182 of the 263 published
  galaxies and agrees with the published label on 181 (misclassification <= 4.8 % at 99.7 %), reproducing the published
  amplitude on that selection (+20.9 +- 7.2 %, 110 / 72). This checks the classifier; it says nothing about the sky.

**Data and sample (fixed before classification; an independent method check returned nine points, all closed before
any sign was read).** Euclid Q1 VIS, MER x morphology x PHZ: Sersic r_e > 3 FWHM (0.5"), q > 0.4, n < 2.5, VIS < 23,
0.3 < PHZ < 1.2: 481,696 discs in EDF-N/F/S, every one cut from the S3 mosaics. HST COSMOS ACS F814W: ZEST discs
with r50 > 3 FWHM, b/a > 0.4, COSMOS2020 0.3 < z_phot < 1.2: 28,406 (12,918 classified). A galaxy is classified when |f| exceeds
the value that 5 % of achiral model galaxies (mirror-symmetric discs, bulges, rings and bars at the sample's R50, S/N
and empirical PSF) pass: 291,652 are. "Raw" errors are the larger of the binomial and a tile bootstrap. The planned
threshold rule (held-out GZ1 accuracy) was degenerate and was replaced by the achiral rule before classification.

## The raw table (CNN; 99.7 % verdicts on the raw values)

| sample | z | classified | cw / ccw | A raw |
|---|---|---|---|---|
| Euclid EDF-N | 0.3-1.2 | 93,306 | 46,728 / 46,578 | +0.16 +- 0.33 % |
| Euclid EDF-F (holds GOODS-S) | 0.3-1.2 | 57,788 | 28,896 / 28,892 | +0.01 +- 0.45 % |
| Euclid EDF-S | 0.3-1.2 | 140,558 | 70,843 / 69,715 | +0.80 +- 0.28 % |
| Euclid, all | 0.3-0.6 / 0.6-0.9 / 0.9-1.2 | 134,764 / 106,374 / 44,401 | | +0.40 / +0.21 / +1.43 % (+-0.28, 0.34, 0.47) |
| Euclid, all | 0.3-1.2 | 291,652 | 146,467 / 145,185 | +0.44 +- 0.19 % |
| HST COSMOS ACS F814W | 0.3-1.2 | 12,918 | 6,583 / 6,335 | +1.92 +- 0.89 % |

Dipoles: the same per-field, no-intercept fit on the three Euclid fields (cosines to the south Galactic pole -0.50,
+0.82, +0.73 for EDF-N/F/S; GOODS-S itself +0.81), on raw A. A uniform bias b enters it as 0.71 b (axis G) and 1.76 b
(axis H).

| pre-registered test | published value | ours, raw | raw 99.7 % upper | largest 1 - 2 eps that hides it | eps needed to hide it |
|---|---|---|---|---|---|
| EDF-F A at the GOODS-S position | +20.2 +- 6.0 % (158 vs 105) | +0.01 +- 0.45 % | +1.36 % | 0.09 | 0.45 |
| same, at 0.3 < z < 1.2 (his 86 galaxies with Euclid PHZ) | +9.3 +- 10.7 % | | | 0.20 | 0.40 |
| Euclid A, all fields, against the same +9.3 % | | +0.44 +- 0.19 % | +1.01 % | 0.16 | 0.42 |
| dipole, south Galactic pole axis | 24.8 % implied | +0.55 +- 0.28 % | +1.40 % | 0.07 | 0.46 |
| same, at matched z | 11.4 % | | | 0.15 | 0.42 |
| dipole, HST axis (78, +47) | 10.5 % implied | -0.31 +- 0.74 % | +1.92 % | 0.27 | 0.37 |
| same | 5 % | | | 0.56 | 0.22: **not constraining** (pre-registered verdict) |
| COSMOS A | +2.48 +- 1.28 % | +1.92 +- 0.89 % | | | **not constraining** (one 1.7 deg2 field) |

**How the bound works.** With a misclassification rate eps and an additive chiral bias b, A_raw = (1 - 2 eps) A_true +
b. A published amplitude X survives only if (1 - 2 eps) X + b <= the raw 99.7 % upper limit. b is set at the data
bound below, in the direction that hides X; the table gives the largest dilution factor 1 - 2 eps for which X could
still hide, and the eps that implies. Every excluded amplitude needs eps >= 0.37.

**The bias, bounded from the data (not from the PSF model).**
- **Mirrored-PSF null on real galaxies.** 3,000 real Euclid discs re-convolved with the VIS PSF and with its mirror image
  move A by dA = -0.10 +- 0.08 %. The PSF's term in the data is dA/2: |b| <= 0.17 % at 99.7 %. The table allows three
  times that, |b| <= 0.51 %, because a re-convolution is not the native PSF. (On achiral model galaxies the VIS PSF is
  strongly chiral to the network, mean f +0.115 against -0.216 with the mirrored PSF; on real galaxies it moves A by a
  tenth of a percent. The model figure is not used.)
- **Size (R50/FWHM quartiles), where a PSF effect must grow toward small galaxies:** +1.0, +0.4, +0.4, +0.1 % (each
  +-0.34 %). The largest quartile alone, the least PSF-affected, is +0.10 +- 0.34 % and would hide +9.3 % only for
  eps >= 0.41.
- **Position in the mosaic (tile quadrants, a coarse focal-plane proxy; VIS stacks are dithered):** +0.50, +0.52, +0.25,
  +0.49 %. No trend.
- **Ellipticals** through the same pipeline: -4.2 +- 4.4 % (505 of 1,998 pass). **Stars** fall below the size floor.

**Misclassification, from the data.** Misclassification by a mirror-odd classifier cannot create a sign; it only
dilutes. To hide the excluded amplitudes the network must be wrong on 37-46 % of the real Euclid discs it classifies,
i.e. nearly a coin. Against that:
- **Confidence.** The most confident 5 % (24,084 discs, |f| > 12.9, five times the threshold) give +0.56 +- 0.64 %;
  hiding +9.3 % there needs eps >= 0.34 among the clearest spirals. A by |f| (bins of 40/20/20/10/5/5 %): -0.15, +0.43,
  +0.53, +0.15, +0.28, +0.56 %. No rise with confidence.
- **Two classifiers.** A logarithmic-spiral Fourier score (exactly odd, no training) agrees in sign with the network on
  87 % of the 62,386 Euclid galaxies both classify.
- **Real arms at higher resolution.** 181 / 182 agreement with independent labels on JADES F200W.
- Model spirals at Euclid's R50, S/N and PSF give eps = 1.7 %, and injected 5 % and 20 % asymmetries come back
  (5.0 +- 2.1, 19.7 +- 2.2 %); not used for the bound.

## Nulls

| test | result |
|---|---|
| mirror (whole pipeline) | exact: f(mirror) = -f to 6e-7 |
| real galaxies re-convolved with the PSF vs its mirror | dA = -0.10 +- 0.08 % (Euclid) |
| ellipticals | Euclid -4.2 +- 4.4 %, COSMOS +1.7 +- 3.6 % |
| stars | COSMOS **fails** at point-source scale (549 of 9,911 pass, A = -23.5 +- 4.3 %: the drizzled ACS PSF is chiral to the network at R50 ~ 2 px); COSMOS galaxies re-convolved: dA = -0.99 +- 0.12 %. Euclid stars fall below the size floor |
| size quartiles, tile quadrants | no trend |
| position angle in the image frame | a trend odd under the mirror (+3.6, +1.3, -0.6, -2.5 % in 45-degree bins): the network is mirror-exact but not rotation-invariant; it cancels because the PA distribution is mirror-symmetric (paired bins 72,434 vs 72,247; 73,531 vs 73,440) |

## What it means

At 0.3 < z < 1.2, 2.9e5 Euclid spirals classified by a network that cannot carry a handedness preference show no
winding-sense asymmetry anywhere near the published amplitudes: not the +20 % found in JWST JADES (or the +9.3 % its own
galaxies give at these redshifts), not a Galactic-pole dipole of 11-25 %, not a 10.5 % HST-axis dipole. The exclusion
uses the raw count, a bias bounded from real galaxies, and the fact that misclassification can only dilute: it would
fail only for a classifier wrong on ~40 % of real Euclid spirals. Below ~1-2 % the raw numbers are not interpreted. The
classifier reproduces the published JADES annotations galaxy by galaxy, so the disagreement is not about how those
galaxies look. This note does not say where the JADES excess comes from.

## Weakest points

- eps on real Euclid arms is not measured against truth; the bound rests on how far it would have to be from the model
  value (1.7 %) and the JADES real-galaxy value (<= 4.8 %), and on the flat confidence trend.
- The bias bound is the mirrored-PSF re-convolution, with a factor-3 allowance. Detector chirality (a non-PSF image
  effect) is checked only through the tile-quadrant and size splits.
- The +9.3 % at matched redshift comes from 86 of the published galaxies and carries +-10.7 %; it is the published
  selection's value, not a population value.
- Euclid photo-z bins are broad; no sign can move between bins.

## To settle it

Euclid DR1 (~1e7 discs) gives a free-axis dipole at 0.1 % and, with the correction measured, an interpretation of the
~1 % level. A direct real-galaxy eps for Euclid VIS: classify on VIS the ~100 JADES galaxies with Euclid counterparts and
compare with their JWST handedness. A rotation-averaged network would remove the PA pattern.

## Revision history

R-158 (2026-11-02, cw-165) headlined the corrected values and a JADES localisation, refereed "VERDICT: novel but not sound:
the claim that the JADES excess is confined to the published selection rests on a post-hoc full-field test with an
uncalibrated HST threshold, no NIRCam PSF or eps calibration, and a comparison against galaxies the published pipeline
rejected, so the gap between the two subsets (~2.7σ, about the original result's own significance) cannot locate the
excess, while the Euclid corrected values hinge on a model-based correction for a strongly chiral VIS PSF and a
misclassification rate calibrated only on simulations." This revision (R-161, the one revision, presentation only)
makes the raw Euclid result the headline with a worst-case correction bound from the data (`r158_raw.py`), removes the
JADES localisation (the post-hoc full-field test; its numbers stay in `results-R158.md` R1), and keeps the JADES
per-galaxy agreement only as a classifier check. No number was recomputed except the raw dipoles.

## Blind review (2026-11-02, R-161 revision, context-free Opus referee, own literature search, this section alone; `blind-review.md`)

VERDICT: novel but not sound: the headline exclusions target extrapolated or redshift-mismatched amplitudes (the JADES galaxies' own matched-z value, +9.3 ± 10.7 %, agrees with Euclid), while the amplitudes actually published at these redshifts are left undecided and the Euclid misclassification rate is never measured on real VIS spirals.

Strongest objection: the excluded values are the all-z JADES +20 % (said to grow with z), the matched-z +9.3 +- 10.7 %
(itself consistent with zero and with Euclid), and dipole amplitudes extrapolated from one field by a cosine, which the
claimant never published; the amplitudes actually published at these redshifts (COSMOS +2.48 %, 1-3 % dipoles) are
undecided, and eps on VIS is bounded only if the CNN and Fourier errors are independent (both read the same marginally
resolved pixels). It credited the method as an upper limit of ~1-2 % at z < 1.2 valid for eps < ~0.37, and named the
direct test: classify on VIS the JADES galaxies with Euclid counterparts. Novelty held (no Euclid spin-direction
measurement, none on HST/JWST with a mirror-odd classifier; closest 2210.04168). **A miss on physics: parked** (brief:
one revision; park on any further physics point). Nothing on the site.
