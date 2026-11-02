# Results (R-158, cosmo-deep-25, 2026-11-02)

Is the spiral winding-sense (spin-direction) asymmetry reported on space-based images at 0.3 < z < 1.2 real? A
bias-controlled test with an exactly mirror-antisymmetric classifier on HST COSMOS, Euclid Q1 VIS and JWST JADES images.
Trajectory check (2026-11-02, in the brief): "partly novel": what is open is an independent, bias-controlled,
pre-registered test of the redshift-growing, axis-aligned asymmetry on a large resolved space-based sample at 0.3<z<1.2.
This file is the worker's own (cloud worker cw-165).

Scripts (all `research/cosmo-discovery/`): `r158_common.py` (IRSA TAP with retries), `r158_counts.py` (the first
all-field async join; abandoned, see F0), `r158_euclid_cat.py` (Euclid catalogue per MER tile), `r158_cosmos_cat.py`
(ZEST x COSMOS2020), `r158_gz1_fetch.py` (training stamps), `r158_cosmos_stamps.py`, `r158_euclid_stamps.py`,
`r158_jades_fetch.py` + `r158_jades_stamps.py` (Shamir's 263 JADES galaxies, positions in `r158_shamir_jades.csv`),
`r158_classify.py` (prep, the antisymmetric CNN, the log-spiral Fourier classifier), `r158_train.py`, `r158_prepall.py`.
Plan: `r158_plan.md`. Data under `data/r158/` (gitignored).

## F0. Access (fetch times; `data/` started empty)

| what | where | size | time |
|---|---|---|---|
| Shamir and rebuttal papers (2011.03723, 2206.00168, 2210.04168, 2403.17271, 2404.06617, 2410.18884, 2502.18781) | arxiv.org/pdf | 48 MB | seconds |
| Galaxy Zoo 1 table 2 (Lintott+2011) | galaxy-zoo-1.s3.amazonaws.com | 20 MB zip | seconds |
| 14,000 GZ1 r-band stamps, 96 px at 0.262" (13,892 valid) | legacysurvey.org cutout service, DR10, 48 threads | 0.4 GB | 10.7 min |
| ZEST morphology (131,532) and COSMOS2020 Classic photo-z at 0.2<z<1.4 (353,621) | IRSA TAP async | | 1.4 min |
| COSMOS ACS F814W tiles (acs_2.0, 569 of 577 hold targets, ~200 MB each), streamed and deleted after cutting | irsa.ipac.caltech.edu/data/COSMOS/images/acs_2.0/I | ~110 GB read | see F2 |
| Euclid Q1 catalogue, MER x morphology x PHZ: one async join over all fields | IRSA TAP async | | killed after 45 min, never finished |
| Euclid Q1 catalogue per MER tile (352 tiles), sync | IRSA TAP | | ~20 s per tile, 16 threads |
| JADES GOODS-S DR2 F200W SCI extension (25000 x 27700 float32) | mast.stsci.edu Download API, 32 byte ranges, 12 threads | 2.77 GB | 51 s |

archive.stsci.edu and skyserver.sdss.org fail TLS through the session proxy (unknown issuer / self-signed); MAST's
download API and the Legacy Surveys viewer work.

## F1. Step 0: feasibility (the gate)

**Shamir's anchors** (fixed from the papers before any count):
- HST COSMOS (2011.03723 Table 2): 3,116 cw vs 2,965 ccw: A = +2.48 % +- 1.28 %.
- JWST JADES GOODS-S (2502.18781): 158 cw vs 105 ccw: A = +20.2 % +- 6.0 %; the author ties it to an axis near the
  Galactic pole (cw excess toward the SGP) growing with redshift. Coordinates of all 263 are in his Tables 1-2.
- Dipole axes (2011.03723): HST best fit (78, +47); SDSS z > 0.15 best fit (71, +61). Implied dipole amplitudes at our
  fields: axis H 10.5 % (P1 / cos 0.237); axis G (SGP) 24.8 % (P2 / cos 0.813); the proposal's ~5 %.

**Counts before classification:**

| dataset | cuts | N | 0.3-0.6 | 0.6-0.9 | 0.9-1.2 |
|---|---|---|---|---|---|
| COSMOS ACS F814W (1.7 deg2) | ZEST type 2, r50 > 0.285" (3 FWHM), a/b < 2.5, I < 24, COSMOS2020 0.3<z<1.2 | 28,710 | 9,955 | 10,800 | 7,955 |
| Euclid Q1 VIS (~63 deg2; 14 tiles sampled at gate time) | VIS < 23, Sersic r_e > 0.5" (3 FWHM), q > 0.4, n < 2.5, PHZ 0.3<z<1.2 | ~1,000 per MER tile, ~2.5e5 in Q1 | 469/tile | 339/tile | 192/tile |
| Euclid, VIS < 22.5 | same | ~670 per tile | 355/tile | 223/tile | 91/tile |
| JADES F200W | Shamir's own 263 | 263 | (his z: mostly > 1) | | |

Euclid cutouts run at ~40 per second through IRSA's IBE cutout service, so the Euclid sample is limited by fetch time,
not by the catalogue: a random-order subsample (mirror-invariant) is fetched for a fixed 75 minutes, ~1e5 galaxies.

**Forecast.** sigma_A = 1/sqrt(N_c) / (1 - 2 eps), N_c the classified number (a fraction f_c of the sample) and eps the
misclassification rate. With f_c = 0.25 and eps = 0.1 (GZ1 held-out accuracy at the chosen threshold is >= 0.9):
- COSMOS: N_c ~ 7,200, sigma_A ~ 1.5 %. P1 (+2.48 %) vs 0 separates at ~1.7 sigma only: **P1 itself cannot be confirmed
  or excluded at 3 sigma** by any sample in one field of that size (Shamir's own was 1.9 sigma).
- Euclid (1e5 fetched, three fields 36/22/42 %): sigma_A per field 1.2-1.7 %; EDF-F (which contains GOODS-S): P2's
  +20 % vs 0 at > 10 sigma overall and > 4 sigma in the 0.9-1.2 bin alone.
- Dipole along axis G (four fields, cos phi -0.67, -0.50, +0.82, +0.73): sigma(A_d) ~ 1.1 %: the implied 24.8 % at
  > 20 sigma. Axis H (cos 0.24, 0.40, 0.20, -0.11): sigma(A_d) ~ 3 %: 10.5 % at ~3.5 sigma; 5 % not (1.7 sigma).
- JADES reclassification of Shamir's 263: A to +-6 % on the same galaxies, plus per-galaxy agreement.

**Gate: passes** for P2, the axis-G dipole and the axis-H dipole at its implied 10.5 %; it fails for P1 (the COSMOS
+2.5 % alone) and for a 5 % axis-H dipole, which are reported as intervals without a 3 sigma verdict. Even at f_c = 0.05
the axis-G and P2 tests stay above 3 sigma.

## M1. Method check (blind-review.md 1b; fresh context-free Opus agent, plan alone, prompt verbatim; 2026-11-02)

Verdict, verbatim:

"METHOD: fixable: calibrate ε on the target population (blind mirror-randomised labels or injected chiral models at the sample's R50/S/N/PSF), measure the armless-disc fraction, and fix a floor on (1−2ε) below which a bin is "not constraining"; redo the injection test by true handedness, not classifier output; evaluate P2 and axis G at 0.3<z<1.2 from Shamir's own galaxies in that range (or a pre-fixed growth law); check the P1/axis-H sign convention on P1's own galaxies and confirm training-stamp parity; correct the Bonferroni bar to ~3.69σ and define how exclusions are counted across tests; add a "not constraining" verdict with a pre-computed power estimate; state the dipole fit's intercept and control the HST-vs-Euclid instrument offset; subtract a measured PSF chirality as a signed bias instead of adding it in quadrature"

Strongest point (its words): the dilution correction "is only as good as ε", and ε from bright low-z GZ1 spirals does not
describe barely resolved arms at z ~ 1. All nine closures are in `r158_plan.md` ("Method-check closures"), fixed before
any target galaxy was classified. Points closed by a measurement (eps, power, P2 at matched z) are reported in M2.

## M2. Closures measured, and the power estimate (sign-free; written before any sign was computed)

- **Thresholds** (revised rule, `r158_plan.md`): t = |score| passed by 5 % of achiral simulations. COSMOS (3,000 each
  achiral, mirror-PSF, chiral; PSF = median of 1,730 ACS stars): t_A = 2.89, t_B = 0.42.
- **Chiral injections at COSMOS R50/S/N/PSF** (`r158_thresholds.py`): CNN passes 81 % of injected spirals with 1.1 %
  sign error; Fourier passes 22 % with 1.9 % error. Injected A = 5.0 % -> recovered 6.6 +- 2.1 % (CNN), 8.9 +- 4.2 %
  (Fourier); injected 20.0 % -> 20.8 +- 2.3 % and 22.4 +- 4.5 %. The dilution model recovers both.
- **PSF term on achiral models (COSMOS):** mean f = -0.012 +- 0.024 (mirrored PSF +0.019); A among passing achiral
  models -5 +- 8 %; Fourier mean S = -0.006 +- 0.003 (mirrored +0.002). No PSF chirality detected for the CNN; the
  Fourier score carries a ~2 sigma mean offset, small against its threshold (0.42).
- **GZ1 held-out (clean):** CNN 97.6 % correct on all, Fourier 76 %; at COSMOS S/N (median 137) CNN 90.5 %, at JADES
  S/N 97.2 %.
- **P2 at matched redshift:** 113 of Shamir's 263 JADES galaxies match a Euclid Q1 object (VIS < 23, r_e > 0.4") within
  1"; 86 have PHZ 0.3 < z < 1.2: 47 cw vs 39 ccw by his labels, A = +9.3 % +- 10.7 %. That is P2's value for test (ii)
  at 0.3 < z < 1.2 (>= 30 matched, so the fallback growth law is not used); axis G's implied A_d at matched z =
  +9.3 / 0.813 = 11.4 %.
- **Power, COSMOS** (counts of |f| > t only): 28,406 discs prepped, 12,918 classified (45 %); armless bound
  w = 0.05 x 28,406 / 12,918 = 0.11; eps = 0.064; 1 - 2 eps = 0.87 (above the 0.3 floor in every bin). Expected
  sigma(A_corr) = 1.0 % overall; 1.7 / 1.6 / 2.0 % in the three bins. So test (i) (P1 = +2.48 %) separates from zero at
  ~2.5 sigma at best: it is "not constraining" unless A_COSMOS comes out below ~-0.5 % or above ~5.5 %. Euclid's power
  is written in M3 before its classification.

## M3. Euclid closures and power (sign-free; written before the Euclid signs were read)

- **Euclid stamps** (F0 addendum): IRSA's IBE cutout service reset connections under load (~10/s), so stamps were cut
  from the MER BGSUB VIS mosaics on IRSA's S3 bucket by HTTP byte-range row stripes (`r158_euclid_s3cut.py`): all 320
  tiles holding targets, 498,003 stamps cut and prepped in 27.8 min (5 processes x 4 streams). The catalogue came from
  the merged-objects HATS parquet on S3 (114 files, column chunks only, 7.5 min; `r158_euclid_cat_s3.py`); a
  pyarrow.dataset filter scan and the per-tile TAP queries were both abandoned as too slow. Sample after cuts and
  dedup: 495,975 discs (EDF-N 155k, EDF-F 95k, EDF-S 231k at 0.3<z<1.2), 1,998 ellipticals, 800 stars (only 30
  pass prep: Euclid stars are below the 1.5 px R50 floor, so they mostly never enter the pipeline).
- **Thresholds** (achiral 5 %, 2,841 sims, PSF = median of 793 VIS stars): t_A = 2.41, t_B = 0.32. Chiral injections:
  CNN passes 87 % with 1.7 % error; injected 5.0 % -> 5.0 +- 2.1 %, 20.0 % -> 19.7 +- 2.2 %. Fourier passes 30 % with
  9.9 % error, and its injections do not recover (5 % -> -2.3 +- 4.3 %, 20 % -> 12.5 +- 4.6 %): **classifier B is not
  usable on Euclid** and is reported only for agreement.
- **PSF chirality, Euclid (null 4):** mean f on achiral models +0.115 +- 0.024 with the VIS PSF and -0.216 +- 0.024 with
  its mirror: the classifier sees the VIS PSF as chiral at ~10 sigma at galaxy scale. A among passing achiral models
  -7 +- 8 % (PSF), -21 +- 8 % (mirror). Per closure 8 this is a signed bias to subtract; its size on real galaxies comes
  from null 4b. COSMOS with 1,730 individual star kernels instead of the median: mean f +0.027 +- 0.029 vs mirror
  +0.088 +- 0.029, A among passing +17 +- 10 % vs +5 +- 10 %: no PSF chirality detected at galaxy scale in ACS.
- **COSMOS stars null (null 2) fails at point-source scale:** of 9,911 ACS point sources 549 pass t_A (5.5 %), with
  A = -23.5 +- 4.3 % and mean f over all stars -0.199 +- 0.014. The drizzled ACS PSF is chiral to the classifier at
  R50 ~ 2 px. The galaxy sample's A shows no size trend (quartiles of R50/FWHM: +0.3, +1.7, +2.1, +2.7 %, each +-1.8 %).
- **Power, Euclid CNN** (counts of |f| > t only; 1 - 2 eps = 0.87-0.89 everywhere): sigma(A_corr) = 0.37 / 0.47 /
  0.30 % for EDF-N / EDF-F / EDF-S over 0.3-1.2 (0.4-1.3 % per bin), before the PSF term. Tests (ii) and (iii) are
  decisive if the PSF term is controlled to that level; (iv) at 10.5 % is too.

## R1. The table (CNN, the primary classifier; Fourier as the second)

A = (N_cw - N_ccw)/N over galaxies with |f| > t. "Raw" carries the larger of the binomial and tile-bootstrap errors.
"Corrected" = (A - PSF term)/(1 - 2 eps). The PSF term is signed: the armless share w (conservative bound) times the
achiral-model sign bias, plus null 4b's dA/2. Its error is propagated. eps = 0.055-0.068 everywhere, so 1 - 2 eps = 0.86-0.89.
`data/r158/stats_rows.csv` and `stats.json` hold the full set (`r158_stats.py`).

| dataset | z | discs | classified | cw / ccw | A raw | A corrected |
|---|---|---|---|---|---|---|
| COSMOS (HST ACS F814W) | 0.3-0.6 | 9,827 | 4,584 | 2,345 / 2,239 | +2.31 +- 1.50 % | +3.9 +- 2.0 % |
| | 0.6-0.9 | 10,701 | 4,985 | 2,524 / 2,461 | +1.26 +- 1.47 % | +2.7 +- 2.0 % |
| | 0.9-1.2 | 7,869 | 3,345 | 1,711 / 1,634 | +2.30 +- 1.73 % | +4.0 +- 2.3 % |
| | 0.3-1.2 | 28,406 | 12,918 | 6,583 / 6,335 | +1.92 +- 0.89 % | +3.5 +- 1.5 % |
| Euclid EDF-N (VIS) | 0.3-0.6 | 57,059 | 36,342 | 18,233 / 18,109 | +0.34 +- 0.52 % | +1.1 +- 1.0 % |
| | 0.6-0.9 | 61,620 | 37,442 | 18,570 / 18,872 | -0.81 +- 0.52 % | -0.2 +- 1.0 % |
| | 0.9-1.2 | 33,268 | 17,519 | 8,954 / 8,565 | +2.22 +- 0.79 % | +3.4 +- 1.3 % |
| | 0.3-1.2 | 155,317 | 93,306 | 46,728 / 46,578 | +0.16 +- 0.33 % | +0.9 +- 0.9 % |
| Euclid EDF-F (contains GOODS-S) | 0.3-0.6 | 44,589 | 28,035 | 14,043 / 13,992 | +0.18 +- 0.63 % | +0.9 +- 1.0 % |
| | 0.6-0.9 | 32,475 | 20,204 | 10,092 / 10,112 | -0.10 +- 0.75 % | +0.6 +- 1.1 % |
| | 0.9-1.2 | 16,027 | 8,274 | 4,122 / 4,152 | -0.36 +- 1.11 % | +0.4 +- 1.6 % |
| | 0.3-1.2 | 95,136 | 57,788 | 28,896 / 28,892 | +0.01 +- 0.45 % | +0.7 +- 0.9 % |
| Euclid EDF-S | 0.3-0.6 | 112,334 | 70,387 | 35,375 / 35,012 | +0.52 +- 0.39 % | +1.3 +- 0.9 % |
| | 0.6-0.9 | 78,548 | 48,728 | 24,637 / 24,091 | +1.12 +- 0.48 % | +2.0 +- 0.9 % |
| | 0.9-1.2 | 35,698 | 18,608 | 9,441 / 9,167 | +1.47 +- 0.73 % | +2.5 +- 1.3 % |
| | 0.3-1.2 | 231,243 | 140,558 | 70,843 / 69,715 | +0.80 +- 0.28 % | +1.6 +- 0.8 % |
| Euclid, three fields | 0.3-1.2 | 481,696 | 291,652 | 146,467 / 145,185 | +0.44 +- 0.19 % | +1.2 +- 0.8 % |
| JADES F200W, Shamir's 263 (his positions) | his 0-3 | 262 | 182 | 110 / 72 | +20.9 +- 7.2 % (his labels on the same 182: 111 / 71) | per-galaxy agreement 181/182 |
| JADES F200W, whole mosaic, our selection (added after the row above) | all | 2,259 | 997 | 514 / 483 | +3.1 +- 3.2 % | |
| ... of which in Shamir's list | | | 163 | 99 / 64 | +21.5 +- 7.6 % | agreement 162/163 |
| ... of which not in Shamir's list | | | 834 | 415 / 419 | -0.5 +- 3.5 % | |

Classifier B (Fourier) raw A, 0.3-1.2: COSMOS +2.5 +- 2.7 % (1,383); EDF-N -0.54 +- 0.60 %; EDF-F -0.94 +- 0.78 %;
EDF-S +0.39 +- 0.50 %; Euclid all -0.17 +- 0.34 % (84,465). Where both classify, the CNN and B agree on 94 % (COSMOS,
958) and 87 % (Euclid, 62,386). B fails its Euclid injections (M3), so it is not used for verdicts.

**Pre-registered tests (99.7 % intervals, corrected values):**

| test | claim | ours | verdict |
|---|---|---|---|
| (i) COSMOS A vs P1 | +2.48 % | +3.5 +- 1.5 % (raw +1.9 +- 0.9 %) | not constraining (interval holds 0 and the claim) |
| (ii) EDF-F A vs P2 | +20.2 % (all z); +9.3 % at 0.3<z<1.2 (his own 86) | +0.7 +- 0.9 % | both excluded |
| (iii) axis G dipole, Euclid fields only (no intercept) | A_d = 24.8 % (11.4 % at matched z) | A_d = +0.95 +- 0.73 % | both excluded; consistent with 0 |
| (iv) axis H dipole | A_d = 10.5 % (5 %) | A_d = +1.4 +- 1.9 % | 10.5 % excluded; 5 % not constraining |
| axis S dipole (reported) | | A_d = -0.1 +- 1.2 % | |

All four fields with a free HST-minus-Euclid offset give the same A_d, because the offset absorbs the single COSMOS
point; that offset is +4.1 +- 1.5 % on axis G. Free axis (Euclid + COSMOS, 2-degree grid): best |A_d|/sigma = 2.85 at
(284, -27). Its global p = 0.023 from 300 noise realisations, so it is not significant. Redshift growth, dA/dz from the three bins
(corrected): EDF-N +2.9 +- 2.6 %, EDF-F -0.8 +- 3.0 %, EDF-S +2.1 +- 2.5 %, Euclid all +1.9 +- 2.2 %, COSMOS -0.1 +- 5.0 %
per unit z. No bin departs from zero at the 3.69 sigma bin-level bar. The largest departure is Euclid's 0.9-1.2 bin: raw
+1.43 +- 0.47 % (3.0 sigma binomial), corrected +2.5 +- 1.1 % (2.3 sigma).

## R2. Nulls

1. **Mirror:** the whole classifier on mirrored inputs returns exactly -f: max |f + f_mirror| = 0 to 6e-7, and Fourier
   to 5e-15, on every dataset. prep is mirror-equivariant to float rounding (max difference 0.0 on the tests).
2. **Stars:** COSMOS fails at point-source scale. 549 of 9,911 ACS point sources pass t_A (5.5 %), with A = -23.5 +- 4.3 %;
   mean f over all stars is -0.199 +- 0.014. Euclid stars fall below prep's R50 floor (30 of 800 enter, none pass), so
   point sources never reach the Euclid classifier.
3. **Ellipticals:** COSMOS +1.7 +- 3.6 % (765 of 2,801 pass); Euclid -4.2 +- 4.4 % (505 of 1,998). Both pass.
4. **Achiral models with the real PSF:** COSMOS shows no chirality, with either the median PSF or 1,730 individual star
   kernels (M2, M3). Euclid VIS is chiral to the classifier at galaxy scale: mean f +0.115 +- 0.024 against -0.216 +- 0.024
   with the mirrored PSF. Among models that pass t, A = -7 +- 8 %.
   4b. **Real galaxies re-convolved** with the PSF and with its mirror (3,000 per dataset): dA = -0.10 +- 0.08 % (Euclid)
   and -0.99 +- 0.12 % (COSMOS). The PSF term in the data is dA/2: under 0.1 % (Euclid) and about 0.5 % (COSMOS), subtracted.
5. **Size:** no trend. Euclid R50/FWHM quartiles give +1.0, +0.4, +0.4, +0.1 % (each +-0.34 %; Q1 - Q4 = 0.9 +- 0.5 %).
   COSMOS gives +0.3, +1.7, +2.1, +2.7 % (each +-1.7 %). **Tile quadrant:** no trend (Euclid +0.50, +0.52, +0.25, +0.49 %;
   COSMOS +2.4, +3.2, -0.6, +2.6 %).
   **Position angle (image frame): a trend, odd under the mirror.** Euclid: PA 0-45 +3.6 %, 45-90 +1.3 %, 90-135 -0.6 %,
   135-180 -2.5 %. COSMOS: +7.0, +1.2, +2.4, -2.9 %. JADES field: the same shape. The CNN is exactly mirror-antisymmetric
   but not rotation-invariant, so its error rate can differ between PA theta and -theta. That makes A(theta) odd under
   theta -> 180 - theta, and it cancels over a PA distribution that is itself mirror-symmetric. The mirror-paired bins
   hold equal counts: Euclid 72,434 vs 72,247 and 73,531 vs 73,440. This is a classifier property, not a sky signal. It
   is what the exact antisymmetry protects against, and it is listed as a finding about the method.
6. **Injections by true handedness** (M2, M3): 5 % and 20 % recovered in COSMOS and Euclid.
7. **JADES per galaxy:** our CNN classifies 182 of Shamir's 263 and agrees with his label on 181.

## F2. Finding check (blind-review.md 2; fresh context-free Opus agent, the section alone, prompt verbatim; 2026-11-02)

VERDICT: novel but not sound: the claim that the JADES excess is confined to the published selection rests on a post-hoc full-field test with an uncalibrated HST threshold, no NIRCam PSF or eps calibration, and a comparison against galaxies the published pipeline rejected, so the gap between the two subsets (~2.7σ, about the original result's own significance) cannot locate the excess, while the Euclid corrected values hinge on a model-based correction for a strongly chiral VIS PSF and a misclassification rate calibrated only on simulations.

Strongest objection (its words): the JADES full-field comparison is "against galaxies the published pipeline rejected",
which "are not a sample of the same population", and the eps is "calibrated only on simulations"; it asked for A as a
function of classifier confidence. **Parked** (physics point; the rule allows one revision for presentation only).
Post-verdict diagnostic, recorded in finding.md: Euclid A by |f|: -0.15, +0.43, +0.53, +0.15, +0.28, +0.56 % (bins of
40/20/20/10/5/5 % of the sample, +-0.23 to 0.64 %); the most confident 5 % alone exclude +20 % and +9.3 % without any
correction.

## V1. Revision (R-161, cosmo-deep-25b, 2026-11-02): the raw result restated, with the worst-case correction bound

Brief `briefs/r163_brief.md` (cloud worker cw-168). Presentation only: no refetch, no new classification; every number below
comes from R1/R2 above. Script: `r158_raw.py` (raw dipoles with the same per-field no-intercept fit as `r158_stats.py`; the
dilution each excluded amplitude would need). Changes to the `finding.md` section: the headline is the raw Euclid result;
the JADES full-field localisation is removed from the findings (the 181/182 per-galaxy agreement stays, as a classifier
check); COSMOS stays "not constraining"; the ~1 % raw excesses are not interpreted.

Output of `python3 r158_raw.py`:

```
bias bound: null 4b PSF term -0.05 +- 0.04 %, 99.7 % |b| <= 0.17 %, allowance B = 0.51 %
JADES real-galaxy eps (1/182 disagree): 99.7 % upper 0.048
cos(GOODS-S, SGP) = 0.813
raw dipole axis G (SGP): A_d = +0.55 +- 0.28 %  cos {'EDF-N': np.float64(-0.499), 'EDF-F': np.float64(0.816), 'EDF-S': np.float64(0.733)}  uniform-bias projection +0.71
raw dipole axis H (HST best fit): A_d = -0.31 +- 0.74 %  cos {'EDF-N': np.float64(0.397), 'EDF-F': np.float64(0.201), 'EDF-S': np.float64(-0.114)}  uniform-bias projection +1.76
raw dipole axis S (SDSS z>0.15): A_d = -0.52 +- 0.45 %  cos {'EDF-N': np.float64(0.613), 'EDF-F': np.float64(-0.005), 'EDF-S': np.float64(-0.337)}  uniform-bias projection +0.27

| test | published | ours, raw | raw 99.7 % upper | largest 1 - 2 eps that hides it (bias at the bound) | eps needed |
|---|---|---|---|---|---|
| EDF-F A | P2, all z = 20.2 % | +0.01 +- 0.45 % | +1.36 % | 0.093 | 0.45 |
| EDF-F A | P2 at 0.3<z<1.2 = 9.3 % | +0.01 +- 0.45 % | +1.36 % | 0.201 | 0.40 |
| Euclid, three fields, A | P2 at 0.3<z<1.2 = 9.3 % | +0.44 +- 0.19 % | +1.01 % | 0.163 | 0.42 |
| dipole, axis G | implied = 24.8 % | +0.55 +- 0.28 % | +1.40 % | 0.071 | 0.46 |
| dipole, axis G | at matched z = 11.4 % | +0.55 +- 0.28 % | +1.40 % | 0.154 | 0.42 |
| dipole, axis H | implied = 10.5 % | -0.31 +- 0.74 % | +1.92 % | 0.268 | 0.37 |
| dipole, axis H | proposal = 5 % | -0.31 +- 0.74 % | +1.92 % | 0.563 | 0.22 |
most confident 5 % (|f| > 12.9): +0.56 +- 0.64 %: hides 9.3 % only if 1-2eps <= 0.322 (eps >= 0.34)
largest R50/FWHM quartile: +0.10 +- 0.34 %: hides 9.3 % only if 1-2eps <= 0.175 (eps >= 0.41)
```

- **Raw dipoles** (new numbers; R-158 gave only the corrected ones): axis G +0.55 +- 0.28 %, axis H -0.31 +- 0.74 %, axis S
  -0.52 +- 0.45 %. A uniform bias enters the no-intercept fit as 0.71 b (G) and 1.76 b (H); the bound uses that.
- **Bias bound from the data:** null 4b's PSF term -0.05 +- 0.04 % gives |b| <= 0.17 % at 99.7 %; the table allows 3x
  (0.51 %), since the re-convolution is not the native PSF. Size quartiles (largest +0.10 +- 0.34 %) and tile quadrants
  (+0.25 to +0.52 %) show no trend.
- **Misclassification only dilutes** (mirror-odd classifier): every excluded amplitude needs eps >= 0.37 on real Euclid discs;
  model spirals give 1.7 %, JADES real galaxies <= 4.8 % (1/182, one-sided 99.87 % Clopper-Pearson), the most confident 5 %
  would need eps >= 0.34. The 5 % HST-axis dipole stays "not constraining" (pre-registered verdict, from corrected values),
  although the raw upper limit (+1.92 %) is below it: hiding it needs only eps >= 0.22.

## V2. Finding check on the revised section (R-161; fresh context-free Opus agent, section alone, prompt verbatim, hints as R-158's)

VERDICT: novel but not sound: the headline exclusions target extrapolated or redshift-mismatched amplitudes (the JADES galaxies' own matched-z value, +9.3 ± 10.7 %, agrees with Euclid), while the amplitudes actually published at these redshifts are left undecided and the Euclid misclassification rate is never measured on real VIS spirals.

Strongest objection (its words): the exclusions refute "amplitudes nobody claimed at this depth and redshift"; the CNN-Fourier
87 % agreement bounds eps (~0.1) only if their errors are independent. **Parked** (physics; the one revision is spent). What it
credits: an upper limit of ~1-2 % on |A| at 0.3 < z < 1.2 in Euclid Q1, valid for eps < ~0.37. Path: classify on VIS the JADES
galaxies with Euclid counterparts (113 matched, 86 at 0.3 < z < 1.2) for a real-galaxy eps; frame as an upper limit at matched
redshift and selection.
