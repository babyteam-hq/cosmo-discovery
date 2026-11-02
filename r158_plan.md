# Analysis plan: spiral winding-sense asymmetry at 0.3 < z < 1.2 on space-based images

Fixed before any target galaxy is classified. Nothing below changes after the first target classification except
where a numbered rule says how.

## Question and claim being tested

Is the excess of one spiral winding sense reported on space-based images real? The claims under test, with their
published numbers (fixed now):

- **P1, HST COSMOS (arXiv:2011.03723, Table 2).** ACS F814W, manual annotation: 3,116 clockwise vs 2,965
  counterclockwise, A = (N_cw - N_ccw)/(N_cw + N_ccw) = +2.48 % +- 1.28 % (binomial). Redshifts mostly z < 1.2.
- **P2, JWST JADES GOODS-S (arXiv:2502.18781).** 158 clockwise vs 105 counterclockwise (Ganalyzer), A = +20.2 % +- 6.0 %,
  at (53.12, -27.81); redshifts spread over 0-3. The author attributes it to a cosmological axis near the Galactic pole with
  a clockwise excess toward the south Galactic pole, growing with redshift.
- **P3, dipole along fixed axes.** A(n) = A_d cos(phi), phi the angle between the field and the axis. Axes fixed now:
  - axis H: (RA, Dec) = (78, +47), the HST best fit (2011.03723);
  - axis S: (71, +61), the SDSS z > 0.15 best fit quoted there;
  - axis G: the south Galactic pole (12.86, -27.13), clockwise excess toward it (2502.18781).
  The amplitudes these claims imply at our fields (cos phi from the field centres):

  | field (centre) | cos phi, axis H | axis S | axis G |
  |---|---|---|---|
  | COSMOS (150.12, +2.20) | 0.237 | 0.125 | -0.671 |
  | Euclid EDF-North (269.73, +66.02) | 0.397 | 0.613 | -0.499 |
  | Euclid EDF-Fornax (52.93, -28.09) | 0.201 | -0.005 | 0.816 |
  | Euclid EDF-South (61.24, -48.42) | -0.114 | -0.337 | 0.733 |

  Implied A_d: axis H 10.5 % (P1's +2.48 % / 0.237); axis G 24.8 % (P2's +20.2 % / 0.813). The proposal's own HST anchor
  is a ~5 % dipole; it is tested too.

Sign convention: "clockwise" is Galaxy Zoo 1's clockwise as seen on the sky, the naming Shamir uses when he compares with
GZ1 (2502.18781 calls GZ1's excess "counterclockwise", as GZ1 does). The mapping is checked per galaxy on Shamir's own
JADES list (sample 3), where available.

## Samples (cuts fixed)

1. **COSMOS, HST ACS F814W** (IRSA acs_2.0 tiles, 0.03"/px, PSF FWHM ~0.095"). ZEST morphology (Scarlata et al. 2007)
   matched at 0.5" to COSMOS2020 Classic LePhare photo-z. Cuts: 0.3 < z_phot < 1.2; ZEST type 2 (disc); r50 > 9.5 px
   (0.285" = 3 PSF FWHM); ACS elongation a/b < 2.5 (b/a > 0.4); I_AB < 24 (ZEST limit); junkflag 0. N = 28,710 before
   classification (9,955 / 10,800 / 7,955 in the three redshift bins).
2. **Euclid Q1 VIS** (MER BGSUB mosaics via IRSA cutouts, 0.1"/px, PSF FWHM ~0.17"). MER catalogue x MER morphology
   (single-Sersic VIS fit) x PHZ. Cuts: 0.3 < phz_median < 1.2; VIS total flux > 2.29 uJy (VIS < 23); Sersic VIS
   radius > 0.5" (3 PSF FWHM); Sersic axis ratio > 0.4; Sersic index < 2.5; point_like_prob < 0.5; spurious_prob < 0.5.
   Sersic radius < 3" (failed fits). Fields EDF-N, EDF-F, EDF-S (LDN1641 excluded: Galactic plane). Counts in
   `results-R158.md` F1 (~2.5e5 pass). Stamps come from IRSA's per-object cutout service (~40/s), so the catalogue is put in
   a random order (seeded) and fetched for a fixed 75 minutes; the fetched set is a random subsample, which is
   mirror-invariant. Stamp side = clip(8 r_e, 6", 20").
3. **JADES GOODS-S, JWST NIRCam F200W** (DR2 mosaic, 0.03"/px, already fetched: all 263 positions inside): Shamir's 263
   galaxies at his published coordinates (2502.18781 Tables 1-2; 105 ccw, 158 cw), 9" stamps, reclassified with our
   pipeline. Reported: A over his 263 (all, and those our classifier classifies), and the per-galaxy agreement with
   his labels (a binomial test against 50 %). This is the one place his annotation step is tested galaxy by galaxy.

No selection uses a handedness-sensitive quantity: every cut is on mirror-invariant catalogue quantities (sizes, axis
ratios, Sersic index, fluxes, photo-z). Zoobot's Euclid vote fractions are not used for selection (they come from a
network whose mirror symmetry is not guaranteed).

## Image handling (exactly mirror-equivariant)

- Stamps are flipped to sky parity (det CD < 0) from their own WCS.
- prep: border-median background, Gaussian-smoothed 2.5 sigma segmentation, the segment at the centre kept and other
  segments zeroed, flux-weighted R50 from the stamp, bilinear resampling onto a 64 x 64 grid spanning +-3 R50 about the
  stamp centre (the catalogue position), arcsinh stretch. The grid is symmetric about the centre, so
  prep(mirror(stamp)) = mirror(prep(stamp)) to float rounding (checked: max difference 0.0).
- Objects with stamp-measured R50 < 1.5 px of the original image are dropped (mirror-invariant).

## Classifiers

- **A (primary): reflection-antisymmetric CNN.** g is a 6-conv-layer CNN (12-48 channels, batch norm); the classifier is
  f(x) = [g(x) - g(mirror(x))]/2, so f(mirror x) = -f(x) exactly for any weights (checked on an untrained net:
  -4.97e-5 vs +4.97e-5). Rotation invariance is learned (rotation augmentation); exactness needs only the reflection. Trained on Galaxy Zoo 1 (Lintott et al. 2011, table 2) spirals with
  P_CW > 0.8 or P_ACW > 0.8 and >= 15 votes, 7,000 of each, Legacy Surveys DR10 r-band stamps (0.262"/px), with
  random rotation, +-15 % scale, Gaussian blur 0-3 grid px and noise (to reach z ~ 1 space-telescope resolution and
  S/N), and mirror-with-label-flip. Training-label chirality bias (GZ1's known anticlockwise excess) cannot enter:
  f is odd under mirror for every weight setting, so a mirror-symmetric population has E[A] = 0 exactly.
  A galaxy is classified when |f| > t_A; t_A is fixed on the held-out GZ1 set (degraded to each dataset's resolution)
  as the smallest |f| where accuracy >= 90 %, before any target galaxy is classified.
- **B (second): logarithmic-spiral Fourier power** (2DFFT-style, Davis et al. 2012). On the same 64 x 64 image, in the
  annulus 0.5-2.8 R50 with (ln r, theta) weights, S = (P+ - P-)/(P+ + P-), with P+- the power at pitch +-p (p = 1..20)
  summed over m = 1, 2, 3. Exactly odd under mirror. Classified when |S| > t_B, t_B fixed the same way (accuracy >= 85 %).
  Its sign convention is fixed by its correlation with the GZ1 labels on the training set.
- Agreement of A and B on the target galaxies is reported (fraction same sign, Cohen's kappa).

## Statistic

- **Per dataset/field and redshift bin** (0.3-0.6, 0.6-0.9, 0.9-1.2, and all): A = (N_cw - N_ccw)/(N_cw + N_ccw) of the
  galaxies classifier A classifies; binomial 68/95 % intervals and a block bootstrap over image tiles (ACS tiles, MER
  tiles; 1,000 resamples). The larger of the two is quoted.
- **Classification-error dilution.** A_obs = (1 - 2 eps) A_true, where eps counts both misread spirals and
  non-spirals classified by chance; the second is measured by the rate at which the achiral simulations (null 4) pass
  the threshold, with their sign random (a sign bias there is the PSF term, not dilution). eps per dataset from the GZ1 held-out set degraded to
  that dataset's median R50/PSF and S/N, cross-checked against the A-B agreement rate (for independent errors
  agreement = (1-eps_A)(1-eps_B) + eps_A eps_B). Both raw and corrected values are reported; the comparison with P1/P2/P3
  uses the corrected value and its interval, with eps's own spread propagated.
- **Dipole amplitude A_d** along axes H, S, G: weighted least squares of the per-field A (corrected) on cos(phi), per
  redshift bin and all; interval from the per-field intervals. With four fields, a free-axis fit is reported only as a
  secondary with its look-elsewhere: the maximum |A_d|/sigma over a 1-degree HEALPix-like grid of axes, calibrated on
  1,000 sign-shuffled catalogues (labels permuted within fields keeps no dipole; global p from the shuffles).
- **Pre-registered tests and decision.** Primary: (i) COSMOS A vs P1 (+2.48 %) and vs 0; (ii) EDF-F A (covers GOODS-S)
  vs P2 (+20.2 %) and vs 0; (iii) A_d on axis G vs 24.8 % and vs 0; (iv) A_d on axis H vs 10.5 % and 5 %. A claim is
  "excluded" when its value is outside our 99.7 % interval, "consistent" otherwise. Three redshift bins x four primary
  tests: any bin-level departure from zero is judged against a 12-trial Bonferroni bar (2.87 sigma for 3 sigma global).
- **Redshift growth.** Slope dA/dz from the three bins per field, vs 0.

## Nulls and chirality tests (each must pass, or the dataset's result is reported as not separated from it)

1. **Mirror:** the whole pipeline on mirrored stamps returns exactly -A (to float rounding) and the same classified set.
2. **Stars:** point sources (COSMOS acs_class_star > 0.95, 19 < I < 23; Euclid point_like_prob > 0.95, 19 < VIS < 23)
   through prep + both classifiers: the classified fraction and A, which must be consistent with 0 (they test PSF
   chirality directly).
3. **Ellipticals:** COSMOS ZEST type 1; Euclid Sersic n > 3.5, axis ratio > 0.6: A consistent with 0.
4. **Simulated achiral galaxies convolved with the real PSF:** exponential discs + bulges + an m = 2 bar or ring (each
   model mirror-symmetric about its own axes), at the sample's R50 and axis-ratio distributions and random position
   angles, convolved with the empirical PSF (median of bright unsaturated stars from the same images), with noise from
   real empty-sky stamps, through prep + both classifiers. The mean f and A must be consistent with 0; a nonzero value
   is a measured PSF-induced chirality, in bins of R50/FWHM.
4b. **Real galaxies re-convolved with the PSF and with the mirrored PSF:** A(stamp x PSF) - A(stamp x mirror(PSF)) on
   the classified galaxies measures how strongly a PSF chirality of the instrument's size moves A; half of it bounds
   the PSF term already in the data.
5. **A vs R50/FWHM** (quartiles), **vs position on the tile** (COSMOS: ACS tile x/y quadrants; Euclid: MER tile
   quadrants; MER mosaics are stacks of dithers, so detector position is not available per galaxy: stated, not
   tested), and **vs galaxy position angle relative to the pixel grid** (four bins): no trend beyond the errors.
6. **Injection:** a known asymmetry (A = 5 % and 20 %) injected by mirroring a random subset of classified stamps and
   re-running: recovered A within errors, confirming the dilution model and the intervals.
7. **Classifier B agreement:** A from classifier B (and from galaxies both classify, with the same sign) reported
   alongside; a disagreement beyond errors is reported as such.

## Error budget and dependencies

- Statistical: binomial / tile bootstrap. Systematic: eps (dilution) from the degraded GZ1 set; the PSF-chirality
  floor from null 4 (its interval is added in quadrature if it is not zero); photo-z errors only move galaxies between
  bins (bins are reported with their photo-z 68 % widths; they do not change the sign of A).
- Parameters the result depends on: t_A and t_B (fixed on GZ1 before classification; a +-25 % variation of each is
  reported), the redshift-bin edges (fixed above), the axes and amplitudes (fixed above).

## Method-check closures (2026-11-02, before any target galaxy is classified)

The method check returned "fixable" (verdict verbatim in `results-R158.md` M1). Each point, as now fixed:

1. **eps on the target population.** Chiral model galaxies (exponential disc modulated by an m = 2 logarithmic spiral,
   arm contrast 0.15-0.8, pitch 10-35 deg, bulge, 0-5 random clumps, q 0.4-1, random PA; true handedness known) are
   drawn at each dataset's own R50 and S/N distribution (from prep, sign-free), convolved with the dataset's empirical
   PSF, noised, and classified (`r158_sims.py`, "chiral"). eps_spiral is measured per dataset in the same R50 and S/N
   mix. The armless-disc term is bounded conservatively: r_a (the achiral pass rate, null 4) x N_sample / N_classified,
   i.e. as if every unclassified galaxy were armless. eps = eps_spiral (1 - w) + 0.5 w with w that bound. **Floor:** a
   bin with (1 - 2 eps) < 0.3 is reported as "not constraining". The sign mapping of the network to the model
   handedness comes from the brightest third of the injections; GZ1's label convention gives the mapping to "clockwise".
   No human labels are used.
2. **Injection by true handedness.** The 5 % and 20 % injections are done on the chiral models (true handedness drawn
   with P(+) = (1 + A_in)/2), not by mirroring classified stamps. Recovered A / (1 - 2 eps) must match A_in.
3. **P2 at matched redshift.** Shamir's 263 JADES galaxies are matched (1") to the Euclid Q1 PHZ catalogue (EDF-F
   covers GOODS-S). P2's target value for test (ii) becomes his own A over his galaxies with 0.3 < z < 1.2, with its
   binomial interval; the 0-3 figure (+20.2 %) is kept as the headline comparison. If fewer than 30 match, a growth law
   fixed now is used instead: A(z) linear through his HST value (+2.48 % at the COSMOS median z ~ 0.7) and his JADES
   value (+20.2 % at z ~ 1.5, the middle of his histogram), evaluated at each bin centre. Axis G's implied A_d is scaled
   the same way.
4. **Sign convention.** P1's galaxy list is not published, so tests (i) and (iv) are reported as conditional on
   Shamir's naming matching GZ1's (which he states when he compares with GZ1). The P2 mapping is tested per galaxy.
   Training-stamp parity checked: Legacy Surveys cutouts have CD1_1 < 0, CD2_2 > 0 (det < 0), the same sky parity the
   target stamps are flipped to.
5. **Trials.** The Bonferroni bar is 3.69 sigma per trial for 3 sigma global over 12 bin-level trials. The four primary
   tests (i)-(iv) each get a verdict at 99.7 % on its own interval, and the report says how many tests were made;
   bin-level departures from zero are judged at 3.69 sigma.
6. **Three verdicts:** "excluded" (claim outside our 99.7 % interval, interval excludes it), "consistent" (claim inside
   and 0 outside, or claim inside and the interval narrower than |claim|), "not constraining" (interval contains both
   0 and the claim). The power estimate (expected interval per test) is recomputed with the measured eps and the
   fetched counts before classification and written to `results-R158.md` M2.
7. **Dipole fit.** No intercept. Primary fit: the three Euclid fields alone (one instrument); COSMOS is an independent
   check against the Euclid A_d. Secondary: all four fields with a free HST-minus-Euclid offset.
8. **PSF term.** A nonzero null-4 or 4b result is subtracted as a signed bias per R50/FWHM bin (4b's response on real
   galaxies), with its error propagated, not added in quadrature.
9. **Limitations stated, not tested:** one median PSF per dataset (no field variation across ACS tiles or Euclid stacks;
   null 2 covers point-source scale only); Euclid detector position is not available in MER mosaics.

## Threshold rule revised before classification (2026-11-02, after the GZ1 and simulation runs, no target sign seen)

The GZ1 rule ("smallest |f| with held-out accuracy >= 90 %") is degenerate: the GZ1 set holds only clear spirals, so
accuracy is >= 90 % at every |f| (t_A = 0.004 clean, 0.007 at COSMOS S/N) and every target galaxy, armless or not,
would be classified, making the closure-1 armless bound w = 1 (nothing constraining). Replaced, for both classifiers and
every dataset, by: **t = the |score| above which 5 % of the achiral simulations (null 4, that dataset's PSF, R50 and
S/N) pass.** It uses only simulations. Everything else in the closures stands (eps from the chiral injections at the
same t; the armless bound w = 0.05 x N_sample / N_classified). COSMOS: t_A = 2.89 (chiral pass 81 %, chiral error
1.1 %), t_B = 0.42 (chiral pass 22 %, error 1.9 %). Classifier B's sign is fixed as -S (its bright-injection agreement
is 65 %, GZ1 agreement 76 %); the CNN's sign maps +f to GZ1 clockwise and to the model handedness the injections call +1.
