# The near-infrared background traced by DESI DR1 spectroscopic slices in SPHEREx Quick Release images: b_I d(νI_ν)/dz is detected at 3-17σ in every one of 18 channel groups from 0.75 to 2.42 µm, at 4-9 nW m^-2 sr^-1 per unit redshift and nearly flat over 0.1 < z < 2.1, with no Hα, [OIII] or Pa-α line detected and values 4-25x below the CIBER x DESI Legacy Survey measurement at 1.8 µm and z < 0.6

## Claim

Cross-correlating 22,322 SPHEREx QR2 wide-survey spectral images (detectors 1-3, 0.75-2.42 µm; zodiacal model
subtracted, catalogue sources masked) with 3.0 million DESI DR1 spectroscopic galaxies and quasars in 32 slices of
Δz = 0.1 over 1,252 deg², the clustering-redshift estimator returns the bias-weighted background intensity per unit
redshift, b_I d(νI_ν)/dz, per channel group (6 groups of 2-3 native channels per detector) and slice:

- detected in every group, 3-17σ per group (all 32 slices), at 4-9 nW m^-2 sr^-1 per unit z: e.g. 1.81 µm 6.0 ± 0.5,
  1.94 µm 6.1 ± 0.4, 1.41 µm 7.7 ± 0.9, 0.96 µm 6.7 ± 1.3 (inverse-variance means over slices, jackknife σ);
- nearly flat in wavelength and in redshift over 0.1 < z < 2.1; at λ > 1.8 µm the z < 0.5 slices are lower (3.5-4.8)
  than z > 0.5 (6.4-7.3); the lowest slice (0.05-0.15) is consistent with zero in every group;
- no line detected: the joint Hα amplitude statistic over 16 groups (χ² 31.2, p = 0.013) is matched by the same
  statistic on redshift-shuffled galaxies (32.0), which carry no line; [OIII] and Pa-α are null; per-group 1σ limits on
  a line component are 3-14 nW m^-2 sr^-1 per unit z;
- at 1.8 µm the values at z 0-0.2, 0.2-0.4, 0.4-0.6 (3.4 ± 2.7, 5.4 ± 1.2, 6.9 ± 1.5) are 4-25x below CIBER x DESI-LS
  (arXiv:2608.12116: 140 ± 40, 63 ± 15, 35 ± 8; 3.4-3.9σ per bin), and consistent with its z > 0.6 limits; ours lie in
  the band of the integrated-galaxy-light models that paper shows (5-20).

The claim is the table (results-R138.md, "Result"); nothing is claimed about galaxy evolution.

## Data and method

- **Images.** SPHEREx QR2 L2 spectral images (IRSA, `nasa-irsa-spherex` S3), every D1-D3 image centred in the field or
  one NSIDE-32 pixel around it (32,005), processed in a fixed random order until a pre-set stop: 22,322 (70 %, ~2.3
  visits per channel per point). Per image: flags (bad bits and the pipeline SOURCE bit, which masks every Reference
  Catalog point source with a magnitude-dependent radius; ~31 % of pixels), pixels > 5 robust σ above the channel
  strip's median dilated by 12"; the pipeline ZODI model subtracted; each channel strip's (0.2° x 3.5°) median of good
  pixels subtracted per exposure; 8x8-pixel block means binned to NSIDE 2048 per channel (17 log-uniform channels per
  detector, R ~ 43) and per half-survey (split at the median exposure time). Map level: Gaia DR3 G < 10 stars masked
  (r = 2' x 10^(0.2(10-G))).
- **Galaxies.** DESI DR1 NGC LSS clustering catalogues (BGS_ANY 0.05-0.45, LRG 0.4-1.1, ELG_LOPnotqso 0.8-1.6, QSO
  0.8-2.1, Δz = 0.1) and random file 0 of each, with WEIGHT; field: 1,252 deg² where BGS, LRG and ELG are all ≥ 80 %
  complete and |b| > 40°.
- **Estimator (arXiv:1810.00885).** w_Ir(θ, z) = mean residual intensity at separation θ around galaxies minus around
  randoms; w̄_Ir = ∫ W(θ) w_Ir dθ over an annulus fixed in comoving distance, θ = max(1 cMpc, 3.5') to 10 cMpc, with
  W ∝ θ^-0.8; b_I dI/dz = w̄_Ir / (b_r w̄_m), w̄_m from CAMB halofit P(k) averaged over the same pixel separations;
  b_r from the slice's own w̄_gg (Hamilton estimator, same annulus) against Limber for the slice's n(z) (b_r: BGS 1.0-1.4,
  LRG 1.9-2.3, ELG 1.1-1.8, QSO 1.8-3.0). Pair sums on 1.72' gnomonic grids, one per NSIDE-8 patch (51), by FFT ring
  convolution.
- **Strip-filter transfer, measured.** Each slice's own galaxy overdensity map is passed through the identical
  per-(image, channel) mean subtraction over a random 3,000-image subset and the estimator run on it filtered and
  unfiltered; the ratio (0.44-0.76) divides w̄. QSO slices take the transfer of the same-z ELG slice (their own
  template is shot noise).
- **Errors:** delete-one jackknife over the 51 patches; the half-survey difference scatters at 0.54-0.88 of it.
- **Nulls:** zodi model alone through the same pipeline (χ² 24-40 / 32 per group: pass); half-survey difference (21-53
  / 32: pass); redshift-shuffled galaxies flat in z within each tracer (18-46 / 28: pass); Gaia 12 < G < 16 stars in
  place of galaxies: not zero (χ² 40-269 / 32), i.e. light associated with stars at 3.5'-1° (PSF wings, ghosts, star
  clustering); splitting the patches at the median star count leaves the galaxy signal no higher where stars are dense
  (low-star half 6-14, high-star half 3-9).
- **Line model:** per group, a quadratic in z plus Hα 0.6563, [OIII] 0.5007, Pa-α 1.8751 µm amplitudes times the
  fraction of each slice whose line wavelength falls in the group; σ from jackknife refits; calibrated against the
  shuffled galaxies.
- **Fixed in advance** (results-R138.md): field, images, channels, groups, slices, masks, annulus, estimator, nulls,
  line model, comparison. Changed after an interim run on 9,743 images (all recorded there with reasons): a pair-map
  indexing bug fixed, the Gaia mask radius law (a counting error: as written it masked 88 % of the map), the QSO
  transfer, the shuffled null's expectation, and an added stellar-density split.

## What it means

The near-infrared background that is spatially correlated with spectroscopic galaxies is measured channel by channel
across 0.75-2.4 µm and 0.1 < z < 2.1 from real SPHEREx data, at 4-9 nW m^-2 sr^-1 per unit z after the pipeline's
catalogue source mask. At the precision reached (σ ~ 1-3 per slice-group, ~0.4-1.6 per group), the spectral
structure in z that Hα, [OIII] or Pa-α would add is not seen; a line is bounded at about the continuum level in its
crossing slice. At 1.8 µm and z < 0.6 the amplitude is an order of magnitude below CIBER x DESI-LS and inside the IGL
model band that work found too low; at z > 0.6 the two agree. The SPHEREx source mask removes all catalogue point
sources, deeper than CIBER's, and the scales and tracers differ, so the table does not decide whether the CIBER excess
is light that SPHEREx masks or something else.

## Weakest points

- The SOURCE mask's depth is the pipeline's (Reference Catalog v0.73, spherex_class ≤ 7), not chosen here; b_I dI/dz is
  that of the background left after it, and the CIBER comparison inherits that difference.
- The per-exposure strip median removes ~25-55 % of the signal; the transfer is measured on the galaxies' own
  overdensity (the same angular shape to first order), not on the intensity field itself.
- Stars give a non-zero null; the split by star count argues it does not carry the galaxy signal, but a galaxy-star
  cross-correlation was not measured directly.
- The line test's jackknife errors are not calibrated at the ~2σ level (shuffled galaxies fail it equally), so the
  line limits are indicative.
- CIBER's values are read off its figure (±10 %); our 1.1 µm-exact group is the filter-edge one (σ ~ 40-100); the
  adjacent 1.14 µm group is used beside it.
- Only D1-D3, and 70 % of the images; D4 (Pa-α at z > 0.3) not processed.

## To settle it

D4-D6 and the remaining 30 % of the images; a CIBER-depth source mask on the SPHEREx images (Legacy Survey to CIBER's
cut) to put both on one mask; the galaxy-star cross-correlation; QR3 and the second-year survey for the line search.
