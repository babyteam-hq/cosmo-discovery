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
