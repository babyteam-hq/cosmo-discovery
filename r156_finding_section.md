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
