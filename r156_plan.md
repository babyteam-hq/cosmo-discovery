# Analysis plan: the Doppler-magnification dipole of galaxy sizes in DESI DR1 BGS at 0.15 < z < 0.45

Fixed before any dipole is computed.

## Question and claim

At fixed observed redshift, a galaxy's line-of-sight peculiar velocity v changes its apparent size by the Doppler
convergence kappa_v = (c/(H r) - 1) v/c (H conformal, r comoving distance; Bonvin et al. 2017, arXiv:1610.05946).
Galaxies behind an overdensity fall toward us and look smaller; those in front look larger. GR predicts an odd (dipole)
cross-correlation between galaxy sizes and galaxy number counts. The claim will be the dipole amplitude A_D (A_D = 1 is
linear-theory GR with Planck 2018 f sigma_8) per redshift bin and combined, with its 68 % interval, against zero and
against one.

## Samples

- Density tracers: DESI DR1 BGS_BRIGHT LSS clustering catalogues v1.5 (NGC + SGC, 7,180 deg2, 3.96 M galaxies), all
  galaxies with 0.08 < z < 0.56, weight WEIGHT; randoms: the first random file per cap (18.7 M points, 4.7x the data),
  same z range, weight WEIGHT.
- Size galaxies: the same catalogue's galaxies at 0.15 < z < 0.45 matched by TARGETID to the DESI DR1 LS-DR9 Tractor
  photometry VAC (`tractorphot`), kept if Tractor TYPE is REX/EXP/DEV/SER, SHAPE_R > 0, SHAPE_R_IVAR > 0,
  FRACMASKED_R < 0.1, MASKBITS free of BRIGHT/GALAXY/CLUSTER (bits 1, 12, 13), and sigma_lambda < 0.5.
  (97 % of the clustering galaxies pass on the first 20 % of the sky.)
- Redshift bins: 0.15-0.25, 0.25-0.35, 0.35-0.45 (1.41 M, 0.92 M, 0.37 M clustering galaxies).

## Size estimator

- lambda = ln R_phys, R_phys = SHAPE_R x D_A(z_obs) (Planck 2018 background).
- kappa_hat = lambda - <lambda>(z, g-r, r-z | type): per Tractor type, least squares on a cubic in z, quadratics in the
  colours, colour-colour and colour-z cross terms, 4-sigma clipped; **no magnitude term** (at fixed z and colour, sizes
  scale as R ~ L^0.49 in this sample, so a magnitude term removes the magnification signal: its response is 0.03).
  Residuals winsorised at +-3 sigma_tot, then the weighted mean removed in z slices of 0.01 per cap.
- Weights u_j = 1/(sigma_int,type^2 + sigma_lambda^2). sigma_int = 0.32-0.45 by type; the effective per-galaxy noise
  is 0.42 (measured).
- **Response eta(z)**, the mean kappa_hat per unit true convergence. A true kappa scales sizes by 1+kappa and fluxes by
  1+2kappa; the sample is cut at r < 19.5, so magnification also admits intrinsically fainter, smaller galaxies:
  eta = 1 + 2.1715 d<kappa_hat>_u / d m_lim, the second term measured on the data in z slices of 0.02 as the slope of the
  weighted mean kappa_hat of galaxies with r < m_cut against m_cut over 19.3-19.5. Measured on 20 % of the sky:
  eta = 0.50 (z 0.15-0.25), 0.24 (0.25-0.35), -0.05 (0.35-0.45). Systematic: the slope from 19.4-19.5 and from
  19.2-19.5. The template carries eta(z) galaxy by galaxy.

## Statistic

- For every size galaxy j and every density point i with 10 < d < 150 Mpc/h, separation s = x_j - x_i and the
  midpoint line of sight, mu = cos(s, x_i + x_j) (mu > 0: size galaxy behind):
  Q_j(d, mu) = sum_gal w_i - alpha sum_ran w_i, R_j(d, mu) = alpha sum_ran w_i; 14 bins of 10 Mpc/h, 10 bins in mu.
- xi(d, mu) = sum_j u_j kappa_hat_j Q_j / sum_j u_j R_j; multipoles xi_l = (2l+1)/2 sum_mu xi(d, mu) Int_bin L_l dmu.
  Normalising per (d, mu) bin removes the leak of the even part into the dipole from the n(z) gradient across a shell.
- The measurement is xi_1(d) at 20 <= d < 150 Mpc/h (13 bins). The same machinery with kappa = 1 gives the galaxy
  autocorrelation xi_gg(d, mu) (for b).

## Model

xi_1(d) = A_D eta T_D(d) + T_rel(d) for each bin, with
- T_D: Bonvin et al. eq. 28 (distant observer), f(z) = Omega_m(z)^0.55, Planck 2018 linear P(k) (CAMB), averaged over
  the bin's size galaxies with their pair weights sum_mu R_j and eta(z_j); b from the monopole of xi_gg at 30-80 Mpc/h
  against the linear Kaiser monopole (b^2 + 2bf/3 + f^2/5) xi_m, per bin.
- T_rel: the number counts' own Doppler term correlated with any density dependence of sizes. Delta contains
  C(z) v/c, C = -Hdot/H^2 - 2/(rH) + 5s(1/(rH) - 1) - E/H (E the luminosity-function evolution at fixed L, from the
  measured n(z) and s); sizes depend on density as kappa_int = b_s delta. Then T_rel = -C b_s/((1/(rH) - 1)(b + 3f/5))
  x T_D(eta = 1). b_s is measured from the monopole ratio xi_0(Delta kappa)/xi_0(gg) (it absorbs physical
  size-environment dependence and density-dependent blending alike); s(z) = dlog10 N/dm at r = 19.5 is measured per
  z slice (0.30 at z = 0.16 rising to 1.03 at z = 0.44). T_rel is subtracted; its uncertainty (b_s's error and +-50 %
  on C) is propagated into A_D.
- Gravitational lensing: < 3 % of the Doppler dipole at z < 0.45 and d < 150 (Bonvin et al. Fig. 3); ignored and
  stated. Wide-angle corrections to the distant-observer dipole: (d/r)^2, about 2-7 % at d = 150; stated.

## Fit and covariance

- Per bin, A_D by generalised least squares on the 13 d bins with the jackknife covariance (200 equal-area regions over
  both caps, HEALPix-based; Hartlap factor applied), T_rel subtracted first.
- Combined: the inverse-variance combination of the three bin A_D with their jackknife 3x3 covariance.
- Cross-checks of the covariance: 300 size shuffles (kappa_hat permuted among size galaxies within z slices of 0.01 and
  caps; pure noise) and the analytic Gaussian covariance (`r156_theory.covariance`).
- Small scales: the fit starts at 20 Mpc/h; the result at d >= 30 Mpc/h is reported as a robustness check (finger-of-
  God velocities put a galaxy's redshift position and its Doppler size shift in step, a non-linear piece the linear
  template omits).

## Nulls (each must give chi^2 with p > 0.01 against zero, dipole, 20-150 Mpc/h)

1. Size shuffles: the distribution defines the noise; the data's chi^2 is compared to it.
2. Randomised redshifts: the pair counts redone with each size galaxy's z replaced by another's (same cap), kappa_hat
   kept: dipole consistent with zero.
3. PSF: kappa replaced by the standardised residual of ln PSFSIZE_R (galaxy's own PSF FWHM): dipole consistent with
   zero; and the regression slope of kappa_hat on ln PSFSIZE_R times that PSF dipole bounds the PSF contamination.
4. Stars: a PSF map from Gaia-confirmed DESI-observed stars (Tractor PSF type, parallax or proper motion > 5 sigma):
   the HEALPix (nside 128) mean of FIBERFLUX_R/FLUX_R, standardised, at each size galaxy's position: dipole consistent
   with zero.
5. Blending: kappa replaced by the standardised FRACFLUX_R (flux fraction from neighbours): dipole consistent with zero;
   and A_D on galaxies with FRACFLUX_R < 0.05 consistent with the full sample.
6. Octupole xi_3: predicted -(2f/5) ... nu_3 (about a tenth of the dipole): consistent with the prediction and with
   zero.
7. The 0.35-0.45 bin, where eta ~ 0 makes the GR prediction for the size dipole ~ 0: consistent with zero (a dipole
   there would be a systematic, not Doppler magnification).
8. Splits: NGC vs SGC, and BASS/MzLS (north) vs DECaLS (south) photometry: A_D consistent.

## Error budget and dependences

- Statistical: size noise (dominant at d < 50), cosmic variance (low z, large d), shot noise. Forecast with the
  measured sigma and eta: S/N 4.4, 1.3, ~0 per bin, 4.6 combined (`r156_forecast.py`, galaxy-normalised noise).
- The result depends on: eta(z) (measured; systematic from the m_cut range), b (measured), f sigma_8 (Planck GR: this
  is what A_D tests), b_s and C(z) through T_rel (measured, +-50 % on C), the distant-observer template (stated).
