# Results (R-156, cosmo-deep-24, 2026-11-02)

First measurement of Doppler magnification, the odd (dipole) cross-correlation of galaxy sizes with galaxy number
counts, at 0.15 < z < 0.45: DESI DR1 BGS_BRIGHT x Legacy Surveys DR9 Tractor sizes. Trajectory check (2026-11-01,
`open-problems/round4/verdicts.md` r2): "partly novel": the claim is restricted to z = 0.15-0.45 and the raw-size
estimator. This file is the worker's own.

Scripts: `r156_fetch_tractor.py` (stream the DR9 photometry VAC), `r156_sample.py` (size sample), `r156_kappa.py`
(size residual), `r156_theory.py` (template, covariance, octupole, C(z)), `r156_forecast.py` (step 0), `r156_acs.py`
(Tractor size response vs HST/ACS), `r156_pairs.py` (numba pair counter, per-galaxy shell sums), `r156_mock.py`
(lognormal light-cone mocks with velocities), `r156_measure.py` (multipoles, fits, nulls). Plan: `r156_plan.md`.
Data under `data/r156/` (gitignored).

## F0. Access (fetch times; `data/` started empty)

| what | where | size | time |
|---|---|---|---|
| BGS_BRIGHT NGC/SGC clustering catalogues + n(z), LSS v1.5 | data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5 | 463 MB | < 1 min |
| BGS_BRIGHT random file 0, NGC + SGC (18.7 M points) | same | 2.4 GB | ~1 min |
| LS-DR9 Tractor photometry of DESI DR1 targets (`tractorphot`, 104 files), streamed, BGS rows + Gaia stars kept | data.desi.lbl.gov/public/dr1/vac/dr1/lsdr9-photometry/iron/v1.1/observed-targets/tractorphot | 46 GB read, 1.2 GB kept | 29.6 min (6 streams, ~26 MB/s) |
| HST/ACS COSMOS morphology (ZEST/GIM2D), acs_mag_auto < 21.5 | IRSA TAP `cosmos_morph_zurich_1` | 23,546 rows | ~1 min (one TLS reset, retried) |
| arXiv 1610.05946, 2511.17858, 2512.03230 | arxiv.org/pdf | | seconds |

The LSS "full" catalogue carries SHAPE_R but not SHAPE_R_IVAR, PSFSIZE or FRACFLUX; the tractorphot VAC carries all of
them, keyed by TARGETID (every one of the 3.96 M clustering galaxies matched).

## F1. Step 0: feasibility (the gate)

**Counts** (DESI DR1 BGS_BRIGHT clustering, NGC + SGC, 7,180 deg2; size cuts in `r156_sample.py`):

| z | clustering galaxies | with good Tractor sizes | per-galaxy noise of kappa_hat | eta (pair-weighted) |
|---|---|---|---|---|
| 0.15-0.25 | 1,409,665 | 1,363,238 (96.7 %) | 0.422 | 0.49 |
| 0.25-0.35 | 920,571 | 885,168 (96.2 %) | 0.424 | 0.22 |
| 0.35-0.45 | 369,459 | 344,914 (93.4 %) | 0.427 | -0.06 |

- **Size scatter at fixed colour and redshift** (the estimator's noise): 0.32-0.45 by Tractor type, 0.42 effective
  with inverse-variance weights; Tractor's own size error is 0.004-0.007, negligible. 2511.17858 assumed 0.39-0.48
  from cosmoDC2: the data agree.
- **The response eta, the new number.** A true convergence kappa scales sizes by 1+kappa and fluxes by 1+2kappa. At
  fixed z and colour BGS sizes scale as R ~ L^0.49, constant surface brightness, which is the direction magnification
  moves a galaxy. So at the r = 19.5 flux limit, magnification admits intrinsically fainter, smaller galaxies, and they
  cancel most of the size increase: eta = 0.57 at z = 0.16 falling to 0.18 at z = 0.30, 0.01 at z = 0.36 and -0.2 at
  z = 0.44 (`data/r156/eta.npz`; in F2 point 2). A magnitude term in the size relation would remove the signal
  outright (its response is 0.03). 2511.17858 assumed eta = 1 ("assuming no selection effects").
- **Forecast** (`r156_forecast.py`; Bonvin et al. eq. 28 template, Gaussian covariance derived per galaxy, the real
  n(z), b = 1.34/D, 20-150 Mpc/h, full covariance across separations): with eta = 1, S/N = 9.0, 5.4, 2.5 per bin.
  With the measured eta: **4.3, 1.2, 0.1; 4.5 sigma combined against A_D = 0. Above 3: the gate passes, on the
  0.15-0.25 bin alone.**
- Why 2511.17858 forecast ~11-13 in every bin: (i) its eq. 14 normalises the size-noise term by the number of grid
  cells V/l_p^3, not the number of galaxies; at BGS densities (nbar l_p^3 = 0.5, 0.17, 0.05 for l_p = 4 Mpc/h) that
  is optimistic by 1/sqrt(nbar l_p^3) in S/N. Reproducing their normalisation here gives 12.0, 12.1, 11.2 per bin,
  their flat numbers. (ii) eta = 1. (iii) Its total S/N looks like a quadrature sum over separation bins, which are
  strongly correlated (the diagonal-only sum here is 16.6, 9.7, 4.3 against 9.0, 5.4, 2.5 with the full covariance).

**What could erase or fake the dipole** (estimates; each has a null or a model term in the plan):
- PSF size and seeing: an angular field, so with a midpoint line of sight it is even in mu up to (d/r)^2. The data
  show d kappa_hat / d ln PSFSIZE_R = +0.04, 0.00, -0.04 per bin (PSF spread 0.155 in ln), so a PSF field enters
  kappa_hat at < 0.01 rms. Nulls 3-4.
- Blending in dense regions: angular (even) for projections; 3D companions act like a size-density dependence b_s,
  which enters the dipole through the counts' own Doppler term (T_rel) and O(d/r) terms; b_s is measured from the
  monopole, the leak from mocks. Null 5.
- Size-magnitude selection: the eta above; measured, not assumed.
- Lensing convergence: < 3 % of the dipole at z < 0.45, d < 150 (Bonvin et al. Fig. 3); its own response is the
  same eta.
- Wide-angle and evolution: with the midpoint line of sight and normalisation per (d, mu) bin, the n(z)-gradient leak
  of the even part (which with a shell-total normalisation would be ~3 (d/3)(dln n/dr) xi_0, comparable to the
  signal at 20 Mpc/h) is removed; the b_s-induced O(d/r) odd terms are calibrated on mocks (F2 point 3).

## F2. Step 1: the method check (blind-review.md 1b)

Plan alone in `r156_plan.md`; the method-check prompt verbatim; a fresh context-free Opus agent; 2026-11-02.
Verdict, verbatim:

> METHOD: fixable: calibrate the Tractor size response (injection/Obiwan or HSC comparison) instead of assuming 1; derive eta by synthetic magnification through every size/flux-dependent cut (TYPE PSF boundary, sigma_lambda, FRACMASKED, fibre-magnitude, redshift success, clipping, weights) and propagate eta's statistical error jointly with xi_1 via jackknife; add the O(d/r) b_s-induced odd terms (midpoint wide-angle with unequal RSD, evolution of b, b_s, eta across the pair) to the model; run end-to-end injection-recovery on velocity-bearing BGS mocks (A_D=1 and A_D=0) and take or validate the covariance there, since 6-degree jackknife regions underestimate cosmic variance at 150 Mpc/h; add depth/extinction/stellar-density substitute-kappa nulls, convert every systematic null into an A_D contamination bound, and write out the octupole prediction

Closures (each done before any dipole was computed):

1. **Tractor size response** (`r156_acs.py`): 311 BGS galaxies (z 0.05-0.5) matched within 1" to HST/ACS F814W sizes
   in COSMOS (ZEST GIM2D half-light radius and r50; ACS PSF 0.1" vs median SHAPE_R 1.4"). OLS slope of the analysis's
   kappa_hat on the ACS size residual: 0.886 +- 0.043 (GIM2D), 0.938 +- 0.040 (r50); both ACS radii are noisy
   (their ln ratio scatters 0.13 robust), which attenuates OLS. Instrumental variables (each ACS radius instrumenting
   the other, independent ACS errors) remove that: **a = 0.966 +- 0.045 and 0.969**. Adopted a = 0.97 +- 0.05.
   Cross-check on the full sample: under "sizes regress toward the PSF", a = 1 - d kappa_hat/d ln PSF = 0.96-1.04.
   This also covers Tractor type switching (the population regression crosses type boundaries); the PSF type itself
   is 0.03 % of the sample.
2. **eta through every cut** (`r156_measure.eta_full`): flux limit (data edge derivative, m_cut 19.3-19.5;
   systematic from 19.4-19.5 and 19.2-19.5: +-0.02-0.04); sigma_lambda < 0.5: no galaxy lies within 1 % of the cut
   (Tractor sigma_lambda ~ 0.005), contribution 0.000; FRACMASKED_R: not moved by magnification at first order;
   fibre-magnitude cut: 0.01 % of galaxies within 0.5 mag of r_fib = 22.9, contribution < 0.001; winsorisation: x
   0.9966 (fraction unclipped); weights: sigma_lambda << sigma_int, contribution < 0.001; redshift success: BGS_BRIGHT
   success ~ 99 % (DESI DR1), so |d eta| <= 2 (1-p) |Delta kappa_hat| <= 0.02, taken as a +-0.02 systematic;
   Tractor response a above. eta(z) = a x 0.9966 + flux edge. The statistical error of eta is propagated by computing
   eta in each jackknife region with xi_1.
3. **b_s-induced O(d/r) odd terms**: computed end to end on the mocks (point 4): a mark kappa = D G (matter contrast,
   unit b_s, no velocity) run through the same pair counts and estimator gives the full leak per unit b_s (midpoint
   wide-angle with unequal RSD, n(z), the evolution of b across the pair); it is scaled by the measured b_s and
   subtracted, with its mock-to-mock scatter as an error. eta's evolution across the pair is in the template galaxy
   by galaxy.
4. **End-to-end mocks and the covariance**: lognormal light cones with linear velocities on the real footprint and
   n(z) (`r156_mock.py`); the Doppler mark eta(z) x kappa_v recovers A_D (target 1) and the velocity-free mark recovers
   A_D = 0. Covariance, fixed now: C = C_shuffle (300 permutations of kappa_hat within z slices of 0.01 and caps: the
   exact size-noise term on the real density field) + the analytic cosmic-variance term (2511.17858 eq. 13 with the
   measured b and eta) + the mock scatter of the b_s leak; the 200-region jackknife is a cross-check only.
5. **Nulls**: substitute marks added: GALDEPTH_R, PSFDEPTH_R, E(B-V) (from MW_TRANSMISSION_R), stellar density (Gaia
   stars, HEALPix nside 128). Every substitute null X becomes an A_D bound: dA_X = beta_X x A_D(X), beta_X the
   regression slope of kappa_hat on the standardised X. Octupole prediction written out: xi_3(d) = -eta f (H/c)
   (1 - c/(Hr)) (2f/5) nu_3(d), nu_3 = (1/2pi^2) Int dk k P(k) j3(kd) (`r156_theory.template_oct`), positive, ~10-15 %
   of |xi_1|.
6. Stated: b_s checked across separation bins; the +-50 % on C(z) is an assumption; A_D tests f sigma_8 given the
   measured b sigma_8 (b from the monopole with Planck sigma_8).

## F3. Mocks (`r156_mock.py`; 3 lognormal light cones, 512^3 in 3 Gpc/h, linear velocities, real footprint and n(z))

Pair counts per mock: 40 % of size galaxies, all mock randoms (2.35x), 12-14 min each. Doppler mark eta(z) kappa_v,
noise-free; velocity-free mark D G (unit b_s).

| z | A_D recovered (target 1), 3 mocks, +-0.035 each | b recovered vs input (20-40 Mpc/h) | velocity-free mark, per unit b_s (jackknife error) |
|---|---|---|---|
| 0.15-0.25 | 0.953, 0.971, 1.040 | 1.49, 1.43, 1.23 vs 1.49 | -35, +36, +14 (+-23-34) |
| 0.25-0.35 | 0.960, 0.936, 0.877 | 1.49, 1.71, 1.88 vs 1.57 | +71, +14, +23 (+-75-90) |
| 0.35-0.45 | 0.926, 0.995, 0.984 | 1.79, 1.54, 1.80 vs 1.65 | +108, -661, +720 (+-440-490) |

- **Caught by the mocks:** with b fitted at the pre-registered 30-80 Mpc/h the mocks returned A_D = 1.37-1.47: their
  monopole carries an offset of ~-0.005 beyond ~45 Mpc/h (ratio to input 3.0 at 15-35, falling to <1 by 55), so b came
  out at half its input. b is now fitted at 20-40 Mpc/h, where the mocks recover it; the table is with that fix.
- The jackknife errors match the mock-to-mock scatter (Doppler mark 0.035-0.04 vs 0.047; velocity-free 23-34 vs 37 in
  bin 1), so the 200-region jackknife does not underestimate the variance at these separations here.

## F4. The measurement (`r156_run.py`, `r156_refit.py`; `data/r156/run.json`, `refit.json`)

Measured per bin: b (20-40 Mpc/h) = 1.14, 1.39, 1.73 (30-80: 0.72, 2.08, 2.37); b_s = 0.0257 +- 0.0046,
0.0338 +- 0.0059, 0.0460 +- 0.0104 (flat at 0.022-0.026 for d = 15-45 in bin 1); C(z) = -3.2, 2.2, 7.5; s(z) at
r = 19.5 = 0.29 (z 0.16) to 1.04 (z 0.44); eta = 0.469, 0.191, -0.083 (stat +-0.002-0.006, syst +-0.06).

**The pre-registered covariance failed.** Size shuffles + analytic cosmic variance + a rank-1 leak term, with the
3-mock mean leak subtracted, gave chi2/dof = 57.2/12, 80.1/12, 34.6/12 (A_D = 1.19 +- 0.25, -1.02 +- 1.05). Two causes:
the subtracted "leak" is the mean of three per-mock values that scatter about zero (mock noise of up to 5e-4 at large d,
larger than the signal), and the shuffles destroy the spatially coherent b_s delta part of the sizes, whose cosmic
variance is then missing. Adopted instead (stated as a departure): the mock-validated jackknife covariance, T_rel
subtracted with its error, no mock-mean leak subtraction, the leak as an analytic systematic +-(0.1|A_D| + 0.1).

| z | A_D | stat | eta | leak | b | fit chi2 | vs zero | d >= 30 | 20-60 |
|---|---|---|---|---|---|---|---|---|---|
| 0.15-0.25 | **0.99** | 0.28 | 0.12 | 0.20 | 0.09 | 12.8/12 (p 0.39) | chi2 25.5/13 (p 0.02) | 1.31 +- 0.31 | 0.87 +- 0.29 |
| 0.25-0.35 | 0.59 | 0.95 | 0.18 | 0.16 | 0.03 | 7.9/12 | 8.3/13 | -0.01 +- 1.15 | 0.71 +- 0.99 |
| 0.35-0.45 | (eta ~ 0: GR predicts ~0) | | | | | | 13.8/13 (p 0.39) | | |

**Combined: A_D = 0.93 +- 0.35 (total), 2.7 sigma from zero, 0.2 sigma from GR.** Bin 1 alone: 3.5 sigma on
statistics, 2.5 sigma total. Sensitivities: planned b range 1.65; planned covariance 1.19 +- 0.25.

## F5. Nulls (`r156_nulls.py`; `data/r156/nulls.json`; jackknife covariance)

| null | bin 1 | bin 2 | bin 3 |
|---|---|---|---|
| randomised redshifts (pair counts redone, `--zshuffle 5`) | A 0.10 +- 0.22, chi2 15.2/13 | -0.53 +- 0.83, 12.9/13 | -6.5 +- 9.2, 13.5/13 |
| octupole vs prediction | A_oct 2.5 +- 2.7 (chi2 vs 0 11.0/13) | -5.4 +- 11.1 | |
| ln PSFSIZE_R: chi2 vs 0 / bound dA_D | 8.0 (p 0.84) / -0.015 +- 0.020 | 9.7 / -0.002 | 11.8 / 0.22 +- 0.43 |
| stars' PSF map (3.61 M Gaia-confirmed stars, G 15-20) | 8.1 (p 0.84) / 0.002 | 28.7 (p 0.007) / -0.001 | 26.6 (p 0.014) / -0.004 |
| FRACFLUX_R (blending) | 14.2 (p 0.36) / -0.10 +- 0.05 | 17.7 / -0.16 +- 0.19 | 7.5 / -0.4 +- 2.2 |
| ln GALDEPTH_R | 11.1 / 0.003 | 6.3 / -0.009 | 6.3 / 0.30 |
| ln PSFDEPTH_R | 9.5 / 0.001 | 5.7 / -0.011 | 6.3 / 0.39 |
| E(B-V) | 20.0 (p 0.09) / 0.006 | 23.3 (p 0.04) / 0.021 | 18.9 / 0.32 |
| stellar density | 23.4 (p 0.04) / 0.000 | 15.2 / 0.007 | 8.8 / -0.29 |

Splits (A_D, jackknife): NGC 1.31 +- 0.33 / SGC -0.08 +- 0.38 (bin 1; **2.8 sigma apart**), bin 2 0.26 +- 1.12 /
0.96 +- 1.43; north (BASS/MzLS) 1.93 +- 0.47 / south (DECaLS) 0.66 +- 0.32 (2.2 sigma); FRACFLUX_R < 0.05:
0.96 +- 0.31.

## F6. The finding check (blind-review.md 2)

Section alone in `r156_finding_section.md` (also appended to `finding.md`); the after-writing prompt verbatim; a fresh
context-free Opus agent with hints (Doppler lensing measurements, DESI DR1 PV 2512.03230, the bright-faint dipole,
size magnification, forecasts 2511.17858 and 1610.05946); 2026-11-02. Verdict, verbatim:

> VERDICT: novel but not sound: the amplitude depends on an η calibration no mock tests end to end, the signal appears only in NGC/BASS-MzLS imaging (2.8σ cap discrepancy), and the b_s leak and covariance are validated only by three mocks after choices made once the first results were in.

Prior work the referee found: 1610.05946, 1401.3694, 1810.12793, 2511.17858 (forecasts; "not yet detected"),
1512.03918 and 1709.07855 (bright-faint dipole, BOSS, 2.7-2.8 sigma), 1504.02662 (FP size residuals vs density),
2512.03230 and 1706.05205 (PV surveys), 0904.4702 and 2609.09151 (size bias / magnification response under flux cuts).
**Parked** on the physics points (one revision is for presentation only). What stands: the first real-data size-
number-count Doppler dipole, A_D = 0.99 +- 0.39 at z = 0.15-0.25 and 0.93 +- 0.35 combined (2.7 sigma); and the
measured flux-limit response of size magnification, eta = 0.57 at z = 0.16 -> 0 at z = 0.36 (2511.17858 assumes 1).
