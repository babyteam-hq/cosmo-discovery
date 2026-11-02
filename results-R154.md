# Results (R-154, cosmo-deep-22, 2026-11-02)

Does the SZ mass bias keep falling at z > 1? CMB-lensing M_lens/M_SZ of ACT DR6 SZ clusters in redshift bins. Every number
comes from `act_szlens.py`, which imports R-122's `erass1_act.py` and R-121's `erass1_cmblens.py`. The fixed plan is in
`plan-R154.md`. Data are under `data/` and logs under `logs/`, both gitignored. This file is the worker's own.

## Rows

### D0. Data access

- The ACT DR6 SZ cluster catalogue v1.0 (arXiv:2507.21459) is public on LAMBDA, which this cloud now reaches (200; it was 403
  for R-122). It has 10,040 clusters, 1,180 of them at z > 1, all optically confirmed. Files: `DR6_cluster-catalog_v1.0.fits`
  (26 MB) and the Legacy subsample; the Nemo selection-function products are a 7.4 GB tgz.
- Fetched from NERSC (`act_dr6_fetch.sh`): the DR6 lensing baseline, tonly, polonly and cibdeproj κ alm, the baseline mask, the
  6' SZ-cluster mask, N_L, and the 400 baseline simulations streamed four in parallel (`act_dr6_simstream.sh`). Also the eRASS1
  catalogues, for the forecast's calibration.
- Scripts: `act_dr6_fetch.sh`, `act_dr6_simstream.sh`.
- Statistic: the whole plan's data are reachable.
- What it means: the catalogue's masses are `A10_M500c`, Arnaud et al. 2010 scaling with the Eddington correction from a Tinker
  prior, the mass behind 2512.18935's 1-b. The default `M500c` uses 0.6× the A10 normalisation, chosen to match weak lensing.

### D1. What the release did at the SZ clusters

- **Source.** Qu et al. 2024 §III.6, §VI.3.1 and App. J (read 2026-11-02).
- **What it says.** The baseline map **subtracts** beam-convolved UPP tSZ templates of every cluster the lensing team's own Nemo
  run detected at SNR > 5, from the single-frequency maps before the coadd. Inpainting those clusters was the pre-unblinding
  method. It was dropped because "small inpainting artifacts at the inpainted cluster locations can correlate easily with the
  true lensing field" in simulations. The switch was motivated by "an entirely different project focused on cluster mass
  calibration".
- **The release mask.** The 6' SZ-cluster mask in `other/` is the old inpainting list.
- **Statistic.** None.
- **What it means.** R-122 called the clusters inside that mask "inpainted". For the baseline map they are template-subtracted.
  This is a loose end for R-122's finding.
  - ACT's own cluster mass calibration on this map is in progress, so the scoop risk is confirmed.

### D2. Footprint counts

- R-122's selection on the DR6 lensing mask: clean within 30', ≥ 80 % clean within 2°.
- Script: `act_szlens.py counts` (`logs/counts.log`).
- 6,263 of the 10,040 clusters are kept. 2,824 of them are inside the SZ-cluster mask, which holds 2,650 of the 3,620 kept
  clusters at SNR > 5.

| z bin | kept | in the SZ mask | median A10_M500c |
|---|---|---|---|
| 0-0.5 | 2,459 | 1,105 | 1.96e14 |
| 0.5-0.8 | 2,125 | 1,004 | 2.08e14 |
| 0.8-1.0 | 922 | 401 | 1.93e14 |
| 1.0-2.0 | 757 | 314 | 1.76e14 |
| 1.2-2.0 | 327 | 145 | — |

- What it means: there are 757 clusters at z > 1 for the stack, against about 480 in 2406.14754's MCMF z > 1 stack on the same
  map.

### F0. Feasibility gate (step 0)

- **Method.** A Fisher calculation over R-122's 9 annuli, clusters treated as independent, with the release N_L plus C_L^κκ.
  - Calibrated by the same calculation on R-122's eRASS1 cosmology sample (1,907 clusters): Fisher σ_A1h is 0.093 against the
    realised 400-sim 0.110, a factor of 1.18.
  - The 2512.18935 law used is R = 1.56((1+z)/1.45)^1.8. The bins' mean z are 0.33, 0.64, 0.89 and 1.23.
- Scripts: `act_szlens.py forecast` (`logs/forecast.log`); the separation calculation is inline in the worker's scratchpad.
- **Calibrated σ(A_1h), two-halo free, A10 template, all clusters:**

| | bin 1 | bin 2 | bin 3 | bin 4 |
|---|---|---|---|---|
| σ(A_1h) | 0.120 | 0.155 | 0.269 | 0.344 |
| dlnA/dlnM | 0.80 | 0.84 | 0.85 | 0.87 |

- **Separation, if 1-b = 0.74 everywhere:** σ(ζ) = 0.51, so ζ = 1.8 against ζ = 0 is separated at **3.5σ**. The z > 1 bin alone
  separates the extrapolated law (R = 3.4) from 1.35 at 3.0σ.
- **If the 2512.18935 law is true:** 5.3σ.
- **Outside the SZ mask only:** 2.2σ (2.0σ from the z > 1 bin alone).
- **Anchor.** 2406.14754's 4.3σ z > 1 stack, with a 5' aperture and about 480 clusters, scales to about 5.4σ for 757 clusters,
  i.e. σ(ln R) ≈ 0.21 at z > 1 and about 4σ separation.
- **What it means.** The gate passes, at 3.0-3.5σ, but only on the full sample. The clean half outside the SZ mask alone does
  not reach 3σ, so the response at the template-subtracted clusters must be calibrated, not dropped. The forecast is
  noise-dominated: the separation at z > 1 comes from the z > 1 bin's 757 clusters.

### M1. Method check (blind-review.md 1b)

- The plan (`plan-R154.md`) alone, the method-check prompt verbatim, a fresh context-free Opus agent.
- The referee's verdict, verbatim:

  > METHOD: fixable: injection-calibrate the reconstruction (NFW κ plus measured-y tSZ, with template subtraction where the
  > release did it) for the response at SZ-mask positions, tSZ leakage and estimator non-linearity per bin, instead of relying
  > on the noisier polonly difference and the confounded in/out-mask and SNR splits; separate CIB from tSZ with a
  > measured-amplitude CIB injection rather than cibdeproj−baseline; model optical-confirmation completeness in the selection
  > mock; extend the Eddington mock to 1-b ≈ 0.25, apply it self-consistently in R, vary scatter and σ8, and run the real
  > profile fit on the mock stacks; correct and budget purity and photo-z/bin migration from the catalogue's false-detection
  > rates and photo-z posteriors; calibrate σ_off per bin from SZ-optical offsets and budget it; cross-match DR5/DR6 A10 masses
  > and report ζ with η fixed at −0.18; include every z-dependent systematic in the ζ interval

- **What it means.** Every point must be closed before the data map is stacked at the clusters. The simulations are stacked
  meanwhile; they carry no cluster signal. The heaviest point is the injection calibration of the reconstruction, which needs a
  map-level pipeline the release does not provide.

## Rows I1-I5: closing the method check (before any data stack)

### I1. DR5/DR6 mass definition, miscentring calibration, photo-z

- **Masses.** DR6 A10_M500c against DR5 M500c (the mass behind 2512.18935's 1-b) for the 3,331 clusters cross-matched within
  1.5' and Δz < 0.05.
  - Median ratio: 1.030, 1.018, 1.013 and 1.017 (bins 1-4).
  - Mean ln ratio: 0.010 ± 0.004, 0.002 ± 0.004, 0.000 ± 0.007, 0.017 ± 0.011.
  - The uncorrected masses agree to 1 %.
- **Miscentring.** SZ-optical offsets, Rayleigh σ from the median offset: 0.38, 0.34, 0.32 and 0.30' (rms-based 0.49, 0.40,
  0.38 and 0.37'). The optical position sources are RM, Vis-BCG, MaDCoWS2, WaZP and others.
- **Photo-z.** Photometric fraction 0.43, 0.58, 0.78 and 0.90. Median phot σ_z 0.010, 0.021, 0.033 and 0.059. Spectroscopic n:
  2,294, 1,413, 310 and 116.
- Script: inline in the worker's scratchpad (`xm.py`); the numbers are in this row.
- **What it means.**
  - The DR6 A10 mass is the DR5 mass within 1-3 %, with no redshift trend. No correction.
  - Miscentring: σ_off = 0.5' is an upper bound, since the offsets include optical-centre errors. The systematic is the shift to
    0'.
  - Photo-z errors are well below the bin widths.

### I2. Purity (the DR6 paper, §2.3 and §3.2)

- **Source.** Noise false positives from signal-free simulations: 11 expected at q > 5.5 in the Legacy footprint. The chance
  association with optical clusters is 1-6 %. The confirmed fraction tracks 1 − F_false, with a small excess at q ≈ 4 that the
  paper attributes to spurious matches.
- **Statistic.** None.
- **What it means.** The contamination of the confirmed catalogue is a few percent, concentrated at q < 5. The budget is a
  one-sided 0-4 % (0-6 % at z > 1). Each percent of contamination lowers R by about 1.2 %.

### I3. Injection calibration of the reconstruction (method check point 1)

- Script: `act_szlens_inject.py calib/run` (`logs/inject_*.log`). The pipeline is described in the plan addendum. Paired
  realisations: 100 (√C_l^tSZ hardening, the release's choice), 60 (1.5' UPP hardening) and 150 (3' UPP hardening, two bins).
- **Gaussian-κ normalisation.** The response is 0.96-1.06 in L bins 300-2000 (8 realisations). The emulated N_0 at
  noise_scale 1.6 is about the release's tonly N_L.
- **Lensing response, √C_l^tSZ hardening.** 1.003 ± 0.036 and 0.995 ± 0.034 (z 0.33; R 1.35 and 3.4); 1.007 ± 0.039 and
  1.003 ± 0.037 (z 0.89); z 0.64 and 1.23 in `logs/inject_clhard.log`.
  - With 1.5' hardening: 1.00-1.01.
  - With 3' hardening: 0.955-0.957 at every z and R.
  - Unhardened: 0.93-0.99, falling with R, which is the estimator's non-linearity.
- **tSZ leakage, in units of the bin's one-halo κ amplitude at R = 1, temperature estimator.** Values for bins 1-4
  (z 0.33, 0.64, 0.89, 1.23):

| case | bin 1 | bin 2 | bin 3 | bin 4 |
|---|---|---|---|---|
| unhardened | −1.12 | −0.79 | −0.62 | −0.55 |
| hardened, √C_l^tSZ, UPP truth | +0.09 | +0.14 | +0.14 | +0.15 |
| same, Planck-2013 GNFW truth | +0.05 | +0.12 | +0.13 | +0.15 |
| hardened, 1.5' UPP | +0.06 | +0.01 | +0.03 | — |
| hardened, 3' UPP | +0.06 | +0.12 | — | — |
| template-subtracted, then hardened, UPP truth | −0.006 | +0.005 | +0.000 | −0.003 |
| template-subtracted, then hardened, P13 truth | 0.000 | +0.013 | +0.014 | +0.016 |

  - At 4-4.5 × 10^14, unsubtracted: ±0.18.
- **What it means.**
  - Without hardening, the temperature estimator would lose 55-110 % of the cluster κ.
  - With the release's hardening, the clusters outside the SZ mask leak +0.09 to +0.15. The template-subtracted ones leak
    ≤ 0.016.
  - The linear response is 1.00 within ±0.04, except under the 3' hardening variant (−4.5 %).
  - In the baseline map the temperature leakage is weighted by w_T ≈ 0.9 and by the outside-mask fraction (0.55).

### I4. CIB, measured and injected (method check point 2)

- Scripts: `act_szlens.py cib` (`logs/cib.log`); `act_szlens_inject.py cib` (`logs/inject_cib_clhard.log`).
- **Planck 857 GHz aperture flux per cluster** (the random null-mean of 334 ± 21 mJy removed):

| | bin 1 | bin 2 | bin 3 | bin 4 |
|---|---|---|---|---|
| 857 GHz (mJy) | 220 ± 61 | 360 ± 65 | 224 ± 99 | 336 ± 109 |
| 545/857, after the A10 tSZ at 545 | 0.32 | 0.36 | 0.44 | 0.56 |
| matching modified blackbody | β 1.5-2.0, T 20-25 K | same | same | same |
| 150 GHz flux (mJy) | 0.9-1.5 | 2.1-3.2 | 1.8-2.6 | 4.3-5.6 |
| CIB leakage, √C_l^tSZ-hardened | −0.03 to −0.05 | −0.07 to −0.10 | −0.06 to −0.08 | −0.12 to −0.14 |
| CIB leakage, template-subtracted | ≤ 0.005 | ≤ 0.005 | ≤ 0.005 | ≤ 0.005 |
| CIB leakage, unhardened | +0.13 to +0.48 | | | |

  The template-subtracted case is small because the subtraction's amplitude fit absorbs the CIB.
- **What it means.**
  - The CIB is detected in the stacks at S/N 2-6, and its SED is pinned by the measured 545/857 ratio.
  - Its leakage after hardening is negative and grows with z. It opposes the tSZ leakage, so the net leakage outside the mask
    is +0.02 to +0.07.
  - Because it is z-dependent, it enters ζ as a variant.

### I5. The selection (Eddington) mock (method check points 3, 4)

- Script: `act_szlens.py selmock` (`logs/selmock.log`).
- **B = R_fit/R_true at σ_int 0.2 (the catalogue's) and σ8 0.80:**

| 1-b | bin 1 | bin 2 | bin 3 | bin 4 |
|---|---|---|---|---|
| 0.9 | 1.07 | 1.05 | 1.07 | 1.05 |
| 0.75 | 1.045 | 1.045 | 1.045 | 1.045 |
| 0.6 | 1.03 | 1.02 | 1.01 | 1.01 |
| 0.4 | 1.00 | 0.99 | 0.98 | 0.97 |
| 0.25 | 0.95 | 0.92 | 0.83 | 0.88 |

- **Variants at 1-b = 0.75 / 0.4:**
  - σ_int 0.1: +5 % / ;
  - σ_int 0.3: −7 % to −10 % / −7 % to −13 %;
  - σ8 0.77 or 0.83: ±1.5 %;
  - optical confirmation at z > 1 (M_50 1-1.5e14): +1 % to +4 %.
- **What it means.**
  - The catalogue's Eddington-corrected A10 masses carry a residual selection bias of +5 % on R at 1-b ≈ 0.75, falling below 1
    toward 1-b ≈ 0.3. It is corrected self-consistently per bin.
  - The scatter assumption is the largest selection systematic (−7 % to +5 %, larger at low 1-b). It is nearly z-independent at
    fixed 1-b, so it mostly moves the normalisation, not ζ.

## Rows V0-V1: validation of the fit (before the data stack)

### V0. Injection into a simulated stack

- **Method.** R = 1.5 × ((1+z)/1.45)^ζ injected (1h + 2h, posterior-averaged) into simulation 6's stack, and noiseless. Script:
  `act_szlens_fit.py fit --variant simtest --nocal --trend`.
- **ζ = 0 with noise:** R per bin 1.51, 1.53, 1.46 and 2.04, ζ = 0.30 [−0.25, 0.80].
- **ζ = 1.8 with noise (simulation 18):** ζ = 1.25 [0.95, 1.65], 95 % [0.55, 2.05].
- **ζ = 1.8 noiseless:** R = 1.29, 1.87, 2.43 and 3.26 against truth 1.28, 1.87, 2.42 and 3.26. ζ = 1.85 [1.40, 2.20]; Δχ² at
  ζ = 0 is 18.1.
- **What it means.** The fit is unbiased, and σ(ζ) ≈ 0.4-0.5 as forecast. A true ζ = 1.8 would be separated from 0 at ~4σ.

## Rows R1-R6: the measurement

### R1. The measurement (the plan's statistic)

- Baseline map; all 6,263 clusters; 400-simulation covariance (joint over bins, Hartlap); method-check calibrations applied
  (response, leakage, selection B). Script: `act_szlens_fit.py fit --trend` (`logs/fit_final.log`).

| z bin | n | ⟨z⟩ | A_1h at R = 1 (2h free) | R = M_lens/M_SZ, 68 % | 95 % | 1-b |
|---|---|---|---|---|---|---|
| 0-0.5 | 2,459 | 0.33 | 1.06 ± 0.12 (8.8σ) | 0.99 [0.85, 1.14] | [0.72, 1.31] | 1.01 [0.88, 1.18] |
| 0.5-0.8 | 2,125 | 0.64 | 0.83 ± 0.14 (6.0σ) | 0.75 [0.60, 0.90] | [0.47, 1.05] | 1.34 [1.11, 1.67] |
| 0.8-1.0 | 922 | 0.89 | 0.76 ± 0.25 (3.0σ) | 0.68 [0.42, 0.96] | [0.20, 1.27] | 1.46 [1.04, 2.36] |
| 1.0-2.0 | 757 | 1.23 | 1.34 ± 0.33 (4.0σ) | 1.32 [0.93, 1.81] | [0.59, 2.32] | 0.76 [0.55, 1.07] |

- **The trend, statistical:** ζ = −0.57 [−1.40, 0.27], 95 % [−2.32, 1.01], R_0 (z = 0.45) 0.96, χ² 28.3 (36 − 6 dof). Δχ² is
  0.43 at ζ = 0 and 8.97 at ζ = 1.8.
- **What it means.**
  - No fall of 1-b with z. The extrapolated 2512.18935 law (ζ = 1.8) is 3.0σ from the best fit, statistically; R5 gives the
    systematics.
  - The normalisation, M_lens ≈ M_A10 at z < 1, lies above the weak-lensing 1-b ≈ 0.65-0.8.

### R2. Nulls

Scratchpad `nulls.py`, `logs/nulls_final.log`.

- **The simulations at the cluster positions:** mean A_1h −0.00 ± 0.006 to 0.016 per bin.
- **40 rotated baseline maps:** mean A_1h +0.03, −0.03, +0.02 and −0.07, each ±0.02-0.05. Their sd (0.10, 0.15, 0.28, 0.32)
  matches the simulations' (0.12, 0.14, 0.25, 0.33).
- **The curl reconstruction on the clusters:** A_1h −0.08 ± 0.12, 0.04 ± 0.14, 0.01 ± 0.25 and 0.09 ± 0.33 with the simulation
  covariance; −0.01, −0.04, −0.17 and 0.20 with its own 40 rotations.
- **The random stack:** +5 × 10^-4 in the central annulus, against a signal of 4-5 × 10^-2.
- **What it means.** Every null passes.

### R3. Map variants (40 rotations each; A_1h at R = 1, variant against baseline on the same covariance)

| map | bin 1 | bin 2 | bin 3 | bin 4 |
|---|---|---|---|---|
| polonly (tSZ- and CIB-free) | 0.75 ± 0.25 / 0.88 | 0.63 ± 0.31 / 0.68 | 1.45 ± 0.59 / 0.89 | 1.23 ± 0.77 / 1.21 |
| tonly | 1.07 / 1.05 | 0.82 / 0.88 | 0.87 / 1.04 | 1.25 / 1.28 |
| cibdeproj | 0.92 / 1.01 | 0.87 / 0.94 | 0.84 / 0.90 | 1.45 / 1.64 |

- **What it means.** The polarization-only map, free of every temperature foreground, agrees with the baseline in every bin.
  The temperature-only map agrees. The CIB-deprojected map runs 5-12 % lower (as in R-122), with no trend in z.

### R4. Diagnostic subsamples (listed in the plan) and matched comparisons

- Scripts: `fit --sample ... [--trend]` (`logs/fit_samples.log`, `logs/fit_split_trends.log`); scratchpad `matched.py`
  (`logs/matched_final.log`).
- **R per bin (bins 1-4) and the trend:**

| subsample | bin 1 | bin 2 | bin 3 | bin 4 | ζ | Δχ²(1.8) |
|---|---|---|---|---|---|---|
| q > 5 (3,620) | 1.28 | 0.91 | 0.40 | 1.39 | −1.41 [−2.37, −0.55] | 15.7 |
| q 4-5 (2,643) | 0.50 | 0.40 | 1.14 | 1.31 | 2.57 [0.90, 4.0] | 0.2 |
| inside the SZ mask (2,824) | 1.24 | 0.89 | 0.49 | 1.38 | −1.63 | 16.4 |
| outside (3,439) | 0.65 | 0.58 | 0.88 | 1.43 | 1.09 [−0.24, 2.41] | 0.3 |
| q > 5 outside the mask (uncorrected) | 1.22 | 0.68 | 0.29 | 1.74 | — | — |
| spectroscopic z (uncorrected) | 1.31 | 0.49 | 0.89 | 3.1 (n 59) | — | — |
| A10 mass < 3e14 (uncorrected) | 1.02 | 0.62 | 0.52 | 1.26 | — | — |
| A10 mass > 3e14 (uncorrected) | 1.14 | 1.37 | 1.94 | 4.99 (n 35) | — | — |

- **Matched to 2512.18935** (q > 5.5, z 0.15-0.43 / 0.43-0.70, uncorrected; B ≈ 1.04): 1-b = 0.66 [0.56, 0.79] and 0.96
  [0.79, 1.19], against their 0.81 (+0.08/−0.06) and 0.58 (+0.06/−0.05).
- **Matched to HSC** (q > 5, their three bins): M_lens/M_SZ = 1.43 [1.21, 1.65], 0.95 [0.76, 1.16] and 0.88 [0.65, 1.11],
  against their 1.3 ± 0.2, 1.6 ± 0.2 and 1.6 ± 0.3.
- **The selection mock split by q** (`logs/selmock_qsplit.log`): B(q 4-5) and B(q > 5) agree within ±7 % in every bin and every
  case.
- **What it means.**
  - The SZ-significant half (and the template-subtracted half, which largely overlaps it) excludes ζ = 1.8 at ~4σ.
  - The low-significance half has R ≈ 0.4-0.5 at z < 0.8, where its expected R is ~1.3. That is a 2.8σ split in bin 1 and 2σ in
    bin 2. The selection mock does not produce it, and contamination at the DR6 paper's 1-6 % chance-association rate is far too
    small to.
  - This low-z deficit pushes that half's ζ positive; it is unexplained.
  - In the high-z bins the halves agree.
  - Matched to the published selections, CMB lensing gives no decline of 1-b from z ≈ 0.3 to 0.6-0.9: 2.6σ against 2512.18935's
    second bin and 2.3σ against HSC's second bin. The low-z bins agree with both.

### R5. Systematics (`act_szlens_fit.py sys`, 400 simulations; `act_szlens_budget.py _final`; `logs/sys_final.log`)

Every z-dependent systematic is run as its own variant of the whole fit. Two-sided pairs are counted as half their difference; single variants as the full shift.

central: R [0.991 0.746 0.685 1.322], zeta -0.57 [-1.40, 0.27] 95 % [-2.32, 1.01]
| systematic | dlnR bin 1 | bin 2 | bin 3 | bin 4 | d zeta |
|---|---|---|---|---|---|
| tSZ profile (UPP vs Planck-13 truth) | 0.009 | 0.000 | 0.009 | 0.010 | 0.03 |
| tSZ mass scaling alpha | 0.000 | 0.000 | 0.009 | 0.000 | 0.01 |
| hardening profile (leakage) | 0.032 | 0.028 | 0.032 | 0.021 | 0.01 |
| CIB SED and shape | 0.009 | 0.009 | 0.018 | 0.010 | 0.02 |
| CIB flux error | 0.009 | 0.009 | 0.028 | 0.019 | 0.05 |
| response (3' hardening) | 0.055 | 0.055 | 0.046 | 0.065 | 0.01 |
| response error | 0.028 | 0.037 | 0.037 | 0.039 | 0.05 |
| selection: scatter 0.1 / 0.3 | 0.056 | 0.062 | 0.081 | 0.102 | 0.11 |
| selection: sigma8 | 0.012 | 0.012 | 0.012 | 0.000 | 0.02 |
| selection: optical confirmation | 0.000 | 0.000 | 0.000 | 0.039 | 0.06 |
| contamination 0-4 % (6 % at z > 1) | 0.048 | 0.048 | 0.048 | 0.073 | 0.04 |
| miscentring 0.5' -> 0 | 0.009 | 0.009 | 0.028 | 0.010 | 0.05 |
| concentration x0.7 / x1.3 | 0.134 | 0.185 | 0.185 | 0.226 | 0.15 |
| two-halo: inner 30' only | 0.131 | 0.102 | 0.092 | 0.186 | 0.07 |
| two-halo: amplitude fixed at 1 | 0.055 | 0.194 | 0.185 | 0.377 | 0.43 |
| mass dependence eta = -0.18 | 0.000 | 0.000 | 0.000 | 0.000 | 0.01 |
| **total (quadrature)** | 0.222 | 0.307 | 0.304 | 0.502 | 0.49 |
zeta with systematics: -0.57 [-1.53, 0.40] 68 %, [-2.57, 1.28] 95 %
zeta = 1.8 is 2.4 sigma above (upper 68 % half-width with systematics); statistical-only Delta chi2 at 1.8: 8.97
bin 1: R 0.991, systematic on ln R 0.222
bin 2: R 0.746, systematic on ln R 0.307
bin 3: R 0.685, systematic on ln R 0.304
bin 4: R 1.322, systematic on ln R 0.502

- The variants' raw outputs, ζ [68 %]:
  - two-halo amplitude fixed at 1: R = 0.94, 0.91, 0.82, 0.91, ζ −0.13 [−0.63, 0.20];
  - fixed and inner 30' only: ζ −0.22 [−0.66, 0.16];
  - c × 0.7: ζ −0.68; c × 1.3: ζ −0.39;
  - scatter 0.3: ζ −0.37;
  - every other variant: ζ −0.63 to −0.51.
- **What it means.**
  - The trend is stable: every variant gives ζ between −0.68 and −0.13. Each one's 68 % upper edge is ≤ 0.43, against 1.8.
  - The normalisation is not stable. The concentration and the two-halo treatment move ln R by 0.13-0.38. The per-bin R is
    reported with that systematic, and the claim is the trend.

### R6. The finding check (blind-review.md, step 2)

- **What was checked.** The section (`finding-R154-section.md`, now at the end of `finding.md`) alone, the after-writing prompt
  verbatim, a fresh context-free Opus agent.
- **Hints given:** CMB-lensing masses of SZ clusters, z splits and evolution, weak- and CMB-lensing calibrations of the
  evolution.
- **The referee's verdict, verbatim:**

  > VERDICT: novel but not sound: the two SZ-significance halves give mutually inconsistent ζ (−1.4 vs +2.6), split almost
  > exactly along template-subtracted versus emulation-calibrated hardened tSZ treatment, so the headline ζ averages over an
  > unexplained systematic larger than the effect tested.

- **What it means.**
  - Novelty is confirmed: no CMB-lensing ζ exists in the literature it found.
  - The miss is on a physics point, so the section is parked, not revised.
  - One post-verdict diagnostic, recorded and not used to revise: the tSZ-free polonly map gives the faint half at z < 0.8
    A_1h = −0.20 ± 0.37 (baseline 0.49 ± 0.19). The q > 5 half gives 1.26 ± 0.21 against 1.21 ± 0.12. The deficit is not
    hardened-tSZ leakage; it belongs to the faint catalogue entries. Scratchpad `polhalf.py`.
  - The next step is a property study of those clusters. If they are mis-associations, the clean (q > 5) sample's
    ζ = −1.4 [−2.4, −0.6] stands as the measurement.
