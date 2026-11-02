# Results (R-148, cosmo-deep-18, 2026-11-02)

The Lyman-alpha forest x Compton-y cross-correlation at z = 2-3.6 (DESI DR1 forests x ACT DR6 y), the proposal
`open-problems/round2/p6-lya-y.md` (trajectory check: "VERDICT: novel: no prior work answers this question"). This file is
the worker's own. **Stopped at the step-0 feasibility gate (first hour): the forecast on the CIB-deprojected map is
4 x 10^-4 sigma, and the upper limit it would give sits ~6,000x above the prediction.** Scripts: `lya_y_model.py`
(model), `lya_y_prep.py` (per-forest reduction), `lya_y_forecast.py` (noise and forecast); data under `data/lya/`,
`data/act_y/` (gitignored); logs `logs/{model,prep_lya,forecast}.log`.

## F0. Data access (all reachable; fetch times)

| what | where | size, time |
|---|---|---|
| DESI DR1 Lya deltas, `delta-lya-0-0` (picca v9, DESI Y1 KP6 baseline: rest frame 1040-1205 A, obs 3600-5772 A, 0.8 A pixels; DLA mask CNN+GP conf > 0.5, NHI > 20.3; BAL mask from the catalogue; `DELTA_BLIND`, BLINDING = desi_y1) | data.desi.lbl.gov `public/dr1/vac/dr1/lya-deltas/v1.0/delta-lya-0-0/Delta/` (1,028 HEALPix files) | 5.3 GB, 2.5 min (8 parallel curls, ~35 MB/s) |
| The CIII red-side control (`delta-ciii-0-0`, rest 1600-1850 A at the same observed 3600-5772 A) and Lyb deltas | same VAC | present, not fetched (the gate stopped first) |
| DESI DR1 Lya correlations VAC (cf, dmat, covariance) | `.../lya-correlations/v1.0/` | present, not fetched |
| ACT DR6 NILC y (standard), CIB-deprojected (beta 1.7, T 10.7 K), beam, wide GAL070 mask | portal.nersc.gov `act/dr6_nilc/ymaps_20230220/` | 7.1 GB, 5 min |
| Planck PR4 NILC y; its CIB-deprojected version | Zenodo records 7947597, 18876143 | reachable, not fetched (not needed for the decision) |
| Websky tSZ and CIB (`tsz_8192.fits` 3.2 GB, `cib_nu0145/0217.fits` 0.8 GB) | mocks.cita.utoronto.ca/data/websky/v0.0 (listing 403, files 200) | reachable, not fetched |

## F1. The model (stated before any cross-correlation)

xi_Fy(theta; z) = b_F(z) (sigma_T / m_e c^2) a <bP_e>(z) chi^-2 Int l dl/2pi J0(l theta) B_l P_m(l/chi, z): Limber, only
k_par ~ 0 enters because y is a projection (so the forest's Kaiser factor is 1).
- <bP_e>(z), proper: halo model with the Battaglia et al. 2012 AGN-feedback profile (P_e = 0.518 P_th), Tinker 2008 mass
  function, Tinker 2010 bias, M200c 10^11-10^16 M_sun, to 4 R200c (colossus, Planck 2018). **0.458, 0.442, 0.414, 0.380,
  0.360 meV cm^-3 at z = 2.1, 2.5, 2.9, 3.3, 3.5** (M > 10^10: +1 to +8 %). Check: 0.145 at z = 0, 0.29 at z = 0.5, 0.40 at
  z = 1; the z = 0 value is ~1.5-2x below the measured 0.2-0.3 (Chiang et al. 2020), so the model is not generous.
- b_F(z) = -0.1078 ((1+z)/3.33)^2.9 (the DESI DR1 Lya auto-correlation value at z_eff = 2.33, arXiv:2404.03001, as quoted
  from the paper, not refitted from the VAC; the decision below is insensitive to it by three orders of magnitude).
- P_m: CAMB halofit. B_l: ACT `ilc_beam.txt`; Planck 10' Gaussian.
- Model xi on the ACT beam (units 10^-10), theta bin centres 1.4, 2.4, 3.9, 6.1, 8.7, 12.2, 17.3, 24.5, 36.7, 52' :
  z 2.0-2.5: -4.37 -3.25 -2.46 -1.75 -1.27 -0.87 -0.55 -0.31 -0.13 -0.05; z 2.5-3.0: -3.83 ... -0.03; z 3.0-3.6: -3.27 ... -0.02.

## F2. The continuum projection removes the signal (`lya_y_prep.py`, `logs/prep_lya.log`)

The estimator of the brief, xi_b(theta) = Sum_q A_qb ybar_q(theta) / Sum_q W_qb, with A_qb = Sum_{i in b} w_i (P delta)_i and
P picca's projection (each forest's weighted mean and slope in log lambda removed, picca's cf default, which is how DESI
models its own continuum-fitting distortion). Because y is one number per sightline at each theta, only each forest's
weighted mean delta in a bin enters, and that mean is what the quasar continuum fit (a_q + b_q log lambda) absorbs. The model
template run through the same projection, forest by forest, on the real pixels and weights:

| | z 2.0-2.5 | z 2.5-3.0 | z 3.0-3.6 |
|---|---|---|---|
| forests (all DR1) | 375,966 | 184,036 | 60,015 |
| pixels (all DR1) | 146.6 M | 66.8 M | 21.4 M |
| forests / pixels in the ACT wide mask (> 0.9) | 221,872 / 86.6 M | 108,634 / 39.3 M | 35,395 / 12.7 M |
| survival of the model after projection, per theta bin | +0.0000 to +0.0001, -0.0004 at 52' | 0.0000 to +0.0003 | -0.0002 to +0.0011 |

The survival is ~10^-4 at every angle, of either sign: the smooth y template along a forest is a constant plus a slope in
log lambda to that precision. Without the projection, the raw deltas carry the same loss (the continuum fit itself absorbs
the mean; their per-forest weighted means, rms 0.018, are fit residuals, not the k_par = 0 mode).

## F3. The forecast (`lya_y_forecast.py`, `logs/forecast.log`)

Noise: the annulus-mean y around 6,000 forest positions rotated by +97 deg in RA (a null placement; no forest-y pairing
read), sigma_ybar = 34, 29, 21, 17, 14, 11, 9.5, 7.4, 5.8, 4.9 x 10^-7 (CIB-deprojected) and 31 ... 4.4 x 10^-7 (standard)
from 1.4' to 52'. Forecast noise on xi: sqrt(Sum A^2) sigma_ybar / Sum W, sightlines independent (optimistic beyond ~5',
where neighbouring annuli overlap). S/N summed in quadrature over the ten theta bins (an upper bound; the bins are correlated):

| case | map | z 2.0-2.5 | z 2.5-3.0 | z 3.0-3.6 |
|---|---|---|---|---|
| **DESI DR1 deltas as released (projected)** | **ACT CIB-deprojected** | **0.00035** | **0.00007** | **0.00011** |
| same | ACT standard | 0.00040 | 0.00008 | 0.00013 |
| hypothetical: continuum known exactly, no projection | ACT CIB-deprojected | 5.96 | 2.11 | 0.87 |
| hypothetical: continuum known to 5 % per sightline | ACT CIB-deprojected | 1.74 | 1.03 | 0.49 |
| hypothetical: continuum known to 10 % (Croft+2005's sigma(tau_eff) = 0.1) | ACT CIB-deprojected | 0.90 | 0.57 | 0.28 |

**Gate.** Below 2.5 sigma on the CIB-deprojected map in every bin, by four orders of magnitude. As an upper limit: the 95 %
limit on the amplitude relative to the prediction would be ~2/0.00035 = **~5,700x the prediction** at z = 2.0-2.5 (~3 x 10^4 in
the other bins); the brief's condition (a limit below the predictions by an informative factor) fails. **Stopped and
reported per step 0.** Planck (10' beam, higher noise) and the CIB-leakage forecast were not computed: the CIB leakage is a
k_par ~ 0 correlation too and is projected away by the same factor, so neither can change the decision. What survives the
projection is anything tied to the individual sightline rather than to the large-scale field (the quasar host's own tSZ,
dust and radio emission correlating with its continuum-fit residuals): after projection a measured signal would be
host systematics, not forest pressure.

## F4. What it means, and what would make it feasible

- The proposal's precision estimate ("~100x the forests of Croft et al.; a few-sigma detection plausible") missed that the
  forest x projected-field two-point function lives entirely in each sightline's mean transmission, which any per-quasar
  continuum fit removes. Croft, Banday & Hernquist 2005 (astro-ph/0512380, read in full) used exactly that mean (tau_eff over
  Lya-Lyb) with a polynomial continuum and assumed the continuum error is noise (sigma(tau_eff) = 0.1), not that it erases the
  mode; Doux et al. 2016 (arXiv:1607.03625) say it outright for CMB lensing ("the mean Lyman-alpha transmission is strongly
  affected by continuum fitting ... For this reason, we instead correlate the CMB lensing convergence with the small-scale
  Lyman-alpha power spectrum") and measure the squeezed bispectrum instead.
- Even a continuum known to 5 % per sightline (a red-side PCA prediction at its best) gives <= 1.7 sigma (upper bound) at
  z 2-2.5 with all of DESI DR1 x ACT DR6, and re-extracting the forests from the DR1 spectra is a multi-TB job.
- The feasible version is the bispectrum: the small-scale forest power (or delta^2) on each sightline against y, as Doux et
  al. did for lensing. It responds to the large-scale pressure through the forest's power response, needs a response model
  (simulations), and is a different proposal: it needs its own trajectory check before any work.
