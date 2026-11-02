# Results (R-150, cosmo-deep-19, 2026-11-02)

The circumgalactic dust reddening profile around DESI DR1 galaxies at z ~ 0.4-1.6, from the colours of DESI DR1 z > 2
quasars behind them (`open-problems/round2/p7-dust-z1.md`; trajectory check "partly novel": a galaxy-selected stacked
reddening profile against impact parameter at z ~ 1-1.5, and its change since z ~ 0.3, is the open part). This file is
the worker's own. Template: Ménard, Scranton, Fukugita & Richards 2010 (arXiv:0902.4240), read in full first: the
estimator is the pair-weighted mean quasar colour excess at separation theta minus its random-position value
(their <E_ab(phi) delta_g(phi + theta)>), errors by resampling, nulls with stars as backgrounds and a Galactic-dust
split, magnification separated as the achromatic part of the per-band brightness shifts. Their result, used as the
reference: <E(g-i)>(theta) = (1.4 +- 0.1) x 10^-3 (theta/1')^(-0.84 +- 0.05) for SDSS i < 21 galaxies at z_eff = 0.36.

Scripts: `dust_z1_prep.py` (slim catalogues), `dust_z1_common.py` (everything fixed in M1), `dust_z1_feas.py` (step 0),
`dust_z1_stack.py` (the measurement). Data under `data/dust/` (gitignored); logs `logs/dust_*.log`.

## F0. Data access (fetch times; `data/` started empty)

| what | where | size, time |
|---|---|---|
| DESI DR1 LSS clustering catalogues v1.5: ELG_LOPnotqso and LRG, NGC + SGC, data and random file 0 | data.desi.lbl.gov `public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5/` | 1.3 GB |
| DESI DR1 QSO full catalogue (all main-survey QSO targets: Z_not4clus, Legacy Surveys DR9 g r z W1 W2 fluxes, MW transmission, Redrock COEFF) | same, `QSO_full.dat.fits` | 1.95 GB |
| DESI DR1 Mg II absorber VAC v1.0 | `public/dr1/vac/dr1/mgii-absorber/v1.0/MgII-Absorbers-DR1.fits` | 66 MB |
| all of the above in parallel | | **88 s** (04:40:06-04:41:34 UTC) |
| desitarget (GitHub, for the DR9 QSO random forest `rf_model_dr9_final.npz` and `myRF.py`) | github.com/desihub/desitarget | 50 MB, 20 s |

## F1. Counts (`dust_z1_feas.py`, `logs/dust_feas.log`)

Background quasars: main-survey QSO targets (DESI_TARGET bit 2) with a good DR1 redshift, z > 2, 16.5 < r < 23, W1 at > 3
sigma: **447,716** (N 84,295, S 363,421; median z 2.41, median r 21.47). Pairs are kept only where z_q > z_f + 0.5.
Pairs at 20 kpc - 3 Mpc (proper, at z_f):

| foreground bin | foregrounds | pairs | quasars with a pair | pairs < 131 kpc | pairs < 458 kpc |
|---|---|---|---|---|---|
| ELG 0.8-1.1 | 1,016,365 | 1,849,262 | 376,292 | 694 | 18,311 |
| ELG 1.1-1.3 | 654,829 | 1,081,861 | 326,667 | 417 | 10,306 |
| ELG 1.3-1.6 | 760,878 | 1,191,746 | 336,853 | 436 | 11,034 |
| LRG 0.4-0.6 | 506,911 | 1,474,929 | 386,187 | 538 | 18,386 |
| LRG 0.6-0.8 | 771,894 | 1,646,627 | 393,090 | 643 | 18,608 |
| LRG 0.8-1.1 | 859,822 | 1,537,867 | 390,949 | 579 | 16,217 |

The DR1 spectroscopic samples are sparse (~250 deg^-2 per tracer against Ménard's 1.8 arcmin^-2 photometric galaxies), so
the pairs inside a halo are few: ~10^3 within 130 kpc per bin. The signal-to-noise lives at 0.5-3 Mpc.

## F2. Per-quasar estimator noise (from the data, no pairing)

Colour residuals (each quasar's g-r, r-z, z-W1 minus the median of quasars in its 0.02 redshift bin, imaging system and
r-magnitude quartile; clipped at +-1 mag): robust sigma 0.171, 0.177, 0.370 mag (std 0.209, 0.214, 0.397). The GLS
projection onto the SMC reddening vector at z_f (E(colour)/A_V,rest = 1.04, 0.62, 1.28 for g-r, r-z, z-W1 at z_f = 1.2)
gives per-quasar sigma(A_V,rest) = **0.17 / 0.18 / 0.23 mag** for ELG z_f = 1.45 / 1.2 / 0.95, and 0.23 / 0.27 / 0.30 for
LRG z_f = 0.95 / 0.7 / 0.5 (lower z_f puts the bands at redder rest wavelengths, where the curve is flatter). It is
intrinsic quasar colour scatter, not photometric noise (brightest to faintest quartile: 0.15 to 0.20 at z_f = 1.2).

## F3. Forecast (pair counts from positions, per-quasar noise from F2; quasars independent, foreground clustering included
through each quasar's summed pair weight)

The reference profile: Ménard's eq. 28 converted to rest-frame A_V at the galaxy, A_V,rest = E(g-i)_obs / [k(4686/1.36) -
k(7480/1.36)], theta -> proper r_p at z = 0.36: A_V,rest(100 kpc) = **5.8 x 10^-3 (MW, O'Donnell R_V = 3.1)** or 4.6 x 10^-3
(SMC bar), slope -0.84 (their 100 h^-1 kpc = 148 kpc: 4.1 x 10^-3, their eq. 30, recovered). (a) no evolution: the same
rest-frame A_V at the same proper r_p; (b) dust-to-gas ~ gas-phase metallicity, -0.2 dex from z = 0.36 to z ~ 1.2 at fixed
stellar mass (mass-metallicity evolution), gas columns at fixed proper r_p unchanged: 0.63 x (a).

| bin | S/N (a) | S/N (b) |
|---|---|---|
| ELG 0.8-1.1 / 1.1-1.3 / 1.3-1.6 | 2.3 / 2.3 / 2.6 | 1.4 / 1.5 / 1.6 |
| LRG 0.4-0.6 / 0.6-0.8 / 0.8-1.1 | 1.7 / 2.0 / 2.3 | 1.1 / 1.2 / 1.4 |
| **ELG combined** | **4.1** | **2.6** |
| **all six** | **5.4** | **3.4** |

(A first pass of this table evaluated the reference at r_p x h and read 7.5 / 4.7; fixed before anything else was run,
log `logs/dust_feas_bugged.log` kept.)

## F4. Selection loss (DESI's own targeting rerun on reddened quasars)

The DR9 QSO random forest (desitarget `rf_model_dr9_final.npz`, 500 trees; the north/south/DES probability thresholds and
the W1, W2 < 22.3 cuts of `isQSO_randomforest`), run on 40,000 random z > 2 quasars from F1: all 40,000 pass at their
catalogue photometry (the selection is reproduced). Each quasar's dereddened fluxes are then reddened by SMC dust at z_f
and the forest rerun (`logs/dust_selection.log`):

| z_f | quasars lost per mag of A_V,rest (A_V = 0.01 / 0.03 / 0.1) | measured / true mean A_V among the survivors |
|---|---|---|
| 0.7 | 0.85 / - / - | 0.900 |
| 0.95 | 1.02 / 0.96 / 0.91 | 0.881 / 0.897 / 0.913 |
| 1.2 | 1.20 / 1.11 / 1.09 | 0.889 / 0.901 / 0.917 |
| 1.45 | 1.34 / 1.23 / 1.21 | 0.885 / 0.901 / 0.918 |

At the expected A_V (10^-3 - 10^-2) the targeting loses ~1 % of quasars per 0.01 mag and the measured mean reddening is
**0.88-0.89 of the true one**: a 11-12 % attenuation, linear in A_V at the level that matters, smaller than the signal and
modelled. (It counts losses only: objects that a reddening would carry into the selection are not in this catalogue. For
reddening, entering means a non-quasar-like object becoming quasar-like, which the RF's training makes rare.)

**Gate: passed.** Forecast 5.4 sigma (a) and 3.4 sigma (b) over all six bins, above 3; the selection bias is 11 % of the
signal and modelled. The margin is thin: each bin alone is a ~2 sigma profile, so the claim will be the combined
amplitude ratio to z ~ 0.36 with its interval, per-bin values with theirs.

## M1. Fixed before any stack (2026-11-02 05:15 UTC; no foreground-quasar colour pairing has been read)

- **Foregrounds.** DESI DR1 LSS v1.5 clustering catalogues (NGC + SGC), weight WEIGHT. ELG_LOPnotqso in z_f 0.8-1.1, 1.1-1.3,
  1.3-1.6; LRG in 0.4-0.6, 0.6-0.8, 0.8-1.1 (the LRG 0.4-0.6 bin bridges to Ménard's z ~ 0.36). Randoms: random file 0 of
  each (carrying the data redshifts), weight WEIGHT.
- **Backgrounds.** As F1: main QSO targets, good z, z > 2, 16.5 < r < 23, W1 > 3 sigma; per pair z_q > z_f + 0.5.
- **Estimators.** Primary: the per-quasar A_V,rest of F2 (GLS over g-r, r-z, z-W1 residuals, SMC-bar curve at the bin's
  central z_f, covariance per r-quartile from the data), each quasar weighted 1/sigma^2 of its quartile. Secondary: the same
  with the MW O'Donnell R_V = 3.1 curve; each colour alone converted to A_V. The spectroscopic continuum slope (Redrock QSO
  PCA reconstruction from COEFF against the median at the quasar's redshift) is secondary and is run only if the template
  reconstruction can be validated; if not, it is reported as not run.
- **Impact-parameter bins.** 8 log bins, 20-3000 proper kpc (edges 20, 37, 70, 131, 245, 458, 857, 1604, 3000), r_p from the
  angle and D_A(z_f) (Planck 2018).
- **Statistic.** A(r_p) = [pair-weighted mean A_hat over foreground pairs] - [the same over random-foreground pairs]
  (Ménard's <E delta_g>), weights w_g w_q, then divided by the selection response of F4 (0.885 at z_f >= 0.95, 0.90 at
  z_f < 0.95; +-0.02 systematic, the spread in F4).
- **Mg II.** Pairs whose quasar has a DR1 Mg II absorber within |dv| < 1000 km/s of z_f are flagged; every number is
  reported with all pairs (primary) and with flagged pairs removed.
- **Magnification.** From the same pairs: per-band mean magnitude shifts (g, r, z, W1, same residualisation as the colours,
  pairs minus randoms) fitted in each r_p bin as dm_b = -m_mu + A_V k_b(z_f), m_mu the achromatic (lensing) term; reported
  alongside, a cross-check on the colour A_V, not used to correct it (colours are achromatic-free; the colour-magnitude trend
  is removed by the r-quartile residualisation).
- **Errors.** Jackknife over 100 sky regions (equal-quasar-count groups of HEALPix nside-16 pixels); the full 8x8
  covariance per bin, Hartlap-corrected.
- **Fit.** A_V(r_p) = A_150 (r_p / 150 kpc)^-gamma over 20-3000 kpc; (i) gamma fixed at Ménard's 0.84, the amplitude ratio
  rho = A_150 / A_150,Ménard with the reference converted under the same curve as the estimator (SMC: 3.5 x 10^-3 at 150 kpc)
  and Ménard's 7 % amplitude error added in quadrature; (ii) gamma free. Per bin and combined (ELG, LRG, all), the combined
  rho by inverse-variance weighting over bins with the jackknife cross-covariance between bins.
- **The claim** is the table: per bin, A_V at r_p = 50, 150, 500, 1500 kpc (the power-law fit) and in the 8 bins, A_150 and
  gamma with intervals, and rho with its interval. No statement about outflows or dust physics beyond what rho says.
- **Nulls** (each in the same statistic, chi^2 against zero with the jackknife covariance): (1) foregrounds behind the
  quasar: ELGs at z_f > z_q + 0.1 around quasars at 0.6 < z_q < 1.5 (same residualisation; tests blending and photometric
  systematics, which do not know the redshift order); (2) random foreground positions: random file 0 as "foregrounds"
  against a second, disjoint half of it as the randoms; (3) foreground positions rotated by +2 deg in RA (footprint kept,
  pairing destroyed); (4) stars as backgrounds: QSO targets whose DESI spectrum is a star (SPECTYPE STAR, same targeting and
  photometry), colours residualised by imaging system and r quartile, the same estimator; (5) Galactic dust: the quasars
  split at the median SFD E(B-V), the two halves' rho compared.

## D1. The stack (`dust_z1_stack.py main`, `dust_z1_table.py`; `logs/dust_stack_main.log`, `logs/dust_table.log`)

Primary estimator (SMC, GLS over three colours, randoms subtracted, divided by the F4 response), jackknife errors (100
regions). rho is the fixed-slope (0.84) amplitude at 150 kpc over Ménard's under the same curve (3.31 x 10^-3), Ménard's 7 %
included. "chi2 vs 0" is the 8-bin detection test.

| bin | A_150 (10^-3) | rho | chi2 vs 0 (p) | free slope: A_150 (10^-3), gamma |
|---|---|---|---|---|
| ELG 0.8-1.1 | -2.95 +- 2.11 | -0.89 +- 0.64 | 8.9/8 (0.35) | -4.5 [-7.7, -1.2], 1.09 [0.59, 1.65] |
| ELG 1.1-1.3 | -1.58 +- 1.92 | -0.48 +- 0.58 | 3.4/8 (0.91) | unconstrained |
| ELG 1.3-1.6 | 2.13 +- 1.51 | 0.64 +- 0.46 | 15.6/8 (0.049) | 3.4 [0.6, 6.2], 1.13 [0.37, 1.67] |
| LRG 0.4-0.6 | 7.70 +- 2.41 | 2.33 +- 0.75 | 17.2/8 (0.028) | 6.2 [2.3, 10.6], 0.72 [0.32, 0.97] |
| LRG 0.6-0.8 | 9.88 +- 2.24 | 2.99 +- 0.71 | 30.3/8 (2e-4) | 16.8 [13.4, 20.1], 1.29 [1.13, 1.48] |
| LRG 0.8-1.1 | 8.79 +- 2.14 | 2.66 +- 0.68 | 20.9/8 (0.008) | 7.7 [4.5, 11.1], 0.75 [0.51, 0.96] |
| **ELG combined** | 0.13 +- 0.96 | **0.04 +- 0.29** (95 % < 0.61) | | |
| **LRG combined** | 8.8 +- 1.6 | **2.65 +- 0.48** | | |
| all six | | 0.79 +- 0.26 | | |

Combined values use the jackknife cross-correlation between bins (|r| <= 0.27). Fixed-slope A_V,rest at 50 / 150 / 500 /
1500 kpc and the eight bin values are in `logs/dust_table.log`.

Secondary estimators, combined rho (ELG / LRG): MW O'Donnell curve 0.12 +- 0.25 / 2.47 +- 0.44; Mg II-flagged pairs removed
-0.15 +- 0.27 / 2.61 +- 0.47 (Mg II flags change nothing beyond noise); g-r alone -0.31 +- 0.34 / 2.51 +- 0.65; r-z alone
0.00 +- 0.54 / 2.13 +- 0.58; z-W1 alone 0.58 +- 0.50 / 3.71 +- 0.89; **Redrock spectroscopic continuum slope** (6000-6800 vs
8800-9600 A of the model, validated: its residual correlates 0.63 with the photometric r-z residual) -0.52 +- 0.44 / 1.71 +-
0.58. Galactic-dust split at the median SFD E(B-V): low half -0.28 +- 0.41 / 2.64 +- 0.58, high half -0.17 +- 0.38 / 2.08
+- 0.66 (consistent).

Magnification (M1 cross-check; per-band brightness, dm_b = -m_mu + A_V k_b with the 0.84 profile): achromatic m_mu(150 kpc)
= 0.010 +- 0.005, 0.003 +- 0.006, 0.013 +- 0.006 mag (ELG) and 0.032 +- 0.006, 0.030 +- 0.005, 0.033 +- 0.006 mag (LRG):
the quasars behind LRGs are brighter by ~3 % at 150 kpc, ~3x the ELG value, as the halo masses would have it. A_V from the
brightness alone: ELG (0.5, -2.4, -0.6) +- 1.5-1.8 x 10^-3, LRG (3.4, 4.8, 5.0) +- 2.1-2.4 x 10^-3, the same sign pattern.

## D2. Post-M1 diagnostic: magnification leaking through the colour-luminosity relation (`dust_z1_stack.py quartiles`)

Not in M1; added after D1 because the per-quasar A_hat rises toward fainter r inside each r-quartile (dA/dr = 0.03-0.05 in
the brightest quartile to 0.12-0.29 in the faintest), so a lensed quasar, intrinsically fainter at a given observed r, looks
redder by dA/dr x (its induced magnitude shift). For the LRGs with m_mu ~ 0.03 mag this could be a large part of the signal.
Dust reddens every quartile equally; the leak scales with each quartile's slope. Per quartile fixed-shape amplitudes, fitted
as alpha_k = a_dust + s_k d (4x4 jackknife covariance):

| bin | quartile amplitudes (10^-3), bright to faint | a_dust -> rho | d (mag) |
|---|---|---|---|
| ELG 0.8-1.1 | -8.5, -3.0, -3.3, 0.7 (+- 3.2-3.8) | -2.91 +- 1.22 | 0.054 +- 0.033 |
| ELG 1.1-1.3 | -0.2, -6.8, 0.2, -1.4 | (see log) | |
| ELG 1.3-1.6 | 4.2, -2.2, -0.4, 4.3 | 0.34 +- 1.01 | 0.000 +- 0.047 |
| LRG 0.4-0.6 | 3.8, 6.6, -0.3, 10.7 | 0.97 +- 1.51 | 0.011 +- 0.031 |
| LRG 0.6-0.8 | 12.1, 5.4, 13.1, 3.4 | 4.24 +- 1.26 | -0.035 +- 0.026 |
| LRG 0.8-1.1 | 6.9, 10.0, 10.2, 7.5 | 2.42 +- 1.13 | 0.004 +- 0.035 |

The LRG quartile amplitudes are flat while the slopes rise ~4x: no leak is seen (d consistent with 0 in every LRG bin), and
the leak-marginalised LRG rho is **2.70 +- 0.73** (bins combined as independent). ELG: -0.94 +- 0.64. The leak does not
carry the LRG signal; marginalising it costs ~50 % in the error.

## D3. Nulls (`dust_z1_stack.py nulls`; `logs/dust_stack_nulls.log`)

chi2 of the 8 bins against zero with the jackknife covariance; the fixed-shape amplitude in units of the reference.

| null | chi2 (p) | amplitude / ref | r_p > 131 kpc |
|---|---|---|---|
| (1) ELGs behind the quasar (489,325 quasars at 0.6-1.5) | 9.4/8 (0.31) | -0.20 +- 0.31 | -0.45 +- 0.35 |
| (2) random positions as foregrounds, ELG | 11.6/8 (0.17) | -0.39 +- 0.22 | |
| (2) the same, LRG | 5.5/8 (0.71) | -0.20 +- 0.32 | |
| (3) foregrounds rotated +2 deg in RA, ELG | 5.7/8 (0.68) | 0.14 +- 0.26 | |
| (3) the same, LRG | 6.2/8 (0.63) | 0.49 +- 0.36 | 0.83 +- 0.44 |
| (4) stars as backgrounds (58,458 QSO targets with a stellar spectrum), ELG | 10.5/8 (0.23) | -0.77 +- 0.61 | -0.37 +- 0.67 |
| (4) the same, LRG | 11.1/8 (0.20) | 1.80 +- 0.74 | 1.70 +- 0.91 |
| (5) Galactic dust split | see D1 | consistent | |

All pass the pre-set chi2 test. The two LRG nulls with a sign (stars 1.8 +- 0.7, rotated 0.5-0.8 +- 0.4) both lean positive at
~2 sigma: an LRG-specific systematic of up to ~1 in rho units (LRG density correlated with something that reddens any
background) is not excluded and is carried as a systematic on the LRG amplitude. The ELG nulls are clean (all |amp| < 0.8 at
< 1.8 sigma); the behind-the-quasar null, which tests blending and imaging systematics directly, is -0.20 +- 0.31.

Run times: stack 7 min, nulls 8 min, slope 7 min, quartiles 7 min (4 cores).

## D4. Finding check (2026-11-02, fresh context-free Opus, after-writing prompt verbatim, the section alone; hints: 0902.4240, 1411.3333, 0706.0898, 1204.1978, 2503.04098, 2606.24159, anyone measuring halo dust reddening at z > 0.5)

VERDICT: novel but not sound: the LRG detection and the ELG–LRG contrast are not separated from an LRG-tied systematic (stellar null 1.8 ± 0.7, strongest in z−W1) and a magnification leak of comparable size, and the 0.25–3 Mpc scales cannot show halo dust.

Novelty held (no galaxy-selected quasar reddening at z_f > 0.6; 2405.16484, chromatic CMASS magnification at z ~ 0.55, overlaps
the LRG 0.4-0.6 bin). The referee calls the ELG limit "more robust". The objection is physics (the LRG systematic), so the
result is **parked** per the brief, without a revision. What stands: ELG rho = 0.04 +- 0.29 (< 0.61, 95 %) at z 0.8-1.6. The
path to a pass: claim the ELG limit alone, run the LRG behind-the-quasar null, the stellar null per colour with Gaia/Legacy point
sources, quote Ménard's error at 25 %, and compare amplitudes after scaling out clustering bias (the S/N is at 0.25-3 Mpc).

# R-152 (cosmo-deep-19b, 2026-11-02, cloud worker cw-159): the ELG limit alone, its two systematics, matched comoving scales

The one revision of R-150 (brief: the headline becomes the ELG result alone, stated as large-scale reddening; the stellar
null and the magnification leak computed for the ELG sample and propagated; the LRG rows moved to "not separated from
systematics"; the limit compared at matched comoving scales with Ménard's profile and with the Mg II dust budget of
arXiv:1204.1978). Script `dust_z1_elg.py` (prep, prepstars, main, nowise, stars, table); log `logs/r152_table.log`. Nothing in
R-150's M1 estimator changed. What R-152 added is labelled post-M1 below.

## E0. Data (refetched; `data/` is gitignored and the container was new)

ELG_LOPnotqso LSS v1.5 (data and random file 0, NGC + SGC) and `QSO_full.dat.fits`, as F0: ~3.2 GB in ~3 min. New:
`zpix-main-dark.fits` (11.3 GB; the last 6.3 GB as 8 parallel byte ranges, ~1 min) and `zpix-main-backup.fits` (1.6 GB), from
`public/dr1/spectro/redux/iron/zcatalog/v1/`. Backup-program stars carry Gaia photometry only (no Legacy Surveys g, r, z,
W1) and are unusable. The dark program gives **327,205 PSF stars** (SPECTYPE STAR, ZWARN 0, ZCAT_PRIMARY), and 286,333 of
them pass the quasars' photometric cuts. Galactic transmission for them is computed from EBV with the coefficients read
off QSO_full (g 3.214, r 2.165, z 1.211, W1 0.184; the same N and S).

## E1. The ELG stack, reproduced, and the ratio at matched comoving scales

The three ELG bins reproduce D1 to the last digit (A_150 = −2.95 ± 2.11, −1.58 ± 1.92, 2.13 ± 1.51 × 10⁻³). R-150's combination
is also reproduced (ρ = 0.045 ± 0.288). R-150 took its cross-bin correlations from an unweighted least-squares amplitude
per jackknife sample. Here they come from the jackknife of the quoted GLS amplitude itself, the consistent form, which gives
**ρ = −0.09 ± 0.26 (proper r_p)**. The difference is half a sigma. Ménard's amplitude error is now 25 % (arXiv:2503.04098), which
changes nothing at ρ ≈ 0.

**Matched comoving scales (post-M1).** Both profiles are r^−0.84 power laws. At a comoving separation R the ratio is
ρ_c = ρ × [(1 + z_eff)/(1 + 0.36)]^0.84, with z_eff the pair- and GLS-weighted foreground redshift of each bin.

| bin | z_eff | ρ (proper) | factor | ρ_c | information share r_p < 131 / 131-458 / > 458 kpc |
|---|---|---|---|---|---|
| ELG 0.8-1.1 | 0.951 | −0.89 ± 0.64 | 1.354 | −1.21 ± 0.87 | 0.09 / 0.14 / 0.77 |
| ELG 1.1-1.3 | 1.198 | −0.48 ± 0.58 | 1.497 | −0.72 ± 0.87 | 0.09 / 0.06 / 0.85 |
| ELG 1.3-1.6 | 1.430 | 0.64 ± 0.46 | 1.628 | 1.04 ± 0.75 | 0.06 / 0.19 / 0.76 |
| **combined** (z_eff 1.22) | | −0.09 ± 0.26 | | **−0.19 ± 0.40** | |

95 % upper limits on ρ_c: 0.47 (Gaussian), **0.67 with a flat ρ_c ≥ 0 prior**. The sensitivity (1.645σ about zero) is 0.66.
More than three quarters of the information comes from r_p > 458 proper kpc (> 1 comoving Mpc), so this is two-halo-dominated
reddening.

**Halo scale.** R_200c of a 10^12 M_⊙ halo at z = 1.2 is 133 proper kpc (91 kpc at 10^11.5). The three bins at 20-131 kpc
hold 1,546 pairs. With the fixed shape they give ρ = −0.62 ± 0.99. A 3σ detection of ρ = 1 inside 131 kpc needs ~1.4 × 10⁴
pairs there, and of ρ = 0.3 ~1.5 × 10⁵: 10-100 times DR1's spectroscopic density. Halo dust at z ~ 1 is out of reach of this
sample.

**WISE-free estimator (post-M1, from the saved pair weights; `nowise`).** GLS over g−r and r−z only: ρ = −0.24 ± 0.31,
**ρ_c = −0.43 ± 0.47 (< 0.68)**. Single colours, ρ_c: g−r −0.46 ± 0.50, r−z −0.26 ± 0.90, z−W1 0.67 ± 0.76.

## E2. Magnification for the ELG sample, and the leak it implies (post-M1)

The achromatic shift from the per-band brightness fit (D1's model), with jackknife errors on the shift itself: m_μ(150 kpc)
= 0.0101 ± 0.0040, 0.0032 ± 0.0052, 0.0132 ± 0.0045 mag brighter. The within-quartile colour-magnitude slope of Â, weighted
over each bin's pairs, is 0.100, 0.063 and 0.064 mag of A_V per mag. The predicted leak is ΔÂ = slope × m_μ, divided by the
targeting response. A lensed quasar is intrinsically fainter at its observed r, so it looks redder.

| bin | leak in ρ (proper) |
|---|---|
| ELG 0.8-1.1 / 1.1-1.3 / 1.3-1.6 | +0.35 ± 0.14 / +0.07 ± 0.11 / +0.29 ± 0.10 |
| **combined** | **+0.23 ± 0.07 (comoving +0.34 ± 0.10)** |

The leak is positive: it adds apparent reddening, so subtracting it can only lower the limit (leak-subtracted ρ_c =
−0.49 ± 0.44, < 0.60). The shift is 0.34 in ρ_c against an error of 0.40, **less than the error**. The headline keeps the leak
unsubtracted, the conservative choice. R-150's quartile marginalisation (D2) gave −0.94 ± 0.64 for the ELGs, also lower.

## E3. Stars as backgrounds around the ELGs (0.8 < z_f < 1.6, z_f = 1.2 response; post-M1 except the first row)

Each star sample is residualised in cells of imaging system × r quartile × Redrock stellar subtype. The subtype is
spectroscopic, so it removes the stellar locus without touching a tested colour. "Own" is R-150's form (the stars' own GLS
estimator). "Through the quasar projection" puts each colour's mean shift at ELG positions through the quasars'
weight-averaged GLS projection vector: the systematic a common per-band photometric offset would put into the quasars' Â.
Amplitudes are in units of the reference (ρ, proper).

| star sample | N | median r | own: χ² (p), amp | per colour as A_V: g−r / r−z / z−W1 | through the quasar projection |
|---|---|---|---|---|---|
| R-150's: QSO targets with a stellar spectrum (QSO_full) | 58,458 | 20.97 | 10.5/8 (0.23), −0.77 ± 0.61 | −1.19 ± 0.65 / +0.41 ± 0.90 / −2.04 ± 1.19 | −1.28 ± 0.66 |
| all DR1 dark-program PSF stars | 286,333 | 18.48 | 27.7/8 (0.0005), −0.29 ± 0.09 | −1.17 ± 0.22 / −0.75 ± 0.32 / −2.05 ± 0.67 | −1.23 ± 0.28 |
| of which QSO-target stars (DESI_TARGET bit 2) | 78,737 | 20.93 | 27.8/8, −1.74 ± 0.43 | −1.51 ± 0.47 / +0.34 ± 0.71 / **−12.4 ± 1.3** | −3.81 ± 0.54 |
| LRG-target stars (bit 0 only) | 21,008 | 21.84 | 32.2/8, −6.65 ± 1.50 | +5.3 ± 2.1 / +13.8 ± 4.2 / **−28.1 ± 3.2** | −2.44 ± 1.60 |
| ELG-target stars (bit 1 only) | 847 | 21.05 | 8.8/7, +4.0 ± 3.9 | (noise) | +2.1 ± 5.8 |
| **not selected by the QSO, ELG or LRG cuts** (209,103 MWS_ANY / 205,251 STD_FAINT / 104,426 STD_BRIGHT bits: mostly spectrophotometric standards, chosen from Gaia and g, r, z, no WISE) | 185,625 | 17.98 | **7.8/8 (0.45), −0.02 ± 0.09** | −0.09 ± 0.09 / −0.23 ± 0.12 / +0.33 ± 0.24 | **−0.03 ± 0.10** |

The null fails only for stars that entered DESI through the W1-based QSO and LRG target cuts, and it fails mostly in z−W1,
the colour those cuts lean on. The profile is flat from 0.36 to 2.4 Mpc (QSO-target stars: −3.9, −3.6, −3.2, −1.5 ×
10⁻³ in the outer four bins), not r^−0.84. That is the signature of a selection effect (a star near the QSO-selection boundary
is admitted with a W1 that depends on the local WISE depth, which the ELG density weakly traces), not of a photometric offset
around galaxies.

The stars that no W1 cut chose are the photometric null: −0.03 ± 0.10, three times tighter than the measurement, χ² 7.8/8.

Two checks on the quasar side show the quasars do not inherit the target-star offset:
- R-150's null 1 (z 0.6-1.5 quasars, same targeting and photometry, ELGs behind them) is −0.20 ± 0.31, against −3.8 ± 0.5 if
  they did.
- The WISE-free estimator (E1) agrees with the full one.

**Propagation (comoving, ×1.51):**

| | ρ_c | 95 % upper (ρ_c ≥ 0) |
|---|---|---|
| measured | −0.19 ± 0.40 | 0.67 |
| photometric systematic (standards' null, −0.04 ± 0.15) subtracted, error in quadrature | **−0.14 ± 0.43** | **0.75** |
| and the leak (+0.34 ± 0.10) also subtracted | −0.49 ± 0.44 | 0.60 |
| R-150's QSO-target-star null taken literally (−5.7) | +5.6 | excluded by null 1 and the WISE-free estimator |

Neither systematic moves the limit by more than its error (0.05 and 0.34 against 0.40). **The stated limit includes the
photometric systematic and leaves the leak unsubtracted: ρ_c = −0.14 ± 0.43, < 0.75 at 95 %.**

## E4. Comparisons at matched comoving scales (post-M1)

**Ménard's profile.** ρ_c < 0.75 (95 %) is the ratio itself.

**What an unchanged galaxy-correlated dust density would give.** In the two-halo regime, with dust traced like matter, linear
bias and the same dust bias at both epochs, the rest-frame A_V at fixed comoving R scales as Ω_d b_d b_g D(z)² (1 + z)² (the
proper dust column).
- Inputs: b_ELG D = 0.84 ± 15 % (DESI's fiducial, arXiv:1611.00036), b = 1.1 ± 0.2 for Ménard's i < 21 sample (assumed),
  D(1.22)/D(0.36) = 0.675.
- That gives S = **1.65 [1.30, 2.10]**: ρ_c ≈ 1.65 if Ω_d b_d were the same at z = 1.22 as at 0.36.
- Measured over S: Ω_d b_d (z ≈ 1.2) / Ω_d b_d (z ≈ 0.36) = −0.09 [−0.35, 0.17], **< 0.48 (95 %, positive prior)**.
- An unchanged Ω_d b_d is 2.9σ away.

**The Mg II dust budget (Ménard & Fukugita 2012, arXiv:1204.1978).**
- Ω_dust(Mg II) = 2.3 ± 0.2 × 10⁻⁶ at z ≈ 0.5, scaling as (1 + z)^−1 (their fit): 1.56 × 10⁻⁶ at z = 1.22. They call it a lower
  limit, low by ~20 % from weak absorbers and ~20 % from obscured quasars.
- Ménard et al. 2010's own budget for the dust behind their profile is Ω_dust ≈ 5 × 10⁻⁶ (halos + IGM, z ~ 0.3; 2.1 × 10⁻⁶ inside
  R_vir).
- If the galaxy-correlated dust at z = 1.22 were the Mg II budget alone, ρ_c = 0.31 S = **0.51 [0.41, 0.65]**.
- Measured over that prediction: −0.27 [−1.12, 0.54], **< 1.56 (95 %)**. The limit is consistent with the absorber-selected
  budget and does not reach it; it is 1.5 times above it.

## E5. The LRG rows: not separated from systematics (moved; no ELG-LRG contrast claimed)

| foregrounds | ρ (proper, D1) | stellar null beside it (R-150, QSO-target stars) | WISE-free and spectroscopic estimators |
|---|---|---|---|
| LRG 0.4-0.6 / 0.6-0.8 / 0.8-1.1 | 2.33 ± 0.75 / 2.99 ± 0.71 / 2.66 ± 0.68 | | |
| LRG combined | 2.65 ± 0.48 | **+1.80 ± 0.74** (rotated foregrounds +0.49 ± 0.36) | g−r 2.51 ± 0.65, r−z 2.13 ± 0.58, fibre-spectrum slope 1.71 ± 0.58 |

Not separated from an LRG-tied systematic. E3 shows that the QSO-target-star null is a target-selection null, not a
photometric one, so the +1.8 does not by itself measure a photometric systematic around LRGs. The standards' null, the
LRG behind-the-quasar null and an LRG magnification-leak propagation were not run (out of this revision's scope). Until
they are, the LRG amplitude is not a result. The ELG-LRG contrast of R-150 (3.55 ± 0.93) is withdrawn.

## E6. Finding check (2026-11-02, fresh context-free Opus, after-writing prompt verbatim, the revised section alone; hints as D4)

VERDICT: novel but not sound: the limit propagates the systematic from a stellar null chosen after the fact that, by design, avoids the W1-dependent QSO targeting effect (−3.8 in ρ around ELGs) to which the background quasars are themselves exposed, so the 95 % limit and the 2.9σ no-evolution exclusion are overstated.

The referee objects that the propagated systematic should be the quasar-side control, not the standards. Propagating null 1
(−0.20 ± 0.31 proper, comoving −0.30 ± 0.47) instead gives ρ_c = 0.12 ± 0.62, < 1.29 at 95 % (ρ_c ≥ 0); against S = 1.65 that is
2.1σ. Novelty held (no galaxy-selected stack with foregrounds at z > 0.5; 2508.15490 and 2601.15961 added to the list). **Parked**:
the brief's last pass, and the objection is physics (which systematic the quasars carry), not presentation.
