# Results (R-138, cosmo-deep-14, 2026-11-01)

The bias-weighted near-infrared background intensity per unit redshift, b dI/dz, per SPHEREx spectral channel, from
cross-correlating masked, zodi-subtracted SPHEREx QR2 spectral images with DESI DR1 spectroscopic galaxies in redshift
slices (the clustering-redshift intensity estimator of arXiv:1810.00885; notation of arXiv:2109.10914, b_I d(νI_ν)/dz).
This file is the worker's own. Scripts: `spherex_desi_prep.py` (DESI DR1 compaction and footprint), `spherex_reduce.py`
(one L2 image to HEALPix-binned residual surface brightness), `spherex_xcorr.py` (maps, estimator, nulls, jackknife).
Data under `data/desi_dr1/` and `data/spherex/` (gitignored).

## Trajectory check (run in the dispatching session, 2026-11-01, context-free Opus; not rerun here)

"VERDICT: partly novel: the SPHEREx×galaxy-slice b dI/dz tomography observable and method are forecast in arXiv:2109.10914
(method demonstrated on GALEX×SDSS in arXiv:1810.00885), and a measured broadband NIR b dI/dz at 1.1/1.8 µm for z<1 exists
from CIBER×DESI Legacy/HSC (arXiv:2608.12116); what remains open is the first measurement with real SPHEREx data, per
spectral channel against DESI spectroscopic slices, separating H-alpha/[OIII]/Pa-alpha line emission from the continuum"

## Step 0: feasibility (2026-11-01, 21:20-21:55Z)

- **SPHEREx QR2 reachable.** IRSA TAP (`spherex.obscore`, `spherex.plane`, `spherex.artifact`) answers; the L2 images
  are public in the S3 bucket `nasa-irsa-spherex` (`qr2/level2/...`), 71.6 MB each (IMAGE, FLAGS, VARIANCE, ZODI, PSF,
  WCS-WAVE), ~24 MB/s per stream. QR2 holds ~182k wide-survey images per detector (D1-D6) and ~45k per detector in the
  two deep fields; QR3 started September 2026 (~32k per detector so far). Each image carries the pipeline zodiacal model
  (ZODI extension, "ModifiedKelsallModelWithHPFT") and a per-pixel wavelength (tabular 'W' WCS).
- **DESI DR1 reachable.** `data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5/`: NGC clustering
  catalogues (BGS_ANY 4.08 M, LRG 1.48 M, ELG_LOPnotqso 1.82 M, QSO 0.79 M) and random file 0 of each, ~20 MB/s.
- **Overlap.** The deep fields do not help: the south one is outside DESI and the north one is covered by a handful of
  DESI DR1 tiles; and 2109.10914 (Fig. "clus_spec", dashed vs solid) forecasts lower S/N in the deep fields than wide,
  sample variance winning. The wide survey covers all of DESI DR1. 1,252 deg² of DESI DR1 NGC has BGS, LRG and ELG all at
  >= 80 % of full random density and |b| > 40 deg; over it (padded by one NSIDE-32 pixel) QR2 has 10.5-10.7k images per
  detector, ~3.3 visits per channel per point.
- **Throughput.** Reduction (download + masking + binning) runs at ~1.8 images/s on this 4-core container with 12
  workers (network-bound): ~6,500 images/h. All of D1-D3 over the field is 32,005 images (~5 h); the budget allows ~2.3 h.
- **Signal-to-noise estimate.** Signal: b dI/dz ~ 3e-3 MJy/sr per unit z at 1 µm (νI_ν ~ 10 nW m^-2 sr^-1 of IGL,
  b_I ~ 1.5); with b_g ~ 1.5 and the annulus-weighted w̄_m ~ 0.01 per unit z (w_p ~ 25 Mpc/h at 2 Mpc/h over c/H ~
  2,600 Mpc/h), w̄_Ir ~ 4e-5 MJy/sr per slice per channel. Noise: the 49" block residual rms is ~0.04 MJy/sr (sky
  structure, not white noise); an annulus of 1-10 comoving Mpc averages ~10^3 blocks and ~3x10^3 independent annuli fit in
  1,250 deg², so σ(w̄) ~ 1-2e-5 MJy/sr per channel per slice from the map side, below the galaxy count side for slices
  with >10^4 galaxies. Expected S/N ~ 2-4 per channel per slice, ~5-10 per group of 3 channels. Feasible; the gate passes
  on the condition that the reduced half of the images covers the field (it does: the image order is random).

## Fixed before any cross-correlation (2026-11-01, 21:58Z)

- **Field:** the 1,252 deg² above (`data/spherex/region32.npy`, NSIDE-32 RING pixels: DESI DR1 BGS_ANY, LRG and
  ELG_LOPnotqso random density each >= 0.8 of its 90th percentile, |b| > 40 deg; RA 148-250, Dec -8 to +39). Galaxies
  and randoms are those of the DR1 NGC clustering catalogues inside it, weighted by WEIGHT.
- **Images:** QR2 wide survey (`spherex_qr2`), detectors D1-D3 (0.75-2.42 µm), every image whose centre falls in the
  field or one NSIDE-32 pixel around it, processed in a fixed random order (seed 138, `data/spherex/list_d123.csv`).
  Reduction stops at 00:15Z; the analysis uses every image reduced by then (a random subset by construction). D4-D6 are
  not processed (budget); Pa-alpha therefore enters only through D3 (z < 0.29).
- **Channels:** 17 log-uniform bins per detector between the QR2 band edges (D1 0.7455-1.1167, D2 1.1031-1.6476, D3
  1.6317-2.4237 µm), channel = the block centre's wavelength; 51 channels, R ~ 43. Channel groups if per-channel S/N is
  low: 6 per detector, channels {0-2, 3-5, 6-8, 9-11, 12-14, 15-16}.
- **Redshift slices:** dz = 0.1. BGS_ANY 0.05-0.45 (4), LRG 0.4-1.1 (7), ELG_LOPnotqso 0.8-1.6 (8), QSO 0.8-2.1 (13).
- **Masks:** SPHEREx FLAGS bad bits {0,1,2,4,5,6,7,9,10,11,15,17,19} and the SOURCE bit 21 (the pipeline's mask of known
  catalogue sources; ~35 % of pixels), then pixels > 5 robust σ above their channel strip's median, dilated by 2 pixels
  (12"); at map level, every NSIDE-2048 pixel within r(G) = 3' x 10^(0.2(12-G)) (min 1.7') of a Gaia DR3 star with G < 12
  (if the Gaia fetch fails, this mask is dropped and the row says so). Deviation from the brief: the Legacy Survey source
  mask is replaced by the SPHEREx SOURCE flag (a fixed catalogue mask applied by the pipeline); its flux limit is the
  pipeline's, not chosen here. b dI/dz is therefore that of the background left after this mask.
- **Zodi:** the ZODI extension is subtracted per pixel; then each channel strip's median of good pixels is subtracted per
  exposure (zodi-model error, offsets). That filter suppresses correlations on scales comparable to the strip (0.2 x 3.5
  deg); its transfer function is measured, not modelled: each slice's own galaxy overdensity map is passed through the
  identical per-(image, channel) footprint mean-subtraction and the estimator is run on it filtered and unfiltered; the
  ratio corrects w̄ per slice (and is reported).
- **Angular scales:** an annulus fixed in comoving Mpc per slice: θ_min = max(θ(1 cMpc), 3.5'), θ_max = θ(10 cMpc) at the
  slice centre (Planck 2018 flat ΛCDM, h = 0.674), weight W(θ) ∝ θ^-0.8 normalised over the annulus (1810.00885).
- **Estimator (1810.00885 eq. w_theta, wbar, wbar_to_dJdz):** w_Ir(θ, z) = <ΔI>_galaxies(θ) - <ΔI>_randoms(θ), mean
  residual intensity at separation θ around weighted galaxies minus around randoms, from pixel pair sums on gnomonic
  1.7'-pixel patches (one per NSIDE-8 pixel of the field, with a θ_max margin of intensity around it); w̄_Ir = ∫W w_Ir dθ;
  b_I dI_ν/dz = w̄_Ir / (b_r w̄_m) with w̄_m from eq. w_m (CAMB, halofit nonlinear P(k)), the same W and annulus.
  b_r from the slice's own clustering: w̄_gg (Landy-Szalay from the same pixel pairs) = b_r² w̄_m^gg, w̄_m^gg Limber for
  the slice's own n(z). Reported in I_ν (MJy/sr per unit z) and in ν I_ν (nW m^-2 sr^-1 per unit z).
- **Errors:** delete-one jackknife over the NSIDE-8 patches; checked against the scatter of the half-survey difference.
- **Nulls (each must be consistent with zero, χ² over all slices of a channel group):** (a) redshift-shuffled galaxies
  (galaxy redshifts permuted within each tracer); (b) Gaia DR3 stars (12 < G < 16) in place of galaxies; (c) the ZODI model
  alone run through the same pipeline; (d) half-survey difference: images split at the median exposure time, (A - B)/2.
- **Line model:** the continuum is smooth in z per channel (a quadratic in z per channel group across all tracers'
  slices); a line is the excess in the slices where Hα 0.6563, [OIII] 0.5007 (+Hβ 0.4861 blended at R~43), Pa-α 1.8751
  µm falls in the channel, λ_obs = λ_rest (1+z), for channel groups. Line significance = the excess over the quadratic
  in the crossing slices, with jackknife covariance.
- **CIBER comparison:** the channels containing 1.1 and 1.8 µm (D1 ch 16 / D2 ch 0; D3 ch 4) against
  arXiv:2608.12116's b_I dI/dz (its Fig. A2h_vs_z, Δz = 0.2 bins, z < 1), as ν I_ν.
- **Claim:** the table per channel group and slice with intervals; line detections with significance; the CIBER
  comparison. A null is reported as limits per channel group. No sentence about galaxy evolution beyond the table.

### Amendment (2026-11-01, 22:20Z, before any cross-correlation was computed or looked at)

The Gaia map-level mask as fixed above (G < 12, r = 3' x 10^(0.2(12-G))) removes 2.89 M of the 3.29 M map pixels:
162k stars with G < 12 lie in the star query box, and every one masks at least the ~10 NSIDE-2048 pixels of a 3' disc.
That was a counting error in the pre-registration, found from the mask fraction alone. Replaced by: Gaia DR3 G < 10,
r = 2' x 10^(0.2(10-G)) (min one pixel); stars 10 < G < 12 are left to the pipeline SOURCE flag, which already masks
catalogue sources at the native 6" pixel. Nothing else changes.

### Interim pipeline checks and second amendment (2026-11-01, 22:50-23:40Z; recorded before the final run)

The pipeline was run end to end on the 9,743 images reduced by 22:50Z, to catch bugs before the stop. It caught one:
the pair stage filled the partner galaxy/random count maps interleaved (G0, R0, G1, ...) while the field list reads
G0..G31, R0..R31, so every slice but the first measured its clustering against the wrong map (b_r ~ 0-0.2). Fixed
(60c078e2); the biases then come out 1.0-1.4 (BGS), 1.9-2.3 (LRG), 1.1-1.8 (ELG), 1.8-3.0 (QSO). The container also
restarted at ~23:25Z; the reduction resumed from the images already done (each image is its own file). The interim
results were seen (signal and nulls). Changes made after seeing them, each for the reason given, none to a choice that
sets the signal's amplitude:
- **QSO transfer:** the strip-filter transfer from the QSO slices' own template is noise (7-14k QSOs per slice; values
  -0.5 to 0.5 ± 0.1-0.4). The transfer depends on the annulus and the filter, not the tracer, so each QSO slice takes
  the transfer of the ELG slice with the same z range (QSO z > 1.6: ELG 1.5-1.6).
- **Shuffled null, expectation corrected:** shuffled galaxies keep their sky positions, so each shuffled slice is a
  random subset of its tracer and carries that tracer's projected (n(z)-averaged) signal; its expectation is not zero.
  The null is now: per tracer, the shuffled b dI/dz is flat in z (χ² about the tracer's own mean), and the line model
  fitted to it (with one offset per tracer) returns no lines. The interim "χ² against zero" column is kept and reads as
  the projected detection.
- **Added robustness check (not pre-registered):** the star null (12 < G < 16) is not zero: around stars the residual
  intensity at the annulus scales is comparable, in the same units, to the galaxy signal (star-associated light: PSF
  wings, ghosts, the clustering of stars). It biases the galaxy estimate only through a galaxy-star correlation, so the
  patches are split at the median star count and the estimate is repeated on each half.
- **Line significance calibration:** the joint Hα χ² of the line model is compared with the same statistic on the
  shuffled galaxies, which have no line by construction; a line counts as detected only if it beats both zero and the
  shuffled distribution.

## Result (final run, 22,322 images; R-138, 2026-11-02)

(The lines above this heading were lost from this file by the R-138 commit c520220, which wrote the report over it; restored verbatim from 586a3ad by R-143.)

Patches (jackknife regions): 51.

Per slice: tracer, z range, N galaxies in the field, clustering bias b_r (jackknife), strip-filter transfer.

| slice | tracer | z | N_gal | b_r | transfer |
|---|---|---|---|---|---|
| 0 | BGS_ANY | 0.05-0.15 | 251375 | 1.05 ± 0.04 | 0.43 ± 0.02 |
| 1 | BGS_ANY | 0.15-0.25 | 346031 | 1.18 ± 0.03 | 0.62 ± 0.02 |
| 2 | BGS_ANY | 0.25-0.35 | 284884 | 1.31 ± 0.03 | 0.63 ± 0.02 |
| 3 | BGS_ANY | 0.35-0.45 | 149285 | 1.42 ± 0.03 | 0.63 ± 0.02 |
| 4 | LRG | 0.40-0.50 | 51108 | 1.86 ± 0.04 | 0.63 ± 0.02 |
| 5 | LRG | 0.50-0.60 | 70438 | 2.01 ± 0.03 | 0.66 ± 0.02 |
| 6 | LRG | 0.60-0.70 | 86031 | 2.13 ± 0.04 | 0.71 ± 0.01 |
| 7 | LRG | 0.70-0.80 | 104646 | 2.19 ± 0.04 | 0.73 ± 0.01 |
| 8 | LRG | 0.80-0.90 | 109927 | 2.31 ± 0.04 | 0.76 ± 0.01 |
| 9 | LRG | 0.90-1.00 | 71006 | 2.33 ± 0.04 | 0.74 ± 0.01 |
| 10 | LRG | 1.00-1.10 | 31533 | 2.34 ± 0.06 | 0.63 ± 0.02 |
| 11 | ELG_LOPnotqso | 0.80-0.90 | 105596 | 1.08 ± 0.04 | 0.43 ± 0.05 |
| 12 | ELG_LOPnotqso | 0.90-1.00 | 123233 | 1.12 ± 0.04 | 0.51 ± 0.02 |
| 13 | ELG_LOPnotqso | 1.00-1.10 | 122021 | 1.28 ± 0.03 | 0.57 ± 0.02 |
| 14 | ELG_LOPnotqso | 1.10-1.20 | 113896 | 1.35 ± 0.03 | 0.60 ± 0.02 |
| 15 | ELG_LOPnotqso | 1.20-1.30 | 109486 | 1.47 ± 0.05 | 0.61 ± 0.02 |
| 16 | ELG_LOPnotqso | 1.30-1.40 | 103296 | 1.48 ± 0.05 | 0.57 ± 0.04 |
| 17 | ELG_LOPnotqso | 1.40-1.50 | 92358 | 1.63 ± 0.04 | 0.60 ± 0.02 |
| 18 | ELG_LOPnotqso | 1.50-1.60 | 62443 | 1.80 ± 0.08 | 0.55 ± 0.03 |
| 19 | QSO | 0.80-0.90 | 7223 | 1.80 ± 0.20 | 0.43 ± 0.05 |
| 20 | QSO | 0.90-1.00 | 8560 | 1.54 ± 0.21 | 0.51 ± 0.02 |
| 21 | QSO | 1.00-1.10 | 10137 | 1.77 ± 0.17 | 0.57 ± 0.02 |
| 22 | QSO | 1.10-1.20 | 11166 | 1.98 ± 0.16 | 0.60 ± 0.02 |
| 23 | QSO | 1.20-1.30 | 12452 | 2.14 ± 0.12 | 0.61 ± 0.02 |
| 24 | QSO | 1.30-1.40 | 13000 | 2.21 ± 0.17 | 0.57 ± 0.04 |
| 25 | QSO | 1.40-1.50 | 13519 | 2.32 ± 0.13 | 0.60 ± 0.02 |
| 26 | QSO | 1.50-1.60 | 14164 | 2.47 ± 0.15 | 0.55 ± 0.03 |
| 27 | QSO | 1.60-1.70 | 13859 | 2.51 ± 0.17 | 0.55 ± 0.03 |
| 28 | QSO | 1.70-1.80 | 13742 | 2.51 ± 0.21 | 0.55 ± 0.03 |
| 29 | QSO | 1.80-1.90 | 13191 | 2.71 ± 0.21 | 0.55 ± 0.03 |
| 30 | QSO | 1.90-2.00 | 12064 | 2.91 ± 0.21 | 0.55 ± 0.03 |
| 31 | QSO | 2.00-2.10 | 10644 | 3.04 ± 0.19 | 0.55 ± 0.03 |

b_I d(νI_ν)/dz per channel group (columns, central µm) and slice (rows), nW m^-2 sr^-1 per unit z, value ± jackknife σ:

| slice | 0.773 | 0.830 | 0.891 | 0.957 | 1.028 | 1.090 | 1.143 | 1.227 | 1.317 | 1.413 | 1.517 | 1.609 | 1.690 | 1.812 | 1.943 | 2.083 | 2.234 | 2.368 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BGS 0.05-0.15 | -8±16 | -10±14 | -10±12 | -5±10 | -3±8 | 19±100 | 3±5 | 5±7 | 3±6 | -1±3 | 0±4 | -1±4 | 5±3 | 3±3 | 1±1 | -0±2 | -1±2 | 1±1 |
| BGS 0.15-0.25 | 12±7 | 12±6 | 15±6 | 16±6 | 12±6 | 101±66 | 14±5 | 16±6 | 16±5 | 13±4 | 12±3 | 12±3 | 8±3 | 6±2 | 5±1 | 6±2 | 5±2 | 4±1 |
| BGS 0.25-0.35 | 6±9 | -0±7 | 4±5 | 8±5 | 4±5 | -37±63 | 9±4 | 10±4 | 11±4 | 10±4 | 4±3 | 4±3 | 6±3 | 5±2 | 4±1 | 3±1 | 2±2 | 3±1 |
| BGS 0.35-0.45 | 11±11 | 9±8 | 8±7 | 4±6 | 4±7 | 10±80 | 8±5 | 7±5 | 7±4 | 3±5 | 4±4 | 3±5 | 4±3 | 4±2 | 4±1 | 5±2 | 6±2 | 5±1 |
| LRG 0.40-0.50 | 6±7 | 7±7 | 11±7 | 11±6 | 7±5 | -67±68 | 13±5 | 12±4 | 11±4 | 11±3 | 7±4 | 7±4 | 5±3 | 5±2 | 7±1 | 7±2 | 7±2 | 6±2 |
| LRG 0.50-0.60 | 11±9 | 17±8 | 22±6 | 20±5 | 17±5 | 87±53 | 8±7 | 11±5 | 16±4 | 14±3 | 12±3 | 12±2 | 8±2 | 10±2 | 9±1 | 7±2 | 6±2 | 7±2 |
| LRG 0.60-0.70 | 2±7 | 1±5 | 9±5 | 13±5 | 10±4 | 75±53 | 13±3 | 8±4 | 12±3 | 15±4 | 9±3 | 9±2 | 8±2 | 8±2 | 8±1 | 8±1 | 9±1 | 8±1 |
| LRG 0.70-0.80 | 11±7 | 11±6 | 11±5 | 11±4 | 8±5 | 37±58 | 9±4 | 8±3 | 9±3 | 13±3 | 13±3 | 10±3 | 6±3 | 7±2 | 6±1 | 6±2 | 6±2 | 8±1 |
| LRG 0.80-0.90 | -5±8 | -11±9 | 0±5 | 1±5 | -4±5 | 0±46 | 7±4 | 9±3 | 10±3 | 5±3 | 3±3 | 2±3 | 4±3 | 5±2 | 8±1 | 10±2 | 9±2 | 8±2 |
| LRG 0.90-1.00 | 5±8 | 2±7 | 2±5 | 5±4 | -2±4 | 151±57 | 6±3 | 7±4 | 10±3 | 6±3 | 5±4 | 8±4 | 5±2 | 6±2 | 8±1 | 9±2 | 9±2 | 8±1 |
| LRG 1.00-1.10 | 4±8 | 9±10 | 10±8 | 11±8 | -1±8 | 117±87 | 8±6 | 13±6 | 9±7 | 10±4 | 8±4 | 7±4 | 11±4 | 11±2 | 11±1 | 12±2 | 12±3 | 8±2 |
| ELG 0.80-0.90 | -18±14 | -24±14 | -11±12 | -13±12 | -20±13 | 18±141 | 11±13 | 24±14 | 13±9 | 7±6 | 4±8 | -1±7 | 4±7 | 7±4 | 9±4 | 6±4 | 6±4 | 6±4 |
| ELG 0.90-1.00 | 6±14 | 0±13 | -11±13 | -3±11 | -9±12 | 29±126 | -3±8 | 17±8 | 15±8 | 3±5 | 3±5 | 1±6 | 14±6 | 12±4 | 5±2 | 6±3 | 4±4 | 5±2 |
| ELG 1.00-1.10 | 2±13 | 1±14 | -4±10 | -1±10 | -14±10 | 77±75 | 3±8 | 3±8 | -5±10 | -0±7 | 4±5 | 7±7 | 2±5 | 5±4 | 6±2 | 10±3 | 11±3 | 7±2 |
| ELG 1.10-1.20 | 15±12 | 9±11 | 9±10 | 9±9 | 21±10 | -78±102 | 15±6 | 5±10 | 11±6 | 11±6 | 4±5 | 6±5 | 7±4 | 7±3 | 6±2 | 9±3 | 8±4 | 6±2 |
| ELG 1.20-1.30 | 12±8 | 7±11 | 4±10 | 7±10 | 5±10 | -138±78 | 8±6 | 3±10 | 8±6 | 10±6 | 7±5 | 4±5 | 6±4 | 3±3 | 7±2 | 8±3 | 10±3 | 8±2 |
| ELG 1.30-1.40 | 4±8 | -27±15 | -13±11 | -7±11 | -22±12 | 39±115 | -4±8 | 0±6 | -4±7 | -3±9 | -3±5 | -2±6 | -0±6 | -1±5 | 5±2 | 11±3 | 10±4 | 6±3 |
| ELG 1.40-1.50 | -18±11 | 1±14 | -8±8 | -9±10 | -12±10 | 21±115 | -7±9 | 2±8 | 6±7 | 3±5 | 1±5 | -7±6 | 3±5 | 3±3 | 3±2 | -0±3 | 4±4 | 2±2 |
| ELG 1.50-1.60 | 2±13 | 22±17 | 7±13 | 1±13 | 11±11 | 15±136 | 2±8 | 6±9 | 11±8 | 10±9 | 5±7 | -1±7 | 9±6 | 8±5 | 6±3 | 2±4 | 5±5 | 6±3 |
| QSO 0.80-0.90 | -3±35 | 2±20 | 18±16 | 7±16 | 1±18 | 299±226 | 6±14 | 41±23 | 31±14 | 12±17 | 19±13 | 21±10 | 0±7 | 4±5 | 15±4 | 11±5 | 14±5 | 13±4 |
| QSO 0.90-1.00 | 38±25 | 14±21 | -21±24 | 4±18 | -3±20 | 403±206 | 28±12 | 29±12 | 35±22 | 2±16 | -0±15 | 6±11 | 4±7 | 6±6 | 11±4 | 17±7 | 9±6 | 11±5 |
| QSO 1.00-1.10 | 28±23 | 10±23 | 2±20 | 20±21 | 7±20 | -0±149 | 20±12 | 16±14 | 6±13 | 21±12 | 11±11 | 15±12 | 4±7 | 8±6 | 12±3 | 16±4 | 13±4 | 11±3 |
| QSO 1.10-1.20 | 5±15 | 18±18 | 13±15 | 11±15 | -3±14 | 36±154 | 10±11 | 15±14 | 18±9 | 28±18 | 3±8 | 6±6 | 7±7 | 7±5 | 5±2 | 5±3 | 5±4 | 3±3 |
| QSO 1.20-1.30 | -6±12 | -8±12 | -13±13 | -5±11 | -5±14 | -13±136 | -2±9 | -9±12 | -3±6 | -4±7 | 3±7 | 3±7 | 5±6 | 2±4 | 3±3 | 5±4 | 3±4 | 2±3 |
| QSO 1.30-1.40 | 37±26 | -15±16 | 4±12 | 3±13 | 0±17 | 128±160 | 1±9 | 0±14 | -13±18 | 3±8 | 4±8 | 8±8 | -11±6 | -3±4 | 4±3 | 10±4 | 13±5 | 7±3 |
| QSO 1.40-1.50 | 4±15 | 14±17 | -11±11 | -8±10 | -15±12 | 76±133 | -24±9 | -7±9 | -13±7 | 5±14 | -5±6 | 1±6 | 4±4 | 3±3 | -0±3 | -5±4 | -2±4 | 0±3 |
| QSO 1.50-1.60 | 10±15 | 23±19 | 10±14 | 10±14 | 8±13 | 11±180 | 8±9 | 14±12 | -6±9 | 8±9 | 5±8 | 2±8 | 10±7 | 7±5 | 5±3 | -0±4 | 3±4 | 6±3 |
| QSO 1.60-1.70 | 14±14 | 18±18 | 9±17 | 17±16 | 14±14 | 156±146 | 4±10 | -4±12 | 14±9 | 5±9 | 17±10 | 10±8 | 5±6 | 3±4 | 3±3 | 6±4 | 3±4 | 6±3 |
| QSO 1.70-1.80 | 16±12 | 10±19 | 12±13 | 3±13 | -1±16 | 52±139 | 4±8 | -30±27 | -6±11 | -2±8 | 2±6 | 4±6 | 0±6 | 6±3 | 3±3 | 1±4 | 2±4 | 4±3 |
| QSO 1.80-1.90 | -5±17 | -8±17 | -3±12 | -6±10 | 5±12 | 6±152 | -9±8 | 2±7 | 0±8 | 3±5 | -1±9 | 4±7 | 6±6 | 8±5 | -0±3 | -3±4 | -4±4 | 1±4 |
| QSO 1.90-2.00 | 22±10 | 20±13 | 22±11 | 18±11 | 22±13 | 139±139 | 8±6 | -1±11 | 7±6 | 10±5 | 3±7 | -1±6 | 5±5 | 6±5 | 7±2 | 5±3 | 5±3 | 0±3 |
| QSO 2.00-2.10 | -4±24 | 5±19 | 8±18 | -6±17 | 1±18 | 219±158 | 8±12 | -8±15 | -8±14 | 3±10 | 8±10 | 14±12 | 1±8 | 0±6 | 2±4 | -2±4 | -3±4 | 3±4 |

Nulls: χ² of each null over the 32 slices per group against zero with jackknife σ (dof 32; p from χ²):

| group µm | signal χ² | shuffled χ² | stars χ² | zodi χ² | half-diff χ² | half-diff σ / jackknife σ |
|---|---|---|---|---|---|---|
| 0.773 | 33/32 (p=4.3e-01) | 45/32 (p=5.8e-02) | 12/32 (p=1.0e+00) | 35/32 (p=3.1e-01) | 21/32 (p=9.4e-01) | 0.70 |
| 0.830 | 36/32 (p=2.8e-01) | 40/32 (p=1.4e-01) | 21/32 (p=9.3e-01) | 34/32 (p=3.6e-01) | 30/32 (p=5.6e-01) | 0.76 |
| 0.891 | 52/32 (p=1.5e-02) | 78/32 (p=8.9e-06) | 21/32 (p=9.2e-01) | 24/32 (p=8.5e-01) | 21/32 (p=9.3e-01) | 0.54 |
| 0.957 | 55/32 (p=6.7e-03) | 119/32 (p=6.9e-12) | 53/32 (p=1.1e-02) | 29/32 (p=6.3e-01) | 33/32 (p=4.3e-01) | 0.59 |
| 1.028 | 48/32 (p=3.6e-02) | 72/32 (p=7.3e-05) | 71/32 (p=8.1e-05) | 32/32 (p=4.5e-01) | 32/32 (p=4.9e-01) | 0.75 |
| 1.090 | 33/32 (p=4.0e-01) | 48/32 (p=3.6e-02) | 269/32 (p=3.0e-39) | 36/32 (p=2.9e-01) | 37/32 (p=2.5e-01) | 1.20 |
| 1.143 | 94/32 (p=5.7e-08) | 165/32 (p=7.7e-20) | 158/32 (p=1.2e-18) | 36/32 (p=3.0e-01) | 27/32 (p=7.0e-01) | 0.70 |
| 1.227 | 85/32 (p=1.1e-06) | 161/32 (p=3.6e-19) | 117/32 (p=1.1e-11) | 26/32 (p=7.5e-01) | 46/32 (p=5.3e-02) | 0.88 |
| 1.317 | 134/32 (p=1.9e-14) | 238/32 (p=2.7e-33) | 51/32 (p=1.7e-02) | 30/32 (p=5.6e-01) | 28/32 (p=6.6e-01) | 0.78 |
| 1.413 | 116/32 (p=2.1e-11) | 282/32 (p=7.0e-42) | 69/32 (p=1.6e-04) | 36/32 (p=2.9e-01) | 26/32 (p=7.5e-01) | 0.69 |
| 1.517 | 93/32 (p=7.7e-08) | 234/32 (p=1.4e-32) | 147/32 (p=1.2e-16) | 25/32 (p=8.0e-01) | 37/32 (p=2.5e-01) | 0.81 |
| 1.609 | 97/32 (p=1.9e-08) | 208/32 (p=9.4e-28) | 160/32 (p=6.1e-19) | 35/32 (p=3.3e-01) | 32/32 (p=4.8e-01) | 0.77 |
| 1.690 | 97/32 (p=1.7e-08) | 196/32 (p=2.1e-25) | 114/32 (p=3.7e-11) | 32/32 (p=4.6e-01) | 53/32 (p=1.2e-02) | 0.86 |
| 1.812 | 175/32 (p=1.0e-21) | 358/32 (p=1.0e-56) | 186/32 (p=1.5e-23) | 28/32 (p=6.6e-01) | 44/32 (p=8.4e-02) | 0.80 |
| 1.943 | 568/32 (p=2.7e-99) | 1135/32 (p=6.4e-218) | 206/32 (p=2.7e-27) | 26/32 (p=7.7e-01) | 42/32 (p=1.2e-01) | 0.74 |
| 2.083 | 294/32 (p=4.0e-44) | 558/32 (p=3.2e-97) | 131/32 (p=7.3e-14) | 33/32 (p=4.0e-01) | 33/32 (p=4.4e-01) | 0.70 |
| 2.234 | 254/32 (p=1.8e-36) | 430/32 (p=3.0e-71) | 106/32 (p=8.0e-10) | 31/32 (p=5.0e-01) | 34/32 (p=3.7e-01) | 0.70 |
| 2.368 | 348/32 (p=7.5e-55) | 646/32 (p=1.5e-115) | 40/32 (p=1.7e-01) | 40/32 (p=1.5e-01) | 28/32 (p=6.5e-01) | 0.66 |

Shuffled-redshift null (amended): χ² of the shuffled b dI/dz about each tracer's own mean (dof = slices - 4), per group:

| group µm | shuffled flatness χ²/dof (p) |
|---|---|
| 0.773 | 21/28 (p=8.1e-01) |
| 0.830 | 24/28 (p=6.8e-01) |
| 0.891 | 21/28 (p=8.2e-01) |
| 0.957 | 21/28 (p=8.2e-01) |
| 1.028 | 29/28 (p=4.0e-01) |
| 1.090 | 29/28 (p=3.9e-01) |
| 1.143 | 23/28 (p=7.5e-01) |
| 1.227 | 19/28 (p=9.0e-01) |
| 1.317 | 18/28 (p=9.2e-01) |
| 1.413 | 18/28 (p=9.3e-01) |
| 1.517 | 30/28 (p=3.8e-01) |
| 1.609 | 30/28 (p=3.6e-01) |
| 1.690 | 25/28 (p=6.5e-01) |
| 1.812 | 33/28 (p=2.3e-01) |
| 1.943 | 46/28 (p=1.7e-02) |
| 2.083 | 41/28 (p=5.1e-02) |
| 2.234 | 34/28 (p=2.1e-01) |
| 2.368 | 32/28 (p=2.8e-01) |

Line model per group: continuum quadratic in z + line amplitudes × the fraction of each slice whose line λ_obs falls in the group; amplitudes in nW m^-2 sr^-1 per unit z, σ from jackknife refits:

| group µm | Halpha (A ± σ, S/N) | OIII (A ± σ, S/N) | Paalpha (A ± σ, S/N) | continuum χ²/dof |
|---|---|---|---|---|
| 0.773 | 8 ± 14 (0.6) | 8 ± 10 (0.8) | — | 21/27 |
| 0.830 | -13 ± 13 (-1.1) | -4 ± 5 (-0.7) | — | 26/27 |
| 0.891 | -12 ± 12 (-0.9) | 4 ± 6 (0.6) | — | 31/27 |
| 0.957 | -1 ± 7 (-0.2) | -9 ± 7 (-1.3) | — | 24/27 |
| 1.028 | 18 ± 7 (2.4) | -4 ± 7 (-0.6) | — | 26/27 |
| 1.090 | 45 ± 64 (0.7) | -236 ± 121 (-1.9) | — | 20/27 |
| 1.143 | 2 ± 4 (0.5) | -3 ± 6 (-0.5) | — | 32/27 |
| 1.227 | 2 ± 4 (0.5) | -5 ± 7 (-0.6) | — | 17/27 |
| 1.317 | 1 ± 4 (0.3) | 7 ± 8 (0.9) | — | 32/27 |
| 1.413 | 3 ± 6 (0.5) | -1 ± 7 (-0.1) | — | 33/27 |
| 1.517 | -3 ± 5 (-0.7) | 11 ± 10 (1.2) | — | 26/27 |
| 1.609 | -9 ± 5 (-1.6) | — | — | 26/28 |
| 1.690 | 5 ± 5 (1.1) | — | — | 18/28 |
| 1.812 | 4 ± 4 (0.9) | — | — | 24/28 |
| 1.943 | 11 ± 3 (4.0) | — | -5 ± 7 (-0.8) | 36/27 |
| 2.083 | 4 ± 13 (0.3) | — | -3 ± 3 (-1.2) | 39/27 |
| 2.234 | — | — | 4 ± 3 (1.2) | 27/28 |
| 2.368 | — | — | -2 ± 3 (-0.7) | 25/28 |

Halpha: 16 groups, χ² against zero 31.2 (p=1.26e-02); inverse-variance mean amplitude 4 ± 1.

OIII: 11 groups, χ² against zero 10.1 (p=5.24e-01); inverse-variance mean amplitude -1 ± 2.

Paalpha: 4 groups, χ² against zero 4.0 (p=4.03e-01); inverse-variance mean amplitude -1 ± 2.

Halpha_shuf: 16 groups, χ² against zero 32.0 (p=1.00e-02); inverse-variance mean amplitude 1 ± 1.

OIII_shuf: 11 groups, χ² against zero 6.4 (p=8.46e-01); inverse-variance mean amplitude -0 ± 1.

Paalpha_shuf: 4 groups, χ² against zero 10.2 (p=3.72e-02); inverse-variance mean amplitude -1 ± 1.

Summary per group: inverse-variance mean of b_I d(νI_ν)/dz over slices (σ from the jackknife of the mean), all slices, z < 0.5 and z > 0.5; and the stellar-density split (patches below/above the median Gaia 12<G<16 count; added after the interim look):

| group µm | all z | z < 0.5 | z ≥ 0.5 | low-star half | high-star half |
|---|---|---|---|---|---|
| 0.773 | 5.9 ± 1.6 | 7.7 ± 3.0 | 5.4 ± 1.7 | 14.2 ± 3.7 | 6.4 ± 1.3 |
| 0.830 | 4.8 ± 1.7 | 6.3 ± 3.4 | 4.2 ± 2.0 | 9.7 ± 4.9 | 3.8 ± 1.7 |
| 0.891 | 5.7 ± 1.5 | 7.7 ± 3.0 | 5.0 ± 1.9 | 12.3 ± 3.6 | 4.1 ± 1.5 |
| 0.957 | 6.7 ± 1.3 | 8.3 ± 2.9 | 6.2 ± 1.4 | 13.6 ± 3.4 | 5.3 ± 1.2 |
| 1.028 | 3.7 ± 1.3 | 6.0 ± 2.7 | 2.9 ± 1.7 | 6.2 ± 2.5 | 2.8 ± 1.3 |
| 1.090 | 42.2 ± 19.1 | 2.6 ± 33.3 | 54.0 ± 21.0 | 31.4 ± 47.6 | 47.7 ± 19.9 |
| 1.143 | 7.5 ± 1.2 | 9.3 ± 2.6 | 6.8 ± 1.2 | 9.3 ± 3.6 | 6.5 ± 1.0 |
| 1.227 | 8.3 ± 1.2 | 10.5 ± 2.4 | 7.5 ± 1.4 | 10.8 ± 2.2 | 7.2 ± 1.2 |
| 1.317 | 9.0 ± 1.4 | 9.8 ± 2.6 | 8.7 ± 1.3 | 8.5 ± 1.8 | 8.9 ± 1.6 |
| 1.413 | 7.7 ± 0.9 | 7.1 ± 1.8 | 7.8 ± 0.8 | 9.4 ± 1.9 | 7.4 ± 1.0 |
| 1.517 | 6.4 ± 0.8 | 6.1 ± 1.5 | 6.5 ± 1.0 | 9.4 ± 1.5 | 5.7 ± 0.8 |
| 1.609 | 6.2 ± 0.9 | 5.9 ± 1.4 | 6.4 ± 1.1 | 6.9 ± 1.8 | 6.0 ± 1.0 |
| 1.690 | 5.9 ± 0.7 | 5.8 ± 1.4 | 5.9 ± 0.9 | 7.1 ± 2.6 | 5.2 ± 0.7 |
| 1.812 | 6.0 ± 0.5 | 4.8 ± 1.0 | 6.4 ± 0.7 | 6.8 ± 1.7 | 5.4 ± 0.5 |
| 1.943 | 6.1 ± 0.4 | 4.2 ± 0.5 | 7.1 ± 0.4 | 6.3 ± 0.7 | 6.1 ± 0.4 |
| 2.083 | 6.1 ± 0.5 | 3.8 ± 0.7 | 7.2 ± 0.7 | 5.2 ± 0.6 | 6.2 ± 0.6 |
| 2.234 | 6.2 ± 0.5 | 3.9 ± 0.7 | 7.3 ± 0.7 | 5.4 ± 0.9 | 6.5 ± 0.6 |
| 2.368 | 5.7 ± 0.3 | 3.5 ± 0.5 | 6.8 ± 0.5 | 6.4 ± 0.7 | 5.4 ± 0.3 |

CIBER comparison (arXiv:2608.12116, CIBER × DESI-LS, read off its Fig. a2h_vs_redshift); ours: inverse-variance mean of the slices whose centre lies in the CIBER bin, σ from jackknife of that mean:

| λ | z | CIBER | SPHEREx group (µm) | this work |
|---|---|---|---|---|
| 1.1 | 0.0-0.2 | 100 ± 30 | 1.090 | 19 ± 100 |
| 1.1 | 0.2-0.4 | 90 ± 18 | 1.090 | 28 ± 50 |
| 1.1 | 0.4-0.6 | 23 ± 9 | 1.090 | 24 ± 40 |
| 1.1 | 0.6-0.8 | < 5.5 (95%) | 1.090 | 58 ± 41 |
| 1.1 | 0.8-1.0 | < 4.0 (95%) | 1.090 | 69 ± 41 |
| 1.8 | 0.0-0.2 | 140 ± 40 | 1.812 | 3 ± 3 |
| 1.8 | 0.2-0.4 | 63 ± 15 | 1.812 | 5 ± 1 |
| 1.8 | 0.4-0.6 | 35 ± 8 | 1.812 | 7 ± 2 |
| 1.8 | 0.6-0.8 | < 14 (95%) | 1.812 | 7 ± 1 |
| 1.8 | 0.8-1.0 | < 7.5 (95%) | 1.812 | 6 ± 2 |
CIBER comparison at the adjacent well-covered 1.1 µm-side groups (the 1.090 µm group is D1's last two channels, at the
filter edge, with ~10x the noise): 1.143 µm (D2 channels 0-2): z 0-0.2 2.6 ± 4.6, 0.2-0.4 10.9 ± 3.3, 0.4-0.6 10.4 ±
3.4, 0.6-0.8 12.0 ± 2.3, 0.8-1.0 6.8 ± 2.8; 1.028 µm (D1 12-14): -2.9 ± 8.4, 7.7 ± 3.5, 10.1 ± 3.8, 9.1 ± 3.2, -4.0 ± 3.4.

Significance of the all-slice mean per group (value/σ): 0.77 µm 3.8, 0.83 2.8, 0.89 3.8, 0.96 5.2, 1.03 2.8, 1.09 2.2,
1.14 6.2, 1.23 6.6, 1.32 6.5, 1.41 8.4, 1.52 8.4, 1.61 7.1, 1.69 7.9, 1.81 11.8, 1.94 17.4, 2.08 11.6, 2.23 12.0,
2.37 16.5.

### What the table says

- **Detected:** b_I d(νI_ν)/dz in every channel group from 0.75 to 2.42 µm, 3-17σ per group over the 32 slices, at
  4-9 nW m^-2 sr^-1 per unit z, nearly flat in wavelength and in redshift over 0.1 < z < 2.1; at λ > 1.8 µm the z < 0.5
  slices sit lower (3.5-4.8) than z > 0.5 (6.4-7.3). The lowest slice (BGS 0.05-0.15) is consistent with zero in every
  group.
- **Lines: not detected.** The joint Hα χ² over 16 groups (31.2, p = 0.013) equals the same statistic on shuffled
  galaxies, which carry no line (32.0, p = 0.010), so the line test's jackknife errors are not calibrated at that level
  and no line counts; [OIII] (p = 0.52) and Pa-α (p = 0.40) are null. Per-group 1σ line limits: 3-14 nW m^-2 sr^-1 per
  unit z (table), i.e. a line no stronger than the continuum in the crossing slice.
- **Nulls:** zodi model alone, χ² 24-40 / 32 in every group (pass); half-survey difference 21-53 / 32 (pass; its scatter
  is 0.54-0.88 of the jackknife σ, so the jackknife includes sample variance); shuffled galaxies flat per tracer, 18-46 / 28
  (pass, 1.94 µm at p = 0.017 the worst). Stars: residual intensity around 12 < G < 16 stars at the same annuli is not
  zero (χ² 40-269 / 32): star-associated light. The patch split by star count does not move the galaxy signal up where
  stars are dense (low-star half 6-14, high-star half 3-9 nW m^-2 sr^-1 per unit z), so that light does not carry it.
- **CIBER:** at 1.8 µm our values (3.4 ± 2.7, 5.4 ± 1.2, 6.9 ± 1.5 at z 0-0.2, 0.2-0.4, 0.4-0.6) are 4-25x below
  CIBER x DESI-LS (140 ± 40, 63 ± 15, 35 ± 8), 3.4-3.9σ in each bin on CIBER's error alone; at 0.6-1.0 both are
  consistent (CIBER limits < 14, < 7.5; ours 7.2 ± 1.3, 6.4 ± 1.5). At 1.14 µm the same pattern (2.6-10.9 vs 100, 90,
  23). Our values sit inside the IGL-model band CIBER drew (5-20 for log M200 = 12-14). The two measurements differ in
  source mask (SPHEREx pipeline catalogue mask of every Reference Catalog point source vs CIBER's brighter cut), in scale
  (our 1-10 comoving Mpc annulus vs their two-halo amplitude at ℓ < 2000), in tracer (spectroscopic vs photometric
  slices) and in field (1,250 deg² vs 5 fields of 4 deg²); the table cannot say which of these carries the factor.

## Blind review (finding check, 2026-11-02, context-free Opus, gate prompt verbatim, the section alone; `blind-review.md`)

"VERDICT: novel but not sound: the amplitude and redshift shape are set by a catalogue source mask and clipping that remove galaxy-correlated light, which the transfer correction does not model, as the factor-of-two star-density split and the failed stellar null show, so the headline amplitudes and the CIBER ratio are not supported."

Strongest objection (the referee's): the SPHEREx SOURCE mask (every catalogue point source, ~31 % of pixels) and the
5σ clip remove light where the galaxies are, by an amount that depends on redshift and star density; the measured
transfer models only the per-strip mean subtraction, so the amplitude and its z shape are set by the mask. It also
reads the star-density split (low-star half ~2x the high-star half at 0.77-1.03 µm) as a sign the amplitude is not under
control, not as reassurance, and notes the jackknife χ² tests use a 51-patch covariance. Novelty held: no prior SPHEREx
x galaxy measurement found (closest: arXiv:2109.10914 forecast, 1810.00885 and 2504.05384 method, 2608.12116 CIBER).
**A miss on physics: parked (brief: one revision on presentation, park on physics).** What stands: a detected
SPHEREx x DESI correlation in all 18 groups with clean zodi, half-difference and shuffled nulls, and no line detected;
what does not: the amplitude as b_I dI/dz of the background, its z shape, the CIBER ratio.

## R-143: the source-mask transfer by injection (cosmo-deep-14b, 2026-11-02; fixed at 01:30Z, before any transfer was computed)

Why: the finding check above. The fix is a measurement: paint mock galaxy-correlated light into the real SPHEREx images,
run the identical chain (FLAGS bad bits + SOURCE bit, 5σ strip clip dilated 2 px, ZODI subtracted, strip median, 49"
blocks to NSIDE 2048), and measure what the estimator gets back. Scripts: `spherex_field.py` (the field and image list
rebuilt from the rule above), `spherex_gaia_fetch.sh` (Gaia DR3 G < 16, IRSA TAP), `spherex_inject.py` (reduction with
injection), `spherex_transfer.py` (maps, pair sums, transfer, corrected table).

- **Data rebuilt** (`data/` was empty): DESI DR1 NGC clustering catalogues and random file 0 (5 min); the field from its
  rule gives 373 NSIDE-32 pixels, 1,252.2 deg², RA 147.7-250.3, Dec -8.4 to 38.7 (R-138: 1,252 deg²); the D1-D3 image
  list from IRSA TAP (D1 10,562, D2 10,726 images, as R-138); Gaia DR3 G < 16 over the field box.
- **Mock** (`spherex_inject.py` docstring has the full definition). Every DESI galaxy inside its tracer's slice range
  carries flux F = WEIGHT × A(z) Δz / (K(z) n_ts) (n_ts the tracer-slice's weighted density in the field, K(z) the number
  of tracers covering z), so the mock's mean d(νI_ν)/dz is A(z), grey in νI_ν. Two components: compact (σ = 1 native
  pixel at the galaxy: the galaxy's own light, which the SOURCE mask can remove) and smooth (σ = 1 comoving Mpc around the
  galaxy: unresolved light correlated with it). Variants: A, A(z) = 6 nW m^-2 sr^-1 per unit z, half compact half smooth
  (fiducial); B, A × 5 (clip non-linearity); C, the smooth half alone; D, the compact half alone; E, A(z) = 6 ((1+z)/1.8)^3,
  half/half (rising shape). The real image is reduced as in R-138; real + A and real + B each with their own clip; C, D, E
  through A's mask with their own strip medians. Kept: the chain's output minus the real image's output on the same mask
  (the recovered mock), the mock on good pixels without the strip median (the mask alone), and the mock over every pixel
  of every block (the truth: no mask, clip or filter).
- **Transfer.** T_k(s, g) = w̄[recovered mock] / w̄[truth], both through R-138's estimator (same galaxy and random centres
  per slice, same annulus, bins and θ^-0.8 weights; w̄ = Σ_bins weight × (⟨I⟩_galaxies − ⟨I⟩_randoms)), per slice s and
  channel group g, jackknife over the 51 NSIDE-8 patches. b_r and w̄_m cancel in the ratio. Split into the mask alone
  (mock on good pixels, no strip median, over truth) and the strip filter (the rest).
- **Correction.** R-138's values were divided by the strip-filter transfer tr(s) alone. Corrected value = R-138 value ×
  tr(s) / T_A(s, g), with T_A's jackknife error and R-138's propagated in quadrature. The composition is not known, so
  the corrected value is bracketed by T_C (all the light smooth: the light outside catalogued sources) and T_D (all of it
  in the galaxies' own compact light). **Rule fixed now:** the amplitude is claimed as a measurement of the light outside
  catalogued sources (R-138 × tr / T_C) if T_A is stable (B/A and E/A within 15 % per group, averaged over slices); the
  total galaxy-correlated light is stated as a range from T_C to T_D, and as a lower bound only if T_D / T_C < 0.5 (then
  the compact light the mask removes is not constrained).
- **Images:** the rebuilt list in its own random order (seed 138), reduced from the top until 04:15Z; the transfer and the
  nulls below use every image reduced by then (the R-138 values are corrected with it; they are not recomputed).
- **Stellar null, star mask fixed:** a fixed-flux block mask of every Gaia DR3 star with G < 16 (block centre within
  1' × 10^(0.2 (12 − G)) + 12"), on top of R-138's masks; the 12 < G < 16 star null and the galaxy signal recomputed on
  this subset with and without it. Pass: χ² against zero with p > 0.01 in every group. Also the galaxy-star
  cross-correlation (DESI galaxies around the same stars, minus randoms, same annuli) and the contamination it implies,
  w̄_gs × N_* × S (N_* the stars expected inside θ_max of a point, S the star null value).
- **Star-density split, predicted:** T_A in the low-star and high-star halves of the patches (R-138's split at the median
  12 < G < 16 count); the predicted low/high ratio of the signal is T_A(low)/T_A(high), compared with R-138's observed one.
- **CIBER** at 1.1 µm (groups 1.090 and 1.143) and 1.8 µm (1.812), corrected as above, same z bins.
- **Added 02:00Z, after a code test on 300 images and 2 patches (numbers not usable, too few patches), before the final
  transfer:** dividing by a noisy per-cell T is unstable. If T_C's jackknife σ exceeds 0.2 in more than a quarter of the
  (slice, group) cells, the correction uses T per slice averaged (inverse variance) over the six groups of each detector:
  T is set by the image geometry, that detector's masks and the slice's annulus, not by the channel within a detector.

## R-143 result (2026-11-02, 04:15-04:45Z; the injection run, 14,508 images)

Reduction with injection 01:34-04:15Z (2 h 41 min, ~5,400 images/h on 4 cores with 10 workers: network-bound, the five extra chains cost ~2 CPU-s per image); maps, pair sums (51 patches, 11 s each once BLAS was single-threaded and catalogues preselected per patch) and the estimator 04:15-04:45Z. Output of `spherex_transfer.py report` and `spherex_r143_extra.py`, verbatim:

Images reduced with injection: 14508 (D1/D2/D3: [4836, 4812, 4860]); patches 51; star-mask block fraction 0.31.

Full-chain transfer T = w̄[recovered mock] / w̄[truth] per channel group, inverse-variance mean over the 32 slices (jackknife σ); the mask alone (mock on good pixels, no strip median) beside it; B is divided by 5:

| group µm | T_A (fiducial) | T_B (×5) | T_C (smooth) | T_D (compact) | T_E (rising z) | mask alone A | mask alone C | mask alone D | T_A, star mask |
|---|---|---|---|---|---|---|---|---|---|
| 0.773 | 0.509 ± 0.018 | 0.349 ± 0.012 | 0.686 ± 0.025 | 0.300 ± 0.013 | 0.360 ± 0.014 | 0.72 ± 0.01 | 0.98 ± 0.00 | 0.45 ± 0.02 | 0.515 ± 0.018 |
| 0.830 | 0.487 ± 0.018 | 0.339 ± 0.012 | 0.670 ± 0.024 | 0.277 ± 0.014 | 0.352 ± 0.014 | 0.70 ± 0.01 | 0.97 ± 0.01 | 0.42 ± 0.02 | 0.489 ± 0.018 |
| 0.891 | 0.474 ± 0.018 | 0.327 ± 0.012 | 0.649 ± 0.026 | 0.268 ± 0.014 | 0.343 ± 0.016 | 0.70 ± 0.01 | 0.97 ± 0.01 | 0.42 ± 0.02 | 0.476 ± 0.019 |
| 0.957 | 0.470 ± 0.018 | 0.325 ± 0.012 | 0.652 ± 0.026 | 0.251 ± 0.014 | 0.341 ± 0.015 | 0.70 ± 0.01 | 0.98 ± 0.00 | 0.40 ± 0.02 | 0.472 ± 0.018 |
| 1.028 | 0.458 ± 0.017 | 0.328 ± 0.012 | 0.638 ± 0.023 | 0.242 ± 0.012 | 0.349 ± 0.013 | 0.69 ± 0.01 | 0.98 ± 0.00 | 0.40 ± 0.02 | 0.460 ± 0.017 |
| 1.090 | 0.546 ± 0.018 | 0.418 ± 0.015 | 0.653 ± 0.026 | 0.385 ± 0.014 | 0.459 ± 0.020 | 0.80 ± 0.01 | 0.99 ± 0.01 | 0.60 ± 0.02 | 0.543 ± 0.019 |
| 1.143 | 0.443 ± 0.018 | 0.315 ± 0.012 | 0.639 ± 0.025 | 0.214 ± 0.011 | 0.347 ± 0.014 | 0.66 ± 0.01 | 0.98 ± 0.01 | 0.34 ± 0.02 | 0.444 ± 0.018 |
| 1.227 | 0.435 ± 0.017 | 0.314 ± 0.012 | 0.643 ± 0.023 | 0.197 ± 0.010 | 0.355 ± 0.014 | 0.64 ± 0.01 | 0.97 ± 0.01 | 0.30 ± 0.02 | 0.435 ± 0.016 |
| 1.317 | 0.424 ± 0.018 | 0.310 ± 0.013 | 0.638 ± 0.027 | 0.179 ± 0.011 | 0.333 ± 0.016 | 0.63 ± 0.01 | 0.98 ± 0.01 | 0.27 ± 0.02 | 0.426 ± 0.019 |
| 1.413 | 0.416 ± 0.017 | 0.317 ± 0.012 | 0.650 ± 0.026 | 0.166 ± 0.011 | 0.340 ± 0.013 | 0.62 ± 0.01 | 0.98 ± 0.01 | 0.24 ± 0.02 | 0.415 ± 0.018 |
| 1.517 | 0.403 ± 0.017 | 0.316 ± 0.012 | 0.646 ± 0.027 | 0.143 ± 0.010 | 0.342 ± 0.015 | 0.60 ± 0.01 | 0.97 ± 0.01 | 0.21 ± 0.02 | 0.402 ± 0.018 |
| 1.609 | 0.384 ± 0.018 | 0.305 ± 0.013 | 0.628 ± 0.028 | 0.118 ± 0.010 | 0.335 ± 0.015 | 0.58 ± 0.01 | 0.97 ± 0.01 | 0.18 ± 0.02 | 0.383 ± 0.018 |
| 1.690 | 0.378 ± 0.014 | 0.319 ± 0.011 | 0.651 ± 0.025 | 0.067 ± 0.006 | 0.347 ± 0.011 | 0.56 ± 0.01 | 0.97 ± 0.01 | 0.13 ± 0.01 | 0.377 ± 0.015 |
| 1.812 | 0.348 ± 0.014 | 0.299 ± 0.011 | 0.616 ± 0.025 | 0.047 ± 0.005 | 0.336 ± 0.010 | 0.54 ± 0.01 | 0.96 ± 0.00 | 0.11 ± 0.01 | 0.349 ± 0.014 |
| 1.943 | 0.356 ± 0.014 | 0.309 ± 0.011 | 0.635 ± 0.026 | 0.045 ± 0.006 | 0.342 ± 0.010 | 0.53 ± 0.01 | 0.96 ± 0.01 | 0.10 ± 0.01 | 0.352 ± 0.015 |
| 2.083 | 0.354 ± 0.014 | 0.317 ± 0.011 | 0.646 ± 0.026 | 0.039 ± 0.005 | 0.342 ± 0.011 | 0.53 ± 0.01 | 0.96 ± 0.01 | 0.09 ± 0.01 | 0.356 ± 0.015 |
| 2.234 | 0.349 ± 0.011 | 0.314 ± 0.009 | 0.649 ± 0.021 | 0.031 ± 0.006 | 0.340 ± 0.008 | 0.51 ± 0.01 | 0.96 ± 0.01 | 0.07 ± 0.01 | 0.353 ± 0.012 |
| 2.368 | 0.353 ± 0.013 | 0.319 ± 0.011 | 0.655 ± 0.025 | 0.031 ± 0.006 | 0.335 ± 0.011 | 0.51 ± 0.01 | 0.96 ± 0.01 | 0.06 ± 0.01 | 0.355 ± 0.013 |

T per slice (rows), inverse-variance mean over the 6 groups of each detector (columns D1/D2/D3), for A, C, D:

| slice | T_A D1 | T_A D2 | T_A D3 | T_C D1 | T_C D2 | T_C D3 | T_D D1 | T_D D2 | T_D D3 |
|---|---|---|---|---|---|---|---|---|---|
| BGS 0.05-0.15 | 0.19±0.07 | 0.20±0.06 | 0.16±0.05 | 0.31±0.11 | 0.30±0.11 | 0.26±0.11 | 0.06±0.03 | 0.09±0.03 | 0.05±0.02 |
| BGS 0.15-0.25 | 0.64±0.17 | 0.55±0.13 | 0.54±0.14 | 1.02±0.27 | 0.92±0.23 | 1.00±0.27 | 0.23±0.06 | 0.15±0.04 | 0.05±0.03 |
| BGS 0.25-0.35 | 0.44±0.11 | 0.44±0.11 | 0.37±0.10 | 0.63±0.17 | 0.59±0.16 | 0.64±0.17 | 0.22±0.05 | 0.26±0.05 | 0.06±0.02 |
| BGS 0.35-0.45 | 0.50±0.08 | 0.47±0.07 | 0.37±0.05 | 0.69±0.11 | 0.67±0.11 | 0.66±0.10 | 0.29±0.04 | 0.24±0.03 | 0.04±0.02 |
| LRG 0.40-0.50 | 0.53±0.08 | 0.47±0.07 | 0.43±0.06 | 0.78±0.10 | 0.72±0.10 | 0.77±0.10 | 0.14±0.04 | 0.08±0.03 | -0.07±0.03 |
| LRG 0.50-0.60 | 0.53±0.07 | 0.44±0.06 | 0.37±0.05 | 0.71±0.09 | 0.71±0.09 | 0.70±0.09 | 0.29±0.04 | 0.10±0.02 | -0.04±0.02 |
| LRG 0.60-0.70 | 0.75±0.07 | 0.59±0.06 | 0.44±0.04 | 0.90±0.08 | 0.90±0.09 | 0.85±0.08 | 0.57±0.06 | 0.18±0.03 | -0.04±0.02 |
| LRG 0.70-0.80 | 0.66±0.06 | 0.47±0.05 | 0.35±0.03 | 0.71±0.07 | 0.66±0.07 | 0.70±0.07 | 0.61±0.06 | 0.23±0.03 | -0.02±0.02 |
| LRG 0.80-0.90 | 0.67±0.11 | 0.63±0.10 | 0.54±0.09 | 0.88±0.16 | 0.95±0.17 | 0.89±0.16 | 0.44±0.07 | 0.30±0.04 | 0.17±0.03 |
| LRG 0.90-1.00 | 0.63±0.06 | 0.57±0.05 | 0.49±0.05 | 0.82±0.08 | 0.83±0.08 | 0.83±0.09 | 0.43±0.05 | 0.30±0.04 | 0.14±0.02 |
| LRG 1.00-1.10 | 0.56±0.06 | 0.43±0.05 | 0.40±0.05 | 0.71±0.08 | 0.65±0.08 | 0.74±0.09 | 0.39±0.06 | 0.19±0.04 | 0.04±0.02 |
| ELG 0.80-0.90 | 0.52±0.08 | 0.50±0.08 | 0.42±0.06 | 0.67±0.12 | 0.67±0.12 | 0.62±0.10 | 0.35±0.08 | 0.33±0.08 | 0.21±0.05 |
| ELG 0.90-1.00 | 0.66±0.09 | 0.64±0.09 | 0.54±0.07 | 0.86±0.12 | 0.88±0.12 | 0.89±0.13 | 0.47±0.07 | 0.41±0.08 | 0.20±0.04 |
| ELG 1.00-1.10 | 0.57±0.09 | 0.39±0.08 | 0.32±0.07 | 0.52±0.11 | 0.44±0.12 | 0.51±0.12 | 0.63±0.15 | 0.31±0.09 | 0.13±0.04 |
| ELG 1.10-1.20 | 0.59±0.11 | 0.55±0.10 | 0.48±0.08 | 0.74±0.14 | 0.82±0.15 | 0.73±0.12 | 0.41±0.12 | 0.26±0.08 | 0.19±0.05 |
| ELG 1.20-1.30 | 0.46±0.06 | 0.43±0.06 | 0.35±0.04 | 0.65±0.09 | 0.64±0.10 | 0.58±0.08 | 0.27±0.04 | 0.23±0.04 | 0.13±0.03 |
| ELG 1.30-1.40 | 0.52±0.06 | 0.41±0.06 | 0.33±0.06 | 0.70±0.11 | 0.62±0.13 | 0.59±0.13 | 0.33±0.05 | 0.18±0.05 | 0.06±0.03 |
| ELG 1.40-1.50 | 0.52±0.06 | 0.41±0.05 | 0.39±0.05 | 0.64±0.09 | 0.60±0.08 | 0.66±0.09 | 0.39±0.08 | 0.24±0.06 | 0.14±0.05 |
| ELG 1.50-1.60 | 0.53±0.07 | 0.44±0.06 | 0.36±0.05 | 0.68±0.10 | 0.66±0.10 | 0.64±0.10 | 0.37±0.07 | 0.21±0.05 | 0.07±0.04 |
| QSO 0.80-0.90 | 0.49±0.17 | 0.57±0.20 | 0.61±0.25 | 0.65±0.13 | 0.71±0.16 | 0.70±0.16 | -0.03±0.37 | 0.10±0.21 | 0.19±0.31 |
| QSO 0.90-1.00 | 0.69±0.18 | 0.57±0.17 | 0.53±0.14 | 0.71±0.19 | 0.62±0.22 | 0.70±0.19 | 0.60±0.25 | 0.46±0.18 | 0.24±0.14 |
| QSO 1.00-1.10 | 0.45±0.12 | 0.32±0.10 | 0.31±0.08 | 0.56±0.14 | 0.50±0.15 | 0.61±0.12 | 0.31±0.14 | 0.07±0.10 | -0.05±0.08 |
| QSO 1.10-1.20 | 0.56±0.13 | 0.56±0.14 | 0.51±0.12 | 0.79±0.16 | 0.90±0.18 | 0.81±0.16 | 0.20±0.12 | 0.04±0.10 | 0.05±0.08 |
| QSO 1.20-1.30 | 0.42±0.11 | 0.39±0.13 | 0.36±0.08 | 0.56±0.16 | 0.52±0.18 | 0.53±0.13 | 0.23±0.13 | 0.21±0.13 | 0.15±0.06 |
| QSO 1.30-1.40 | 0.15±0.09 | 0.16±0.08 | 0.18±0.08 | 0.26±0.13 | 0.32±0.13 | 0.33±0.14 | 0.00±0.10 | -0.01±0.08 | 0.00±0.06 |
| QSO 1.40-1.50 | 0.50±0.15 | 0.44±0.13 | 0.36±0.13 | 0.71±0.21 | 0.67±0.22 | 0.60±0.26 | 0.29±0.16 | 0.23±0.11 | 0.13±0.08 |
| QSO 1.50-1.60 | 0.33±0.08 | 0.18±0.09 | 0.19±0.07 | 0.42±0.12 | 0.29±0.14 | 0.38±0.12 | 0.23±0.08 | 0.07±0.07 | -0.04±0.07 |
| QSO 1.60-1.70 | 0.32±0.07 | 0.29±0.07 | 0.22±0.06 | 0.55±0.11 | 0.56±0.13 | 0.48±0.11 | 0.09±0.07 | 0.02±0.06 | -0.04±0.04 |
| QSO 1.70-1.80 | 0.25±0.10 | 0.23±0.08 | 0.22±0.07 | 0.43±0.11 | 0.46±0.11 | 0.47±0.11 | 0.06±0.11 | -0.02±0.08 | -0.05±0.06 |
| QSO 1.80-1.90 | 0.23±0.08 | 0.20±0.08 | 0.23±0.07 | 0.54±0.10 | 0.50±0.11 | 0.52±0.11 | -0.11±0.13 | -0.11±0.09 | -0.14±0.09 |
| QSO 1.90-2.00 | 0.25±0.09 | 0.22±0.11 | 0.25±0.06 | 0.49±0.14 | 0.43±0.19 | 0.44±0.13 | -0.04±0.10 | -0.07±0.14 | 0.03±0.07 |
| QSO 2.00-2.10 | 0.26±0.08 | 0.23±0.07 | 0.24±0.06 | 0.43±0.12 | 0.42±0.13 | 0.46±0.12 | 0.07±0.06 | 0.02±0.05 | 0.00±0.04 |

Fraction of (slice, group) cells with σ(T_C) > 0.2: 0.12; the correction uses T per cell.

Clustering bias of each slice on this subset (b_r, jackknife): BGS 0.05 1.05±0.04, BGS 0.15 1.18±0.03, BGS 0.25 1.31±0.03, BGS 0.35 1.42±0.03, LRG 0.40 1.86±0.04, LRG 0.50 2.01±0.03, LRG 0.60 2.13±0.04, LRG 0.70 2.19±0.04, LRG 0.80 2.31±0.04, LRG 0.90 2.33±0.04, LRG 1.00 2.34±0.06, ELG 0.80 1.08±0.04, ELG 0.90 1.12±0.04, ELG 1.00 1.28±0.03, ELG 1.10 1.35±0.03, ELG 1.20 1.47±0.05, ELG 1.30 1.48±0.05, ELG 1.40 1.63±0.04, ELG 1.50 1.80±0.08, QSO 0.80 1.80±0.20, QSO 0.90 1.54±0.21, QSO 1.00 1.77±0.17, QSO 1.10 1.98±0.16, QSO 1.20 2.14±0.12, QSO 1.30 2.21±0.17, QSO 1.40 2.32±0.13, QSO 1.50 2.47±0.15, QSO 1.60 2.51±0.17, QSO 1.70 2.51±0.21, QSO 1.80 2.71±0.21, QSO 1.90 2.91±0.21, QSO 2.00 3.04±0.19.

Corrected b_I d(νI_ν)/dz per group (nW m^-2 sr^-1 per unit z), inverse-variance mean over slices (jackknife σ, which carries the transfer's own error): uncorrected (w̄ / b_r w̄_m, no transfer), corrected by T_A (fiducial), by T_C (all smooth: the light outside catalogued sources), by T_D (all compact); z < 0.5 and z ≥ 0.5 for T_C:

| group µm | uncorrected | /T_A | /T_C | /T_D | /T_C z < 0.5 | /T_C z ≥ 0.5 | S/N (T_C) |
|---|---|---|---|---|---|---|---|
| 0.773 | 3.96 ± 0.89 | 7.3 ± 1.7 | 5.5 ± 1.3 | 10.3 ± 2.4 | 6.9 ± 2.2 | 4.9 ± 1.5 | 4.2 |
| 0.830 | 2.78 ± 0.98 | 5.0 ± 1.7 | 3.9 ± 1.4 | 6.8 ± 2.4 | 5.8 ± 2.2 | 2.8 ± 1.6 | 2.8 |
| 0.891 | 3.93 ± 0.86 | 7.4 ± 1.6 | 5.6 ± 1.3 | 9.5 ± 2.3 | 7.6 ± 2.4 | 4.8 ± 1.7 | 4.4 |
| 0.957 | 4.30 ± 0.73 | 8.1 ± 1.7 | 6.3 ± 1.2 | 9.9 ± 2.6 | 7.7 ± 2.6 | 5.7 ± 1.5 | 5.2 |
| 1.028 | 3.21 ± 0.90 | 5.8 ± 2.0 | 4.5 ± 1.4 | 6.9 ± 3.3 | 5.6 ± 2.8 | 4.0 ± 2.0 | 3.1 |
| 1.090 | 30.26 ± 12.16 | 55.2 ± 22.3 | 47.8 ± 18.1 | 60.5 ± 27.6 | 40.8 ± 28.5 | 50.8 ± 23.7 | 2.6 |
| 1.143 | 3.78 ± 0.77 | 8.2 ± 1.7 | 6.0 ± 1.2 | 11.0 ± 2.6 | 8.1 ± 2.1 | 5.4 ± 1.3 | 5.0 |
| 1.227 | 4.21 ± 0.84 | 9.4 ± 1.7 | 6.6 ± 1.2 | 13.4 ± 2.8 | 9.6 ± 2.0 | 5.8 ± 1.4 | 5.4 |
| 1.317 | 4.61 ± 0.69 | 10.6 ± 1.6 | 7.3 ± 1.1 | 16.5 ± 2.8 | 7.0 ± 1.7 | 7.3 ± 1.3 | 6.5 |
| 1.413 | 4.20 ± 0.53 | 9.4 ± 1.2 | 6.4 ± 0.8 | 14.3 ± 2.2 | 7.1 ± 1.6 | 6.2 ± 0.9 | 7.8 |
| 1.517 | 3.70 ± 0.54 | 7.6 ± 1.3 | 5.1 ± 0.9 | 10.7 ± 2.9 | 7.3 ± 1.7 | 4.5 ± 1.1 | 5.8 |
| 1.609 | 3.28 ± 0.58 | 7.6 ± 1.5 | 5.0 ± 1.0 | 12.1 ± 3.9 | 6.7 ± 1.6 | 4.5 ± 1.2 | 5.1 |
| 1.690 | 3.76 ± 0.48 | 8.6 ± 1.1 | 5.4 ± 0.7 | 13.3 ± 3.8 | 3.9 ± 1.6 | 5.9 ± 0.9 | 8.1 |
| 1.812 | 3.18 ± 0.30 | 8.0 ± 0.9 | 5.0 ± 0.5 | 11.1 ± 3.1 | 4.8 ± 1.1 | 5.1 ± 0.7 | 9.1 |
| 1.943 | 3.45 ± 0.20 | 9.5 ± 0.7 | 5.6 ± 0.4 | 18.0 ± 3.2 | 4.8 ± 1.0 | 6.0 ± 0.5 | 15.4 |
| 2.083 | 3.93 ± 0.33 | 10.8 ± 1.1 | 6.2 ± 0.6 | 13.8 ± 4.9 | 4.8 ± 1.2 | 6.7 ± 0.7 | 10.2 |
| 2.234 | 3.66 ± 0.32 | 9.5 ± 1.0 | 5.3 ± 0.5 | 14.3 ± 5.9 | 2.9 ± 0.9 | 6.3 ± 0.7 | 10.2 |
| 2.368 | 3.05 ± 0.23 | 7.9 ± 0.7 | 4.5 ± 0.4 | 11.4 ± 3.5 | 2.9 ± 0.7 | 5.5 ± 0.5 | 11.4 |

Per slice b_I d(νI_ν)/dz corrected by T_C (the light outside catalogued sources), nW m^-2 sr^-1 per unit z, value ± jackknife σ:

| slice | 0.773 | 0.830 | 0.891 | 0.957 | 1.028 | 1.090 | 1.143 | 1.227 | 1.317 | 1.413 | 1.517 | 1.609 | 1.690 | 1.812 | 1.943 | 2.083 | 2.234 | 2.368 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BGS 0.05-0.15 | 1±11 | 5±13 | 4±15 | 1±11 | 3±12 | 55±126 | 2±7 | 1±7 | 2±5 | -1±4 | -7±6 | -6±7 | 9±8 | 6±6 | -5±4 | -3±3 | -1±3 | 1±2 |
| BGS 0.15-0.25 | 10±4 | 10±5 | 12±5 | 11±4 | 5±4 | 120±49 | 9±4 | 9±4 | 8±4 | 8±3 | 8±3 | 7±3 | 4±2 | 5±2 | 5±2 | 4±1 | 2±1 | 2±1 |
| BGS 0.25-0.35 | 6±4 | 7±4 | 5±4 | 9±4 | 8±5 | -19±62 | 10±4 | 11±4 | 11±3 | 10±3 | 7±3 | 9±4 | 4±3 | 4±2 | 4±2 | 8±6 | 2±2 | 3±1 |
| BGS 0.35-0.45 | 8±5 | 5±6 | 7±6 | 7±6 | 3±7 | 9±91 | 7±4 | 9±5 | 4±4 | 6±4 | 9±4 | 6±4 | 2±4 | 4±2 | 6±2 | 7±3 | 5±2 | 4±1 |
| LRG 0.40-0.50 | 6±4 | 3±4 | 7±4 | 5±4 | 5±5 | -28±71 | 8±4 | 12±4 | 6±3 | 8±2 | 10±3 | 9±4 | 5±3 | 6±2 | 6±1 | 6±2 | 6±2 | 5±2 |
| LRG 0.50-0.60 | 18±5 | 20±6 | 19±5 | 22±5 | 17±6 | 59±58 | 15±3 | 13±4 | 13±4 | 15±3 | 12±3 | 9±3 | 8±3 | 8±2 | 7±2 | 9±3 | 6±2 | 7±2 |
| LRG 0.60-0.70 | 6±3 | 2±4 | 5±3 | 10±4 | 13±4 | 100±52 | 9±2 | 7±3 | 8±3 | 9±2 | 7±2 | 7±3 | 7±2 | 9±4 | 8±2 | 7±2 | 6±1 | 6±1 |
| LRG 0.70-0.80 | 15±5 | 12±5 | 13±5 | 10±5 | 5±6 | 72±72 | 8±4 | 9±3 | 11±3 | 8±3 | 11±3 | 8±4 | 8±4 | 8±2 | 5±2 | 7±2 | 8±2 | 9±2 |
| LRG 0.80-0.90 | 2±3 | -0±4 | 1±5 | 1±5 | -1±5 | -30±50 | 4±3 | 6±3 | 6±3 | 6±2 | 2±2 | 4±2 | 5±2 | 5±2 | 7±2 | 9±2 | 10±3 | 6±2 |
| LRG 0.90-1.00 | 5±4 | 5±5 | 4±4 | 5±4 | -1±4 | 141±76 | 2±3 | 3±3 | 7±3 | 7±2 | 4±3 | 5±4 | 4±3 | 4±2 | 7±1 | 8±2 | 8±2 | 6±1 |
| LRG 1.00-1.10 | 4±5 | 0±7 | 7±6 | 9±7 | -1±8 | 48±91 | 3±5 | 10±5 | 13±5 | 8±3 | 8±4 | 8±5 | 12±5 | 12±5 | 11±3 | 9±2 | 9±2 | 7±2 |
| ELG 0.80-0.90 | -12±9 | -8±8 | -10±10 | -10±8 | -12±11 | -118±124 | 9±7 | 13±7 | 6±5 | 3±4 | -3±5 | 3±6 | 9±6 | 6±5 | 7±2 | 8±3 | 9±3 | 3±2 |
| ELG 0.90-1.00 | -7±7 | 5±6 | -1±7 | -3±7 | -3±8 | 90±107 | -1±4 | 0±4 | 2±4 | 2±2 | -0±3 | 1±4 | 6±4 | 5±2 | 4±2 | 1±2 | 2±2 | 3±2 |
| ELG 1.00-1.10 | 1±13 | -14±14 | -3±14 | -8±11 | -6±13 | -16±128 | -7±9 | -3±8 | 9±12 | 6±6 | 10±7 | 9±8 | 8±9 | 4±4 | 13±7 | 10±3 | 9±4 | 8±4 |
| ELG 1.10-1.20 | 9±7 | 1±8 | 8±11 | 6±8 | 12±9 | -12±109 | 10±6 | 11±6 | 10±5 | 7±4 | -0±4 | 2±4 | 3±4 | 3±3 | 6±3 | 6±4 | 6±3 | 5±3 |
| ELG 1.20-1.30 | 6±7 | 0±10 | 11±9 | 5±8 | 2±9 | 115±104 | 4±6 | 7±5 | 8±6 | 10±5 | 10±6 | 3±4 | 4±6 | 5±4 | 9±4 | 11±4 | 9±4 | 5±3 |
| ELG 1.30-1.40 | 8±7 | -7±7 | -8±10 | -7±9 | -5±15 | 120±108 | 2±5 | 5±7 | 7±5 | 5±5 | -4±6 | -2±5 | 2±5 | -2±5 | 5±3 | 11±5 | 7±4 | 6±4 |
| ELG 1.40-1.50 | -15±11 | -9±9 | -12±7 | -12±9 | -13±13 | 138±97 | 5±7 | 3±6 | -2±5 | -3±4 | -2±7 | -7±8 | 2±4 | 2±3 | 4±4 | -1±3 | -2±4 | 2±2 |
| ELG 1.50-1.60 | 6±11 | 12±10 | 13±10 | 5±9 | 3±12 | 235±160 | 5±5 | 10±7 | 8±6 | 2±6 | 3±6 | 4±7 | 10±6 | 8±5 | 5±3 | 7±6 | 3±4 | 3±3 |
| QSO 0.80-0.90 | 11±10 | 4±14 | 11±11 | 17±12 | 2±13 | -93±215 | 8±9 | 14±7 | 15±8 | 12±7 | 9±7 | 22±11 | 3±6 | 6±6 | 7±4 | 4±5 | 9±5 | 7±4 |
| QSO 0.90-1.00 | 3±12 | 12±14 | 4±14 | 5±12 | -3±18 | 183±210 | 10±8 | 20±10 | 18±13 | 10±7 | 5±7 | 16±16 | 3±5 | -1±6 | 5±5 | 17±7 | 5±5 | 14±6 |
| QSO 1.00-1.10 | 14±19 | -19±19 | 17±26 | 15±23 | 6±17 | -267±220 | 7±11 | 11±16 | 20±16 | 17±12 | 31±17 | 19±18 | 13±10 | 5±10 | 8±3 | 15±6 | 15±5 | 8±4 |
| QSO 1.10-1.20 | -3±10 | -2±14 | -6±12 | 1±11 | 9±12 | -27±113 | 10±8 | 11±8 | 7±7 | 7±5 | 1±5 | 1±5 | 7±6 | -1±4 | -0±5 | 6±3 | 3±3 | 3±3 |
| QSO 1.20-1.30 | -13±11 | -28±17 | -14±14 | -17±15 | 3±17 | -185±233 | -3±8 | -11±9 | -16±18 | -11±13 | 10±8 | -5±8 | 5±7 | -3±8 | -1±5 | 3±5 | 6±5 | 3±3 |
| QSO 1.30-1.40 | 62±104 | -29±25 | 3±21 | -3±40 | 61±49 | -185±395 | 8±14 | 8±17 | 5±17 | 7±15 | -10±16 | 4±23 | -12±22 | -12±9 | 10±8 | 6±11 | 28±57 | 14±9 |
| QSO 1.40-1.50 | -13±11 | -12±9 | -4±10 | -12±13 | -10±12 | 20±160 | -9±11 | -10±6 | -4±6 | -7±7 | -5±5 | 1±6 | 8±7 | 1±5 | 4±9 | 0±6 | 0±3 | -2±3 |
| QSO 1.50-1.60 | 19±15 | 29±23 | 3±19 | 8±17 | -6±25 | 145±321 | 10±11 | 3±16 | 12±21 | 7±16 | -16±17 | -3±35 | 29±17 | 20±12 | 13±7 | 4±13 | 6±6 | 8±4 |
| QSO 1.60-1.70 | 10±13 | 14±16 | 18±18 | 8±15 | 19±17 | -61±221 | 3±9 | -7±7 | 9±8 | 8±7 | 6±8 | 12±7 | 12±9 | 9±6 | 5±4 | 4±5 | 7±4 | 7±4 |
| QSO 1.70-1.80 | 18±13 | 18±17 | 2±14 | 17±18 | 1±24 | -84±643 | -2±8 | -8±10 | 1±11 | -1±8 | 7±8 | -5±17 | -0±9 | 6±6 | 4±4 | 5±6 | 3±5 | 10±6 |
| QSO 1.80-1.90 | -12±13 | -5±15 | -8±15 | 7±11 | 3±19 | 158±219 | -5±9 | 1±8 | -0±8 | 8±6 | -1±8 | 7±8 | 10±6 | 4±6 | 3±3 | 2±5 | -3±6 | 0±5 |
| QSO 1.90-2.00 | 10±11 | 17±14 | 16±17 | 24±18 | 2±16 | 31±295 | 5±13 | 7±12 | 8±20 | 11±8 | -2±10 | -10±9 | -1±8 | 3±6 | 4±7 | 5±6 | 3±6 | 2±5 |
| QSO 2.00-2.10 | -18±23 | 17±22 | 7±19 | -9±23 | -6±33 | 110±377 | -8±11 | 1±11 | -5±12 | -3±13 | 19±12 | 25±19 | -1±9 | -1±11 | 2±4 | 5±5 | 0±6 | 9±5 |

Stellar null (12 < G < 16 stars as centres, R-138's estimator and units with T_C), χ² against zero over 32 slices (jackknife σ), without and with the fixed-flux star mask; the galaxy signal with the star mask; the galaxy-star cross-correlation w̄_gs per slice and the contamination it implies, max over slices of |w̄_gs N_* S| / σ(signal):

| group µm | star null χ² (R-138 masks) | star null χ² (star mask) | signal /T_C, star mask | max contamination / σ |
|---|---|---|---|---|
| 0.773 | 16/32 (p=9.9e-01) | 41/32 (p=1.4e-01) | 5.1 ± 1.2 | 29.80 |
| 0.830 | 18/32 (p=9.8e-01) | 41/32 (p=1.3e-01) | 3.6 ± 1.2 | 21.18 |
| 0.891 | 48/32 (p=3.2e-02) | 26/32 (p=7.7e-01) | 5.1 ± 1.2 | 14.11 |
| 0.957 | 71/32 (p=7.8e-05) | 23/32 (p=8.7e-01) | 5.7 ± 1.1 | 11.94 |
| 1.028 | 58/32 (p=3.7e-03) | 50/32 (p=2.2e-02) | 3.7 ± 1.4 | 36.18 |
| 1.090 | 64/32 (p=7.5e-04) | 57/32 (p=4.2e-03) | 46.4 ± 18.1 | 50.17 |
| 1.143 | 131/32 (p=7.4e-14) | 60/32 (p=1.9e-03) | 6.0 ± 1.1 | 1.66 |
| 1.227 | 174/32 (p=2.0e-21) | 90/32 (p=1.8e-07) | 6.5 ± 1.2 | 7.48 |
| 1.317 | 85/32 (p=9.9e-07) | 33/32 (p=4.2e-01) | 7.1 ± 1.1 | 3.31 |
| 1.413 | 109/32 (p=3.0e-10) | 44/32 (p=7.0e-02) | 6.1 ± 0.8 | 7.00 |
| 1.517 | 106/32 (p=8.2e-10) | 45/32 (p=5.8e-02) | 4.9 ± 0.9 | 9.52 |
| 1.609 | 164/32 (p=1.4e-19) | 107/32 (p=4.6e-10) | 4.9 ± 1.0 | 8.96 |
| 1.690 | 52/32 (p=1.5e-02) | 19/32 (p=9.7e-01) | 5.2 ± 0.6 | 6.60 |
| 1.812 | 51/32 (p=1.8e-02) | 14/32 (p=1.0e+00) | 4.9 ± 0.5 | 8.06 |
| 1.943 | 47/32 (p=4.0e-02) | 21/32 (p=9.3e-01) | 5.5 ± 0.4 | 8.82 |
| 2.083 | 37/32 (p=2.4e-01) | 19/32 (p=9.6e-01) | 6.1 ± 0.6 | 5.50 |
| 2.234 | 20/32 (p=9.6e-01) | 38/32 (p=2.1e-01) | 5.3 ± 0.5 | 12.84 |
| 2.368 | 18/32 (p=9.8e-01) | 19/32 (p=9.6e-01) | 4.5 ± 0.4 | 14.86 |

w̄_gs per slice (galaxies around 12 < G < 16 stars, Hamilton estimator, same annulus): BGS 0.05 +0.0094±0.0090, BGS 0.15 -0.0007±0.0042, BGS 0.25 +0.0036±0.0045, BGS 0.35 -0.0079±0.0040, LRG 0.40 -0.0080±0.0050, LRG 0.50 +0.0018±0.0037, LRG 0.60 +0.0064±0.0053, LRG 0.70 -0.0042±0.0041, LRG 0.80 +0.0070±0.0036, LRG 0.90 +0.0011±0.0045, LRG 1.00 -0.0047±0.0049, ELG 0.80 +0.0034±0.0028, ELG 0.90 +0.0000±0.0027, ELG 1.00 -0.0042±0.0027, ELG 1.10 -0.0005±0.0030, ELG 1.20 -0.0007±0.0033, ELG 1.30 -0.0024±0.0034, ELG 1.40 +0.0013±0.0044, ELG 1.50 -0.0088±0.0056, QSO 0.80 +0.0034±0.0083, QSO 0.90 +0.0069±0.0061, QSO 1.00 -0.0127±0.0051, QSO 1.10 +0.0013±0.0061, QSO 1.20 -0.0036±0.0056, QSO 1.30 +0.0064±0.0063, QSO 1.40 -0.0060±0.0074, QSO 1.50 +0.0010±0.0072, QSO 1.60 -0.0035±0.0067, QSO 1.70 -0.0030±0.0083, QSO 1.80 -0.0137±0.0076, QSO 1.90 +0.0110±0.0089, QSO 2.00 +0.0097±0.0067. Stars per sr 1.23e+06.

Star-density split (patches at the median 12 < G < 16 count): T_A per half and the predicted low/high ratio T_A(low)/T_A(high); the observed low/high ratio of the uncorrected signal on this subset, and R-138's:

| group µm | T_A low-star | T_A high-star | predicted ratio | observed ratio (this subset) | R-138 ratio |
|---|---|---|---|---|---|
| 0.773 | 0.435 ± 0.067 | 0.536 ± 0.015 | 0.81 ± 0.13 | 1.31 ± 0.59 | 2.22 |
| 0.830 | 0.390 ± 0.048 | 0.521 ± 0.015 | 0.75 ± 0.09 | 4.18 ± 2.81 | 2.55 |
| 0.891 | 0.384 ± 0.050 | 0.504 ± 0.014 | 0.76 ± 0.10 | 2.39 ± 0.93 | 3.00 |
| 0.957 | 0.367 ± 0.052 | 0.503 ± 0.015 | 0.73 ± 0.11 | 1.95 ± 0.65 | 2.57 |
| 1.028 | 0.350 ± 0.053 | 0.494 ± 0.014 | 0.71 ± 0.11 | 1.57 ± 1.09 | 2.21 |
| 1.090 | 0.462 ± 0.063 | 0.570 ± 0.016 | 0.81 ± 0.11 | 1.13 ± 0.70 | 0.66 |
| 1.143 | 0.344 ± 0.043 | 0.466 ± 0.017 | 0.74 ± 0.10 | 1.51 ± 0.56 | 1.43 |
| 1.227 | 0.370 ± 0.043 | 0.456 ± 0.017 | 0.81 ± 0.10 | 1.27 ± 0.50 | 1.50 |
| 1.317 | 0.388 ± 0.042 | 0.440 ± 0.018 | 0.88 ± 0.10 | 1.32 ± 0.37 | 0.96 |
| 1.413 | 0.382 ± 0.047 | 0.435 ± 0.016 | 0.88 ± 0.11 | 1.28 ± 0.30 | 1.27 |
| 1.517 | 0.369 ± 0.035 | 0.424 ± 0.018 | 0.87 ± 0.09 | 1.23 ± 0.37 | 1.65 |
| 1.609 | 0.355 ± 0.037 | 0.403 ± 0.018 | 0.88 ± 0.10 | 0.76 ± 0.40 | 1.15 |
| 1.690 | 0.348 ± 0.048 | 0.402 ± 0.012 | 0.87 ± 0.12 | 1.63 ± 0.46 | 1.37 |
| 1.812 | 0.313 ± 0.037 | 0.372 ± 0.013 | 0.84 ± 0.10 | 1.54 ± 0.46 | 1.26 |
| 1.943 | 0.312 ± 0.031 | 0.375 ± 0.013 | 0.83 ± 0.09 | 0.99 ± 0.17 | 1.03 |
| 2.083 | 0.369 ± 0.042 | 0.366 ± 0.012 | 1.01 ± 0.12 | 0.60 ± 0.11 | 0.84 |
| 2.234 | 0.304 ± 0.040 | 0.370 ± 0.009 | 0.82 ± 0.11 | 0.66 ± 0.11 | 0.83 |
| 2.368 | 0.322 ± 0.037 | 0.368 ± 0.009 | 0.87 ± 0.10 | 0.89 ± 0.15 | 1.19 |

CIBER comparison (arXiv:2608.12116, read off its Fig. a2h_vs_redshift as R-138): this subset, inverse-variance mean of the slices whose centre lies in the CIBER bin, corrected by T_C and by T_D (the range for the total):

| λ | z | CIBER | group µm | /T_C | /T_D |
|---|---|---|---|---|---|
| 1.1 | 0.0-0.2 | 100 ± 30 | 1.090 | 54.8 ± 126.5 | 454.3 ± 2916.0 |
| 1.1 | 0.2-0.4 | 90 ± 18 | 1.090 | 65.7 ± 41.1 | 61.2 ± 214.1 |
| 1.1 | 0.4-0.6 | 23 ± 9 | 1.090 | 21.3 ± 41.5 | 70.3 ± 96.2 |
| 1.1 | 0.6-0.8 | < 5.5 (95%) | 1.090 | 90.2 ± 37.0 | 101.0 ± 43.9 |
| 1.1 | 0.8-1.0 | < 4.0 (95%) | 1.090 | 19.1 ± 40.6 | 29.8 ± 54.4 |
| 1.1 | 0.0-0.2 | 100 ± 30 | 1.143 | 2.2 ± 7.2 | 7.5 ± 24.8 |
| 1.1 | 0.2-0.4 | 90 ± 18 | 1.143 | 9.3 ± 2.6 | 21.6 ± 6.7 |
| 1.1 | 0.4-0.6 | 23 ± 9 | 1.143 | 10.7 ± 2.7 | 30.7 ± 9.9 |
| 1.1 | 0.6-0.8 | < 5.5 (95%) | 1.143 | 8.4 ± 1.8 | 16.4 ± 4.2 |
| 1.1 | 0.8-1.0 | < 4.0 (95%) | 1.143 | 3.2 ± 2.5 | 6.3 ± 5.2 |
| 1.8 | 0.0-0.2 | 140 ± 40 | 1.812 | 5.9 ± 6.0 | 40.0 ± 56.0 |
| 1.8 | 0.2-0.4 | 63 ± 15 | 1.812 | 4.4 ± 1.4 | 38.7 ± 18.2 |
| 1.8 | 0.4-0.6 | 35 ± 8 | 1.812 | 6.4 ± 1.6 | 4.5 ± 28.0 |
| 1.8 | 0.6-0.8 | < 14 (95%) | 1.812 | 8.6 ± 2.0 | -266.6 ± 220.8 |
| 1.8 | 0.8-1.0 | < 7.5 (95%) | 1.812 | 4.6 ± 1.3 | 12.6 ± 4.6 |

Transfer stability (the rule fixed at 01:30Z) and the extra null statistics:

| group µm | T_B/T_A | T_E/T_A | stable (both within 15 %) | T_D/T_C |
|---|---|---|---|---|
| 0.773 | 0.68 | 0.71 | no | 0.44 |
| 0.830 | 0.70 | 0.72 | no | 0.41 |
| 0.891 | 0.69 | 0.72 | no | 0.41 |
| 0.957 | 0.69 | 0.73 | no | 0.39 |
| 1.028 | 0.72 | 0.76 | no | 0.38 |
| 1.090 | 0.76 | 0.84 | no | 0.59 |
| 1.143 | 0.71 | 0.78 | no | 0.33 |
| 1.227 | 0.72 | 0.82 | no | 0.31 |
| 1.317 | 0.73 | 0.79 | no | 0.28 |
| 1.413 | 0.76 | 0.82 | no | 0.26 |
| 1.517 | 0.78 | 0.85 | no | 0.22 |
| 1.609 | 0.79 | 0.87 | no | 0.19 |
| 1.690 | 0.84 | 0.92 | no | 0.10 |
| 1.812 | 0.86 | 0.97 | yes | 0.08 |
| 1.943 | 0.87 | 0.96 | yes | 0.07 |
| 2.083 | 0.89 | 0.96 | yes | 0.06 |
| 2.234 | 0.90 | 0.97 | yes | 0.05 |
| 2.368 | 0.90 | 0.95 | yes | 0.05 |

w̄_gs against zero over 32 slices: χ² = 39.6/32 (p = 0.17); inverse-variance mean -0.0006 ± 0.0008

Star mask minus R-138 masks, T_C-corrected signal (nW m^-2 sr^-1 per unit z), inverse-variance mean over slices, jackknife of the difference; χ² of the per-slice differences against zero:
0.773: -0.37 ± 0.30; χ² 22/32 (p=0.90)
0.830: -0.32 ± 0.32; χ² 25/32 (p=0.79)
0.891: -0.43 ± 0.29; χ² 25/32 (p=0.79)
0.957: -0.44 ± 0.32; χ² 30/32 (p=0.57)
1.028: -0.68 ± 0.35; χ² 39/32 (p=0.18)
1.090: -2.27 ± 2.28; χ² 22/32 (p=0.91)
1.143: -0.03 ± 0.13; χ² 12/32 (p=1.00)
1.227: -0.09 ± 0.13; χ² 20/32 (p=0.96)
1.317: -0.10 ± 0.13; χ² 18/32 (p=0.97)
1.413: -0.21 ± 0.16; χ² 25/32 (p=0.81)
1.517: -0.07 ± 0.13; χ² 29/32 (p=0.61)
1.609: -0.08 ± 0.10; χ² 30/32 (p=0.54)
1.690: -0.13 ± 0.08; χ² 22/32 (p=0.91)
1.812: -0.10 ± 0.07; χ² 25/32 (p=0.79)
1.943: -0.08 ± 0.05; χ² 31/32 (p=0.52)
2.083: -0.05 ± 0.05; χ² 31/32 (p=0.51)
2.234: -0.04 ± 0.05; χ² 34/32 (p=0.38)
2.368: -0.03 ± 0.05; χ² 24/32 (p=0.84)

### What the R-143 table says

- **The transfer, measured.** For smooth galaxy-correlated light (σ = 1 cMpc) the SOURCE mask and clip cost almost
  nothing (the mask alone keeps 0.96-0.99) and the strip median removes a third: T_C = 0.62-0.69 in every group, flat in
  wavelength, 0.26-0.31 in the lowest slice (BGS 0.05-0.15) and 0.4-0.9 elsewhere (jackknife σ 0.02-0.03 per group). For
  the galaxies' own compact light the SOURCE mask is what matters: the mask alone keeps 0.45 at 0.77 µm falling to 0.06 at
  2.37 µm (the Reference Catalog is deeper at longer wavelength), and T_D = 0.30 → 0.03. The fiducial half/half mock
  has T_A = 0.51 → 0.35. So R-138's referee was right in kind: the chain removes galaxy-correlated light the strip
  transfer did not model, by a factor that depends on wavelength and on what the light is.
- **Stability rule (fixed 01:30Z): fails below 1.8 µm, holds at 1.81-2.37 µm.** At five times the amplitude the 5σ clip
  starts removing the compact mock (T_B/T_A = 0.68-0.79 below 1.6 µm), and the rising-z mock changes T by 15-29 % there
  (T_E/T_A = 0.71-0.87); at 1.81-2.37 µm both stay within 15 % (0.86-0.97), because the SOURCE mask has already removed
  the compact light there. T_D/T_C < 0.5 in every group but the filter-edge 1.090 µm one (0.59).
- **Hence, by the rule:** at 1.81-2.37 µm the amplitude is a measurement of the galaxy-correlated light outside
  catalogued sources, b_I d(νI_ν)/dz = 5.0 ± 0.5, 5.6 ± 0.4, 6.2 ± 0.6, 5.3 ± 0.5, 4.5 ± 0.4 nW m^-2 sr^-1 per unit z
  (9-15σ; the jackknife carries the transfer's own error); below 1.8 µm the same quantity (3.9-7.3, 2.6-8.1σ; the filter-edge 1.090 µm group 48 ± 18) is stated
  as a bound, not a measurement. The total galaxy-correlated light (including what the catalogue mask removes) is only
  bounded from below, by the same numbers; /T_D (7-18) shows how far above it could lie and is not an upper bound.
- **Shape:** after the per-slice correction the z dependence survives: at λ > 1.8 µm the z < 0.5 slices sit at 2.9-4.8
  and z ≥ 0.5 at 5.1-6.7; the lowest slice is consistent with zero in every group (its T_C is 0.26-0.31, so its
  corrected σ is large). Detection per group (T_C-corrected, this subset): 2.6-15σ, all 18 groups > 2.5σ.
- **Stellar null, star mask fixed:** with R-138's masks the 12 < G < 16 star null fails in 11 of 18 groups on this subset
  (χ² up to 174/32); with the fixed-flux G < 16 mask (31 % of blocks) it passes (p > 0.01) in 14 of 18 and fails at
  1.090, 1.143, 1.227 and 1.609 µm (p = 4e-3 to 5e-10), the detector-edge groups of D1/D2. The galaxy signal does not
  move when the star mask is added: the change is -0.03 to -0.68 nW per unit z, within 2σ in every group (χ² of the
  per-slice differences p = 0.18-1.00). The galaxy-star cross-correlation is zero: w̄_gs = -0.0006 ± 0.0008, χ² 39.6/32
  (p = 0.17). The pre-registered contamination column (|w̄_gs N_* S| / σ, 1.7-50) is not informative: it multiplies
  w̄_gs's noise (±0.003-0.009 per slice) by the ~2,000 stars inside θ_max and assumes the star profile is flat out to θ_max;
  the star-mask test above is the direct one.
- **Star-density split, predicted:** the transfer is lower, not higher, in the low-star patches (T_A low/high =
  0.71-1.01), so the transfer predicts a low/high signal ratio ≤ 1. Observed on this subset 1.3-4.2 ± 0.6-2.8 at 0.77-1.03 µm
  (R-138: 2.2-3.0), 1.9σ above the prediction at best in any one group; ~1 at λ > 1.9 µm. The transfer does not explain
  the split; with 25 patches per half and jackknife σ of 30-70 % on the ratio, the split is not significant either. Open.
- **CIBER:** the light outside catalogued sources at 1.81 µm (4.4 ± 1.4, 6.4 ± 1.6 at z 0.2-0.4, 0.4-0.6) is below CIBER
  x DESI-LS (63 ± 15, 35 ± 8), but it is a lower bound on the total, and dividing by T_D instead (38.7 ± 18.2, 4.5 ±
  28.0) reaches CIBER's values: the SPHEREx catalogue mask removes enough light to account for the difference. R-138's
  "4-25x below CIBER" is withdrawn as a tension. At 1.143 µm, z 0.6-0.8, the lower bound 8.4 ± 1.8 sits above CIBER's
  95 % limit (< 5.5), 1.6σ; at 1.8 µm z > 0.6 both agree.
- **Not redone here:** the zodi, half-survey and shuffled nulls and the line model (R-138 above; the line limits scale
  by 1/T like the continuum). The χ² tests use diagonal jackknife σ, so no Hartlap factor applies; the 51-patch
  covariance is not inverted anywhere.

## Blind review of the R-143 revision (finding check, 2026-11-02, context-free Opus, gate prompt verbatim, the section alone (`briefs/r143_finding_section.md`); `blind-review.md`)

"VERDICT: novel but not sound: the reported 1.8–2.4 µm amplitudes depend on the unvaried compact/smooth split and smooth-light profile of the injected mock, which the pre-set stability test never probes, so "light outside catalogued sources" is not a mock-independent measurement."

Strongest objection: the corrected signal is T_C·S_smooth + T_D·S_compact divided by T_C, so leaked compact light (profile and PSF wings past the mask radius, faint satellites) with transfers between T_D and T_C enters the 1.8-2.4 µm amplitude at a level set by the mock's half/half split and 1-cMpc kernel, which the A/B and A/E tests hold fixed; it also notes the strip median's scale dependence (the z < 0.5 deficit sits where T is lowest), the diagonal-only χ², the star mask chosen with the null in view, and the CIBER reconciliation resting on 39 ± 18 and 5 ± 28. Novelty held: no SPHEREx x galaxy measurement found (closest: 2608.12116 CIBER, 2109.10914 forecast, 1810.00885 method, 2608.03415, 2305.08847, 2601.21111). What survives, per the referee: the detection in all 18 groups (2.6-15σ). **A miss on physics: parked (brief: one revision; park on any further physics point).**
