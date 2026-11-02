# Results (R-130, cosmo-deep-10, 2026-11-01)

Archival strongly lensed supernovae behind the known galaxy-scale strong lenses, by forced photometry on the public ZTF
difference images at every lens position, below the alert threshold, with an injection-calibrated efficiency and a rate
per lens class against the per-lens forecasts. This file is the worker's own. Scripts: `ztf_lsn_targets.py` (the fixed
lens list, `lsn_targets.csv`, committed), `ztf_lsn_fetch.py` (IRSA IBE difference-image cutouts; R-127's access code),
`ztf_lsn_phot.py` (forced photometry at the lens centre and at four control positions), `ztf_lsn_detect.py` (the
statistic, the null from the controls, the candidates), `ztf_lsn_inject.py` (the injection set and the efficiencies),
`ztf_lsn_rate.py` (the table). Data under `data/ztf_lsn/<id>/` (gitignored): `epochs.csv`, `cut/*_diff.fits`, `phot.csv`.

## Trajectory check (run in the dispatching session, 2026-11-01, context-free Opus; not rerun here)

"VERDICT: partly novel: archival per-lens lensed-SN searches of the same known-lens catalogues exist (Sheu et al.
arXiv:2301.03578 on 5,807 Legacy Survey lenses in DECam epochs; Magee et al. arXiv:2303.15439 cross-matching ZTF transients
to known lenses, null; Townsend et al. arXiv:2405.18589 ZTF archive candidates); what remains open is ZTF/ATLAS forced
photometry at every lens position below the alert threshold, with injection-recovery efficiency and a per-lens rate upper
limit compared against the Goldstein/Wojtak/Arendse forecasts".

## Row 1: what Sheu et al. 2023 (arXiv:2301.03578) contains, read in full (the LaTeX source)

A targeted search of 5,807 lenses and candidates (Huang et al. 2020, 2021, Storfer et al. 2022 and eight smaller lists) in
every DECam g, r, z exposure of Legacy Surveys DR9 (DECaLS region, Dec < 32; BASS/MzLS excluded), with difference imaging
(Bramich 2008 and SFFT) and a detection grouping. It contains: (i) a **forecast** of what its own pipeline should find,
from Shu et al. 2018's rate formulation (SFR-based Ia and CC rates at a source redshift drawn as the lens redshift times a
truncated N(2, 0.5); per-image magnifications lognormal with ln mu ~ N(1.5, 0.35), mean 4.77; 2 or 4 images, 0.7/0.3;
delays N(36, 4) d; SALT3 Ia with M_B = -19.25 and Nugent/Vincenzi CC templates), run 1,000 times on the real DECam
exposure times: 2.90 ± 1.68 lensed Ia and 4.53 ± 2.12 lensed CC with three or more sub-detections, 6.58 ± 2.57 and
9.35 ± 3.08 with two or more; (ii) a **pipeline test** on 32 known DES SNe Ia (100 % recovered, photometry consistent with
DES), which is a bright-end check, not an efficiency curve: no injections into the lens cutouts, no efficiency per peak
magnitude or per lens class; (iii) **seven lensed-SN candidates** (one grade A, DESI-344.6252-48.8977, two grade B) plus
one unlensed SN, each with four light-curve postulations fitted (unlensed/lensed Ia/CC) and the magnification each lensed
postulation requires; no rate or rate limit is derived from the candidates against the forecast (the paper presents the
candidates and the forecast side by side, with the candidates needing follow-up). Of the eight transients, seven are at
Dec < -30; one, DESI-035.1374+00.4676 (grade C; best fit an unlensed Ia at z = 0.27), is in the ZTF footprint and is a
check target here (class X).

**So the open part is as the trajectory check said, with one narrowing:** an efficiency-calibrated per-lens rate limit
exists nowhere (Sheu's forecast is not calibrated on its own candidates' efficiency, and its candidates are unconfirmed),
so the result here is the ZTF sub-threshold search with its injection efficiency and the comparison with the per-lens
forecast (Shu et al. 2018, which gives N_Ia and N_cc per year for each of its 128 lenses) and with the blind-survey ZTF
forecasts (Goldstein et al. 2019 arXiv:1809.10147: 1.23 Ia and ~7 CC per year; Wojtak et al. 2019 arXiv:1903.07687: 2.1
Ia and 4.8 CC per year by the magnification method at r < 20.6; Sainz de Murieta et al. 2023 / Arendse et al. 2024
arXiv:2312.04621: ~0.2 Ia per year with the observing logs, as Townsend et al. quote them), which are for the whole
survey, not per known lens, and are reported beside the table as context.

## Reachability from the cloud (2026-11-01)

IRSA IBE answers (image search ~15 s per position; difference-image cutouts at ~16 files/s over 24 keep-alive
connections; 48 connections draw 503s). The ZTF forced-photometry service (ztfweb) needs an account (401), so the forced
photometry here is done on the difference-image cutouts, as in R-127. The ATLAS forced-photometry server answers its
front page but queues jobs only for registered accounts: **no ATLAS; ZTF only**. VizieR via IDIA serves SLACS (Auger et
al. 2009), Huang et al. 2020 and 2021; the Storfer et al. 2022 table came from the project site's Google Drive CSV
(`uc?export=download`); Shu et al. 2018's Table 1 from the arXiv LaTeX source. sncosmo installs and downloads its
templates (SALT3, Nugent, Hsiao). MasterLens is down; the BELLS and SL2S VizieR tables are not at IDIA (Shu's compilation
carries the BELLS positions, so BELLS is in).

## Fixed before any light curve was extracted (2026-11-01, 07:40 UTC; `lsn_targets.csv` committed with this file)

**Lens list and classes** (`ztf_lsn_targets.py`): class **S**, the 128 spectroscopic galaxy-galaxy lenses of Shu et al.
2018 (SLACS 63, S4TM 40, BELLS 25; all at Dec > -30), each with its forecast N_Ia, N_cc per year (110 of them; 18 have no
[O II] rate); class **L**, 311 grade-A Legacy Surveys candidates (Huang 2021 197, Huang 2020 60, Storfer 2022 54) at
Dec > -30, not within 3" of a class-S lens, duplicates merged; class **X**, Sheu's in-footprint candidate. Order: X, then
S (by N_cc + N_Ia descending, no-rate lenses last) and L (by lens r magnitude ascending) interleaved one for one. The
sample is the prefix of that order the download budget allows (about 4 h at 16 files/s; r band only, ~1,000-2,500 epochs
per lens), reported as such. g band is fetched only for candidates.

**Window:** every public ZTF r-band epoch IRSA serves for the position (2018-03 to the current public release), after
the cuts: INFOBITS != 0, SEEING > 4", MAGLIM < 19.0, cutouts with more than half their pixels undefined.

**Photometry** (`ztf_lsn_phot.py`): on each 30" difference cutout, a Gaussian PSF with FWHM = the header SEEING at the
catalogue lens position, flux and constant background by weighted least squares in a 15 × 15 pixel box, weights from the
clipped sky variance of the cutout outside 3" of the centre, pixels off by more than 5 sigma masked and the fit repeated;
flux to microJansky with MAGZP. The same at four **control positions**, 10" N, S, E and W of the lens, in the same
cutouts. Per (position, field, ccd, quadrant) the median flux is the baseline (a transient that fell in the reference
coadd leaves a constant negative offset, which this absorbs) and the formal errors are scaled by the robust standard
deviation of (flux - baseline)/error over the light curve, 3-sigma clipped, floored at 1 (the excess noise of galaxy-core
subtractions).

**Detection statistic** (`ztf_lsn_detect.py`): windows of 15 days starting every 5 days; in each window with at least two
epochs, S_w = sum(f/sigma^2) / sqrt(sum(1/sigma^2)) of the baseline-subtracted fluxes (the S/N of the window's weighted
mean); **T = max_w S_w**, one-sided. The null is the T distribution of the 4 × N control light curves. A lens is a
**candidate** when its T exceeds both 5.0 and the 99.5th percentile of the control T's (about 0.5 % × N false candidates
expected, every one vetted and reported). A **detection** is a candidate that survives vetting: at least two adjacent
windows positive at S_w > 3 (a lensed SN lasts longer than 15 observed days), the stacked difference image of the window
centred within 2.5" of the lens, a g-band light curve (fetched for the candidate) consistent with a transient of SN
colour, no persistent variability at the position (a lensed quasar or AGN shows in every season), no known transient
or asteroid at the time, no residual at the four controls in the same window (a bad subtraction shows everywhere).

**Injection set** (`ztf_lsn_inject.py`): per lens, 200 lensed SNe drawn as the forecasts draw them: type Ia with
probability N_Ia/(N_Ia + N_cc) for the lens (class S; the class-S mean 0.11 for L and X), else CC split IIP 55.8 %, Ic
17.0 %, IIb 12.4 %, Ib 9.0 %, IIL 3.3 %, IIn 2.4 % (Sheu's table, from Shu); templates SALT3 (x1 ~ N(0, 1), c ~
N(0, 0.1), M_B = -19.25), nugent-sn2p, nugent-sn1bc (Ib, Ic), nugent-sn2l, nugent-sn2n, with the peak M_B and scatter of
Sheu's table (IIP -16.75 ± 0.98, Ic -17.66 ± 1.18, IIb -16.99 ± 0.92 on the Ib/c template, Ib -17.45 ± 1.12, IIL -17.98
± 0.86, IIn -18.53 ± 1.36); source redshift z_s = the lens's (class S) or z_l times a N(2, 0.5) truncated below at 1
(L, X, as Sheu); 2 or 4 images (0.7/0.3) with per-image magnifications lognormal ln mu ~ N(1.5, 0.35), rescaled for class
S so their sum is Shu's total mu; relative delays uniform on 0-30 d (galaxy-scale lenses); image positions on a circle of
radius 1.2" about the lens centre at random position angles (the Einstein radii of SLACS lenses are 0.7-1.8"); explosion
time uniform over [first epoch - 100 d, last epoch + 100 d] (so seasonal gaps are in the efficiency); no host dust, no
Milky Way extinction (the sample is at |b| > 20). Each image's r-band light curve enters the measured light curve
through the measurement operator's response to a point source at that offset in that epoch's seeing (the forced
photometry is linear in the image, so this equals a pixel-level injection and re-measurement, minus the 5-sigma clip),
added to the real baseline-subtracted fluxes, and the statistic is recomputed. **Efficiency** per lens = detected
(T above the candidate threshold, with the two-window shape cut) / 200; reported per class, per type, and per peak
observed r magnitude of the summed images in 0.5 mag bins from 18.5 to 23.

**The rate** (`ztf_lsn_rate.py`): expected detections per class E = sum over lenses of (N_Ia + N_cc) × T_i × eps_i, where
T_i is the injection window in years and eps_i the lens's efficiency; class S uses Shu's per-lens rates (class median for
the 18 without); class L has no per-lens forecast, so the SLACS-like medians (N_Ia 0.004, N_cc 0.030 per year) are the
proxy and the BELLS-like medians the high bracket, stated as such. The claim per class is N lenses, the efficiency, the
candidates with their vetting, E, and the detections/E ratio with its Poisson 95 % interval (upper limit 3.00/E for zero).
Nothing about lens-population physics beyond what the ratio says.

## Rows (filled as the work lands)

**Row 2, the check target (class X), 08:10 UTC.** Sheu's transient in DESI-035.1374+00.4676 was seen in DECam exposures
of 2017-11-16 to 2017-12-24 (MJD 58073-58111, the paper's detection figure), four months before ZTF's first epoch, and
about 10" from the lens galaxy (their postulation is an unlensed Ia at z = 0.27 in a neighbour). So ZTF cannot test it;
the ZTF light curve at the lens (379 clean r epochs, 2018-07 to 2026-01) has T = 2.7, nothing. Kept in the table as its
own class, with no forecast weight.

**Row 3, validation on real supernovae and the PSF calibration (08:20-08:40 UTC; before any lens candidate above
threshold was examined, none being above it in the first 15 targets).** Eight ALeRCE SNe Ia (lc_classifier probability
> 0.8, short-lived, alert peak r 18.4-20.5; `val_targets.csv`) were fetched and measured exactly as the lenses. The alert
exposures are partnership ones that the public IBE does not serve (nearest public epoch 1-5 d away), so the comparison is
with the nearest public epoch within 2 d. The statistic finds the SNe at their peaks (ZTF18adjlalc T = 15.8 at +1 d,
ZTF24abdooqu S = 16.3 at peak, ZTF20aaxyglx 6.6, ZTF21absjikb 4.5) but the fluxes came out 0.1-0.9 mag fainter than the
alert PSF magnitudes. Cause, found with the pipeline's own difference-image PSF stamps (`diffimgpsf.fits`,
`ztf_lsn_psfcal.py`): the difference image's PSF is the science and reference PSFs convolved to a common one, broader than
the science seeing (moment FWHM 3.5" at SEEING 2.2" for ZTF24abdooqu), so the Gaussian of FWHM = SEEING with a free
background recovers kappa0 = 0.77 (0.50-0.83) of a central point source and kappa12 = 0.50 (0.21-0.82) of one 1.2" off
centre, instead of the Gaussian's own 1.00 and 0.63. **Amendment to the fixed protocol (a calibration, not a selection
choice):** for every target the PSF stamp of every 10th r epoch is fetched, kappa0 and kappa12 measured per stamp, and
carried to every epoch by a linear fit on SEEING per (field, ccd, quadrant); the injections use kappa12(epoch) as the
response of an image 1.2" off centre (`ztf_lsn_inject.py`, `psf_response`), which lowers every efficiency relative to the
Gaussian response (the injections run before this amendment were discarded); candidate magnitudes are reported corrected by
kappa0. The data side of the statistic (measured fluxes and errors) is unchanged: the estimator is sub-optimal for the
broad PSF, and that inefficiency is in the efficiency, where it belongs. With kappa0 applied, ZTF24abdooqu's epochs near
peak read r 20.2-20.5 against the alerts' 20.1-20.3.

The validation table after the calibration (`data/ztf_lsn/validation.csv`; kappa0 per target 0.68-0.81, kappa12
0.37-0.55, PSF moment FWHM 3.4-3.7" at SEEING 1.8-2.4"): the three SNe with many alert detections agree with the alerts'
PSF magnitudes to within 0.3 mag at the nearest public epoch (ZTF24abdooqu +0.19 ± 0.08 over 5 epochs, ZTF18adjlalc
−0.31 ± 0.08 over 3, ZTF25acfxoib +0.23 over 2, the sign being mine minus alert), and the statistic finds them at their
peaks (T = 17.6, 15.8, 11.9; ZTF24abdooqu's best window is 87 d after the alert peak, in its denser late sampling, with
S = 16.3 at the peak itself). ZTF20aazpora, ZTF21absjikb and ZTF20aaxyglx (3, 7 and 8 r-band alerts; the last flat at
r 20.5-20.7 for 117 d, not a SN Ia light curve; the first with its best window 1,751 d before the alert) come out 0.77 mag
fainter: they are variable nuclei, whose flux relative to my seven-year median is not the alerts' flux relative to the
reference, so they do not calibrate the photometry. ZTF21aafxvah has 231 clean epochs of 833 (INFOBITS 2097152 on most)
and ZTF19aaseomv a negative S at the alert peak (−5.8: its reference holds the SN); both are noted, neither used. So the
forced photometry is right to ~0.2-0.3 mag on real SNe after the PSF calibration, and the statistic recovers r ~ 20 SNe
at T = 12-18; the efficiency does not depend on this check (it uses the stamps' kappa12 directly).

**Row 4, candidates in the first 122 targets (09:00-09:20 UTC; the full-sample row supersedes the counts).** Lens T:
median 3.08, 99.5th percentile 12.2, 10 of 114 above 5.0; control T (456 light curves): median 2.67, 99.5th percentile
4.07, max 5.38 (one above 5, a two-epoch window). So the lens positions have a heavier tail than blank sky even after the
error scaling, and the threshold is the fixed 5.0. The ten, vetted (`ztf_lsn_vet.py`, g band fetched for the top four;
ALeRCE cone search 6" for all four: no ZTF alert object at any of them):

| lens | T | T2 | window (JD) | r light curve (15-d means, microJansky, measured) | g | verdict |
|---|---|---|---|---|---|---|
| SLACS J1250+0523 (z_l 0.23, z_s 0.80, mu 10) | 13.7 | 13.7 | 2461110 (2026-03) | −1.6, 7.0, **13.6**, 13.5, 5.0, 0.5, 1.6 (−30 d to +90 d); controls 0 ± 0.7; centroid 0.4" from the lens; earlier seasons' max S ≤ 3.7 | −0.3, 8.6, **10.6**, 6.4, 0.4: peaks earlier, fades faster | **a real transient the alert stream missed**: kappa-corrected peak r ≈ 20.8, g−r ≈ +0.3 at the r peak reddening to ≈ +0.8 at +20 d, ~45 d above half maximum. At the lens redshift that is M_r ≈ −19.5 with a SN Ia's colour evolution and duration: an unlensed SN (Ia-like) in the lens galaxy, or a nuclear flare. Not a lensed SN of the z = 0.80 source: a Ia there has g in the rest-frame far UV (observed g−r > 2) and peaks near r 24.5, a magnification of ~30 against the model's 10; a CC SN would need mu > 50. **Candidate, vetted as not lensed.** |
| SLACS J1621+3931 (0.24, 0.60, 8) | 10.9 | 10.9 | 2458990 (2020-05) | season-long: 0.5, −1.0, **4.0**, 3.4, 3.1, 3.3, 2.7; season means −1.7, −0.8, +2.4, +1.9, +0.6, +0.8 (2018-2024) | +0.7 | slow, red, lasting a season and more: nuclear variability of the lens galaxy, not a SN. Vetted out (persistent variability). |
| S4TM J0753+3416 (0.14, 0.96, 27) | 9.2 | 2.3 | 2458845 | 203 epochs in one 15-d window (a high-cadence run) at 1.9 ± 0.2, the next window −4.3 ± 0.8 | | the shape cut fails; a high-cadence night's correlated residual at a bright core. Vetted out. |
| S4TM J0757+1956 (0.12, 0.83, 9) | 7.0 | 5.4 | 2461080 (2026-02) | 1.9, 11.6 ± 11, **9.2**, 3.4, 1.9, 2.2, −1.2; controls ±1; centroid 0.3" | 4.7, 1.7, 3.3 | a short red bump (kappa-corrected r ≈ 21.1, g−r ≈ +0.7), 15-30 d. At z_l = 0.12 M_r ≈ −17.6: a faint SN or a flare in the lens galaxy; not a lensed SN of the z = 0.83 source at mu ~ 9 (needs mu > 100). **Candidate, vetted as not lensed** (and marginal: T2 = 5.4). |
| DESI-266.4075+39.8584 (L) | 6.0 | 3.3 | 2460440 | 4.1, 1.6, **6.0**, 1.7, 2.9, 1.0, 1.6 | | single window; fails the shape cut. |
| S4TM J1031+3026 (0.17, 0.75) | 5.7 | 3.0 | 2458870 | −1.0, 0.7, **3.2**, −0.4, 1.6 | | single window; fails. |
| DESI-189.9885+12.6693 (L) | 5.3 | 3.8 | 2459340 | 0.8, −1.5, **2.2**, 0.8, 1.5, 1.4 | | single window at 2 microJansky; fails. |
| DESI-341.0206+27.9877 (L) | 5.3 | 2.0 | 2458350 | −0.1, −0.6, **4.1**, −3.0, 0.7 (140 epochs in all) | | single window; fails. |
| SLACS J1153+4612 (0.08, 0.88) | 5.1 | 2.2 | 2458565 | 2.8, 0.9, **7.1**, 1.9, 1.6, 1.9, 2.7; control W 2.9 ± 1.0 in the window | | single window with a control also high; fails. |
| SLACS J1432+6317 (0.12, 0.66) | 5.1 | 4.0 | 2461010 (2025-11) | **12.1** ± 2.4, 6.3, 11.4 ± 6, −1.7, 9.9 ± 3.5 (5, 3, 2, 3, 5 epochs); season mean +3.2 ± 0.5 after four quiet seasons | | sparse and noisy; T2 = 4.0 < 5: fails the cut as fixed, but the season-level rise is noted as a loose end. |

The control null is not blank sky: at 190 targets (744 controls) its 99.5th percentile is 5.28 and its maximum 9.9, the
latter 10" south of DESI-353.5287-07.7712 at JD 2459140, where the alert stream holds a single-detection object,
ZTF20acjkvsr (MJD 59137): a real transient in a neighbour, found by six sub-threshold epochs here and once by the alerts.
So the control-derived threshold is conservative (it carries real transients at the 0.5 % level), and the statistic does
what it is for. The band filter was tightened to r epochs only after the g-band fetch for the four top candidates was
found to feed g epochs into their statistic and injections (09:55 UTC; those four injections rerun).

The lensed-SN vetting applied to the survivors (the brief's "colour" and "lens model" items, made explicit here): a
candidate counts as a lensed SN only if its kappa-corrected peak and g−r are consistent with a SN at the source redshift
magnified within the lens model's range (Shu's total mu times 0.3-3); both survivors are transients in the lens galaxies
themselves, which forced photometry at the lens centre cannot separate from a lensed image by position alone, and which
colour and brightness separate here.

**Photometry check (first 15 targets):** at the four control positions the formal errors are right (robust scatter of
flux/error 0.90-1.05); at the lens centres the scale factor is 1.0-1.4 (galaxy-core subtraction residuals), carried into
the errors as fixed. 43 % of the epochs carry INFOBITS != 0 (mostly bit 22, 4194304); at the controls those epochs look as
clean as the rest (scatter 0.90, median MAGLIM 20.9 against 20.5), so the pre-registered cut costs depth without buying
purity; it stands for the claim, and `ztf_lsn_detect.py --keep-infobits` is reported as a robustness row. The 15-day
window error at the lens is 2.1-2.6 microJansky (5 sigma: r ~ 21.2-21.5) for an equatorial target with ~380 clean epochs.


## R-135 (the resume, 2026-11-01 19:15-23:40 UTC)

**Row 5, the full sample and the table (23:30 UTC).** The container's data did not survive R-130's idle, so the cutouts
were fetched again (same scripts, IRSA IBE at 15-22 files/s over three or four fetchers); the fetch stopped at 22:15 with
the **contiguous prefix of ranks 0-184** of the fixed order: the class X check, 92 class S, 92 class L. Nine of those have no
usable photometry (the difference images are masked at the position: constant fill values, or no cutout), and four more
have fewer than 20 clean epochs or 8 PSF stamps, so **168 lenses carry an efficiency** (85 S, 82 L, 1 X). Window: JD
2458200-2461300 (2018-03 to 2026-08); injection window per lens 8.4-8.6 years.

Candidate threshold from the 704 control light curves: T > **5.31** (their 99.5th percentile; maximum 9.89, the neighbour
transient ZTF20acjkvsr of row 4). Lens T: median 3.19, 99.5th percentile 11.3. The ten candidates, all vetted:

| lens | class | T | T2 (shape cut) | window | verdict |
|---|---|---|---|---|---|
| SLACS J1250+0523 | S | 13.75 | 13.75 | 2026-03 | row 4: a transient in the lens galaxy (r ≈ 20.8, SN-Ia-like colour and duration at z_l = 0.23); not a lensed SN of the z = 0.80 source (needs mu ~ 30 for a Ia, > 50 for a CC, against 10) |
| SLACS J1621+3931 | S | 10.88 | 10.88 | 2020-05 | row 4: season-long nuclear variability; out (persistent variability) |
| S4TM J0753+3416 | S | 9.18 | 2.31 | 2019-12 | fails the shape cut (one high-cadence night's residual) |
| S4TM J0757+1956 | S | 6.98 | 5.42 | 2026-02 | row 4: a short red bump in the lens galaxy (r ≈ 21.1, M_r ≈ −17.6 at z_l); not lensed (needs mu > 100) |
| SLACS J0037-0942 | S | 6.24 | 0.92 | 2024-09 | fails the shape cut (one 10-epoch window, the neighbours one epoch each) |
| DESI-266.4075+39.8584 | L | 5.99 | 3.32 | 2024-05 | fails the shape cut |
| DESI-234.8707+16.8379 | L | 5.89 | 5.64 | 2020-05 | passes the shape cut; new here. Season means +1.6, +1.6, **+2.8, +2.6**, −2.1, −1.2, −1.4, −1.2 microJansky (2018-2025): a two-season plateau then a step down, controls quiet, centroid 0.5": persistent variability of the lens galaxy (an AGN-like step), not a SN; out |
| S4TM J1031+3026 | S | 5.75 | 2.96 | 2020-01 | fails the shape cut |
| S4TM J1101+1523 | S | 5.58 | 2.84 | 2023-10 | fails the shape cut; season means drift −1.7 to +2.2 to −1.1 (slow variability) |
| SLACS J0959+0410 | S | 5.55 | 2.51 | 2023-04 | fails the shape cut (4 epochs) |

**Detections: none.** (R-130's 122-target sub-threshold names DESI-189.9885+12.6693, DESI-341.0206+27.9877, SLACS
J1153+4612 and SLACS J1432+6317 sit below the final 5.31 threshold.)

**Two amendments to the injections, both found on the full sample before the table was read, both lowering E:**
(1) a recovery counts only in windows from 60 d before the injected SN's model t0 to 200 observed days after its last
image (the fixed text says "detected (T above the candidate threshold, with the two-window shape cut)", which, applied to a
lens whose own light curve passes the cut, credits every injection; the untimed efficiencies were 0.075 (L) and 0.059 (S),
kept in `inject_summary_untimed.csv`); (2) where the lens's own light curve already passes the cut inside that span, the
injection counts as missed (no efficiency is credited during a real event; this touches the four lenses whose T2 > 5.31).

**The table** (`ztf_lsn_rate.py`; E = sum over lenses of (N_Ia eps_Ia + N_cc eps_cc) × years):

| class | N | epochs/lens | 15-d depth (5σ, r) | eff all | eff Ia | eff CC | eff, peak r < 21 | candidates | detections | forecast to occur | E | detections/E, 95 % |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S (Shu 2018 per-lens rates) | 85 | 791 | 21.58 | 0.024 | 0.110 | 0.013 | 0.43 | 8 | 0 | 90.1 | 0.589 | < 5.1 |
| L (SLACS-like proxy) | 82 | 639 | 21.65 | 0.064 | 0.219 | 0.045 | 0.56 | 2 | 0 | 23.5 | 1.573 | < 1.9 |
| L (BELLS-like bracket) | 82 | | | | | | | | 0 | | 7.48 | < 0.40 |
| X (Sheu's candidate) | 1 | 377 | 21.12 | 0.025 | | | 0.21 | 0 | 0 | 0.3 | 0.007 | (no weight) |

Class S by survey: SLACS 44 lenses E = 0.523, S4TM 24 E = 0.246, BELLS 17 E = 0.021 (BELLS holds 46 of the 90 forecast
events, at z_s 1-1.5, where the injected peaks are r 26-29). Ia carries 0.34 of S's 0.59. Efficiency per summed peak r
magnitude (class S / L): 18.5-19 0.55/0.51, 19-19.5 0.47/0.69, 19.5-20 0.50/0.72, 20-20.5 0.39/0.55, 20.5-21 0.24/0.36,
21-21.5 0.12/0.15, 21.5-22 0.03/0.06, 22-22.5 0.01/0.01. Forecast-free: lensed SNe with summed peak r < 21 occur at
< 0.020 per lens-year (S, 147.6 effective lens-years) and < 0.0086 (L, 347 effective lens-years), 95 %.

So: the per-lens forecast for the 85 Shu lenses predicts 90 lensed SNe in the window and 0.59 detections at the measured
efficiency; seeing none limits the rate to < 5.1 times the forecast. The class L limit is < 1.9 times the SLACS-like proxy
and excludes the BELLS-like rates for the grade-A candidates (< 0.40 of them), with the proxy's assumption carrying that.

**Row 6, the finding check (2026-11-01 ~23:40 UTC; fresh context-free Opus, the section alone, `blind-review.md`'s prompt
verbatim):** "VERDICT: novel but not sound: the quoted efficiencies come from injections that skip the colour, magnification and variability vetting that every real candidate had to pass, and the class-L limits that carry the headline come from an assumed source population rather than from the data."
A miss on a method point; parked with the path to pass recorded under the section in `finding.md`. No revision pass.


## R-141 (cosmo-deep-10b, the one revision, 2026-11-02 from 00:00 UTC)

The finding check's two points (R-135 row 6) are both addressed here without new data: (1) the candidate vetting as code,
`ztf_lsn_vetting.py`, applied unchanged to every real light curve and to every injected one (`ztf_lsn_inject_vetted.py`;
R-135's injection set drawn identically, same seeds, same kappa12 response); (2) the headline becomes the forecast-free
limit, with the forecast comparisons as secondary rows labelled with the source population each assumes.

**Row 7, the vetting as code, fixed before any injection was vetted (00:50 UTC).** A light curve is a detection when all
of these pass, in order (the cuts R-135 applied by hand, rows 4-5, each made a number):

| stage | cut | real-data input | injection input |
|---|---|---|---|
| 0 detection | T2 > 5.31 (R-135's threshold from its 704 controls): a window above threshold with a neighbour at S > 3; for injections, R-135's two amendments (time-matched; no credit during a real event) | windows of the lens light curve | the same, with the injected flux |
| 1 controls quiet | no control position (10" N, S, E, W) has S > 3 in the detection window | the real controls | the real controls (the injection does not touch them) |
| 2 centroid | window centroid within 2.5" of the lens | the stacked window image (`ztf_lsn_vet.stack`) | flux-weighted mean of the image positions (1.2", random position angles) plus centroid noise PSF sigma / S per axis (the injections are in flux, not pixels) |
| 3 no persistent variability | no season (365.25 d from JD 2458000, >= 5 epochs, outside [t_w − 100, t_w + 250]) deviating from the median season mean by more than VAR_AMP at more than VAR_SIG sigma | r light curve | r light curve with the injection |
| 4 colour and magnification | some SN hypothesis (7 types: SALT3 Ia x1 = c = 0, M_B −19.25; IIP, Ic, IIb, Ib, IIL, IIn at their mean M_B, as the injections draw them) at the source redshift (class S: measured; L, X: z_l × 1.0-3.0 in steps of 0.25) fitted to the r epochs in [t_w − 60, t_w + 120] (amplitude and explosion time, kappa0-corrected as the candidates' magnitudes were) requires a total magnification within 0.3-3 × the lens model's (S: Shu's mu; L, X: the injection draw's median total, 10.6), and the g epochs in the span are consistent with it (chi2 p > 0.001, the amplitude's error included; no g epochs = untested) | real g photometry where fetched | synthetic g at the lens's real clean g epochs (IBE rows): the SN's ZTF-g flux × the kappa12-on-seeing fit at the g seeing, plus Gaussian noise of the r band's error at that limiting magnitude (median e / (F_5σ/5) of the lens's r epochs, times its error scale) |
| 5 no alert object at the time | (ALeRCE, by hand) | done in row 4 | none by construction |

VAR_AMP and VAR_SIG are the one place where R-135's call ("persistent variability") had no number; they are set on the four
real candidates that passed stages 0-2 in R-135 (the two variables must fail, the two lens-galaxy transients must pass) and
frozen in row 8 before the injections are vetted, with their false-rejection rate on the control positions reported.

**Row 8, the vetting run on the real candidates, and its frozen numbers (00:55 UTC, before any injection was vetted).**
The three R-135 candidates that passed stages 0-2 and lie in ranks 55-89 were refetched first (r band in full; g band only
in JD 2460990-2461260 for the two transients; DESI-234.8707+16.8379, rank 126, was not refetched: IRSA throttled, see the
data note below). Through `ztf_lsn_vetting.vet`:

| lens | T (R-135) | stage reached | numbers |
|---|---|---|---|
| SLACS J1621+3931 | 10.78 (10.88) | fails 3 (variability) | a season 2.3 microJansky off the median of the other seasons at 12 sigma |
| SLACS J1250+0523 | 13.61 (13.75) | fails 4 (colour) under the peak form; **passes** under the chi2 form | Ia at z_s 0.80 needs mu = 26.6 ± 1.7 (inside 3-30), but g over the r peak is 7.1 ± 1.3 microJansky against 1.6 predicted (+4.3 sigma); IIn needs 24.7, g 7.1 against 11.6 (−3.4 sigma); every CC type but IIn needs mu > 140 |
| S4TM J0757+1956 | 7.15 (6.98), T2 5.53 | **passes every stage** (both forms) | Ia at z_s 0.83 needs mu = 17.0 ± 2.4 (inside 2.7-27); g over the r peak 1.8 ± 1.0 against about 1 predicted |
| DESI-234.8707+16.8379 | (5.89) | fails 3 by R-135's recorded season means (+2.8, +2.6 against a median near −1.2, errors ~0.3) | not rerun |

**Frozen: VAR_AMP = 1.0 microJansky, VAR_SIG = 5** (J1621 at 2.3/12 sigma and DESI-234 at ≈ 4/>10 sigma fail; the two
transients have no significant season outside their event span, 0.0; on the first 11 lenses' 44 control light curves the
largest season deviation is 2.5 sigma, so the cut removes no blank sky).

**Amendment to row 7's colour test (made here, before any injection was vetted, and applied identically to the
injections):** row 7 froze "g epochs in the span consistent with the hypothesis, chi2 p > 0.001" over [t_w − 60, t_w + 120].
On J1250 that test passes a lensed Ia (p = 0.76) although g at the r peak is 4.3 sigma too bright for it: the 58 g epochs
outside the peak dilute the chi2. R-135's call was made on the colour at the r peak ("g−r at the r peak"), so the coded test
is now that: the weighted mean g flux over the r peak windows [t_w, t_w + 30) within 3 sigma of A × the hypothesis's g
(A's error included). The first-frozen chi2 form is kept as a robustness row (`--colour chi2`); under it J1250 also
survives. Also fixed here: the synthetic g errors carry a factor 1.15 (g's e/F_5σ over r's at the lens centre on the two
lenses with real g: 1.18, 1.12).

**What the coded vetting changes in R-135's candidate list:** S4TM J0757+1956 was "vetted as not lensed" in R-135 on "needs
mu > 100"; that holds for every core-collapse type (100-260) but not for a SN Ia at the source redshift (peak r 24.25 at
z = 0.83 for M_B = −19.25, so mu ≈ 17 for its kappa0-corrected r ≈ 21.1, inside Shu's mu = 9 × 0.3-3), whose g at the
peak and duration (29 d above half maximum at z = 0.83, against the bump's 15-30 d) are also consistent. It is equally a
faint transient in the lens galaxy (M_r ≈ −17.6 at z_l = 0.12); ZTF photometry does not separate the two. **It counts as
a detection in every limit below** (the conservative reading; the k = 0 numbers are given beside).

**Data note (R-141).** The container's data did not survive; IRSA was rate-limited throughout (503s at any concurrency above
~4; 4.4-4.7 files/s in total against R-135's 15-22), so instead of the 168 lenses the refetch ran 00:08-03:36 UTC (about
3.5 h) over one rank-ordered queue with the fixed epoch cut applied before download (`ztf_lsn_refetch.py`), and the sample
is the contiguous prefix of ranks 0-65: 66 targets, 5 masked, 7 with < 20 clean epochs or < 8 stamps, **54 with an
efficiency** (27 S, 26 L, 1 X); 53,000 cutouts and stamps. The three candidates of rows 4-5 at ranks 83-89 were fetched out
of order for the vetting check only (row 8) and are not in the sample.

**Row 9, the vetted table (03:40 UTC; `ztf_lsn_rate_vetted.py --stop 66 --detections S4TMJ0757+1956`).** Candidates in
the sample above 5.31: S4TM J0757+1956 (passes every stage: counted as a detection) and DESI-266.4075+39.8584 (fails the
shape cut). Control 99.5th percentile in this sample 4.08 (216 light curves).

| class | N | eff before / after vetting (all injected) | eff at summed peak r < 21, before / after | X<21 lens-yr, before / after | detections | R95 per lens-yr, peak r < 21 (with eff uncertainty) |
|---|---|---|---|---|---|---|
| L | 26 | 0.083 / 0.046 | 0.62 / 0.34 | 113.5 / 79.7 | 0 | < 0.045 (< 0.038 without; < 0.030 for r < 20) |
| S | 27 | 0.0028 / 0.0019 | 0.38 / 0.19 (21 injections) | 12.0 / 5.6 | 1 | < 1.7 (< 0.95 if k = 0) |
| X | 1 | 0.030 / 0.030 | | 2.3 | 0 | (no weight) |

By stage, class L at r < 21: detected 0.624, controls 0.619, centroid 0.619, variability 0.351, colour/magnification 0.344.
Per summed peak r (L, before → after): 18.5 0.94→0.38, 19.0 0.90→0.50, 19.5 0.82→0.50, 20.0 0.62→0.40, 20.5 0.31→0.19,
21.0 0.12→0.08, 21.5 0.03→0.01. Secondary: S Shu per-lens 67.5 occur, E 0.095 → 0.065, < 73 × forecast; L SLACS-like
E 0.65 → 0.37, < 8.2; L BELLS-like E 3.09 → 1.73, < 1.7 (no longer excluded). Robustness (L, X<21, R95 without unc): chi2
colour 76.1, < 0.039; wide magnification 95.0, < 0.032; all images at 2.0" (kappa20 from the stamps) 22.7, < 0.13.

**Row 10, the finding check (03:45 UTC; fresh context-free Opus, the revised section alone, `blind-review.md`'s prompt
verbatim, hints as R-135's):** "VERDICT: novel but not sound: the headline limit rests on putting every lensed image 1.2" from the centre, an offset the efficiency is very sensitive to (all at 2" makes the limit three times weaker) and whose quoted 18 % uncertainty disagrees with the note's own 2σ figure, so the limit of fewer than 0.045 per lens-year and its error budget are not supported."
A miss on a physics point (the image offset). Novelty held again (closest: LCOLSS arXiv:2111.01680; Sheu et al.
arXiv:2301.03578; Magee et al. arXiv:2303.15439; Townsend et al. arXiv:2405.18589). This was the last pass: **parked.** One
erratum, not a rewording: the offset systematic used in the computation is 0.36 (half the relative exposure change, 0.715),
not the "18 %" the section stated; the limit printed (0.045) is the one computed with 0.36.
