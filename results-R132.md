# R-132 (cosmo-deep-12): the one-halo tSZ signal of group catalogues at 10^12.5-10^13.5 M_sun on the ACT DR6 y map

Trajectory check (2026-11-01, context-free Opus, in the dispatching session, verbatim): "VERDICT: partly novel: Planck-resolution
stacks on spectroscopic and DESI DR9 groups already measured the group Y-M down to about 10^12.5 M_sun (Lim et al.
arXiv:1710.06856; arXiv:2201.12591), and ACT DR6 stacks exist only for galaxies or LRGs (arXiv:2607.08721, 2609.08939, 2608.00164)
with no group-catalogue Y-M; the ACT-resolution one-halo Y-M amplitude and slope on Tully 2015/2MRS/DESI DR9 groups below
10^13.5 M_sun remains open". Not rerun.

## Pre-registration (written 2026-11-01 10:00 UTC, before any stack; `groups_ysz.py`)

**Reachability (this container).** NERSC ACT DR6 NILC y maps (`portal.nersc.gov/project/act/dr6_nilc/ymaps_20230220`: the
ACT+Planck ILC y map, the CIB-deprojected variants, the four masks, the beam) answer; the Planck Legacy Archive answers for the
MILCA y map and the R2.01 masks by `product-action?MAP.MAP_ID=`; the IDIA VizieR mirror answers (CDS, CfA, NAO, China-VO, IUCAA
are cut); HEASARC's NVSS and SUMSS tdat files answer; `data.desi.lbl.gov/public` answers. **The SJTU group-catalogue server
(`gax.sjtu.edu.cn`) is cut** from the cloud (connection reset over http and https), so the Yang et al. SDSS DR7 groups are not
reachable; the DESI DR9 Yang et al. 2021 groups are, as the DESI DR1 value-added catalogue `gfinder` v1.0
(`DESIDR9.y1.v1_group.fits`, 5.6 GB: IGRP, RICH, GRP_RA, GRP_DEC, GRP_Z, GRP_LOGM = log10(M_h h/M_sun), GRP_LOGL; MAG_Z < 21,
spec-z where DESI DR1 has one, photo-z otherwise).

**Catalogues and floors (fixed).**
1. *DESI DR9 Yang groups* (primary). z in [0.1, 0.6) (photo-z for most; the halo masses lean on the group luminosity, whose
   completeness at MAG_Z < 21 falls off past z ~ 0.6). Halo mass M_h (180 x mean density, Yang et al. 2021, Planck18 h = 0.674)
   converted to M500c with an NFW profile and the Dutton & Maccio 2014 concentration; bins in log10 M_h [M_sun]: 12.5-12.75 and
   12.75-13.0 (**below the referee's floor; reported, flagged, not in the headline fit**), 13.0-13.25, 13.25-13.5, 13.5-13.75,
   13.75-14.0, 14.0-14.5 (the cluster anchor). Two redshift halves [0.1, 0.3) and [0.3, 0.6) and the whole range. All
   richnesses in the primary sample; RICH >= 3 as the check.
2. *Tully 2015 2MRS groups* (VizieR J/AJ/149/171 table3; 25,474 nests): 1,500 < <Vcmba> < 10,000 km/s, M_lum (luminosity mass,
   ~M200c; converted to M500c) bins 12.5-13.0, 13.0-13.5, 13.5-14.0, >= 14.0; Nmb >= 2 (a single-galaxy "group" has no group
   luminosity). Low-z, resolved: theta500 of 10-40'.
3. *Saulder et al. 2016 SDSS DR12 groups* (VizieR J/A+A/596/A14 grlist_s; 229,893 groups; mock-calibrated luminosity masses
   logMtot, ~M200c) as the stand-in for the unreachable Yang SDSS catalogue: z in [0.02, 0.11), the same four mass bins, Ntot >= 2.
4. Optional, if time allows after 1-3 are done: Zou et al. 2021 DESI DR9 clusters (VizieR J/ApJS/253/56; M500 from L_1Mpc,
   12.5 <= log M500 <= 15.2) in the same bins, as a second mass calibration on the same imaging.
Mass scatter propagated into the template: lognormal with sigma = 0.25 dex for every catalogue (0.15 and 0.35 as checks).

**Footprint.** ACT: the release's wide GAL070 mask >= 0.99, not inpainted, not a subtracted source (pixel level, R-131's
`clean_bool`); a group enters if >= 99 % of its 15' disc is clean. Planck: the R2.01 GAL60 mask and the point-source mask.

**Aperture and template (fixed).** Annulus means of the y map around each group's position, edges 0, 1, 2, 3, 4, 6, 8, 10, 15,
20, 30, 45, 60' (ACT; the DESI fit uses the annuli to 30', the low-z catalogues all 12) and 0, 2.5, 5, 7.5, 10, 15, 20, 30, 40,
60' (Planck), computed as FFT convolutions of the masked map with the annulus kernels in declination strips of the CAR grid
and sampled at the nearest pixel (ACT 0.5', Planck reprojected to 2'); rotation nulls and randoms are the same lookups shifted
in RA. Template per group: the Arnaud et al. 2010 UPP at the group's R500 (from its M500c and z), truncated at 5 R500,
normalised to the Planck 2013 XX reference Y500 at the group's mass (E^-2/3 D_A^2 Y500 / 1e-4 Mpc^2 = 10^-0.19
[(0.8 M500)/6e14]^1.79; the same reference as R-131, so A reads directly as "observed over the Planck cluster relation"),
convolved with the release beam (ACT: `ilc_beam.txt`; Planck: 10' Gaussian and the nside-2048 pixel window), averaged over the
mass posterior, then averaged over the bin's groups. Fixed-aperture check: the 0-3' disc minus the 3-6' annulus (ACT).

**Model and statistic (fixed).** Per bin and map, y_k = A T1h_k + B T2h_k + c over the annuli k, fitted by GLS with the
covariance of the bin mean from the per-group scatter (N groups; the groups' profiles are sampled independently of each
other), validated against 23 RA-rotation nulls (every position shifted by 15 n degrees; where the nulls' scatter of A exceeds
the GLS error the error is inflated to it). T1h: the bin-averaged template above. T2h: the neighbour model, the y map painted
from every catalogue group with log M_h >= 12.5 at its reference template (Y-M at the reference relation), stacked at the bin's
positions with the same annuli, minus each group's own template; B is free (B = 1 is "the neighbours sit on the reference
relation"). c: the local background. Statistic: Y500 = A <Y_ref E^-2/3> over the bin, in Mpc^2 (D_A^2 Y500 E^-2/3); the claim
is Y500 per bin with its interval, and per catalogue and map the power law log10 Y500 = log10 A_Y + B_Y log10(M500c / 10^13.5)
fitted over the headline bins (log M_h >= 13.0 for DESI; all four for Tully and Saulder) with the bins' covariance, plus the
ratio of the ACT and Planck amplitudes at 10^13.5 with its interval. Self-similar reference: B_Y = 5/3.

**Radio sources (fixed).** NVSS (dec > -40) and SUMSS (dec < -40; S_843 scaled to 1.4 GHz with alpha = -0.7): a group is
"radio" if a source with S_1.4 > 5 mJy lies within 3' of its centre. The headline is the radio-masked sample; the unmasked
stack is reported beside it in every table, with the fraction masked per bin. The release's inpainted and subtracted-source
pixel masks apply always.

**Foreground check.** The CIB-deprojected map (`deproj_cib_1.7_10.7`) stacked identically; the shift in A per bin reported.

**Nulls.** Mask-matched random positions (A consistent with 0); the 23 rotations; the two-halo check is B against 1 and the
stack on the neighbours alone.

**Planck comparison.** The same groups, templates and model on the Planck MILCA map at 10'; and the literature at Planck
resolution: Lim et al. 2018 Table 1 (Y~500 per log M500 bin, 10^-6 arcmin^2 at d_A = 500 Mpc; 2.85 +- 1.17/sqrt(19020) at
13.5) and Chen, Zhang & Yang 2022 (arXiv:2201.12591: E^-3/2 D_A^2 Y / (A 1e-4 Mpc^2) = (M_L / 10^14.5 M_sun/h)^alpha with
A = 1.13 +- 0.13, 1.2 +- 0.09, 1.02 +- 0.08, 0.65 +- 0.08 and alpha = 1.81 +- 0.06, 2.03 +- 0.05, 2.1 +- 0.05, 2.12 +- 0.08 in
z 0-0.2, 0.2-0.4, 0.4-0.6, 0.6-1; richness >= 5; Y to R_vir with a KS profile).

**What the prior ACT work did that this must do better** (read before computing): Battaglia et al. 2026 (arXiv:2607.08721)
split Planck LBGs by an NVSS match within 1' with no flux cut, compensated aperture 3' and no two-halo term, stellar-mass
bins only above 10^11.1 ("noise dominated" below); Moore et al. 2026 (arXiv:2609.08939) masked targets with an NVSS source
inside the 2.1-2.7' aperture, found the two lowest-mass samples negative with the radio hosts in and the one- and two-halo
terms degenerate inside 4' ("the individual one-halo and two-halo profiles varied wildly"), and assigned halo masses from
simulations; Siegel et al. 2026 (arXiv:2608.00164) fitted dust and radio per bin from the single-frequency maps at 1.6' and
forward-modelled the two-halo term with FLAMINGO, declining to convert to Y500. Here: halo masses from the catalogues, the
UPP at the catalogue R500 (so Y500 is the statistic), the two-halo term from the catalogue's own positions with a free
amplitude and annuli to 30-60', and the radio cut with a stated flux limit, with and without.

**Planck-era numbers for the comparison (read from the papers' tables, 2026-11-01).** Lim et al. 2018 Table 1 (Y~500 = Y500
E^-2/3 (d_A/500 Mpc)^2 in 10^-6 arcmin^2; scatter among groups, N): log M500 12.9: 0.0675 (0.330, 37,344); 13.1: 0.144
(0.450, 32,063); 13.3: 0.735 (0.600, 25,744); 13.5: 2.85 (1.17, 19,020); 13.7: 6.60 (1.77, 12,500); 13.9: 30.3 (3.00,
6,203); 14.1: 59.7 (7.50, 2,163). 1 x 10^-6 arcmin^2 at 500 Mpc is 2.116 x 10^-8 Mpc^2, so Lim's 13.5 bin is 6.0 x 10^-8
Mpc^2 against the Planck 2013 reference (1-b = 0.8) of 22.3 x 10^-8: 0.27 of the reference. Chen, Zhang & Yang 2022 Table 1
(Yang DESI groups on MILCA, KS-profile amplitude A1 relative to their model with S/N): [0, 0.2) x [13, 13.5): N 60,269, A1
0.364 at S/N 2.2; [0.2, 0.4) x [13, 13.4): N 252,970, A1 0.382 at 2.4; [0.4, 0.6) x [13, 13.4): N 225,638, A1 0.399 at 1.9;
[0.2, 0.4) x [13.4, 13.7): 145,052, 0.249 at 3.7; [0.2, 0.4) x [13.7, 14): 67,347, 0.27 at 7.7. They drop log M_L < 13 (and
< 13.5 at z > 0.6) for "relatively high chi2 and significant detection of the background term"; the one-halo term at
10^13-10^13.5 is a 2σ Planck measurement in each redshift bin.

**Deviation from the pre-registration, recorded before the fits (10:35 UTC).** The DESI neighbour map is painted from the
groups with log M_h >= 13.0 at z < 0.9 (3-4 M) rather than >= 12.5 (23 M): the direct painting loop would not finish in the
budget. By n(M) M^1.79 the unpainted 12.5-13.0 half-dex carries about 0.4 of the painted 13.0-13.5 half-dex's y, so the
model's amplitude is ~25 % low at fixed B; B is free, and the shape (the groups' clustering times the beam) is what is used.

## Results

R-136 (resume, cloud worker cw-143, 2026-11-01 19:25 UTC on): R-132 left no stack; its container's `data/` was gone, so every
input was refetched (row D0). Code changes before any fit, all to match the pre-registration: Tully's 2MASS `_RAJ2000` columns,
the Nmb/Ntot >= 2 cut and the Tully 1,500-10,000 km/s window (the code had neither), and the ACT inpainted/subtracted masks,
whose R-131 bool caches are True = untouched (the code inverted them, which would have masked 99.8 % of the map).

| row | what | where | result |
|---|---|---|---|
| D0 | Data access (R-136, all in under 20 min) | NERSC `dr6_nilc/ymaps_20230220/` (y map, `deprojections/cib/..._1.7_10.7`, `deprojections/cib_cibdBeta/..._1.7_10.7`, beam, the three masks as in R-131, converted to bool and deleted); PLA MILCA + masks R2.01; DESI DR1 VAC `gfinder/v1.0/DESIDR9.y1.v1_group.fits` (5.6 GB, 99.6 M groups, columns extracted and deleted); HEASARC `tdat_files/heasarc_{nvss,sumss}.tdat.gz`; VizieR (IDIA mirror) J/AJ/149/171/table3 and J/A+A/596/A14/grlist_s | DESI: 22.95 M groups with log M_h >= 12.5 in the dec band, radio-flagged fraction 0.20 (NVSS/SUMSS > 5 mJy within 3'; about 0.17 of that is chance coincidence at the NVSS density). Tully: 2,421 in the dec band and window; Saulder: 29,582 |
| D1 | **Planck MILCA, DESI groups, pre-registered fit** (`fit --map milca --cat desi`) | `logs/fit` (this row) | Footprint: 3.97 M groups. Radio-masked A (obs/ref) by log M_h bin 12.5-14.5: 6.5 ± 2.3, 1.3 ± 1.3, −0.23 ± 0.61, −0.13 ± 0.23, 0.16 ± 0.19, **0.33 ± 0.08**, **0.28 ± 0.03**; neighbour amplitude B = 0.97-1.19 ± 0.12-0.14 in every bin (the catalogued neighbours on the reference relation describe the large-scale signal); χ²/dof 3-14/4. Power law over the five headline bins: log10 A_Y(10^13.5) = −7.06 ± 0.15, B_Y = 2.16 ± 0.36. At Planck's 10' the one-halo term is detected only above 10^13.5 (M500c ≳ 10^13.4), at 0.28-0.33 of the Planck 2013 reference: the level of Lim et al. 2018 (0.27 at 10^13.5) and Chen, Zhang & Yang 2022 (A1 0.25-0.40 on the same Yang DESI groups). The Planck like-for-like reproduces the literature. |
| D2 | **ACT DR6 standard ILC map, DESI groups, pre-registered fit** (`fit --map act --cat desi`) | `logs/fit_act_desi.log` | **The pre-registered model fails**: χ²/dof 120-640/7 in every bin; A radio-masked −7.1 ± 0.3, −4.2 ± 0.3, −1.31 ± 0.10, −0.39 ± 0.06, −0.10 ± 0.02, 0.12 ± 0.02, 0.218 ± 0.006; B negative (−1.8 to −4.1) below 10^13.5. Two causes, read off the profiles: (i) **the central annulus (0-1') is negative at −7 to −9 × 10^-8 in every bin up to log M_h 13.5**, with the 1-3' annuli near zero: a compact negative source at the group centre; it is −0.4 to +0.2 × 10^-8 on the CIB-deprojected map (D3), so it is dust emission of the member galaxies leaking into the standard ILC y; (ii) **the annulus means fall with radius at the 400 k mask-matched randoms too** (−0.15 at 0-1' to −5.5 × 10^-8 at 20-30'), while the 23 RA rotations of the same groups (not 15'-clean-selected) are flat at −3.5 × 10^-8: requiring a clean 15' disc puts the centre away from the subtracted and inpainted sources, whose negative wings sit in the outer annuli. The pre-registered constant c cannot absorb a radial trend, so the neighbour term takes it with a negative B. |
| D3 | ACT CIB-deprojected (1.7, 10.7 K) map, DESI, pre-registered fit (`--map act_cib`) | `logs/fit_act_cib_desi.log` | The central dip is gone (0-1' annulus −0.4, −0.6, −1.6, +0.2 × 10^-8 in the four lowest bins) but the radial trend (ii) remains and the fit fails the same way (χ²/dof 190-1230/7; B −2.5 to −3.2 below 10^13). The pre-registered "shift in A per bin" from the CIB deprojection: +2.6, +1.9, +0.8, +0.5, +0.21, +0.09, +0.04 in the seven bins, many times the errors below 10^13.5. |
| D4 | Saulder 2016 SDSS groups and Tully 2015 2MRS groups, pre-registered fit (ACT, ACT CIB, MILCA) | this row | **Saulder**: on ACT the 0-1' annulus is −26, −29, −21 × 10^-8 (4σ) in the three bins below 10^14 with the radio mask applied, −9 to −18 on the CIB-deprojected map; A(ACT CIB) = −0.50 ± 0.52, −0.23 ± 0.22, −0.10 ± 0.09, and **0.134 ± 0.032 in the cluster anchor (log M_tot 14.0-15.5, 167 groups)**; MILCA gives 0.18 ± 0.06 there. Saulder's luminosity masses put the clusters at 0.13-0.18 of the reference, below D1's 0.28 for DESI: their mass scale is not the one the template assumes (the column is "total mass", mock-calibrated), so Saulder is reported, not fitted for a relation. **Tully**: 86, 184, 68 groups pass the 15'-clean footprint (the 2MRS groups are 10-40' wide); A(ACT CIB) = −0.01 ± 1.06, 0.60 ± 0.27, 0.13 ± 0.16; no bin above 2.3σ. Neither low-z catalogue constrains the relation at this footprint. |

**Deviation, recorded 20:35 UTC before any revised fit (the D2-D3 diagnosis forces it).** The pre-registered model is kept
(A T1 + B T2 + c over the same annuli, the same templates, the same radio cut, the same rotation-null inflation) with three
changes: (1) the mean profile of 3 M mask-matched randoms (the same 15'-clean cut and dec band; `randoms` command) is subtracted
from every stack and its covariance added to the bin's; (2) the headline map is the CIB-deprojected (β = 1.7, T = 10.7 K) ACT
map, with the CIB + dβ moment-deprojected map (same SED) as its SED check and the standard map reported beside it; the
standard map's group centres are dust-dominated (D2 i) at every mass below 10^13.5; (3) the quoted error is the GLS error times
max(1, rotation inflation, sqrt(χ²/dof)). The Planck comparison gets the same random subtraction. The claim is whatever this
gives, including a null; the pre-registered numbers (D1-D4) stay as recorded.

| row | what | where | result |
|---|---|---|---|
| D5 | The revised fit (deviation above), 3 M randoms, ACT CIB-deprojected, DESI | `logs/fit_act_cib_desi_rsub.log` (first version) | Subtracting the all-random profile halves the χ² but leaves a **step at exactly 3'** in the radio-masked stacks (1-3' annuli 1-2 × 10^-8 above 3-6'), absent from the unmasked stacks: dropping groups with an NVSS/SUMSS > 5 mJy source inside 3' empties the inner annuli of the (negative-y) radio sources and leaves the outer ones at the field density. **Second deviation (21:00 UTC, after seeing this step, before any further fit): the radio-masked stacks subtract the randoms with the same radio condition (no > 5 mJy source within 3'); the unmasked stacks subtract all randoms.** |
| D6 | **ACT, the revised fit with radio-matched randoms** (`fit --rsub`), DESI, the UPP at the catalogue centre | `logs/fit_{act_cib,act,act_cibdb}_desi_rsub.log` | χ²/dof now 14-25/7 below 10^13.5, 35-90/7 above (errors inflated by sqrt(χ²/dof)); B = 0.1-0.35 ± 0.15-0.3 below 10^13.5. Radio-masked A by log M_h bin (12.5 … 14.5) on the **CIB-deprojected map: −0.69 ± 0.84, −0.11 ± 0.44, −0.24 ± 0.16, −0.01 ± 0.09, 0.03 ± 0.06, 0.157 ± 0.032, 0.235 ± 0.022**; standard map −4.5, −2.6, −1.13 ± 0.18, −0.43 ± 0.15, −0.13 ± 0.09, 0.09 ± 0.05, 0.21 ± 0.03; CIB+dβ map −3.1, −2.3, −0.93 ± 0.36, −0.34 ± 0.22, −0.05 ± 0.12, 0.15 ± 0.07, 0.25 ± 0.03. Fixed aperture (0-3' minus 3-6', CIB): −0.42 ± 0.32, −0.12 ± 0.16, 0.12 ± 0.11, 0.23 ± 0.06, 0.26 ± 0.04 from 10^13.0. Power law over the five headline bins (CIB): log10 A_Y(10^13.5) = −7.40 ± 0.13, **B_Y = 2.78 ± 0.31** (χ² 4.1/3); RICH >= 3: 2.79 ± 0.32; σ_M 0.15 / 0.35 dex: 2.80 ± 0.34 / 2.75 ± 0.27 (A moves ± 25 %, the nulls stay nulls). Unmasked (CIB): −0.30, −0.13, −0.03, 0.09, 0.21 from 10^13.0: radio hosts pull A down by 0.03-0.12. |
| D7 | Planck MILCA, the revised fit (radio-matched randoms), DESI | `logs/fit_milca_desi_rsub.log` | A from 10^13.0: 0.60 ± 0.79, 0.19 ± 0.26, 0.29 ± 0.20, **0.38 ± 0.08, 0.295 ± 0.034**; B = 0.92-1.05 ± 0.12-0.20; χ² 2-9/4; B_Y = 1.78 ± 0.25, log10 A_Y(10^13.5) = −6.87 ± 0.10. ACT/Planck per bin (no central term): 0.41 ± 0.12 (13.75-14), 0.80 ± 0.12 (14-14.5); at the power-law pivot 0.30 (+0.13 −0.09). |
| D8 | **The profile at the catalogue centre is not the UPP** (post-hoc diagnosis) | `logs/fit_*_rsub_psrc.log`, `fit_act_cib_desi_rsub_x*.log` | Even at 14-14.5 the stack is broader than the template (1-2'/0-1' = 0.59 vs 0.48). Two one-parameter descriptions: **(a) a beam-shaped central term S of free amplitude per bin** (`--psrc`): on the CIB map S = −14.6 ± 4.4, −22.8 ± 4.7, −25.9 ± 4.7, −45.6 ± 5.1 × 10^-8 (0-1' annulus units) from 10^13.25 and A = 0.76 ± 0.39, 0.52 ± 0.19, 0.39 ± 0.09, 0.36 ± 0.03 (B = 0.6-1.1, χ² 4.5-20/6); B_Y = 1.66 ± 0.20; standard map A 1.26, 0.72, 0.46, 0.40 with S −32 to −67; CIB+dβ A 1.43, 0.89, 0.69, 0.45 with S −38 to −83. **S grows with redshift at fixed mass** (14-14.5: −30.5 ± 7.3 at <z> 0.22, −56.8 ± 6.4 at 0.47 on CIB), the opposite of a fixed-luminosity radio source (flux ∝ D_L^-2, 5× brighter in the low-z half), and differs between maps by 20-40 (± 5-13): a beam-scale deficit at the catalogue centre with a partly non-tSZ spectrum. **(b) Gaussian miscentering with rms x R500 per axis**, one x per catalogue (`--xoff`): total χ² over the headline bins 214 (x = 0), 136 (0.25), 166 (0.5), 417 (0.75), against 65 for (a); A at 14-14.5 goes 0.24 → 0.30 → 0.44. Neither is pre-registered; (a) is preferred by Δχ² = 149 for five parameters. |
| D9 | Saulder and Tully | — | Not refitted for a relation (D4: their mass scales put the clusters at 0.13-0.18 of the reference, and Tully has < 200 groups per bin). |
| D10 | **Common footprint** (`fit --common`; the 1.75 M groups clean in both maps, the 357 k randoms clean in both; the random sets share positions, seed 136) | `logs/fit_{milca,act_cib}_desi_common_rsub*.log` | Planck A from 10^13.0: −0.56 ± 1.45, −0.28 ± 0.66, 0.79 ± 0.32, 0.49 ± 0.11, 0.36 ± 0.05 (B_Y 1.65 ± 0.29); ACT CIB, UPP at centre: −0.18 ± 0.14, −0.15 ± 0.11, 0.01 ± 0.08, 0.16 ± 0.03, 0.23 ± 0.02 (B_Y 2.83 ± 0.36); ACT CIB with the central term: −0.37 ± 0.81, −0.01 ± 0.57, 0.62 ± 0.24, 0.40 ± 0.06, 0.34 ± 0.02 (B_Y 1.72 ± 0.21). **ACT/Planck at 14-14.5: 0.63 ± 0.09 (UPP at centre), 0.93 ± 0.13 (central term)**; 13.75-14: 0.32 ± 0.10, 0.83 ± 0.22. (A first run without restricting the randoms to both footprints gave unstable Planck fits, B −8.6 to +5.6; superseded.) |
| D11 | **Finding check, pass 1** (fresh context-free Opus, the section alone, `blind-review.md` prompt verbatim; hints: the prior works of the trajectory check) | `finding.md`, the group section | "VERDICT: novel but not sound: the beam-scale central "deficit" that drives the ACT bracket and claim (4) is most plausibly flux-limited radio/dust residual contamination (already documented in ACT DR6×DESI stacks), the redshift test cannot exclude it, and without frequency-resolved separation the ACT one-halo amplitude below 10^14 M_sun is unconstrained." Strongest objection: the central term is most plausibly flux-limited radio and dust residual (arXiv:2607.08721, 2609.08939, 2502.08850), and the redshift test cannot exclude it, because a flux cut admits more luminous sources at higher z. **A physics point: parked after one pass.** Novelty held (no ACT DR6 tSZ Y-M on the Yang DR9 groups). |
