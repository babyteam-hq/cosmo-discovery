# Analysis plan: does the Giant Arc reappear on independent DESI DR1 sightlines?

## Question
The Giant Arc (GA) is an arc ~1 Gpc long made of Mg II absorbers at z ~ 0.80, found by single-linkage clustering of absorbers
in SDSS quasar spectra (claimed delta rho/rho = 1.3 +- 0.3 inside its convex hull of member spheres, CHMS). One reading says
it is a real overdensity. The other says it is a pattern the linking algorithm picks out of an ordinary (clustered) point
field. A real overdensity must also show up as an excess of absorbers on quasar sightlines SDSS never observed. A chance
pattern must show no more than ordinary absorber clustering gives.

## Template (frozen before any DESI data was opened)
Zhu & Ménard SDSS DR7 + DR12 "trimmed" Mg II catalogues, merged (a DR7 absorber within 1.5" and |dz| < 0.005 of a DR12 one is
dropped). Slice z = 0.802 +- 0.060 (4,752 absorbers). Single linkage at 95 Mpc in comoving coordinates (flat LCDM, Om 0.27,
h 0.7, the discovery paper's choice). The GA is the group containing the absorber nearest (RA, Dec) = (216, 39): 46 members,
RA 205.0-224.6, Dec 36.2-42.2, z 0.745-0.860, member W2796 median 0.82 A. Template volumes: the union of spheres of radius
R = 47.5, 70 and 95 Mpc around the 46 members, and the CHMS (convex hull of spheres of radius half the mean minimum-spanning-tree
edge length around each member).

## Sightlines and absorbers
- Sightlines: DESI DR1 main-survey dark-time quasars (the quasar catalogue the DESI Mg II absorber catalogue searched), one
  per object; exposure depth TSNR2_QSO; Legacy Surveys g and r fluxes; redshift z_q.
- Independence: every DESI quasar within 1.5" of an SDSS DR16Q quasar or of any quasar in the Zhu & Ménard searched lists is
  removed (322,606 of 1,556,507). 1,233,901 independent sightlines remain. No absorber on an SDSS sightline enters.
- Absorbers: DESI DR1 Mg II absorber catalogue (270,529; the physically impossible z_abs > z_q + 5000 km/s are a separate
  file), attached to sightlines by target id.
- Path per sightline: z_abs from max(0.30, 1250 (1 + z_q)/2796 - 1) (the catalogue does not search the Lya forest; dN/dz at
  z ~ 0.8 drops to exactly zero for z_q > 3.4) to (1 + z_q)(1 - 5000 km/s / c) - 1 (proximate zone removed).
- Path inside a template: exact ray-sphere (or ray-convex-hull) intersection in comoving distance, intervals merged per
  sightline and clipped to that sightline's allowed path. The template's spheres extend in redshift beyond the slice; the
  path is used wherever it lies (z ~ 0.70-0.90).

## EW matching
DESI EW_2796 is rest-frame. On 21,789 absorbers present in both catalogues (DESI re-observations of SDSS sightlines), DESI
EW / Zhu & Ménard EW = 1.03 (median; 1.02-1.05 in EW bins 0.3-0.6, 0.6-1, 1-2, 2-5 A; 1.03 in the slice). Cuts W2796 > 0.3 A
and > 0.6 A on the Zhu & Ménard scale = DESI EW_2796 > 0.309 and > 0.618 A, applied to both the counts and the dN/dz.

## Expected count and completeness
Detection completeness in DESI DR1 depends overwhelmingly on the quasar's brightness: all-sky dN/dz (W > 0.3, z = 0.8) runs
0.00004 / 0.001 / 0.011 / 0.049 / 0.280 across g-flux quintiles, against 0.052 / 0.063 / 0.090 across TSNR2 tertiles.
Removing SDSS quasars removes the brighter ones inside the SDSS footprint, so the GA's independent sightlines are fainter
than the sky average. Therefore dN/dz(z) is measured per stratum (TSNR2 tertile x 12 g-flux quantile bins, finest at the
bright end: percentiles 30, 50, 60, 70, 77, 83, 88, 92, 95, 97.5, 99 x 4 z_q bins split at 1.5, 2.2, 2.8; 144 strata), in z
bins of 0.01 over 0.55-1.10 smoothed with a 5-bin boxcar, on all 1.23M independent sightlines over the whole sky. A
sightline's expected count is the integral of its own stratum's dN/dz over its path inside the template. Check: over the
whole sky in the slice the stratified expectation is 7,844 against 7,864 observed.

## Statistic
delta = N_obs / N_exp - 1, per template (R = 47.5, 70, 95 Mpc; CHMS) and EW cut (0.3, 0.6 A). Under Poisson counts with mean
(1 + delta) e_i on sightline i, the maximum-likelihood delta is sum n_i / sum e_i - 1, so the stratified ratio is the
optimal TSNR2- and flux-dependent weighting. A literal TSNR2-weighted ratio, sum TSNR2_i n_i / sum TSNR2_i e_i - 1, is
reported as a variant.

**Primary cell (fixed in advance): R = 95 Mpc, W > 0.3 A.** It has the largest expected count and the best separation of
the two readings. The other seven cells are reported in full as secondary.

Expected counts at the GA (W > 0.3): 4.13 (R = 47.5), 9.64 (70), 16.39 (95), 8.86 (CHMS), through 3,534 / 4,807 / 6,282 /
3,823 sightlines.

## Nulls
1. Template shift (coverage, completeness and large-scale-structure variance): the 46-member template is rotated rigidly on
   the sky to a random centre (drawn from independent DESI sightline positions) with a random spin, redshifts kept. A shift
   is kept only inside the SDSS footprint (> 30 % of DESI quasars within 1 deg of every 5th shifted member are SDSS quasars;
   the GA itself has 39-53 %), so the SDSS exclusion thins the sightlines as it does on the GA. It is dropped if a shifted
   sphere centre lies within 2R of a GA member (hull: centre within 2x the hull radius), and kept only if its N_exp is
   within 0.75-1.33x the GA's, so the scatter is the GA's own at its depth. 3,000 trials. Feasibility (600 trials, ~370
   kept): mean -0.045 / -0.052 / -0.073 / -0.072, sd 0.50 / 0.32 / 0.24 / 0.34 (R = 47.5 / 70 / 95 / CHMS).
2. Clustering baseline (what a chance pattern built from real absorbers gives): spheres of radius R on random sets of 46
   Zhu & Ménard slice absorbers (GA members excluded; only absorbers with at least one independent DESI sightline through
   their sphere), 1,000 sets. Feasibility: mean +0.092 / +0.013 / -0.010 (R = 47.5 / 70 / 95), so the clustering excess
   above random positions is +0.14 / +0.07 / +0.06.
3. Other linkage groups (secondary): every other single-linkage group from the same slice and recipe with >= 20 members,
   evaluated the same way: whether the algorithm's groups in general carry an excess on new sightlines.

## Claim
Delta = delta_GA - mean(clustering null), per cell, with a 68 % and a 95 % interval from the template-shift null's
percentiles around delta_GA (its spread holds the Poisson, coverage and LSS variance at the GA's depth) combined in quadrature
with the clustering-null mean's standard error. Significance of Delta against zero and against the real-structure prediction.

## Predictions
- Real structure with the published amplitude: Delta ~ +1.3 in the CHMS cell, and >= +1.3 in the sphere cells (spheres hug
  the members more tightly than the hull). The published amplitude is measured on its own selected members, so a real
  structure's true amplitude can be lower (selection bias); any Delta clearly above zero but below 1.3 means "real but
  weaker than claimed".
- Chance pattern: Delta ~ 0 (delta_GA ~ the clustering baseline).
- Separation at the GA's depth: 5.4 sigma (R = 95), 4.1 (70), 3.8 (CHMS), 2.6 (47.5).

## Known limits
- One DESI release; the GA's sightlines are fainter than average, which raises the variance (reflected in the shift null).
- The Big Ring and Giant Ring are not tested (no reproducible member list).
- The rough count seen at proposal time (flat dN/dz, partial exclusion) was biased by the completeness gradient above.
