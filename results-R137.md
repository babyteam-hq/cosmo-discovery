# Results (R-137, cosmo-deep-13, 2026-11-01)

The mass of eRASS1 superclusters outside their member clusters, from the ACT DR6 lensing convergence integrated over each
supercluster's footprint. Scripts: `sc_members.py` (the member lists), `sc_lens.py` (apertures, member model, control, nulls,
simulations, the fit). It reuses `erass1_cmblens.py`'s halo model (Duffy concentrations, Tinker bias, CAMB linear P(k) at the
eRASS1 σ8 = 0.88, Σ_crit for a CMB source), `erass1_act.py`'s band (20 ≤ L ≤ 2048, cosine taper), mask threshold (≥ 0.99) and
posterior averaging (lognormal of the catalogue's M500_L/H width, five Gauss-Hermite nodes), and `act_dr6_fetch.sh` /
`act_dr6_simstream.sh`. Data under `data/` (gitignored); logs in `logs/` (gitignored). This file is the worker's own.

Trajectory check (2026-11-01, context-free Opus, in the dispatching session), verbatim: "VERDICT: novel: no prior work answers
this question".

## Fixed before any stack (2026-11-01 21:37 UTC, after the footprint counts, before any map was integrated in an aperture)

**Members.** The released supercluster table (`erass1sc.fits`, 1,338 systems) has each system's centre (the mean of its members'
X-ray centroids), redshift (their mean), multiplicity and mtot (= 1.46 Σ M500 / 0.39: M200 ≈ 1.46 M500 summed over members, then
divided by Chon et al. 2014's simulated cluster-mass fraction 0.39), but not the members. The member table (Liu et al. 2024,
Table 2) is at the CDS only, and the CDS resets every connection from this cloud. So `sc_members.py` rebuilds the members with the
paper's recipe (the cluster sample of its Sect. 2.2, 8,910 clusters here against the paper's 8,862; friends-of-friends on 3D
comoving distances with l = (n f)^-1/3, f = 10, l ≤ 50 Mpc, n in z ± 0.005 and in five exposure bins), then for every catalogued
system searches the nearest rebuilt groups with up to three members removed and three primary clusters added for the member
set that reproduces the catalogue's centre (to 7.6e-5 deg in RA and Dec, the 4-decimal rounding), mean redshift (to 7.6e-5) and
mtot (to 0.06 %). A system is used only if exactly one such set exists. **991 of the 1,338 are reproduced** (679 of 818 pairs, 312
of 520 with ≥ 3 members), none ambiguously. The rest are not used; the reproduced set is biased to fewer members (60 % of the
rich systems against 83 % of the pairs).

**Sample.** Reproduced systems at 0.05 < z < 0.6 whose two apertures are each ≥ 80 % clean in the ACT DR6 baseline mask (nside
512, clean = mean of the nside-4096 mask ≥ 0.99): **365 systems: 246 pairs and 119 with ≥ 3 members** (of 531 catalogued systems
in the z range with a clean centre, 210 with ≥ 3; 384 and 124 of those reproduced). Samples: ≥ 2 members (365), ≥ 3 members (119);
pairs, ≥ 3 at z < 0.2 and z > 0.2, and ≥ 5 members as information subsets.

**Apertures.** 'hull': the members' convex hull in the comoving tangent plane at the catalogue centre, dilated by 5 comoving Mpc
(a pair's hull is the segment); 'circle' (the check): a disc of 20 comoving Mpc radius around the catalogue centre. The
integral is over the aperture's clean pixels only, for the data and the model alike.

**Member model.** Each member a Navarro-Frenk-White halo truncated at r200c (Duffy et al. 2008 c(M), M500c → M200c) plus its
linear two-halo term (Tinker et al. 2010 bias × linear P(k), σ8 = 0.88, the eRASS1 cosmology), κ(L) exact in the flat sky,
multiplied by the map's band filter and back-transformed to κ(θ), averaged over the member's lognormal mass posterior, summed over
the aperture's clean pixels. K_m = Σ_members (1h + 2h). As a check, the 'all detected clusters' model: the same for every
primary-catalogue cluster with PCONT < 0.5 at any redshift whose profile reaches the aperture.

**Statistic.** K_d = Σ_clean pixels κ Ω_pix on the band-filtered (20-2048) map at nside 512. Excess
**E = Σ_s w_s (K_d,s − K_m,s) / Σ_s w_s K_m,s**, w_s = K_m,s / var(K_d,s), the variance from the first 100 of the release's
simulated reconstructions (no clusters in them, so the weights see only noise). E is the mass in the aperture outside the
members' one- and two-halo terms, in units of those terms.

**Control.** Random groupings of isolated eRASS1 clusters (the supercluster sample selection, member of no rebuilt system, in the
clean footprint: 1,574 clusters): each member replaced by an isolated cluster with |Δz| < 0.02 and |Δ ln M500| < 0.15, the
supercluster's aperture carried with it so that the member's position falls on the replacement (offsets in comoving Mpc,
re-projected at the replacement's distance), each shifted aperture ≥ 80 % clean. Summing the shifted apertures superposes N
isolated clusters at the members' positions: it carries each cluster's own one- and two-halo term and no supercluster selection.
E_control is the same E on that superposition, the same weights; ten realisations, averaged. **The result is E − E_control, per
sample and aperture, with its interval.** As a check on the control, the same with a population pool (any cluster of the
supercluster sample selection, members of other systems included).

**Errors and nulls.** Covariance and null distribution of E, E_control and E − E_control from the release's 400 simulated
reconstructions (the same pipeline). Random-position null: each supercluster's aperture (orientation kept) at a random clean
position, 40 realisations, normalised by the supercluster's own K_m. Rotated null: the whole aperture set turned in RA by
30-330° in 15° steps (21 rotations), the mask re-applied. Member-mass systematic: the eRASS1 weak-lensing mass calibration
uncertainty propagated as one correlated scale on every member mass (b_M: E − E_control scales as 1/b_M, since the control is
built on the same masses); stated with the calibration's quoted width. tSZ and CIB: the baseline map is profile-hardened; the
polarization-only map (no tSZ or CIB) and the CIB-deprojected map as checks, errors from their random-position nulls.

**Planck.** The same apertures on Planck PR3's tSZ-deprojected reconstruction (rotated to equatorial; ACT ∩ Planck clean pixels;
errors from its random-position and rotated nulls). Planck PR4 (the brief's like-for-like check) is distributed only as a GitHub
release (github.com/carronj/planck_PR4_lensing), and github.com returns 403 to this cloud; PR3 is the substitute.

**Claim rule.** The claim is the table: per sample and aperture, N, E, E_control, E − E_control with its 68 % interval. If
E − E_control is below 3σ, the result is its 95 % upper limit (and lower limit), written as such.

| # | test | script | statistic | significance | what it means |
|---|---|---|---|---|---|
| D0 | Data access. eRASS1 primary v3.2, cosmology v1.1, the optical follow-up catalogue (Kluge et al. 2024) and the supercluster table from erosita.mpe.mpg.de (110 MB, under a minute); ACT DR6 lensing v1 from NERSC (baseline, polonly, cibdeproj κ alm, 96 MB each; the baseline mask, 805 MB SZ-cluster mask; seconds to two minutes each); the 400 baseline simulations streamed one at a time (`act_dr6_simstream.sh`, 160 MB each, about 5 per minute with two consumers); Planck PR3 tSZ-deprojected lensing from the PLA (320 MB, under a minute). Not reachable: the CDS (the supercluster member table; connection reset), github.com (Planck PR4; 403) | `sc_members.py`, `sc_lens.py` | — | — | everything but the member table and PR4 in hand within 20 minutes of the start; the members rebuilt (above) |
| D1 | Members rebuilt | `sc_members.py`, `sc_members.py refine 3 3` (`logs/refine33.log`) | 991 of 1,338 catalogued systems reproduced exactly (centre, redshift, multiplicity, mtot); FoF alone 780, the constrained search 211 more; 0 ambiguous | — | the members used are the catalogue's, to its rounding, for three quarters of it |
| D2 | Footprint and sample | `sc_lens.py prep` (`logs/sc_prep.log`) | 365 systems (246 pairs, 119 with ≥ 3); isolated pool 1,574; 72,782 apertures (the two apertures, 10 control and 10 population realisations, 40 random-position and 21 rotated nulls); median clean hull 1.0 deg², circle 3.7 deg²; in the member model the one-halo term is 32 % of the hull integral and 15 % of the circle's (the members' two-halo terms dominate the footprint) | — | — |
| D3 | **A look before the weights were fixed (a deviation, stated)**: `fit` was run once with weights from the 28 simulations then processed (the pre-registration says the first 100), on the baseline map | `sc_lens.py fit` | ≥ 2 members, hull: E − E_control = +0.32 ± 0.38; circle −0.29 ± 0.43; ≥ 3 members: +0.24 ± 0.59, −0.03 ± 0.69 | < 1σ | the precision is ±0.4-0.7 on the excess in units of the member model, so the result will be an interval with no detection; the statistic, samples and apertures were not changed after this look. After it (and so as information rows only), two aperture-size variants were added: hull margin 2.5 / circle 10 cMpc and hull margin 10 / circle 30 cMpc |
| D4 | **The measurement** (the statistic fixed above; baseline map; weights frozen from simulations 0001-0100, `weights_var100.npz`; errors the spread over the 400 simulations; member model) | `sc_lens.py fit` (`data/erass1_sc/fit_baseline.txt`, `logs/fit_base_baseline.log`) | **≥ 2 members (365): hull E = +0.46 ± 0.24, E_control = +0.17 ± 0.28, E − E_control = +0.29 ± 0.36, 95 % [−0.42, +1.00]; circle −0.25 ± 0.25, +0.06 ± 0.32, −0.31 ± 0.40, [−1.08, +0.46]. ≥ 3 members (119): hull +0.30 ± 0.34, +0.03 ± 0.38, +0.27 ± 0.49, [−0.69, +1.23]; circle −0.18 ± 0.38, −0.13 ± 0.43, −0.05 ± 0.56, [−1.15, +1.05]**. The aperture integral itself: K_d/K_m = 1.46 (hull, 6.2σ above zero), 0.75 (circle, 3.0σ); ≥ 3: 1.30 (3.9σ), 0.82 (2.2σ). Information subsets (hull / circle E − E_control): pairs +0.32 ± 0.47 / −0.55 ± 0.48; ≥ 3 at z < 0.2 (49) +1.14 ± 1.22 / +0.74 ± 1.38; ≥ 3 at z > 0.2 (70) +0.08 ± 0.53 / −0.23 ± 0.61; ≥ 5 (23) +0.11 ± 0.84 / −0.20 ± 0.97 | every row < 1.2σ | no excess of supercluster lensing beyond the members' one- and two-halo terms relative to matched isolated clusters; the 95 % upper limits are 1.0 (hull) and 0.5 (circle) times the member model for ≥ 2 members, 1.2 and 1.1 for ≥ 3 |
| D5 | Derived: the members' one-halo (r < r200c) share of the measured aperture integral (weighted, the model's 1h / (1 + E), its range from E's error) | same | ≥ 2, hull 0.24 [0.20, 0.28] (control superposition 0.29 [0.23, 0.38]); circle 0.22 [0.16, 0.33] (control 0.15 [0.12, 0.21]); ≥ 3, hull 0.22 [0.17, 0.29] (0.27 [0.20, 0.42]); the model alone puts it at 0.35 (hull) and 0.16 (circle) | — | about three quarters of the band-filtered lensing inside the 5-Mpc hull lies outside the members' r200c, for superclusters and for the isolated-cluster superposition alike; most of that is the members' own two-halo term in the model |
| D6 | Nulls: the simulations' mean of each E; random-position (40) and rotated (21 turns; 5-7 keep half the sample) apertures on the data map, normalised by each supercluster's own K_m | same | sims mean of E: 0.000 ± 0.001 in every row; random null, ≥ 2 hull +0.05 ± 0.23, circle +0.01 ± 0.27; ≥ 3 +0.02 ± 0.35, −0.06 ± 0.39; rotated +0.14 ± 0.17 (6), +0.01 ± 0.27 (6); ≥ 3 +0.10 ± 0.23, −0.15 ± 0.35 | random and rotated nulls within 1σ of zero; their spreads equal the simulation errors (0.23 vs 0.24 hull, 0.27 vs 0.25 circle) | the error is the noise, and the map has no mean offset at random apertures; no null-mean correction is applied |
| D7 | Control check: the population pool (any cluster of the supercluster sample selection) in place of the isolated pool | same | E_pop: ≥ 2 hull +0.23, circle +0.22; E − E_pop +0.24 ± 0.27, −0.48 ± 0.29; ≥ 3 +0.14 ± 0.43, −0.55 ± 0.49 | < 1.7σ | the population control sits 0.06-0.16 above the isolated one (isolated clusters have, as selected, less company), within the errors; the mean cluster carries +0.2 ± 0.3 beyond the linear-two-halo model in these apertures, the direction of R-122's free two-halo amplitude 1.2-1.5 |
| D8 | The 'all detected clusters' model (the members' 1h + 2h plus the one-halo term of every other primary cluster with PCONT < 0.5, any z, in the aperture) | same | ≥ 2 hull E − E_control = +0.43 ± 0.30, circle −0.10 ± 0.30; ≥ 3 +0.44 ± 0.41, +0.15 ± 0.43 | ≤ 1.4σ | removing the other detected clusters' halos changes nothing outside the errors |
| D9 | Foregrounds: polarization-only map (no tSZ, no CIB) and CIB-deprojected map, same apertures and weights; errors the baseline simulations' scaled by the random-null spread ratio | `fit --map polonly / cibdeproj` | polonly: ≥ 2 hull −0.28 ± 0.84, circle −0.26 ± 0.60; ≥ 3 −0.58 ± 0.91, +0.39 ± 0.81. cibdeproj: +0.31 ± 0.34, −0.19 ± 0.34; +0.34 ± 0.49, −0.16 ± 0.52 | < 1σ | the baseline (profile-hardened) and CIB-deprojected maps agree to 0.02-0.12; polarization alone agrees within its 0.6-0.9; no tSZ or CIB leakage is visible at this precision, which bounds it at roughly the statistical error |
| D10 | Planck PR3 tSZ-deprojected map (rotated to equatorial), apertures wholly inside the Planck mask (65 %: 207 superclusters, 61 with ≥ 3); errors scaled as in D9 | `fit --map planck` | ≥ 2 hull +0.45 ± 0.99, circle +0.62 ± 0.73; ≥ 3 +0.30 ± 1.22, +1.10 ± 1.36 | < 1σ | the same null at two to three times the error; PR4 not reachable (github.com 403) |
| D11 | Aperture-size variants, defined after the first look (information only): hull margin 2.5 / circle 10 cMpc, and 10 / 30 cMpc; their own footprint selection, weights from their own first 100 simulations | `SC_VAR=2.5_10`, `10_30` | 2.5/10: ≥ 2 hull −0.02 ± 0.33, circle +0.49 ± 0.48; ≥ 3 +0.23 ± 0.47, +0.57 ± 0.72. 10/30: ≥ 2 +0.11 ± 0.37, −0.09 ± 0.41; ≥ 3 +0.25 ± 0.49, +0.23 ± 0.54 | < 1σ | every aperture size from 2.5 to 30 cMpc gives the same null; the sign changes between apertures are within the errors |
| D12 | Member-mass systematic: one scale b_M on every member mass; the model's 1h term scales as b_M and its 2h term as the bias, b ∝ M^0.34-0.46 (the grid's Tinker slope at 1-3 × 10^14, z 0.1-0.3), so E − E_control scales as 1/(h b_M + (1 − h) b_M^0.4), h the 1h share (0.35 hull, 0.16 circle) | arithmetic on D4 | b_M = 0.9 / 1.1: × 1.07 / × 0.94 (hull), × 1.05 / × 0.95 (circle); b_M = 0.70 (R-122's ACT value): × 1.24 / × 1.19. 95 % intervals with ± 10 % on b_M: ≥ 2 hull [−0.45, +1.07], circle [−1.14, +0.49]; ≥ 3 [−0.74, +1.31], [−1.21, +1.11] | — | the mass calibration moves the interval by 5-7 %, a fifth of its width; at R-122's lower mass scale the upper limits rise by a fifth (hull ≥ 2: +1.24) |
| D13 | The finding section through the blind gate (the section alone in a file, the prompt of `blind-review.md` verbatim, fresh context-free Opus), first pass | Agent tool, Opus | referee's verdict (verbatim): "VERDICT: novel but not sound: the central statistic E − E_control is compared with zero, but no ΛCDM null has been computed for it, because the overlapping two-halo terms and the isolated (underdense) control mean zero isn't the expected value, so the bounds can't be read as limits on mass outside the member clusters." Closest prior work it found: arXiv:2203.11306 (Planck κ on DES Y3 superclusters, 8.4σ), 2105.11936 (DESI LS superstructures × Planck κ), 2409.04535, 1805.04555, 2310.04387, 2404.17345 | a miss | physics (the comparison's ΛCDM expectation needs a light-cone with the eRASS1 selection), so parked per the brief, no second pass; the title corrected to "bounds" after the verdict. What stands: the aperture integral of 365 eRASS1 superclusters is detected at 6.2σ (hull) in ACT DR6, and supercluster and isolated-cluster superpositions agree within ±0.36-0.56 of the member model in every pre-registered row |
