# Results (R-140, cosmo-deep-13b, 2026-11-01/02)

The ΛCDM null for R-137's statistic (E − E_control: the ACT DR6 convergence integrated over eRASS1 superclusters' apertures
beyond their members' one- and two-halo model, against matched isolated-cluster superpositions), from a simulated light-cone.
R-137's gate verdict, verbatim: "VERDICT: novel but not sound: the central statistic E − E_control is compared with zero, but no
ΛCDM null has been computed for it, because the overlapping two-halo terms and the isolated (underdense) control mean zero isn't
the expected value, so the bounds can't be read as limits on mass outside the member clusters." This file computes that
expectation. Script: `sc_null.py` (imports `sc_lens.py`'s apertures, member model, control and statistic unchanged, and
`sc_members.py`'s finder). Data under `data/` (gitignored); logs in `logs/` (gitignored). This file is the worker's own.

## Fixed before running (2026-11-02 00:00 UTC, after the selection was matched and before any simulated map was integrated)

**Simulation.** Websky (Stein et al. 2020, arXiv:2001.08787; Ωm 0.31, h 0.68, σ8 0.81), from its CITA data site
(mocks.cita.utoronto.ca/data/websky/v0.0): the full-sky CMB-lensing convergence `kap.fits` (nside 4096, all structure to the
last-scattering surface) and the halo light-cone `halos.pksc` (862,923,142 halos sorted by mass; only the first 30 M are read,
over HTTP ranges: M200m ≥ 1.6 × 10^13 M☉, complete for every mass the selection reaches). Reachable from the cloud (the
directory listing is 403, the files 200). Websky's own frame is read as equatorial, as ACT's Websky foreground reconstructions do.

**Patches.** 25 copies of the simulated sky: patch 0 the Websky frame as is, patches 1-24 random rotations of the whole sky
(halos and convergence together; `scipy` `Rotation.random(1000 + p)`). Each is put through the real footprint: the eRASS1
footprint of the finder's input sample (the nside-16 pixels holding at least one of its 8,910 clusters, 15,376 deg²) and the ACT
DR6 baseline mask (nside 512, clean = ≥ 0.99), exactly as R-137. Each patch is then a data-sized sample; the spread over patches
is the cosmic variance of the statistic for a sample of the data's size (two random footprints overlap by ~10 % of their area).

**Simulated clusters.** Halo mass M200m = 4/3 π ρ_m R_TH³ (Websky's definition) → M500c through NFW with Duffy et al. 2008's
c200m(M, z). Observed redshift: cosmological + line-of-sight velocity + a Gaussian error with σ drawn from the eRASS1 finder
sample's BEST_ZERR at that redshift (0.02 bins). Selection: probability Φ((ln M500c − ln M_c(z)) / s), with M_c(z) per 0.02 in
observed z set so the expected number density equals the eRASS1 finder sample's (8,910 clusters in 15,376 deg²), and s = 0.6
chosen from {0, 0.2, ..., 1.0} as the value that best matches the sample's median M500 per bin at 0.05 < z < 0.6
(`logs/select.log`: a sharp threshold, s = 0, puts the simulated medians 1.3-1.6 times the data's; s = 0.6 brings them to
0.97-1.02 at 0.15 < z < 0.6, 0.9-1.2 at z < 0.15). Full sky: 23,680 selected against the data's 23,905 full-sky equivalent.
The released eRASS1 selection function (Clerc et al. 2024) is a function of count rate and exposure, which would need a
simulated X-ray luminosity for each halo; the matched mass-redshift selection is the brief's fallback, stated. Each simulated
cluster's model mass is its true M500c, with a posterior width drawn from the catalogue's at that redshift (the model is
averaged over it exactly as R-137 averages the catalogue's).

**Finder.** R-137's rebuild of Liu et al. 2024's friends-of-friends (`sc_members.py`): 3D comoving positions (H0 70, Ωm 0.3),
l = (n f)^-1/3 with n from the clusters in z ± 0.005 over the catalogue's hemisphere-area convention, f = 10, l ≤ 50 Mpc, the
'min' rule (R-137's best reproduction). One exposure bin instead of R-137's best five: the simulated selection does not depend on
exposure, so the five bins would have equal densities (the data's median l is 50 Mpc, the cap, either way). Every system with
≥ 2 members is used (the data used the 991 of 1,338 catalogued systems whose members were reproduced exactly; that set is biased
to pairs, so a multiplicity-matched subsample of the simulated systems is an information row).

**Sample, apertures, member model, control, statistic, weights.** R-137's, imported unchanged from `sc_lens.py`: systems at
0.05 < z < 0.6 with both apertures ≥ 80 % clean; the hull (5 cMpc margin) and the 20-cMpc circle; NFW (r200c) + linear two-halo
member model at the eRASS1 cosmology (σ8 0.88), posterior-averaged, band-filtered 20-2048 exactly as the map; ten control
realisations from the isolated pool (selected clusters in no system, centre clean) matched in |Δz| < 0.02 and |Δ ln M| < 0.15;
E = Σ w (K_d − K_m) / Σ w K_m, w = K_m / var(K_d), the variance from the ACT DR6 simulated reconstructions 0001-0100 on each
patch's own apertures (R-137's rule). Samples: ≥ 2 members and ≥ 3 members (the claim's); pairs and the multiplicity-matched
subsample (information).

**Maps.** The simulated convergence: `kap.fits` → alm to L = 2048 → R-137's band filter → nside 1024 in the Websky frame, read
in each patch's frame by bilinear interpolation at the nside-512 pixel centres (R⁻¹ applied to each pixel's direction).
Noise: the ACT DR6 noise realisations n_i = recon_i − input_i (the release's baseline simulated reconstructions minus their
input convergence, `sim_inputs/kappa_alm`), i = 1-100, band-filtered identically, added to the noiseless integrals (the
statistic is linear in the map, so this is the same as adding the maps). Foregrounds (patch 0 only, where the frames agree):
ACT's own foreground-only reconstruction from Websky's tSZ, CIB and radio skies through the DR6 baseline pipeline
(`websky/kappa_fg_freqcoadd_prh.fits` in the ACT DR6 lensing release) added as a third term.

**The null and the comparison (the claim rule).** Per sample and aperture: the ΛCDM expectation ⟨E − E_control⟩ = the mean
over the 25 noiseless patches; its distribution for a data-like sample = E − E_control over the 25 × 100 (patch, noise)
pairs, quoted as the mean, the standard deviation (split into the patch-to-patch and the noise parts) and the 2.5-97.5 %
range. The comparison: Δ = (E − E_control)_observed (R-137, D4) − ⟨E − E_control⟩_ΛCDM, with σ = that distribution's
standard deviation; the observed value's percentile in the distribution; the 95 % interval Δ ± 1.96 σ, which is the
outside-member excess relative to ΛCDM in units of the member model. The same for E and E_control separately.
**The null is meaningful if** the simulated patches' sample sizes are within ± 30 % of the data's (365 and 119) and the
selected medians within 15 % of the data's at 0.15 < z < 0.6 (met by the selection, above); otherwise the result is that it is not.
Known difference, stated before: Websky's σ8 = 0.81 against the member model's 0.88, so the model's two-halo term overstates
the simulation's by ~(0.88/0.81)² ≈ 1.18 at fixed bias; E and E_control share that offset, and E − E_control carries it only
through the supercluster-to-isolated difference in two-halo share.

| # | test | script | statistic | significance | what it means |
|---|---|---|---|---|---|
| D0 | Data access. Websky `kap.fits` (805 MB, a minute) and the first 30 M halos of `halos.pksc` (1.2 GB of 34.5 GB by HTTP range, a minute) from mocks.cita.utoronto.ca; the eRASS1 primary, cosmology, optical and supercluster tables from erosita.mpe.mpg.de (seconds); from the ACT DR6 lensing release at NERSC, the baseline mask, simulations 0001-0100 with their input κ (`sim_inputs/kappa_alm`, 316 MB a pair, ~10 pairs a minute with three streams) and ACT's own Websky foreground-only reconstruction (`websky/kappa_fg_freqcoadd_prh.fits`) | `sc_null.py halos`, `noise` | — | — | everything in hand within 15 minutes; no substitute needed |
| D1 | Simulated selection (`logs/select.log`) | `sc_null.py select` | s = 0.6; 23,680 selected over the full sky (data: 23,905 full-sky equivalent); per 0.02 in z the density matches by construction; median M500 sim/data 0.91 (z 0.05), 1.23 (0.11), 1.01 (0.19), 0.98 (0.29), 0.93 (0.39), 0.95 (0.49), 0.94 (0.59). A sharp threshold (s = 0) gives 1.05-1.40 | — | the selection is matched in number and in median mass to ≤ 7 % at 0.15 < z < 0.6 (the pre-registered 15 %), 9-23 % at z < 0.15 |
| D2 | Patches (`logs/prep_*.log`) | `sc_null.py prep --patch 0..24` | per patch: 8,818 clusters in the eRASS1 footprint (data 8,910), 1,254 FoF systems at all z (catalogue 1,338), sample 354 (311-398; data 365): 224 pairs (246), 129 with ≥ 3 (107-149; 119), 31 with ≥ 5 (23); isolated pool 1,686 (1,574); 0 failed control draws; the model's one-halo share 0.29 (hull), 0.16 (circle) against the data's 0.35, 0.16 | — | the simulated samples match the data's sizes to ± 15 % (pre-registered ± 30 %): the null is meaningful by the rule fixed before. The simulated systems are a little richer (no member-reproduction cut); the multiplicity-matched rows test that |
| D3 | Patch independence | inline (`sc_null.py` geometry) | footprint (eRASS1 finder footprint ∩ ACT clean) 4,865 deg²; the 25 patches' union covers 95 % of the sky, each point ~3.1 times; mean pairwise overlap 12 % (max 56 %) | — | ~8 independent patches: the error on the ΛCDM mean uses N_eff = 8, not 25 (fixed after this count, before the comparison was read) |
| D4 | **The ΛCDM null** (noiseless Websky through the pipeline; mean over 25 patches; ± its error with N_eff 8; then one data-sized sample's distribution over 25 × 100 (patch, ACT noise) pairs: sd, split, 95 % range) | `sc_null.py integrate`, `fit` (`data/sc_null/fit_null.txt`) | **≥ 2, hull: E +0.007, E_control −0.100, E − E_control = +0.107 ± 0.081; one sample sd 0.363 (patches 0.230, noise 0.278), 95 % [−0.59, +0.83]. ≥ 2, circle: −0.135, −0.255, +0.120 ± 0.089; sd 0.368 (0.252, 0.266), [−0.56, +0.87]. ≥ 3, hull: −0.111, −0.240, +0.129 ± 0.122; sd 0.495 (0.345, 0.353), [−0.81, +1.12]. ≥ 3, circle: −0.139, −0.343, +0.204 ± 0.166; sd 0.603 (0.470, 0.380), [−0.96, +1.35].** Information: pairs +0.080 ± 0.061 (hull), +0.047 ± 0.108 (circle); multiplicity-matched mix (R-137's 246 : 96 : 23) +0.060 ± 0.078, +0.080 ± 0.082 | null ≠ 0 at 1.3σ (≥ 2 hull) to 1.3σ (circle); every row positive | ΛCDM does not predict zero: superclusters carry +0.1 to +0.2 of their member model more than matched isolated-cluster superpositions, the isolated control sitting lower (E_control −0.10 to −0.34 against E +0.01 to −0.14), as the R-137 referee argued; but the offset is a third of one sample's spread. One data-sized sample scatters by ± 0.36-0.60 under ΛCDM, the same as R-137's error (0.36-0.56): its Gaussian-simulation errors were right in size, the non-Gaussian sample variance (0.23-0.47) taking the place of the Gaussian κ in its simulations |
| D5 | **The comparison**: R-137's observed E − E_control (D4 there) minus the ΛCDM mean; σ = one sample's sd ⊕ the mean's error | same | **≥ 2, hull: Δ = +0.18 ± 0.37 (+0.5σ; percentile 69), 95 % [−0.55, +0.91]; ≥ 2, circle: −0.43 ± 0.38 (−1.1σ; 12), [−1.17, +0.31]; ≥ 3, hull: +0.14 ± 0.51 (+0.3σ; 61), [−0.86, +1.14]; ≥ 3, circle: −0.25 ± 0.63 (−0.4σ; 35), [−1.48, +0.97].** Pairs: +0.24 ± 0.45, −0.60 ± 0.46. E and E_control separately: observed E +0.46 against ΛCDM +0.01 ± 0.22 (+2.0σ, hull), −0.25 against −0.14 ± 0.26 (−0.4σ, circle); E_control +0.17 against −0.10 ± 0.31 (+0.9σ), +0.06 against −0.26 ± 0.34 (+0.9σ) | every row within 1.3σ of the ΛCDM expectation | the eRASS1 superclusters carry no mass outside their members beyond what ΛCDM puts there: relative to ΛCDM, the excess is below 0.91 (hull) and 0.31 (circle) of the member model for ≥ 2 members, 1.14 and 0.97 for ≥ 3 (95 %); a deficit larger than 0.55 (hull) and 1.17 (circle) is excluded. The limits move down by 0.1-0.2 against R-137's (read against zero) because the null is positive. The data's hull E sits 2.0σ above the simulation's, its E_control 0.9σ: a common offset (masses or σ8) that the difference removes |
| D6 | Hull against circle (the referee's 'sign flip') | inline on `fit_null.npz` | ΛCDM: D_hull − D_circle = −0.00 ± 0.33 (≥ 2; the two correlate at 0.60), observed +0.60, 1.9σ (percentile 97); ≥ 3: −0.07 ± 0.46, observed +0.32, 0.8σ | 1.9σ | the null has no aperture-dependent offset (+0.11 and +0.12): the data's hull-circle difference is a 2σ fluctuation, not the null's shape; it is stated, not used |
| D7 | Foregrounds (patch 0, ACT's Websky tSZ + CIB + radio foreground-only reconstruction through the DR6 baseline pipeline, added) | `fit` | shift in E − E_control: −0.005 (≥ 2 hull), +0.013 (circle), +0.046 / +0.049 (≥ 3), −0.073 / −0.020 (pairs); in E alone −0.02 to −0.11 | ≤ 0.2σ | residual foregrounds after profile hardening move the difference by ≤ 0.07, a fifth of its error or less; they lower E slightly (tSZ-induced deficit) in both E and E_control |
| D8 | The model's two-halo term at Websky's own σ8 (× (0.81/0.88)² = 0.847, information) | `SC_NULL_K2SCALE=0.847 fit` | E − E_control null: +0.117 (≥ 2 hull), +0.137 (circle), +0.144, +0.233 (≥ 3); E rises to +0.125 / −0.007 | — | the null is insensitive to the σ8 mismatch (shifts ≤ 0.03); E itself moves by +0.12 |
| D9 | Selection variant: a sharp threshold (s = 0, the brief's literal 'mass-redshift threshold'; medians 1.05-1.40 × the data's), own patches and noise (information, defined after D4) | `SC_NULL_VAR=s0` | null: +0.20 ± 0.09 (≥ 2 hull), +0.16 ± 0.11 (circle), +0.26 ± 0.15, +0.21 ± 0.17 (≥ 3); Δ: +0.09 ± 0.36, −0.47 ± 0.39, +0.01 ± 0.52, −0.26 ± 0.58 | ≤ 1.2σ | a sharper, more massive selection raises the null by ~0.1 (the isolated control falls further: −0.17 to −0.40); the comparison moves by ≤ 0.13 and stays within 1.2σ |
| D10 | Mass variant: the model and the control matching use a noisy mass, M500 × exp(σ_post N(0, 1)), as the catalogue's are estimates (information, defined after D4) | `SC_NULL_VAR=mnoise` | null: +0.11 ± 0.06, +0.14 ± 0.08, +0.12 ± 0.11, +0.17 ± 0.15; Δ: +0.18 ± 0.33, −0.45 ± 0.35, +0.15 ± 0.46, −0.22 ± 0.57 | ≤ 1.3σ | mass-matching noise in the control does not move the null (≤ 0.03) |
| D11 | The revised R-137 section of `finding.md` (the ΛCDM null and the comparison added) through the blind gate: the section alone in a file, the prompt of `blind-review.md` verbatim, fresh context-free Opus, hints as R-137's plus "anyone who computed supercluster or cluster-pair lensing expectations in simulations"; pass 2 of the finding | Agent tool, Opus | referee's verdict (verbatim): "VERDICT: the findings are confirmed novel and sound if the calculations hold" | **a pass** | strongest objection: an environment-dependent bias in the count-rate masses is not in the simulated null (sized by the referee at ~0.1, a quarter of σ). Corrections after the verdict, none to a number: the title no longer says superclusters "carry no mass" and quotes the expected sensitivity ± 0.75; "first lensing measurement on the eRASS1 superclusters" narrowed (arXiv:2404.17345 weak-lensed eRASS1 supercluster members against isolated clusters); two weakest points added. Closest prior work it found: arXiv:2203.11306, 2404.18278, 2105.11936, 1709.02543, 2507.08561, 1402.3302, 1702.08485, 2107.05523, 2409.04535, 2404.17345 |
