# Results (R-147, cosmo-deep-16b, 2026-11-02)

Question: does type-dependent photo-z isolation leakage explain the early/late split in the KiDS lensing radial
acceleration relation (RAR)? Brouwer et al. 2021 (arXiv:2106.11677) find log⟨g_obs,E/g_obs,L⟩ = 0.17 dex (Sérsic split) and
0.27 dex (u−r split) over the full g_bar range, 0.14 and 0.19 dex at g_bar > 10⁻¹³ m s⁻² (R < 0.3 Mpc/h70). They assumed photo-z
errors affect the isolation of both types equally (their Sec. 5.3, third caveat). Their App. A found, in MICE, ~30 % more ESD at
R > 0.3 Mpc/h70 for the photo-z-isolated sample than for the truly isolated one, for the whole sample. Mistele et al. 2024
(arXiv:2310.15248, Sec. 4.2 and App. C.5) find the ETG RAR moves with R_isol and the LTG RAR does not; with 4 Mpc/h70 isolation
and their ETG masses (Q_ETG = 1.4, hot gas) the split disappears.

Trajectory verdict on file (brief, 2026-11-01): "partly novel": contamination by type against spectroscopic truth, and the
spurious offset it predicts against 0.14-0.19 dex, are open. Scripts: `r147_leak.py` (this run). Data under `data/r147/`
(gitignored): the BGS group VAC (Williams et al. 2026, 290 MB) and the DESI DR1 LSS `BGS_BRIGHT_full_HPmapcut` catalogue (5.2 GB,
for the Legacy Survey morphology).

## Fixed before computing (written 03:45 UTC, before any contamination or bias number of this run)

**Data.** DESI DR1 BGS group VAC (`bgs-groups-dr1-v1.0.fits`): r < 19.5, 2.23 M galaxies, fastspecfit LOGMSTAR (h = 1, converted
with +2 log10(1/0.7)), rest-frame G_R, QUIESCENT, group id IGRP, IS_SAT, group halo mass M_HALO (taken as M_200m in M☉/h,
converted to h70 by /0.7; one value per group; each group has exactly one central). Flat ΛCDM, H0 = 70, Ωm = 0.3; physical
(proper) distances, as Brouwer.

**Lenses.** Spectroscopic redshift (Z_ASSIGNED_FLAG 0 or −1), 0.05 < z_spec < 0.3. In the emulated photo-z view, Brouwer's
cuts: 0.1 < z_ph < 0.5 and log M*_ph < 11.0, in his four stellar-mass bins [8.5, 10.3, 10.6, 10.8, 11.0] (log M☉ h70⁻²).

**Neighbours.** Every galaxy in the VAC (spectroscopic or assigned redshift; KiDS sees all of them photometrically). The truth
view uses the catalogue redshift. Variant: spectroscopic neighbours only (R-145's choice).

**Emulated photo-z (KiDS-bright error model).** Primary: Bilicki et al. 2021 (arXiv:2101.06010) published fit to KiDS-bright
minus GAMA, a generalized Lorentzian N(Δz) ∝ (1 + Δz²/(2as²))^−a with a = 2.613, s = 0.0149, Δz = δz/(1+z), i.e. a Student t
with ν = 2a − 1 = 4.23 and scale 0.0166 (its wings are the outliers; std 0.023). Every galaxy gets an independent draw.
Stellar masses rescale by D_L(z_ph)²/D_L(z)² and get 0.12 dex Gaussian scatter (Brouwer's MICE emulation, App. A).
Variant: Brouwer's Gaussian σ_z = 0.02(1+z). Five seeds each.

**Isolation criteria (photo-z view).** No galaxy with M*_ph > 0.1 M*_ph,lens within a 3D physical distance R_isol (projected
separation at D_A(z_ph,lens), line-of-sight distance from the photo-z comoving distances / (1 + z_ph,lens)). Brouwer: R_isol =
3 Mpc/h70. Mistele: 4 Mpc/h70 (with h70 = 1 here).

**Contamination statistic (truth view).** For each photo-z-isolated lens: does a real neighbour with M* > 0.1 M*_lens (catalogue
masses) lie within projected R and |Δv| < 500 km/s? Reported at R = 0.3, 0.5, 1, 2, 3, 4 Mpc, per type, per mass bin, per
criterion, as the fraction of photo-z-isolated lenses.

**ESD bias model.** The contaminating halos of a lens are the distinct groups of its contaminating neighbours (within R_isol,
500 km/s, > 10 %). Each adds to the lens's own profile:
- (a) another group's halo: an NFW of mass M_HALO (Duffy et al. 2008 c(M, z), M_200m), centred on that group's central, at its
  projected offset d, azimuthally averaged about the lens: ΔΣ_off(R | M, d).
- (b) the lens is a satellite and its own group holds a contaminant: the host halo, NFW(M_HALO) at the offset to the central.
- (c) the lens is a central and its own group holds a contaminant: its halo is the group halo instead of an isolated one,
  ΔΣ_NFW(M_HALO) − ΔΣ_NFW(M_iso(M*)), centred, where M_iso(M*) is the catalogue's median M_HALO of N_SAT = 0 centrals at the
  same stellar mass (0.1 dex bins); floored at 0.
The lens's own profile is the same for both types (the null being tested is "both types share one RAR"): ΔΣ_own(g_bar) =
g_M16(g_bar) / 4G, McGaugh et al. 2016, g† = 1.2 × 10⁻¹⁰ m s⁻², the relation the isolated KiDS LTGs follow. g_obs = 4GΔΣ (Brouwer's
SIS conversion, his Eq. 9) in his 15 log bins 10⁻¹⁵ < g_bar < 5 × 10⁻¹² m s⁻². Each lens enters bin k at R_i,k = (G M_gal,i /
g_k)^½, M_gal = M*_ph (1 + f_cold), log f_cold = −0.69 log M* + 6.63 (Brouwer Eq. 12). Lens weight per bin ∝ (R_i,k / D_A)² Σ_crit⁻²
(sources per annulus times Brouwer's W_ls), Σ_crit with a KiDS-like source n(z) ∝ z² exp[−(z/0.5)^1.5] and z_s > z_l + 0.2.
Bias per type and bin: b_k = Σ_i w ΔΣ_cont,i(R_i,k) / Σ_i w ΔΣ_own(g_k); predicted log g_obs bias = log10(1 + b_k).

**Statistic and claim.** Types mass-matched as Brouwer did (both reweighted to the overlap of their M*_ph histograms in 0.05 dex
bins; his random removal in expectation). Predicted spurious offset δ = log10( mean_k (1 + b_E,k)/(1 + b_L,k) ), exactly his
log10⟨g_obs,E/g_obs,L⟩, over (i) bins with centre g_bar > 10⁻¹³ (vs 0.14 Sérsic / 0.19 colour) and (ii) all 15 bins (vs 0.17 /
0.27). The claim is the table: per criterion and split, δ and δ / published with interval, plus b_k against radius and g_bar.

**Type splits.** Brouwer's u−r > 2.5 (observer frame) is not in DESI; the DESI equivalents: (1) colour, rest-frame G_R > 0.75
(the valley of the G_R histogram at 9.8 < log M* < 10.5, 0.05 < z < 0.3, read before any lensing number); (2) QUIESCENT
(fastspecfit); (3) morphology, Legacy Survey MORPHTYPE from the LSS catalogue matched by TARGETID: DEV → bulge (n = 4), EXP and REX
→ disc (n = 1), SER by its index if the catalogue carries it (n > 2 bulge), else dropped; PSF dropped. Colour is compared with
Brouwer's u−r numbers (0.19 / 0.27), morphology with his Sérsic numbers (0.14 / 0.17); QUIESCENT is shown against both.

**Intervals.** Per seed, a fresh photo-z draw and a fresh halo-mass draw (each group's log M_HALO + N(0, 0.2), plus one global
shift N(0, 0.15) dex for the mass definition and h convention); 100 bootstrap resamples over lenses per seed. Interval: the 16-84
and 2.5-97.5 percentiles of the pooled seed × bootstrap distribution.

**Checks.** (i) True isolation: the spectroscopically isolated sample has no contaminant, so zero bias by construction (stated,
and run through the code). (ii) Zero scatter (z_ph = z_spec, no mass scatter) through the same pipeline: what is left is the 3D
criterion's own redshift-space leak. (iii) Gaussian 0.02(1+z) variant. (iv) Depth: DESI BGS stops at r = 19.5, KiDS-bright at
r = 20, so BGS misses neighbours between 19.5 and 20 that KiDS could have used to reject a lens and that contribute mass. Stated
systematic, bounded by restricting lenses to r < 17.0, for which every > 10 % neighbour (≤ 2.5 mag fainter) is inside BGS. (v)
Spectroscopic-only neighbours.

**What each outcome would mean (written before computing).**
- Fraction ≳ 0.7 at g_bar > 10⁻¹³: type-dependent leakage alone produces most of the published small-radius offset; Brouwer's
  equal-contamination assumption would carry the split.
- 0.2-0.7: leakage explains part; the rest must come from the baryonic masses (Mistele's ETG M/L and hot gas) or be real.
- ≲ 0.2 with the interval reaching 0: leakage does not explain the small-radius offset; the split there is set by masses or is
  real. The radius dependence then says where leakage does matter (large R, small g_bar), which is where Mistele's ETG R_isol
  dependence sits. A large prediction over the full range but a small one at g_bar > 10⁻¹³ is the "partial" reading for the
  0.27 / 0.17 dex numbers and the "null" one for 0.19 / 0.14.
- Prior expectation (geometry, before computing): an offset halo at projected distance d adds little ΔΣ at R < d, so the bias
  should grow with R and be small inside 0.3 Mpc unless the lens is a satellite of a massive host.

## Run (04:00-04:20 UTC; `r147_leak.py`, `r147_table.py`, `r147_morph.py`)

- 2,229,727 galaxies with a redshift; 778,963 lens candidates per seed (median z_spec 0.20, median log M*_ph 10.61).
  Photo-z isolation passes 46 % of them at 3 Mpc and 22 % at 4 Mpc (about 360,000 and 168,000 lenses). Each run takes 65 s.
- Morphology: the LSS catalogue has MORPHTYPE but no Sérsic index, and 65 % of its rows are SER. As fixed, SER is dropped, so
  the morphology split uses DEV against EXP/REX: about 40 % of the lenses, 42,000 bulge and 110,000 disc lenses at 3 Mpc.
- Check (i) holds: no lens that is isolated spectroscopically carries a contaminating halo. The code asserts this for every
  seed, split and criterion.
- Mass matching keeps 147,000 ETG and 210,000 LTG lenses per seed (colour split, 3 Mpc).

### Contamination: fraction of photo-z-isolated lenses with a real > 10 % neighbour within R and 500 km/s

Primary error model, colour split, mean over 5 seeds. The seed spread is 0.005-0.026. Full tables for every split and
criterion are in `data/r147/table_lorentz.md`, regenerated by `r147_table.py`.

| Criterion | log M*_ph | ETG 0.3 / 0.5 / 1 / 3 Mpc | LTG 0.3 / 0.5 / 1 / 3 Mpc |
|---|---|---|---|
| 3 Mpc | 8.5-10.3 | 0.39 / 0.54 / 0.73 / 0.95 | 0.22 / 0.35 / 0.57 / 0.91 |
| 3 Mpc | 10.3-10.6 | 0.32 / 0.46 / 0.66 / 0.93 | 0.19 / 0.30 / 0.52 / 0.88 |
| 3 Mpc | 10.6-10.8 | 0.29 / 0.42 / 0.63 / 0.92 | 0.18 / 0.29 / 0.50 / 0.87 |
| 3 Mpc | 10.8-11.0 | 0.28 / 0.40 / 0.61 / 0.91 | 0.18 / 0.28 / 0.50 / 0.87 |
| 3 Mpc | all | 0.30 / 0.43 / 0.63 / 0.92 | 0.20 / 0.31 / 0.53 / 0.89 |
| 4 Mpc | all | 0.27 / 0.39 / 0.59 / 0.90 | 0.18 / 0.28 / 0.49 / 0.86 |

QUIESCENT gives the same picture, 0.31 / 0.44 against 0.20 / 0.32 at 0.3 / 0.5 Mpc. Morphology narrows the gap, 0.27 / 0.39
against 0.20 / 0.32. What drives the gap: 28 % of photo-z-isolated red lenses are group satellites with a > 10 % neighbour
inside 3 Mpc, against 12 % of blue ones. The 4 % and 3 % that are centrals with such a satellite add little.

### Predicted log g_obs bias against g_bar and radius

Colour split, 3 Mpc. R is the bin radius for M_gal = 10^10.6.

| g_bar (m/s²) | R (Mpc) | ETG | LTG | ETG − LTG |
|---|---|---|---|---|
| 3.8e-12 | 0.04 | 0.023 | 0.013 | 0.009 |
| 1.2e-12 | 0.07 | 0.047 | 0.024 | 0.023 |
| 3.9e-13 | 0.12 | 0.088 | 0.040 | 0.048 |
| 1.25e-13 | 0.21 | 0.153 | 0.065 | 0.088 |
| 4.0e-14 | 0.37 | 0.244 | 0.103 | 0.141 |
| 1.3e-14 | 0.66 | 0.350 | 0.159 | 0.191 |
| 4.1e-15 | 1.16 | 0.445 | 0.227 | 0.218 |
| 1.3e-15 | 2.04 | 0.560 | 0.370 | 0.190 |

- Component (b), the host halo of satellite lenses, carries 70-90 % of the ETG bias at every radius (87 % at 1.25e-13). Component (a),
  other groups' halos, is negative just inside its offset and positive beyond, so it nets to ≤ 0.02 inside 1 Mpc.
- At 4 Mpc the ETG bias is 0.005-0.02 dex lower at 10⁻¹⁴ < g_bar < 10⁻¹³, and up to 0.09 dex lower at the lowest g_bar.

### The claim: predicted spurious offset log10⟨g_obs,E/g_obs,L⟩ against the published one

Primary error model (Bilicki Lorentzian), 5 seeds × 100 bootstrap resamples, with the halo-mass draws inside.

| Split (published) | Criterion | δ at g_bar > 10⁻¹³ (68 %; 95 %) | fraction of 0.19 / 0.14 | δ over all 15 bins (68 %; 95 %) | fraction of 0.27 / 0.17 |
|---|---|---|---|---|---|
| colour G_R (u−r: 0.19 / 0.27) | 3 Mpc | 0.041 (0.038-0.046; 0.036-0.048) | 0.22 (0.20-0.24) | 0.121 (0.111-0.136; 0.109-0.140) | 0.45 (0.41-0.50) |
| colour G_R | 4 Mpc | 0.040 (0.036-0.045; 0.034-0.048) | 0.21 (0.19-0.24) | 0.112 (0.104-0.126; 0.100-0.131) | 0.41 (0.39-0.47) |
| QUIESCENT (vs u−r) | 3 Mpc | 0.048 (0.044-0.054; 0.043-0.057) | 0.26 (0.23-0.29) | 0.130 (0.120-0.143; 0.117-0.148) | 0.48 (0.44-0.53) |
| QUIESCENT | 4 Mpc | 0.047 (0.042-0.051; 0.039-0.056) | 0.25 (0.22-0.27) | 0.122 (0.111-0.135; 0.107-0.140) | 0.45 (0.41-0.50) |
| morphology DEV vs EXP/REX (Sérsic: 0.14 / 0.17) | 3 Mpc | 0.031 (0.029-0.035; 0.027-0.039) | 0.22 (0.20-0.25) | 0.086 (0.080-0.091; 0.076-0.096) | 0.51 (0.47-0.54) |
| morphology | 4 Mpc | 0.032 (0.028-0.037; 0.024-0.040) | 0.23 (0.20-0.27) | 0.083 (0.077-0.091; 0.072-0.096) | 0.49 (0.45-0.53) |

The published values are Brouwer's, all at 3 Mpc. Mistele published no offset at 4 Mpc with Brouwer's masses, so the 4 Mpc
rows compare against the same numbers. They show the prediction barely depends on the criterion: the photo-z line-of-sight
error, σ_z c/H ≈ 100 Mpc, dwarfs the 1 Mpc difference between the two radii. Per Brouwer mass bin, δ at g_bar > 10⁻¹³ is
0.038-0.045 (colour, 3 Mpc), so no single bin carries it.

### Checks and systematics (colour split, 3 Mpc; δ at g_bar > 10⁻¹³ / all bins, fraction of 0.19 / 0.27)

| Variant | δ | Fraction |
|---|---|---|
| Primary (Lorentzian, all neighbours) | 0.041 / 0.121 | 0.22 / 0.45 |
| Gaussian 0.02(1+z), Brouwer's MICE emulation | 0.043 / 0.128 | 0.23 / 0.47 |
| Spectroscopic neighbours only (2 seeds) | 0.040 / 0.119 | 0.21 / 0.44 |
| Zero scatter (z_ph = z_spec): the 3D criterion's own redshift-space leak | 0.012 / 0.044 | 0.06 / 0.16 |
| Depth: lenses r < 17.0 (~1,300 lenses; every > 10 % neighbour inside BGS) | 0.059 (68 % 0.041-0.073) / 0.099 | 0.31 (0.22-0.38) / 0.37 |
| Depth, QUIESCENT, r < 17.0 | 0.066 (0.050-0.084) / 0.106 | 0.35 (0.26-0.44) / 0.39 |

**Amplitude check against Brouwer's App. A.** Over all photo-z-isolated lenses with no type split, the contamination adds
1.05 times the own profile at g_bar < 7.6 × 10⁻¹⁴ (his R > 0.3 Mpc/h70) and 0.15 above it. His MICE test found the
photo-z-isolated ESD about 0.30 above the truly isolated ESD. The denominators differ:
- Ours is the McGaugh one-galaxy profile.
- His is the truly isolated total ESD, which also holds everything beyond R_isol (the 2-halo term). Around 10^10.6 galaxies
  that term dominates at R ≳ 1 Mpc.

So the two agree only if the environment beyond R_isol is about 2.5 times the own profile there. That is plausible, but this
run cannot check it. Two consequences follow:
- The full-range fractions (0.41-0.51) are upper estimates in that sense.
- At R < 0.3 Mpc the 1-halo term dominates and the comparison does not apply.

### Reading, against the outcomes written in advance

- **g_bar > 10⁻¹³ (R < 0.3 Mpc).** The prediction is 0.21-0.26 of the published offset for every split and criterion. The
  95 % range is 0.19-0.30, and 0.31-0.35 for the bright, depth-complete lenses. That sits at the low edge of "partial": leakage
  accounts for about a fifth to a third of the 0.14 / 0.19 dex. The other 0.10-0.15 dex has to come from the masses (Mistele's
  ETG M/L and hot gas) or be a real difference.
- **Full range.** The prediction is 0.41-0.51 of 0.17 / 0.27 dex, the upper estimate above: partial.
- **Radius.** The predicted ETG−LTG offset grows from 0.01 dex at R = 0.04 Mpc to 0.09 at 0.2 Mpc and 0.19-0.22 at 0.7-2 Mpc.
  It is set by red lenses being twice as often hidden satellites.
- **Mistele's criterion.** At 4 Mpc the leakage barely moves (≤ 0.02 dex at 10⁻¹⁴ < g_bar < 10⁻¹³). Leakage of > 10 %
  neighbours therefore does not predict the change in the ETG RAR between 3 and 4 Mpc that Mistele saw. Within this model,
  the agreement of the two types at 4 Mpc has to come from their masses, from statistics, or from something the model leaves
  out (neighbours under 10 %, or the environment beyond R_isol).

## Finding check (gate), 2026-11-02, context-free Opus, the section alone (`briefs/r147_finding_section.md`)

VERDICT: novel but not sound: the spurious-offset fractions depend on an uncalibrated NFW contamination amplitude that exceeds Brouwer's MICE benchmark by a factor of about 3.5 (1.05 against 0.30 of the profile), divided by a reference profile that is the same for both types, so the claimed share of the published offset is not supported.

- **Novelty held.** The referee found no type-resolved photo-z isolation contamination in print. Closest: arXiv:2106.11677,
  2310.15248, 1509.06762, 2104.05355, 2406.09685, 2609.14085.
- **The miss is a physics point,** so it is parked as the brief says, with no revision pass. Two other objections came with it:
  - The truth window (500 km/s, about 7 Mpc along the line of sight) is wider than the 3 Mpc 3D criterion. That is the 0.012
    dex zero-scatter residual.
  - BGS fibre incompleteness in dense regions.
- **What stands:** the contamination and satellite fractions by type.
- **What is open:** the share of the published offset, bracketed at 0.06-0.26 (g_bar > 10⁻¹³) between the MICE-scaled and the
  computed amplitude.

# R-149 (cosmo-deep-16c, 2026-11-02, cloud worker cw-156): the contamination amplitude calibrated on a mock

Brief: the R-147 gate failed on the amplitude (an NFW contamination added to a common McGaugh profile, 1.05 of it at
R > 0.3 Mpc against Brouwer's MICE 0.30). This revision computes the amplitude on a public halo-based mock, against truly
isolated lenses of the same type. Scripts: `r149_mock.py` (mock isolation, contamination, ESD from the mock's halos),
`r149_maglim.py` (the BGS-calibrated magnitude model for the mock), `r149_desi_ci.py` (the DESI leakage interval),
`r149_table.py`. Data under `data/r149/` (gitignored).

## R-149 fixed before computing (written 04:43 UTC, before any mock or interval number)

**Mock: UniverseMachine DR1 on Bolshoi-Planck** (Behroozi et al. 2019, arXiv:1806.07893), snapshot a = 0.827654 (z = 0.208,
the median BGS lens redshift), `halos.as.arizona.edu/UniverseMachine/DR1/SFR/sfr_catalog_0.827654.bin` (1.7 GB, fetched in
2 min 20 s; the BGS group VAC again in 20 s). Box 250 Mpc/h, Planck cosmology (h = 0.678, Ωm = 0.307), 13.9 M halos and
subhalos (Rockstar), each with position, peculiar velocity, BN98 virial mass and radius (current mass for subhalos), central
or satellite (UPID), and a galaxy with stellar mass and SFR.
- Why this mock: it is public without registration and carries what the test needs at once: halos (so the ESD comes from
  the mock's own mass distribution, satellites in their true hosts), stellar masses down to ~10^8.5, and a star-forming /
  quenched split whose dependence on halo and central/satellite status is fit to SDSS quenched fractions, including those of
  isolated centrals and the clustering of quenched and star-forming galaxies. That dependence is what sets the
  type-dependent leakage.
- Not used: MICE (Brouwer's mock) needs a CosmoHub account; the AbacusSummit halos and light cones are public
  (data.desi.lbl.gov/public/cosmosim) but carry no galaxies, so colours would have to come from an HOD, which is the thing
  under test; the DESI AbacusSummit BGS mocks with colours are not in the public tree.
- Not a light cone: one snapshot, observed at virtual redshifts z_v (below), with no evolution across them.

**Mock galaxies.** True stellar mass = UM `obs_sm` (with UM's observational scatter, the analogue of fastspecfit masses).
Type: quenched if obs_sfr / obs_sm < 10⁻¹¹ yr⁻¹ (the analogue of DESI QUIESCENT and of the colour split; compared with both
DESI rows). Units: physical Mpc at the snapshot (comoving Mpc/h × a / h), halo masses M☉ (÷ h). Line of sight: each box
axis in turn (3 lines of sight) × 2 photo-z seeds = 6 realisations.

**Magnitudes (selection).** Each mock galaxy gets an apparent r magnitude from the BGS itself: at each z_v, BGS galaxies in
|z − z_v| < 0.02 give m_r − DM(z) = a_T + b_T (log M* − 10.5) + N(0, σ_T) per type T (DESI QUIESCENT), fit by maximum
likelihood with the truncation at r = 19.5 (`r149_maglim.py`). Selections: BGS-like r < 19.5 (the comparison with DESI) and
KiDS-bright-like r < 20 (the prediction for Brouwer). Lenses and neighbours both need the limit.

**Virtual redshifts.** BGS-like: z_v = 0.10, 0.20, 0.28, weighted by the fraction of R-147's photo-z-isolated lenses in
0.05-0.15, 0.15-0.25, 0.25-0.30. KiDS-like: z_v = 0.15, 0.25, 0.35 weighted 1/4, 1/2, 1/4 (KiDS-bright lenses, mean z 0.25).

**Photo-z view.** Same as R-147: Bilicki's Lorentzian (a = 2.613, s = 0.0149) in Δz/(1+z_v), per galaxy; the line-of-sight
position moves by c δz / H(z_v) (to physical), periodic in the box; z_ph = z_v + (1+z_v)(Δz + v_los/c). Photo masses:
log M* + 2 log10(D_L(z_ph)/D_L(z_obs)) + N(0, 0.12). Lens cuts in the photo view as Brouwer: 0.1 < z_ph < 0.5,
8.5 < log M*_ph < 11. Photo-z isolation: no neighbour with M*_ph > 0.1 M*_ph,lens within 3 Mpc in 3D (photo positions);
the 4 Mpc criterion as a variant.

**Truly isolated (the reference).** The same criterion with true masses in true redshift space (real position + v_los/H).
This is Brouwer's App. A reference. Both samples are binned in g_bar with photo masses (what KiDS would do), so the ratio
isolates the isolation effect; mass scatter is the same for both types.

**ESD from the mock's halos.** For each lens, ΔΣ(R_k) = Σ over distinct halos (UPID = −1) with M_vir > 10^11 M☉ whose
centre lies within projected 8 Mpc and ±30 Mpc along the line of sight in real space, each as an NFW truncated at its R_vir
(catalogue mass and radius, Dutton & Macciò 2014 c_vir(M, z)), azimuthally averaged about the lens at its projected offset;
plus, when the lens is a satellite, its own subhalo (current mass and radius, truncated NFW, centred); plus the lens's baryons
as a point mass M_gal = M*(1 + f_cold) (Brouwer Eq. 12). The lens's own host halo enters with offset 0 if it is the central,
whatever its mass. No amplitude is assumed: masses and positions are the mock's.

**Amplitude statistic, per type.** A_T,k = ⟨ΔΣ⟩(photo-z-isolated, T) / ⟨ΔΣ⟩(truly isolated, T), in Brouwer's 15 g_bar bins,
with R-147's lens weights (R_k² / D_A² Σ_crit⁻² at z_v); each sample mass-matched between types on M*_ph as Brouwer did.
Brouwer's benchmark (amplitude check): A over all lenses (no type split) in fixed radial bins, averaged over 0.3-3 Mpc,
against his MICE ~1.30.

**Predicted spurious offset.** δ = log10 mean_k (A_E,k / A_L,k), over g_bar > 10⁻¹³ (vs 0.14 / 0.19) and all 15 bins (vs
0.17 / 0.27); also per bin against radius. Equivalent to Brouwer's log⟨g_E/g_L⟩ measured with photo-z isolation minus the same
with true isolation, so a real type difference in the truly isolated profiles drops out (the referee's second point).

**Check against the data, and the rescaling.** The mock's contamination statistic exactly as R-147 (a real neighbour with
M* > 0.1 M*_lens within projected R and |Δv| < 500 km/s, fraction of photo-z-isolated lenses, mass-matched, BGS-like
selection) and its hidden-satellite fraction, against DESI's (colour 0.30 / 0.43 vs 0.20 / 0.31 at 0.3 / 0.5 Mpc; QUIESCENT
0.31 / 0.44 vs 0.20 / 0.32; satellites 28 % vs 12 %). If they differ, the rescaled row multiplies each type's excess
(A_T,k − 1) by f_T,DESI / f_T,mock at 0.5 Mpc (QUIESCENT row, the same type definition); variants: 0.3 Mpc, the satellite
fraction.

**Intervals.** Mock: 16 jackknife regions (4 × 4 in the plane of the sky) per realisation, the jackknife variance plus the
scatter over the 6 realisations; quoted as 68 % and 95 % (± 1 and 2 σ). DESI leakage (the headline): R-147's pipeline rerun
for 5 seeds, the ETG − LTG difference and ratio of the fractions with a sky jackknife (≈ 50 regions of equal lens count in
RA, Dec) plus the seed spread.

**Variants.** 4 Mpc criterion; untruncated NFW (R-147's profile); halos down to 10^10.5 M☉; r < 19.5 vs r < 20; per z_v.

**What each outcome would mean (written before computing).**
- If the mock's all-lens A at 0.3-3 Mpc is near Brouwer's 1.3, the mock reproduces his benchmark and its type split stands
  in for his. If far above or below, say so; the type ratio is still the quantity, but the benchmark is not reproduced.
- δ / published ≳ 0.7 at g_bar > 10⁻¹³: leakage carries the small-radius split; 0.2-0.7: part; ≲ 0.2 with the interval
  reaching 0: it does not. Prior (from R-147 scaled to MICE): 0.06-0.26.
- If the mock's fractions match DESI's within their intervals, the raw and rescaled rows agree; if not, the rescaled row is
  the claim and the mismatch is reported.

## R-149 run (04:50-05:35 UTC; `r149_maglim.py`, `r149_mock.py`, `r149_desi_ci.py`, `r149_table.py`)

- Fetch: the UM snapshot (1.7 GB) in 2 min 20 s, the group VAC (290 MB) in 20 s. R-147's pipeline rerun for 5 seeds (65 s
  each). Mock: 6 virtual redshifts × 3 lines of sight × 2 seeds = 36 runs of ~2.5 min, plus 3 variant runs. Full tables:
  `data/r149/table_r3.md`, `table_r4.md` (regenerated by `r149_table.py`).
- Magnitude model (BGS, truncated fit at z = 0.20): m_r − DM = −20.42 − 1.99 (log M* − 10.5), scatter 0.35 (quiescent);
  −20.72 − 1.51 (log M* − 10.5), scatter 0.42 (star-forming). The other z_v differ by ≤ 0.1 mag.
- Mock at z_v = 0.20, r < 19.5: 131,000 galaxies, 108,000 lens candidates, 37 % quenched; photo-z isolation passes 53 % (3 Mpc)
  and 26 % (4 Mpc), true isolation 18 % and 9 %. DESI: 46 % and 22 % pass photo-z isolation.
- Implementation choices (none changes a definition above): the ESD stacks of the four samples (photo-z / truly isolated ×
  type) are all reweighted to one common M*_ph shape (the bin-wise minimum of their normalised histograms), so A_T compares
  like with like; the benchmark reweights the truly isolated to the photo-z-isolated shape. The DESI fractions here are
  mass-matched between types (R-147's table was not), hence 0.33 / 0.46 rather than 0.30 / 0.43 for red.
- Two variant runs were killed for memory when run side by side. They were rerun alone at z_v = 0.20, r < 19.5.

### Row 1 (the headline): the measured type-dependent leakage, DESI DR1 BGS, 3 Mpc, with interval

5 seeds, sky jackknife (50 regions) plus seed scatter; types mass-matched.

| split | statistic | ETG | LTG | ETG − LTG (68 %; 95 %) | ETG / LTG (68 %; 95 %) |
|---|---|---|---|---|---|
| colour G_R | neighbour within 0.3 Mpc | 0.319 ± 0.002 | 0.192 ± 0.001 | 0.128 (0.126-0.130; 0.124-0.132) | 1.67 (1.65-1.68; 1.64-1.69) |
| colour | within 0.5 Mpc | 0.453 ± 0.002 | 0.303 ± 0.002 | 0.149 (0.147-0.151; 0.145-0.153) | 1.49 (1.48-1.50; 1.48-1.51) |
| colour | within 1 Mpc | 0.651 ± 0.002 | 0.517 ± 0.002 | 0.134 (0.132-0.136; 0.129-0.139) | 1.26 (1.25-1.26; 1.25-1.27) |
| colour | hidden satellite | 0.290 ± 0.002 | 0.113 ± 0.001 | 0.177 (0.176-0.179; 0.174-0.181) | 2.57 (2.55-2.60; 2.53-2.62) |
| QUIESCENT | within 0.3 Mpc | 0.327 ± 0.002 | 0.200 ± 0.001 | 0.127 (0.124-0.129; 0.122-0.131) | 1.63 (1.62-1.65; 1.61-1.66) |
| QUIESCENT | within 0.5 Mpc | 0.461 ± 0.002 | 0.314 ± 0.002 | 0.147 (0.145-0.149; 0.143-0.152) | 1.47 (1.46-1.48; 1.45-1.49) |
| QUIESCENT | hidden satellite | 0.310 ± 0.002 | 0.118 ± 0.001 | 0.192 (0.190-0.194; 0.188-0.196) | 2.62 (2.60-2.65; 2.57-2.68) |

The interval is statistical (photo-z draws and sky). The emulation's systematics (BGS depth, fibre incompleteness, the error
model) are larger and are not in it; R-147's variants moved the fractions by ≤ 0.03.

### Check against the data: the mock's leakage (QUIESCENT-like split), 3 Mpc, 6 realisations, jackknife + scatter

| statistic | mock Q | mock SF | mock Q − SF | DESI Q | DESI SF | DESI Q − SF |
|---|---|---|---|---|---|---|
| within 0.3 Mpc, BGS-like r < 19.5 | 0.369 ± 0.004 | 0.220 ± 0.004 | 0.149 ± 0.006 | 0.327 | 0.200 | 0.127 |
| within 0.5 Mpc, BGS-like | 0.484 ± 0.005 | 0.322 ± 0.005 | 0.162 ± 0.007 | 0.461 | 0.314 | 0.147 |
| hidden satellite, BGS-like | 0.295 ± 0.005 | 0.140 ± 0.003 | 0.155 ± 0.006 | 0.310 | 0.118 | 0.192 |
| within 0.3 / 0.5 Mpc, KiDS-like r < 20 | 0.358 / 0.472 | 0.211 / 0.311 | 0.146 / 0.161 | | | |
| hidden satellite, KiDS-like | 0.289 | 0.135 | 0.153 | | | |

The mock reproduces DESI's type-dependent leakage: within 0.02-0.04 of each fraction and 0.02 of the 0.5 Mpc difference.
Formally the 0.3 / 0.5 Mpc fractions differ by more than the (statistical) intervals, so the rescaled rows are computed as
fixed: factors f_DESI / f_mock (Q, SF) of 0.95, 0.98 at 0.5 Mpc; 0.89, 0.91 at 0.3 Mpc; 1.05, 0.84 for hidden satellites.

### Amplitude: Brouwer's App. A benchmark (all lenses, ESD of photo-z-isolated / truly isolated, 3 Mpc)

| R (Mpc) | 0.03 | 0.05 | 0.08 | 0.12 | 0.19 | 0.30 | 0.48 | 0.75 | 1.19 | 1.89 | 3.00 | mean 0.3-3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BGS-like | 1.01 | 1.02 | 1.05 | 1.10 | 1.15 | 1.22 | 1.37 | 1.55 | 1.80 | 2.16 | 3.01 | 1.85 ± 0.11 |
| KiDS-like | 1.01 | 1.03 | 1.05 | 1.09 | 1.14 | 1.21 | 1.35 | 1.52 | 1.76 | 2.10 | 2.89 | 1.81 ± 0.07 |

The mock does not reproduce Brouwer's MICE ~1.30: its excess at 0.3-3 Mpc is 0.8, about 2.7 times his. It reaches 1.3 only
at 0.4 Mpc. Brouwer's number may refer to a narrower radial range or to MICE's own galaxy model; this run cannot tell which.
It matters little for the claim: the claim is the ratio of two types' amplitudes, and those agree at every radius, so scaling
both excesses by a common factor leaves δ near zero.

### Row 2: the predicted spurious offset, δ = log10 mean_k (A_Q,k / A_SF,k), per-type references (3 Mpc)

Mean; 68 %; 95 % (jackknife + realisation scatter). Fraction of 0.19 (u−r) / 0.14 (Sérsic) at g_bar > 10⁻¹³; of 0.27 / 0.17
over all bins.

| selection | amplitude | δ at g_bar > 10⁻¹³ | fraction | δ, all 15 bins | fraction |
|---|---|---|---|---|---|
| KiDS-like (the prediction for Brouwer) | mock | −0.001 (−0.005 to 0.002; −0.008 to 0.006) | −0.01 (−0.03 to 0.01) / −0.01 | −0.008 (−0.020 to 0.003; −0.031 to 0.014) | −0.03 (−0.07 to 0.01) / −0.05 |
| KiDS-like | rescaled to DESI, 0.5 Mpc | −0.002 (−0.005 to 0.001; −0.008 to 0.005) | −0.01 (−0.03 to 0.01) / −0.01 | −0.010 (−0.021 to 0.001; −0.032 to 0.011) | −0.04 (−0.08 to 0.00) / −0.06 |
| KiDS-like | rescaled, 0.3 Mpc | −0.002 (−0.005 to 0.001; −0.008 to 0.005) | −0.01 / −0.01 | −0.010 (−0.020 to 0.000) | −0.04 / −0.06 |
| KiDS-like | rescaled, hidden satellites | 0.004 (0.001 to 0.007; −0.003 to 0.011) | 0.02 (0.00 to 0.04) / 0.03 | 0.010 (−0.001 to 0.021; −0.013 to 0.032) | 0.04 (−0.01 to 0.08) / 0.06 |
| BGS-like | mock | 0.001 (−0.004 to 0.006; −0.008 to 0.011) | 0.01 (−0.02 to 0.03) / 0.01 | −0.004 (−0.015 to 0.007; −0.026 to 0.018) | −0.01 (−0.05 to 0.03) / −0.02 |
| BGS-like | rescaled, 0.5 Mpc | 0.001 (−0.004 to 0.005) | 0.00 / 0.00 | −0.006 (−0.016 to 0.005) | −0.02 / −0.03 |
| BGS-like | rescaled, hidden satellites | 0.007 (0.002 to 0.011; −0.003 to 0.016) | 0.03 (0.01 to 0.06) / 0.05 | 0.015 (0.004 to 0.026; −0.007 to 0.038) | 0.06 (0.02 to 0.10) / 0.09 |
| KiDS-like, 4 Mpc | mock | 0.006 (−0.002 to 0.013; −0.009 to 0.020) | 0.03 (−0.01 to 0.07) / 0.04 | 0.000 (−0.014 to 0.015) | 0.00 / 0.00 |
| KiDS-like, 4 Mpc | rescaled, hidden satellites | 0.013 (0.005 to 0.021; −0.003 to 0.029) | 0.07 (0.03 to 0.11) / 0.09 | 0.018 (0.002 to 0.033) | 0.07 / 0.11 |

**Radius dependence** (KiDS-like, 3 Mpc; R at M* = 10^10.6). log10(A_Q/A_SF) is −0.005 to +0.002 (± 0.002-0.008) at
0.04-0.23 Mpc, −0.002 to +0.003 (± 0.01-0.02) at 0.3-0.7 Mpc, and turns negative beyond 1 Mpc: −0.005, −0.019, −0.028, −0.076
(± 0.03-0.06) at 0.96, 1.27, 1.69, 2.25 Mpc. Each type's own amplitude grows from 1.01-1.02 at 0.04 Mpc to 1.12 at 0.23 Mpc,
1.30 at 0.54 Mpc and 2.1-2.5 at 2.2 Mpc.

**Variants** (z_v = 0.20, r < 19.5, one realisation; δ at g_bar > 10⁻¹³ / all bins): baseline 0.003 / −0.011; untruncated
NFW 0.003 / −0.005; halos down to 10^10.5 M☉ 0.003 / −0.007.

**Diagnostic (the referee's two points separated).** The mock's own excess ESD put over a common reference (the truly isolated
star-forming profile, R-147's construction): δ = 0.010 (0.006-0.015) at g_bar > 10⁻¹³ (0.06 / 0.07 of the published) and
0.082 (0.071-0.093) over all bins (KiDS-like, 3 Mpc); 0.023 / 0.090 at 4 Mpc. With the mock's amplitude, R-147's
construction gives a quarter of R-147's small-radius number, matching its MICE-scaled estimate (0.01). The rest goes once each
type has its own reference.

**Not pre-registered, an observation.** The truly isolated mock lenses themselves differ by type: log10(ESD_Q/ESD_SF) at fixed
M*_ph is 0.09-0.21 at g_bar > 10⁻¹³ and 0.25-0.42 below. Brouwer's statistic on them gives 0.14 (g_bar > 10⁻¹³) and 0.26 (all
bins). This is the ΛCDM halo-mass gap at fixed stellar mass as UniverseMachine sets it. It is why the fractional leakage
cancels: red lenses gain more contaminating mass, but their own halos are heavier in about the same proportion.

### Reading, against the outcomes written in advance

- **Benchmark:** not reproduced (1.8 against 1.3); stated. The type ratio is still the quantity, and it is insensitive to a
  common amplitude.
- **δ / published at g_bar > 10⁻¹³:** −0.01 (68 %: −0.03 to 0.01) for the pre-registered KiDS-like prediction, ≤ 0.07 in
  every variant, with every interval reaching 0 or within 0.01 of it. That is the "≲ 0.2 with the interval reaching 0"
  outcome: **type-dependent leakage does not explain the small-radius split.** Over all bins it is −0.04 to +0.07.
- **Data check:** the mock's fractions match DESI's to 0.02-0.04, and the rescaled rows agree with the raw ones.
- R-147's 0.21-0.26 was set by the common reference, not by the amplitude alone.

## R-149 finding check (gate), 2026-11-02, context-free Opus, the section alone (`briefs/r149_finding_section.md`)

VERDICT: novel but not sound: the lensing calibration is relative to each type's own isolated halo in a mock that already builds in the early/late halo gap being tested, and under the null of equal halos the note's own common-reference check leaves roughly half of the full-range offset attributable to leakage, so "not a leakage artefact" does not follow.

- **Novelty held.** The referee found no measurement of isolation leakage by type against spectroscopic truth, and no mock
  calibration of its effect on the early/late ratio. Closest: arXiv:2106.11677 (App. A, whole sample; the caveat is in its
  Sec. 5.4, not 5.3 as written), 2310.15248 (qualitative), 2408.05013 (UniverseMachine photo-z mocks for HSC, no type split, no
  isolation), 1107.4093 and 1507.04301 (halo-model satellite fractions by colour).
- **The miss is a physics point, so it is parked** (the brief's last pass). The objection: per-type references make the test
  circular. The mock already has the halo gap at fixed M* (0.14 / 0.26 dex) whose reality is in question. Under the strict null
  of equal halos, the common-reference row is the right one, and it leaves 0.08-0.09 dex over all bins (about half of 0.17).
  - Other points raised: the mock's all-lens amplitude is 2.7 times MICE's; fibre incompleteness and group-finder assignments
    could bias the DESI gap low; the photo-z model ignores colour; there is no Sérsic counterpart; the hidden-satellite rescaling
    changes the sign; and the intervals are statistical only.
- **What stands:**
  - Row 1, the measured type-dependent leakage: red 0.32 / 0.45 against blue 0.19 / 0.30 within 0.3 / 0.5 Mpc, hidden satellites
    0.29 against 0.11 (statistical intervals).
  - The mock reproduces it to 0.02-0.04.
  - At g_bar > 10⁻¹³ the leakage effect is small under either reference: ≤ 0.09 of the published offset with per-type references,
    and 0.06-0.17 with a common reference.
- **What is open:** the full-range share. It is 0 with per-type references and 0.3-0.6 under the equal-halo null. Which applies
  depends on whether isolated early types really have heavier halos at fixed M*, which is the question itself.
