# Results (R-146, cosmo-deep-17, 2026-11-02)

CMB lensing of z ~ 5 Lyman-break galaxies (HSC GOLDRUSH r-dropouts) with ACT DR6 and Planck PR3, and with the dropouts'
auto-spectrum the amplitude of matter fluctuations at z ~ 5 relative to Planck ΛCDM; the z ~ 4 g-dropouts (Miyatake et al. 2022,
arXiv:2103.15862) through the same pipeline beside it. Every number from `z5_cat.py` (catalogue), `z5_lens.py` (step 0 `feas`,
the measurement; usage in its docstring), theory from `qso_z4_theory.py` (R-139). Data under `data/` (gitignored); logs in
`data/z5/`. The trajectory check ("VERDICT: novel: no prior work answers this question", 2026-11-01) is in the brief; not rerun.

## What Miyatake et al. 2022 did (read in full before fixing anything)

g-dropouts only (no r-dropouts): 1,643,863 HSC PDR2 Wide g-dropouts with total i < 25.5 (1,473,106 after the Planck mask, 270
deg²), GOLDRUSH colour cuts, ⟨z⟩ = 3.8. Statistic: the stacked real-space κ profile around the LBGs on Planck PR3 MV (2′ Gaussian),
random-point mean field subtracted, jackknife covariance (262 regions): 5.1σ. Interlopers: f_cont ≈ 0.25 (Ono et al. 2018), their
κ profile modelled from HSC galaxy-galaxy lensing (significance against the interloper signal 3.5σ). CIB: a SMICA-noSZ temperature
stack shows no gradient. **Magnification bias: not modelled.** Halo mass 2.9 (+9.5 −2.5) × 10¹¹ h⁻¹ M_⊙ (NFW + Tinker bias,
θ = 6-20′); with the projected clustering (1.5-6.5′): σ8 = 0.46 (+0.25 −0.20), b = 6.5 (+5.2 −4.2), Ω_m unconstrained; 1.4σ below
Planck.

## Step 0: feasibility (02:50 UTC, before any κ map was crossed with any galaxy map)

Fetch: GOLDRUSH Zenodo 5512721 (3.3 GB zip; 5 min with 12 parallel range requests, a single stream ran at ~1 MB/s), ACT DR6
baseline/CIB-deprojected/curl alm and masks (2 min), Planck PR3 `COM_Lensing_4096_R3.00` (1 min). ACT simulations streamed from
02:44 (two streams after a stalled transfer; `act_dr6_fetch.sh` now aborts a transfer below 1 MB/s for 20 s).

**The r-dropout count, verified:** GOLDRUSH (Harikane et al. 2022, arXiv:2108.01090, Table 3) holds 139,359 r-dropouts in all
layers, **119,651 in the Wide layer** (`masked_riz_W_*`): not "a few hundred thousand". 1,692,323 Wide g-dropouts, 1,569,187 at
i < 25.5. The published catalogue is the S18A (≈ PDR2) one; no PDR3 dropout catalogue is public. GOLDRUSH ships no randoms: the
galaxy window is built from the catalogue's own holes (below).

| sample | in the galaxy window | × ACT DR6 overlap | × Planck PR3 overlap |
|---|---|---|---|
| r-dropouts (z ~ 5, all) | 115,950 on 352 deg² (effective) | 101,005 on 302 deg² (334/deg²) | 103,677 on 320 deg² |
| g-dropouts (z ~ 4, i < 25.5) | 1,550,740 | 1,354,372 on 302 deg² | 1,391,462 on 320 deg² |

ACT covers W02-W05 (HECTOMAP W06, Dec +43°, is outside); Planck covers all fields less its point-source/cluster holes.

Magnification and interlopers at these depths (from the published LF and Table 4, before any cross-spectrum): the sample-mean
intrinsic count slope s = 1.23 (r) and 0.97 (g), so 5s − 2 = 4.1 and 2.8: **magnification is not small**: at z ~ 5 the
(5s − 2) C^{κμ} term is as large as the clustering term for b ~ 4. Mean interloper fraction 0.21 (r), 0.10 (g).

Gaussian forecast of the cross-power significance (Limber, Planck 2018 ΛCDM, halofit; n(z) ∝ C(z) dV/dz from the published
completeness; clustering diluted by (1 − f_int); magnification with the s above; ACT 100 ≤ L < 1000 with the release N_L, Planck
100 ≤ L < 400 with the PR3 N_L):

| | b = 3 | b = 4 | b = 5 | b = 6 |
|---|---|---|---|---|
| r-drops × ACT | 9.5 | 9.6 | 9.5 | 9.4 |
| r-drops × Planck | 5.8 | 5.7 | 5.5 | 5.3 |
| g-drops × ACT | 13.0 | 12.3 | 11.7 | 11.2 |
| g-drops × Planck | 6.9 | 6.3 | 6.0 | 5.7 |

**Gate: passed** (z ~ 5 forecast ≥ 5σ on either map; the brief's stop line was 3σ). The forecast is flat in b because the
magnification term carries much of the signal; the bias therefore leans on s, which the fit carries as a prior (below). The z ~ 4
× Planck forecast (6-7σ in C_ℓ) is the scale against which Miyatake's 5.1σ (a stacked profile, θ ≤ 20′) is checked.

## Fixed before the first cross-spectrum (2026-11-02 02:55 UTC)

**Samples.** GOLDRUSH S18A Wide-layer masked catalogues (`masked_{riz,gri}_W_W0{2..7}.cat`), the published colour cuts already
applied by the catalogue (r-drops: r − i > 1.2, i − z < 0.7, r − i > 1.5(i − z) + 1.0, g < 2σ, z S/N > 5; g-drops: g − r > 1.0,
r − i < 1.0, g − r > 1.5(r − i) + 0.8, i S/N > 5). Magnitudes: the catalogue's aperture-corrected 2″ `convolvedflux_0_20` (its
total magnitudes) less `a_*`. *r*: every r-dropout. *g*: i < 25.5 (Miyatake's cut). Randoms nulls: Poisson points on the galaxy
window with each field's count.

**Masks.** Galaxy window at nside 1024: the fraction of each pixel's four nside-2048 sub-pixels holding ≥ 1 g-dropout (any
magnitude; ~4.8 per sub-pixel, so an empty sub-pixel is a mask hole or edge), pixels below 0.5 dropped; δ = N / (n̄_field w) − 1
with n̄ per HSC field. The same window for r and g (the masks are the same bright-star/inputcount masks). κ windows: ACT baseline
mask squared (degraded from 4096), Planck PR3 binary mask rotated to equatorial; both restricted to the HSC region (2° around the
fields). "Planck outside ACT" = W06 only.

**Multipoles.** Cross: ACT 100 ≤ L < 1000 (edges 100, 200, 300, 400, 550, 700, 850, 1000), Planck 100 ≤ L < 400 (100, 150, 200,
300, 400). Auto: 100 ≤ ℓ < 1000, ACT's edges. The fields are 6-10° across; ℓ < 100 is not used.

**Estimator.** Pseudo-C_ℓ with NaMaster's mode coupling (nside 1024, NaMaster lmax 2000; κ alm to L = 3000 before the map), galaxy
pixel window in the model. Auto: Poisson noise subtracted at the coupled level (mean over the sky of w / n̄) and decoupled. ACT
transfer (reconstruction × input over input × input, the first 60 release simulations, the same windows) applied per band.

**Model.** Planck 2018 ΛCDM, halofit, Limber (`qso_z4_theory.py`). LBG n(z) ∝ C(z) dV/dz with C(z) read off Harikane et al. 2022
Fig. 4 (⟨z⟩ ≈ 4.9 and 3.8). With A the amplitude of matter fluctuations at the LBG redshifts relative to Planck ΛCDM (A² scales
every term whose kernel is the LBG n(z)):

    C^{κg} = A² (1 − f) b C^{κm}_LBG + f b_int C^{κm}_int + (5s − 2) C^{κμ}
    C^{gg} = A² (1 − f)² b² C^{mm}_LBG + f² b_int² C^{mm}_int + 2 f b_int (5s − 2) C^{mμ}_int + 2 A² (1 − f) b (5s − 2) C^{gμ}_LBG + (5s − 2)² C^{μμ}

Nuisances with Gaussian priors, marginalised by Monte Carlo over 400 prior draws: f = the count-weighted Table 4 (Wide) fraction
(0.21 r, 0.10 g) ± 0.10 (the COSMOS2015 cross-check in the same paper allows ~30 % at the faint end), clipped to [0, 0.6];
interlopers at z = 0.75 ± 0.25 (r; the 4000 Å break between r and i) and 0.45 ± 0.20 (g), b_int = 1.5 ± 0.5 (r), 1.3 ± 0.4 (g);
s = 1.23 ± 0.15 (r), 0.97 ± 0.15 (g): s is the sample-mean intrinsic count slope of the published z ~ 5 / z ~ 4 UV LF (galaxy +
AGN double power laws, Harikane et al. Table 6) at each object's magnitude (the selection depends on flux, not on colour, which
lensing leaves unchanged: R-139's amendment); its spread over the DPL/Schechter forms, z ± 0.3 and m ± 0.1 is ≤ 0.08, and 0.15
is the floor for selection effects beyond total flux. Variant (systematic, not in the headline): n(z) shifted by +0.1 (the
spectroscopic redshifts sit higher than the simulation).

**Fit.** (1) Cross alone: b with A = 1. (2) Cross + auto: (b, A) on a grid, flat priors b ∈ [0, 20], A ∈ [0, 3]; A's marginal
posterior median and 68 % interval; σ8(z) = A σ8^{Planck}(z). Covariance: cross from the 400 ACT simulated reconstructions × the
real galaxy maps (Planck: 400 Gaussian realisations of the PR3 C + N); auto and cross-auto blocks from 400 correlated Gaussian
κ-galaxy realisations with Poisson sampling at the observed densities on the same windows (theory at b = 4 (r), 3 (g), A = 1);
Hartlap-corrected. Combined lensing: ACT on W02-W05 ⊕ Planck on W06 (no galaxy under two maps). Significance: √(χ²_null − χ²_best)
of the cross. Claim: the table of N, cross S/N, b (cross alone), b and A from cross + auto with intervals, for z ~ 5 and z ~ 4.

**Nulls and consistency (each with a rule fixed now).** (1) ACT curl × galaxies (covariance from the baseline sims). (2)
Randoms × κ. (3) ACT vs Planck on the same sky (Planck inside ACT): Δb with independent-sim errors. (4) Per-field: b (cross alone)
in W02, W03, W04, W05 on ACT; **if the fields' χ² of consistency has PTE < 0.05, or (3) differs at > 2σ, the combined b and A errors
are inflated by √(χ²/dof) of that test (R-139's gate lesson).** (5) CIB: the ACT CIB-deprojected reconstruction in place of the
baseline; |Δb| beyond its expected scatter (σ²_cib − σ²_base) is added to the error in quadrature. (6) The z ~ 4 reproduction:
g-drops × Planck PR3 passes if its cross S/N is 5.1 ± 1.5 (a different statistic over θ ≤ 20′, so a 1.5 window, not 1) and its
amplitude A(z ~ 4) is within 2σ of Miyatake's σ8 = 0.46 (+0.25 −0.20), i.e. A = 0.57 (+0.31 −0.25) for σ8^{Planck} = 0.811.

**Amendment (03:20 UTC; after the data spectra existed, before any fit or covariance; disclosed as post-data).** The fixed model
puts (5s − 2) on the whole sample, but interlopers dominate the bright, steep end (f = 0.45 at z = 23 where s ≈ 1.4), and their
counts are magnified only by matter below z ~ 0.75. The consistent coefficient on C^{κμ}_LBG is the LBG-weighted mean
⟨(1 − f(m))(5 s(m) − 2)⟩ over the objects: **3.20 (r), 2.52 (g)**, against (5s̄ − 2) = 4.13 and 2.85 in the fixed form (interloper
self-magnification neglected: z ~ 0.75 sources, C^{κμ} ~ 10 % of the z ~ 5 one, 5s_int − 2 ≈ 0 for faint red galaxies). The prior
moves with the nuisances: a = E1 − δf E2 + 5 δs E3 − 5 δf δs (E1 = ⟨(1−f)(5s−2)⟩, E2 = ⟨5s−2⟩, E3 = ⟨1−f⟩; δf ~ N(0, 0.10),
δs ~ N(0, 0.15)). What prompted it: a look at the first data cross-spectra against the fiducial model, which sat below the
magnification term alone; the error was found by re-deriving the term, not by tuning to the data. **Both forms are reported**
(`fit --model fixed` and `--model amended`); the amended one is the headline. ACT transfer from 30 input pairs (preliminary):
0.90-0.97 per band (the window's m⁴/m² loss, 0.93 averaged), not the cause of the low cross.

## Rows

| # | what | result |
|---|---|---|
| 1 | Data access (above): GOLDRUSH 5 min (parallel ranges), ACT data + masks 2 min, Planck PR3 1 min; ACT simulations streamed in the background | — |
| 2 | ACT simulations: 400 baseline reconstructions streamed 02:44-03:32 (48 min on two streams; one stalled transfer until the stall abort) and stored on the HSC region as float16 (`simstore`), 60 inputs (03:04); NaMaster 3.0.1 from pip after `apt-get install libcfitsio-dev libgsl-dev libfftw3-dev`. Galaxy window 397 deg² (352 effective): ~15 % above the paper's 307.9 deg² because a sub-pixel with one g-dropout counts as covered (holes smaller than 1.7′ are missed) | — |
| 3 | ACT transfer (reconstruction × input / input × input, 60 pairs, these windows): 0.920, 0.921, 0.922, 0.929, 0.960, 0.960, 0.947, 0.936 (± 0.013-0.020), applied to the model | the window's m⁴/m² edge loss (0.93 averaged) |
| 4 | **The cross alone** (`fit`, amended model; A = 1; nuisances marginalised; 400 ACT sims, 400 Planck Gaussian sims): r × ACT S/N 3.3 (χ²_null 16.0/8, PTE 0.04), b = −1.6 ± 2.4; ACT ⊕ Planck W06 S/N 3.1, b = −2.1 ± 2.4. g × ACT S/N 5.6 (χ²_null 36.6/8, PTE 1.4 × 10⁻⁵), b = 0.35 ± 1.0. Planck PR3 (L < 400): r S/N 1.3, g 1.8. The fixed model (5s − 2 on all): r −3.4 ± 2.7, g −0.05 ± 1.1 | the first CMB-lensing cross-power of z ~ 5 galaxies (3.3σ); at both redshifts the cross is what the LF-predicted magnification alone gives, with no room for the clustering term at A = 1 |
| 5 | Nulls and consistency: curl × r χ² 2.8/8 (PTE 0.95), × g 3.6/8 (0.89); randoms × ACT PTE 0.86 (r), 0.96 (g); × Planck 0.72, 0.07. ACT vs Planck on the ACT sky: Δb = 4.7 ± 4.1 (r), 1.7 ± 1.4 (g). **Per field (ACT): r W02 3.8 ± 6.0, W03 4.8 ± 3.2, W04 −9.3 ± 2.9, W05 −4.9 ± 3.2, χ² 12.4/3 (PTE 0.01); g −2.3 ± 1.7, 2.2 ± 1.2, −2.2 ± 1.2, 1.1 ± 1.1, χ² 9.7/3 (PTE 0.02)**: the rule inflates the errors ×2.03 (r), ×1.80 (g). **CIB-deprojected ACT: r b = −4.2 ± 2.1 (Δ −2.4), g −0.9 ± 0.7 (Δ −1.1), beyond the expected scatter (0.45, 0.15)**: 1.9 and 0.9 added in quadrature | curl and randoms pass; the fields disagree (W04 low in both samples, on the same sky); the CIB-deprojected map moves b down by ~1σ |
| 6 | z ~ 4 reproduction (rule fixed before: g × Planck S/N 5.1 ± 1.5): **S/N 1.8 at 100 ≤ L < 400: fails.** Post hoc, Planck to L < 1000 (ACT's bands, 300 Gaussian sims, `pvariant`): S/N 2.8; in Miyatake's own model (b C^{κm}, no magnification) b = 2.1 ± 0.8 (their halo mass 2.9 × 10¹¹ h⁻¹ M_⊙ implies b ≈ 3; their joint b = 6.5 +5.2 −4.2). r × Planck to L < 1000: S/N 2.5 | the amplitude reproduces within 1σ; the significance does not (their stack used θ ≤ 20′ on Planck to L = 2048, a jackknife covariance, and no magnification term) |
| 7 | **The joint result** (`fit`, amended; 400 joint κ-galaxy sims; `data/z5/fit_amended.md`): auto-spectra (shot noise subtracted) 9.6, 1.4, 1.4, 0.77, 0.60, 0.42, 0.40 × 10⁻⁶ (r), 3.6, 0.74, 0.75, 0.40, 0.34, 0.24, 0.23 × 10⁻⁶ (g); lowest band 3.8× / 2.6× the b = 4 / 3 model. Auto alone (A = 1): b A = 6.18 (−1.00 +1.10), χ² 7.4/6 (r); 4.17 (−0.47 +0.55), 3.3/6 (g). Cross alone: b A² = −2.05 ± 2.4 (r), 0.20 ± 1.0 (g). **Ratio A = −0.33 (−0.41 +0.38) (r), 0.05 (−0.25 +0.24) (g); with the inflation and CIB term ±0.8 and ±0.5: A = 1 at 1.6σ and 2.0σ (3.5σ, 4.0σ statistical).** Grid fit (flat b ∈ [0, 20]): A = 0.41 (−0.08 +0.11), 0.29 (−0.06 +0.10), b at the prior edge. At A = 1, b = b_auto, the cross needs a = −0.5 ± 1.0 (r), −0.8 ± 0.6 (g) (LF: 3.20, 2.52), s_eff ≈ 0.25. Fixed form of a: b A² = −3.8 (r), −0.15 (g); ratio A = −0.65, −0.04. The joint Gaussian sims' ACT cross σ is 1.4-2.8× the release sims' (their N_L convention or galaxy model), so only their correlations are used | a shortfall of lensing per unit clustering at both redshifts; not read as growth (see the finding's "What it means") |
| 8 | The finding section through the blind gate (the section alone in a file, the `blind-review.md` prompt verbatim, fresh context-free Opus, own literature search; hints: Miyatake et al. 2022, the forecasts 1904.13378, 2401.13166, 2511.22243, anyone who measured growth or σ8 at z ≥ 4, and magnification of dropout samples) | referee's verdict (verbatim): **"VERDICT: novel but not sound: the headline A values are a ratio of a numerator dominated by an assumed magnification coefficient that was never measured, over an auto-spectrum the note itself shows is contaminated by systematics, from a pipeline that failed its pre-set reproduction check, so the claimed 2σ shortfall is not supported."** Novelty granted for the z ~ 5 cross (closest arXiv:2103.15862; 1904.13378, 2401.13166, 2511.22243; 1204.2830; 1511.01352); the z ~ 4 amplitude is prior work. Strongest objection: an unvalidated, post hoc magnification coefficient over a systematics-contaminated auto-spectrum, from a pipeline that failed its pre-set reproduction | **a miss; PARKED** (physics: the magnification response and the auto-spectrum systematics, not a presentation point; brief: park on a physics point) |
