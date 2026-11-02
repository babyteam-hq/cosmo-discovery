# Results (R-139, cosmo-deep-15, 2026-11-01)

The linear bias of quasars at 3.5 < z < 5 from the cross-power of their overdensity with CMB lensing convergence (ACT DR6,
Planck), and with it the discrimination between the quasar auto-correlation answer (b ~ 13; Shen et al. 2007, astro-ph/0702214)
and the quasar-galaxy cross-correlation answer (b ~ 3-6; He et al. 2018, arXiv:1704.08461). Every number from `qso_z4_feas.py`
(step 0) and `qso_z4_lens.py` (the measurement; usage in its docstring); theory in `qso_z4_theory.py`. Data under `data/`
(gitignored); logs in `data/qso_z4/`. The trajectory check ("VERDICT: novel: no prior work answers this question", 2026-11-01)
is in the brief and in `LOG.md`; not rerun.

## Step 0: feasibility (22:43 UTC, before any κ map was crossed with any quasar map)

Refetch: ACT DR6 baseline alm + mask (1 min), SDSS DR16Q (1.7 GB, 1 min), DESI DR1 `QSO_full_HPmapcut` (2.0 GB, 5 min),
Planck PR3 `COM_Lensing_4096_R3.00` from the PLA (720 MB, 1 min). **Planck PR4 lensing is not reachable from this cloud**: not on
the PLA, 404 at NERSC (`planck2020/PR4_lensing`, both `project` and `cfs` paths, with and without `data/`), and the GitHub release
assets 404 through the proxy (as R-129 found). The Planck leg uses PR3 (MV, 8 ≤ L ≤ 400, its own conservative range).

| sample | 3.5<z<5 | in footprint, \|b\|>20° | × ACT DR6 overlap | × Planck PR3 overlap |
|---|---|---|---|---|
| SDSS DR16Q (all quasars, `Z`) | 9,761 | 9,435 (9,984 deg²) | 3,060 on 3,359 deg² | 9,153 on 9,653 deg² |
| DESI DR1 (`QSO_full_HPmapcut`, good z) | 14,402 | 14,141 (9,838 deg²) | 6,200 on 3,972 deg² | 13,520 on 9,453 deg² |

3,403 DR16Q high-z quasars have a DESI high-z counterpart within 2". Of DR16Q's 9,761 only 1,615 are BOSS CORE (uniform XDQSO)
targets; 5,788 are BONUS: the DR16Q sample is not homogeneously selected (below).

Gaussian forecast of the cross-power significance (Limber, Planck 2018 ΛCDM, halofit; the spectroscopic n(z); magnification
included with s = 0.4 / 0.8, i.e. 5s−2 = 0 / 2; ACT 40 ≤ L ≤ 1000 with the release N_L, Planck 8 ≤ L ≤ 400 with the PR3 N_L):

| | b = 3 | b = 6 | b = 13 |
|---|---|---|---|
| DR16Q × ACT | 1.0 / 1.8 | 1.9 / 2.8 | 4.1 / 4.9 |
| DR16Q × Planck | 1.1 / 2.2 | 2.1 / 3.2 | 4.5 / 5.5 |
| DESI × ACT | 1.4 / 2.6 | 2.8 / 4.0 | 5.8 / 6.9 |
| DESI × Planck | 1.3 / 2.7 | 2.6 / 3.9 | 5.4 / 6.7 |

**Gate: passed in its second form.** b = 13 forecasts 4-7σ in every survey-map pair; b = 3 does not (1-3σ). The measurement
cannot detect a b ~ 3-4 signal at high significance, but it separates the two published answers: the expected Δχ² between b = 13
and b = 4 is (S/N₁₃ − S/N₄)² ≈ 9-16 per pair before combination. Continued. The referee's caution holds: magnification is
as large as the clustering term per unit (C^{κμ} ≈ 1.7 C^{κm} at L = 100 for unit 5s−2), so s matters at the ±1 level in b.

## Fixed before the first cross-spectrum (2026-11-01 22:50 UTC)

**Samples.**
- *DESI DR1*: `QSO_full_HPmapcut.dat.fits` (LSS v1.5), good QSO redshift (`Z_not4clus` finite), 3.5 < z < 5. Data weight
  `WEIGHT_ZFAIL / FRACZ_TILELOCID`; the residual completeness by `NTILE` (the weighted all-z good-QSO count per target: 0.644,
  0.687, 0.677, 0.668, 0.663, 0.659, 0.654 for NTILE 1-7, measured before any cross-spectrum) goes on the randoms. Footprint
  and vetoes from the DR1 clustering randoms `QSO_{NGC,SGC}_0_clustering.ran.fits` (angular positions only; their `NTILE`).
- *SDSS DR16Q*: every quasar with 3.5 < `Z` < 5 (`Z` is the visual-inspection redshift for 92 %). No public randoms cover
  DR16Q's heterogeneous selection; the footprint is the plate coverage traced by all 750k DR16Q quasars (nside-64 pixels holding
  > 30 % of the mean count), unit weights. Its selection varies with chunk and imaging (BONUS); this adds variance to the cross
  power and biases it only through selection correlated with κ (dust, stars), tested by the E(B−V) variant below.
- *Combined*: the joint fit of the DESI and DR16Q cross-spectra with their joint covariance from the same κ simulations (the
  sims crossed with the fixed quasar maps carry the 3,403 shared quasars' correlation).
- *End-to-end check*: DESI DR1 at 2.2 < z < 3.5 through the same pipeline, against the published DESI DR1 × Planck PR4 bias.

**Masks.** Quasar maps at nside 512, weight = the footprint coverage fraction (pixels below 0.5 dropped), δ = N_w / (n̄ f) − 1
per region (DESI N/S photometry, DR16Q NGC/SGC each with its own n̄); |b| > 20°; E(B−V)_SFD < 0.15 is *not* applied in the
baseline (variant). ACT: the release baseline mask degraded to nside 512, squared as the κ-field window (the quadratic estimator
weights by the mask twice); Planck PR3: its mask rotated to equatorial, 1° apodized (as `kappa_templates.py`).

**Multipoles.** ACT 40 ≤ L < 1000 in bins with edges 40, 70, 100, 150, 200, 300, 400, 550, 700, 850, 1000; Planck 8 ≤ L < 400 with
edges 8, 40, 70, 100, 150, 200, 300, 400. Variant: ACT L < 600.

**Estimator.** Pseudo-C_ℓ with the mode-coupling matrix (NaMaster 3.0.1, built here), the quasar map's pixel window deconvolved;
the ACT transfer function (reconstructed × input over input × input) measured on the release simulations with the same
windows and applied if it differs from 1 by more than 2 %.

**Model.** C^{κg}(L) = b C^{κm}(L) + (5s − 2) C^{κμ}(L), Limber, Planck 2018 ΛCDM, halofit P(k), the sample's spectroscopic
n(z), binned through NaMaster's bandpower windows. s = d log₁₀ N(<m)/dm of the sample's own counts at its flux limit (DESI:
dereddened r = 23.0; DR16Q: dereddened r = 21.85, the BOSS/eBOSS r limit), over the faintest 0.25 mag; uncertainty the
Poisson error plus the spread over Δm = 0.15-0.5, floor 0.1; propagated as a Gaussian prior on 5s − 2 (the model is linear in
b and 5s − 2, so the marginalisation is exact: C → C + σ²_{5s−2} t tᵀ with t the binned C^{κμ}).

**Fit.** b free (flat), χ² with the simulation covariance (Hartlap-corrected): ACT from the release's baseline simulated
reconstructions; Planck PR3 from 400 Gaussian κ realisations of the PR3 C^{κκ}+N^{κκ} (`nlkk.dat` column 3; the PR3 simulations
are not on the PLA), masked identically. Claim: b with its 68 % interval per survey and map, the combined value, the cross-power
significance (√χ² of the null against the best fit's Δχ²), and the likelihood ratio L(b = 13)/L(b = 4) with s marginalised.

**Nulls.** ACT curl reconstruction × quasars (covariance from the curl simulations); the quasar randoms (a random subsample
made into a δ map with the same mask) × κ; ACT − Planck on their overlap (difference of b; the covariance from independent
sims, conservative because the shared κ signal cancels in the difference); the 2.2 < z < 3.5 end-to-end check.

**Amendment (22:50 UTC, still before any cross-spectrum): the magnification slope.** The fixed estimator (the cumulative count
slope at the flux limit) gave s = 0.008 (DESI, r = 23.0) and 0.046 (DR16Q, r = 21.85): neither sample is cut by a hard flux limit
at those magnitudes (DESI's high-z counts turn over at r ≈ 21, DR16Q's at r ≈ 20.4, a selection-completeness roll-off), and for a
selection that depends on magnitude only (colours are unchanged by lensing) the response is the sample-weighted *intrinsic*
differential count slope, s_eff = ⟨d log₁₀ n_int/dm⟩ over the selected objects (integrate the selection by parts: the
completeness P(m) drops out). Replaced by: s(r) from a quadratic fit of log₁₀ n(r) of DESI's own 3.5 < z < 5 quasars over
18.0 ≤ r ≤ 20.75 (taken as complete), and s = 0.25 beyond r = 20.75 (the z ~ 4 faint-end LF slope α = −1.3 to −2.0, i.e. s =
0.12-0.40, the range carried as the error); s_eff the mean over each sample's objects; error = half the faint-end range, the fit-range
shifts (20.5/21.0, 17.5/18.5) and a 0.05 floor. Result (`prep`): **DESI s = 0.33 ± 0.10 (5s−2 = −0.35 ± 0.49; 60 % of the sample
beyond r = 20.75), DR16Q s = 0.40 ± 0.07 (5s−2 = 0.00 ± 0.34), DESI 2.2 < z < 3.5 s = 0.30 ± 0.12.** NaMaster bins carry an overflow
bin to L = 1500 (the field band limit), dropped from every output.

## Rows

| # | what | result |
|---|---|---|
| 1 | Data access and refetch (step 0 above): 9 min for everything but the simulations; ACT baseline simulations streamed four at a time (400 × 160 MB, 85 min with the consumer, `act_dr6_simstream.sh`), 60 sim inputs for the transfer, 100 curl simulations. NaMaster 3.0.1 builds here once `apt-get install libcfitsio-dev libgsl-dev libfftw3-dev` runs (R-129 found it unbuildable; the headers are installable after `apt-get update`). The ACT data alm hold NaN at L = 0, 1 (set to zero) | — |
| 2 | Magnification slope (amended before any cross-spectrum, above): DESI s = 0.33 ± 0.10, DR16Q 0.40 ± 0.07, DESI 2.2 < z < 3.5 0.30 ± 0.12 (the published DESI DR1 × PR4 fitted values for 2.1-3.5 are 0.19-0.22 ± 0.05-0.06, arXiv:2506.22416). At z ≈ 3.8 C^{κμ}/C^{κm} = 1.6 at L = 100, 1.05 at 300, 0.95 at 700 per unit 5s−2 | consistent with the published DESI value at 2-3.5 within 1σ |
| 3 | ACT transfer (reconstructed × input over input × input, 60 sims, the DESI window): 0.97-0.985 ± 0.004-0.010 in every band; applied (> 2 % in 7 of 10 bands) | — |
| 4 | **End-to-end check**, DESI 2.2 < z < 3.5 (326,980 weighted quasars): b = 2.98 ± 0.53 (ACT), 3.37 ± 0.53 (Planck PR3), 3.23 ± 0.48 (ACT ⊕ Planck outside ACT), χ² PTE 0.16-0.83, S/N 6.6-8.8; ACT and Planck on the same sky 2.98 ± 0.53 vs 2.91 ± 0.67. Published DESI DR1 × Planck PR4 (arXiv:2506.22416 Table/Sec. 6): 2.53 ± 0.23 (2.1-2.5), 2.93 ± 0.36 (2.5-3.5); their linear-theory fit 2.86 and 3.43 | **passes**: the pipeline recovers the published bias at the redshifts below the target. DESI × PR4's top bin ends at z = 3.5 (2.5 ≤ z ≤ 3.5), so it carries no z > 3.5 information |
| 5 | **The measurement, 3.5 < z < 5** (`fit`, 400 ACT sims, 400 Planck Gaussian sims, Hartlap-corrected; s marginalised): see the table under the finding section. Per map: DESI × ACT b = 5.45 ± 2.66 (S/N 2.0), DR16Q × ACT −0.91 ± 3.64, DESI × Planck 4.49 ± 2.70, DR16Q × Planck 3.48 ± 3.07; Planck outside ACT: DESI 11.6 ± 3.7, DR16Q 9.0 ± 3.5. **Combined (DESI + DR16Q, joint covariance, ACT ⊕ Planck outside ACT, every quasar under one lensing map once): b = 6.21 ± 1.92, χ² 26.5/33 (PTE 0.78), cross-power S/N 3.2; L(b=13)/L(b=4) = 3.7 × 10⁻³.** ACT L < 600: 5.41 ± 1.98, L ratio 8.4 × 10⁻⁴. DESI alone (ACT ⊕ Planck outside ACT) 7.55 ± 2.23, L ratio 0.18; DR16Q alone 4.19 ± 2.58, 3.0 × 10⁻³ | the first CMB-lensing bias of z > 3.5 quasars; 3.2σ; b = 13 sits 3.5σ above the combined value; DESI alone does not separate the two answers |
| 6 | Nulls. Randoms × κ: DESI b_null = −0.56 ± 2.34 (ACT), −2.1 ± 2.5 (Planck), −6.2 ± 3.3 (Planck outside ACT); DR16Q −2.1 ± 3.1, 0.8 ± 3.1, 2.2 ± 3.9; χ² PTE 0.07-0.72. E(B−V) < 0.15: no change (≤ 0.3 in b; the cut removes < 0.2 % of the footprint at \|b\| > 20°). **Curl × DESI χ² 3.7/10 (PTE 0.96); curl × DR16Q χ² 20.5/10 (PTE 0.03), carried by L = 550-1000 (−3.4σ, −2.0σ, −1.7σ in the top three bands; covariance from 100 curl sims)**; curl × DESI 2.2-3.5 PTE 0.20. ACT − Planck on the ACT footprint: DESI 5.45 ± 2.66 vs −3.29 ± 3.78 (Δ = 8.7 ± 4.6, 1.9σ); DR16Q −0.91 ± 3.64 vs −6.29 ± 5.23 (0.8σ). Planck inside vs outside ACT for DESI: −3.3 ± 3.8 vs 11.6 ± 3.7 (2.8σ between independent skies) | one null marginal (DR16Q curl at high L; the L < 600 variant, which drops those bands, moves the combined b by −0.8); the Planck legs scatter more than their errors between the two skies |
| 7 | Post hoc (not in the plan): DR16Q at i < 20.2 (Shen et al. 2007's limit, 4,763 quasars, s = 0.56 ± 0.05 from the hard cut) on Planck PR3: b = 3.4 ± 4.1 (L ratio 0.07); outside ACT 10.0 ± 5.0, inside −6.6 ± 7.3 | uninformative: the bright sample alone cannot separate 4 from 13 |
| 8 | The finding section through the blind gate (the section alone in a file, the `blind-review.md` prompt verbatim, fresh context-free Opus, own literature search; hints: the DESI DR1 × PR4, Quaia, eBOSS quasar-lensing works, Shen et al. 2007, He et al. 2018) | referee's verdict (verbatim): **"VERDICT: novel but not sound: the headline b = 6.2 ± 1.9 and likelihood ratio of 270 come from averaging mutually inconsistent sky and lensing-map legs (Planck inside vs outside ACT differ by 2.4–2.8σ, and the Planck legs outside ACT alone favour b ≈ 13) without inflating the error, and the ratio compares samples of different luminosity to Shen et al.'s bright quasars, so the preference for the cross-correlation answer is not supported."** Novelty granted: no CMB-lensing quasar bias above z = 3.5 found (closest: arXiv:2506.22416, 2512.17865, 2507.08798, 1809.04196, 1207.4543, 1902.06955, 2011.01234, 2410.23372; clustering at z > 5, arXiv:2602.02778). Strongest objection: the legs disagree by sky region and the combination does not carry that scatter | **a miss; PARKED** (physics: a sky-dependent difference between the legs, not a presentation point; brief: park on a physics point). The referee's test, computed once for the record: DESI + DR16Q on the ACT sky 3.30 ± 2.42 vs on Planck outside ACT 10.40 ± 2.86, Δ = 7.1 ± 3.7 (1.9σ, independent skies; χ² 3.6 for 1 dof); scaling the combined error by √3.6 gives b = 6.2 ± 3.6, b = 13 at 1.9σ, L(13)/L(4) = 0.21: not a decisive preference |
