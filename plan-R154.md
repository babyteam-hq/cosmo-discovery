# Analysis plan: CMB-lensing masses of ACT DR6 SZ clusters in redshift bins, and the redshift evolution of the SZ mass bias

Fixed 2026-11-02, before any stack of the data map at the cluster positions.

## Question and claim

Galaxy weak lensing of ACT DR5 SZ clusters (arXiv:2512.18935, DES Y3) finds the SZ mass bias 1-b = M_SZ/M_true falling with
redshift: 0.81 (+0.08/-0.06) at 0.15 < z < 0.43, 0.58 (+0.06/-0.05) at 0.43 < z < 0.70. Their fit is
M_true/M_SZ = A_mass (M_SZ/3e14)^eta ((1+z)/1.45)^zeta, with zeta = 1.8 (+0.6/-0.5), so zeta > 0 at 99.95 %. Here M_SZ is the
ACT catalogue mass from the Arnaud et al. 2010 (A10, UPP) scaling. Galaxy lensing stops near z ~ 0.7-1. CMB lensing has its
source plane at z = 1100, so it can carry the test past z = 1.

The claim will be this table, with intervals:
- per redshift bin, R = M_lens/M_SZ (and 1-b = 1/R);
- the evolution exponent zeta, fitted with the same form as 2512.18935;
- both compared with 2512.18935's bins and exponent, and with HSC (arXiv:2407.08201: M_WL/M_UPP = 1.3 ± 0.2, 1.6 ± 0.2,
  1.6 ± 0.3 at mean z 0.32, 0.57, 0.88).

## Data

- **Catalogue.** ACT DR6 SZ cluster catalogue v1.0 (arXiv:2507.21459; LAMBDA `DR6_cluster-catalog_v1.0.fits`). It has 10,040
  clusters at SNR > 4, all optically confirmed with redshifts (4,133 spectroscopic, 5,907 photometric), 1,180 of them at z > 1.
- **Mass M_SZ.** The column `A10_M500c`: the A10 scaling relation, flat ΛCDM Ωm = 0.3, H0 = 70. It is the catalogue's posterior
  mass given the measured y and z, with the Tinker 2008 mass function at σ8 = 0.80 as the prior, so it is Eddington-corrected.
  The assumed intrinsic scatter is 0.2 in ln. This is the mass definition behind 2512.18935's 1-b.
  - The catalogue's default `M500c` uses a normalisation 0.6 times A10's. It is reported as a secondary ratio.
  - Each cluster's mass posterior is taken as lognormal, its width from the errPlus/errMinus columns, averaged with 5
    Gauss-Hermite nodes (R-122's posterior averaging).
- **Lensing.** The ACT DR6 lensing release v1 (Qu et al. 2024, arXiv:2304.05202; Madhavacheril et al. 2024, arXiv:2304.05203).
  - The baseline map is a profile-hardened, temperature-plus-polarization quadratic estimator. In it, the clusters detected at
    SNR > 5 by the lensing team's own Nemo run had beam-convolved UPP tSZ templates subtracted from the single-frequency maps
    before the reconstruction. Before unblinding, these clusters were inpainted instead; the release keeps that 6' mask under
    `other/`.
  - Variants used: tonly, polonly (tSZ-free), cibdeproj (CIB-deprojected, not profile-hardened) and curl (a null).
  - Simulations: the 400 baseline simulated reconstructions and their input κ fields.
  - Band: 20 ≤ L ≤ 2048, with a cosine taper.

## Sample (fixed)

- **Footprint selection (R-122's).** A cluster is kept if it is clean (mask ≥ 0.99) everywhere within 30' and ≥ 80 % clean
  within 2°. This keeps 6,263 clusters.
  - 2,824 of the kept clusters lie inside the release's 6' SZ-cluster mask, the list whose tSZ templates were subtracted. This
    includes 2,650 of the 3,620 kept clusters at SNR > 5.
  - All 6,263 form the primary sample.
  - The two halves (in the SZ-cluster mask, so tSZ-subtracted; outside it, so tSZ only profile-hardened) are each reported per
    bin as a consistency test.
- **Redshift bins.**

| bin | z range | kept | in the SZ mask | median A10_M500c |
|---|---|---|---|---|
| 1 | 0 – 0.5 | 2,459 | 1,105 | 1.96e14 |
| 2 | 0.5 – 0.8 | 2,125 | 1,004 | 2.08e14 |
| 3 | 0.8 – 1.0 | 922 | 401 | 1.93e14 |
| 4 | 1.0 – 2.0 | 757 | 314 | 1.76e14 |

- **Diagnostic subsamples** (reported, not the claim): fixed_SNR > 5, fixed_SNR 4-5, spectroscopic redshifts only, z > 1.2.

## Model and statistic

- **Predicted κ per cluster.** An NFW one-halo term plus a linear two-halo term, computed in the catalogue's cosmology.
  - One-halo: Duffy et al. 2008 c(M, z).
  - Two-halo: Tinker 2010 bias, σ8 = 0.81.
  - The prediction is averaged over the cluster's mass posterior and filtered exactly as the map: band, taper and HEALPix pixel
    windows, with nside 1024 for the annuli inside 30' and 512 outside.
- **Stacked data vector.** The mean of the per-cluster masked annulus means in 9 annuli: 0-5, 5-10, 10-15, 15-20, 20-30, 30-50,
  50-80, 80-110 and 110-150'.
- **Null-mean correction.** The mask-matched random stack (100,000 random positions passing the same selection) is subtracted.
- **Miscentring.** A Rayleigh offset of the SZ position with σ_off = 0.5', applied as exp(-L²σ_off²/2) on the one-halo κ(L).
  Varied to 0' and 1.0'. As a check, optical centres are used where the catalogue gives them.
- **Per-bin statistic.** R_b = M_lens/M_SZ.
  - Every cluster's mass in the bin is multiplied by R_b and its profile rebuilt; the two-halo amplitude of the bin is free and
    profiled out.
  - The 68 % and 95 % intervals come from Δχ² < 1 and < 3.84 on a grid in R.
  - Also reported: the one-halo amplitude A_1h at R = 1, with A_2h free.
- **Trend statistic.** ln R(z) = ln R_0 + ζ ln((1+z)/1.45).
  - It is fitted to all four bins jointly: each cluster's mass is scaled by R(z_i) at its own redshift, with one free two-halo
    amplitude per bin, profiled.
  - The χ² is over the 4 × 9 stacked vector with the joint covariance from the 400 simulations.
  - The interval on ζ comes from the profile Δχ² over R_0.
  - Compared with: ζ = 0 (constant bias); ζ = 1.8 (2512.18935); and the 2512.18935 law extrapolated to the z > 1 bin.
- **Mass dependence.** The 2512.18935 law has η = −0.18 ± 0.15. Its effect on R at the bins' median masses (2512.18935's pivot
  is 3 × 10^14) is reported. η is not fitted.

## Covariance and nulls

- **Covariance.** From the 400 baseline simulations, stacked at the real cluster positions with the same selection and filter;
  Hartlap-corrected; joint across bins. It is checked against 100 rotated-map stacks.
- **Nulls.**
  - The random-position stack.
  - 100 rotated-map stacks: the map and mask turned together in RA by 25-335°, the clusters fixed, the selection re-applied. The
    fitted R and ζ on each rotation must scatter about 0, consistent with the simulations.
  - The curl reconstruction stacked on the clusters: amplitude consistent with 0.
  - The simulations' mean at the cluster positions: consistent with 0.
- **Normalisation.** R-122's check of the release map against its simulation inputs found T(L) = 0.98-1.02 over 20-2048. It is
  reused, and repeated at the SZ-mask positions: output against input κ of the simulations, stacked at the cluster positions in
  the SZ mask and outside it. This tests whether the reconstruction responds to lensing differently where the templates were
  subtracted. (Whether the simulations themselves had templates subtracted is not stated; the test is reported for what it is.)

## Systematics and how each is bounded

1. **tSZ leakage into the reconstruction.**
   - Baseline: profile hardening, plus template subtraction for the SNR > 5 clusters.
   - Bounds, per bin: polonly (tSZ- and CIB-free, about 2× the noise), tonly, and cibdeproj, which has no hardening and is the
     tSZ-sensitive extreme (R-122 found it 0.14 low).
   - Also per bin: the in-mask (template-subtracted) half against the out-of-mask half.
   - The baseline − polonly difference per bin, with its error from rotations on both maps, is the tSZ-leakage bound. Its
     redshift dependence is carried into ζ as a systematic.
2. **CIB from cluster members at z > 1.** The cibdeproj − baseline difference per bin; and polonly.
3. **The SZ selection (Eddington) bias.** The catalogue's posterior mass is conditioned on the measured y, so selecting on y does
   not bias it if the prior and scaling relation are right. Assuming A10 with b = 0 when the truth has 1-b ≈ 0.6-0.8 evaluates
   the mass-function prior at the wrong mass.
   - **Mock.** Draw Tinker 2008 halos over z in 0-2, with this model's σ(M). Assign y from A10 with 1-b ∈ {0.6, 0.75, 0.9}, a
     log-normal scatter of 0.2, and Gaussian measurement noise at the catalogue's own per-cluster noise. The noise comes from the
     empirical err_y_c distribution per redshift bin, so the noise and completeness are the catalogue's own.
   - Select q > 4. Compute the A10 posterior the catalogue's way (b = 0, scatter 0.2, Tinker prior).
   - Measure the expected ⟨M_true^0.87⟩^(1/0.87)/⟨M_post⟩ per bin against the input 1/(1-b). The exponent 0.87 is the stack's
     dlnA/dlnM.
   - The per-bin offset from the input is a correction applied to R with its spread as a systematic. The uncorrected masses
     (`A10_M500cUncorr`) are the extreme.
4. **Purity and redshift errors.**
   - False or misassigned clusters dilute R. The fixed_SNR > 5 against 4-5 split per bin bounds the low-SNR contamination, and
     the spectroscopic subsample bounds catastrophic photo-z errors.
   - A contamination fraction f lowers R by (1-f)^(1/0.87). No correction is applied; the split is reported.
5. **Profile model.**
   - c(M) × 0.7 and × 1.3.
   - The two-halo term with halo exclusion inside 1.5 r200c.
   - The fit restricted to the annuli inside 30'.
   - The spread of these on R per bin and on ζ is the model systematic, added in quadrature to the intervals.
6. **Quadratic-estimator non-linearity for massive clusters.** A mass split per bin (A10_M500c above and below 3 × 10^14); the
   L ≤ 2048 band.

## Forecast (feasibility; done before the plan was fixed)

- **Method.** A Fisher calculation over the 9 annuli, clusters independent, with the release N_L plus the lensing signal power.
  It is calibrated by the same calculation on R-122's eRASS1 cosmology sample (1,907 clusters) against its realised 400-sim error
  (0.110); the factor is 1.18.
- **σ(A_1h), two-halo free, A10 template.** 0.12, 0.155, 0.27 and 0.34 in bins 1-4; dlnA/dlnM = 0.80-0.87.
- **If 1-b = 0.74 everywhere:** σ(ln R) = 0.12, 0.14, 0.24, 0.31 and σ(ζ) = 0.51, so ζ = 1.8 against 0 is separated at 3.5σ.
  The z > 1 bin alone separates the extrapolated 2512.18935 law (R = 3.4) from R = 1.35 at 3.0σ.
- **If the 2512.18935 law is true:** σ(ζ) = 0.34, a 5.3σ separation.
- **Outside the SZ mask only:** σ(ζ) = 0.81, so 2.2σ.
- **Anchor.** arXiv:2406.14754 detected a z > 1 stack of ACT-DR5 MCMF clusters on this map at 4.3σ, with a 5' aperture,
  600 < L < 2100 and about 480 clusters. Scaled by √N to 757 clusters, that is ~5.4σ before any gain from the profile fit, i.e.
  σ(ln R) ≈ 0.21 at z > 1.

## Error budget and what decides

- Statistical errors come from the 400 simulations.
- The systematics from items 3 and 5 are added in quadrature per bin.
- Item 1's redshift-dependent part goes into ζ: half the spread of ζ across baseline/tonly/cibdeproj, plus the polonly
  difference if significant.
- The pre-stated test: is ζ = 0 inside the 95 % interval, and is the 2512.18935 extrapolation (ζ = 1.8) inside it?

Trials: every variant above is computed and reported. The claim is the baseline primary-sample table.

## Addendum (2026-11-02, after the method check, still before any data stack at the clusters): how each point is closed

The method check returned "fixable". Every point it named is closed below. The scripts are `act_szlens_inject.py`,
`act_szlens.py selmock` / `cib`, and `act_szlens_fit.py`, which applies the corrections; rows I1-I5 of `results-R154.md` give
the numbers.

**1. Reconstruction response, tSZ leakage and estimator non-linearity: injection.**
- **The pipeline.** The release's reconstruction is not runnable here, so it is emulated in flat sky:
  - a 256' box at 0.5' pixels;
  - CAMB CMB, lensed by remapping through the cluster's NFW deflection;
  - A10/UPP tSZ at 90 and 150 GHz from the cluster's A10 Y500 at M_SZ;
  - ACT beams and noise, scaled so the emulated N_0 matches the release's tonly N_L;
  - the release's subtraction for the in-mask clusters (UPP template at the nearest of 15 Nemo scales, amplitude
    matched-filter-fitted per frequency on the noisy map);
  - the inverse-variance coadd;
  - the Hu-Okamoto TT estimator over CMB 600-3000, profile-hardened with the release's choice of profile, u(l) = √C_l^tSZ
    (here a Tinker × A10 halo model);
  - the stack's filter and its inner five annuli;
  - the amplitude fitted with the real fit's per-cluster covariance.
- **Normalisation.** Checked on Gaussian κ: the response is 1 at L > 300 within ±5 %.
- **Lensing response.** Paired realisations, 100-150 pairs, at each bin's median M_SZ, at R = 1.35 and 3.4, plus two
  high-mass cases (4-4.5 × 10^14).
- **tSZ leakage.** The mean bias is deterministic, QE(s, s). It is computed for UPP truth and for Planck-2013-GNFW truth (a
  profile the template does not match), hardened and not. The template-subtracted residual comes from Monte Carlo
  (40-60 realisations).
- **Variants.** Hardening profiles of 1.5' and 3' UPP instead of √C_l^tSZ.
- **Into the baseline map.** The temperature leakage enters with w_T = ⟨N_L^MV/N_L^T⟩ over L = 500-2000, from the release N_L
  curves; polarization carries none.
- **The correction.** Per bin, the leakage is weighted by the bin's in-mask (template-subtracted) fraction. The data are
  corrected as d → (d − l_b m1_b(R=1)) / r_b.
- **Systematic.** Each variant's shift, added in quadrature.

**2. CIB separated from tSZ.**
- **The amplitude is measured.** The mean Planck PR3 857 GHz aperture flux of each bin's clusters (0-8' minus a 15-25' ring;
  the null-mean from 20,000 randoms), and the 545/857 ratio, after removing the A10 tSZ increment at 545.
- **The SED.** A modified blackbody whose 545/857 ratio matches the measured one (β 1.5-2, T 20-25 K). The extrapolation to 90
  and 150 GHz gives the bin's fluxes.
- **The injection.** A point source and a 1' Gaussian at the centre, with the cluster's tSZ, through the hardened estimator, and
  through the subtraction (whose amplitude fit absorbs part of the CIB).
- **Systematic.** The SED and shape spread, and the 857 flux error.
- cibdeproj − baseline is no longer used as the CIB bound. It is reported only.

**3. Optical confirmation in the selection.** A variant of the selection mock confirms clusters at z > 1 with probability
1/(1 + (M_50/M_true)^4), with M_50 = 1.0 and 1.5 × 10^14. Its ratio to the central B is a systematic.

**4. The Eddington / selection mock, extended and self-consistent.**
- The grid is 1-b ∈ {0.25, 0.4, 0.6, 0.75, 0.9}, with the scatter varied over 0.1, 0.2 and 0.3 and σ8 over 0.77, 0.80 and 0.83.
- Each bin's own y-M slope (1.22-1.47) is measured from fixed_y_c against A10_M500cUncorr. The noise per mock halo is drawn from
  the catalogue's own per-cluster k = q/M^slope.
- The mock is run through the real fit: posterior nodes and the two-halo term free.
- B(R_true) = R_fit/R_true is inverted per bin (R_fit = B(R_true) R_true). The variants' ratios are systematics, one-sided
  where they are.

**5. Purity and photo-z.**
- **Contamination.** The DR6 paper gives a 1-6 % chance association with optical clusters, and noise false positives that are
  negligible above q = 5.5. The central correction is 0. The one-sided systematic allows up to 4 % (6 % at z > 1); it lowers R
  by (1-f)^(1/0.85).
- **Photo-z.** The catalogue's photo-z errors are 0.010, 0.021, 0.033 and 0.059 (medians, bins 1-4). The κ kernel varies
  slowly with z, and the effect enters through the angular size and Σ_crit. Bin migration across a bin edge is ≲ σ_z/Δz_bin of
  the clusters near the edge.
- Both are reported. The q > 5 subsample (higher purity) and the spectroscopic subsample are stated as checks, not as the
  budget.

**6. Miscentring, calibrated.** The SZ-optical offsets in the catalogue have a Rayleigh σ of 0.30-0.49' (bins 4 to 1); these
include optical-centre errors, so they are an upper bound. The central value is σ_off = 0.5'; the systematic is the shift to
0'.

**7. Mass definition.**
- The DR6 A10_M500c against DR5 M500c (the 2512.18935 mass) for 3,331 cross-matched clusters: median ratio 1.03, 1.02, 1.01 and
  1.02 in the four bins, mean ln ratio 0.000-0.017. No correction is applied.
- ζ is also reported with η = −0.18 (M_SZ/3e14)^η.

**8. Error assembly.**
- **Per bin:** the statistical interval from the 400 simulations, with these systematics added in quadrature on ln R:
  - leakage variants;
  - the response error and its hardening variant;
  - the selection variants;
  - contamination;
  - c(M) × 0.7/1.3;
  - the inner-30' fit (two-halo shape);
  - miscentring.
- **ζ:** each z-dependent systematic is applied as its own variant of the whole fit. The quadrature sum of the shifts in ζ is
  added to ζ's interval. The test (ζ = 0 and ζ = 1.8 against the 95 % interval) uses the combined interval.
