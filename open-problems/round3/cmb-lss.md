# Round 3: CMB and LSS growth, quick-kill survivors

Checked 2026-11-01 with WebSearch (extended), the arXiv API (export.arxiv.org) and arXiv's
site search. Two candidates survive. A "not found" means not found by the queries listed.

**Killed or dropped before listing:**
- tau from WMAP+LFI (Cosmoglobe): WMAP+LFI tau without HFI already exists (Planck 2015
  0.071±0.013; Planck 2018 0.069±0.011; BeyondPlanck 2011.05609). A rerun moves the error
  from about 0.011 to 0.010 and does not separate 0.055 from 0.09.
- DESI DR1 spectroscopic LRG x ACT DR6 lensing: done (2407.04606; 2505.20656 gives S8=0.808±0.017).
- PR4 parity asymmetry in EE/TE: EE at ℓ≤30 is noise-dominated in Planck, so the ratio cannot
  carry a p~1% anomaly (failed check 3).
- T x phi lack-of-correlation (S1/2): reconstruction noise at L<10 swamps the ISW-lensing term
  (failed check 3).
- DESI DR1 ELG x ACT lensing: not found, but σ(S8)≈0.04-0.05 at z≈1.3. Magnification
  (α poorly known for ELGs, about 15% of the signal) is as large as that error, and unWISE x ACT
  already gives about 0.02 there (failed check 3).
- ACT DR6 aberration/Doppler: overlaps "sky dipoles" on the skip list.

---

## 1. Do SZ cluster masses really fall further below true masses at z>1? CMB-lensing masses of ACT DR6 clusters

**Question.** Does the mass bias of ACT SZ clusters, 1-b = M_SZ/M_true, keep falling with
redshift beyond z≈0.8, where shear calibration stops? That decline is what makes the cluster
counts at high z disagree.

**Shape.** Both a disagreement and a regime extension.
- (a) ACT DR5 x DES Y3 (2512.18935): 1-b falls with z at 99.95% confidence (power-law exponent
  -1.8 +0.5/-0.6; high-z bin 0.58, low-z bin 0.81). HSC Y3 (PRD 110, 103006, 2024) finds
  M_WL/M_SZ = 1.3±0.2, 1.6±0.2 and 1.6±0.3 up to z≈1.18. KiDS (2304.10219) finds no mass trend
  and has no z split. Hydrostatic-bias models usually expect weak evolution.
- (b) The step in dispute is the shear calibration at high z: few sources lie behind z≈1
  clusters, and their photo-z decides the result. 2608.00164 notes that WL calibrations reach
  only z≲0.8. CMB lensing removes that step, because its source plane is fixed at z=1100.
- Prior CMB-cluster lensing covers other samples only: SPT-SZ 513 clusters at 3.1σ (1412.7521;
  2111.07491), SPT-3G on DES redMaPPer at about 8% (2404.02153), and ACT on 677 MaDCoWS
  infrared-selected clusters at 4.2σ (2009.07772). No CMB-lensing mass exists for ACT
  SZ-selected clusters, and none for any SZ sample split at z>1.

**Public data.**
- ACT DR6 SZ cluster catalogue (2507.21459): 10,040 clusters, 3,747 in the Legacy sample,
  1,180 at z>1, median z=0.58, median M500c=2.8e14 Msun. Catalogue and selection function:
  https://lambda.gsfc.nasa.gov/product/act/actadv_prod_table.html (Nemo: https://nemo-sz.readthedocs.io)
- ACT DR6.02 f090/f150 maps: https://lambda.gsfc.nasa.gov/product/act/act_dr6.02/act_dr6.02_maps_prod_table.html
  (paper 2503.14451). mnms noise-simulation code and models: 2303.04180, at NERSC dr6.02.
- Planck PR4 frequency maps and a tSZ+CIB-deprojected gradient map: https://pla.esac.esa.int

**Method.**
1. Use a gradient-inversion / maximum-likelihood CMB-cluster lensing estimator (1705.00411;
   foreground-immune variant 2305.06326). The gradient leg comes from a Planck+ACT ILC with tSZ
   and CIB deprojected, filtered to ℓ<2000; the small-scale leg is the ACT DR6 f090+f150 coadd.
2. Stack it on DR6 clusters in z bins (z<0.5, 0.5-0.8, 0.8-1.0, >1.0) and fit an NFW κ profile
   with a c(M) prior, giving M_lens/M_SZ per bin (M_SZ is the catalogue's uncalibrated UPP mass).
3. Errors come from stacks at 10^3 random positions with the same mask and depth. Nulls:
   rotated-gradient (curl) stacks, and f090-only vs f150-only.
4. Calibrate bias with injected NFW lensing plus tSZ/CIB profiles in lensed-CMB sims with
   mnms noise. Fit 1-b(z) = (1-b0)((1+z)/1.5)^ζ and compare ζ with 2512.18935.

**Overlap.** Complete by construction: the clusters were detected in the same ACT DR6 maps. There
are 10,040 clusters over 16,293 deg² (12,590 deg² clean), 1,180 of them at z>1. Planck covers
the full sky for the gradient leg.

**Feasibility.**
- **Scaling.** The anchor is MaDCoWS: 677 clusters at z≈1 with M500c≈1.7e14 gave 4.2σ. Scaling
  by mass (about 2.3e14 for DR6 clusters at z>1), by √N (1,180) and by about 1.3 for DR6 depth
  gives about 10σ. A cost of about 1.3 for foreground cleaning brings that to 7σ, so 13-20% on
  M_lens/M_SZ at z>1. The full sample reaches about 4-9%. As a cross-check, SPT-3G reached 8%
  on 8,865 lower-mass clusters with about 3x lower noise.
- **What could erase it.**
  - tSZ in the lensing estimator: its central signal is far larger than the lensing dipole. A
    tSZ-free gradient leg removes the leading-order bias, and an injection test checks the rest.
  - CIB from dusty members at z>1, at 150 GHz: deproject it in the gradient leg. In the
    small-scale leg it adds variance, not first-order bias.
  - Miscentering: under 0.5' for SZ peaks.
  - The assumed c(M): a few-% effect on the mass.
  - SZ selection coupling to the lensing dipole: odd in orientation, so small. Test it with
    injection sims passed through Nemo selection.
  None of these is comparable to the predicted evolution.

**Expected precision and whether it decides.** At z>1, M_true/M_SZ is about 1.35 with no
evolution. The 2512.18935 power law extrapolated to z≈1.2 gives about 2.5-2.7, and HSC gives
1.6±0.3. With 13-20% errors the measurement tells no-evolution from the power-law extrapolation
apart at about 3-4σ, but HSC-like values from no-evolution at only about 1σ. A second result
falls out: the low-z bins give an independent check of 1-b≈0.75-0.8 at about 8-10% precision.

**Risk.** ACT has not published a CMB-lensing calibration of the DR6 catalogue, but likely
plans one, so scoop risk is moderate to high. 2602.03917 already fits DR6 counts and needs a
"significant residual mass bias", so the calibration is wanted.

**Quick-kill queries.**
- WebSearch "ACT DR6 Sunyaev-Zel'dovich cluster number counts cosmological constraints S8 2026":
  catalogue and fNL paper only.
- arXiv API all:ACT AND cluster AND DR6 AND (counts OR abundance): 2602.12238 (review) and
  2602.03917 (fNL, WL-calibrated); no S8 and no CMB-lensing masses.
- WebSearch "CMB lensing mass calibration SZ-selected clusters ACT DR6 stacked ... hydrostatic":
  WL only (2512.18935, 2304.10219, 2606.17984 Planck x HSC) plus a 2018 redMaPPer CMB-lensing
  study.
- arXiv API abs:"cluster lensing" AND CMB AND SZ: SPT only (1412.7521, 2111.07491), plus
  forecasts.
- arXiv API abs:"CMB-cluster lensing": 2404.02153 (SPT-3G x DES redMaPPer), 2401.08910 (method)
  and older work; no ACT SZ sample.
- arXiv search "cluster lensing" ACT DR6: no results.
- WebSearch "ACT SZ clusters stacked CMB lensing masses redshift evolution ... z>1 2025 2026":
  no CMB-lensing z>1 calibration of ACT clusters.

---

## 2. E-mode hemispherical power asymmetry: does half-mission coherence survive disjoint-detector (PR4 A/B) splits?

**Question.** Is the E-mode local-power dipole real CMB signal, aligned with the temperature
HPA axis? Or does its coherence come from systematics shared by the PR3 half-mission maps,
which use the same detectors?

**Shape.** A disagreement.
- (a) 2607.28774 (PR3 Commander, half-mission split, ℓ=15-40): the two dipole directions lie
  8.7° apart, p<0.33% from 300 sims, near the T axis. The coherence appears only in EE and
  vanishes without leakage correction.
- (b) On PR4, 2306.14880 finds the E-mode variance dipole at only 2.8% (SEVEM), strongly
  mask-dependent, with "no clear conclusion". 2504.05597 finds no significant clustering of
  E-mode power directions.
- The steps that differ: the release (PR3 vs PR4), the statistic, and the split. Half-missions
  share each detector's bandpass mismatch and ADC/transfer systematics; PR4 A/B detector sets
  do not. Running the same statistic on PR4 A/B removes both confounds.

**Public data.**
- PR4 SEVEM and Commander CMB maps with A/B detector splits: https://pla.esac.esa.int
- PR4 SEVEM: 600 A/B sims, with separate signal-only and noise-only sets. PR4 Commander: 400
  A/B Q/U sims. Both are at NERSC, via the PLA wiki:
  https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/Simulation_data
- PR3 half-mission maps, to reproduce the published number first.

**Method.**
1. Reproduce the 2607.28774 PR3 result first: local EE bandpowers in overlapping discs,
   pseudo-C_ℓ E→B deconvolution, dipole fit at ℓ=15-40.
2. Apply the identical pipeline to PR4 SEVEM A and B, then Commander A and B. Measure the
   angular separation between the A and B directions, each split's dipole amplitude with its
   interval, and the offset from the T axis.
3. Get p-values from the 600 SEVEM and 400 Commander A/B sims run through the same pipeline.
4. Nulls and checks: an (A-B)/2 noise map (should give no coherent dipole), three masks,
   ℓ=15-25 vs 25-40, and the cross-release direction PR3 vs PR4.

**Overlap.** Same sky. The PR4 common polarisation mask gives f_sky≈0.7-0.78. The sims
outnumber PR3's 300 by 1.3-2x, so p-values down to about 0.2% are resolvable.

**Feasibility.**
- **Statistical power.** Noise in each A/B split is close to that of a PR3 half-mission (half
  the detectors for the full time, with NPIPE's lower noise), so the power matches the PR3
  test. If the PR3 coherence is CMB sky, the sky term is common to A and B and the separation
  stays small (≲20°). If it is a detector-shared systematic, A and B decorrelate and the
  separation looks isotropic (median about 60-90°).
- **What could erase it.**
  - NPIPE's low-ℓ polarisation transfer function: it suppresses mainly ℓ<10, below the band,
    and the sims include it.
  - Residual Galactic foregrounds: these are common to A and B, so they could fake coherence.
    Compare SEVEM with Commander and vary the masks.
  - Leakage-correction artefacts: these are identical in the data and the sims, so the
    sim-calibrated p-value absorbs them.
- **Signal-to-noise.** Comparable to the PR3 run, which reached p<0.33%. The outcome splits
  cleanly either way.

**Expected precision and whether it decides.** A/B separation with a p-value at the 0.2-1%
level, plus per-split dipole amplitudes. If it stays under about 15° near the T axis with
p≲1% in both SEVEM and Commander, the E-mode HPA survives the removal of shared detectors: a
CMB-level signal and a sharp LiteBIRD target. If it is consistent with isotropy while PR3
reproduces, the PR3 coherence is detector-shared. That decides between 2607.28774 and the PR4
nulls.

**Overlap with our skip list.** This is a CMB power-modulation dipole, not a count or velocity
dipole. Check this reading against "sky dipoles" before taking it on.

**Quick-kill queries.**
- WebSearch "hemispherical power asymmetry E-mode polarization Planck PR4 NPIPE split detector
  2026": 2306.14880, 2607.28774 (PR3 half-mission), 2601.13830 and 2411.15786 (T only); no PR4
  split test.
- arXiv API abs:asymmetry AND polarization AND (PR4 OR NPIPE): 2306.14880 and 2504.05597; no
  A/B coherence.
- arXiv search "dipolar modulation E-mode polarization" (newest first): only 2607.28774 in
  2024-26, and no follow-up.
- WebSearch on PR4 A/B sims: confirmed SEVEM 600 and Commander 400 A/B sims.
