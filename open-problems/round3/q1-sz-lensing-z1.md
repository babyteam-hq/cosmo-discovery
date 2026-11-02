# Proposed test: Do SZ cluster masses really fall further below true masses at z>1? CMB-lensing masses of ACT DR6 clusters

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
