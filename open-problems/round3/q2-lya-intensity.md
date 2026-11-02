# Proposed test: Mean Lyα surface brightness at z 2-3.5 from DESI DR1 Lyα forest x fibre-residual cross-correlation

**Question.** Is the diffuse Lyα emission seen in quasar-centred stacks cosmic, from star-forming galaxies,
or is it powered by the quasars themselves? Measure b_I<I>, the bias-weighted mean Lyα intensity, with an
interval at z ≈ 2.2, 2.6 and 3.1. The tracer is the Lyα forest, which has no quasar-specific term.

**Shape: disagreement.** Three published answers disagree:
- Quasar x Lyα-intensity cross-correlations give a large signal. Croft et al. 2016 (arXiv:1504.04088), BOSS.
  Croft et al. 2018 (1806.06050) is about 50% lower and calls it quasar-dominated within 15 Mpc/h. Lin et al.
  2022 (2207.10682), eBOSS DR16, finds a luminosity density of 6.6 (+3.3/-3.1) x 10^40 erg/s/cMpc^3, 8-9x the
  narrowband LAE estimates, which they put down to star-forming galaxies, not quasars.
- HETDEX LAE x intensity (2510.11427) measures b_g b_I <I> = 6.7±3.1, 11.7±1.4 and 8.3±1.5 x 10^-22
  erg/s/cm^2/arcsec^2/Å at z = 2.1, 2.6 and 3.2. That is "significantly smaller" than all the quasar-cross
  values. The authors suggest the quasar crosses are dominated by quasar-affected small-scale physics.
- Lyα forest x intensity on BOSS DR12 (Croft 2018) gives a 95% upper limit μ < 1.2 x 10^-22 in the same
  units. The HETDEX mid-z bin sits 2.6σ above that limit.

**The disputed step** is whether the quasar-centred signal includes emission sourced by the quasar
(fluorescence, or leakage from neighbouring quasars) on top of the galaxy-sourced field. The forest cross
contains only the galaxy-sourced term, b_F b_I <I>. Measuring it, and the quasar cross on the same fibres,
removes that step. BOSS couldn't, because it detected nothing in the forest cross.

**Public data.**
- DESI DR1 coadded spectra and redrock fits (template coefficients) for LRG, ELG and sky fibres:
  https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/ (docs: https://data.desi.lbl.gov/doc/releases/dr1/)
- DESI DR1 quasar catalogue and Lyα forest spectra (the same release; deltas rebuilt with the public picca
  code, https://github.com/igmhub/picca)
- DESI DR1 LSS catalogues: https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/
- Comparison values: HETDEX 2510.11427 (table), Croft 2018, Lin 2022

**Method.**
1. Residual spectra are coadd minus redrock best-fit model, for DR1 LRG and ELG fibres at 3600-5550 Å
   (Lyα z 1.96-3.56). Median-filter out the broadband, and mask sky lines and each galaxy's own emission lines.
2. Cross-correlate the residual pixels with DR1 forest δ_F, picca-style, in (r_⊥, r_∥) to 30 Mpc/h.
   Use the forest distortion matrix, and drop every pair whose two fibres share an exposure, which kills
   sky-model and CCD-crosstalk terms at r_∥ = 0. Do the same with DR1 quasars, with the ±5-fibre exclusion.
3. Fit b_F(1+β_F μ²) b_I <I> ξ_lin, with b_F and β_F fixed to DESI DR1 Lyα values. Fit the quasar cross with
   an added r^-2 quasar-sourced term.
4. Nulls: sky fibres x forest (must be zero), and forest x residuals at 5600-9800 Å (no Lyα at the forest's
   redshifts). Errors come from jackknife over HEALPix, checked on the public DESI Lyα mocks (LyaCoLoRe; see
   2607.27412) with injected intensity.

**Overlap.** It is one survey, so overlap is complete over the DR1 dark footprint. Counts: about 420k Lyα
forests (DR1 Lyα BAO, 2404.03001), 2.14M LRG and 2.43M ELG redshifts (DR1 BAO samples), and about 0.7M quasars
at z > 1.96 as tracers. BOSS DR12 had about 1.2-1.4M LRG fibres and about 160k forests, so DESI DR1 has about
3.4x the emission fibres and 2.7x the forests: about 9x the pair count. Fibre-pair loss from the
same-exposure cut is estimated at 30-50% at r_⊥ < 15 Mpc/h (0.2° at z = 2.5), since DR1 dark tiles cover
most of the footprint 2-4 times. **Verify** this on the tile files first.

**Feasibility.**
- **Erasers.**
  - Correlated sky-subtraction residuals (DESI 2404.03003) and crosstalk land exactly at r_∥ = 0, where
    the signal peaks. The same-exposure cut removes them, and the sky-fibre null tests that it did.
  - The forest continuum fit removes line-of-sight modes, which the distortion matrix models.
  - The fibre area is smaller: 1.77 vs 3.14 arcsec², about 1.3-1.8x the per-fibre surface-brightness
    noise variance of BOSS.
  - Interpreting the result needs b_F and β_F, which DESI measures to a few percent. Comparing with HETDEX
    needs b_g of the LAEs (1.7-2.0), a 15% systematic, which is stated.
- **Precision (scaled, not simulated).**
  - BOSS's 95% limit of 1.2e-22 implies σ ≈ 0.6e-22. DESI's σ ≈ 0.6e-22 x sqrt(1.5/9) x 1/sqrt(0.6)
    (pair loss) ≈ 0.3e-22, about 2x better than BOSS.
  - A HETDEX-level intensity would be seen at about 4-6σ. A quasar-cross-level intensity (several times
    larger) would be seen at >10σ.
  - A null at about 2σ = 0.6e-22 excludes the HETDEX z = 2.6 value at about 4σ.

**Tells the answers apart:** yes. A forest detection at the HETDEX level means the quasar crosses are
inflated by quasar-sourced emission, and fitting the quasar-specific term on the same fibres sizes that
inflation. A null below 0.6e-22 means HETDEX's mid-z intensity is too high.

**Quick-kill queries (7):**
1. WebSearch "Lyman-alpha intensity mapping DESI quasar cross-correlation fiber residual spectra": BOSS/eBOSS
   only, no DESI.
2. arXiv API `all:"intensity mapping" AND Lyman AND quasar`: 1504.04088, 1806.06050, 2207.10682, 2510.11427;
   no DESI.
3. arXiv search "Lyman-alpha intensity mapping DESI": only 2605.26518 and 2504.03848, neither a measurement.
4. arXiv API `DESI AND "Lyman-alpha emission" AND residual`: 0 results.
5. arXiv API `abs:"Lyman-alpha forest" AND intensity AND emission`: no forest x Lyα-emission measurement
   since 2018.
6. WebSearch on the quasar-luminosity dependence and the BOSS/eBOSS vs HETDEX discrepancy (2025-26): no
   resolution paper. BOSS found no quasar-luminosity dependence.
7. Author check (Renard, Lyman): unrelated (PAUS quasar LF, 2407.19020).
