# Round 3: high-z galaxies, reionisation and black holes, quick-kill survivors

Checked 2026-11-01. Searches used WebSearch (extended), the arXiv API (export.arxiv.org, which answered this
session) and arXiv's own search. The precision numbers are back-of-envelope scalings from the nearest
published measurement, and say so. Pre-registration should redo them on mocks before any compute.

**Result: two survivors.** One is a disagreement about Lyα intensity at z 2-3.5. The other is a
disagreement about how much quasars cluster below 100 kpc.

## 1. Mean Lyα surface brightness at z 2-3.5 from DESI DR1 Lyα forest x fibre-residual cross-correlation

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

## 2. Quasar clustering excess at 20-150 kpc from DESI DR1, against its own Mpc-scale clustering

**Question.** What is the ratio of the projected quasar correlation w_p(r_p) at r_p ≈ 20-150 kpc (proper) to
the extrapolation of the same sample's 0.5-20 Mpc/h power law? Measure it in 3-4 redshift bins over
0.8 < z < 3, with an interval.

**Shape: disagreement.**
- Hennawi et al. 2006 (astro-ph/0504535), spectroscopic binaries: an order-of-magnitude excess at
  < 40 h^-1 kpc, about 30x near 10 h^-1 kpc.
- Myers et al. 2007 (astro-ph/0612191), photometric pairs: only a slight excess, at most 4.3x. They argue
  that surveys biased to z > 2 inflate it.
- Eftekharzadeh et al. 2017 (MNRAS 468, 77; 47 binaries): A = 24.1 ± 3.6 at 26.6 h^-1 kpc and a steep
  γ = 1.97 ± 0.03, against 1.5-1.8 on Mpc scales.

The disputed step is the completeness and selection of heterogeneous, follow-up-built pair samples. DESI DR1
has a homogeneous spectroscopic quasar sample with fibre-assignment probabilities (PIP bitweights from 128
alt-MTL realisations), and measures the large-scale amplitude from the same objects, which removes it.
Recent work counts pairs but doesn't measure clustering:
- Jing et al. 2505.03103: 1,020 confirmed DESI DR1 pairs within 110 kpc, pair fraction 6.2 x 10^-4.
- 2504.06415: SDSS x Gaia at 3-30 kpc, a pair fraction compared with simulations, not with the clustering
  extrapolation.
- 2603.03042: Gaia pair search IV.

**Public data.**
- DESI DR1 LSS v1.5pip catalogues (BITWEIGHTS, angular upweights), QSO data and randoms:
  https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/ (the "pip" directories; method
  2411.12025, Appendix C of 2411.12020)
- DESI DR1 QSO target catalogue (for angular pair counts of targets) and spectra:
  https://data.desi.lbl.gov/public/dr1/
- Legacy Surveys DR9/DR10 images, to vet lensed pairs: https://www.legacysurvey.org/
- The Jing et al. 2505.03103 pair list, as a cross-check

**Method.**
1. Take DR1 QSOs at 0.8 < z < 3.0 with |Δv| < 2000 km/s (π_max ≈ 20 h^-1 Mpc). Compute w_p(r_p) from 15 kpc to
   20 Mpc/h with PIP + angular-upweighted pair counts. Restrict to θ > 2.5" (Legacy Survey deblending) and set
   the smallest bin by that.
2. Remove lensed pairs: flag matching spectra with |Δv| < 100 km/s and θ < 5", and carry the uncertainty as a
   systematic, using the 58/1,220 lens-candidate rate from 2505.03103.
3. Fit r_0 and γ on 0.5-20 Mpc/h. Report the excess E(r_p) = w_p,obs / w_p,PL in each redshift bin.
4. Errors come from jackknife over HEALPix, checked on DR1 alt-MTL realisations run on mock quasar
   catalogues (the AbacusSummit DR1 mocks) with no small-scale excess injected. That null also tests whether
   the fibre-assignment correction leaves a bias.

**Overlap and counts.** One survey: about 1.6M DR1 quasars (0.86M in the 0.8 < z < 2.1 BAO sample). There are
1,020 confirmed plus 142 candidate pairs at < 110 kpc (2505.03103).
- Expected pairs at r_p < 100 h^-1 kpc under a power-law extrapolation (r_0 ≈ 6 h^-1 Mpc, γ = 1.8,
  n ≈ 2 x 10^-5 h^3 Mpc^-3, π_max 20) ≈ 20 random pairs x ξ̄ ≈ 30, so about 600.
- That is the same order as what's observed, so the ratio is measurable and not swamped. Hennawi's 10x would
  mean about 6,000.

**Feasibility.**
- **Erasers.**
  - Fibre collisions. DESI positioners can't place two fibres a few arcsec apart on one tile, so close
    pairs come only from multi-pass coverage. Pairs with zero assignment probability are recovered by angular
    upweighting from the target catalogue, which is the method's stated limit (2411.12025).
  - Single-pass areas give zero-probability pairs. Restrict to the area covered by 2+ dark tiles and quote
    the area kept.
  - Lensed-pair contamination (about 5% of pairs).
  - Target-selection blending below about 2".
  - Assumed coefficient: γ, from the same sample's Mpc fit. Carry its error through.
- **Precision.** Poisson on about 1,000 pairs gives about 3%. Completeness modelling is estimated at 10-20%
  (from the null on the alt-MTL mocks). So E(50-100 kpc) comes to about ±15-25% per redshift bin.

**Tells the answers apart:** yes. E ≈ 10 (Hennawi), E ≤ 4 (Myers) and E ≈ 1 differ by much more than ±25%.
The redshift bins test Eftekharzadeh's claim that the kpc-scale amplitude tracks the Mpc-scale evolution.

**Quick-kill queries (6):**
1. WebSearch "small-scale quasar clustering excess binary quasars DESI pairs kpc scales 2025": the SDSS-era
   papers, the 2505.03103 pair catalogue and a DESI-LS lensed/dual discovery paper. No DESI clustering
   measurement.
2. arXiv API `quasar AND pairs AND DESI AND clustering`: only 2512.18281 (transverse proximity, not
   quasar clustering).
3. arXiv API `abs:quasar AND ("small-scale clustering" OR "binary quasars" OR "quasar pairs")`, newest 50:
   pair catalogues, Gaia searches (2504.17777, 2603.03042, 2607.00093) and the SDSS x Gaia pair fraction
   (2504.06415, no clustering extrapolation). Also Trinity VII (2506.05612, model vs compiled data). No
   DESI w_p at < 150 kpc.
4. Abstract checks of 2505.03103 and 2504.06415: neither computes ξ or the excess over the extrapolation.

## Killed or dropped

| candidate | shape | why dropped |
|---|---|---|
| Halo mass vs BH mass from weak lensing of AGN | disagreement on the M_BH-M_halo link | Done: UNIONS, 36k AGN in M_BH bins (2402.10740) |
| DESI Lyα P1D extended to z 4.4-5 (reionisation relics / UVB fluctuations) | extension of DR1 P1D (2505.07974, 2505.09493; z 2.2-4.2) | Public DR1 has few z > 4.5 quasars (why DR1 stops at 4.2). Relic signal is degenerate with T0 and the mean flux at DESI resolution, so feasibility fails |
| Mean free path at z 4.6-5.3 from DESI DR1 stacks | extension of 2411.15838 (3.2-4.6) | That paper already used all of DR1 (Y1); beyond z = 4.6 needs DR2 (not public) |
| Periodic-quasar (SMBH binary) candidates tested with ZTF | disagreement on whether CRTS/PTF candidates are real | Crowded: 2206.11497 (ZTF, 17-yr baselines), 1803.05448, the PS1 MDS persistence tests |
| Patchy τ(n̂) from ACT DR6 | extension | Feasibility: expected limits far above model amplitudes (see `09-high-z.md`) |
| Galaxy x Lyα-forest cross at z ≈ 3.8 (HSC g-dropouts x DESI forest) | extension | Photo-z width (Δχ ≈ 500 Mpc/h) dilutes the line-of-sight forest signal about 50x, and no dispute to decide; dropped on feasibility |
