# Proposed test: Quasar clustering excess at 20-150 kpc from DESI DR1, against its own Mpc-scale clustering

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
