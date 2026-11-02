# Proposed test: CMB lensing around z~5 Lyman-break galaxies (HSC r-dropouts) with ACT DR6

- **Question:** Does the CMB-lensing x galaxy cross-power of z~5 r-dropouts
  match LCDM's sigma8(z~5) and the clustering-derived halo mass, i.e. is matter
  clustering at z~5 as LCDM predicts?
- **Shape:** regime extension (redshift). Prior work stops at z~4:
  Miyatake et al. 2022, arXiv:2103.15862 (1.5M HSC g-dropouts x Planck lensing,
  ~5 sigma); UNIONS LBGs at z~2.5 x Planck (2607.09943); forecasts only for
  z=3-4 (1904.13378) and fNL (2511.22243).
- **Data:**
  - GOLDRUSH HSC dropout catalogues z~4-7 (Harikane et al. 2022,
    arXiv:2108.01090), Zenodo https://zenodo.org/record/5512721 ("masked_"
    files, with masks).
  - ACT DR6 lensing maps plus their noise simulations:
    https://lambda.gsfc.nasa.gov/product/act/actadv_dr6_lensing_maps_info.html
  - Planck PR4 lensing with its simulations, for the HSC fields outside the
    ACT footprint and for a joint estimate.
- **Method:**
  1. Build r-dropout (z~5) and g-dropout (z~4, to recover the Miyatake result)
     overdensity maps on the HSC-Wide footprint, with random masks.
  2. Measure C_ell^{kappa g} with a pseudo-C_ell (NaMaster) on ACT DR6 and
     Planck PR4 kappa, for ell between 100 and 1500.
  3. Get errors from the released lensing simulations crossed with the real
     galaxies, and run nulls on curl-mode and rotated catalogues.
  4. Fit A_lens x b, with b fixed from the GOLDRUSH angular clustering, so
     that the result is sigma8(z~5)/sigma8_LCDM; model the low-z interloper
     fraction with the published contamination estimates.
- **Expected precision:** Miyatake et al. got ~5 sigma with 1.5M z~4 galaxies
  on Planck. A few hundred thousand r-dropouts with ACT DR6 noise (lower per
  mode than Planck at ell~1000 over the equatorial HSC fields) should give
  roughly 3-5 sigma, so the amplitude is known to ~20-30%. That's a first
  matter-clustering check at z~5 (it can't resolve percent-level S8 shifts).
  The catalogue counts per dropout band still need checking: this estimate
  assumes a few hundred thousand r-dropouts.
- **Quick kill (6 queries):**
  - WebSearch, "CMB lensing x LBG z~5/z~6 HSC r-/i-dropout": only 2103.15862
    (z~4).
  - WebSearch, "LBG CMB lensing z~5 ACT DR6": none.
  - arXiv API, all:dropout AND "CMB lensing": 2511.22243 and 1904.13378,
    both forecasts.
  - arXiv search, abstract "Lyman break" "CMB lensing": 5 hits, none at z>=4.5.
  - arXiv search, abstract "CMB lensing" "Hyper Suprime-Cam": 6 hits, none at
    z>=4.5.
  - arXiv search, abstract "dropout lensing z~5": none.
