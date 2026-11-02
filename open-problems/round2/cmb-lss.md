# Round 2: CMB and LSS growth, quick-kill survivors

Checked 2026-11-01. The arXiv API (export.arxiv.org) returned HTTP 429 after the
first four calls, so the remaining arXiv checks used arXiv's own advanced
abstract search (arxiv.org/search/advanced, 2018-2026, newest first), which
covers the same index. Query counts below include both kinds plus WebSearch
(extended).

**Killed in the quick kill:** Lyman-alpha forest x CMB lensing on DESI DR1
(done: 2405.14988, DESI Y1 quasars x Planck PR3 bispectrum; also 2609.36949,
forest weak lensing validated on ACT DR6). CMB lensing x DESI ELGs at z 1.1-1.6
was dropped without a full check: DESI DR1 ELG x lensing is too obvious a
collaboration target, and Quaia already covers z~1.5. Supervoid ISW on DESI DR1
was not taken further (DESI LRG void lensing exists, and void work is on the
skip list).

---

## 1. CMB lensing around z~5 Lyman-break galaxies (HSC r-dropouts) with ACT DR6

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

## 2. Thermal energy around quasar hosts: feedback excess or not?

- **Question:** What is the mean thermal energy (integrated Compton-y) of
  quasar host haloes at z~1-3: ~1e62 erg (feedback-dominated) or ~5e60 erg
  (consistent with gravitational heating)?
- **Shape:** disagreement.
  - High: Ruan et al. 2015, arXiv:1502.01723. Planck y stack on 26.7k SDSS
    quasars, E_th ~ 1.1e62 erg at z~1.5, read as quasar feedback.
  - Low: Soergel et al. 2017, arXiv:1612.06296 (~5e60 erg, 1.6 sigma) and
    Crichton et al. 2016, arXiv:1510.05656 (ACT and Herschel, ~6.2e60 erg,
    3-4 sigma, dust-model dependent).
  - Related: Hall et al. 2019, arXiv:1907.11731 (ACT, Herschel and VLA, 3.8
    sigma at z>1.9, dust-model dependent); Verdier et al. 2016,
    arXiv:1509.07306 (tSZ only at z 2.5-4).

  The disagreement is driven by the dust model and by the Planck beam, which
  mixes in clustered haloes.
- **Data:**
  - DESI DR1 quasar catalogue and LSS randoms, https://data.desi.lbl.gov/doc/releases/dr1/
  - ACT DR6.02 needlet-ILC y maps with CIB-deprojected variants (several beta
    and T values, and moment deprojection):
    https://lambda.gsfc.nasa.gov/product/act/act_dr6.02/act_dr6.02_nilc_prod_table.html
  - ACT DR6 + Planck y map:
    https://lambda.gsfc.nasa.gov/product/act/actadv_dr6_compton_maps_info.html
  - Planck 353/545/857 GHz maps (PLA) for the dust SED.
- **Method:**
  1. Stack the standard and CIB-deprojected ACT y maps on DESI DR1 quasars
     (radio-loud ones cut by a FIRST/VLASS match), in z and L_bol bins, with
     compensated apertures of 1-3 arcmin.
  2. Errors from random-position stacks and a bootstrap over sky patches;
     nulls on randoms and on the difference between deprojection variants.
  3. Fit tSZ, dust (modified blackbody), and a 2-halo term fixed by quasar
     clustering bias.
  4. Convert Y to E_th and compare with the ~1e62 and ~5e60 erg answers, and
     with the gravitational expectation for M_h ~ 10^12.5 Msun.
- **Expected precision:** Hall et al. got 3.8 sigma from ~110k SDSS quasars on
  older ACT maps. DESI DR1 gives several times more quasars overlapping ACT
  DR6's ~1/3 of the sky, at lower noise, so the statistical error should be
  ~10 sigma. The two answers differ by ~20x, so they're easily told apart if
  the result is stable across CIB-deprojection variants. Dust is the limiting
  systematic, and the variants measure it directly.
- **Quick kill (5 queries):**
  - WebSearch, "tSZ stacking quasars z>3 DESI ACT DR6": nothing on quasars.
    The ACT x DESI tSZ stacks are of LRGs, BGS and clusters (2502.08850,
    2609.08939).
  - WebSearch, Hall/Crichton/Ruan disagreement: confirms the dust-driven
    split and finds no resolution.
  - WebSearch, "quasar tSZ Compton-y stacking eBOSS/DESI 2021-2025": only the
    2015-2019 papers.
  - arXiv search, abstract "quasar Sunyaev thermal" 2018-2026: Hall 2019, an
    ALMA stack of 10 quasars (2310.09445) and simulations, with no ACT DR6 or
    DESI quasar stack.
  - arXiv search, abstract "quasars SZ stack" 2018-2026: none.

## 3. Bias-weighted gas pressure <b P_e>(z) beyond z = 1.2

- **Question:** How does the cosmic mean thermal pressure of gas, <b P_e>,
  evolve from z = 1 to z ~ 3, where halo-model predictions and simulations
  (Magneticum, 2105.15043) begin to diverge?
- **Shape:** regime extension (redshift).
  - Chiang et al. 2020, arXiv:2006.14650: Planck and SDSS, direct constraints
    to z~1.
  - DES Y3 x SPT/Planck, arXiv:2210.08633: up to z = 0.97, which they call the
    highest so far.
  - ACT x WISE-DESI galaxies, 2508.09514: up to z = 1.2, as a stacked CGM
    measurement, not <b P_e>.
  - Simulations: FLAMINGO 2605.11083, z <= 1.
- **Data:** the same y maps as candidate 2 (ACT DR6.02 NILC with CIB
  deprojection, and ACT+Planck), plus Planck PR4 y. DESI DR1 ELGs (z 1.1-1.6)
  and quasars (z 0.8-3.5) with randoms:
  https://data.desi.lbl.gov/doc/releases/dr1/
- **Method:**
  1. Cross-power C_ell^{yg} for DESI ELG and quasar slices of width 0.2-0.4 in
     z, on large scales (ell < 1000), where the 2-halo term dominates.
  2. Measure the tracer bias b_g from DESI's own auto-clustering, or from
     CMB-lensing cross-correlation.
  3. Dust: fit with CIB-deprojected y variants and treat the spread between
     them as a systematic.
  4. Errors from Gaussian simulations with the measured spectra plus jackknife;
     nulls on randoms and on the deprojection differences. Report <b P_e>
     against z from 1.0 to 3.
- **Expected precision:** DES Y3 reached ~10-20% per bin at z<1. At z 1-2, the
  DESI ELG and quasar densities with ACT DR6 noise should give ~20-30% per bin,
  and ~40% at z 2-3, where the quasars are sparse. That's enough to test the
  halo-model trend against the Magneticum divergence, which is tens of percent
  by z~2. The CIB is the dominant risk at z>1.5.
- **Quick kill (3 queries):**
  - arXiv search, abstract "thermal history Sunyaev pressure redshift"
    2019-2026: 4 hits, none measuring above z~1.
  - WebSearch, "bPe z>1 quasars ELG Compton-y tomography": highest is DES Y3
    at z=0.97.
  - arXiv search, abstract "Compton-y cross-correlation redshift pressure"
    2022-2026: 3 hits, none above z=1.2.
  - Risk: ACT and DESI members are probably working on this (several ACT x
    DESI tSZ papers in 2025-26). Chiang 2020's sky spectra reach z~3 but give
    no pressure constraint above z~1.
