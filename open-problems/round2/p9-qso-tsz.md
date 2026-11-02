# Proposed test: Thermal energy around quasar hosts: feedback excess or not?

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
