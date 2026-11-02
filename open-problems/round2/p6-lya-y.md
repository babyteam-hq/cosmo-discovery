# Proposed test: Lyman-alpha forest x Compton-y at z = 2-3.5 (thermal pressure of the high-z IGM and its halos)

- **Question.** What is the bias-weighted thermal pressure of the gas traced by the Lya forest at z ~ 2-3.5, measured
  as the cross-correlation of forest flux fluctuations with a Compton-y map?
- **Shape.** Regime extension. The only measurement is Croft, Banday & Hernquist (astro-ph/0512380): 3,488 SDSS DR3
  forests x WMAP-1, a "preliminary" upper limit Omega_b,ionised < 0.8. Bias-weighted pressure <bP_e> tomography with
  galaxies/quasars detects only to z ~ 1 and gives limits to z ~ 3 (Chiang et al. 2020, ApJ 903, abb403; DES Y3
  2210.08633; the FLAMINGO comparison 2605.11083 covers z = 0.1-1). Public data now reach ~100x more forests and
  arcminute y maps.
- **Data.** DESI DR1 Lya forest catalogue and quasars (https://data.desi.lbl.gov/doc/releases/dr1/); eBOSS DR16 Lya
  (https://www.sdss4.org/dr17/); ACT DR6 Compton-y maps including CIB-deprojected variants (Coulton et al. 2307.01258,
  https://lambda.gsfc.nasa.gov/product/act/); Planck PR4 NILC y (2305.10193; Planck legacy archive https://pla.esac.esa.int).
  Websky simulated y + CIB skies for contamination (https://mocks.cita.utoronto.ca/data/websky/).
- **Method.**
  1. delta_F per pixel from DESI DR1 continuum fits, in z bins 2.0-2.5, 2.5-3.0, 3.0-3.6.
  2. Cross-correlate with y (pixel stack vs angular separation, 1'-60'), on standard and CIB-deprojected maps.
  3. Errors and nulls: rotated/flipped y maps, forests shuffled between sightlines, ACT/Planck difference maps.
  4. Model as b_F x <b P_e>(z), with b_F from DESI's own Lya auto-correlation; compare with hydro predictions.
- **Precision.** Not forecast yet. Croft et al. sat within ~10x of the predicted signal with WMAP-1 and 3.5k forests;
  DESI DR1 has ~100x the forests and ACT DR6 far lower small-scale noise, so a few-sigma detection, or a limit below
  the hydro predictions, is plausible. A first measurement, not a tie-breaker. Main risk: CIB at z ~ 2 traces the same
  structure; the deprojected maps and the frequency dependence of the stack are the check.
- **Quick kill (5 queries).** WebSearch (extended) "Lyman-alpha forest cross-correlation thermal SZ Compton-y map":
  only astro-ph/0512380 and a simulation paper (2006.09742). WebSearch (extended) for <bP_e> at z > 3: none beyond
  Chiang 2020's limits. arXiv API: abs:Lyman AND abs:Sunyaev AND abs:forest -> only 0512380 measures it;
  abs:Lyman AND abs:forest AND (Compton-y OR "thermal SZ" OR tSZ OR y-map) -> only 0512380;
  abs:"Lyman-alpha" AND abs:Sunyaev AND abs:cross -> only 0512380.
