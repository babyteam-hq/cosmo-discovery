# Round 2: clusters/baryons and small-scale dark matter (quick-kill survivors)

Checked 2026-11-01. Candidates killed in the quick kill (prior work already answers them):
Faraday-rotation CGM magnetism around galaxies (SPICE-RACS z~1, 2605.16924; MIGHTEE 2308.11391; LOFAR 2302.06617);
the Milky Way's declining rotation curve tested by a third tracer (DESI Y3 K giants 2606.18849, streams, review 2608.10189);
Gaia wide binaries vs subhalos/MACHOs (2209.08100; Bootes I 2509.04555);
kSZ vs X-ray gas fractions (2509.10455 and 2512.02954 already combine them).

Search note: arXiv's API and search page (and Semantic Scholar) rate-limited this machine's curl (HTTP 429) after a
burst, so the arXiv queries counted below went through WebFetch on `export.arxiv.org/api/query`.

---

## 1. Lyman-alpha forest x Compton-y at z = 2-3.5 (thermal pressure of the high-z IGM and its halos)

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

## 2. tSZ around quasars at z = 3.5-5: is the quasar-stack excess feedback or host halo?

- **Question.** How much thermal energy surrounds quasars at z > 3.5, where host-halo gravitational heating is
  predicted to be small?
- **Shape.** Both. Extension: quasar tSZ stacks stop at z = 3.5 (Ruan et al. 1502.01723, Planck, mean z ~ 1.5;
  Crichton et al. 1510.05656 and Hall et al. 1907.11731, ACT, 0.3 < z < 3.5, detection at z > 1.9). Disagreement:
  Ruan/Crichton read the signal as quasar feedback; Cen & Safarzadeh (1509.07306) and Soergel et al. (1703.00920)
  say host-halo gas explains most of it. At z > 3.5 the halo term falls steeply while feedback energy need not, so
  the redshift dependence of the stacked thermal energy separates the two.
- **Data.** DESI DR1 quasars at z > 3.5 (https://data.desi.lbl.gov/doc/releases/dr1/) plus SDSS DR16Q; ACT DR6
  f090/f150/f220 maps and DR6 y maps (LAMBDA, above); Planck PR4 HFI maps for the dust SED; Websky halos and CIB for
  leakage nulls. Host halo mass from published quasar clustering (an input, not reanalysed).
- **Method.**
  1. Aperture-photometry stacks on DESI quasars in z bins to 5, with Hall et al.'s bins repeated at z < 3.5.
  2. Fit tSZ + modified-blackbody dust per bin with ACT DR6 + Planck HFI, as Hall et al. did.
  3. Errors from random-position stacks and bootstrap; nulls from Websky quasar-mass halos with no feedback.
  4. Compare E_thermal(z) with the halo-only expectation from the quasar halo mass.
- **Precision.** Hall et al. got ~3-4 sigma per bin at z ~ 2-3.5 from ~110k SDSS quasars on older ACT maps. DESI DR1
  has fewer quasars at z > 3.5 (order 10^4-10^5) but ACT DR6 is deeper: expect ~2-3 sigma, enough to tell a feedback
  energy of ~10^60 erg from a ~10x smaller halo-only expectation only if dust is controlled. Dust leakage is the risk.
- **Quick kill (5 queries).** WebSearch (extended) "quasar host tSZ stacking z>3 ACT DR6 DESI": nothing beyond
  z = 3.5; ACT x DESI tSZ work is on LRGs/BGS (2502.08850, 2609.08939). arXiv API: (Compton-y OR tSZ) AND quasars AND
  redshift -> highest is z <= 3.5 (1907.11731) and ALMA z 2.2-2.7 (2310.09445); quasar AND (thermal Sunyaev OR tSZ)
  AND stack -> nothing after 2015; quasars AND SZ AND (DESI OR "ACT DR6") -> kSZ only (2604.04867); quasars AND SZ AND
  feedback -> no z > 3.5 measurement.

## 3. Circumgalactic dust reddening around galaxies at z ~ 0.8-1.6

- **Question.** How much dust sits in the halos of galaxies at z ~ 1-1.5 (the reddening profile of background
  sources against impact parameter), and how does it compare with z ~ 0.3?
- **Shape.** Regime extension (redshift). Galaxy-selected statistical reddening exists only at low z: Menard et al.
  0902.4240 (SDSS, z ~ 0.3), Peek et al. 1411.3333 (z ~ 0.05), 2503.04098 (redMaGiC, WISExSuperCOSMOS), KiDS-DR5
  2609.34615 (z < 0.4), Rubin DP1 2606.24159 (4.6 deg^2, photo-z). At z > 1 only absorber-selected reddening exists
  (Mg II absorbers, e.g. 1603.06434 and DESI's Mg II extinction paper) or two protoclusters with JWST (2601.15961).
- **Data.** DESI DR1 ELGs (0.8 < z < 1.6) and LRGs to z ~ 1 as foregrounds; DESI DR1 quasars at z > 2 with spectra as
  backgrounds (per-object continuum slope) (https://data.desi.lbl.gov/doc/releases/dr1/); Legacy Surveys DR10
  photometry and photo-z for a larger photometric foreground sample (https://www.legacysurvey.org/dr10/); SFD/Planck
  Galactic dust maps for systematics.
- **Method.**
  1. For each background quasar, a reddening estimator: continuum slope against the DESI composite at its redshift,
     plus g-r, r-z, z-W1 colour residuals.
  2. Stack against foreground ELG/LRG impact parameter, 20 kpc-3 Mpc; flag pairs with Mg II at the foreground z.
  3. Nulls: foregrounds behind the quasar, random positions, stars as backgrounds; errors by jackknife.
  4. Fit A_V(r_perp) and compare with Menard's z ~ 0.3 profile; magnification from the same pairs as a cross-check.
- **Precision.** ~10^6 spectroscopic pairs within 1 Mpc, ~10^7 with photometric foregrounds; Menard got ~20 sigma
  from ~85k quasars x 24M galaxies, so a several-sigma profile at z ~ 1.2 looks reachable. Risk: DESI quasar targeting
  is colour-selected, so reddened quasars drop out; model it with the DESI selection.
- **Quick kill (5 queries).** WebSearch (extended) "circumgalactic dust reddening background quasars foreground
  galaxies z>1 DESI": none galaxy-selected at z > 1. arXiv API: abs:"circumgalactic dust" -> low z only;
  reddening AND circumgalactic AND dust AND quasars -> nothing at z > 1; reddening AND ("foreground galaxies" OR
  "around galaxies" OR circumgalactic) AND (quasars OR "background galaxies") -> low z or absorber-selected only;
  dust AND extinction AND (ELG OR z~1) AND halo AND quasars -> only 0902.4240 (z ~ 0.3).

## 4. Subhalo abundance from the density power spectra of several thin streams together

- **Question.** What relative subhalo abundance (10^6-10^9 Msun) do the density fluctuations of several thin streams
  imply jointly, and is each stream's power at its predicted epicyclic scale or broadband?
- **Shape.** Extension (sample) that bears on a disagreement. Density-power subhalo constraints exist only for GD-1
  and Pal 5 (Banik et al. 1911.02662; Banik & Bovy 1809.09640). Ibata et al. 2002.01488 say GD-1's power is
  epicyclic, no subhalos needed; Banik et al. reject that. Streams on different orbits predict different epicyclic
  scales but share one subhalo population, so a joint fit separates the readings. Phoenix power spectra exist only as
  a survey-systematics demonstration (2510.07511); the 14-stream paper 2608.19321 is a forecast.
- **Data.** DES DR2 (https://des.ncsa.illinois.edu/releases/dr2) and DELVE (https://datalab.noirlab.edu/delve/)
  photometry; Gaia DR3 proper motions; S5 DR1 radial velocities (https://s5collab.github.io/); Patrick et al.'s
  uniform DES density models of 13 streams (2206.04241); STREAMFINDER DR3 atlas (2311.17202). Simulation: galpy
  streampepperdf (public) for subhalo-perturbed mocks.
- **Method.**
  1. Members and 1D density for ATLAS-Aliqa Uma, Phoenix, Jhelum, Indus, Turranburra (+ GD-1) from matched filter and
     proper motions, depth-corrected with the 2510.07511 selection maps.
  2. Density power spectrum per stream; null = the same filter on off-stream fields (catches scanning/depth patterns).
  3. Forward-model each stream over subhalo abundance (streampepperdf) and with an epicycle-only model; ABC as Banik.
  4. Joint posterior on the relative abundance; per-stream test of power at the epicyclic scale against broadband.
- **Precision.** GD-1 alone gave abundance ~0.7 (about a factor of 2; < 0.9 at 95%). Four or five more streams of
  similar depth could tighten that toward ~30-50%, and the epicyclic-scale test could separate the readings at
  ~2-3 sigma if baryonic perturbers (bar, GMCs, Sagittarius for Jhelum) are modelled; they ruined Pal 5 as a probe.
- **Quick kill (4 queries).** WebSearch (extended) x2: no multi-stream density-power subhalo measurement; Phoenix
  limits only from 1606.04946; the uniform DES density paper 2206.04241 gives no subhalo constraint. arXiv API:
  stream AND subhalo AND "power spectrum" -> only GD-1 (incl. DESI GD-1, 2606.15863) and Pal 5; stream AND
  (Jhelum OR Phoenix OR ATLAS OR Orphan OR Indus OR Tucana) AND (subhalo OR gaps OR "density variations") ->
  single-stream descriptions (Indus peaks read as epicycles, 2603.09305), the 14-stream forecast, a review
  (2610.22202) and 2011-12 gap counting (1109.6022, 1201.1347); no joint power-spectrum measurement. The single-stream
  GD-1 kinematic constraint 2510.02247 is noted.
