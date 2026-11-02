# Proposed test: Circumgalactic dust reddening around galaxies at z ~ 0.8-1.6

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
