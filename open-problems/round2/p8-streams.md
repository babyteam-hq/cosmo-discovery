# Proposed test: Subhalo abundance from the density power spectra of several thin streams together

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
