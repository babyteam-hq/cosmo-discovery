# Round 5, Euclid B: day-one DR1-Foundation measurements (search and reasoning only; no code run, no data downloaded)

Harvested 2026-11-02. Scope: DR1-Foundation (Nov 2026: ~1,900 deg² wide, ~500 deg² north and ~1,400 deg² south, calibrated
VIS/NISP images, MER catalogues, PHZ photo-z, SIR/SPE slitless spectra; LE2 only, no LE3 clustering, shear or cluster
catalogues until mid-2027). Sources: https://www.cosmos.esa.int/web/euclid/dr1-timeline,
https://www.cosmos.esa.int/web/euclid/euclid-dr1-expected-coverage, the DR1 content memo
(https://www.cosmos.esa.int/documents/10647/12245842/DR1-memo-20250416.pdf/, read only through search snippets: the PDF
could not be rendered without running code).

**Consortium check.** No public list of the Consortium's DR1 key papers was found (one extended search; the ESA and EC pages
give none). As a proxy we read the full list of 41 Q1 special-issue titles (https://www.cosmos.esa.int/web/euclid/q1-papers)
and assumed every Q1 topic continues into DR1: morphology, dwarfs, bars, AGN and red quasars, LRDs at z<4 (photometric),
bright z>6 galaxies, strong lensing (galaxy, double-source-plane, cluster), clusters, cosmic web, ultracool dwarfs, known
transients, secondary nuclei, radio. Also the Q1 ICL fraction (2610.29684) and Q1 quasar spectroscopy (2512.08803,
2609.20947). Anything on those topics as a census was dropped. Both survivors use a step the Q1 papers did not take
(slitless spectra of LRD candidates; clustering), but they sit next to a Consortium topic, so the scoop risk is stated.

**Spectroscopic depth (it limits both survivors).** SIR extracts spectra only for H_E ≤ 22.5 (Copin et al., 2503.15307), and
Q1 quasar identification is secure to about H_E ≲ 21.3 (2512.08803). LRD-like sources have very high Hα equivalent widths
(several hundred Å rest frame), so their Hα reaches the 2e-16 erg/s/cm² line limit near H_E ≈ 22.5–23. Dusty
star-forming interlopers (EW ~100 Å) do not reach it.

---

## 1. Are Euclid's z < 2 "little red dots" really LRDs? A slitless-spectroscopy census of the photometric candidates (settles a disagreement)

**Question.** The Q1 photometric selection finds 3,341 LRD candidates over 63 deg². Their number density *rises* from high z
down to z = 1.5–2.5. JWST slitless and spectroscopic samples instead find a decline after cosmic noon, with LRDs ≲ 3 % of
AGN at z ~ 2 and classical AGN dominating above M5100 ≈ −22.5. That bright regime is exactly where Euclid selects.
In DR1, with H_E ≤ 22.5 and 0.84 < z_phot < 1.88 (Hα in the red grism), what fraction of candidates shows (a) broad Hα,
(b) only narrow lines, (c) no line? Combining this with stacked hot-dust and X-ray weakness separates LRD-like objects from
reddened quasars. What is the spectroscopically corrected density of LRD-like broad-line sources at M5100 < −22.5?

**Shape.** A disagreement settled by removing the disputed step: Euclid's photometric-only purity (the V-shape from broad
bands, with no spectra or IRAC across the wide survey). Prior work:
- Bisigello et al., Q1, 2503.15323: the candidates. The paper calls spectroscopic follow-up "pivotal".
- Kapoor, Matthee et al. (EIGER), 2607.00084: 5 LRDs at z = 1.55–3.18, ~7e-6 cMpc⁻³ at z = 1.9–2.5. LRDs ≲ 3 % of AGN, and
  classical AGN dominate brighter than M5100 ≈ −22.5.
- Lisiecki et al. (VIPERS), 2609.10319: 14 LRD-like sources at 0.5 < z < 1.75 with spectra; "rapid decline".
- Zhang et al., 2606.02773: abundance falls at z < 3.
- 2608.05267: SDSS local analogues, 0.08 % of broad-line objects.
- Method precedent for broad-line identification in slitless spectra: 2505.14025 (3D-HST G141).

**Data.**
- Euclid DR1-Foundation: MER catalogue (VIS, Y, J, H, compactness), PHZ, SIR 1D spectra and SPE redshifts and line
  measurements (https://www.cosmos.esa.int/web/euclid/euclid-dr1; archive https://eas.esac.esa.int/).
- unWISE/CatWISE W1/W2 and AllWISE W3 for the red side and hot dust (https://catwise.github.io/,
  https://irsa.ipac.caltech.edu/Missions/wise.html); IRAC is absent in the wide survey, so WISE replaces it.
- eRASS1 (https://erosita.mpe.mpg.de/dr1/) for stacked X-ray weakness, south only.
- Gaia DR3 (https://gea.esac.esa.int/archive/) to remove stars with proper motions.
- DESI DR1 (https://data.desi.lbl.gov/doc/releases/dr1/) for redshift checks and Mg II/Hβ in the north.

**Method (four lines).**
1. Reproduce the Bisigello selection on DR1 MER and PHZ (frozen before looking at spectra), with W1/W2 in place of IRAC.
   Keep H_E ≤ 22.5 and 0.84 < z_phot < 1.88.
2. Classify each SIR spectrum with a pre-registered fit: broad Hα (FWHM > 1200 km/s after the ~660 km/s resolution), narrow
   only, or none. Completeness against flux and FWHM comes from DESI/SDSS broad-line quasars and DESI ELGs in the same
   footprint, which serve as known-answer controls.
3. Split the broad-line objects into LRD-like and reddened quasars with stacks: W3/W2 (hot dust) and eRASS1 counts (X-ray
   weakness), against matched blue and red DR16Q/DESI quasars.
4. Turn the fractions into a corrected density by z and M5100, with binomial plus completeness errors. Nulls: the same
   pipeline on non-V-shaped compact sources of matched magnitude, and on random positions for the stacks.

**Overlap and feasibility.**
- Euclid alone for the census: the full ~1,900 deg². From Q1's 53 candidates per deg², DR1 holds ~1e5 candidates.
- The fraction with H_E ≤ 22.5 and z_phot in 0.84–1.88 is unknown. If it is 2–5 %, that gives ~2,000–5,000 spectra.
  **Feasibility gate (hour 1):** ≥ 300 such candidates, from the public Q1 catalogue scaled up or from the first DR1 tiles.
- DESI DR1 controls lie in the ~500 deg² north and part of the south. eRASS1 overlaps only the southern DR1 fields in the
  eROSITA-DE half; size this on day one from the coverage map.
- Everything is catalogue level, with no image processing: days, not weeks.

**Expected separation and anchor.**
- Hypothesis B (the photometric density is real): most candidates are LRDs, and ~80 % of LRDs show broad Balmer lines
  (2503.15323 intro), so the broad-line, X-ray- and hot-dust-weak fraction is ≳ 40 %.
- Hypothesis A (the JWST/EIGER/VIPERS decline): at M5100 < −22.5 LRDs are rare, so that fraction is ≲ 10 %, with the broad-line
  objects being reddened quasars (hot dust, X-ray normal) and the rest interlopers.
- With N ≥ 300, the binomial σ ≲ 3 %: a separation of ≳ 10σ, and ≥ 3σ even at N = 50.
- Anchor: EIGER's 5-object density, with a Poisson error of ~45 %, and VIPERS's 14 objects.
- What could erase the signal:
  - SIR contamination from overlapping spectra, which mimics broad wings. Use the multi-dither decontamination flags and
    the control quasars.
  - Hα falling outside the window because of photo-z errors for V-shaped SEDs. Report the "no line" class separately and
    bound it with DESI redshifts in the north.
  - At R ~ 450, Balmer absorption and the line shapes of LRDs cannot be resolved. The class therefore rests on the broad
    line plus the stacks, not on line shape.

**Quick-kill queries and results.**
1. arXiv API `abs:"little red dots" AND abs:Euclid` (9 entries). Only 2503.15323 (photometric) and 2604.16178 (Euclid+IRAC
   COSMOS z>4, no spectra) use Euclid data. No Euclid spectra of LRDs.
2. arXiv API `abs:"little red dots" AND abs:clustering` (77 results; the 37 most recent read). Clustering appears only for JWST
   samples at z > 2.3 (2606.09721, 2606.02773, 2610.23806).
3. arXiv API on LRDs at low z and number density: EIGER 2607.00084, VIPERS 2609.10319, SDSS 2608.05267. None mentions
   Euclid Q1.
4. WebSearch (extended) "Euclid NISP slitless spectra little red dot candidates broad H-alpha Q1 confirmation": none.
5. WebSearch (extended) "little red dots clustering halo mass z<4 Euclid Q1 angular correlation": none.
6. The Q1 quasar spectroscopy papers (2512.08803, 2609.20947) select with Gaia/WISE/DESI and do not mention LRDs.

**Why the Consortium is unlikely to publish it first.**
- Medium risk, stated plainly: the Q1 LRD team wrote that this follow-up is needed.
- But the Q1 LRD work is photometric, and the Q1 quasar spectroscopy works from Gaia/WISE preselection.
- The test needs non-Euclid partners: DESI as known-answer controls, eRASS1 and WISE stacks.
- It answers a question posed by JWST groups (EIGER, VIPERS) that the Consortium's AGN papers have not engaged.
- Consortium DR1 science papers follow internal review cycles of months. A catalogue-level public analysis can come out
  within days of release.
- **Precondition:** this requires that DR1-Foundation ships SIR/SPE spectra for H_E ≤ 22.5 as Q1 did (the ESA timeline says
  "spectra" are in Foundation).

---

## 2. Halo mass of bright LRD-like sources at z ≈ 1–1.8 from cross-correlation with Euclid galaxies (regime extension; depends on #1)

**Question.** Where do the survivors of #1, LRD-like broad-line sources at z ≈ 1–1.8 and M5100 < −22.5, live? In the
~10^11.3 M_sun halos measured for LRDs at z ~ 3.5, or in the ~10^12.5 M_sun halos of quasars at the same z and luminosity?
The answer comes from their cross-correlation bias with Euclid photometric galaxies in the same z slice, relative to Euclid
spectroscopic quasars in the same footprint.

**Shape and prior work.** A regime extension in redshift and luminosity.
- LRD halo masses exist only from JWST at z > 2.3:
  - Zhang et al. 2606.02773 (98 spectroscopic LRDs): 10^10.1 M_sun at z~7.5, rising to 10^11.3 at z~3.5, "approaching normal
    galaxies".
  - NEXUS 2606.09721: "several × 10^11", large errors, 2.3 < z < 7.4.
  - 2603.15736: assembly-bias predictions.
- Quasar bias at z ~ 1.5 is b ≈ 2.4–2.6 from eBOSS/DESI clustering.

**Data.** As #1, plus PHZ galaxies (H_E < 24) in Δz = 0.2 slices. The comparison quasars are Euclid SPE broad-line quasars,
and DESI DR1 QSOs in the northern DR1 footprint (https://data.desi.lbl.gov/doc/releases/dr1/).

**Method (four lines).**
1. Tracers: the spectroscopic LRD-like set from #1 (N_spec), plus the photometric candidates weighted by #1's purity by
   magnitude.
2. Compute w_cross(θ) with PHZ galaxies at 0.5–10 Mpc comoving (Landy–Szalay, the MER mask, randoms from the depth map), in
   two z bins.
3. Fit the same estimator to Euclid/DESI quasars matched in z and M5100. Report b_LRD/b_QSO, which cancels the galaxy bias
   and photo-z dilution to first order, and convert to M_h with a fixed Tinker et al. relation.
4. Errors from a 200-patch jackknife and from mocks with the photo-z scatter. Null: mismatched redshift slices, which should
   give zero cross-power.

**Overlap and feasibility.**
- Euclid-only for galaxies (~1,900 deg²). The quasar comparison uses Euclid SPE quasars over the full footprint: Q1's
  ~5,500 per 63 deg² scales to ~1.5e5.
- Scaling from quasar-galaxy cross-correlations of ~1e5 objects (b to ~1–2 %), N_spec ≈ 2,500 gives σ_b ≈ 0.2–0.3.

**Expected separation and anchor.**
- b(10^11.5 M_sun, z = 1.4) ≈ 1.6 against b(10^12.5) ≈ 2.5, a difference of ≈ 0.9.
- That is ≥ 3σ when σ_b ≤ 0.3. **Borderline:** it holds only if #1's gate gives N_spec ≳ 2,500, or a purity calibration to
  ±5 % for the photometric sample (then N ~ 1e4–1e5 and σ_b is set by the purity error).
- Anchor: 2606.02773's z ~ 3.5 point (10^11.3) and its trend.
- What could erase the signal:
  - Interlopers at the same redshift (dusty star-forming galaxies, b ≈ 2–2.5) pull the bias toward the quasar value.
    Hence the spectroscopic purity from #1.
  - Slitless confusion around close neighbours, which suppresses small-scale pairs. Fit only above 0.5 Mpc.

**Quick-kill queries and results.** Queries 1, 2 and 5 of #1: no LRD clustering below z = 2.3, and none with Euclid.
WebSearch "little red dots clustering halo mass z<4 Euclid Q1 angular correlation" found only 2603.15736 (theory) and
2503.15323 (no clustering).

**Why the Consortium is unlikely to publish it first.** Medium risk. Q1 did not measure clustering for LRDs, and the
Consortium's clustering products and analyses are LE3 work for mid-2027. A catalogue-level angular cross-correlation needs
only Foundation products. It depends on #1, so the two would ship together.

---

## Dropped this round (one line each)

| Idea | Why dropped |
|---|---|
| UDG abundance against eRASS1 cluster mass slope | Answered many times (1706.02704 1.11±0.07; 2210.00009 0.89±0.04; 2210.14994). Euclid dwarfs/LSB is a Q1 Consortium topic. Same question, new data. |
| Obscured vs unobscured WISE quasar halo mass, with Euclid dN/dz | The dN/dz step was already tested (Petter et al. 2023, ApJ acb7ef: only an implausible z ≈ 3.5 configuration erases the gap; 1705.05306). |
| Strong-lensing overconcentration from an eRASS1-selected sample | Strong-lensing clusters are a Q1 Consortium series (Bergamini et al.; Klein et al. on ICM-selected clusters). |
| ICL/BCG light against eRASS1 mass with space-based sky subtraction | Q1 ICL fraction published (2610.29684); the Consortium has an ICL group. |
| Size magnification by eRASS1 clusters (shear-free S8 check) | Mass precision about 5–10 % against a ~15 % gap: under 3σ. |
| Lensed fraction of z > 6 quasars at 0.1–1″ | Only ~30–50 known z > 5.7 quasars in the footprint. |
| Dual-quasar fraction at 0.1–0.4″ (a PTA amplitude link) | ~30–50 pairs expected. Separation-slope test under 3σ. |
| X-ray stack of Euclid LRD candidates on eRASS1 alone | Stacked depth ~2e-16 against the expected ~3e-16 erg/s/cm²: ~2σ. Kept only as a classifier inside #1. |
| [O III] doublet α variation from NISP quasar spectra | Dispersion-scale calibration ~1e-3 against claims of ~1e-6. |
| Cosmic-string lens pairs in VIS | Reaches only Gμ ≳ 1e-8, already excluded by the CMB. |
| SN Ia host near-IR step and hostless intracluster SNe (ZTF × Euclid) | ~100–150 and ~20–30 SNe in the footprint: under 3σ. |
| Optical EBL (IGL) from VIS counts against New Horizons | 3 nW/m²/sr gap with LORRI's ±1.65 error: < 2σ even with a perfect IGL. |
