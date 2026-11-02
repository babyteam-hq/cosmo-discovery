# Proposed test: Are Euclid's z < 2 "little red dots" really LRDs? A slitless-spectroscopy census of the photometric candidates (settles a disagreement)

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
