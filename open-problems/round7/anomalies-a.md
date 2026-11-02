# Round 7, harvester A: unretested anomaly claims in LSS, clusters, CMB secondaries and lensing (2026-11-02)

Brief: published claims at >= 3 sigma (2018-2026) that nobody has re-tested on independent data, where public data now
allow a re-test that separates "claim real" from "chance" by >= 3 sigma. Search and reasoning only; no code run, no data
downloaded. Skipped as taken: the Giant Arc, the CMB cooling family (incl. the filament / hot-void follow-ups), spiral spin
asymmetry, quasar polarisation alignment, sky dipoles, plus everything killed in rounds 1-6 (parity 4PCF, Ho'oleilana,
supervoid ISW, Huge-LQG, LoTSS kpc-jet alignment, the lensing RAR, the SZ hydrostatic-bias redshift trend, ...).

**One survivor.** Most claims in these four areas are either already re-tested on other data, single objects (El Gordo,
the Bullet), comparisons with simulations that new observations cannot settle (galaxy-galaxy strong lensing), or under
3 sigma on what is public (the Giant Ring, the Big Ring, the GRB walls). Kills are listed at the end.

---

## A1. Quasar jet alignment regions (QJARs): do parsec-scale jets in four >5 sigma regions stay aligned on new VLBI sources?

**(1) Claim.** Mandarakas, Blinov, Casadio, Pelgrims, Kiehlmann, Pavlidou & Tassis 2021, "Local alignments of parsec-scale
AGN radiojets", A&A 653, A123, **arXiv:2106.06546**. Jet position angles (PAs) of 7,290 AGN with "good" visual jets in the
Astrogeo VLBI image database (version of 2020-06-10), redshifts from OCARS. Four regions in (RA, Dec, z), found by a grid
search and then optimised (particle swarm over centre, radius 0.5-2x, z range), whose PAs deviate from random at > 5 sigma:

| QJAR | RA | Dec | radius | z range | S.L. | mean PA |
|---|---|---|---|---|---|---|
| 1 | 36.18 | +25.95 | 14.82 deg | 0.02-1.53 | 5.54 sigma | -16.9 deg |
| 2 | 168.75 | +47.54 | 18.11 deg | 0.36-1.89 | 6.74 sigma | -17.7 deg |
| 3 | 180.15 | +23.56 | 16.07 deg | 0.29-0.65 | 5.72 sigma | -25.6 deg |
| 4 | 201.83 | +50.97 | 18.27 deg | 0.07-0.35 | 5.38 sigma | +37.6 deg |

Statistic: a parallel-transported dispersion (S-test) D over the sources in the region, S.L. = 1 - Phi((D - <D>_MC)/sigma_MC)
from 10^4 random-PA shuffles. The authors call the regions possibly "some of the largest structures known" (~400 Mpc).
Caveats in the paper itself: no look-elsewhere correction for the grid plus optimisation; two of six initial regions
dropped for beam-driven N-S PAs at low declination; for two of the initial six "systematics cannot be fully discarded".
**Rebuttals:** none specific. A global test (Pelgrims/Blinov et al. 2020, arXiv:2002.00982, 6,388 PAs) finds no
*global* alignment, which doesn't address local regions. Jung et al. 2025 (arXiv:2502.03730) cite it as part of an
"inconclusive" question and study kpc-scale LOFAR jets vs the cosmic web, which is a different observable.

**(2) No independent-data re-test (searches).**
- WebSearch (extended): "parsec-scale radio jet alignment regions Mandarakas 2021 follow-up re-test new VLBI data": only
  the claim and the 2020 global paper.
- WebSearch (extended), citation style: "\"Mandarakas et al. 2021\" jet alignment regions large-scale structure": citing
  works are Jung et al. 2025 (LOFAR kpc jets, cosmic web) and a theory paper (arXiv:2108.05011); no re-test of the regions.
- WebSearch (extended): "\"QJAR\" OR \"QJARs\" ...": only the claim paper.
- arXiv API `all:QJAR OR all:QJARs OR abs:"jet alignment regions"`: 0 results.
- arXiv API `abs:jet AND abs:alignment AND (abs:VLBI OR abs:parsec-scale) AND (...)`, sorted by date: 21 entries since
  2021-06; none tests large-scale alignment of pc-scale PAs (2411.16951 and 2411.09099 use Plavin's PAs for astrometry
  and jet-host alignment).
- Overlap with our work: R-160 tests radio *polarisation* angles of quasars (SPICE-RACS). Three of the four QJARs overlap
  optical/radio polarisation-alignment regions, so the two tests are complementary, not the same: this one uses jet
  directions from VLBI images.

**(3) Independent public data and overlap.**
- New sources: Astrogeo held images of 15,235 sources at the claim's cut-off (2020-06-10). The RFC now has 22,839 objects
  (rfc_2026b), and rfc_2025d lists images for 22,184 sources, so about 7,000 sources were imaged first after the cut-off.
  None of them can be in the claim sample. They are fainter, so expect fewer good jets than the 2020 rate of 48 %; at
  25-35 % that is about 1,800-2,500 new PAs.
- New redshifts: sources the claim had to drop for lack of an OCARS redshift (unknown number; OCARS covered roughly half
  the VLBI sky then) that now have one from DESI DR1, SDSS DR19, or the current OCARS. These are also independent of the
  3D claim. DESI DR1 dark time covers QJAR2-4 (NGC, Dec +5 to +70) and part of QJAR1 (SGC north).
- Cross-check PAs (not independent sources): Plavin, Kovalev & Pushkarev 2022 (arXiv:2203.13750, ApJS 260, 4) give
  automated, visibility-based jet directions for 9,220 AGN (VizieR J/ApJS/260/4). They let us check our PA recipe against
  the claim era on the sources both cover.
- Footprint per region: each QJAR covers about 650-950 deg2 of the northern VLBA sky. The gate is the count of new sources
  with a good PA and a redshift in each region's (RA, Dec, z) volume.

**(4) Frozen template (from the claim paper, fixed before looking).** The four regions exactly as in Table 1: centre,
radius, z range. Each region's claimed axis is its published mean PA. Statistic per region: the V-test (Rayleigh test
with known direction) on doubled, parallel-transported PAs against the published mean PA,
V_k = sqrt(2 n'_k) * <cos 2(PA_i - PA_k)>. The headline is the four regions combined with Stouffer weights
(sum V_k / 2). The paper's own S-test is a secondary statistic (direction free, 10^4 shuffles). Null checks:
- the same V-test on new sources in the same declination bands outside the regions;
- the V-test against each image's restoring-beam PA;
- per-region rotation by 90 deg.

**(5) Expected separation (anchored on the claim).** The S-test D behaves like the mean resultant length of the doubled
angles. Under random PAs, <D> ~ 0.886/sqrt(n) and sigma_D ~ 0.463/sqrt(n), so the published S.L. gives
R*sqrt(n) = 0.886 + 0.463 Z: 3.45, 4.01, 3.53 and 3.38 for QJAR1-4. On n' new sources per region with the claimed
alignment, the expected V is sqrt(2) R sqrt(n) sqrt(n'/n), which is 4.9, 5.7, 5.0 and 4.8 times sqrt(f), with
f = n'/n. Combined, sqrt(f) x 10.2 sigma; under chance it is 0 +/- 1.
- **>= 3 sigma for f >= 0.09.** The plausible f of 0.25-0.6 gives 5-8 sigma.
- Optimisation inflates the claimed R (winner's curse). Even at half the published amplitude, f >= 0.35 still separates
  the two readings at >= 3 sigma.
- The gate: count n'_k (new good PAs with a redshift inside each region). Stop if sum-weighted f < 0.35.
- What could erase it:
  - The beam. Three of the four mean PAs lie within 26 deg of N-S, the direction of VLBA beam elongation. The beam-PA
    control and the outside-region control at matched declination test this.
  - Fainter new sources may have shorter, noisier jets. Plavin's ~10 deg PA errors cost only cos(20 deg) = 0.94.
  - PA definitions that differ from the claim's visual inner-jet PA. Calibrate on the overlap with Plavin.

**Data URLs.**
- Astrogeo VLBI image database: http://astrogeo.org/vlbi_images/ (FITS images and uv data, with dates)
- RFC catalogue, with the first-observation epoch per source: https://astrogeo.org/rfc/ (rfc_2026b: https://astrogeo.org/sol/rfc/rfc_2026b/)
- Plavin et al. 2022 jet directions: https://cdsarc.cds.unistra.fr/viz-bin/cat/J/ApJS/260/4
- OCARS redshifts: http://www.gao.spb.ru/english/as/ac_vlbi/ocars.txt
- DESI DR1 redshifts: https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/zcatalog/
- SDSS DR19 spectra: https://www.sdss.org/dr19/

**Method in four lines.**
1. Take Astrogeo/RFC sources whose first VLBI image postdates 2020-06-10, plus older-imaged sources whose only redshift
   postdates OCARS 2020 (DESI DR1, SDSS DR19). Match them to redshifts within 1".
2. Measure each jet PA with one frozen automated recipe (the core-to-brightest-jet-component direction on
   circular-beam-restored images, checked against Plavin 2022 on the shared sources). Keep the "good" class only.
3. Select the sources inside each QJAR's frozen (RA, Dec, radius, z range) and compute the V-test against the published
   mean PA. Combine the four regions; run the paper's S-test with 10^4 shuffles as a secondary.
4. Run the controls: matched-declination sources outside the regions, the beam-PA V-test, and a 90 deg rotation. Report
   the combined V against the predicted sqrt(f) x 10.2, and the published amplitude as confirmed or excluded.

---

## Screened and dropped (in-scope claims; one line each)

| Claim | Why dropped | Queries |
|---|---|---|
| Giant Ring (Lopez & Clowes 2026, arXiv:2604.17534): elliptical shells at 3.0-4.7 sigma, 2D power at 320 Mpc at 3.5 sigma | **Fails 5.** The template is frozen (the paper gives the ellipse a, b, theta, centre and 40 Mpc shell; z 0.722-0.862), but the predicted GA+NA shell's 94 absorbers vs ~66 expected (delta ~ 0.42) scale to ~31 expected on DESI DR1 independent sightlines (the R-155 GA ratio of 8.9/19). That predicts ~13 +/- 5.6 excess, **~2.3 sigma**, and its southern part is the GA that R-155 already excluded. The visually fit shell (4.0 sigma) is post hoc and reaches ~2.9 sigma. Retry at DESI DR2. | arXiv abs + HTML |
| Big Ring (arXiv:2402.07591, up to 5.2 sigma) | **Fails 4.** There's still no member table (absorbers appear only as IDs 1-51 in figures). The centre is given only as tangent-plane (0, 240) Mpc, and the alpha hull (delta = 0.75, 21.8e6 Mpc^3) isn't published. The anchored separation would be ~3.1 sigma even with a frozen hull (62 vs ~35 expected, scaled to DESI). | arXiv HTML; R-155 note |
| Galaxy-galaxy strong-lensing excess in clusters (Meneghetti 2020, arXiv:2009.04471; 2023 "persistent") | The anomaly is observed vs simulated; new clusters only re-measure the observed side, which isn't disputed. Hard to make a frozen statistic on public lens models of other clusters. | 1 WebSearch |
| Galaxy rotation coherent with neighbour motion (Lee et al. 2019, arXiv:1901.02151, 3.5 sigma < 0.8 Mpc; 1908.10972, 2.8 sigma at 1-6 Mpc) | The several-Mpc claim is < 3 sigma. The < 1 Mpc signal is expected in N-body (Kim et al. 2022) and was already searched for on MaNGA/SAMI (arXiv:2310.10723). | 1 WebSearch |
| Hercules-Corona Borealis Great Wall / Giant GRB Ring | The claims date from 2014-2015, outside the window. GRBs with redshifts since 2020 add only a few per region in 1.6 < z < 2.1, and northern follow-up bias dominates. | reasoning |
| El Gordo 6 sigma infall velocity (Asencio 2021) | A single object; it can't be re-tested on independent data. | reasoning |
| Cluster/S8, A_L lensing excess, birefringence | Re-tested by ACT/SPT/eROSITA/DES (outside "unretested"). | 2 WebSearch |

Search budget: 8 WebSearch (6 extended), 5 arXiv API queries, 7 fetches. Generic "anomaly + sigma" queries returned
mostly H0/S8 tensions; the QJARs came from an arXiv API query on alignment + sigma + large-scale.
