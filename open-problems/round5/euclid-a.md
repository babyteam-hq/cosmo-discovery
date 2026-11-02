# Round 5, harvester euclid-a: Euclid DR1-Foundation, day-one measurements (2026-11-02; search and reasoning only, no code, no data)

## What DR1-Foundation is, and where it lies

- **Contents** (https://www.cosmos.esa.int/web/euclid/dr1-timeline, update of 15 June 2026): raw data, calibrated VIS/NISP
  images, catalogues and NISP spectra over about 1,900 deg², in November 2026. The clustering and weak-lensing products
  (calibrated shear, clustering catalogues with randoms) follow in mid-2027. The DR1 content memo
  (https://www.cosmos.esa.int/documents/10647/12245842/DR1-memo_250624.pdf) covers wide-survey data to 21 May 2025 and requires a
  ground-based counterpart (DES in the south, UNIONS in the north). I could not extract the memo's text (the PDF streams are
  compressed), so the item-by-item product list is unchecked.
- **Footprint**, read from the published DR1 input-coverage map
  (https://www.cosmos.esa.int/web/euclid/euclid-dr1-expected-coverage, `cgv_map_dr1input_equ_rgb_small.png`; 2,100 deg² input,
  about 1,900 deg² expected). By eye there is one large southern block at RA about 1h-6h, Dec about -20° to -62° (around EDF-S
  and EDF-Fornax, inside DES). That is roughly 1,300-1,500 deg². The rest is northern strips at Dec about +55° to +75°, at
  RA about 7h-11h and 16h-19h (the NEP/EDF-N side), roughly 400-600 deg². These are by-eye estimates from a figure. The HEALPix
  file `cgv_map_dr1input_o13.fits.gz` gives exact numbers at the gate.
- **Partner overlap**, inferred from the footprint:
  - ACT DR6 (Dec -63° to +23°) covers almost all of the southern block, about 1,200-1,400 deg² after its Galactic mask. The
    ACT DR6 maps paper (arXiv:2503.14451) gives 8,800 deg² of ACT x Euclid for the full survey.
  - eRASS1 (the western Galactic half) covers the whole southern block. SPT-3G also overlaps it.
  - DESI DR1 overlaps little: only the northern strips where DESI reaches Dec > 55°.
  - LoTSS-DR2 overlaps part of the northern strip at RA 7h-11h.

## The Consortium's announced DR1 Key Projects (the drop list)

From https://www.euclid-ec.org/about/contact/dr1-key-project-leads/ (318 planned key papers, per the NASA APAC Euclid update,
2025). They cover a lot:

- **Clusters:** cosmology, scaling relations, ICL, protoclusters, WL detection, SIDM (CL-1 to CL-8).
- **CMB cross-correlations:** CMB lensing x Euclid observables (CMBX-2), CMB lensing + photometric probes (CMBX-3), tSZ x
  GC/WL/clusters (CMBX-4), and kSZ peculiar velocities (CMBX-5).
- **Strong lensing:** galaxy-scale lenses, lensed quasars, cluster SL, lensed sub-mm galaxies (SL-1/2/4/5).
- **Cosmology and theory:** voids (GC-7), photo-z clustering (GC-8), extended models and isotropy/homogeneity (TH-1/2).
- **Galaxy evolution:** morphology and mergers (JL-2), SHMR (JL-3), cosmic web (GEV-1), quenching (GEV-2), dense environments
  (GEV-5), clustering (GEV-6), AGN demography (GEV-7).
- **High redshift and the local universe:** z > 7 galaxies and quasars (PU-1/3), the CIB (PU-7), LSB (LU-1), SBF distances
  (LU-4), Fornax (LU-3), Magellanic Clouds (MW-3).

Every candidate below was checked against this list.

## Candidates killed in this harvest (one line each)

- **Arc / Einstein-radius statistics of X-ray-selected (eRASS1) clusters**, the overconcentration question: the Consortium
  already runs an arc finder on Q1 (ARTEMIDE, arXiv:2511.03064), has Q1 ICL in eRASS1 clusters (arXiv:2610.29683), and has
  KPs SL-4 and CL-3. Scoop risk too high.
- **Splashback at z > 0.7:** done to z = 1.65 (MaDCoWS2, arXiv:2602.06434) and for eFEDS (arXiv:2301.03626).
- **SAGA vs ELVES quenched satellites** by Euclid LSB depth: the dispute is mostly resolved by matched definitions
  (arXiv:2210.03748; ELVES arXiv:2203.00014; SAGA III, Mao et al. 2024, doi:10.3847/1538-4357/ad64c4, at 1σ). The residual is
  below 10^6.75 Msun, which Euclid can't reach at the SAGA host distances.
- **Size vs environment at fixed mass at z ~ 1-1.5** (disputed: arXiv:1307.0003 larger vs 2212.06853 null): covered by GEV-5.
- **DES/HSC source n(z) high-z tail by clustering-z with NISP Hα:** HSC isn't in the DR1 footprint, and DESI already reaches
  these redshifts.
- **CMB Cold Spot supervoid at z > 0.9:** the spot's centre (RA 3.25h, Dec -19.5°) lies just outside the block's northern edge
  by eye, and the question is widely treated as closed (Mackenzie et al. 2017, 2CSz; Kovács et al. 2022, DES).
- **eRASS1 S8 via Euclid size-magnification masses:** the precision can't settle 0.86 vs ~0.80. **LoTSS-DR2 n(z) from Euclid
  hosts:** the north overlap is small and the S8 offset it would test is ~1σ.

---

## Survivor 1: spiral spin-direction (handedness) asymmetry at 0.3 < z < 1.2 with a mirror-equivariant classifier on Euclid VIS

**Question.** Is there a cosmological excess of one spiral winding sense, rising with redshift and aligned with a dipole
axis, as Shamir claims? Or is the published excess an artefact of classifier or annotator chirality bias and of unresolved
arms in ground-based images?

**Shape: a disagreement, settled by removing the disputed step.**

- **The claims:**
  - SDSS z < 0.3 (arXiv:1207.5464).
  - DES, 739k galaxies (arXiv:2208.00912).
  - DESI Legacy South, 8.1e5 galaxies, a 4.66σ dipole near the CMB Cold Spot (arXiv:2106.07118).
  - HST COSMOS, 8.7k galaxies, 2.8σ (arXiv:2011.03723).
  - HSC DR3 (arXiv:2410.15269).
  - JWST JADES, 158 clockwise vs 105 counter-clockwise, about 20 % ± 6 % (Shamir, MNRAS 2025; arXiv:2403.17271 for the earlier
    34-galaxy version).
- **The rebuttals:**
  - A reflection-equivariant CE-ResNet on SDSS/DESI images drops a 7σ human-label asymmetry below 1.8σ (Jia et al.,
    arXiv:2210.04168).
  - Patel & Desmond (arXiv:2404.06617) and Stiskalek & Desmond (arXiv:2410.18884) find isotropy.
  - Iye et al. 2021 find a 0.29σ dipole.
  - Shamir rebuts each (arXiv:2404.13864, 2411.08723, 2208.00893, 2309.03418).
- **What remains open:**
  - The redshift-growing part of the claim (HST, JWST) has only been tested on small space-based samples with Shamir's own
    pipeline.
  - The equivariant-classifier rebuttal stops at low z with ground-based images.
- **What Euclid removes:** Euclid VIS (0.16" FWHM) resolves spiral arms to z ~ 1 over about 1,900 deg². That takes away both
  disputed steps at once: unresolved arms, and an asymmetric classifier (by construction, with a mirrored-input null).

**DR1 products and partner data.**

- **Euclid DR1 MER VIS mosaics and the MER catalogue** (positions, VIS/NISP photometry, sizes, morphology) and PHZ photo-z (if
  PHZ is in Foundation; it was in Q1). Euclid Science Archive: https://eas.esac.esa.int/ ; DR1 page:
  https://www.cosmos.esa.int/web/euclid/euclid-dr1
- **Gaia DR3** for star masking: https://gea.esac.esa.int/archive/
- **DESI DR1 spectroscopic redshifts** for the northern strips (optional): https://data.desi.lbl.gov/doc/releases/dr1/
- **Shamir's published axes and catalogues** (DES, DESI Legacy, HSC) for the a-priori predictions.

**Method, in four lines.**

1. Select disc galaxies with VIS r_e > 3 PSF FWHM, 0.3 < z_phot < 1.2 and VIS < 23.
2. Classify winding sense with a reflection-equivariant network (the CE-ResNet design of arXiv:2210.04168) and with Shamir's
   own Ganalyzer pipeline. Each is run on the images and on their mirror images; the mirrored run must return the exact
   negative asymmetry.
3. Measure A = (N_cw - N_ccw)/N per patch (south block, north strips), in redshift bins, and as an amplitude along Shamir's
   published axes. The axes and predicted amplitudes are fixed before computing.
4. Errors and nulls:
   - Errors from jackknife over tiles, random reassignment, and injection of a known asymmetry into mirrored subsets.
   - Ellipticals and point-like sources pushed through the classifier must give A = 0.
   - A must not depend on focal-plane position, PSF orientation or galaxy size/PSF.

**Overlap and feasibility.**

- No partner overlap is required: the test is internal to Euclid, and Gaia is all-sky.
- **Sample size:** Galaxy Zoo Euclid on Q1 classified a few hundred thousand galaxies over 63 deg². Scaling to DR1 gives
  ~10^7 resolved galaxies; at a conservative 2-5 % with a confidently classified winding sense, that is ~2-5 x 10^5 spirals at
  0.3 < z < 1.2. To be checked on the DR1 catalogue at the gate.
- **Cost:** classification is a GPU-day; the whole measurement takes days.

**Expected separation and anchor.**

- **Anchors:**
  - HST COSMOS, 2.8σ on 8.7k galaxies, implies a dipole amplitude of about 5 % (σ_A ≈ sqrt(3/N) ≈ 1.9 %).
  - JWST JADES gives about 20 % ± 6 % (3.3σ).
  - Shamir's DESI-South axis lies near the Cold Spot, adjacent to the Euclid southern block, so the block sits near a claimed
    pole.
- **Euclid precision:** with N ~ 2 x 10^5, σ_A per patch is about 0.25-0.4 %. A 5 % claim projected onto the south block
  (cos θ ~ 0.7-0.9) is excluded or confirmed at ≥ 9σ. A residual 1.5 % amplitude still separates at ≥ 4σ.
- **What could erase it:**
  - A chirality in the Euclid VIS PSF. Mirroring the image also mirrors the PSF, so the mirror null does not catch it. It is
    tested by A vs size/PSF ratio and by focal-plane position.
  - Photo-z errors smear the redshift trend but not the sign.
  - The patchy two-region footprint weakens a free dipole fit, but not the per-patch test against an a-priori axis.

**Quick-kill queries and results.**

1. Web (extended): "Shamir galaxy spin direction asymmetry Euclid Q1": no Euclid analysis; returned the HSC, DES, JWST and
   rebuttal papers.
2. arXiv API: `all:galaxy AND all:spin AND all:direction AND all:asymmetry` (32 hits): none use Euclid. Rebuttals and
   classifier-bias work as above.
3. arXiv API: `(handedness OR chirality OR "spin direction" OR clockwise) AND Euclid`: one irrelevant hit (neutrino wakes).
4. arXiv API: `abs:spiral AND ("spin direction(s)" OR handedness OR clockwise)`, newest 30: nothing on Euclid, DESI or JWST
   since 2024. The 2026 SDSS filament-handedness paper (arXiv:2609.07792) finds parity symmetry, at low z only.
5. Web: JADES spin-direction rebuttal: no published rebuttal of the JWST 158/105 claim; status "no consensus".

**Why the Consortium is unlikely to publish it first.**

- It's not a KP. TH-2 (homogeneity/isotropy) is about counts, clustering and dipoles, and JL-2 is about morphology and mergers.
- The topic belongs to one contrarian author, and the Consortium has no programme reason to adjudicate it.
- The real scoop risk is Shamir running his own pipeline on DR1. That would not remove the disputed classifier step, so an
  equivariant, mirror-nulled measurement stays the settling one.

---

## Survivor 2: host-halo mass of Euclid-selected little red dots at 1.5 < z < 3.5 from ACT DR6 CMB lensing (plus angular clustering)

**Question.** Euclid's wide-area V-shaped compact sources ("LRD candidates") live in halos of one of two masses:

- about 10^11-10^11.5 Msun, like JWST LRDs (clustering/overdensity masses 10^10.1-10^11.3 at z ~ 7.5-3.5, arXiv:2606.02773;
  several 10^11/h, NEXUS arXiv:2606.09721; 10^11.3 for a z = 7.3 LRD, Nature Astronomy, doi:10.1038/s41550-025-02660-1);
- about 10^12.5-10^13, like obscured or red quasars (the quasar x CMB-lensing analogue: ACT, arXiv:1207.4543).

That answer decides whether Euclid's bright z < 4 sample is the same population as JWST's LRDs. The disagreement is live:

- Bisigello et al. (Euclid Q1, arXiv:2503.15323) find 3,341 candidates in 63 deg² at z = 0.33-3.6, with number density
  rising to z = 1.5-2.5.
- JWST spectroscopy finds a "reduced incidence of LRDs at z < 3" (arXiv:2606.02773).

The halo mass also separates LRD formation channels (DCBH vs SIDM vs PBH; assembly-bias predictions in arXiv:2603.15736 and
2512.18000).

**Shape.**

- **Regime extension.** LRD halo masses exist only from JWST pencil beams (≲ 1 deg², N ≲ 100). This is the first large-scale
  (linear-bias) halo mass, and the first with CMB lensing, which is free of the photo-z clustering systematics and of the
  stellar contamination that biases an auto-correlation.
- **Disagreement.** Euclid-selected and JWST-selected "LRDs" disagree in redshift trend. The lensing mass removes the disputed
  step: whether the Euclid colour/compactness selection picks the same objects.

**DR1 products and partner data.**

- **Euclid DR1** MER catalogue (VIS, Y, J, H photometry; compactness/size), PHZ photo-z, and the MER mosaics for the
  Bisigello et al. selection (reproduced on DR1 exactly as published); NISP spectra of the bright subset, as a purity and
  redshift check. https://eas.esac.esa.int/
- **DES DR2 griz** (in the MER catalogue as the external photometry): https://des.ncsa.illinois.edu/releases/dr2
- **ACT DR6 lensing maps, baseline and CIB-deprojected variants, with simulations:**
  https://lambda.gsfc.nasa.gov/product/act/actadv_prod_table.html
- **Planck PR4 lensing** (cross-check): https://github.com/carronj/planck_PR4_lensing ; PLA: https://pla.esac.esa.int/
- **Gaia DR3** for star rejection: https://gea.esac.esa.int/archive/
- Optional: SPT-3G lensing over the southern block.

**Method, in four lines.**

1. Reproduce the Q1 LRD selection on DR1. Add Gaia, VIS point-source and NISP-spectra star/brown-dwarf rejection, then split
   in photo-z: 1.5-2.5 and 2.5-3.5.
2. Measure the angular cross-power C_L^{κg} with the ACT DR6 baseline and CIB-deprojected maps (40 ≤ L ≤ 1000, the
   release's range), with the mask coupling corrected.
3. Fit b with the photo-z n(z) (widths marginalised), including the magnification term from the measured count slope s.
   Fit the auto C_ℓ^{gg} separately as a consistency check.
4. Errors from the ACT DR6 lensing simulations crossed with the real catalogue. Nulls:
   - Curl and rotated κ.
   - A Gaia-star sample.
   - The same pipeline on a control of non-V-shaped red galaxies matched in magnitude and photo-z.

**Overlap and feasibility.**

- **Area:** the southern block x ACT DR6 is about 1,200-1,400 deg² (estimated from the coverage figure; exact from the
  HEALPix map at the gate).
- **Sample:** at the Q1 density (53 deg⁻²), about 6-7 x 10^4 candidates. With 30 % lost to purity cuts, about 4-5 x 10^4.
- **Feasibility scaling** is anchored on our own R-139 forecast (`results-R139.md`; DESI DR1 3.5 < z < 5 quasars x ACT DR6):
  6.2k objects on 4,000 deg² give 1.4σ at b = 3. Scaling by sqrt(N) (x 2.7-3.2) and by the stronger lensing kernel at
  z ~ 2-3 (x ~1.3-1.5) gives:
  - S/N ≈ 4-5 for b ≈ 2.3 (10^11.3 halos at z ~ 2.5);
  - S/N ≈ 8-10 for b ≈ 4.5 (10^12.5).
  To be redone with the real n(z) at the gate.

**Expected separation and anchor.**

- **Nearest published anchors:**
  - JWST LRD halo masses (10^11.3 at z ~ 3.5, arXiv:2606.02773; large errors).
  - ACT x SDSS quasar lensing, 3.8σ (arXiv:1207.4543).
- **Separation:**
  - b = 2.3 vs 4.5 is a Δ(S/N) of about 4-5. Expected separation ≈ 4-5σ before contamination.
  - With 30 % unclustered stellar contamination left in, both b values dilute by 0.7 and the separation falls to about 3-3.5σ.
    So the stellar fraction must be bounded first.
- **What could erase it:**
  - Low-z dusty interlopers (b ~ 1.5 at z ~ 0.5-1, with a different lensing kernel). These are bounded by the NISP-spectra
    subset and the control sample.
  - The magnification term for a steep-count faint sample (as in R-139, comparable to clustering per unit 5s-2). It must be
    measured from the counts, not assumed.
  - CIB leakage into κ, handled with the CIB-deprojected ACT variant.

**Quick-kill queries and results.**

1. Web (extended): "Euclid Q1 little red dots clustering halo mass bright LRD wide area": only the Bisigello Q1 paper, which
   names large-scale clustering as future work, and the theory paper arXiv:2603.15736. No measurement.
2. arXiv API: `abs:"little red dots" AND (clustering OR "halo mass" OR lensing OR bias)`, 101 hits, the newest ~60 read:
   - every observed halo mass is JWST, ≲ 1 deg², N ≲ 100 (2606.02773, 2606.09721);
   - wide-area entries are only forecasts (2609.06926, 2609.26553) or simulations (2608.02739);
   - no CMB lensing, no Euclid.
3. Web (extended): `"little red dots" "CMB lensing" ACT Planck halo mass Euclid`: none. Closest analogues are ACT x quasar
   lensing (arXiv:1207.4543) and DESI LRG x ACT/Planck (arXiv:2407.04607).

**Why the Consortium is unlikely to publish it first.**

- No KP names LRDs.
- CMBX-2/3 are cosmology-sample programmes, timed to the mid-2027 clustering products.
- GEV-7 (AGN demography) is the nearest home, but its natural tools (clustering with randoms, galaxy-galaxy lensing) wait for
  the 2027 products.
- The Q1 LRD team are Consortium members, though, and could do this. Moderate scoop risk, so speed matters: the measurement
  is days of work on day-one catalogues plus public ACT maps.

---

## Notes for the referee and the gate

- Survivor 1 is a foundations/isotropy test against a contrarian literature. Its value is a clean, pre-registered null or
  detection. The referee may judge it low payoff.
- Survivor 2 sits next to R-139 in method (CMB-lensing bias of a rare population; parked there for precision). The difference
  is ~10x more objects and a better lensing kernel. Its gate must first bound stellar contamination and the magnification
  slope.
- Both footprint numbers are by-eye readings of the ESA coverage figure. The first gate step is intersecting the HEALPix map
  with the ACT DR6 mask.
