# Proposed test: host-halo mass of Euclid-selected little red dots at 1.5 < z < 3.5 from ACT DR6 CMB lensing (plus angular clustering)

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
