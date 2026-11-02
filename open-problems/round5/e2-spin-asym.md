# Proposed test: spiral spin-direction (handedness) asymmetry at 0.3 < z < 1.2 with a mirror-equivariant classifier on Euclid VIS

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
