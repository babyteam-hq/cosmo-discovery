# Proposed test: Large-scale alignment of radio quasar polarisation angles, re-tested with broadband, Faraday-derotated angles

**Question.** Are the intrinsic (RM-derotated) polarisation angles of radio quasars coherently aligned over ~150-500 Mpc,
as claimed from JVAS/CLASS 8.4 GHz, or randomly oriented?

**Shape: disagreement.** On the same JVAS/CLASS data:
- Tiwari & Jain 2013 (arXiv:1201.5180, 4,400 sources with P > 1 mJy) found alignment above 3 sigma at about 150 Mpc
  (n_v about 10, S_D statistic, sigma_random = 0.005), and 5.4-6 sigma in other procedures.
- Pelgrims & Hutsemekers 2015 (arXiv:1503.03482, 4,155 sources) found significant alignment in quasar-only samples.
- Joshi et al. 2007 found none.
- Pelgrims & Hutsemekers 2016 (arXiv:1604.03937) tie the angles to the axes of large quasar groups.

The disputed step is instrumental: angle biases in a single-band survey (Tiwari & Jain's own von Mises test gives a
non-uniform PA distribution at 4.7 sigma, and NVSS shows PA clustering at multiples of 45 deg). An independent
instrument, frequency and calibration, with derotated angles, removes it. The arXiv:2207.05765 review lists the effect as
an open challenge to the cosmological principle.

**Data.**
- SPICE-RACS DR2: derotated angles and RMs (arXiv:2605.16917, CASDA).
- Quasar redshifts: Milliquas v8 and DESI DR1 QSO (URLs as in R6-1).
- JVAS/CLASS polarisation (Jackson et al. 2007), VizieR J/MNRAS/376/371, for the same-source comparison.
- VLASS Quick Look, https://cirada.ca/, for spectral index (flat-spectrum selection).

**Method (four lines).**
1. Select SPICE-RACS sources matched to spectroscopic quasars, with |b| >= 20 deg and an angle error below 15 deg
   (sigma_RM x lambda0^2). Make two subsets: flat-spectrum (alpha from RACS-VLASS) and all.
2. Compute the S_D / Z_D statistics with parallel transport in 2D and 3D, at the published n_v (2-2500) and scales.
3. Get significance from three nulls:
   - angle shuffles among sources inside each ASKAP tile and beam position, which keeps any instrumental pattern;
   - a beam-position angle regression;
   - mocks with leakage-like angle patterns injected.
4. Same-source check: for JVAS/CLASS sources also in SPICE-RACS, compare the 8.4 GHz angles with the derotated
   0.9 GHz angles, and rerun the JVAS statistic on the derotated angles.

**Overlap and feasibility.**
- **Footprint.** SPICE-RACS (dec < +49) covers the whole of the quasar sky the claim uses below dec +49.
- **Same-source subset.** JVAS/CLASS lies at 0 < dec < 75, so the shared strip is dec 0-49. That holds about 2,500
  JVAS/CLASS sources, about 1,000 of them polarised in SPICE-RACS. This is a 1-2 sigma instrument check only.
- **Main test.** About 5,000-15,000 SPICE-RACS sources match spectroscopic quasars (estimate for the gate;
  arXiv:2605.16924 already matched 2,483 DESI DR1 quasars after strict cuts).
- **Angle errors.** Per-source derotation error is about 2 rad m^-2 x 0.1 m^2, roughly 11 deg, small against the
  uniform 52 deg rms.

**Expected separation and anchor.**
- **Anchor:** the JVAS effect is above 3 sigma at N = 4,400, sigma(S_D) = 0.005, so Delta S_D is about 0.015.
- **Separation:** at N = 10,000, sigma(S_D) is about 0.0033, so the effect at its claimed amplitude would show at about
  **4.5 sigma**. N = 5,000 gives about 3.2 sigma.
- **What could erase it:**
  - Depolarisation selects against compact flat-spectrum cores at 0.9 GHz, so the sample may be lobe-dominated, a
    different population from JVAS. Report the flat-spectrum subset separately.
  - ASKAP leakage makes angle patterns on beam and tile scales (about 1-6 deg, i.e. tens to hundreds of Mpc at
    z ~ 1). The within-tile shuffle null controls this, at the cost of power on those scales.
- **Note:** a physical prior for real alignment is weak. The value is in settling a twenty-year anomaly with an
  independent instrument.

**Quick kill (5 queries).**
- arXiv API, polarization AND alignment AND (radio OR quasars) AND "large scale": 45 hits. The latest radio tests are
  2015-2016 (1503.03482, 1604.03937); 2009.01637 covers LQGs at 2.4 sigma.
- arXiv API, polarization AND alignments AND (quasars OR blazars OR "radio sources"): no test since 2017.
- arXiv API, ("polarization angles" OR "polarization vectors") AND alignment AND radio: only ELAIS-N1 fields and single
  sources.
- WebSearch (extended), SPICE-RACS polarisation-angle alignment: none.
- WebSearch (extended), NVSS/VLASS/LoTSS polarisation alignment 2020-2026: none. The LOFAR 2020 paper (Osinga,
  2008.10947) is about jet position angles, a different question, which round 4 already answered.

**Overlap with our work:** sky anisotropy in counts and dipoles (tried) is a different statistic, and this candidate
uses no counts.

---
