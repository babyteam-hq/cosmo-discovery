# Proposed test: CMB "cooling" around nearby galaxies, filaments and voids: an independent test (disagreement)

**Question.** Is there an achromatic CMB temperature decrement within the halos of nearby large spirals (z < 0.02), and an
excess in nearby voids, beyond ΛCDM plus known foregrounds?

**The two answers.**
- Claimed by one group: Luparello et al. 2022 (arXiv:2206.14217), Hansen et al. 2023 (arXiv:2305.00268), Garcia Lambas et al.
  2023 (arXiv:2310.13755, which explains the Cold Spot as this foreground), Hansen et al. 2024 (arXiv:2411.15307: 5.7 sigma,
  3-4 sigma after look-elsewhere), Dominguez Feldman et al. 2025 (arXiv:2506.08833, filaments, published in PRD) and Hansen
  et al. 2025 (arXiv:2506.08832: hot voids at z < 0.03, ~6.5 sigma, accepted in A&A). All of these use Planck PR3
  SMICA/SEVEM with 2MRS tracers and 100-10,000 simulations.
- ΛCDM with known foregrounds: no such signal. The tSZ of these halos is ~1-2 uK and chromatic (2506.08833 itself puts its
  100-217 GHz difference down to tSZ). Nobody outside the group has published a test, confirmation or rebuttal.

Strictly, only one side is a published measurement. This is a claimed >5 sigma detection set against the standard model,
never checked outside the group.

**Public data.**
- Planck PR4 (NPIPE) maps and their 600 end-to-end simulations: https://portal.nersc.gov/project/cmb/planck2020/ ; PR3 component maps and FFP10 sims: https://pla.esac.esa.int/
- ACT DR6 maps (f090/f150/f220), for sub-degree halo profiles and a second instrument: https://lambda.gsfc.nasa.gov/product/act/act_dr6.02/
- Tracers: 2MRS (http://tdc-www.harvard.edu/2mrs/) to reproduce the claim; new, independent tracers from the Siena Galaxy Atlas 2020 (sizes of ~380k large galaxies, https://www.legacysurvey.org/sga/sga2020/) with DESI DR1 BGS / 6dF redshifts; DESIVAST DR1 voids (https://data.desi.lbl.gov/doc/releases/dr1/vac/desivast/).

**Method.**
1. Freeze the claimed profile and cuts of 2411.15307 (spirals with R > 8.5 kpc, z < 0.015, pixels inside the halo) as the hypothesis, before looking at the data.
2. Apply them to galaxies **not in 2MRS** (SGA with DESI/6dF redshifts) and to DESIVAST voids, so the sky positions are new.
3. Stack ΔT profiles on Planck PR4 component-separated and single-frequency maps and on ACT DR6 f090/f150/f220; fit an amplitude A (claimed profile = 1) and the SED (achromatic vs tSZ vs dust).
4. Errors from the 600 PR4 sims (plus many Gaussian ΛCDM realisations), with the same masks; null: rotated tracer positions.

**Expected precision.** The CMB sky is shared, so a second instrument tests only processing artefacts; the independent test
is the new tracer positions. If SGA roughly triples the large-spiral sample at z < 0.02 outside 2MRS (still to be checked),
a real signal gives A = 1 at about 3 sigma, and chance or selection gives A < 0.3 at 95%. Frequency differences cancel the
primary CMB, so the SED test is limited by noise and foreground residuals, not by cosmic variance: well under 1 uK.
**Can tell apart: yes** for the halo claim; the void claim is weaker, because the local voids are few and large.

**Quick kill (6 queries, plus 3 papers read).**
- arXiv "CMB cooling galaxies foreground spiral": only the group's 2411.15307.
- arXiv "CMB temperature nearby galaxies 2MRS": only the group's 2310.13755 and 2206.14217.
- arXiv "CMB temperature decrement nearby galaxies Hansen Luparello": only 2310.13755.
- arXiv "local voids CMB temperature excess sign change ISW": only 2506.08832 itself.
- WebSearch "CMB temperature decrement around nearby large spiral galaxies Luparello Hansen independent test rebuttal foreground ACT SPT": group papers only; no rebuttal, no ACT/SPT test.
- WebSearch "independent test CMB cooling nearby galaxies filaments ... ACT DR6 OR Planck PR4 OR no evidence 2025 2026": no outside test.
- Read: the 2506.08833 HTML (PR3 SMICA, SEVEM frequency check, 100-1000 sims) and the 2411.15307 and 2506.08832 abstracts.

---
