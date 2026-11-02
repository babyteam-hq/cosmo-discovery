# Proposed test: Does the Giant Arc reappear on independent sightlines? DESI DR1 Mg II absorbers inside the arc (disagreement)

**Question.** The Giant Arc (GA) is a ~1 Gpc arc of 44 Mg II absorbers at z = 0.80, found in SDSS quasar spectra. Is it
a real overdensity of absorbing galaxies, or a chance pattern that the linking algorithm picked out of Poisson noise?
A real structure must show up again as an excess of absorbers on *different* quasar sightlines through the same volume.
A pattern picked out of noise must not.

**Shape: disagreement.**
- *Real structure:* Lopez, Clowes & Williger 2022 (arXiv:2201.06875) give δρ/ρ = 1.3 ± 0.3, with significance ~4.5σ (SLHC/CHMS),
  3.0σ (Cuzick–Edwards) and 4.8σ (power spectrum). The same group's Big Ring (arXiv:2402.07591) and Giant Ring
  (arXiv:2604.17534) papers also back this side, as does their rebuttal of the simulation argument (arXiv:2504.14940, which
  finds no Gpc structures in FLAMINGO-10K).
- *Artefact:* Sawala et al. 2025 (arXiv:2502.03515) and Sawala & Teeriaho 2025 (arXiv:2505.11072) find the same kind of
  pattern in uniform random points and in FLAMINGO-10K halos. They argue the overdensity is "an algorithmic artefact",
  "typical of similar point patterns and unrelated to the underlying matter density".
- *The step they dispute:* whether the point pattern traces real matter. Both sides argue it with simulations of the
  *same* SDSS points. Neither has looked at independent sightlines: the "independent" checks so far are visual overlays
  of DR16Q quasars and DESI photometric clusters. DESI DR1 measures the absorber density on new sightlines, which settles
  the step directly.

**Public data.**
- DESI DR1 Mg II absorber VAC: 270,529 absorbers from 1.47M quasars,
  https://data.desi.lbl.gov/public/dr1/vac/dr1/mgii-absorber/v1.0/MgII-Absorbers-DR1.fits
- DESI DR1 QSO sightlines (clustering catalogue, plus `QSO_full` for z > 2.1):
  https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5/QSO_NGC_clustering.dat.fits
- Zhu & Ménard SDSS DR7/DR12 Mg II catalogues (the GA's own data), and the searched-quasar lists:
  https://www.guangtunbenzhu.com/jhu-sdss-metal-absorber-catalog (files `Trimmed_SDSS_DR7_107.fits`,
  `Trimmed_BOSS_DR12_107.fits`, `Searched_Quasars_*`)
- SDSS DR16Q, for the strict sightline exclusion: https://data.sdss.org/sas/dr16/eboss/qso/DR16Q/

**Method in four lines.**
1. Rebuild the GA from the Z&M DR7+DR12 trimmed catalogues using the paper's own recipe: z = 0.802 ± 0.060, single
   linkage at 95 Mpc, Ωm = 0.27, h = 0.7. Freeze the members (and the Big Ring and Giant Ring templates) before touching DESI.
2. Keep only DESI DR1 sightlines that are not SDSS DR16Q quasars (1.5″ match). For each, compute the comoving path inside
   the template (spheres of radius R around the members; R = 47.5, 70 and 95 Mpc, plus the CHMS hull) with the
   proximate zone removed. The expected count is that path times dN/dz, with dN/dz measured on the same sightlines over
   the whole sky.
3. Report δ_GA = N_obs/N_exp − 1 with matched EW cuts (W2796 > 0.3 Å and > 0.6 Å), weighting by sightline TSNR2.
4. Calibrate against two nulls. Template shifts across the DESI footprint give the sightline and coverage variance.
   Spheres placed on random sets of 46 SDSS slice absorbers give the baseline from ordinary absorber clustering, i.e.
   what a chance pattern built from real absorbers predicts. Quote δ_GA minus that baseline, with its interval.

**Overlap and feasibility (counted from the files).**
- The rebuilt GA is a 46-member group (the paper has 44; my DR7/DR12 deduplication differs slightly). It spans RA
  205.0–224.6°, Dec 36.2–42.2°, centred (216.0, 39.0), next to PSZ2 G069.39+68.05 at (215.4, 38.4) as in the paper.
- DESI DR1 coverage of the arc band (RA 202–229, Dec 35.5–41.5, 127 deg²): 77% of 0.5° cells hold ≥ 2 DR1 absorbers at
  any z. The band has 184 DR1 absorbers in the slice, and coverage is patchy at Dec > 40.
- Of 767,096 DESI NGC sightlines with z_q > 0.9, 735,460 are not in the Z&M searched lists. Using DR16Q will remove more.
- dN/dz on those sightlines in the slice is 0.074, from 6,554 absorbers.
- Sightlines through the template: 5,711 (R = 47.5), 7,795 (R = 70), 10,104 (R = 95). Expected null counts: 17.8, 39.1, 62.1.
- Null scatter of δ: σ = 0.29 for template shifts (39 shifts, median expected count 12.9) and 0.29 for the clustering
  null (40 sets, median expected 11.4). Both scale as roughly N_exp^−½, so σ ≈ 0.23 at R = 47.5 and about 0.18 at
  R = 95. The clustering-null mean is +0.07: this is what ordinary clustering alone gives.
- *What could erase or fake it:* (a) DESI re-observations of SDSS sightlines carry the same absorbers (removed with the
  DR16Q match). (b) Detection completeness varies with sightline S/N across the sky (the template-shift null and TSNR2
  weighting absorb this). (c) The EW floor differs between SDSS and DESI (matched cuts). (d) Selection bias: a real
  structure's true δ is below the δ of its own selected members. That is what the test measures, not a nuisance.
  Magnification and dust are negligible on 50–100 Mpc paths at z = 0.8. No assumed coefficient enters.

**Expected precision; does it decide?** σ(δ_GA) ≈ 0.2. If the GA is real with δ = 1.3 ± 0.3, the arc should show
δ ≈ 1.3 + 0.07. If it is a chance pattern, δ ≈ 0.07 ± 0.2. The two differ by more than 5σ, and even δ = 0.5 differs at
about 2.5σ. **Yes, it decides.** The unblinded feasibility count gives δ = −0.16 (15 observed vs 17.8 expected) at
R = 47.5, −0.28 (28 vs 39.1) at R = 70 and −0.29 (44 vs 62.1) at R = 95. That favours the artefact reading, but it is not
the measurement: SDSS exclusion is only partial, there is no EW matching and no TSNR2 weighting. The same pipeline
also runs on the Big Ring and Giant Ring templates (Anand et al. 2021 DR16Q absorbers). The Big Ring sits partly in the
sparse Dec > 40 part of DR1, so check its coverage first.

**Quick-kill queries and results.**
1. WebSearch (extended) "Giant Arc Big Ring Mg II absorbers DESI quasars independent test": found only the Lopez-group
   papers (2201.06875, 2402.07591, 2604.17534), whose independent checks are visual overlays of quasars and DESI
   photometric clusters.
2. arXiv API `all:"Giant Arc"`, sorted by date. 139 hits. I read the 34 newest (back to 2020-05, before the GA was
   found); all others are cluster-lensing arcs. GA papers: 2604.17534, 2505.11072, 2504.14940, 2502.03515, 2409.14894
   (review), 2402.07591, 2201.06875. None uses new sightlines.
3. arXiv API `all:"Big Ring" AND all:absorbers`: 2504.14940 and 2402.07591 only.
4. arXiv API `abs:DESI AND (abs:MgII OR abs:"Mg II") AND abs:absorbers`: 18 entries, all CGM, dust or catalogue work.
   None tests Gpc-scale structure.
5. WebSearch (extended) `"Giant Arc" DESI absorbers confirm OR refute independent sightlines 2026`: nothing.
6. arXiv API ultra-large/uLSS/gigaparsec structures since 2024: 2607.01172 (Sawala, a DESI DR1 *local* anisotropy claim,
   not the GA) and 2509.17696 (Quaia voids and clusters, whose photo-z is too coarse). Nothing on the GA.

---
