# Round 3: transients and foundations, quick-kill survivors

Mapped 2026-11-01. Two shapes were sought: regime extensions and disagreements, each a measurement on public data with
errors from nulls, mocks or simulations. Every candidate got the strict quick kill (WebSearch extended, then the arXiv
API), a footprint check from the real public files, and a physics check. **One survives.** No transient-probe candidate
survived: each one was either already done or too weak to carry the signal (see the table at the end).

The footprint and feasibility numbers below come from downloading the public catalogues and counting. The scripts are
`research/cosmo-discovery/r3_ga_rebuild.py`, `r3_ga_desi.py`, `r3_lqg_parse.py` and `r3_lqg_desi.py`. **Caveat for the
blind review:** the feasibility count for survivor 1 is itself a first unblinded look at the answer (δ ≈ −0.2 ± 0.25).
It is not the measurement. The pre-registered version (below) must fix the template, the sightline exclusion and the EW
cut before computing. Anyone running it should know that this rough count already points to a null.

---

## 1. Does the Giant Arc reappear on independent sightlines? DESI DR1 Mg II absorbers inside the arc (disagreement)

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

## Dropped

| candidate | shape | why dropped | queries |
|---|---|---|---|
| Huge-LQG (U1.27) in DESI DR1 quasars on new sightlines: Clowes 2013 (1211.6256, q = 0.40, main set 1.2) vs Nadathur 2013 (1306.1700) and Park 2015 (1502.03563, chance groups) | disagreement | **Fails 3 at DR1** (power). DR1 covers about 35% of the field (RA 160–170, Dec 8–22): 7.5 vs 21 QSO deg⁻² in the slice. Spheres of R = 33 Mpc on the 73 members hold 47 DESI QSOs not in DR7 against 26.1 expected (δ = 0.80). The clustering baseline (spheres on random DR7 QSOs) is 0.30 ± 0.19; at R = 50 it is 0.37 against 0.13 ± 0.11. That is about 2σ above the chance-pattern prediction, but σ ≈ 0.25 can't separate q = 0.40 from it. Worth rerunning on DESI DR2 (not public yet) | arXiv API "Huge-LQG" OR "large quasar group" (28 entries; no test with new tracers); 2 WebSearches (none with eBOSS or DESI) |
| Cluster infall F_G with DESI DR1 BGS | extension | **Fails 3.** The only work is a forecast (2610.37782, σ = 0.033). Its own abstract names mass-observable scatter as the main degeneracy, and DR1 has no external mass calibration below ~10%. redMaPPer projection effects bias the redshift-space cluster–galaxy signal by 10–20%. There is no disagreement to settle, and the forecast group is placed to do it | 1 WebSearch, 1 arXiv API (infall AND clusters AND gravity: only 1211.1379 SDSS groups, 2610.06841 refracted gravity) |
| Ξ0 (modified GW propagation) from GWTC-4 spectral sirens | extension | Done: 2509.04348 (GWTC-4, Ξ0 = 1.4 +1.0/−0.4), 2605.27227 (GWTC-5, 1.1 +0.6/−0.3) | 1 WebSearch |
| Fraction of BBH mergers with ZTF AGN-flare counterparts, O4 | extension | Done: Zhu et al. (ApJ ae47ea, GWTC-4.0 × ZTF DR23, f = 0.07 +0.24/−0.05), ALeRCE O4a/b (2603.04342), 2603.26972 | 1 WebSearch |
| DES5YR low-z offset (Efstathiou 2408.07175 vs Vincenzi 2501.06664) via the same SNe in ZTF DR2 | disagreement | **Fails 2.** No overlap: DES5YR low-z SNe are CfA/CSP/Foundation (to 2017), and ZTF DR2 starts in 2018 | from round-1 map row 3 |
| Dark-siren H0 with DESI DR1 spec-z vs photo-z | extension | **Fails 3.** σ(H0) is about ±12 from GW numbers, so the spec-z shift is ≪ σ; close to our GW host bias work | map row 2 |
| SN Ia residuals × CMB lensing κ | new | **Fails 3.** S/N ~ 1 (the κ kernel peaks at z ~ 2) | map row 12 |
| Homogeneity scale in DESI DR1 | extension | Done: 2610.38083, 2606.07854 | sweep-4 |
| CMB hemispherical asymmetry in DESI clustering | extension | **Fails 3.** Planck limits A to ≲ 0.01–0.02 at the ℓ that DESI k-modes map to, and DR1 sensitivity is A > 0.04 (1707.06555) | map row 4 |
| Ho'oleilana single BAO shell (2309.00677) in DESI BGS | extension | **Fails 3.** BGS traces the same LSS in the same volume, so it can't tell a real shell from a found pattern; it only adds volume. No published counter-answer except a non-arXiv note | arXiv API (1 entry), 1 WebSearch |
| Radio-galaxy birefringence (polarisation vs jet axis) against the CMB β ≈ 0.3° | new | **Fails 3.** NVSS RM-derotation error is ~25° per source, and the angle calibration of ~1° is bigger than 0.3° | none needed |
| Wide-binary MOND test with DESI DR1 radial velocities | disagreement | **Fails 3.** RV errors of 1–2 km/s swamp an orbital signal of 0.1–0.5 km/s; Gaia DR4 (Dec 2026) isn't out | none needed |

Not searched (gaps, not conclusions): GW231123 lensing claims, kilonova and SLSN candles, Fermi-LAT lensed-blazar
autocorrelation searches.
