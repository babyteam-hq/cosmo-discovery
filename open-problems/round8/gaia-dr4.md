# Round 8, Gaia DR4 harvester: cosmology measurements for the days after DR4

2026-11-02. Search and reasoning only; no code run, no data downloaded. DR4 is set for **2 December 2026, ~12:00 CET**
(66 months; https://www.cosmos.esa.int/web/gaia/release). DR4 contents (https://www.cosmos.esa.int/web/gaia/dr4): epoch
astrometry for every source in the main astrometric solution (~130 billion CCD records, 62 TB, via DataLink); about 312 M
mean radial velocities (`all_source_rvs`; DR3 had 33.8 M) and 6.9 billion epoch RVs; 3.37 M acceleration solutions,
5.9 M NSS orbits; `qso_candidates` 5.37 M; `gaia_crf4` 4.2 M; galaxy surface-brightness profiles; a microlensing catalogue with
lens masses.

**DPAC's own DR4 papers** (https://www.cosmos.esa.int/web/gaia/dr4-papers, 116 titles) take most of the obvious cosmology,
so these are out of bounds:
- Gaia-CRF4 (the glide / secular aberration drift).
- PV "Global Patterns in the Quasar Proper Motions" (glide, spin, the l = 2 terms, the low-frequency GW background).
- "Search for Continuous Gravitational Waves in Astrometric Residuals".
- Solar and planetary light bending (PPN gamma, the Jupiter quadrupole).
- PV "From the Magellanic Clouds to M31 and M33" (LMC/SMC/M31/M33 proper motions).
- PV "The sounding line in the Local Group" (the Cepheid ladder).
- PV "Time delays in lensed quasars".
- PV "High accuracy dynamics of the nearest stars".
- PV "Star clusters: detection, dynamics, and stellar variability" (a partial risk for survivor 1, below).
- The microlensing processing paper.

## Survivor 1: does the leading-tail excess of nearby open clusters survive radial-velocity membership? (MOND vs Newton)

**Question.** Kroupa et al. count more stars in the leading than the trailing tidal tail of every nearby open cluster and read
that as a 9 sigma falsification of Newtonian dynamics, with Milgromian N-body tails matching. All their counts come from
convergent-point (CP) membership, which projects on-sky velocities and doesn't use radial velocities. Does the asymmetry
survive when tail members are confirmed by radial velocity (6D), and when counts are compared within a completeness limit that
is the same for both tails?

**Shape and prior work.** A published disagreement that turns on membership.
- 2210.13472 (Kroupa et al. 2022): first claim. It notes the full tails of Praesepe, Coma Ber and NGC 752 also fit
  Newtonian symmetry.
- 2405.09609 (Kroupa et al. 2024, ApJ): epsilon = (n_l - n_t)/(n_l + n_t), counted from 10 pc to d_max.
  - Hyades: 234/184 (epsilon 0.12, Roeser+19) but 40/38 (0.03, Meingast & Alves).
  - Praesepe: 214/124 (0.27).
  - COIN-Gaia13: 222/77 (0.49).
  - Coma Ber: 44/21 (0.35, Fuernkranz) but 8/5 (Tang).
  - Combined 762/449, epsilon ~ 0.26, "8.99 sigma". The paper says the tail data "need significant improvement" and that
    teams' counts differ substantially.
- 2301.02251: the stochastic asymmetry in Newtonian models, where the Hyades is a 6.7 sigma outlier.
- 2103.12080 (Jerabkova+21): explains the Hyades asymmetry by an encounter with a massive (~1e7 Msun) Galactic lump, i.e. a
  dark-subhalo reading, Newtonian.
- 2411.13675, A&A 693 A127 (2025): MOND N-body runs, where the leading tail holds up to 2x the trailing one.
- 2501.17225 (Risbud, Jadhav, Kroupa 2025): CP tails for 19 of the 21 clusters within 400 pc. DR3 RVs appear only as a
  result (trailing tails at higher RV), never as a membership criterion.
- 2508.15056: catalogues give 38:31 leading:trailing, but tail length, magnitude limits and the Sun's position bias the counts.
- 2607.13747 (2026 review): projection and contamination can make or exaggerate tails; proposes RVs as the validation layer.

No paper counts leading vs trailing members with RV-confirmed membership and a matched completeness limit.

**DR4 products and partners.**
- `gaia_source`: DR4 astrometry.
- `all_source_rvs`: 312 M RVs, roughly to G_RVS ~ 15-16.
- `rvs_epoch_parameters_single` and the NSS tables: flag spectroscopic binaries, whose RV jitter would scatter members out.
- GSP-Spec [M/H]: separates the cluster (Hyades [Fe/H] ~ +0.13) from the Hyades moving group.
- No partner data needed. LAMOST DR10/11 RVs can serve as an external check.

**Method.**
1. For the four clusters in 2405.09609 (Hyades, Praesepe, COIN-Gaia13, Coma Ber) plus the 400 pc set of 2501.17225, select
   stars whose DR4 6D velocity lies within 3 sigma + 1.5 km/s of the velocity a Newtonian tail model predicts at their
   position (no CP projection).
2. Cut to an absolute-magnitude limit that is RV-complete at the farthest point of both tails, and drop RV-variable stars
   (epoch RVs).
3. Count n_l and n_t in the 2405.09609 regions (10 pc to d_max).
4. Compute epsilon per cluster and combined, with three nulls, each run through the same pipeline:
   - (a) field stars shifted in velocity by +-5 km/s;
   - (b) the Kroupa CP members, to reproduce their counts;
   - (c) Newtonian test-particle tails passed through the same selection.

**Feasibility.**
- One ADQL query per cluster, run on release day.
- Hyades tail stars lie 20-220 pc from the Sun. An RV-complete cut near G_RVS 15 keeps M_G <~ 7.5-8 (about 40% of members),
  giving n ~ 170 in the Hyades region.
- Praesepe (186 pc) and COIN-Gaia13 keep a similar fraction, so roughly 400-600 RV-confirmed members in total across the four.
- DR4 RV errors at G_RVS 14-15 are a few km/s (DR3 was ~6 km/s at 14). The tail dispersion is near 1 km/s and the field
  spread near 30 km/s, so the RVs remove most interlopers whose tangential motion matched by chance.

**Expected separation and anchor.**
- The anchor is the combined epsilon = 0.26 of 2405.09609, against Newton's ~0 (models < 1.4 sigma per cluster).
- With N ~ 500 RV-confirmed members, sigma_epsilon ~ 1/sqrt(N) ~ 0.045.
- If the asymmetry is real, epsilon ~ 0.26 sits ~5.8 sigma from 0. If it comes from CP contamination, epsilon ~ 0 sits
  ~5.8 sigma below the anchor. Either way it clears 3 sigma on the combined sample.
- The Hyades alone doesn't clear it: Roeser's 0.12 at sigma 0.077 is 1.6 sigma.
- The 400 pc set of 2501.17225 roughly doubles N.

**What could erase the signal.**
- Hyades moving-group stars share the cluster's velocity, so RVs can't reject them; GSP-Spec [M/H] only partly separates them.
- RV completeness that differs between the near and far tail, if the limit is mis-set.
- Unresolved binaries with RV offsets (lost from both tails alike, but cutting N).
- Newtonian epicyclic overdensities, or the bar and spiral arms (2606.06248), placed asymmetrically inside the d_max window,
  so the window must be fixed in advance and varied.
- The CP null (b) failing to reproduce the published counts would make the comparison ambiguous.

**Quick-kill queries.**
- WebSearch (extended) "open star cluster tidal tails asymmetry MOND Newtonian test Gaia DR3 radial velocities rebuttal": only
  the Kroupa-group papers and the 2026 review.
- WebSearch (extended) "Hyades tidal tails radial velocities leading trailing asymmetry dark matter subhalo encounter
  Jerabkova": no RV-membership asymmetry test.
- WebSearch "tidal tails open clusters asymmetry Milgromian 2025 2026 ... convergent point": 2508.15056, nothing with RVs.
- arXiv API `abs:"tidal tails" AND abs:asymmetr* AND abs:cluster*` (19 entries) and `abs:"tidal tails" AND abs:"open clusters"
  AND (Milgrom* OR MOND OR asymmetr* OR "radial velocit*")`: nothing past 2405.09609 tests the asymmetry with RV membership.
- Fetched 2501.17225, 2607.13747 and 2405.09609.

**Why DPAC is unlikely to publish it first.**
- DPAC papers don't test alternative gravity.
- The star-cluster PV paper may show DR4 tails, possibly with RVs, but a leading/trailing count against MOND isn't a DPAC kind
  of result.
- The real competitor is the Bonn/Prague group (Kroupa, Pflamm-Altenburg, Jadhav), who want the test and will run it on DR4
  quickly. The window is days to weeks.
- Our edge is the RV-membership and matched-completeness design, which their CP pipeline doesn't do.

## Considered and dropped

| Candidate | Why out | Evidence |
|---|---|---|
| Secular extragalactic parallax (CMB-frame dipole in nearby-galaxy proper motions; distance-free H0) | Below 3 sigma. Maloney, Zhang, Croft & Malanchev (2606.00280, DR2-DR3 position differences, 67k galaxies) get a near-field amplitude of -13 +- 102 uas Mpc/yr against an expected 77.8, i.e. sigma ~ 1.3x the signal. Paine+20 (1912.11935) forecast 5-10 sigma at end of mission (10.5 yr). Scaling by T^1.5 to 66 months gives 1.9-3.8 sigma, and correlated peculiar velocities offset the DR2 fit. The Croft group (a small team) is set up to rerun it on DR4. | web (extended) + 2 fetches |
| PBH dark matter from astrometric microlensing in DR4 epoch astrometry | Verma & Rentala (2208.14460) forecast f_PBH ~ 3e-4 at ~10 Msun for exactly 5 yr / 70 epochs, but assume a 1-sigma single-epoch threshold and **negligible background**. Toward the bulge, stellar lenses and astrophysical BHs outnumber f = 3e-4 PBHs by orders of magnitude. Toward the LMC (clean), my estimate gives ~270 events at f = 1, so a limit near 1e-2, no better than OGLE's 20-yr LMC limit (Mroz+24). Processing 62 TB is weeks, not days. | web (extended) + ar5iv |
| Glide vs redshift, quasar l = 2 GW background, LMC term in the glide | DPAC's CRF4 and "Global Patterns" papers. The LMC term (~0.2 uas/yr) is also only ~1.3 sigma at the expected DR4 glide error (~0.13 uas/yr, EDR3's 0.35 scaled by T^1.5). Our DR3 version: sweep-4 row C. | DR4 papers page |
| uHz stochastic GW background from epoch astrometry | DPAC's continuous-GW paper owns the residuals. A quasar-only limit is Omega ~ 1e-2 at 1e-7 Hz (my estimate), useless next to the N_eff bound. With stars, attitude modelling absorbs the signal, which only DPAC can undo. | reasoning |
| M31/M33, LMC/SMC proper motions, timing-argument LG mass, M31 collision | DPAC PV "From the Magellanic Clouds to M31 and M33". The timing argument on EDR3 is already done (Chamberlain+23, Benisty+22, Sawala+25). | DR4 papers page |
| Internal PM dispersions of dSphs (core vs cusp) | Precision. Sculptor/Draco: sigma ~ 22-25 uas/yr against a DR4 per-star error ~ 150-190 uas/yr, so the signal is ~1-2% of the error variance and the errors would need calibrating to better than 1%. | reasoning |
| Quasar time dilation in Gaia epoch light curves | The disagreement is settled: n = 1.28 +0.28/-0.29 (2306.04053) confirmed by 2501.04171. DR4 only tightens n around 1, with no >= 3 sigma separation to win. | web |
| Wide-binary MOND test with DR4 multiplicity screening; Keplerian rotation-curve decline; VPOS; oldest-star ages | Same question already answered several ways on DR3 (prior rounds' verdicts), with many groups queued for DR4. RV errors swamp the wide-binary orbital signal. | prior rounds |
