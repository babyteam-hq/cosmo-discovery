# Round 4, probes C: early universe and PBHs, neutrinos and light relics, GW and transients, gravity and foundations

Harvested 2026-11-01 with search and reasoning only: no code was run and no data downloaded. Every number below comes from a
published paper, cited at the point of use. Rounds 1-3, the area maps and the "Tried by us" list were read first.

**Result: two survivors, both tests of gravity using DESI DR1 BGS.** They share a sample, but they are not the same measurement:
- Survivor 1 is the size-convergence x density dipole.
- Survivor 2 is the number-count dipole.

Survivor 2 is close to round-1 candidate D: it uses the same statistic for a different question (see its overlap flag).
Nothing survived in PBHs, neutrinos, light relics or GW/FRB transients; see the table at the end.

---

## 1. Doppler magnification: first measurement of the size-density dipole in DESI DR1 BGS x Legacy Surveys sizes (regime extension)

**Question.** At z < 0.5, peculiar velocities make galaxies look larger or smaller (Doppler convergence,
kappa_v ~ (1/(aH chi) - 1) v.n/c). GR predicts a specific antisymmetric (dipole) cross-correlation between galaxy sizes and
galaxy number counts. Is that dipole there, with the amplitude A_D = 1 that GR predicts?

**Shape: regime extension.**
- Magnification from galaxy sizes has been measured only in the gravitational-lensing regime. Huff & Graves, arXiv:1111.1070
  (SDSS sizes and FP); Schmidt et al. 2012 (COSMOS); DES Y3 size/flux magnification is cited in 2511.17858 (Elvin-Poole et al. 2023).
- At z < 0.45 the Doppler term dominates kappa (2511.17858). This regime has never been measured.
- Theory and forecasts: Bonvin et al., arXiv:1610.05946; Andrianomena et al., arXiv:1810.12793; Ye et al., arXiv:2511.17858
  (DESI BGS + LSST, S/N >= 10 per z bin); Bacon et al., arXiv:1401.3694. All are forecasts or simulations.

**Public data.**
- DESI DR1 BGS LSS catalogues, randoms and weights: https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/
  About 5.5M reliable BGS redshifts (DESI 2024 III/VI, arXiv:2404.03002).
- Legacy Surveys DR9/DR10 Tractor catalogues, where shape_r is the half-light radius and shape_r_ivar its weight:
  https://www.legacysurvey.org/dr10/ . BGS targets were selected from DR9, so every BGS galaxy has a Tractor size.
- Siena Galaxy Atlas 2020 for the largest galaxies, where Tractor sizes are unreliable: https://www.legacysurvey.org/sga/sga2020/

**Method (four lines).**
1. Use physical size from the *observed* redshift, lambda = ln R_phys(z_obs). Subtract a smooth mean lambda(z) per bin, and
   per morphological type (Tractor SER/DEV/EXP and colour). The function of z only must be fixed before computing. Do **not**
   regress on absolute magnitude: Doppler magnification shifts the inferred luminosity too, and a lambda-M regression would
   remove about 0.6 of the signal.
2. Cross-correlate kappa_hat = lambda - <lambda>(z) with the number-count contrast. Use the Bonvin et al. 2017 dipole
   estimator, with a midpoint line of sight, separations 10-150 Mpc/h, and bins z = 0.05-0.15, 0.15-0.25, 0.25-0.35, 0.35-0.45.
3. Fit A_D against the linear-theory GR template. The template takes f from the same sample's quadrupole, and b and the
   magnification slope from the BGS luminosity function at r = 19.5. Add the wide-angle and intrinsic-size-density terms as
   free nuisances.
4. Errors:
   - covariance from size shuffles within z and type bins (kills the signal, keeps geometry and shot noise) plus spatial jackknife;
   - nulls: the even/odd tracer-swap null, a stars-as-sizes null (PSF leakage), and the same pipeline on redshift-shuffled randoms.
   - Result: A_D with a 68% interval per bin and combined.

**Overlap.**
- 100%: every DR1 BGS galaxy is a Legacy Surveys DR9 source with a Tractor shape.
- DR1 BGS density is about 5.5M over a bright-time footprint of roughly 7,500-9,500 deg2 (DESI 2024 III/VI; check the exact
  effective area in the LSS catalogue header). That is about 580-730 deg^-2.
- 2511.17858's forecast assumed a BGS-like density of 629 deg^-2 (276,578 galaxies in its 440 deg2 cosmoDC2 light cone) and
  5,900 deg2 of overlap. DR1 matches it in density and exceeds it in area.

**Feasibility and precision.**
- *Noise is set by intrinsic size scatter, not image quality.* 2511.17858 finds sigma_kappa,intrinsic = 0.39-0.48 against a
  measurement scatter of 0.03 at LSST Y1 depth.
- Legacy Surveys seeing (~1.2") might raise the measurement term to ~0.1 for r < 19.5 galaxies. It is still in quadrature
  under 0.4, so the LSST-to-Legacy loss is under 5%.
- 2511.17858's S/N per bin is 13.0, 13.0, 12.7 and 10.4 (z = 0.1-0.5) for 5,900 deg2. Scaled to DR1's comparable density
  and area, that is about 10-13 per bin, or above 20 combined.
- Allow a factor 2-3 loss for incompleteness, the smooth-trend removal and nuisance marginalisation, and it is still
  **>= 4 sigma per bin, about 8 sigma combined, against A_D = 0.**
- Anchor: the 2511.17858 simulation forecast (no real-data measurement exists).

*What could erase or fake it.*
- **The angular-size and Malmquist gradient inside a z bin** (more distant galaxies are smaller and, in a flux-limited
  sample, intrinsically larger). This gives a spurious dipole of the same order as the signal (about r/chi). The physical
  size at the observed z and the smooth <lambda>(z) remove it. The removal is a smooth function of z, which cannot absorb a
  pairwise dipole. Size-shuffled mocks with an injected Doppler kappa measure the transfer.
- **The intrinsic size-density correlation** is even in the line of sight and enters the dipole only through wide-angle
  terms, about r/d ~ 0.05 times a small coefficient. It is modelled, and is about 5x below the Doppler dipole at z ~ 0.1.
- **PSF and seeing patterns** are angular and symmetric under a front/back swap, so they cancel in the dipole (as 2511.17858
  argues). The stars-as-sizes null checks this.
- **Fibre assignment** is angular and symmetric too.
- No assumed coefficient sets the detection. The GR amplitude test uses f, b and s measured from the same sample.

**Quick-kill queries.**
1. WebSearch (extended), `"Doppler lensing" OR "Doppler magnification" detection measurement galaxy sizes spectroscopic survey SDSS DESI first measurement`:
   forecasts only (1810.12793, 2511.17858, 1610.05946). The only real-data size magnification found is gravitational
   (Huff & Graves).
2. arXiv API, `abs:"Doppler lensing" OR (abs:Doppler AND abs:convergence AND abs:sizes)`: 25 entries, none on real data
   (simulations 2601.13820, 2205.00458, 2111.05439, 2011.12936; theory and forecasts).
3. arXiv API, `abs:"Doppler magnification"`: 10 entries, all theory or forecasts. The newest is 2511.17858.
4. arXiv site search, "Doppler magnification galaxy": 23 results, none on real survey data. The Jan-2026 review 2601.05320
   lists it as future work.

---

## 2. Relativistic (gravitational-redshift + Doppler) dipole of BGS bright x faint at linear scales in DESI DR1 (regime extension)

**Question.** Does the antisymmetric dipole of the bright x faint galaxy cross-correlation appear at linear scales
(r > 25 Mpc/h) with the amplitude GR predicts, A_rel = 1 vs 0?

**Shape: regime extension.**
- The only measurement on real data is Alam et al. 2017, arXiv:1709.07855: BOSS CMASS, dipole at 2.8 sigma, dominated by
  scales of about 10 Mpc/h, in the nonlinear regime.
- Gaztanaga et al., arXiv:1512.03918 (BOSS), detected only the geometric wide-angle dipole and saw no relativistic signal.
- The Jan-2026 review arXiv:2601.05320 says no Stage-III or early Stage-IV detection exists.
- Forecasts for BGS:
  - Bonvin et al., arXiv:2306.04213 (relativistic N-body light cone): flux-limited 50/50 split gives 9.0 sigma at z = 0.25
    and 6.0 sigma at z = 0.35, over a BGS-like ~13,800 deg2;
  - 2004.08014: >= 4.4 sigma;
  - 2109.06012: S/N 10-20.

**Overlap flag.** Round-1 candidate D (`round1/d-dark-force-dipole.md`, verdict "never measured, beta too weak") used the
same statistic to fit a dark-force beta, with the GR dipole as a nuisance. This candidate measures the GR dipole itself, a
different question whose separation is about 5 sigma rather than below 1. Send it to the referee with D's file attached.
Scoop risk is high, because the forecast authors are DESI members.

**Public data.**
- DESI DR1 BGS LSS catalogues with PIP/fibre-assignment weights and randoms:
  https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/
- The DR1 BGS AbacusSummit mocks under `dr1/mocks/`, for the even-multipole covariance and the null.
- Legacy Surveys DR9 photometry, for the magnification slope s and the evolution bias.

**Method (four lines).**
1. In z bins 0.05-0.15, 0.15-0.25 and 0.25-0.35, split BGS at a fixed apparent-magnitude cut into bright and faint
   (2306.04213 fixes the 50/50 split; freeze the split before computing).
2. Measure the cross-correlation dipole xi_1(r) for 25-150 Mpc/h with the midpoint line of sight, and the even multipoles
   for b_B, b_F and f.
3. Fit A_rel times the linear GR template (Doppler + gravitational redshift; light-cone term) plus the fixed wide-angle term.
   - s_B and s_F come from the number-count slope at the cut; f_evo comes from the measured luminosity function.
   - The detection (A_rel = 0 vs not) does not depend on these coefficients; the GR-amplitude test does.
4. Errors and nulls:
   - covariance from mocks rescaled by the jackknife diagonal (2306.04213's "combined" covariance);
   - nulls: a random 50/50 split (identical populations, no dipole) and the hemisphere and NGC/SGC split.

**Overlap.**
- One dataset. DR1 BGS has about 5.5M redshifts (2404.03002) over roughly 7,500 deg2 effective, against about
  13,800 deg2 and full BGS density in 2306.04213.
- The dipole is shot-noise dominated (2306.04213), so S/N scales roughly as n*sqrt(A):
  0.6 (density) x 0.74 (sqrt of the area ratio), about 0.45.

**Feasibility and precision.**
- Scaled: 9.0 x 0.45 ~ 4 sigma at z = 0.25, and 6.0 x 0.45 ~ 2.7 sigma at z = 0.35. Combined about 5 sigma, plus the
  z = 0.05-0.15 bin, which 2306.04213 expects to be stronger.
- **The expected separation from A_rel = 0 is about 4-6 sigma.** Anchor: Alam et al. 2.8 sigma on CMASS at nonlinear scales.

*What could erase or fake it.*
- Fibre collisions depend only on angle, so they are symmetric under a front/back swap and give no dipole.
- Redshift failures depend on magnitude, not on pair order.
- The colour/flux "Doppler bias" (arXiv:2410.04705) is part of the physical Doppler term, set by the measured s.
- Wide-angle dipoles are of the same order and must be modelled from the measured b and f.
- Nonlinear scales below 25 Mpc/h are cut.
- This is weaker than survivor 1. If only one is dispatched, take survivor 1.

**Quick-kill queries.**
1. WebSearch (extended), "relativistic dipole galaxy cross-correlation bright faint measurement DESI DR1 gravitational
   redshift odd multipoles detection": forecasts only (2306.04213); no DR1 result.
2. WebSearch (extended), "first detection relativistic effects galaxy clustering dipole Bright Galaxy Survey DR1
   measurement 2025 2026": forecasts only; the 2601.05320 review says nothing has been detected.
3. arXiv API, `abs:dipole AND abs:relativistic AND (abs:DESI OR abs:BOSS OR abs:SDSS) AND abs:galaxies`: 8 entries.
   Only BOSS measured anything: 1512.03918 (geometric only) and 1810.05051 (wide-angle odd multipoles).
4. arXiv site search, `"gravitational redshift" dipole clustering`: 7 results. The only real-data one is 1709.07855
   (CMASS, about 10 Mpc/h).

---

## Dropped

| candidate | area | why dropped | queries |
|---|---|---|---|
| LoTSS-DR2 radio-jet orientation alignment, 2D and 3D with DESI spec-z | foundations | **Answered:** Osinga et al. 2020 (2008.10947, LoTSS DR1) found 2D alignment at ~4 deg driven by bright sources and none in 3D. DR2 is the same test on more data. ELAIS-N1 null 2303.00773. The jet-filament papers 2502.03730 and 2512.07985 are a different question | WebSearch x2, arXiv API x2, arXiv site search |
| Superradiance limits on ultralight bosons from GW231123 spins | light relics | **Done:** 2507.20979, 2507.21788 | 1 WebSearch |
| GW231123 strong-lensing interpretation | GW | Crowded (2512.19118, 2512.20890, 2512.19077, 2609.03935, GWTC-4 lensing 2512.16347), and the answer turns on waveform systematics, not on data a small team adds | 1 WebSearch |
| "Phoebe" LMC PBH microlensing claim (May 2026) and its rebuttal | PBH | A check of one event, already rebutted ("Phoebe - a Mirage of a Primordial Black Hole"). Not a measurement on a sample | 1 WebSearch |
| E_G at z > 1 with DESI ELG/QSO x CMB lensing | gravity | **Fails 1/3:** eBOSS QSO x Planck E_G exists. DESI DR1 E_G is done (2507.16098, 2505.20656, gamma = 1.17 +/- 0.11 read as a projection effect). No disagreement to settle | 1 WebSearch |
| PBH Poisson isocurvature from DESI DR1 Lya P1D (map 07 row 6) | PBH | **Fails 3:** the white-noise term is strongest at high k, and DESI's k_max (~0.04 s/km) is below MIKE/HIRES (2508.04767). Only a weaker, independent bound, with no two answers | none needed |
| Spin-2 PNG via DESI DR1 LRG intrinsic alignments (map 07 row 13) | early universe | **Fails 3:** sigma(f_NL^{s=2}) improves only ~1.5-2x on BOSS's +/-280 (2302.02925), with no target value to separate | none needed |
| PTA-amplitude vs local SMBH mass function from a DESI velocity-dispersion function | early universe | **Fails 3:** the disputed step is the massive-end mass function, already set by SDSS at z < 0.2. The result hangs on the assumed M-sigma normalisation | none needed |
| ALP spectral oscillations in XRISM Resolve spectra of NGC 1275 | light relics | **Fails 3:** the limit depends on the assumed cluster magnetic-field model | none needed |

Not searched (gaps, not conclusions): CHIME Outrigger-localised FRB samples (2026), Euclid DR1 if public, JWST
caustic-crossing transient rates as compact-object limits.
