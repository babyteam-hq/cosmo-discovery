# Proposed test: Doppler magnification: first measurement of the size-density dipole in DESI DR1 BGS x Legacy Surveys sizes (regime extension)

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
