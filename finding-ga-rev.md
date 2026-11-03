# On new sightlines, a Giant Arc made real cannot be told from chance: once a real structure is required to reproduce the Arc's member overdensity through the same linkage recipe, it predicts δ ≈ +0.6 (68 %: +0.25 to +1.05) in the Arc's hull on DESI DR1 sightlines SDSS never observed; the hull holds 8 absorbers against 8.86 expected (δ = −0.10), 1.3σ below that prediction and 0.5σ below chance

## Claim

The Giant Arc (GA; Lopez, Clowes & Williger 2022, arXiv:2201.06875) is a ~1 Gpc group of Mg II absorbers at z ≈ 0.80, found by
single-linkage clustering (95 Mpc) of the Zhu & Ménard (Z&M) SDSS DR7 + DR12 absorbers in z = 0.802 ± 0.060, with
δρ/ρ = 1.3 ± 0.3 measured on its own members inside the convex hull of its member spheres (CHMS). Whether it traces real
matter or is a pattern the algorithm picks out of ordinary points (arXiv:2502.03515, 2505.11072) can be tested on quasar
sightlines that played no part in finding it. On 3,822 DESI DR1 quasar sightlines through the frozen CHMS that are in neither
DR16Q nor the Z&M searched lists, the DESI Mg II catalogue holds 8 absorbers (W2796 > 0.3 Å on the Z&M scale) against 8.86
expected from dN/dz measured per stratum of quasar flux, exposure depth and redshift on all 1,233,901 independent sightlines.

A member-measured overdensity is inflated by selection, so the 1.3 is not what a real structure would show on new sightlines.
We measure what it would show by injecting real overdensities into the real SDSS catalogue and re-running the recipe:

| injected shape | GA-like recoveries (of 200 sites per amplitude) | A_true given a GA-like recovery (16/50/84) | predicted δ in the hull (16/50/84) | P(N ≤ 8 \| real) | P(N ≤ 8 \| chance) | likelihood ratio real / chance |
|---|---|---|---|---|---|---|
| **H: uniform in the GA's CHMS (primary)** | 280 (peak 28 % at A_true = 1.0; 0 at 4.0) | 0.5 / 1.3 / 2.0 | **+0.25 / +0.62 / +1.05** | **0.10 (1.3σ)** | 0.32 | **0.37** |
| S: member spheres only (clumpy) | 170 | 1.6 / 3.0 / 4.0 | +0.15 / +0.26 / +0.42 | 0.22 (0.8σ) | 0.24 | 0.91 |
| L: hull of 60.8 Mpc spheres (diffuse) | 200 | 0.5 / 1.0 / 1.3 | +0.22 / +0.56 / +1.05 | 0.12 (1.2σ) | 0.32 | 0.42 |
| H, SDSS placement weighted by quasar g magnitude | 244 | 0.5 / 1.0 / 1.6 | +0.19 / +0.48 / +0.86 | 0.13 (1.1σ) | 0.33 | 0.46 |

The claim is the primary row: a real structure, uniform in the GA's hull, strong enough that the linkage recipe recovers a
GA-like group from it (30-70 members, 0.5-2× the GA's extent, member overdensity within ±0.3 of the GA's own 1.32 under the
same estimator), predicts δ = +0.62 (68 %: +0.25 to +1.05) on the new sightlines. The observed δ = −0.10 lies 1.3σ below it
(P = 0.10) and 0.5σ below chance. The data favour chance by a factor of about 2.7, which is not decisive. **The test cannot
tell a real Giant Arc from chance**: no shape or variant reaches 2σ (0.7-1.4σ), and the clumpy shape is indistinguishable
from chance (ratio 0.91). Measured against the published 1.3 directly, the hull's fitted amplitude, A = −0.23 (95 %: −0.90
to +0.63), excludes it at 3.2σ, but that comparison is not a test of reality, because a real structure recovered by this
recipe does not carry 1.3 on new sightlines.

## Data and method

- **Template and DESI count** (frozen before any DESI file was opened; unchanged from the first version): Z&M DR7 + DR12
  (69,206 after dropping 8,441 duplicates), 4,752 in the slice; single linkage at 95 Mpc (flat Ωm 0.27, h 0.7); the GA group
  has 46 members (the paper has 44); CHMS sphere radius 30.4 Mpc (half the mean MST edge). DESI DR1 main-dark quasars minus
  every one within 1.5″ of a DR16Q or Z&M-searched quasar; path redward of Lyα to 5000 km/s below z_q. The 8.86 is the sum
  over the 3,822 sightlines of each one's path in the hull times the DESI dN/dz of its stratum (TSNR2 tertile × 12 g-flux bins
  × 4 z_q bins); the all-sky check gives 7,844 expected against 7,864 observed.
- **Member overdensity estimator.** All slice absorbers in a group's CHMS over the SDSS expectation in CHMS ∩ slice, minus 1;
  SDSS expectation = path on the 129,410 Z&M searched sightlines (DR7 ∪ DR12; redward of Lyα, 5000 km/s below z_q) times one
  measured slice dN/dz per class (BOSS, DR7-only; the classes reproduce the slice's 4,752). On the GA rebuild it gives
  **1.32** (48 absorbers against 20.7), the published 1.3. Variant model (seven g-magnitude classes from DR16Q PSF g): 1.03.
- **Injection.** The GA's CHMS shape moved rigidly to 200 places in the SDSS footprint (> 25° from the GA, SDSS expectation
  0.75-1.33× the GA's, holding < 10 members of any real slice group of ≥ 40 members; 14 independent patches > 20° apart).
  At each site and A_true ∈ {0, 0.25, …, 3, 4}, Poisson(A_true × SDSS expectation) absorbers are added uniformly per unit SDSS
  path in the shape, into the real catalogue (real clustering, real chance groups), and linkage is re-run. The recovered group
  is the one of ≥ 20 members holding the most injected absorbers (at A_true = 0, the most members inside the shape).
- **What DESI would see.** For each GA-like recovery, two numbers on real DESI data: the share f_rec of the recovered hull's
  DESI expectation lying inside the injected structure, with both rotated back onto the GA's own sightlines (median 0.54: the
  recovered hull is offset from the structure that seeded it); and the chance part bg, DESI's own N_obs/N_exp − 1 in the
  recovered hull at the site (the injected absorbers are not in DESI): +0.04 under "real", +0.15 under "chance" (the real
  slice's own GA-like groups give +0.20). The hull's predicted count is Poisson(8.86 (1 + bg + A_true f_rec)) with a lognormal
  for clustering beyond Poisson (width 0.058, from the GA's shift null: sd 0.341 at 8.86 against Poisson's 0.336); compared at
  count level with the 8 observed. Posterior over A_true: flat prior on the grid, times the share of sites giving a GA-like
  recovery.
- **Process.** The injection plan passed a blind method check as "fixable" (seven points: count-level errors, the chance part
  measured per hypothesis, the DESI/SDSS relative bias, Monte Carlo errors and independent sites, a completeness variant,
  three shapes, neighbouring groups), all closed before any injection ran. The decision rule, fixed with the plan: an
  exclusion only if the real reading's P ≤ 0.023 (2σ) for the primary shape; otherwise the claim is that the test cannot
  distinguish. **Disclosure:** the DESI count (R-155, 2 November) was known when the plan was written.

## Checks

- **Window, prior, EW cut.** GA-like window ±0.5: P = 0.11. Prior 1/A_true: P = 0.15. W > 0.6 Å (8 against 8.40): P = 0.13.
- **Relative bias** r = δ_DESI / δ_SDSS for one structure, scanned: P = 0.28 / 0.19 / 0.14 / 0.10 / 0.08 at r = 0.25 / 0.5 /
  0.75 / 1 / 1.25. No plausible r reaches 2σ.
- **Monte Carlo.** Block bootstrap over the 14 independent patches: P (primary) 0.09-0.11. L(A_true = 4) = 0 for H and L, so
  the grid does not truncate them; for S it does not (L rises to 4.0), which only lowers that shape's prediction further.
- **Merges.** Shape H recoveries draw a median 0 members from a pre-existing real group; S, a median 28 (the clumpy injection
  mostly relabels a real group, which is why its f_rec is 0.06).

## What it means

The gap between "3.2σ below 1.3" and "1.3σ below a real Arc" is the selection shrinkage. To come out of the linkage looking
like the GA, a uniform real structure needs A_true of about 0.5-2. The hull the algorithm draws around it is offset, and only
about half of its DESI path lies in the structure. So the excess a real GA would carry on new sightlines is about +0.6, not
+1.3, and about 9 expected absorbers cannot separate +0.6 from chance (+0.15). The first version's exclusion of the
published amplitude stands as arithmetic, but it does not bear on whether the GA is real. This measurement does not decide
between the two readings.

## Weakest points

- One DESI release and ~9 expected absorbers: the test's power is the limit, not its systematics.
- Injected shapes are uniform or sphere-clumped; a real GA could be filamentary along the arc. Three shapes span the range
  and none reaches 2σ.
- The SDSS completeness is modelled by class, not per spectrum.

## To settle it

DESI DR2 (about 2× the sightlines) brings the count's error from ~0.34 to ~0.24 in δ; the real reading's +0.6 would then sit
about 2σ from chance, if the observed δ held at −0.10. The Big Ring is the next test once its member list is public.

## Blind review (2026-11-03, R-171, context-free Opus referee, own literature search, this section alone; `blind-review.md`)

VERDICT: the findings are confirmed novel and sound if the calculations hold

Strongest objection: the real-Arc prediction rests on uniform or sphere-clumped injections, not the thin filament the
discoverers describe, so "cannot tell" is shown for those shapes, and the DR2 forecast should be redone with a filament
injection. Closest work: Sawala & Teeriaho 2025 (arXiv:2505.11072); also 2502.03515, 2201.06875, 2402.07591, 2604.17534,
2504.14940, 1211.6256, 1306.1700, 1603.03260, 2412.15383. Section unchanged since review. Briefly withdrawn by R-171 under the
plan's rule; published as its interval at `site/research-giant-arc.html` (R-172): `results-R155.md` F8, F9.
