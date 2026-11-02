# cosmo-deep-1 (R-121): candidates for one cross-community confrontation

2026-11-01. Shape required by the brief: a public model, reconstruction, map or catalogue-level prediction released by one
community, tested with an independent public dataset from another community, at a precision estimated before touching data,
where a search shows the two have not been confronted. Parked threads excluded (SN residuals vs dust; CF4 zero points and
bulk flows; NGC 1052; eRASS1 anisotropy; Pantheon+ anisotropy; count dipoles; CHIME DM vs structure; SN clustering bias).

## Killed on data before the gate

- DESI DR1 peculiar velocities (FP and TF) vs the 2M++ (Carrick+2015) velocity field that Pantheon+ uses: the natural
  complement of the CF4 result, but DESI DR1 ships no peculiar-velocity value-added catalogue
  (data.desi.lbl.gov/public/dr1/vac/dr1/ has no `pv/`; checked 2026-11-01). Not sent to a referee.

## The candidates

The thread running through A-C: the eROSITA community's eRASS1 cluster mass scale. The eRASS1 primary catalogue (Bulbul+2024,
arXiv:2402.08452; 12,247 clusters, half sky) gives every cluster an M500 from the count-rate-mass relation calibrated on
DES Y3, KiDS-1000 and HSC-Y3 weak lensing (Grandis+2024 arXiv:2402.08455; Kleinebreil+2024 arXiv:2402.08456), and that mass
scale is what makes the eRASS1 cluster-count S8 = 0.86 ± 0.01 (Ghirardini+2024 arXiv:2402.08458) sit above the cosmic-shear
values (0.76-0.78) and above what Planck SZ counts give with a hydrostatic bias of 1 − b = 0.8. Three independent public
datasets from other communities can test that mass scale, each by stacking on the catalogue: the Planck Compton-y map (CMB
community), the Planck CMB-lensing convergence (CMB community, a different observable), and DESI DR1 spectroscopy (the
galaxy-redshift community). Each is a measurement of the eRASS1 masses in another unit, with the same catalogue, so the
three can be combined and cross-checked. Candidate D is from a different domain, for breadth in the choice.

### A. eRASS1 masses vs the Planck Compton-y map: the Y_SZ-M500 relation at eRASS1's weak-lensing mass scale

- Model under test: eRASS1 M500 per cluster (WL-calibrated count-rate relation), with z and R500, for the 12,247 clusters.
- Independent data: Planck 2015 MILCA and NILC all-sky y maps (nside 2048, 10' FWHM; PLA `COM_CompMap_Compton-SZMap-
  milca-ymaps_2048_R2.00.fits` and the NILC twin, with the release's point-source and Galactic masks).
- Method: per cluster, a matched-filter (or 5R500 aperture) Y_500 with the Arnaud+2010 universal pressure profile at the
  catalogue's R500 and z, convolved with the 10' beam; stack in bins of eRASS1 M500 (and z); fit Y_500 E(z)^-2/3 D_A^2 = A ×
  (M500 / 3e14)^B. Compare with the Planck collaboration's own Y-M relation (Planck 2014 XX, arXiv:1303.5080: 10^-0.19 ×
  (M/6e14)^1.79 × 1e-4 Mpc^2 with 1 − b = 0.8). The ratio of the measured normalisation to Planck's at fixed mass gives the
  mass bias between the eRASS1 WL scale and Planck's hydrostatic-based scale: A/A_Planck = ((1 − b_eff)/0.8)^1.79. Nulls:
  the same filters at random positions inside the eRASS1 footprint and mask (background and noise), and the catalogue
  rotated in longitude; covariance from bootstrap over clusters and a jackknife over sky regions; trials: the bins and the
  two y maps. Decomposition: by redshift (z < 0.1, where eRASS1 is complete to low masses and the beam resolves R500;
  0.1-0.3; > 0.3), by X-ray extent likelihood, by count-rate S/N (Malmquist: at fixed true mass the X-ray-selected clusters
  near the limit are X-ray-bright, and X-ray and SZ scatter correlate), and by whether the cluster is in PSZ2. Joint fit:
  normalisation, slope and the eRASS1 mass-scatter dilution together (the bins are in a noisy mass, so the slope is
  attenuated; the normalisation at the pivot is not, to first order).
- Detectability (before data): Planck's PSZ2 detects a 5e14, z = 0.2 cluster (Y_500 ≈ 1e-3 arcmin^2) at S/N ≈ 5, so the
  per-cluster noise on Y_500 is ~2e-4 arcmin^2. The Planck relation gives a 2e14 cluster at z = 0.2 Y_500 ≈ 2e-4 arcmin^2
  (S/N ≈ 1 per cluster); eRASS1 has ~2,000 clusters per decade-third of mass around 2e14, so each mass bin stacks to S/N
  ≈ 30-40 and the normalisation at the pivot comes out to ~3-4 % statistical. The quantity tested is a 20 % mass offset,
  which moves Y at fixed mass by 1.79 × 20 % = 36 %: a 10σ distinction if systematics are held to 10 % in Y. Known
  systematics: the pressure-profile shape (UPP vs the eRASS1 clusters' own; ~10 %), CIB and radio leakage in the y map
  (few %, differs between MILCA and NILC), the X-ray-selection/SZ scatter correlation (10-20 % at the limit, which is why
  the S/N and z splits are in the plan), and the eRASS1 masses' own Eddington bias (the catalogue's masses are from the
  posterior with their selection model, to be checked in the catalogue paper).
- Closest prior work (to be verified by the referee): Planck 2013 XI (arXiv:1212.4131: Y-M down to 1e13 by stacking SDSS
  locally brightest galaxies); Planck 2015 XXIV (arXiv:1502.01597: SZ counts need 1 − b ≈ 0.6-0.7); Bulbul+2024
  (arXiv:2402.08452) and Ghirardini+2024 (arXiv:2402.08458) for eRASS1; any eRASS1 × Planck or eRASS1 × ACT DR6 tSZ
  stacking paper (2024-2026) is the comparison that would kill or narrow this.

### B. eRASS1 masses vs Planck PR3 CMB lensing: a CMB-lensing mass calibration of the eRASS1 catalogue

- Model under test: the same eRASS1 M500 catalogue.
- Independent data: Planck 2018 lensing convergence, the tSZ-deprojected TT reconstruction made for cluster work
  (PLA `COM_Lensing-Szdeproj_4096_R3.00.tgz`, alm to L = 4096, with the MV reconstruction as the cross-check).
- Method: stacked κ profile at the clusters (and a matched filter with an NFW convergence profile at the catalogue's M500,
  z and a concentration-mass relation), fitted for one amplitude A_κ against the eRASS1-mass prediction; nulls from random
  positions and from the PR3 lensing mean-field; covariance by bootstrap and sky jackknife. Splits by z (the lensing kernel
  suppresses z < 0.1) and by mass.
- Detectability: Zubeldia & Challinor 2019 (arXiv:1904.07887) got 1 − b = 0.71 ± 0.10 from 439 PSZ2 clusters in Planck
  lensing. Scaling by mass (κ ∝ M^0.8) and number: eRASS1's ~8,000 clusters at 0.1 < z < 0.6 with a mean M500 near 1.5-2e14
  give S/N ≈ (1.7/5)^0.8 × sqrt(8000/439) × 7 ≈ 12, so A_κ to ~8-10 %: a 20 % mass offset at ~2σ alone, but
  model-light (no gas physics) and therefore the arbiter if A finds a discrepancy. Systematics: tSZ leakage (the deprojected
  map is the point), the mean field, the concentration-mass relation (±5 %), miscentring (eRASS1 X-ray centres are good).
- Closest prior work: Zubeldia & Challinor 2019 (arXiv:1904.07887); Geach & Peacock 2017 (arXiv:1707.09369, redMaPPer in
  Planck κ); Madhavacheril+2020 (arXiv:2009.07772, ACT × redMaPPer); Raghunathan+2019 (arXiv:1810.10998, SPTpol × DES);
  any ACT DR6 or Planck PR4 lensing stack on eRASS1 clusters.

### C. eRASS1 masses vs DESI DR1 spectroscopy: the velocity-dispersion-mass relation of X-ray clusters

- Model under test: the same eRASS1 M500 catalogue (converted to M200c with a c-M relation).
- Independent data: DESI DR1 BGS_BRIGHT large-scale-structure clustering catalogues (NGC + SGC, 460 MB; RA, Dec, z,
  weights), r < 19.5, complete to z ≈ 0.2 for M_r < −20.
- Method: for each eRASS1 cluster in the DESI DR1 footprint at 0.05 < z < 0.3, members inside R500 (and R200) and
  ±3,000 km/s with iterative 3σ clipping and the biweight scale estimator; clusters with N ≥ 10; stack ln σ_v against ln
  M500 in bins and fit σ = σ_15 (h(z) M200 / 1e15)^α; compare with the dark-matter calibration σ_DM = 1083 (h(z) M200 /
  1e15)^0.336 km/s (Evrard+2008, arXiv:0702241) times the galaxy velocity bias of simulations (b_v = 1.0-1.1, Munari+2013
  arXiv:1301.1682). Nulls: the same apertures at the DESI randoms; interloper correction from the DESI randoms'
  redshift-shuffled field; splits by N, aperture, z.
- Detectability: ~1,500 eRASS1 clusters in the DR1 footprint with N ≥ 10 at 25 % per-cluster σ_v precision give the
  normalisation to ~1 % statistical, but M ∝ σ^3 and the velocity bias is known to ±5 % at best, so the mass scale is
  tested to ±15 %: the 20 % offset at ~1.3σ. A systematics-limited consistency check, not a decisive test, unless A and B
  need a third opinion.
- Closest prior work: Capasso+2019 (arXiv:1812.06094, SPT dynamical masses); Kirkpatrick+2021 (arXiv:2104.09541, SDSS
  redMaPPer velocity dispersions); Kluge+2024 (arXiv:2402.08453, eROMaPPer, spectroscopic cluster redshifts incl. DESI);
  any eRASS1 velocity-dispersion paper with SDSS or DESI.

### D. (other domain) The Planck y map vs the DESI DR1 BGS group catalogue: the tSZ of spectroscopic groups at 1e13-1e14.5

- Model under test: the DESI DR1 `bgs-groups` value-added catalogue's halo masses (group finder, abundance matched).
- Data: the Planck MILCA/NILC y maps as in A.
- Method and detectability: as A, binned by group mass; Planck 2013 XI reached 2e13 with 260,000 SDSS LBGs at ~3σ per bin,
  and the DESI groups (tens of thousands at 1e13-1e14) would give the Y-M relation at 1e13.5-1e14 at ~10 % per bin.
- Closest prior work: Planck 2013 XI (arXiv:1212.4131); Lim+2018 (arXiv:1712.08619, Y-M of SDSS/2MRS groups with Planck);
  the DESI × ACT DR6 tSZ stacks (2024-2025). Expected "partly novel" at best (the measurement exists for SDSS groups; the
  DESI catalogue is new). Listed for breadth; A-C carry the cosmological stake.

## Trajectory verdicts (context-free Opus referees, blind-review.md prompt verbatim, one per candidate)

| candidate | verdict (verbatim) |
|---|---|
| A (eRASS1 masses vs the Planck y map) | VERDICT: partly novel: the eRASS1 WL mass scale vs Planck SZ (hydrostatic-bias-corrected, 1−b≈0.84) masses is already compared cluster-by-cluster for PSZ2–eRASS1 matches in Aymerich et al. arXiv:2509.02068 (and vs ACT DR6, ⟨M_eRASS1⟩=0.99⟨M_ACT⟩, arXiv:2507.21459); a stacked Planck/ACT Compton-y Y_SZ–M500 relation over the full eRASS1 sample, including sub-SZ-detection masses and its mass/z/selection dependence, remains open |
| B (eRASS1 masses vs Planck PR3 CMB lensing) | VERDICT: novel: no prior work answers this question |
| C (eRASS1 masses vs DESI DR1 velocity dispersions) | VERDICT: partly novel: arXiv:2411.09735 (Li+2024) already shows for 22 eRASS1 clusters that kinematic (≈hydrostatic) masses fall well below the eRASS1 WL-calibrated masses (WL ~110% higher), and arXiv:2402.08453 / arXiv:2609.31852 already measure eRASS1 dispersions and eRASS1×DESI DR1 dynamics; what remains open is a direct, large-sample (~10^3 clusters) stacked σ_v–M500 test of the eRASS1 WL-calibrated mass scale with DESI DR1 spectroscopy |
| D (Planck y map vs the DESI DR1 BGS group catalogue) | VERDICT: partly novel: the stacked Planck matched-filter Y_SZ–halo-mass relation for spectroscopic groups at 1e12–1e14.5 is answered by Lim et al. 2018 (arXiv:1710.06856, 2MRS/6dF/SDSS groups), and for DESI DR9 groups by arXiv:2201.12591; what remains open is testing the specific DESI DR1 BGS group catalogue's halo masses against tSZ, which no published work has done |

Design points the referees added, carried into the plan: (A) the eRASS1 scale against SZ masses is known at high mass
(eRASS1 ≈ ACT DR6 masses to 1 %, arXiv:2507.21459; ≈ Planck's DES-calibrated masses, arXiv:2509.02068), so the information
is in the sub-SZ-detection masses and the mass, redshift and selection dependence, and Planck's early X-ray-selected stacks
came out low against models (the "SZ cluster-stacking puzzle"), so the Malmquist and correlated-scatter splits are needed.
(B) ACT DR6 lensing would be the stronger map (its MCMF stack reached 16σ, arXiv:2406.14754) but is unreachable from this
network (LAMBDA); the S/N of 12 is optimistic once tSZ/CIB leakage and the two-halo term are in; the eRASS1 scale has been
checked against optical weak lensing, SZ masses (SPT-3G arXiv:2607.01175, ACT arXiv:2507.21459) and hydrostatic masses
(Li+2024 arXiv:2411.09735) but never against CMB lensing. (C) velocity bias limits every such test to ±15 % in mass; the
eRASS1 release itself carries 1,499 velocity dispersions (Kluge+2024). (D) two-halo dominance below 1e13.5 at Planck's beam
(arXiv:1608.04160, 1706.03753) and dust/radio contamination (arXiv:2609.08939).

## Choice

**B.** The only plain "novel" verdict, a model-light observable (no gas physics between the mass and the signal), and the
cosmological stake is the eRASS1 mass scale itself. Fallback: A, narrowed to its open part (the stacked Y-M below the SZ
detection limit and its mass/z/selection dependence), which shares the catalogue, the footprint and the stacking code.

## Outcome (end of R-121)

B was run to paper depth: `finding.md` sixth section, rows B0-B7 in `results-R121.md`. A's y-map stack was done as B's
corroboration (row B5) rather than as a separate Y-M measurement; the narrowed A (the Y-M relation below the SZ detection limit,
with a pressure-profile matched filter) remains open for a later job, with the y-map, the mask and the per-cluster profiles already
in `data/erass1_cmblens/`. C and D not run.
