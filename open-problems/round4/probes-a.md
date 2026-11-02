# Round 4, probes A: expansion/H0, dark energy, CMB and secondaries, LSS growth

Checked 2026-11-01 with search only: no code was run and no data downloaded. Sources were WebSearch
(extended), the arXiv API (export.arxiv.org, which answered this session) and arXiv site search.
Every number comes from published papers; the scalings are noted as estimates.

**One survivor, and its precision is marginal.** H0 and w0wa are still mined out (see round 3). In
CMB secondaries and growth, most of the measurements left open on public data were done in 2025-26
by the DESI-ACT kSZ and lensing series.

## Killed or dropped

- **Moving-lens effect, ACT DR6 x DESI:** first detection already published (2605.18938, ACT DR6
  NILC x DESI LS LRGs, 3.7-4.8 sigma). BOSS x Planck transverse-velocity dipoles: 2504.02525.
- **Pairwise and velocity kSZ, DESI DR1/DR2 x ACT DR6:** done (2511.23417, 9.3 sigma, LRGs;
  2510.14135, BGS; 2604.19744 and 2604.19745 for LRG, BGS and ELG; 2604.04867).
- **HSC Y3 high-z bin shifts (Delta z3, Delta z4) vs photo-z, decided with DESI clustering
  redshifts:** done (2511.18133 and 2511.18134; bins 3 and 4 move much less than the
  self-calibrated shifts).
- **"Lensing is low" retested with CMB lensing at small scales (fixed source plane):** galaxy-lensing
  systematics are already excluded by Lensing Without Borders (2111.13805). The live dispute is
  galaxy-halo modelling (2211.01744, 2305.09637, 2610.23828), which a CMB-lensing version still
  depends on. ACT DR6 reconstruction noise at L>1500 limits the one-halo regime anyway. Dropped:
  it doesn't remove the disputed step.
- **Supervoid ISW on DESI DR1 spectroscopic voids:** the disputed step is void selection on
  spectroscopic data. That dispute already exists on BOSS (Nadathur & Crittenden: A=1.64+/-0.53;
  Kovacs et al. 1811.07812 and DES Y3 + BOSS: about 4-5), and DESI LS superstructures are done (Hang et
  al. 2021: A=1.52+/-0.72). DESI DR1 would be the same test on a larger, largely overlapping sample.
  It also overlaps "void lensing" on our list.
- **Gravitational redshift in DESI DR1 clusters:** outside the four areas (a gravity test, not
  growth). SDSS versions exist (A&A 2023, 3058 clusters; 2102.11156), and 2601.19861 sets up the
  equivalence-principle version.

---

## 1. Rotational kSZ beyond z = 0.09: cluster gas rotation with ACT DR6 + Planck PR4 and DESI DR1 members

**Question.** Does the rotational-kSZ (rkSZ) dipole, aligned with the spin axis estimated from
member-galaxy redshifts, appear in an independent cluster sample at 0.09 < z < 0.25? What is its
amplitude A_rkSZ (on the same template), with an interval from random-position stacks?

**Shape.** Regime extension (redshift, sample and angular resolution) that also replicates the
first measurement independently.
- First evidence: Goldstein & Hill 2025, arXiv:2512.10951 (PRL 137, 161001). They used 25 clusters
  at 0.02<z<0.09 with at least 50 SDSS members, stacked on Planck PR4 NILC with tSZ deprojected and
  5' smoothing. The dipole is about 25 uK, 3.6 sigma by random stacks, and 2.9 sigma under their
  noise model with cluster-correlated noise. A_rkSZ = 1.05+/-0.32 against a v_rot=1500 km/s template.
- Groups: 2603.14494 (SDSS DR7 groups x Planck), tentative at 2.3 sigma.
- Forecasts only for SO: 2610.39977.
- The authors left ACT out because their clusters are large on the sky and ACT loses large
  scales. They expect clusters at z>~0.5 to need lensing-contamination modelling. Nobody has
  measured rkSZ at z>0.09, with ACT or SPT, or with DESI members.

**Public data.**
- DESI DR1 redshifts (zcatalog; BGS r<19.5, 2 mag deeper than SDSS main):
  https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/zcatalog/
- Cluster centres and masses: the ACT DR6 SZ catalogue (2507.21459; 10,040 clusters),
  https://lambda.gsfc.nasa.gov/product/act/actadv_prod_table.html , plus MCXC-II for X-ray centres.
- ACT DR6.02 f090/f150 maps and mnms noise sims:
  https://lambda.gsfc.nasa.gov/product/act/act_dr6.02/act_dr6.02_maps_prod_table.html
- Planck PR4 (NPIPE) frequency maps, for the large-scale leg and tSZ deprojection:
  https://pla.esac.esa.int

**Method.**
1. Find members (DESI DR1 + SDSS redshifts within 3 R500c and |dv| < 3 sigma_v) and fit the spin
   angle exactly as in 2512.10951: the hemisphere mean-z difference over 5-degree trials, keeping fits with r >= 0.8.
2. Build a tSZ-deprojected ILC from ACT DR6 f090/f150 and Planck PR4 30-545 GHz, with sources
   inpainted. Remove primary CMB with the E-mode template step of 2512.10951. Cut out 2 R500c,
   rotate each cutout to its spin axis and stack, weighted by member count.
3. Get the errors from at least 10^4 stacks at random positions with random angles, on the same
   mask and weights. Nulls: shuffled spin angles, a 545 GHz dust stack, f090 vs f150, and a stack of
   the tSZ-only map.
4. Fit A_rkSZ with the 2512.10951 template (and 2610.39977's shell model), in a z<0.09 bin (new
   clusters only) and a 0.09-0.25 bin, and compare.

**Overlap and feasibility.**
- **Footprints.** DESI DR1 BGS covers 7,473 deg^2 (2404.03000). Using the share of DR1 LRGs that fall
  on ACT DR6 f150 (1.69 M of 2.39 M, 2511.23417), the DESI x ACT overlap is about 70% of the DR1 dark
  footprint, roughly 4,100 deg^2. Planck covers the whole DESI footprint.
- **Cluster count (estimate).** SDSS gave 54 X-ray clusters with at least 50 members at z<0.09, and 25 of
  them pass the spin-fit cut. DESI BGS is about 1.9 mag deeper, so the same 50-member cut reaches
  z of about 0.2. The 0.09<z<0.2 shell holds about 9 times the z<0.09 volume. Applying DR1 fibre
  completeness in cluster cores (about 0.5, assumed) and the ACT overlap fraction (0.55) gives about
  130 clusters before the spin cut and about 60 after it, all independent of the 25. Clusters that
  fell short of 50 SDSS members at z<0.09 add more, at lower mass.
- **Physics.**
  - tSZ is even under the 180-degree flip, so the dipole projection cancels it, and the ILC
    deprojects it too.
  - CMB lensing of the primary gradient and the moving-lens dipole follow the CMB gradient and the
    transverse velocity, not the spin axis set by the member redshifts. They add noise, not bias.
  - CIB from members has no side preference under the velocity-gradient split; the 545 GHz null checks it.
  - Lensing contamination grows at z>0.5 (2512.10951), and we stop at 0.25.
  - No fitted coefficient sets the detection, since the significance comes from random stacks.
    A_rkSZ is relative to an assumed 1500 km/s template, so quote it as a ratio to the
    2512.10951 value, not as a velocity.
- **Resolution.** At z = 0.15, R500 is about 6', which Planck's 5' smoothing washes out and ACT's 1.4-2.1' beam
  resolves.

**Expected separation and its anchor.** The anchor is 3.6 sigma (2.9 sigma with correlated noise)
from 25 clusters, or about 0.6-0.7 sigma per cluster.
- With about 60 independent clusters at equal per-cluster S/N: 3.6 x sqrt(60/25) = 5.6 sigma (4.5
  with the correlated-noise discount).
- If the per-cluster S/N halves (fewer members, so noisier spin angles, and smaller angular size),
  that drops to 2.3-2.8 sigma.
- Higher mass at fixed member cut (signal goes roughly as tau x v_rot, about M^0.7) works the other way.
- Central estimate: about 4 sigma separation from no rotation, with a plausible range of 2.5-5.5.
- At A_rkSZ = 1.05 it measures the amplitude to about 25%, enough to test whether the z<0.09 value
  holds at z of about 0.15.
- **Risk:** it is marginal at the low end. The first gate is to count the real clusters that pass the 50-member and spin-fit cuts
  in DESI DR1 (the 0.5 completeness is assumed). Goldstein and Hill (Hill is in ACT) are the
  obvious people to do it next, so scoop risk is moderate to high.

**Quick-kill queries (5).**
- WebSearch extended "rotational kinematic Sunyaev-Zel'dovich galaxy clusters ACT DR6 detection
  DESI spectroscopic members": only pairwise kSZ (2511.23417) and DESI DR2 x ACT kSZ (2604.19744);
  no rkSZ.
- arXiv API abs:"rotational kinematic" OR abs:"rotational kSZ" (64 hits, 30 read): data results are
  only 2512.10951 (Planck+SDSS) and 2603.14494 (Planck+SDSS groups). 2610.39977 is an SO forecast.
- arXiv API (abs:rotation OR rotating) AND abs:Sunyaev AND (ACT OR Atacama OR SPT): no measurement;
  1909.04690 is a feasibility study for galaxy halos.
- WebSearch extended "rotational kSZ cluster rotation Atacama Cosmology Telescope measurement 2026
  Goldstein Hill follow-up DESI": no ACT rkSZ, no follow-up.
- arXiv site search "rotational kinetic Sunyaev" (19 results, newest first): no ACT, SPT or DESI rkSZ.

## Quick-kill queries for the dropped candidates

- Moving lens: arXiv API abs:"moving lens" (2605.18938 detection).
- HSC bins: WebSearch extended "HSC Y3 cosmic shear redshift bins 3 4 shift clustering redshifts
  DESI" (2511.18133, 2511.18134).
- Lensing is low: WebSearch extended '"lensing is low" CMB lensing small scales ACT DR6 DESI
  one-halo'; arXiv API abs:"lensing is low" (40 hits; 2111.13805, 2610.23828, 2606.23651).
- Supervoid ISW: WebSearch extended "supervoid ISW DESI DR1 spectroscopic voids Planck stacking
  excess" (no DESI DR1 spectroscopic version; BOSS and DESI LS versions exist).
- Gravitational redshift: WebSearch extended "gravitational redshift galaxy clusters DESI
  spectroscopic BCG members"; arXiv API abs:"gravitational redshift" AND cluster AND DESI (0).
