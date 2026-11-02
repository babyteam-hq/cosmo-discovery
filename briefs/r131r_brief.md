# RESUME of R-131 (it failed on 2026-11-01 when a usage limit idled its worker key; nothing was wrong with the analysis)

Your predecessor landed partial work on staging: read `research/cosmo-discovery/results-R131.md`, its `finding.md` section if any, its `### cosmo-deep-*` notes in `LOG.md`, and `git log --oneline origin/staging | grep R-131` first. Continue from the last landed milestone; do not redo landed rows or rerun a gate pass already recorded (a recorded first-pass miss leaves one revision pass if it was a presentation point). Same rows file, same rules. The original brief follows.

# cosmo-deep-11: the stacked Sunyaev-Zeldovich signal of the eRASS1 clusters below the SZ catalogues' mass limit, against their weak-lensing-calibrated masses (R-131)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff at the top), `research/cosmo-discovery/blind-review.md` (the gate), `finding.md`
(the eRASS1 x ACT lensing section is the one pass of this kind: samples and statistic fixed before the fit, the claim its interval;
the eRASS1 x Planck section records the first y stack and the morphology split that twice failed the gate: do not claim a
morphology split here), `deep-1.md` (candidate A and its trajectory verdict), `results-R121.md` and `results-R122.md`. Reuse
`erass1_cmblens.py` (catalogue handling, the Planck y stack `ystack`) and `erass1_act.py` (ACT release handling, sparse annulus
matrices, the 400-simulation covariance machinery). Follow the agents' working rules: worktree on staging; keep-alive (a background loop running
the session keep-alive every 20 min); land after every milestone with
`git pull --rebase origin staging && git push origin HEAD:staging`. Rows in `results-R131.md` (own file). Long jobs with
`nohup python -u` into a log file and an `until` watcher; OMP_NUM_THREADS=1 in workers; disk ~28 GB, memory ~15 GB.

## Trajectory check: already run (R-121, 2026-11-01, context-free Opus), verbatim in `results-R121.md` row A:
"VERDICT: partly novel: the eRASS1 WL mass scale vs Planck SZ (hydrostatic-bias-corrected, 1−b≈0.84) masses is already compared
cluster-by-cluster for PSZ2–eRASS1 matches in Aymerich et al. arXiv:2509.02068 (and vs ACT DR6, ⟨M_eRASS1⟩=0.99⟨M_ACT⟩,
arXiv:2507.21459); a stacked Planck y–M relation for the eRASS1 catalogue below the SZ detection limit ... remains open" (read the
full row). The open part is this job: the stacked Y_SZ-M relation for the eRASS1 clusters at 1-4 x 10^14 M_sun, two decades
below the SZ catalogues, with the clusters' own weak-lensing-calibrated masses, on the Planck and ACT DR6 y maps. Before any
compute, spend 30 minutes confirming nobody published it since (search "eRASS1 stacked Compton-y", "eROSITA clusters tSZ stacking
Planck ACT", "Y-M relation eRASS1 weak lensing masses", 2024-2026; the eRASS1 scaling-relation paper arXiv:2511.14356 and the
eROSITA-SPT/ACT papers arXiv:2402.08455, 2507.21459 are the ones to read in their methods); if it is published, stop and report.

## The question
The eRASS1 masses are calibrated on galaxy shear and are the scale behind S8 = 0.86. The thermal SZ signal Y_SZ is an
independent observable of the hot gas; its relation to mass (Y ∝ M^(5/3) E(z)^(2/3) self-similar) has been measured on the SZ-
selected clusters above about 4 x 10^14 M_sun. Stacking the Planck MILCA/NILC and ACT DR6 Compton-y maps on the eRASS1 clusters
binned in weak-lensing-calibrated mass and redshift gives the Y-M relation down to 1 x 10^14 M_sun: its amplitude at fixed mass
against the SZ-calibrated relations (Planck, ACT, SPT) measures the eRASS1 mass scale with the hydrostatic bias as the one
nuisance, and its slope tests self-similarity at masses where feedback is expected to break it. The claim is the relation: the
amplitude and slope with intervals, and the implied mass-scale ratio with the hydrostatic-bias prior stated.

## Data (public; verify reachability; zenodo.org, LAMBDA data directories and GitHub release assets 403; PLA and NERSC work)
eRASS1 primary and cosmology catalogues (R-121's files); Planck 2015 MILCA and NILC y maps (PLA; R-121 fetched them; files over
~1 GB are cut by the proxy, resume with curl -C -); ACT DR6 Compton-y (the DR6 ILC y maps and masks on NERSC's ACT portal under
dr6.02 or the ymaps release; find the README); the Planck and ACT y-map simulations or noise realisations if released, else
rotated-map and random-position nulls as in R-122.

## Method
1. Fix before any stack, in `results-R131.md`: the samples (the cosmology sample as primary; all PCONT < 0.5 as supplementary),
   the mass-redshift bins (e.g. 5 mass bins x 2 redshift bins over 1-8 x 10^14, 0.1-0.8), the aperture (a matched filter with
   the Arnaud et al. 2010 universal pressure profile at each cluster's R500 from its catalogue mass, convolved with the map beam;
   and a fixed 10' aperture as a check), the statistic (Y_500 per bin, cylindrical to spherical with the UPP), the model
   (Y ∝ A (M/M_pivot)^B E(z)^(2/3) D_A^-2, A and B free), and the comparison relations (Planck 2013 Y-M with 1-b = 0.8;
   ACT DR6 / SPT calibrations as published).
2. Stack with the per-cluster posterior-averaged mass (R-122's lognormal), with the beam and pixel windows in the filter; the
   covariance from rotated-map stacks and random positions (and y-map sims where available); the two-halo and
   cluster-neighbour contamination estimated by the same stack on the catalogue's own neighbours.
3. Fit A and B per map; the mass-scale ratio b_M implied under the Planck relation with the hydrostatic prior stated and varied
   (0.7-0.9); by redshift; the two maps against each other; CIB-deprojected variants as the foreground check; the relativistic
   SZ correction stated.
4. The claim is the table: per bin Y_500 with its interval; A, B with intervals per map; the implied mass ratio with its interval
   and the prior it rests on. No morphology split; no sentence about S8 beyond what the ratio says.
5. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus
   agent, section alone; hints: the prior works above and anyone who stacked tSZ on eROSITA or X-ray-selected clusters). Verdict
   verbatim in `results-R131.md` and under the section. One revision pass on a presentation point; park on a physics point.
   Nothing on the site; a pass is the headline of your report.
6. End: `### cosmo-deep-11` in `LOG.md`. Final report: headline, the table, the null checks, friction.
   Budget 6 h; stop at the budget with everything landed.
