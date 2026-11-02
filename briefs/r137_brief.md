# cosmo-deep-13: the mass of eRASS1 superclusters outside their member clusters, from stacked ACT DR6 CMB lensing (R-137)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff at the top), `research/cosmo-discovery/blind-review.md` (the gate), `finding.md`
(the eRASS1 x ACT lensing section is the pass this job copies in shape: samples and statistic fixed before the fit, the claim its
interval, covariance from the 400 released simulations), `results-R122.md`, and `screen-1/p4-superclusters.md` (the proposal as the
referee saw it). Reuse `erass1_act.py` (ACT DR6 release handling, masks, sparse annulus matrices, the 400-simulation covariance,
`simstream.sh` to stream the simulations one at a time) and `erass1_cmblens.py` (eRASS1 catalogue handling, R-122's lognormal
per-cluster masses). `research/cosmo-discovery/data/` is gitignored and empty in a new container: refetch what you need (NERSC
for ACT; the eRASS1 DR1 site for catalogues) and say in your report how long it took. Follow the agents' working rules: worktree on staging;
keep-alive (a background loop running the session keep-alive every 20 min); land after every milestone
with `git pull --rebase origin staging && git push origin HEAD:staging`. Rows in `results-R137.md` (own file). Long jobs with
`nohup python -u` into a log file and an `until` watcher; OMP_NUM_THREADS=1 in workers; disk ~28 GB, memory ~15 GB.

## Trajectory check: already run (2026-11-01, context-free Opus, in the dispatching session), verbatim:
"VERDICT: novel: no prior work answers this question". The referee's cautions, which this brief builds in: the eRASS1 supercluster
catalogue (Liu et al. 2024, arXiv:2402.08454) has 1,338 systems, 818 pairs and about 520 with >= 3 members, western Galactic
hemisphere only, so the ACT overlap is a few hundred, not ~1,000; the outside-member excess depends on the supercluster aperture and
on how members are modelled; part of it overlaps the two-halo term already fitted in cluster-lensing stacks (arXiv:2404.17345,
arXiv:2310.04387). Closest in method: Lokken et al. arXiv:2409.04535 (oriented ACT/DES stacks around clusters), Tanimura et al.
arXiv:1805.04555 (tSZ between clusters in superclusters). Read those two methods sections before fixing the statistic. Do not
rerun the trajectory check.

## The question
What fraction of a supercluster's total mass lies outside its member clusters? The catalogue's own supercluster "mass" is the sum
of member M500 (weak-lensing calibrated). The ACT DR6 lensing convergence integrated over the supercluster's footprint, minus a model
of its members at their catalogue masses, measures the rest. The claim is that fraction (or the excess in units of the member mass),
with its interval, and how it compares with the same statistic on LCDM expectations (random groupings of clusters with matched
member count and redshift, which carry the two-halo term but no supercluster selection).

## Method
1. Fix before any stack, in `results-R137.md`: the sample (superclusters with >= 2 and with >= 3 members inside the ACT DR6 lensing
   mask, 0.05 < z < 0.6; report the counts), the aperture (the supercluster's member convex hull dilated by a fixed comoving margin,
   and a second, circular aperture of fixed comoving radius as the check), the member model (each member an NFW plus its own
   two-halo term out to the aperture, at its posterior-averaged mass, beam-convolved), the statistic (excess = integrated kappa in the
   aperture minus the member model, divided by the member model's integral), and the control (the same statistic on random groupings
   of non-supercluster eRASS1 clusters matched in member count, masses and redshift). The result is the excess minus the control, in
   both apertures.
2. Covariance from the 400 ACT simulations; random-position and rotated-map nulls; the member-mass systematic (the eRASS1 mass
   calibration's stated uncertainty) propagated as a correlated nuisance; tSZ leakage checked with the CIB/tSZ-deprojected
   (profile-hardened) lensing variant if released, else stated as a limit.
3. The same on Planck PR4 lensing as the like-for-like check.
4. The claim is the table: per sample and aperture, N superclusters, excess, control, difference with its interval. No sentence
   about baryons, filaments or cosmology beyond what the number says. If the sample is too small for a 3 sigma excess, the result
   is the upper limit, written as such.
5. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus
   agent, section alone; hints: the works above and anyone who weighed superclusters or inter-cluster matter with lensing). Verdict
   verbatim in `results-R137.md` and under the section. One revision pass on a presentation point; park on a physics point.
   Nothing on the site; a pass is the headline of your report.
6. End: `### cosmo-deep-13` in `LOG.md`. Final report: headline, the table, the null checks, the data refetch
   time, friction. Budget 6 h; stop at the budget with everything landed.
