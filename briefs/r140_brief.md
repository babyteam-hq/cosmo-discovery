# cosmo-deep-13b: the LCDM null for the eRASS1 supercluster outside-member lensing statistic, from a simulated light-cone (R-140)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff at the top, and `### cosmo-deep-13`), `research/cosmo-discovery/blind-review.md`
(the gate), `results-R137.md` and the R-137 section of `finding.md` (the measurement and its gate verdict), and R-137's scripts.
Follow the agents' working rules: worktree on staging; keep-alive (a background loop running the session keep-alive
every 20 min); land after every milestone with `git pull --rebase origin staging && git push origin HEAD:staging`. Rows in
`results-R140.md` (own file). Long jobs with `nohup python -u` into a log file and an `until` watcher; OMP_NUM_THREADS=1;
disk ~28 GB, memory ~15 GB. Background tasks die at 10 minutes: use nohup for anything longer. `research/cosmo-discovery/data/`
is empty in a new container: refetch what R-137 used (it took ~20 min plus ~70 min for the 400 ACT simulations).

## Why this job
R-137 measured the ACT DR6 lensing excess outside the member clusters of 365 eRASS1 superclusters (E - E_control consistent with
zero; 95 % upper limits 1.0x and 0.5x the member model). The gate, verbatim: "VERDICT: novel but not sound: the central statistic
E − E_control is compared with zero, but no ΛCDM null has been computed for it, because the overlapping two-halo terms and the
isolated (underdense) control mean zero isn't the expected value, so the bounds can't be read as limits on mass outside the member
clusters." Novelty held. This job computes the missing LCDM expectation, so the measurement can be read against it. It is the one
revision R-137's verdict allows: the method, not the data, was incomplete.

## Method
1. A public simulated light-cone with halos and a matching CMB-lensing convergence map: the Websky simulation (halo catalogue and
   kappa map, public at the Websky data site / LAMBDA mirror if reachable; else the AGORA or MDPL2-based light-cones with lensing;
   verify reachability first and state which). Fix before running, in `results-R140.md`: the simulated cluster sample (halos above
   the eRASS1 selection mapped to M500 and redshift: use the eRASS1 selection function where released, else a mass-redshift
   threshold matched to the eRASS1 cosmology sample's distribution, stated), the supercluster finder (R-137's rebuilt FoF with the
   catalogue's linking length), the control (matched isolated-cluster superpositions, exactly as R-137), the apertures, and the
   beam and noise (add ACT DR6 noise realisations from the released simulations to the noiseless kappa).
2. Run R-137's pipeline unchanged on the simulation: E, E_control and E - E_control per sample and aperture, with the
   distribution over sky patches or noise realisations as the interval of the LCDM expectation.
3. The claim: the observed E - E_control (from R-137) against the LCDM expectation, as a difference with its interval, and the
   upper limit on the outside-member mass re-expressed relative to LCDM. If the simulation's selection cannot be matched closely
   enough to make the null meaningful, say so as the result.
4. Revise the R-137 section of `finding.md` to include the null, then the finding check again (after-writing prompt verbatim,
   fresh context-free Opus agent, section alone; hints as R-137's plus "anyone who computed supercluster or cluster-pair lensing
   expectations in simulations"). Verdict verbatim in `results-R140.md` and under the section. This is the last pass: park on any
   further physics point. Nothing on the site; a pass is the headline of your report.
5. End: append to `### cosmo-deep-13` in `LOG.md`. Final report: headline, the simulation used, the null
   table, the comparison, friction. Budget 5 h; stop at the budget with everything landed.
