# cosmo-deep-16b: does type-dependent photo-z isolation leakage explain the early/late split in the lensing radial acceleration relation? (R-147)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff), `research/cosmo-discovery/blind-review.md` (the gate), `finding.md` (passed
sections: samples and statistic fixed before the fit, the claim its interval), `results-R145.md` (R-145's feasibility work and the
by-product this job turns into a result: photo-z-isolated lenses with a real >10 % neighbour, 87-94 % within 3 Mpc; within 0.5 Mpc
38-48 % for quiescent vs 27-30 % for star-forming at fixed M*), R-145's scripts, and `open-problems/round2/p5-isolation-leakage.md`
(the proposal as the referee saw it). No shear catalogue is needed. Follow the agents' working rules: worktree on staging; keep-alive (a background
loop running the session keep-alive every 20 min); land after every milestone. Rows in `results-R147.md`.
Long jobs with `nohup python -u` and an `until` watcher; foreground waits cap at 10 min, bare `sleep` is blocked, match processes
with `ps -eo pid,args`. R-145 fetched the BGS group VAC (290 MB) in 21 s.

## Trajectory check: already run (2026-11-01, context-free Opus, 4 web searches and ~9 fetches), verbatim:
"VERDICT: partly novel: Brouwer et al. arXiv:2106.11677 already measured photo-z isolation contamination for the full sample (about
30% ESD excess beyond 0.3 Mpc/h70 in MICE, with σz = 0.02) but assumed it was the same for both types, and Mistele et al.
arXiv:2310.15248 showed empirically that the early-type lensing RAR depends on the isolation radius and the late-type one does not;
contamination by type against spectroscopic truth (DESI BGS), and the spurious early/late offset in log g_obs it predicts against the
published 0.14–0.19 dex, remain open". It qualifies because it removes the step the two published answers dispute (Brouwer's
equal-contamination assumption). The referee's cautions, built in: the published 0.14-0.19 dex was measured at R < 0.3 Mpc/h70,
where Brouwer judged contamination small, so report the predicted bias as a function of radius, not only one number; and Mistele
attributes part of the difference to early-type mass-to-light, so a partial or null explanation is a plausible outcome: say in
advance what each would mean. Read Brouwer et al. Appendix A and Mistele et al. Sec. 4 / App. C.5 in full first. Do not rerun the
trajectory check.

## Method
1. Fix before computing, in `results-R147.md`: the lens sample (DESI DR1 BGS, 0.05 < z < 0.3, the KiDS stellar-mass bins), the type
   splits (colour u-r or the DESI equivalent, and morphology/Sersic where available, both as the KiDS papers did), the emulated
   photometric redshifts (the KiDS-bright error model: scatter and outlier rate as published, several seeds), both published
   isolation criteria (Brouwer's 3 Mpc/h70 and Mistele's 4 Mpc/h70, 10 % mass ratio), the contamination statistic (fraction with a
   real neighbour within the radius and 500 km/s from spectroscopic redshifts), and the ESD bias model (each contaminating
   neighbour's halo, from the group VAC halo mass, added to the lens's own profile; converted to g_obs with the deprojection the KiDS
   papers used, at their radial bins).
2. Compute per type, per mass bin, per criterion: the contamination fraction and the predicted log g_obs bias as a function of radius
   (and g_bar); then the predicted spurious ETG-LTG offset averaged exactly as Brouwer averaged the published 0.14 / 0.19 dex.
3. Intervals from seeds, bootstrap over lenses and the halo-mass uncertainty; the stated check that true (spectroscopic) isolation
   gives zero bias; the DESI-vs-KiDS-bright depth difference as a stated systematic.
4. The claim is the table: per criterion, the predicted spurious offset against the published offset, as a fraction with its
   interval, and the radius dependence. No sentence about MOND or dark matter beyond what the fraction says.
5. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus agent,
   section alone; hints: the papers above and anyone who quantified isolation contamination in galaxy-galaxy lensing by type).
   Verdict verbatim in `results-R147.md` and under the section. One revision pass on a presentation point; park on a physics point.
   Nothing on the site; a pass is the headline of your report.
6. End: `### cosmo-deep-16` in `LOG.md`. Final report: headline, the table, nulls, friction. Budget 5 h.
