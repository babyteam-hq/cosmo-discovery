# cosmo-deep-23: does the Giant Arc reappear on independent sightlines? DESI DR1 Mg II absorbers inside the arc (R-155)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff), `research/cosmo-discovery/blind-review.md` (the gate, including the method check
1b), `finding.md` (passed sections: samples and statistic fixed before the fit, the claim its interval), and
`open-problems/round3/q5-giant-arc.md` (the proposal as the referee saw it: data URLs, the four-line method, the overlap and null
numbers). Follow the agents' working rules: worktree on staging; keep-alive (a background loop running the session keep-alive
every 20 min); land after every milestone. Rows in `results-R155.md`. Long jobs with `nohup python -u` and an `until`
watcher; foreground waits cap at 10 min, bare `sleep` is blocked, match processes with `ps -eo pid,args`. `data/` starts empty: report
fetch times.

## Trajectory check: already run (2026-11-01, context-free Opus, 2 web searches, 4 API queries, 3 fetches), verbatim:
"VERDICT: novel: no prior work answers this question". The dispute: Lopez, Clowes & Williger (arXiv:2201.06875; Big Ring 2402.07591;
Giant Ring 2604.17534; FLAMINGO-10K 2504.14940) against Sawala et al. (arXiv:2502.03515) and Sawala & Teeriaho (arXiv:2505.11072),
both sides arguing from simulations of the same SDSS points. Read all of them first. Do not rerun the trajectory check.

**Disclosure that must go into the finding:** the harvester that proposed this ran an unblinded rough count before the analysis was
fixed (delta = -0.16 / -0.28 / -0.29 at R = 47.5 / 70 / 95 Mpc, with partial SDSS exclusion, no EW matching, no TSNR2 weighting). The
analysis below is fixed as the proposal states it, before your own first look; say in the finding that a rough count was seen at
proposal time, and that nothing in the plan was changed after it.

## Step 0: feasibility (first hour)
Reproduce the proposal's overlap numbers (sightlines through the template, expected counts, null scatter) from the public files; stop
and report if the separation between "real" (delta ~ 1.3) and "chance" (delta ~ 0.07) falls below 3 sigma.

## Step 1: method check (blind-review.md 1b), before the measurement
The fixed plan (the template rebuild and freeze, the sightline selection and SDSS exclusion, the path computation and proximate cut,
dN/dz on the same sightlines, the EW matching and TSNR2 weighting, the two nulls and the subtraction of the clustering baseline, the
completeness treatment, what delta value each reading predicts) alone in a file; the method-check prompt verbatim in a fresh
context-free Opus agent. Close every point it names before measuring, or stop. Verdict verbatim in `results-R155.md`.

## Method (as the proposal fixes it)
1. Rebuild the Giant Arc from the Zhu & Ménard DR7+DR12 trimmed catalogues with the discovery recipe (z = 0.802 +- 0.060, single
   linkage 95 Mpc, Om 0.27, h 0.7); freeze members and the Big Ring / Giant Ring templates before touching DESI.
2. DESI DR1 sightlines not in SDSS DR16Q (1.5" match); comoving path inside the template (R = 47.5, 70, 95 Mpc, and the CHMS hull),
   proximate zone removed; expected count from dN/dz on the same sightlines over the whole sky.
3. delta_GA = N_obs/N_exp - 1 at matched EW cuts (W2796 > 0.3 and > 0.6 A), TSNR2-weighted.
4. Nulls: template shifts across the DESI footprint (coverage variance); spheres on random sets of 46 SDSS slice absorbers (the
   ordinary-clustering baseline). The claim is delta_GA minus that baseline, with its interval, per R and EW cut; the same for the
   Big Ring and Giant Ring where DR1 coverage allows (check coverage first).
5. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus agent,
   section alone; hints: the works above). Verdict verbatim in `results-R155.md` and under the section. One revision pass on a
   presentation point; park on a physics point. Nothing on the site; a pass is the headline of your report.
6. End: `### cosmo-deep-23` in `LOG.md`. Final report: headline, method-check verdict, the table, nulls,
   friction. Budget 5 h.
