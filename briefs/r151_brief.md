# cosmo-deep-20: is thin-stream density power at the epicyclic scale or broadband? A cross-stream test of the epicycle and subhalo readings (R-151)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff), `research/cosmo-discovery/blind-review.md` (the gate), `finding.md` (passed
sections: samples and statistic fixed before the fit, the claim its interval), `open-problems/round2/p8-streams.md` (the proposal as
the referee saw it) and `open-problems/round2/verdicts.md` row 8. Follow the agents' working rules: worktree on staging; keep-alive (a background loop
running the session keep-alive every 20 min); land after every milestone. Rows in `results-R151.md`. Long
jobs with `nohup python -u` and an `until` watcher; foreground waits cap at 10 min, bare `sleep` is blocked, match processes with
`ps -eo pid,args`. `data/` starts empty: report fetch times.

## Trajectory check: already run (2026-11-01, context-free Opus, 2 web searches and 9 fetches), verbatim:
"VERDICT: partly novel: joint stream density-power subhalo abundance already done for GD-1+Pal 5 (Banik et al. 1911.02662,
n_sub/n_CDM=0.4+0.3-0.2); open: extending it to DES/S5 streams (ATLAS-Aliqa Uma, Phoenix, Jhelum, Indus, Turranburra) and a
cross-stream epicyclic-scale vs broadband test of the Ibata (2002.01488) dispute". Only the second part qualifies here: it settles a
published disagreement (Ibata et al. 2002.01488: GD-1's power is epicyclic; Banik et al.: subhalos) by a test neither side ran.
Extending the abundance to more streams on its own is the same question on new data and is not this job's claim. Read 2002.01488,
1911.02662 and 2603.09305 (Indus epicycles) in full first. Do not rerun the trajectory check.

## Step 0: feasibility gate (first 90 minutes; stop and report if it fails)
For each candidate stream (ATLAS-Aliqa Uma, Phoenix, Jhelum, Indus, Turranburra, GD-1), compute from its published orbit the predicted
epicyclic density-oscillation scale, and from public photometry (DES DR2 / DELVE, Patrick et al. 2206.04241 density models, Gaia DR3
proper motions) the member count and the 1D density noise. Forecast, with galpy streampepperdf mocks at two subhalo abundances and an
epicycle-only model (and the bar / GMC / Sagittarius perturbers where they apply, as the referee warns), the significance with which
the joint data separate "power concentrated at each stream's own epicyclic scale" from "broadband power". If the forecast separation
is below 2.5 sigma, stop and report: that is the answer.

## Method
1. Fix before measuring any power spectrum, in `results-R151.md`: the streams and their member selection (matched filter plus proper
   motions, with the survey selection maps of 2510.07511 for depth), the 1D density estimator, the power-spectrum estimator and
   scales, the epicyclic-scale prediction per stream (from its orbit and progenitor mass), the statistic (the fraction of each stream's
   excess power within a fixed band around its own epicyclic scale, combined across streams), and the two model families (epicycle
   only; epicycle plus subhalos at the CDM-relative abundance as a free parameter), each forward-modelled with the same selection.
2. Nulls: the same filter on off-stream fields (scanning and depth patterns); shuffled stream assignments of the epicyclic scales
   (the statistic must lose its signal when each stream gets another stream's scale).
3. The claim is the table: per stream and combined, the epicyclic-band fraction with its interval, and its likelihood under each
   model family; the subhalo abundance only as a secondary row.
4. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus agent,
   section alone; hints: the works above and anyone who tested epicyclic versus subhalo origins of stream density structure). Verdict
   verbatim in `results-R151.md` and under the section. One revision pass on a presentation point; park on a physics point. Nothing
   on the site.
5. End: `### cosmo-deep-20` in `LOG.md`. Final report: headline, feasibility numbers, the table, nulls,
   friction. Budget 6 h.
