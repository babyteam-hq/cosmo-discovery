# cosmo-deep-16c: calibrate the contamination amplitude on mocks; headline the type-dependent leakage (R-149)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff, `### cosmo-deep-16`), `research/cosmo-discovery/blind-review.md`,
`results-R147.md` (all rows), the R-147 section of `finding.md`, R-147's scripts. Follow the agents' working rules: worktree on staging; keep-alive
(a background loop running the session keep-alive every 20 min); land after every milestone. Rows continue
in `results-R147.md`. Long jobs with `nohup python -u` and an `until` watcher; foreground waits cap at 10 min, bare `sleep` is blocked,
match processes with `ps -eo pid,args`. `data/` starts empty: report fetch times.

## Why this job
R-147's gate, verbatim: "VERDICT: novel but not sound: the spurious-offset fractions depend on an uncalibrated NFW contamination
amplitude that exceeds Brouwer's MICE benchmark by a factor of about 3.5 (1.05 against 0.30 of the profile), divided by a reference
profile that is the same for both types, so the claimed share of the published offset is not supported." Novelty held. What stands
(measured against spectroscopic truth): photo-z isolation leaks unevenly by type (red 30 % / 43 % vs blue 20 % / 31 % within 0.3 /
0.5 Mpc at fixed M*; hidden group satellites 28 % vs 12 %), so Brouwer's equal-contamination assumption fails. This is the one
revision: the method point (the amplitude) is computable on mocks.

## Method
1. Fix before computing, in `results-R147.md`: a public mock with halos, positions, colours and stellar masses close to DESI BGS
   (e.g. AbacusSummit DESI BGS mocks with halo catalogues, or the MICE catalogue Brouwer used if reachable; state which and why), the
   same emulated photo-z isolation, and the contamination ESD computed directly from the mock's halos (no assumed NFW amplitude):
   per type, the stacked ESD of photo-z-isolated lenses divided by that of truly isolated lenses of the same type (not a common
   reference profile, the referee's second point).
2. Check the mock against the data: the mock's type-dependent leakage fractions against R-147's DESI measurement; if they differ,
   rescale with the measured fractions and state it.
3. The headline is the measured type-dependent leakage with its interval (the part that stands). The predicted spurious ETG-LTG
   offset, now from the calibrated amplitude and per-type references, is the second row, as a fraction of the published value with
   its interval and its radius dependence.
4. Revise the section, then the finding check again (after-writing prompt verbatim, fresh context-free Opus agent, section alone;
   hints as R-147's). Last pass: park on any further physics point. Nothing on the site; a pass is the headline of your report.
5. End: append to `### cosmo-deep-16`. Final report: headline, the mock used, the table, friction. Budget 5 h.
