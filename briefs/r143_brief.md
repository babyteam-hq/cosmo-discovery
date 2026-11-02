# cosmo-deep-14b: the SPHEREx x DESI tomography with the source-mask transfer measured by injection (R-143)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff, and `### cosmo-deep-14`), `research/cosmo-discovery/blind-review.md`,
`results-R138.md` (all rows and amendments), the SPHEREx section of `finding.md`, and R-138's scripts (`spherex_*.py`). Follow
the agents' working rules: worktree on staging; keep-alive (a background loop running the session keep-alive every 20
min); land after every milestone. Rows continue in `results-R138.md`. Long jobs with `nohup python -u` and an `until` watcher;
foreground waits cap at 10 min, bare `sleep` is blocked, match processes with `ps -eo pid,args`. `data/` is empty in a new
container: R-138 reduced ~6,500 images/h with 4 cores; refetch and reduce only what the injection test and the rerun need, and say
how long it took.

## Why this job
R-138's finding check, verbatim: "VERDICT: novel but not sound: the amplitude and redshift shape are set by a catalogue source mask
and clipping that remove galaxy-correlated light, which the transfer correction does not model, as the factor-of-two star-density
split and the failed stellar null show, so the headline amplitudes and the CIBER ratio are not supported." Novelty held. The fix
is a measurement, not new data: this is the one revision.

## Method
1. Fix before rerunning, in `results-R138.md`: the injection test. Paint mock galaxy-correlated diffuse light into real SPHEREx
   images (a field built from the DESI galaxies with a stated luminosity-weighting and a smooth unresolved component correlated with
   them, at several amplitudes and redshift shapes), run the identical source mask, clipping and pipeline, and measure the transfer
   (recovered over injected) per channel group and slice. Also: the stellar null redone with the star mask fixed; the star-density
   split predicted by the measured transfer.
2. Correct the measured b_I d(nuI_nu)/dz by the injected transfer; propagate its uncertainty; redo the CIBER comparison at 1.1 and
   1.8 um.
3. The claim is the corrected table with intervals. If the transfer cannot be measured well enough to support amplitudes, the result
   is the detection per channel and the shape with the amplitude stated as a bound: say so.
4. Revise the section in `finding.md`, then the finding check again (after-writing prompt verbatim, fresh context-free Opus agent,
   section alone; hints as R-138's). Verdict verbatim in `results-R138.md` and under the section. Last pass: park on any further
   physics point. Nothing on the site; a pass is the headline of your report.
5. End: append to `### cosmo-deep-14`. Final report: headline, the transfer, the corrected table, the nulls, the CIBER comparison,
   friction. Budget 6 h.
