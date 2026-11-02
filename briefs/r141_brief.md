# cosmo-deep-10b: the lensed-SN rate limit with injections through the full vetting, and the forecast-free limit as the headline (R-141)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff at the top, and `### cosmo-deep-10`), `research/cosmo-discovery/blind-review.md`
(the gate), `results-R130.md` (all rows, including R-135's), the lensed-SN section of `finding.md`, and the scripts R-135 used
(`ztf_lsn_*.py`). Follow the agents' working rules: worktree on staging; keep-alive (a background loop running
the session keep-alive every 20 min); land after every milestone with
`git pull --rebase origin staging && git push origin HEAD:staging`. Rows in `results-R130.md` (continue it). Long jobs with
`nohup python -u` into a log file and an `until` watcher; foreground waits are capped at 10 min and bare `sleep` is blocked; match
processes with `ps -eo pid,args` and an anchored pattern, never `pkill -f`.

**Data:** `research/cosmo-discovery/data/` is empty in a new container. R-135 refetched for ~3 h. Shard the IRSA fetch over 3-4
fetchers on disjoint contiguous rank ranges (15-22 files/s in total) and refetch only the 168 usable lenses and the 704 control
light curves; say how long it took.

## Why this job
R-135's finding check, verbatim: "VERDICT: novel but not sound: the quoted efficiencies come from injections that skip the colour,
magnification and variability vetting that every real candidate had to pass, and the class-L limits that carry the headline come
from an assumed source population rather than from the data." Novelty held. Both points are fixable without new data; this is the
one revision.

## Method
1. Fix before rerunning, in `results-R130.md`: the vetting as code (the exact colour, magnification-consistency, variability and
   shape cuts R-135 applied to real candidates, written as a function), applied identically to every injected light curve; the
   injection set unchanged otherwise.
2. Recompute the efficiency per class and per peak magnitude with the vetting in; recompute the limits.
3. The headline is the forecast-free limit (lensed SNe with summed peak r < 21 per lens-year, 95 %, per class), with the efficiency
   now including the vetting. The forecast comparisons (Shu 2018 per-lens, the SLACS- and BELLS-like brackets) become secondary
   rows, each labelled with the source population it assumes.
4. Revise the section in `finding.md`, then the finding check again (after-writing prompt verbatim, fresh context-free Opus agent,
   section alone; hints as R-135's). Verdict verbatim in `results-R130.md` and under the section. This is the last pass: park on any
   further physics point. Nothing on the site; a pass is the headline of your report.
4. End: append to `### cosmo-deep-10` in `LOG.md`. Final report: headline, the vetted efficiencies, the
   limits, the refetch time, friction. Budget 5 h; stop at the budget with everything landed.
