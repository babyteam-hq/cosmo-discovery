# cosmo-deep-28: do the four quasar jet alignment regions stay aligned on new VLBI sources? (a re-test of arXiv:2106.06546)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Your Request number is the R-N in your prompt; use it for the keep-alive and commit messages.
Read first: `LOG.md` (handoff), `research/cosmo-discovery/blind-review.md` (the gate, including the method check
1b), `finding.md` (passed sections; the Giant Arc section passed with this exact shape: a claim re-tested on independent data against
a template frozen from the claim paper, with nulls; and the withdrawn anisotropy notes: pre-register directions and scales), and
`open-problems/round7/t1-jet-alignment.md` (the proposal as the referee saw it). Follow the agents' working rules: worktree on staging; keep-alive (a
background loop running the session keep-alive every 20 min); land after every milestone. Rows in
`results-jet-alignment.md`. Long jobs with `nohup python -u` and an `until` watcher; foreground waits cap at 10 min, bare `sleep` is
blocked, match processes with `ps -eo pid,args`. `data/` starts empty: report fetch times.

## Trajectory check: already run (2026-11-02, context-free Opus, 3 web searches and 3 fetches), verbatim:
"VERDICT: novel: no prior work answers this question". The claim: Mandarakas et al. 2021 (arXiv:2106.06546, A&A 653, A123), four
(RA, Dec, z) regions of Astrogeo parsec-scale jet position angles deviating from random at > 5 sigma, found by a grid search and
optimisation. Related: 2002.00982 (no global alignment), 2203.13750 (Plavin et al. automated jet directions for 9,220 AGN, no regional
test), 2411.09099, 2411.16951. Read 2106.06546 and 2203.13750 in full first. Do not rerun the trajectory check.

## Step 0: feasibility gate (first hour)
Freeze from 2106.06546 the four regions (centres, radii in angle and redshift) and their published mean PAs and spreads. Build the
independent sample: VLBI jet PAs (Plavin et al. directions and the current Astrogeo / rfc releases) for sources NOT in the 2021
sample (imaged after its cutoff, or newly given a redshift by DESI DR1 / SDSS), inside each frozen region. Count them; forecast the
separation between "alignment real at the published strength" and "random", and at half strength. If below 3 sigma, stop and report.

## Step 1: method check (blind-review.md 1b), before any PA statistic is computed on the new sample
The fixed plan (the frozen regions and PAs, the independent-sample definition, the statistic (the Rayleigh / dispersion statistic of
2106.06546 and the mean-PA difference from the published value, on axial angles), the beam-elongation controls (PA against the
restoring-beam position angle per source; the same statistic on sources outside the regions at matched declination; beam-PA-matched
random draws), the look-elsewhere (none needed for the frozen regions; state it), what each reading predicts) alone in a file; the
method-check prompt verbatim in a fresh context-free Opus agent. Close every point before computing, or stop.

## Method
1. Per frozen region: the alignment statistic and the mean PA on the new sample, the significance from beam-matched random draws, and
   the difference from the published PA; combined over the four regions.
2. The claim is the table: per region, N new, the statistic, the significance, the PA difference with its interval; whether the
   published alignment survives on independent sources. A null is the result if the data say so.
3. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus agent,
   section alone; hints: the works above). Verdict verbatim in the results file and under the section. One revision pass on a
   presentation point; park on a physics point. Nothing on the site; a pass is the headline of your report.
4. End: `### cosmo-deep-28` in `LOG.md`. Final report: headline, feasibility numbers, method-check verdict, the
   table, nulls, friction. Budget 5 h.
