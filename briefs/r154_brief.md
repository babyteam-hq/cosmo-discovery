# cosmo-deep-22: does the SZ mass bias keep falling at z > 1? CMB-lensing M_lens/M_SZ of ACT DR6 clusters in redshift bins (R-154)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff), `research/cosmo-discovery/blind-review.md` (the gate, including the method check
1b), `finding.md` (the eRASS1 x ACT DR6 lensing section passed the gate with this exact machinery: samples and statistic fixed before
the fit, the claim its interval, 400-simulation covariance), `results-R122.md`, and `open-problems/round3/q1-sz-lensing-z1.md` with
its verdict in `open-problems/round3/verdicts.md`. Reuse `erass1_act.py` (ACT DR6 lensing, masks, the 400 simulations, `simstream.sh`,
the matched-filter stack). Follow the agents' working rules: worktree on staging; keep-alive (a background loop running the session keep-alive
every 20 min); land after every milestone. Rows in `results-R154.md`. Long jobs with `nohup python -u` and an
`until` watcher; foreground waits cap at 10 min, bare `sleep` is blocked, match processes with `ps -eo pid,args`. `data/` starts
empty: R-139 streamed the 400 ACT simulations in ~85 min.

## Trajectory check: already run (2026-11-01, context-free Opus, 2 web searches and 11 fetches), verbatim:
"VERDICT: partly novel: CMB-lensing detection and rough mass normalisation of ACT SZ-selected clusters, including a separate 4.3σ z>1
stack, already exist (arXiv:2406.14754, ACT-DR5 MCMF x DR6 lensing), and Planck SZ clusters have a single CMB-lensing
1-b=0.71±0.10 (arXiv:1904.07887); what remains open is a CMB-lensing M_lens/M_SZ in redshift bins, and its evolution exponent, to
test the falling-with-redshift bias of arXiv:2512.18935 at z>1". It qualifies as a test of a published claim (galaxy-lensing 1-b
falling from 0.81 to 0.58 with z, arXiv:2512.18935; milder in HSC, arXiv:2407.08201; SPT's evolution consistent, arXiv:2401.02075)
by a method free of galaxy-lensing source limits at z > 1. Read 2406.14754 (its estimator and its 4.3 sigma z > 1 stack),
2512.18935 and 2401.02075 in full first. Scoop risk is high: work fast. Do not rerun the trajectory check.

## Step 0: feasibility gate (first hour; stop and report if it fails)
Count the ACT DR6 SZ clusters (the DR6 cluster catalogue if public, else DR5) per redshift bin inside the DR6 lensing mask; forecast,
with R-122's matched filter and the 400-simulation noise, the error on M_lens/M_SZ per bin, anchored to 2406.14754's 4.3 sigma at
z > 1 (scale from it, not from MaDCoWS). Then the separation between "1-b keeps falling as in 2512.18935" and "constant 1-b" from the
redshift trend. If it is below 3 sigma, stop and report: that is the answer.

## Step 1: method check (blind-review.md 1b), before any stack
The fixed plan (sample, bins, filter and profile, miscentring and the SZ-selection (Eddington) bias treatment, the two-halo term, the
tSZ-leakage treatment in the lensing reconstruction, the statistic and the trend fit, the nulls) alone in a file, the method-check
prompt verbatim in a fresh context-free Opus agent. Close every point it names before computing, or stop. Verdict verbatim in
`results-R154.md`.

## Method
1. Stack ACT DR6 kappa on the clusters per redshift bin with R-122's matched filter at each cluster's M_SZ; fit M_lens/M_SZ per bin
   and the evolution exponent; tSZ leakage bounded with the profile-hardened / CIB-deprojected lensing variants; the SZ selection's
   Eddington bias modelled from the catalogue's own completeness.
2. Nulls: random positions, rotated maps, curl; covariance from the 400 simulations.
3. The claim is the table: per bin, M_lens/M_SZ (and 1-b) with its interval; the trend exponent; the comparison with 2512.18935's
   bins and SPT. No sentence about S8 beyond what the trend says.
4. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus agent,
   section alone; hints: the works above). Verdict verbatim in `results-R154.md` and under the section. One revision pass on a
   presentation point; park on a physics point. Nothing on the site; a pass is the headline of your report.
5. End: `### cosmo-deep-22` in `LOG.md`. Final report: headline, feasibility numbers, method-check verdict, the
   table, nulls, friction. Budget 6 h.
