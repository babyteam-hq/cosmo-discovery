# cosmo-deep-29: does the leading-tail excess of nearby open clusters survive radial-velocity membership? MOND vs Newton (dispatch when Gaia DR4 is public)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Your Request number is the R-N in your prompt; use it for the keep-alive and commit messages.
Read first: `LOG.md` (handoff), `research/cosmo-discovery/blind-review.md` (the gate, including the method check
1b), `finding.md` (passed sections; the Giant Arc passed by re-testing a claim on independent data with a frozen template), and
`open-problems/round8/g1-tidal-tails.md` (the proposal as the referee saw it). Follow the agents' working rules: worktree on staging; keep-alive (a
background loop running the session keep-alive every 20 min); land after every milestone. Rows in
`results-tidal-tails.md`. Long jobs with `nohup python -u` and an `until` watcher; foreground waits cap at 10 min, bare `sleep` is
blocked, match processes with `ps -eo pid,args`. Work fast: the Kroupa group will redo this on DR4.

## Trajectory check: already run (2026-11-02, context-free Opus, 5 web searches and 5 fetches), verbatim:
"VERDICT: partly novel: Hyades tails mapped with RV-informed membership already show a weaker trailing tail (2012.02061, also
2007.02969), and CP-based leading/trailing counts exist (2405.09609); a completeness-matched, RV-confirmed leading/trailing count across
clusters as a MOND-vs-Newton test remains open". It settles a published disagreement (2210.13472 and 2405.09609, ~9 sigma
leading-tail excess from convergent-point membership, read as MOND; Newtonian chance asymmetry 2301.02251; the membership caveats of
2607.13747) by removing the disputed step: membership without radial velocities. Read 2405.09609, 2301.02251, 2012.02061 and
2607.13747 in full first. Do not rerun the trajectory check.

## Step 0: feasibility gate (first hour)
Confirm Gaia DR4 is public; for the clusters of 2405.09609 (Hyades, Praesepe, Coma Ber, NGC 752, COIN-Gaia13 as available), count the
convergent-point tail candidates with DR4 RVS velocities, the completeness limit (G_RVS) common to both tails, and the members
surviving an RV cut; forecast the separation between the published asymmetry (0.26) and Newtonian symmetry with the chance spread of
2301.02251. Name the moving-group confound (field stars sharing the cluster velocity) and estimate it from off-tail control windows at
the same velocity. If below 3 sigma, stop and report.

## Step 1: method check (blind-review.md 1b), before counting tails
The fixed plan (clusters, the tail windows frozen from 2405.09609, the membership (convergent point plus RV within a stated window of
the cluster's predicted tail velocity), the common completeness limit, the asymmetry statistic, the Newtonian null from N-body or the
2301.02251 distribution, the moving-group control, what MOND and Newton each predict) alone in a file; the method-check prompt verbatim
in a fresh context-free Opus agent. Close every point before counting, or stop.

## Method
1. RV-confirmed leading and trailing counts per cluster to the common limit; the asymmetry per cluster and combined, with intervals.
2. Against the Newtonian chance distribution and the published MOND-like asymmetry, as a likelihood ratio; the moving-group control
   subtracted with its error.
3. The claim is the table. Write the `finding.md` section, gate (finding check verbatim), record; one revision on a presentation point;
   park on physics. Nothing on the site; a pass is the headline.
4. End: `### cosmo-deep-29`. Budget 5 h.
