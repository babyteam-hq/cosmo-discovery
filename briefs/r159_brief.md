# cosmo-deep-26: the host-halo mass of Euclid-selected little red dots at 1.5 < z < 3.5 from ACT DR6 CMB lensing (dispatch when Euclid DR1-Foundation is public; your Request number is the R-N in your prompt, not 159)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff), `research/cosmo-discovery/blind-review.md` (the gate, including the method check
1b), `finding.md` (passed sections; the eRASS1 x ACT DR6 lensing section passed with this machinery), `results-R122.md` and
`results-R139.md` (CMB lensing of a high-z sample: magnification, interlopers and sky-leg consistency were where R-139 and R-146 failed
the gate), and `open-problems/round5/e3-lrd-lensing.md` (the proposal as the referee saw it). Reuse `erass1_act.py` (ACT DR6 lensing,
masks, the 400 simulations). Follow the agents' working rules: worktree on staging; keep-alive (a background loop running the session keep-alive
every 20 min); land after every milestone. Rows in `results-euclid-lrd.md`. Long jobs with `nohup python -u` and an
`until` watcher; foreground waits cap at 10 min, bare `sleep` is blocked, match processes with `ps -eo pid,args`.

## Trajectory check: already run (2026-11-02, context-free Opus, 4 web searches and 6 API calls), verbatim:
"VERDICT: novel: no prior work answers this question". The referee's comparison set, built in: JWST LRD clustering halo masses
(arXiv:2606.02773: ~10^11.3 Msun at z~3.5; 2606.09721 NEXUS; 2610.23806) against obscured/red quasars by CMB lensing and clustering
(arXiv:1411.0527 b 2.57+-0.24 obscured vs 1.89+-0.19 unobscured at z~1; 1511.04469; 2302.00690; 1705.05306; 2409.08314), scaled to
z ~ 2-3 before comparing. Euclid Q1 candidates: arXiv:2503.15323. Do not rerun the trajectory check.

## Step 0: feasibility gate (first hour; stop and report if it fails)
Confirm Euclid DR1-Foundation is public and fetch its catalogues; apply the Q1 LRD colour/compactness selection (2503.15323) at
1.5 < z < 3.5; count candidates inside the ACT DR6 lensing mask; bound stellar (brown dwarf) and low-z interloper contamination from
the DR1 spectra and Gaia; measure the magnification slope s from the candidates' own counts; forecast the cross-power significance and
the separation between the JWST-LRD halo mass (~10^11.3) and the obscured-quasar mass (~10^12.7-13.2) with magnification included.
If below 3 sigma, stop and report.

## Step 1: method check (blind-review.md 1b), before any cross-spectrum
The fixed plan (selection, contamination model and its propagation, dN/dz from photo-z and the DR1 spectra, the magnification term
with the measured slope, the multipole range, the estimator, the bias-to-mass conversion, nulls, ACT vs Planck legs and the error
inflation if they disagree) alone in a file; the method-check prompt verbatim in a fresh context-free Opus agent. Close every point
before computing, or stop. Verdict verbatim in `results-euclid-lrd.md`.

## Method
1. kappa x LRD-candidate cross-power on ACT DR6 (Planck PR4 as the check), plus the candidates' angular auto-correlation; joint fit
   of b with the magnification term and the contamination fraction marginalised.
2. Nulls: curl, randoms, stars, a matched sample of ordinary galaxies at the same z (must reproduce its known bias).
3. The claim: b and the host-halo mass with intervals, and which of the two published reference masses it favours (as a likelihood
   ratio). No sentence about LRD physics beyond what the mass says.
4. Write, gate (finding check verbatim), record; one revision on a presentation point; park on physics. Nothing on the site.
5. End: `### cosmo-deep-26`. Budget 6 h.
