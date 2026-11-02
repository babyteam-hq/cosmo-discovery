# cosmo-deep-27: do radio quasar polarisation angles align at 150-500 Mpc? An independent, Faraday-derotated re-test with SPICE-RACS (R-160)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff), `research/cosmo-discovery/blind-review.md` (the gate, including the method check
1b), `finding.md` (passed sections; the Giant Arc section passed with this shape: a disputed large-scale claim re-tested on independent
data with a frozen template and nulls; and the withdrawn anisotropy notes, the lesson of sky-direction claims: pre-register the axes,
scales and statistic, correct for look-elsewhere), and `open-problems/round6/s3-pol-align.md` (the proposal as the referee saw it).
Follow the agents' working rules: worktree on staging; keep-alive (a background loop running the session keep-alive every
20 min); land after every milestone. Rows in `results-R160.md`. Long jobs with `nohup python -u` and an `until` watcher; foreground
waits cap at 10 min, bare `sleep` is blocked, match processes with `ps -eo pid,args`. `data/` starts empty: report fetch times.

## Trajectory check: already run (2026-11-02, context-free Opus, 5 web searches, 5 API queries, 3 fetches), verbatim:
"VERDICT: partly novel: whether radio quasar polarisation angles align on large scales is already answered with conflicting results on
JVAS/CLASS 8.4 GHz data (alignment at ~150 Mpc: 1201.5180, 1503.03482; none: 0705.2548), and a Faraday-insensitive dual-frequency
test finds isotropy only at ≥800 Mpc (1809.01270); still open is a test at 150-500 Mpc with an independent instrument and broadband
Faraday-derotated angles (SPICE-RACS), plus a check of the same sources against JVAS". It settles a published disagreement by removing
the disputed step (JVAS instrumental-polarisation residuals, 0705.2548). Read 0705.2548, 1201.5180, 1503.03482, 1604.03937,
1809.01270 and 2002.00982 in full first. Do not rerun the trajectory check.

## Step 0: feasibility gate (first 90 minutes; stop and report if it fails)
Fetch SPICE-RACS DR2 (arXiv:2605.16917; CASDA) and match to spectroscopic quasars (SDSS DR16Q, DESI DR1 QSO) and to JVAS/CLASS.
Count matched quasars, and separately the flat-spectrum, core-dominated subset (the population the JVAS claim is about; spectral index
from RACS-low/RACS-mid/NVSS); count the pairs at 150-500 Mpc comoving. Forecast the alignment statistic's power against the published
JVAS alignment amplitude (anchor: 1201.5180 / 1503.03482). If the flat-spectrum subset cannot test the published amplitude at 3 sigma,
stop and report.

## Step 1: method check (blind-review.md 1b), before any statistic is computed
The fixed plan (samples, the Faraday derotation (the catalogue's intrinsic angle chi_0 with its error), the alignment statistic
(Jain et al. / Hutsemekers S_D with parallel transport, at fixed scales 150, 300, 500 Mpc), the instrumental nulls (angles against
position in the ASKAP beam footprint and the PAF tile pattern; a randomised-angle null that keeps positions; shuffled redshifts), the
JVAS same-source check (angle differences SPICE-RACS minus JVAS for matched flat-spectrum sources), the look-elsewhere, what each
reading predicts) alone in a file; the method-check prompt verbatim in a fresh context-free Opus agent. Close every point before
computing, or stop. Verdict verbatim in `results-R160.md`.

## Method
1. The alignment statistic and its significance from position-preserving angle randomisations, per scale, for the flat-spectrum
   subset (headline) and all matched quasars (check); the same on the JVAS angles for the same sources.
2. The claim is the table: per scale, the alignment statistic and its significance, the amplitude bound against the published JVAS
   amplitude, and the same-source angle comparison. A null is the result if the data say so. No sentence about cosmic birefringence
   or new physics beyond what the numbers say.
3. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus agent,
   section alone; hints: the works above). Verdict verbatim in `results-R160.md` and under the section. One revision pass on a
   presentation point; park on a physics point. Nothing on the site; a pass is the headline of your report.
4. End: `### cosmo-deep-27` in `LOG.md`. Final report: headline, feasibility numbers, method-check verdict,
   the table, nulls, friction. Budget 6 h.
