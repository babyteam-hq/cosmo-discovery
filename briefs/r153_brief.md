# cosmo-deep-21: an axion-like-particle decay-line search in SPHEREx images of dwarf spheroidals (R-153)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff), `research/cosmo-discovery/blind-review.md` (the gate, including the new method
check, 1b), `finding.md` (passed sections: samples and statistic fixed before the fit, the claim its interval),
`open-problems/round3/q4-alp-spherex.md` (the proposal as the referee saw it), and `results-R138.md` (R-138 and R-143's SPHEREx QR2
handling: public S3 L2 images, 72 MB each, ZODI and per-pixel wavelength included, ~6,500 images/h reduction on 4 cores; reuse the
code). Follow the agents' working rules: worktree on staging; keep-alive (a background loop running the session keep-alive
every 20 min); land after every milestone. Rows in `results-R153.md`. Long jobs with `nohup python -u` and an `until` watcher;
foreground waits cap at 10 min, bare `sleep` is blocked, match processes with `ps -eo pid,args`. `data/` starts empty: report fetch
times.

## Trajectory check: already run (2026-11-01, context-free Opus, 2 web searches and 4 API fetches), verbatim:
"VERDICT: partly novel: whether SPHEREx dSph data can reach g_agamma of a few 1e-12 over 0.5-3 eV is already forecast
(arXiv:2412.12286), and the JWST order-of-magnitude disagreement is already set out and traced to instrument modelling
(arXiv:2609.10674, vs 2503.11753/2310.15395/2503.14582); what remains open is an actual decay-line search in SPHEREx Quick Release
data with same-exposure ON-OFF, which no paper has done, as a JWST-independent arbiter at about 1-3 eV". It qualifies twice: posed
(forecast) but never measured on data, and it settles a published disagreement by removing the disputed step (JWST instrument
modelling). Read 2412.12286 (method, D-factors, the targets) and 2609.10674 (why the JWST limits disagree) in full first. Caveats from
the referee: MUSE-Faint (2307.07403, 2310.01996) already limits ~2.65-5.3 eV, so the new range is ~1-2.65 eV; below ~1 eV the reach is
weak; the headline reach is a photon-noise forecast, unshown on real data. Do not rerun the trajectory check.

## Step 0: feasibility gate (first 90 minutes; stop and report if it fails)
Fetch the QR2 images covering Draco, Ursa Minor and the other dSphs in the proposal; measure, from the real images, the noise of the
same-exposure ON-minus-OFF spectrum per spectral resolution element against the photon-noise forecast, and the residual from the
dwarfs' own starlight (the stellar continuum and its absorption features through the same ON-OFF). Convert to the expected 95 %
limit on g_agamma per mass with the published D-factors. If the expected limit at 1.5-2.5 eV does not reach below the weaker of the
two disputed JWST limits (so cannot arbitrate), stop and report.

## Step 1: method check (blind-review.md 1b), before any limit is computed
Write the fixed plan (targets, apertures, the ON/OFF definition, the stellar-light model, the line profile including the dSph
velocity dispersion and SPHEREx's spectral response, the D-factors and their uncertainties, the statistic and how the limit is set,
the nulls) alone in a file and run the method-check prompt verbatim in a fresh context-free Opus agent. Close every point it names
before computing, or stop and report. Verdict verbatim in `results-R153.md`.

## Method
1. Same-exposure ON (within the dSph's half-light or D-factor aperture) minus OFF (an annulus in the same exposure) per spectral
   channel, stacked over exposures and dwarfs with D-factor weights; the stellar light modelled from the dwarfs' own photometry and a
   stated old-population spectrum, with its uncertainty propagated.
2. Line search across ~0.5-3 eV with the line width set by the velocity dispersion and the SPHEREx resolution; the limit by a
   profile likelihood with the D-factor uncertainty marginalised.
3. Nulls: the same pipeline on blank fields with no dark matter (fake dwarfs), injected lines recovered at known coupling, the line
   statistic across neighbouring channels (look-elsewhere stated).
4. The claim is the limit curve g_agamma(m_a) with its interval, against the JWST limits and MUSE-Faint; any excess reported with its
   global significance.
5. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus agent,
   section alone; hints: the works above and anyone who searched for eV dark-matter decay lines). Verdict verbatim in
   `results-R153.md` and under the section. One revision pass on a presentation point; park on a physics point. Nothing on the site.
6. End: `### cosmo-deep-21` in `LOG.md`. Final report: headline, feasibility numbers, the method-check verdict,
   the limit curve, nulls, friction. Budget 6 h.
