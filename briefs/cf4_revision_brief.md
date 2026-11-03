# cosmo-cf4-rev: calibrate the null of the Type Ia supernovae vs Cosmicflows-4 test, then one revision (R-167 brief)

Goal (the only one): make the published Cosmicflows-4 paper (`site/research-cf4-supernovae.html`, under revision since
2026-11-03) pass the finding check on a calibrated null, or park it. One revision, then park: no second revision.
Your Request number is the R-N in your prompt; use it for the keep-alive and commit messages.

Read first: `LOG.md` (handoff), `research/cosmo-discovery/blind-review.md` (the gate, including the method
check 1b), `finding.md`'s first section (the CF4 finding, its "To settle it" and its blind reviews), `results.md` (rows "CF4
follow-up 1, 2, 2b", "CF4 axes", "CF4 bulk flow"), the page itself, and the verdict history at the end of
`outreach/authors-cf4-supernovae.md` (seven re-checks, 2026-11-02 to 03, all "novel but not sound", all on the same point). The
scripts are `cf4_*.py`, `ztf_tests.py`, `bulk_flow.py`, `bulkflow_free.py`. Follow the agents' working rules: worktree on staging; keep the key
alive as the worker brief says; land after every milestone. Rows in `results-cf4-rev.md` (own file). Long jobs with
`nohup python -u` and an `until` watcher; foreground waits cap at 10 min, bare `sleep` is blocked. `data/` starts empty: refetch
Pantheon+SH0ES (data release GitHub), the CF4 2023 and CF4++ grids, DES-SN5YR Dovekie, ZTF DR2 lite and the 2M++ field, and report
fetch times.

## What the re-checks found (do not argue with it; close it)
The estimator's null is uncalibrated. `cf4_mocks.py` injects the CF4++ field itself as the truth, so a correct field returns
β ≈ 0.88; but a real reconstruction is a shrunk posterior mean, not the truth, and Pantheon+'s own 2M++ correction, the
best-validated field, returns 0.39 ± 0.18 on the same test. So no β (CF4++ 0.10 to 0.26 non-dipole, 0.33 to 0.81 dipole) can be read
as a deficiency of either field. Also named: the Pantheon+ covariance already carries the 2M++ peculiar-velocity terms; the
Pantheon+ supernovae are Cosmicflows-4 inputs; the CF4 vs 2M++ comparison (1.5σ) treats two fits on the same supernovae as
independent; the bulk-flow direction offset has no uncertainty. Closest prior work, to cite and read in full: Stiskalek et al. 2025
(arXiv:2502.00121), Stiskalek, Desmond & Lavaux 2025 (arXiv:2509.14997), Watkins, Clark & Feldman (arXiv:2512.03168, the CF4 dipole
split), Nusser & Tully 2026 (arXiv:2608.14265), Carrick et al. 2015 (arXiv:1504.04627), Courtois et al. 2023, 2025
(arXiv:2211.16390, 2502.01308).

## Trajectory check: not rerun
The question passed novelty in the finding's third pass (`finding.md`, "Blind review"); the re-checks fault soundness only. Do not rerun the trajectory check.

## Step 0: feasibility gate (first hour)
The fix stands on a CF4-like null: what β does a correct field return through the same estimator, when the field is a
reconstruction of the truth with Cosmicflows-4's errors, sky coverage and shrinkage, not the truth itself? Find which route is open:
1. Public posterior realisations of the CF4 field (CF4++ ensembles or HMC samples: Courtois et al. 2025, Hollinger et al., Valade
   et al.; the CF4 data pages at IP2I and EDD). Each realisation stands in for the truth; the posterior mean is the field under test.
2. Otherwise, a CF4-like mock chain built here: ΛCDM truth velocity fields (constrained or unconstrained Gaussian realisations at
   the CF4 grid, linear theory), a mock CF4 catalogue at the real CF4 galaxy positions with each galaxy's distance error, a Wiener
   filter reconstruction in the same grid, and the posterior mean as the field under test. The same chain for 2M++ (density-based,
   Carrick et al.'s β fit) gives the 2M++ null, which the data's 0.39 must match for the estimator to be trusted.
Forecast for the open route: the calibrated null's spread across realisations, and the separation of CF4++'s observed β from that
null in the sample where it is tightest. If neither route is open within the budget, or the forecast separation is below 2σ in
every sample, stop and report: the paper stays under revision and the default is to withdraw it (Will's call, put to them in chat).

## Step 1: method check (blind-review.md 1b), before computing any calibrated β
The fixed plan alone in a file: the samples (Pantheon+ 618, DES 191 as pipeline check only, ZTF DR2 1,122, and Pantheon+ minus
CF4 inputs from `cf4_notin.py`), the residuals (Pantheon+ refitted with its covariance rebuilt without the 2M++ peculiar-velocity
terms, from the STAT-only matrix and the systematics files; ZTF and DES as released), the statistic (β_dipole, β_rest, joint fit,
GLS), the null (route 1 or 2 above, its number of realisations), the 2M++ field on every sample including ZTF, the paired
CF4 minus 2M++ difference with its error from a joint fit of both fields on the same supernovae and a bootstrap over them, the
bulk-flow direction error from the same bootstrap, and what each reading predicts. Run the method-check prompt verbatim in a fresh
context-free Opus agent. Close every point it names before computing, or stop and report the one you cannot close.

## Method
1. The calibrated null per sample and field: the distribution of β_rest and β_dipole a correct, CF4-like (and 2M++-like)
   reconstruction returns.
2. The data's β per sample and field, as a fraction of its own calibrated null, with the null's spread in the error.
3. The 2M++ check: if 2M++'s observed β is far from its own 2M++-like null, the estimator is not trusted and you stop there: that
   is the result, stated as such.
4. The paired CF4 minus 2M++ difference per sample, ZTF first (the one sample independent of Cosmicflows-4).
5. The claim is the table: per sample and field, observed β, calibrated null, the ratio with its interval, and the paired
   difference; the headline is no stronger than the tightest interval. A null result is a result.
6. Rewrite the finding as a `finding.md` section (a revision of the first section; keep the old one and its verdicts), then the
   finding check (after-writing prompt verbatim, fresh context-free Opus agent, the section alone; hints: the prior works above;
   failure modes: an uncalibrated or self-referential null, shared supernovae treated as independent, covariance double-counting
   2M++, CF4 inputs in the test sample, trials). Verdict verbatim in `results-cf4-rev.md` and under the section.
7. One revision only. A pass: rewrite the page's title, abstract and results from the passed section, drop the status line,
   re-check the revised title and abstract alone with the same prompt, and regenerate `site/research-erass1-act-lensing.html`'s head
   from it only if `briefs/mkpaper.py` still applies (it copies the CF4 head). Anything short of a pass: park; the page keeps its
   "Under revision" line and the report proposes withdrawing it (Will's call). `outreach/authors-cf4-supernovae.md` stays held
   either way until Will decides.
8. End: `### cosmo-cf4-rev` in `LOG.md`, and the methodology page's record (`briefs/paper_method_body.html`,
   regenerated with `briefs/mkpaper_method.py`) updated to the outcome. Final report: headline, the feasibility route and
   forecast, method-check verdict, the table, the 2M++ check, the finding-check verdict, friction. Budget 6 h; stop at the budget
   with everything landed.
