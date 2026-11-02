# cosmo-deep-25: is the spiral spin-direction asymmetry at 0.3 < z < 1.2 real? A bias-controlled test on space-based images (R-158)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff), `research/cosmo-discovery/blind-review.md` (the gate, including the method check
1b), `finding.md` (passed sections; and the withdrawn eRASS1 anisotropy and the sky-dipole misses, the lesson of this neighbourhood:
claims about sky directions fail on selection and look-elsewhere, so pre-register the axes and amplitudes before looking), and
`open-problems/round5/e2-spin-asym.md` (the proposal as the referee saw it). Follow the agents' working rules: worktree on staging; keep-alive (a
background loop running the session keep-alive every 20 min); land after every milestone. Rows in
`results-R158.md`. Long jobs with `nohup python -u` and an `until` watcher; foreground waits cap at 10 min, bare `sleep` is blocked,
match processes with `ps -eo pid,args`. `data/` starts empty: report fetch times. GPU is not available; size the classifier for CPU.

## Trajectory check: already run (2026-11-02, context-free Opus, 4 web searches and 6 API queries), verbatim:
"VERDICT: partly novel: removing classifier and annotator chirality with a reflection-equivariant network is already done at low z on
ground-based SDSS/DESI images (Jia, Zhu & Pen, arXiv:2210.04168), and Shamir's mirror-symmetric pipeline already runs on space data
(arXiv:2206.00168, 2502.18781); what remains open is an independent, bias-controlled, pre-registered test of the redshift-growing,
axis-aligned asymmetry on a large resolved space-based sample at 0.3<z<1.2, and no Euclid spin-direction analysis exists". It settles a
published disagreement (Shamir's space-based claims: HST COSMOS 2011.03723, 2206.00168; JWST 2403.17271, 2502.18781; vs the bias
explanations 0803.3247, 1610.07060, 2210.04168 and the isotropy tests 2404.06617, 2410.18884) by removing the disputed step:
classifier, selection and instrument chirality. Read all of these first. Do not rerun the trajectory check.

## Data, now (no need to wait for Euclid DR1)
HST COSMOS ACS F814W mosaics and the COSMOS2020 catalogue (photometric redshifts), the same images behind Shamir's HST claim; Euclid
Q1 VIS images and catalogues (public since March 2025, ~63 deg^2); JWST JADES / CEERS / COSMOS-Web public imaging where Shamir's JWST
claims sit. Euclid DR1-Foundation (planned November 2026) extends the sample when it appears: note it, do not wait for it.

## Step 0: feasibility gate (first 90 minutes; stop and report if it fails)
Count resolved spirals at 0.3 < z < 1.2 per dataset (size > a stated multiple of the PSF FWHM, axis ratio cut), and forecast the error
on the asymmetry fraction and on a dipole amplitude along Shamir's published axes (fixed now from his papers). Anchor: Shamir's own
HST and JWST amplitudes. If the combined samples cannot detect his claimed amplitude at 3 sigma (or exclude it at 3 sigma), stop and
report.

## Step 1: method check (blind-review.md 1b), before classifying anything
The fixed plan (samples and cuts, the classifier: a reflection-equivariant network trained on a stated public labelled set, validated
by exact sign flip on mirrored images; the PSF-chirality tests: asymmetry against focal-plane position, against size/PSF, and on
simulated symmetric galaxies convolved with the real PSFs; the statistic: global asymmetry and the dipole along the pre-registered
axes, with the look-elsewhere for any free axis; the redshift bins; the nulls) alone in a file; the method-check prompt verbatim in a
fresh context-free Opus agent. Close every point it names before classifying, or stop. Verdict verbatim in `results-R158.md`.

## Method
1. Classify handedness with the equivariant network; Ganalyzer-style or SpArcFiRe as a second classifier; agreement reported.
2. The asymmetry fraction per redshift bin and dataset, and the dipole amplitude along Shamir's axes, with binomial and bootstrap
   intervals; the PSF-chirality and mirror nulls; the same pipeline on stars and on simulated symmetric galaxies.
3. The claim is the table: per dataset and redshift bin, the asymmetry and the dipole amplitude with intervals, against Shamir's
   published values and against zero. No sentence about cosmological parity beyond what the numbers say.
4. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus agent,
   section alone; hints: the works above). Verdict verbatim in `results-R158.md` and under the section. One revision pass on a
   presentation point; park on a physics point. Nothing on the site; a pass is the headline of your report.
5. End: `### cosmo-deep-25` in `LOG.md`. Final report: headline, feasibility numbers, method-check verdict, the
   table, nulls, what Euclid DR1 would add, friction. Budget 6 h.
