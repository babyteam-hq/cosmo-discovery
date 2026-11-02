# cosmo-deep-7: time delays of cluster-scale lensed quasars from ZTF difference imaging, tested against pre-published lens-model predictions (R-127)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff at the top), `research/cosmo-discovery/blind-review.md` (the gate), `finding.md`
(what passed and what missed; the misses were claims stronger than their intervals). Follow the agents' working rules: worktree on staging;
keep-alive (a background loop running the session keep-alive every 20 min); land after every
milestone with `git pull --rebase origin staging && git push origin HEAD:staging`. Rows in `results-R127.md` (own file). Long
jobs with `nohup python -u` into a log file and an `until` watcher; logs to files.

## Trajectory check: already run (2026-11-01, context-free Opus, in the dispatching session), verbatim:
"VERDICT: partly novel: measuring wide-separation lensed-quasar time delays from public ZTF forced-photometry light curves is
already done (Nunez-Pizarro et al. 2026, A&A 706 L12, GraL J1651-0417, ~10" quad; arXiv id not found), and PSF-forward-model
deblending for delays is established on dedicated monitoring (TDCOSMO XVII, arXiv:2504.02932); what remains open is ZTF delays for
cluster/group-scale systems (J1004, J1029, J2222, J1326, COOL J0542 etc.), scene-model forced photometry on ZTF difference images,
and a sample-wide comparison with pre-published cluster-model predictions". The open part is the job; it stands alone because
no ZTF delay exists for any cluster-scale lensed quasar and several have model predictions published before any measurement
(J2222 images D, E, F arXiv:1609.08848; COOL J0542-2125 arXiv:2209.03972; COOL J0335-1927 arXiv:2305.14317; COOL J1153+0755
arXiv:2605.13514; J1326+4806 no published delay). Published delays to validate against: J1004+4112 (arXiv:2206.08597), J2222+2745
(Dahle et al. 2015, arXiv:1505.06187), J1029+2623 (Fohlmeister et al. 2013). Also arXiv:2608.07470 (ZTF delays vs HST models for
J1651). Do not rerun the trajectory check; cite these. Caveat the referee raised: the 2018-2025 baseline limits which long-delay
pairs are measurable; say so per pair before fitting.

## Data (public; verify reachability from the cloud first; zenodo.org, LAMBDA and GitHub release assets 403; IRSA is to be tried)
ZTF public data (IRSA): the forced-photometry service (per-position light curves), the difference-image cutout service
(image cutouts at a position and epoch), and the science-image cutouts; image positions from HST and Legacy Survey DR10 (the
discovery papers' tables); the public cluster lens models (RELICS, SGAS, the papers' lenstool / GLAFIC releases) or their published
delay predictions. If the difference-image cutouts are unreachable at scale, fall back to the forced-photometry service at each
image position and state the blending it implies; the scene model is the better product but not the gate.

## Method
1. Targets: the cluster- and group-scale lensed quasars with image separations >= 3" and r <= 20.5: SDSS J1004+4112, J1029+2623,
   J2222+2745, J1326+4806, COOL J0542-2125, COOL J0335-1927, COOL J1153+0755, and any other wide system in the Lemon et al. and Gaia
   GraL compilations (arXiv:2603.12810 lists 73 systems with delays; take the wide ones without ZTF delays). Fix the target list and
   the pairs to be fitted before any light curve is extracted; mark which pairs the baseline can reach.
2. Photometry: on each ZTF difference image, a fixed-position multi-PSF scene (N point sources at the HST positions, PSF from field
   stars, only fluxes free) in g and r; per-image light curves 2018-2025. Validate on J1004 and J1029 first: recover their published
   delays within errors or stop and report why.
3. Delays: PyCS3-style curve shifting with spline microlensing models and free-knot extrinsic variability; uncertainties from
   simulated light curves with the same cadence, gaps and noise (the standard TDCOSMO procedure); report each pair's delay with its
   interval, or the excluded delay range when no delay is recovered.
4. The test: each new delay against the model prediction published before the measurement, with the prediction's own interval; the
   sample-wide statement is the count of predictions confirmed and rejected, and by how much. H0 is not the claim; a delay set is.
5. The claim is the table of delays with intervals and the prediction comparison; nothing beyond it.
6. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus agent,
   section alone; hints: the prior works above). Verdict verbatim in `results-R127.md` and under the section. One revision pass
   on a presentation point; park on a physics point. Nothing on the site; a pass is the headline of your report.
7. End: `### cosmo-deep-7` in `LOG.md`. Final report: headline, the table (system, pair, delay, interval,
   predicted, confirmed?), the validation on J1004/J1029, friction. Budget 8 h; stop at the budget with everything landed.
