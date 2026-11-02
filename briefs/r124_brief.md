# cosmo-deep-4: NANOGrav 15-year continuous-wave limits on the periodic-quasar binary candidates (R-124)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff at the top), `research/cosmo-discovery/blind-review.md` (the gate; trajectory
check before any compute, finding check after writing), `finding.md` (what passed and what missed; the misses were claims stronger
than their intervals). Follow the agents' working rules: worktree on staging; keep-alive (a background loop running
the session keep-alive every 20 min); land after every milestone with
`git pull --rebase origin staging && git push origin HEAD:staging`. Rows in `results-R124.md` (own file). Long jobs with
`nohup python -u` into a log file and an `until` watcher; logs to files.

## The question
Time-domain surveys have published hundreds of periodic-quasar candidates claimed as supermassive black-hole binaries (Graham et
al. 2015 CRTS, Charisi et al. 2016 PTF, Chen et al. 2020/2024 ZTF and DES, Liu et al. PS1; later compilations). Each candidate's
claimed period, redshift and black-hole mass imply a continuous gravitational-wave strain at a known sky position and frequency.
Targeted continuous-wave searches on pulsar-timing data exist for a handful of famous sources (3C 66B, OJ 287, PKS 2131-021) and
all-sky limits exist, but a systematic pass over the whole candidate catalogue, framed as the population-level fraction of
candidates excluded by the NANOGrav 15-year data (and what chirp mass and distance each surviving one is limited to), appears not
to have been published. The result is an interval either way: "x of N candidates are excluded as binaries at their claimed
parameters (95 %)" or "none are, and the limits sit above the predictions by a factor y".

## Data (public; verify reachability from the cloud first; zenodo.org and GitHub release assets 403; source trees and NANOGrav's
own data server are to be tried)
NANOGrav 15-year data set (timing files, noise models; the `enterprise`, `enterprise_extensions` and PINT packages); the 15-year
continuous-wave paper (arXiv:2306.16222) for the method and its all-sky limits; the periodic-quasar candidate catalogues from
VizieR or the papers' tables; black-hole masses from the catalogues or standard virial estimates where given.

## Method
1. Trajectory check first (blind-review.md's before-computing prompt, verbatim, fresh context-free Opus agent): hints: anyone who
   ran targeted PTA continuous-wave searches on the periodic-quasar candidate catalogues as a population, or compared the
   candidates' implied strains with published PTA upper limits (e.g. Sesana et al. 2018, Xin et al., the NANOGrav 15yr CW paper's
   discussion). If the population comparison against the all-sky limits is already published, the open part is the targeted
   (fixed sky position and frequency) analysis per candidate; decide whether that stands alone before computing.
2. For each candidate: implied strain h0 from its period, redshift and mass (with the mass uncertainty and the unknown mass ratio
   and inclination treated as nuisances, stated), sky position and GW frequency (twice the orbital frequency; state the
   period-doubling assumption).
3. Cheap pass: compare each h0 against the 15-year all-sky or sky-dependent 95 % upper limit at that frequency (the published
   sensitivity map if released; else your own from the noise models). Report the fraction excluded with the mass-ratio and
   inclination nuisances marginalised.
4. Targeted pass on the top candidates (the highest implied strain, say 20 to 50): fixed-position, fixed-frequency Bayesian upper
   limits with `enterprise` on the 15-year data, red noise and the common process included as the 15-year CW paper does; the
   Bayes factor for signal vs noise; a null from the same analysis at scrambled sky positions.
5. The claim is the table of candidates with their limits and exclusion status, plus the population fraction with its error; no
   statement about what the candidates are beyond "excluded as a binary at the claimed parameters" or "not excluded".
6. Write the result as a `finding.md` section (ninth), then the finding check (after-writing prompt verbatim, fresh context-free
   Opus agent, section alone; hints as in step 1). Verdict verbatim in `results-R124.md` and under the section. One revision pass
   on a presentation point; park on a physics point. Nothing on the site; a pass is the headline of your report.
7. End: `### cosmo-deep-4` in `LOG.md`. Final report: headline, the table (candidate, period, z, mass, h0
   implied, h0 limit, excluded?), the fraction, the null, friction. Budget 6 h; stop at the budget with everything landed.
