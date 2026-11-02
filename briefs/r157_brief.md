# Site paper: the Giant Arc on independent sightlines (R-155, passed the blind gate) (R-157)

Goal: publish the fourth finding that passed the blind gate as a standalone paper on the site, previewed on staging.
Read first: the agents' working rules, the design-review rules, the R-155 section of
`research/cosmo-discovery/finding.md` (its gate verdict: "VERDICT: the findings are confirmed novel and sound if the calculations
hold"), `results-R155.md` (all rows, including the method-check verdict F4 and the disclosed rough count at proposal time),
`### cosmo-deep-23` in `LOG.md`, and the last paper as the template: `site/research-erass1-superclusters.html`,
built by `research/cosmo-discovery/briefs/mkpaper_sc.py` from `paper_sc_body.html` (R-142, R-144: wiring, sitemap, JSON-LD, og tags,
references with authors, year, title and arXiv id). Follow the agents' working rules: worktree on staging; land with `git pull --rebase origin staging
&& git push origin HEAD:staging`; run the CI jobs the CI change scoper scopes in (`site`).

## What to write
`site/research-giant-arc.html`, built by a `mkpaper_ga.py` from a `paper_ga_body.html` in the same pattern: a standalone paper. Title,
authors (as the earlier papers credit them), abstract, numbered sections, references; no Hermes header or site styling beyond what the
earlier research pages use. A figure where it reads better: the observed and expected absorber counts per cell and template, or the
amplitude A against the published 1.3 and the chance baseline, as inline SVG from the landed numbers.

Content rules (the gate passed this wording; do not strengthen it):
- The claim is the interval: A = -0.23 (68 % -0.71 to +0.26; 95 % -0.90 to +0.63) on 3,822 DESI DR1 sightlines never searched by
  SDSS; the published delta = 1.3 excluded at 3.2 sigma (3.7-3.8 with the other baselines); a weaker real structure (A <~ 0.6) is not
  excluded.
- State the referee's strongest objection (the published 1.3 is selection-inflated: the test falsifies the published amplitude but
  cannot split a modest real structure from chance without mock-based shrinkage).
- State the method check and its six closed points, and that a rough count was seen at proposal time with the plan unchanged after it.
- Data and methods: Zhu & Ménard SDSS catalogues, the frozen template, DESI DR1 Mg II VAC and sightlines, SDSS DR16Q exclusion, the
  nulls. Prior work with full references: Lopez et al. 2022 (2201.06875), the Big Ring (2402.07591), the Giant Ring (2604.17534),
  FLAMINGO-10K (2504.14940), Sawala et al. 2025 (2502.03515), Sawala & Teeriaho 2025 (2505.11072); look up authors and titles on the
  arXiv API (https://export.arxiv.org/api/query?id_list=...).
- Reproducibility: name the scripts in `research/cosmo-discovery/`.

## Then
1. Wire it as R-142 wired the supercluster paper; the site tests pass.
2. Land, deploy staging (the staging deploy command), check the page on staging.hermesai.dev.
3. Final report: the staging URL, the commit, anything the founders should look at before a Promote. Do not promote production.
   Budget 2 h.
