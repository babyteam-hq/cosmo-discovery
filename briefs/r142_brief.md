# Site paper: the mass of eRASS1 superclusters outside their member clusters (R-137 + R-140, passed the blind gate) (R-142)

Goal: publish the third finding that passed the blind gate as a standalone paper on the site, previewed on staging.
Read first: the agents' working rules, the design-review rules, the R-137 section of
`research/cosmo-discovery/finding.md` as revised by R-140 (its gate verdict: "VERDICT: the findings are confirmed novel and sound if
the calculations hold"), `results-R137.md`, `results-R140.md`, `### cosmo-deep-13` in `LOG.md`, and the last
paper's page and commit as the template: `site/research-erass1-act-lensing.html` and commit 90d15e04 (how it was wired: JSON-LD,
og tags, sitemap, SEO tests). Follow the agents' working rules: worktree on staging; land with `git pull --rebase origin staging && git push origin
HEAD:staging`; run the CI jobs the CI change scoper scopes in (`site` at least) before landing.

## What to write
`site/research-erass1-superclusters.html`: a standalone paper. Title, authors (as the earlier papers credit them), abstract,
numbered sections, references; no Hermes header or site styling beyond what the earlier research pages use. A figure where it reads
better than a table (the null comparison per sample and aperture is the natural one: observed vs the ΛCDM expectation with
intervals), drawn as inline SVG or a static image generated from the landed numbers.

Content rules (the gate passed this wording; do not strengthen it):
- The claim is the interval: E - E_control against the ΛCDM null per sample and aperture, the data within 1.1 sigma of ΛCDM in
  every row, the upper limits relative to ΛCDM (0.9 / 0.3 / 1.1 / 1.0 of the member model) as R-140 states them, and the 6.2 sigma
  detection of footprint lensing.
- State the referee's strongest objection (an environment-dependent count-rate mass bias not in the null, ~0.1, a quarter of
  sigma) and the hull-versus-disc 1.9 sigma difference as stated-not-used.
- Data and methods: eRASS1 superclusters (Liu et al. 2024, arXiv:2402.08454; members rebuilt because CDS was unreachable, 991/1338
  reproduced exactly: say so), ACT DR6 lensing with its 400 simulations, the Websky light-cone null (Stein et al. 2020), Planck PR3
  as the check. Prior work: arXiv:2203.11306, 2409.04535, 1805.04555, 2404.17345, 2310.04387.
- Reproducibility: point to `research/cosmo-discovery/` scripts by name.

## Then
1. Wire it as 90d15e04 wired the last one (sitemap, JSON-LD, og tags, the research index if there is one); the site tests pass.
2. Land, then deploy staging (the staging deploy command), and check the page on staging.hermesai.dev.
3. Final report: the staging URL, the commit, anything the founders should look at before a Promote. Do not promote production.
   Budget 2 h.
