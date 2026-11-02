# cosmo-sweep-6: candidates from the literature's own open follow-ups (R-120)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read `LOG.md` (handoff at the top), `research/cosmo-discovery/README.md`, `blind-review.md`,
`sweep-3.md` and the R-119 section of `results.md` first. Follow the agents' working rules: worktree on staging, land after every row.

## Why this sweep is different
Three sweeps (R-115/116/118/119, ~90 tests) brainstormed anisotropy, dust and flow tests on public catalogues; the literature
check killed most, and the rest were nulls or catalogue systematics that the gate judged "novel but not sound". The candidate
source is the problem. This sweep takes candidates only from papers' own words: astro-ph.CO / GA / HE papers from 2026-04 to
2026-10 (arXiv listings, ADS) whose limitations, discussion or future-work sections name a test they did not do, where public
data this network can reach already suffice (DESI DR1 LSS catalogues, Pantheon+, DES-SN5YR, ZTF DR2, Quaia, CatWISE2020, Gaia
DR3, CHIME Cat 2, eRASS1, Cosmicflows-4, 2MRS, Planck PR3, SDSS/BOSS/eBOSS; zenodo.org, LAMBDA data dirs and GitHub release
assets answer 403 here, so skip anything needing them). Parked threads, do not pick from: SN Hubble residuals vs dust maps;
Cosmicflows-4 zero points and bulk flows; NGC 1052 trail; eRASS1 anisotropy; Pantheon+ sky anisotropy; quasar and galaxy
number-count dipoles. A candidate is a stated gap plus the public data that closes it plus a detectability estimate.

## Method
1. Harvest: at least 40 papers read at the level of their discussion and conclusions (the full text via arxiv.org/pdf or
   HTML, not the abstract), each yielding zero or more "they say it is not done, and it could be done with X" lines. Record
   every line in `sweep-4.md` with the paper's arXiv id and the quoted sentence (under 15 words of quote), then rank by
   detectability x cheapness.
2. Trajectory check (blind-review.md, the before-computing prompt) on each of the top candidates before any computation; a
   candidate whose gap another paper has since closed is dropped with the id.
3. Compute the survivors, best first, one row each in `results-R120.md` (your own file; the main session merges it; this
   avoids the rebase conflicts R-119 reported). Every number from a script in `research/cosmo-discovery/`, usage line in its
   docstring, data under `data/` (gitignored, never run as code, `python -I` for anything that reads a download).
4. A result that is a measurement on data with p < 1e-3 after a trials count and a spatially correlated null where the template
   is smooth gets a `finding.md` section and the finding check (blind-review.md's after-writing prompt, verbatim, one fresh
   context-free Opus agent, the section alone in a file). Record the verdict line verbatim. A miss is a miss: one revision
   pass on the referee's strongest objection if it is a presentation point, none if it is a physics point; then park it and go
   to the next candidate. Nothing goes on the site.
5. Keep-alive: a background loop running the session keep-alive every 20 min for the whole job
   (the worker key's one-hour gauge kills idle sessions). Land (`git pull --rebase origin staging && git push origin HEAD:staging`)
   after every row; a conflict in a shared file keeps both sides.
6. Update `LOG.md` with a `### cosmo-sweep-6` section at the end: rows, verdicts, what to do next.

## Report
Final report: headline (did anything pass the gate), rows added with one line each, verdicts verbatim, the harvest's size
(papers read, candidates found, candidates checked), friction. Budget: 8 h. Stop at the budget with everything landed.
