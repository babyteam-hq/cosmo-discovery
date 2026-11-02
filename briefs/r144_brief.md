# Site paper fix: full references in the supercluster paper (R-144)

Goal: finish the references of `site/research-erass1-superclusters.html` (R-142, the paper for the finding that passed the blind
gate) before a Promote. Follow the agents' working rules: worktree on staging; land with `git pull --rebase origin staging && git push origin
HEAD:staging`; run the CI jobs the CI change scoper scopes in (`site`); then the staging deploy command and check the page.

The page is built by `research/cosmo-discovery/briefs/mkpaper_sc.py` from `research/cosmo-discovery/briefs/paper_sc_body.html`:
edit the body and rebuild, so the two stay in step.

1. Replace three reference items (looked up on the arXiv API):
   - `arXiv:2203.11306. CMB lensing by DES Y3 superclusters, stacked on the Planck lensing map (8.4&sigma;).` →
     `Kovács, A., et al. 2022, Dark Energy Survey Year 3 results: imprints of cosmic voids and superclusters in the Planck CMB lensing map, arXiv:2203.11306.`
   - `arXiv:2310.04387. The two-halo term in cluster-lensing stacks.` →
     `Toscano, F., et al. 2023, A CMB lensing analysis of the extended mass distribution of clusters, arXiv:2310.04387.`
   - `arXiv:2404.17345. Weak lensing of eRASS1 supercluster members against isolated clusters.` →
     `Liu, A., et al. 2024, The SRG/eROSITA All-Sky Survey: exploring halo assembly bias with X-ray selected superclusters, arXiv:2404.17345.`
2. Find arXiv IDs for Clerc et al. 2024 and Kluge et al. 2024 (the eRASS1 papers the page cites without one) on the arXiv API
   (https://export.arxiv.org/api/query?search_query=...) and add them; if one has none, leave it as is and say so.
3. Check every other reference has authors, year, title and ID in the same form; fix any that don't.
4. Final report: the staging URL, the commit, the CI result. Do not promote production. Budget 45 min.
