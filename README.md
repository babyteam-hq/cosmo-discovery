# cosmo-discovery

An AI-agent research programme in observational cosmology, run by HermesAI since 31 October 2026: Claude agents choose
the questions, check the literature, write and run the analyses on public data, and write up the results, with a founder
steering by chat. Everything that went into the published papers is here, the scripts, the numbers, the referee verdicts
and the work log, and so is everything that did not make it.

The methodology, in one page: <https://hermesai.dev/research-methodology.html>.

## The papers

Three working papers stand, not peer reviewed, each passed by the blind gate described below; three more were
published and then withdrawn (after the table).

| Paper | Date | Scripts | Numbers |
|---|---|---|---|
| [ACT DR6 CMB lensing detects eRASS1 clusters at 5σ; their mass rescaling does not yet test the eROSITA weak-lensing calibration](https://hermesai.dev/research-erass1-act-lensing.html) | 1 Nov 2026 | `erass1_act.py` (imports `erass1_cmblens.py`), `act_dr6_fetch.sh`, `act_dr6_simstream.sh` | [`results-R122.md`](results-R122.md) |
| [ACT DR6 CMB lensing over 365 eRASS1 superclusters, against the ΛCDM expectation from the Websky light-cone](https://hermesai.dev/research-erass1-superclusters.html) | 2 Nov 2026 | `sc_members.py`, `sc_lens.py`, `sc_null.py` (with `erass1_cmblens.py`, `erass1_act.py`, `act_dr6_fetch.sh`, `act_dr6_simstream.sh`) | [`results-R137.md`](results-R137.md), [`results-R140.md`](results-R140.md) |
| [On new DESI sightlines the Giant Arc cannot yet be told from chance: a real Arc recovered by the discovery recipe predicts δ = +0.62 (68 per cent: +0.25 to +1.05) in its hull, and DESI DR1 sees 8 against 8.86 expected, 1.3σ below a real Arc and 0.5σ below chance](https://hermesai.dev/research-giant-arc.html) | 2 Nov 2026, revised 3 Nov 2026 | first version: `r155_template.py`, `r155_sightlines.py`, `r155_core.py`, `r155_run.py`, `r155_measure.py`, `r155_common.py`, [`r155_plan.md`](r155_plan.md); the revision, after a re-check found the first version compared against a selection-inflated amplitude: `garev_inject.py`, `garev_infer.py`, [`ga_rev_plan.md`](ga_rev_plan.md), [`finding-ga-rev.md`](finding-ga-rev.md) | [`results-R155.md`](results-R155.md) |

Withdrawn after publication, when the blind referee found them not sound or not novel:
[the eRASS1 cluster anisotropy note](https://hermesai.dev/research-erass1-anisotropy.html) (`erass1_lxt.py`),
[the NGC 1052 bullet-dwarf trail note](https://hermesai.dev/research-ngc1052-trail.html) (`trail_ngc1052.py`,
`trail_roman.py`, `dm_content.py`) and
[the Type Ia supernovae and Cosmicflows-4 note](https://hermesai.dev/research-cf4-supernovae.html), which passed on its
third pass, was found on re-check to have an uncalibrated null, and was withdrawn when its one revision, which built the
calibrated null, found the reference field 2M++ failing it too (first version: `cf4_vs_2mpp.py`, `cf4_decompose.py`,
`cf4_variants.py`, `cf4_2023.py`, `cf4_robust.py`, `cf4_notin.py`, `cf4_des.py`, `ztf_tests.py`, `bulk_flow.py`,
`cf4_mocks.py`, `cf4_referee2.py`, with `cf4_io.py`; [`results.md`](results.md); the revision: `cf4rev_common.py`,
`cf4rev_chain.py`, `cf4rev_fit.py`, `cf4rev_realwf.py`, `cf4rev_table.py`, [`plan-cf4-rev.md`](plan-cf4-rev.md),
[`results-cf4-rev.md`](results-cf4-rev.md), [`finding-cf4-rev.md`](finding-cf4-rev.md)).
Their sections and verdicts are still in [`finding.md`](finding.md).

## How the work is done

- **The blind gate** ([`blind-review.md`](blind-review.md)). A fresh Claude Opus agent with no conversation context
  reads only the text, does its own literature search and answers in a fixed verdict line. It runs three times per
  test: the *trajectory check* before any download (is the question novel?), the *method check* on the pre-registered
  plan before any compute (will the plan support the claim?), and the *finding check* on the written result (is it
  novel and sound?). Only the verbatim line "the findings are confirmed novel and sound if the calculations hold" is a
  pass. Every verdict is recorded verbatim, next to the result it judged.
- **The open-problem map** ([`open-problems.md`](open-problems.md), with the harvests in [`open-problems/`](open-problems)).
  Candidates come from a fixed taxonomy of twelve areas, harvested from primary sources, quick-killed against the
  literature, then put to the trajectory check. The yield log there records each round and what was changed when a
  stage stalled.
- **Dispatch.** Each candidate that survives goes to one agent with a brief (in [`briefs/`](briefs)) and a feasibility
  gate in its first hour. The work log, [`LOG.md`](LOG.md), is the running record of what was dispatched, what came
  back and why.
- **Sweeps and deep dives.** The earlier rounds' ranked candidate lists are `sweep*.md` and `deep-*.md`; a row per test
  is in `results.md` and `results-*.md`.

**The misses are kept on purpose.** Most tests here missed: a "novel but not sound" verdict, a feasibility gate that
stopped the work in its first hour, a question the literature had already answered. Each miss stays in the record as
it was judged, unreworded, because a programme that only shows its passes cannot be checked.

## Reproducing a paper

1. The scripts were written to run from the directory above this tree, with the tree at `research/cosmo-discovery/`,
   and some open their files by that path. Keep that layout:

       mkdir cosmo && cd cosmo
       git clone https://github.com/babyteam-hq/cosmo-discovery research/cosmo-discovery

2. Python 3.11 or later with numpy, scipy, astropy, healpy and matplotlib covers most scripts; each script's
   docstring names anything else it needs (pixell, sncosmo, dustmaps) and its usage line.
3. The data are public and not in this repository: each script reads from `data/` (create it in the tree), and
   [`SCRIPTS.md`](SCRIPTS.md) and the docstrings say where each file comes from: Pantheon+SH0ES
   (github.com/PantheonPlusSH0ES/DataRelease), the Cosmicflows-4 grids (projets.ip2i.in2p3.fr/cosmicflows), DES-SN5YR,
   ZTF SN Ia DR2 (ztfcosmo.in2p3.fr), the eRASS1 catalogues (eROSITA-DE DR1), ACT DR6 lensing
   (portal.nersc.gov/project/act/dr6_lensing_v1, fetched by `act_dr6_fetch.sh`), the Planck 2018 lensing
   reconstruction (Planck Legacy Archive), the Websky light-cone (CITA), the Zhu & Ménard Mg II catalogues, the DESI DR1
   Mg II and quasar catalogues (data.desi.lbl.gov) and SDSS DR16Q.
4. Run the paper's scripts (table above) as their usage lines say; the numbers to compare with are in the results
   file named in the table.

Paths in the notes are as they were in the working repository: `research/cosmo-discovery/` is this tree's root, and
`LOG.md` was the work log. References to the agents' tooling are described in plain words where they are not public.

## Licence

Code (`*.py`, `*.sh`, `*.js`) is under the MIT licence ([`LICENSE`](LICENSE)); the notes, results and figures are under
CC BY 4.0 ([`LICENSE-text`](LICENSE-text)). To cite this work, see [`CITATION.cff`](CITATION.cff). Questions to
team@hermesai.dev.
