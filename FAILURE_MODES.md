# Failure modes of an AI research agent at a blind gate

`gate-failures.csv` lists every gate verdict in this repository that was not a pass: 132 of them, from 31 October to
3 November 2026. Claude agents chose, ran and wrote up tests on public cosmology data. A blind referee (a fresh model with
no context, given only the text, doing its own literature search, `blind-review.md`) checked them at up to four stages:

| stage | when | pass line | rows |
|---|---|---|---|
| trajectory | before any download: is the question already answered? | "VERDICT: novel: no prior work answers this question" | 74 |
| method | after the plan is fixed, before compute: will the plan carry the claim? | "METHOD: sound as planned" | 4 |
| finding | after the write-up: novel and sound? | "VERDICT: the findings are confirmed novel and sound if the calculations hold" | 44 |
| feasibility | the test's own step-0 forecast: can the data reach the line at all? | the forecast clears the brief's threshold | 10 |

The method check started on 2026-11-01, after five finding-stage misses in a row on physics points knowable before compute,
so it has few rows. Feasibility stops are the agent's own forecast, not the referee's; they are in because they stopped a
test, and the referee was never asked. So 122 rows are the blind referee's and 10 are the agent's own. The referee is a fresh
model of the same family as the agent, so their blind spots may overlap, and errors both miss are not in the file.

## The denominator: passes

`gate-passes.csv` lists the passes the same way (one row per verdict, same columns but `mode` and `after`):

| stage | passes | misses (`gate-failures.csv`) | misses / all verdicts |
|---|---|---|---|
| trajectory | 9 | 74 | 89 % |
| method | 0 | 4 | 100 % |
| finding | 4 | 44 | 92 % |
| feasibility | 7 | 10 | 59 % |
| all | 20 | 132 | 87 % |

Read the last column with care: it counts verdicts, not tests. A test can miss several times and then pass (the four
finding-stage passes are the four working papers, each passed after earlier misses), and a "partly novel" trajectory miss
usually narrowed a test that then ran. Feasibility passes count only the step-0 forecasts the files record as passing;
two forecasts that cleared their line without a recorded pass line (`plan-R154.md`, the ELG section of `finding.md`) are
not counted. The four method checks all asked for a fix; none passed as planned.

The record table on the methodology page (https://hermesai.dev/research-methodology.html) counts by Request and note: 48
finding-check verdicts (5 passes, the fifth a re-check of a page's reworded title and abstract) and 6 tests stopped at a step-0 gate. These files count every verdict line in the
record, so their numbers are larger: they add the forecast kills inside the earlier sweeps (GF-034 to GF-036 and GF-103)
and the first review of the NGC 1052 note (GF-032), among others.

## Columns

- `id`: GF-001 to GF-132, by date, then stage.
- `date`: the verdict's date as the file states it (blank where the file gives none).
- `stage`: trajectory, method, finding or feasibility.
- `test`: the Request (R-N) and the row or note, as the source names it.
- `verdict`: the verdict line verbatim, cut at about 300 characters with "...". For a table cell, the verdict cell and the
  prior-work cell, joined by " | ".
- `mode`: one label from the ten below, the referee's main reason.
- `after`: what the files say changed next (killed, narrowed, re-run, parked, withdrawn), in a few words.
- `source`, `line`: the file in this repository and the line holding the verdict, as of tag `v1.2`. `LOG.md` is the work log.
- `trace`: the same place as a link, the file at tag `v1.2` with the line highlighted
  (`https://github.com/babyteam-hq/cosmo-discovery/blob/v1.2/<source>?plain=1#L<line>`). Later versions of a file may move
  the line; the tag keeps it where the row says.

Each row was extracted from the committed files and checked against them by script: the quoted verdict is at the cited
line of the published tree (the export's test checks every row, in both files). Where one verdict is recorded in several files, it has one row, citing the primary record (a `results-*.md` or
`sweep*.md` row, a `finding.md` section, or an `open-problems/round*/verdicts.md` table) over the work log. Summary lines
that only restate a verdict recorded elsewhere have no row.

What it is not: a benchmark with ground truth, or a list of proven errors (a row is a referee's flag). A verdict is one model's judgement. Some misses were later passed on revision
(GF-088), and a "partly novel" trajectory verdict usually narrowed a test rather than killing it. The labels are ours,
assigned by an agent reading each verdict, one label per row, not checked by a second rater; treat them as a starting
taxonomy.

## The modes

Counts are rows (trajectory / method / finding / feasibility).

### partly-novel: 70 (70 / 0 / 0 / 0)
Part of the question was already answered in the literature; the test narrowed to what was open, or was dropped when the
open part was too small. The most common verdict by far: the agent's first framing almost always overlapped published work.
**Example, GF-018** ([sweep-2.md:83](sweep-2.md?plain=1#L83)): the agent proposed testing Migkas's cluster residuals
against the ROSAT exposure map; the referee found "residual vs RASS exposure and anisotropy vs exposure already tested,
null: Migkas+2020 2004.03305 Sec. 6.4.1-6.4.2". The agent's script had repeated a published null.

### not-novel: 6 (4 / 0 / 2 / 0)
The question, or the whole result, is already in the literature. Includes the finding-stage "sound but not novel" and
"neither". **Example, GF-070** ([results-R125.md:12](results-R125.md?plain=1#L12)): CHIME fast-radio-burst dispersion
against the Compton-y map, killed in its first half hour with nothing downloaded: Takahashi et al. (arXiv:2511.02155) had
already done it. **Withdrawn note, GF-032** ([results.md:109](results.md?plain=1#L109)): "VERDICT: neither" on the NGC 1052
trail note, whose main tests were the discovery papers' own and which misdescribed them; the published page became a
withdrawal notice the same day.

### unseparated-systematic: 18 (0 / 1 / 17 / 0)
The signal could be produced by a known systematic or selection effect the analysis did not separate or model. The most
common finding-stage miss. **Example, GF-026** ([finding.md:539](finding.md?plain=1#L539)): the eRASS1 cluster
anisotropy tracked the survey's exposure, but in a half-sky the exposure is nearly collinear with the anisotropy axis, so
the fit cannot tell them apart. Three passes later the claim had narrowed to "neither confirms nor excludes", and the
published note was withdrawn.

### undetectable: 10 (0 / 0 / 0 / 10)
The step-0 forecast says the data cannot reach the threshold, so the test stops before measuring. **Example, GF-127**
([results-R148.md:69](results-R148.md?plain=1#L69)): Lyman-alpha forest x Compton-y. The per-quasar continuum fit removes
the very mode the y map sees; the forecast was about 1e-4 sigma per bin, so nothing was measured.

### unvaried-mock: 7 (0 / 1 / 6 / 0)
A mock, simulation or null choice the result depends on was never varied or calibrated. **Example, GF-116**
([results-R138.md:696](results-R138.md?plain=1#L696)): the SPHEREx x DESI amplitudes "depend on the unvaried compact/smooth
split and smooth-light profile of the injected" mock; the detection stood, the corrected amplitudes did not.

### assumed-coefficient: 5 (0 / 2 / 3 / 0)
A coefficient or amplitude taken from elsewhere instead of measured on the data. **Example, GF-117**
([results-R146.md:130](results-R146.md?plain=1#L130)): the z~5 lensing amplitude was a ratio "dominated by an assumed
magnification coefficient that was never measured". The method check (GF-112, [results-R156.md:75](results-R156.md?plain=1#L75))
was added to catch this before compute, and it did: "calibrate the Tractor size response ... instead of assuming 1".

### misstated-errors: 5 (0 / 0 / 5 / 0)
Error bars or covariance mis-specified, so the significance is wrong. **Example, GF-102**
([results-R139.md:105](results-R139.md?plain=1#L105)): a quasar bias of 6.2 +- 1.9 averaged two sky regions that disagree
at 2.4-2.8 sigma; with the error scaled to cover them it is 6.2 +- 3.6, and the preferred model is no longer decisive.

### overclaimed: 5 (0 / 0 / 5 / 0)
The conclusion or interval says more than the numbers carry. **Example, GF-030** ([results.md:104](results.md?plain=1#L104)):
"a regression slope beta < 1 cannot tell an over-predicted amplitude from structure uncorrelated with the truth"; the note's
"over-predicted three to ten times" was withdrawn and the paper restated, after which it passed.

### look-elsewhere: 3 (0 / 0 / 3 / 0)
A trials factor ignored, or a sample, cut or statistic chosen after seeing the data. **Example, GF-088**
([results-R122.md:29](results-R122.md?plain=1#L29)): the eRASS1 x ACT lensing intervals came "from a sample chosen after
unblinding". The note went back to the samples fixed before the first fit and passed on its second check.

### circular: 3 (0 / 0 / 3 / 0)
The method assumes what it tests. **Example, GF-119** ([results-R147.md:399](results-R147.md?plain=1#L399)): the
leakage calibration used a mock that "already builds in the early/late halo" difference under test. **GF-092**
([results-R126.md:175](results-R126.md?plain=1#L175)): an "injection-free" count of gravitational-wave events equals the
known list by construction.

## What changed after

Of the 44 finding-stage misses, 27 rows end in the note being parked, most after one revision; the rest record a
revision (a null run, a claim narrowed or reworded) before the next pass. Two notes went on to pass after revision: the
Cosmicflows-4 supernova note (GF-021, GF-030) and the eRASS1 x ACT lensing note (GF-088). Three published notes were
withdrawn: NGC 1052 (GF-032), the eRASS1 anisotropy note (GF-026, GF-029, GF-031) and the Cosmicflows-4 supernova note,
whose one revision built the calibrated null the re-check asked for and found the 2M++ reference field failing it too
(GF-132). The `after` column has each one.
