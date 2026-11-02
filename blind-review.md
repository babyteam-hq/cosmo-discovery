# Blind review: the novelty and soundness gate

Will, 2026-10-31: an agent asking itself whether its idea is novel is not evidence. The gate is a fresh Opus with no
conversation context, given only the text, doing its own literature search, and answering in a fixed form. The pass line is
the verbatim verdict "the findings are confirmed novel and sound if the calculations hold". Anything short of it is a miss.

Run it twice per test, not once at the end, and run the first one in the dispatching session before a worker is started on a
brief (R-124 and R-125, 2026-11-01: both killed as not novel in their first half hour; one referee call here would have saved
each). **The trajectory check is the one that saves time**: a finding check catches a
wrong result after a day's compute; the trajectory check catches a wrong direction before the first download (Will,
2026-10-31). No sweep row is run, and no follow-up from a finding's "To settle it" list is started, without its trajectory
verdict recorded in `sweep-2.md` (section "Blind trajectory verdicts") or the next sweep's file. A "not novel" verdict kills
the row; "partly novel" narrows it to what the verdict says is open; a worker that computes without a verdict on file has
skipped the gate.

1. **Before computing** (the trajectory check): give the referee the question, the data, the method and the detectability
   estimate, as a sweep row reads. Ask for novelty only. A "not novel" verdict kills the row before any download.
1b. **After pre-registration, before computing** (the method check; added 2026-11-01 after five 'novel but not sound' misses in a row, each on a physics point knowable before compute: an assumed magnification coefficient, an uncalibrated contamination amplitude, an unvaried mock split, an unseparated systematic, a statistic the data cannot carry). The worker puts its fixed plan (samples, model, statistic, nulls, error budget, the parameters the result will depend on) alone in a file and runs the prompt below in a fresh context-free Opus agent. It fixes every point the referee names before computing, or stops and reports if one cannot be fixed. Recorded verbatim in the results file.
2. **After writing** (the finding check): give the referee the finding section as it stands, alone in a file. Ask for novelty
   and soundness. Record the verdict line verbatim in `results.md` next to the row, and in `finding.md` under the section.

How: the Agent tool, `subagent_type: general-purpose`, `model: opus`, one agent per paper, the paper in its own file so the
referee sees nothing else. The prompt below, unchanged except for the file path and the "closest prior work" hints (which
must name the obvious comparisons, not steer toward the answer). Never tell the referee who wrote it, what earlier referees
said, or that a null would be disappointing.

## Prompt (finding check)

    You are an independent referee for a short research note in observational cosmology. Read the note at <path> (use
    the Read tool). You have no other context about it and should not assume anything about who wrote it or why.

    Your job, in order:
    1. NOVELTY. Do your own literature search (WebSearch, extended mode; fetch arXiv abstracts with WebFetch where
       needed). Do not trust the note's own statements about what is or is not in the literature; verify them. Find the
       closest prior work: <the obvious comparisons for this note, e.g. "anyone who compared X with Y, or tested Z">.
       Cite what you find with arXiv ids.
    2. SOUNDNESS. Assume every number in the note is computed as described. Is the method sound and the conclusion
       supported? Look for <the failure modes that apply: circularity, mis-specified errors, selection effects, the trials
       problem, convention pitfalls, whether the robustness checks are the right ones>. Name the single strongest
       objection a referee would raise.
    3. VERDICT. End with exactly one of these lines, verbatim, and nothing after it:
       "VERDICT: the findings are confirmed novel and sound if the calculations hold"
       "VERDICT: novel but not sound: <one sentence why>"
       "VERDICT: sound but not novel: <closest prior work, arXiv id>"
       "VERDICT: neither: <one sentence why>"
    Be rigorous and adversarial. Give the full reasoning (up to ~600 words) before the verdict line, with the prior works
    you found listed with their arXiv ids.

## Prompt (method check, after pre-registration)

    You are an independent referee. Read the analysis plan at <path>: a measurement that has not been made yet, with its samples,
    model, statistic, null tests and error budget fixed in advance. You have no other context. Assume the data are as described.
    Your only question is soundness: when the numbers come in, will the stated claim be supported by this plan? Name every
    assumption the result will depend on that the plan does not measure or vary (a coefficient taken from elsewhere, an amplitude
    that is assumed, a mock choice never varied, a systematic not separated, a selection not modelled, a statistic the data cannot
    constrain), and the null or calibration that would close each. End with exactly one line, verbatim, and nothing after it:
       "METHOD: sound as planned"
       "METHOD: fixable: <the points to close before computing, semicolon-separated>"
       "METHOD: not measurable as planned: <one sentence why>"

## Prompt (trajectory check, before computing)

    You are an independent referee. Read the proposed test at <path>: a question, the public data, the method and the
    expected precision. You have no other context. Do your own literature search (WebSearch, extended mode; arXiv
    abstracts via WebFetch) for the closest prior work: <the obvious comparisons>. Do not trust the proposal's own
    literature statements. Then end with exactly one line, verbatim, and nothing after it:
       "VERDICT: novel: no prior work answers this question"
       "VERDICT: not novel: <closest prior work, arXiv id, and what it already answers>"
       "VERDICT: partly novel: <what is already answered, arXiv id; what remains open>"
    Give the reasoning and the prior works with arXiv ids before the verdict line.

## Record

Each verdict goes in `results.md`'s row (a "blind review" column entry or a sentence in "what it means") and, for findings,
in `finding.md` under the section as "Blind review (date, model): <verdict line>", followed by the referee's strongest
objection in one sentence. A miss is recorded as a miss; the section is not reworded to pass.
