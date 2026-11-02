# R-123 (cosmo-deep-3): which dataset pays for the DESI DR2 w0wa preference

Cloud worker R-123, 2026-11-01. Script `desi_shapley.py` (run with `data/r123/venv/bin/python`); data in `data/r123/`
(gitignored): CobayaSampler/bao_data and sn_data clones, the Planck 2018 baseline likelihood tarball from the PLA, the
cobaya `packages/` tree. Outputs in `data/r123/out/` (JSON per run, logs).

## Blind trajectory verdict (before any compute; `blind-review.md` prompt verbatim, context-free Opus)

VERDICT: partly novel: dataset dependence of the DESI w0wa preference (SN sample dominant, CMB likelihood minor, DR1 vs
DR2, priors) and drop-one tests are already answered piecemeal (DESI DR2 arXiv:2503.14738; arXiv:2504.18464;
arXiv:2411.11685; arXiv:2509.13318; arXiv:2603.05472); an exact Shapley attribution of Δχ² over dataset coalitions and a
pre-registered specification curve with variance decomposition by analysis choice remain open

Decision (entry rule: "partly novel" is a no unless the open part is a result on its own): go. The open part is the whole
proposed deliverable, not a corner of it. Prior work the referee named, to cite and not repeat: DESI DR2's own drop-one
table (2503.14738, Table 6), Notari, Redi & Tesi (2411.11685: ~3σ without DESI at all, the low-z SNe shared by DES5Y and
Pantheon+), Efstathiou (2408.07175: per-sample Δχ² of DES5Y vs Pantheon+, a DES5Y low-z magnitude offset), the
prior-dependence paper (2509.13318: theory priors take it to 1.3-1.8σ), the evidence grid of Ong, Yallup & Handley
(2603.05472: DR1/DR2 × three Planck likelihoods × four SN samples incl. DES-Dovekie, under which the preference goes
away), 2504.18464 (ACT DR6 swapped in: 2.5-4.2σ set by the SN sample), 2607.28918 (an e-value reanalysis calling the DR2
result specification-dependent). The referee's own expectation: the SN sample, DES-SN5YR's low-z half in particular,
will carry most of the weight. Risk recorded before computing: if that is what comes out, the result is the exact
attribution with errors, not a surprise, and it is written as the table and nothing more.

## Pipeline and the reproduction check (gate: published Δχ² within 0.3σ)

Background-only w0waCDM (photons, CAMB's neutrino convention, massive neutrino by the exact Fermi-Dirac integral,
distances by quadrature; r_d, r_* and z_* from CAMB tabulated on an (ω_b, ω_c, Σm_ν) grid). Validated against CAMB:
D_M to < 1e-5 at z ≤ z_* (3e-5 at Σm_ν = 0.3 eV), θ_* to 3e-6 (the compressed prior's σ(θ_*)/θ_* is 8e-6), r_d and r_* to
1e-5, Ω_m and ω_ν to 1e-4. 0.8 ms per evaluation. Likelihoods as cobaya ships them: DESI DR2 and DR1 BAO Gaussians;
Pantheon+ (z_HD > 0.01, 1590 SNe, STAT+SYS), Union3 (22 bins), DES-SN5YR (1829, covsys + diagonal), DES-Dovekie (1820,
inverse total covariance), each with the absolute magnitude marginalised analytically (as cobaya and DESI do); the
compressed CMB is DESI DR2's (θ_*, ω_b, ω_bc) Gaussian (2503.14738 App. A). Priors as DESI's: w0 ∈ [−3, 1], wa ∈ [−3, 2],
w0 + wa < 0. Minimiser: Nelder-Mead then Powell from 6-8 starts (the DESI-like point, ΛCDM, random in the box), best kept.
Significance from Δχ² with 2 degrees of freedom, as DESI does.

Rows of DESI DR2 Table 6 that a compressed pipeline can reproduce (Δχ²_MAP, significance):

| combination | this pipeline | DESI DR2 | best fit here (w0, wa, Ω_m, h) |
|---|---|---|---|
| DESI | −4.7 (1.7σ) | −4.7 (1.7σ) | −0.18, −2.71, 0.386 |
| DESI + (θ_*, ω_b, ω_bc)_CMB | −8.1 (2.4σ) | −8.0 (2.4σ) | −0.443, −1.65, 0.351, 0.637 (DESI: −0.43 ± 0.22, −1.72 ± 0.64, 0.353, 63.7) |
| DESI + Pantheon+ | −4.9 (1.7σ) | −4.9 (1.7σ) | −0.89, −0.21, 0.305 |
| DESI + Union3 | −10.0 (2.7σ) | −10.1 (2.7σ) | −0.70, −1.00, 0.330 |
| DESI + DESY5 | −13.5 (3.3σ) | −13.6 (3.3σ) | −0.78, −0.72, 0.321 |

Every published row reproduces to 0.1 in Δχ² (0.0-0.1σ). With the compressed CMB in place of the full one the headline
rows come out at −7.6 (2.3σ), −13.8 (3.3σ), −18.0 (3.8σ) for Pantheon+, Union3, DESY5 against DESI's full-CMB −10.7
(2.8σ), −17.4 (3.8σ), −21.0 (4.2σ): the full CMB (ISW, lensing smoothing, CMB lensing) adds 3-4 units of Δχ² that the
(θ_*, ω_b, ω_bc) compression cannot carry, as DESI's own −8.0 / −9.7 / −12.5 sequence (compressed / full without lensing /
full) shows. The compressed grid therefore understates the CMB player; the full-likelihood check on the headline
coalitions (below) measures by how much. DESI + CMB_c + DESY5 with z > 0.1: −3.9 (1.5σ).

## Frozen specification list (written before the curve was computed)

Factors and levels; every combination is run unless marked as collapsing:

1. SN sample: none; Pantheon+; Union3; DES-SN5YR; DES-Dovekie.
2. SN low-z cut: as released; z > 0.1 (collapses when SN = none).
3. BAO release: DESI DR2; DESI DR1.
4. CMB: none; DESI's (θ_*, ω_b, ω_bc) prior; Planck 2018 (R, l_A, ω_b) distance priors (Chen, Huang & Wang 2019, wCDM row).
5. Dark-energy model: CPL (w0, wa), 2 dof; wCDM (w0), 1 dof.
6. Prior: DESI box with w0 + wa < 0; DESI box without the cut; narrow box w0 ∈ [−2, 0], wa ∈ [−2, 1] (collapse for wCDM:
   the cut is moot, the narrow box differs).
7. Neutrino mass: Σm_ν = 0.06 eV fixed; Σm_ν free in [0, 0.5] eV (profiled).
8. Metric: profile Δχ² to a Gaussian-equivalent significance with the model's dof; the parameter-shift distance of
   (w0, wa) = (−1, 0) from the best fit in the Hessian metric at the best fit, same dof.

Not in the grid (cost or reach): the full CMB likelihoods (the headline-coalition check instead), a thawing or binned
w(z), Bayesian evidence. The variance decomposition is a main-effects ANOVA of the significance over the grid (sum of
squares per factor over the total), with the factor-level means listed.

## Rows

### Row 1: the 5-player Shapley game on the compressed likelihoods (`shapley5.json`)

Players: DESI DR2 BAO, CMB_c (DESI's compressed prior), Pantheon+, Union3, DES-SN5YR. 32 coalitions, each with its
ΛCDM and w0waCDM profile minima (4 starts); v(S) = χ²_min(ΛCDM) − χ²_min(w0wa) ≥ 0. Grand coalition (all five, the three
SN samples stacked as if independent, which they are not: they share low-z SNe) Δχ² = 23.9 (4.5σ). Exact Shapley values:

| player | φ (Δχ² units) | share |
|---|---|---|
| DESI DR2 BAO | 9.64 | 40 % |
| CMB_c | 3.08 | 13 % |
| Pantheon+ | −0.63 | −3 % |
| Union3 | 4.18 | 18 % |
| DES-SN5YR | 7.59 | 32 % |

The 3-player games that do not stack SN samples (DESI, CMB_c, one SN sample), read off the same coalition values:

| game | grand Δχ² | DESI | CMB_c | SN |
|---|---|---|---|---|
| with Pantheon+ | 7.6 (2.3σ) | 5.83 (77 %) | 1.55 (20 %) | 0.22 (3 %) |
| with Union3 | 13.8 (3.3σ) | 7.03 (51 %) | 2.05 (15 %) | 4.75 (34 %) |
| with DES-SN5YR | 18.0 (3.8σ) | 7.99 (44 %) | 2.17 (12 %) | 7.85 (44 %) |

Per-dataset χ² gain at the joint best fits (χ²_k at the ΛCDM minimum minus χ²_k at the w0wa minimum; the "sacrifice"
view, which Efstathiou 2408.07175 and DESI give in part): DESI+CMB_c+Pantheon+: DESI +2.7, CMB_c +2.2, SN +2.7;
+Union3: +4.2, +2.8, +6.8; +DES-SN5YR: +4.6, +2.0, +11.4; all five: DESI +6.2, CMB_c +1.0, Pantheon+ +0.5, Union3 +5.5,
DES-SN5YR +10.5. Single-player values: DESI alone 4.7, DES-SN5YR alone 5.9 (at the wa = −3 prior edge), Union3 alone 3.4
(same edge), Pantheon+ alone 0.5, CMB_c alone 0. Pairs without DESI: CMB_c+DES-SN5YR 6.6 (2.1σ), CMB_c+Union3 4.8,
CMB_c+Pantheon+ 1.0. Errors on the shares: the parametric bootstrap, row 2.

### Row 3: the specification curve (`speccurve.json`, `desi_speccurve.png`; run before row 2 finished)

540 specifications from the frozen list, all run (0 failures), compressed likelihoods throughout, 3 starts each.
Significance from the profile Δχ² with the model's dof (2 for CPL, 1 for wCDM). Range 0.01-3.84σ, mean 1.59, sd 0.94.
Main-effects ANOVA of the significance over the 480 specifications that contain an SN sample (where the low-z-cut factor
is defined; the 60 SN-free ones behave the same without that factor):

| factor | share of variance | level means (σ) |
|---|---|---|
| SN low-z cut | 44.5 % | as released 2.27; z > 0.1 0.97 |
| DE model (dof) | 15.8 % | CPL 1.93; wCDM 1.14 |
| SN sample | 6.4 % | Pantheon+ 1.39; Union3 1.86; DES-SN5YR 1.87; DES-Dovekie 1.35 |
| prior box | 2.6 % | DESI 1.54; DESI without w0+wa<0 1.93; narrow 1.54 |
| BAO release | 1.6 % | DR2 1.74; DR1 1.49 |
| CMB form | 0.0 % | none 1.63; (θ_*, ω_b, ω_bc) 1.63; (R, l_A, ω_b) 1.59 |
| neutrino mass | 0.0 % | fixed 1.63; free 1.60 |

Main effects sum to 71 %; the rest is interaction (the cut matters most where DES-SN5YR is the sample; the prior only
where the data are weak). The headline slice (DESI DR2 + CMB_c + SN, CPL, DESI priors, Σm_ν fixed) and what each
single change does to it:

| SN sample | base | SN z > 0.1 | wCDM | DR1 BAO | no CMB | (R, l_A) CMB | Σm_ν free | prior without the cut / narrow |
|---|---|---|---|---|---|---|---|---|
| none | 2.37 | — | 0.72 | 1.96 | 1.66 | 2.12 | 1.91 | 2.37 / 2.37 |
| Pantheon+ | 2.28 | 2.02 | 0.90 | 2.02 | 1.71 | 2.12 | 2.08 | 2.28 / 2.28 |
| Union3 | 3.29 | 2.40 | 0.88 | 2.98 | 2.72 | 3.13 | 3.06 | 3.29 / 3.29 |
| DES-SN5YR | 3.84 | 1.47 | 2.08 | 3.49 | 3.26 | 3.71 | 3.71 | 3.84 / 3.84 |
| DES-Dovekie | 2.84 | 1.18 | 1.24 | 2.62 | 2.21 | 2.68 | 2.65 | 2.84 / 2.84 |

The parameter-shift metric (Hessian distance of (−1, 0) from the best fit) tracks the Δχ² metric within 0.1-0.2σ
wherever the best fit is inside the prior box; it is undefined for the 10 SN-free CPL fits that sit on the wa prior
edge, which are excluded from its ANOVA (same ranking of factors: cut 36 %, dof 11 %, sample 9 %).
Caveat on the cut: for Union3 "z > 0.1" drops its two lowest bins (binned distances at z = 0.05 and 0.10); for the
other samples it is a cut on z_HD.

### Row 2: parametric bootstrap of the Shapley values (`boot_w0wa.json`)

40 realisations: every data vector (DESI DR2 BAO, the three SN samples, the compressed CMB) drawn from its covariance
around the w0wa best fit of the five-player grand coalition (w0 = −0.81, wa = −0.64, Ω_m = 0.31), the 32 coalitions
refitted in both models from the data's best fits (one start each; 44 minutes on 4 cores). Mean ± sd of φ over the
realisations, and the 16-84 % interval of the share φ/Δχ²_grand; the data's value for comparison:

| game | player | data φ | bootstrap φ (mean ± sd) | data share | bootstrap share 16-84 % |
|---|---|---|---|---|---|
| 5-player (grand Δχ² 23.9; bootstrap 23.1 ± 9.1) | DESI DR2 | 9.64 | 6.0 ± 3.2 | 40 % | 17-41 % |
| | CMB_c | 3.08 | 2.2 ± 2.0 | 13 % | 2-15 % |
| | Pantheon+ | −0.63 | 7.4 ± 4.4 | −3 % | 21-45 % |
| | Union3 | 4.18 | 1.7 ± 1.9 | 18 % | −2-19 % |
| | DES-SN5YR | 7.59 | 5.8 ± 3.0 | 32 % | 16-35 % |
| DESI, CMB_c, Pantheon+ (data 7.6; bootstrap 13.5 ± 6.5) | DESI | 5.83 | 4.1 ± 2.4 | 77 % | 15-48 % |
| | CMB_c | 1.55 | 1.2 ± 1.4 | 20 % | 0-16 % |
| | Pantheon+ | 0.22 | 8.1 ± 4.8 | 3 % | 39-78 % |
| DESI, CMB_c, Union3 (data 13.8; bootstrap 5.3 ± 3.5) | DESI | 7.03 | 2.2 ± 2.4 | 51 % | 9-64 % |
| | CMB_c | 2.05 | 0.9 ± 1.4 | 15 % | −9-33 % |
| | Union3 | 4.75 | 2.3 ± 2.1 | 34 % | 8-77 % |
| DESI, CMB_c, DES-SN5YR (data 18.0; bootstrap 11.7 ± 5.3) | DESI | 7.99 | 4.1 ± 2.7 | 44 % | 15-51 % |
| | CMB_c | 2.17 | 1.1 ± 1.2 | 12 % | −1-15 % |
| | DES-SN5YR | 7.85 | 6.5 ± 3.1 | 58 % | 36-81 % |

Reading: the sampling error of a Shapley value is 2-3 units of Δχ² for DESI and the SN player and 1-2 for the CMB, i.e.
30-40 % of the data's values; the shares' 68 % intervals span 30 percentage points. Under the fitted w0wa model the
expected pattern differs from the data's: Pantheon+ (the most precise sample) would carry 8 ± 5 units in its game and
carries 0.2; Union3 would carry 2.3 ± 2.1 and carries 4.8; the three-player grand values of the data sit at +0.9σ
(Pantheon+ low), +2.4σ (Union3 high) and +1.2σ (DES-SN5YR high) of the bootstrap expectation. The one statement the
intervals allow in every game: the compressed CMB's share is small (16 % or less at 84 %). Bootstrap fits use one start
from the data's best fit, so part of the scatter is minimiser noise; the values are upper bounds on the errors. The ΛCDM-truth null
(`boot_lcdm.json`, 40 realisations around the grand coalition's ΛCDM best fit, same procedure): the grand Δχ² is
1.6 ± 1.5 (maximum 6.5 of 40; the three-player games 1.3-1.7 ± 1.5, maxima 5.9-9.1), none within reach of the data's
7.6-23.9; the noise floor on a Shapley value is ±1.0-1.5 units per player (φ means 0.3-1.0, DESI's the largest because
it is the player most often on its own), which is also the size of the minimiser noise in these one-start fits.

### Row 4: the full-CMB check on the headline coalitions (`desi_fullcmb.py`, `data/r123/out/fullcmb/*.json`)

Planck 2018 PR3: plik_lite TT,TE,EE (ℓ ≥ 30), Commander TT (ℓ < 30), SimAll EE (ℓ < 30), the PR3 lensing
reconstruction (conservative range), through cobaya's clik likelihoods (clipy 0.15) and CAMB 2.0.5 with PPF dark energy,
one massive neutrino of 0.06 eV; free (ω_b, ω_c, 100θ_MC, τ, ln 10^10 A_s, n_s, A_planck [Gaussian 1 ± 0.0025], w0, wa),
whitened with the Planck w0wa+BAO covariance; L-BFGS-B (0.1σ then 0.03σ finite differences, capped at 450 + 200
evaluations) and a 250-evaluation Nelder-Mead polish, from the compressed pipeline's best fit; 35-56 minutes per fit at
1.7-2.5 s per likelihood call. The χ² floor is flat to 0.02 over the last 200 evaluations of each fit.

| coalition | χ²(ΛCDM) → χ²(w0wa) | Δχ² here (σ) | DESI DR2 (PR4 CamSpec, PR4+ACT lensing) | compressed grid | w0wa best fit here (w0, wa, H0, Ω_m) |
|---|---|---|---|---|---|
| DESI + CMB, no lensing | 1018.47 → 1008.81 | 9.66 (2.7σ) | 9.7 (2.7σ) | 8.1 (2.4σ) | −0.42, −1.73, 63.7, 0.352 |
| DESI + CMB (PR3 lensing) | 1027.75 → 1017.66 | 10.09 (2.7σ) | 12.5 (3.1σ) | 8.1 | −0.47, −1.58, 64.0, 0.348 |
| DESI + CMB + DES-SN5YR | 2677.11 → 2657.23 | 19.87 (4.1σ) | 21.0 (4.2σ) | 18.0 (3.8σ) | −0.76, −0.80, 66.9, 0.318 |

DESI's own numbers reproduce to 0.1 in Δχ² without lensing and to 1.1 with DES-SN5YR (0.1σ); the lensing row is the
one place the PR3 stand-in falls short: the PR3 reconstruction adds 0.4 to the preference where DESI's PR4+ACT DR6
combination adds 2.8. The compressed (θ_*, ω_b, ω_bc) form loses 1.6-1.9 units of the full primary-CMB preference. Per
dataset at the joint best fits (χ² at the ΛCDM minimum minus at the w0wa minimum): DESI + CMB without lensing: plik_lite
+6.1, Commander TT −1.1, SimAll EE +0.4, DESI +4.2; with lensing: +4.1, −0.7, +1.4, lensing +0.3, DESI +4.7; with
DES-SN5YR: plik_lite +2.5, TT −0.2, EE +0.9, lensing +0.3, DESI +5.0, DES-SN5YR +11.3. So with the full likelihood the
primary CMB is not a bystander that is "sacrificed" by 2 units (the compressed grid's reading) but a gainer of 2-5 units:
the w0wa best fit of DESI + CMB has Ω_m ≈ 0.35 and H0 ≈ 64 and a higher lensing amplitude in the high-ℓ spectra, which
plik_lite prefers (the A_L direction). Best fits match DESI's (−0.43 ± 0.22, −1.72 ± 0.64, 63.7, 0.353 for DESI + CMB
in the compressed form; −0.76 ± 0.07, −0.80 ± 0.3, 67.3, 0.319 with DES-SN5YR). Running at the time of writing: the CMB
alone, CMB + DES-SN5YR, DESI + CMB + Union3 and + Pantheon+ pairs, for the three-player game with the full CMB as the
player (added below if they finish within the budget).

### Blind finding check, pass 1 (2026-11-01, context-free Opus referee, own literature search; `blind-review.md`)

VERDICT: novel but not sound: the Shapley table is computed on a compressed-CMB game that the note's own full-likelihood check shows misattributes the CMB, its ordering claims sit inside errors from a bootstrap that assumes the w0wa signal is real, and the variance decomposition is run on an unbalanced, nonlinear-σ grid without an expected-under-signal baseline.

Strongest objection (referee's words): the attribution is computed for a compressed-CMB game whose grand Δχ² is 2-3
units below DESI's; the full-likelihood check shows this game mis-assigns the CMB's share and the full-CMB game did not
finish; the ordering claims sit inside the bootstrap errors, and those errors assume a real signal. Also raised: the
specification counts do not close (unbalanced grid, order-dependent sums of squares); σ is a nonlinear, dof-dependent
response; a cut removes lever arm under a real signal too, so the 44 % needs an expected-under-signal baseline; no
trials correction; Shapley's additivity axiom has no physical meaning for a profile Δχ² (Banzhaf or Owen would split
the interactions differently); the BBN prior switched on only in CMB-free coalitions; prior edges in weak coalitions.
Prior work the referee added: Karim, Speagle & Hložek 2609.04334 (leave-one-redshift-block-out predictive scoring: DESI
BAO bins, not SNe, drive the tension), Cortês & Liddle 2504.15336, Huang, Cai & Wang 2502.04212, Popovic+ 2511.07517,
Ong+ 2511.10631; the Shapley and specification-curve methods themselves it found nowhere. One revision pass (below): the
presentation points fixed (claims as intervals, balanced decomposition on Δχ², the cut baseline, the prior work cited,
the conventions stated); the full-CMB three-player game added if its fits finish within the budget.

### Row 5: revision pass (`desi_r123_revise.py`, `speccurve_balanced_anova.json`, `cut_baseline.json`)

Balanced sub-grid: 4 SN samples × 2 cuts × 2 BAO × 3 CMB × 2 models × 2 priors (DESI, narrow) × 2 neutrino = 384
specifications (all present), main effects orthogonal. Shares of the variance:

| response | cut | cut × sample | sample | CMB form | BAO | model dof | ν | prior | main effects |
|---|---|---|---|---|---|---|---|---|---|
| Δχ² (CPL only, n = 192) | 55 % | 20 % | 15 % | 6 % | 2 % | — | 0.2 % | 0 % | 79 % |
| σ (CPL only) | 58 % | 19 % | 13 % | 8 % | 2 % | — | 0.3 % | 0 % | 80 % |
| σ (both models, n = 384) | 43 % | 12 % | 6 % | 1 % | 2 % | 17 % | 0 % | 0 % | 68 % |

Level means of Δχ² (CPL): cut released 10.0 vs z > 0.1 3.3; Pantheon+ 4.8, Dovekie 5.0, Union3 8.4, DES-SN5YR 8.4;
CMB none 5.1, θ_* 7.6, R-l_A 7.3; DR2 7.3 vs DR1 6.0. Same order as the unbalanced 540 grid.

Cut baseline: 12 parametric realisations of DESI + CMB_c + SN around that coalition's own w0wa best fit, each fitted
with and without the z > 0.1 cut (one start, from the data's best fits):

| sample | data Δχ² full → cut | expected full | expected cut | expected ratio | data ratio | P(ratio ≤ data) |
|---|---|---|---|---|---|---|
| DES-SN5YR | 18.0 → 3.9 | 19.2 ± 6.2 | 6.7 ± 4.2 | 0.34 ± 0.18 | 0.22 | 0.42 |
| Pantheon+ | 7.6 → 6.3 | 8.4 ± 4.2 | 5.8 ± 4.0 | 0.75 ± 0.37 | 0.83 | 0.58 |
| Union3 | 13.8 → 8.2 | 14.9 ± 9.0 | 11.3 ± 8.4 | 0.71 ± 0.21 | 0.60 | 0.25 |

Decided: the data's drops under the cut are what a real w0wa signal of the fitted size loses with the z < 0.1
supernovae (the lever arm), so the specification curve's top factor is read as "where the information is", not as a
systematic; the earlier row-3 wording ("consistent with a low-z offset") is withdrawn in the finding section.

### Row 6: the full-CMB three-player game (DES-SN5YR) and the Union3 headline (`fullcmb/*.json`)

Full Planck 2018 fits added: CMB alone 1011.80 → 1008.63 (Δχ² 3.17; w0wa best fit w0 = −0.71, wa = −2.22, H0 = 74.4,
the geometric degeneracy), CMB + DES-SN5YR 2655.98 → 2647.79 (8.19, 2.4σ; −0.76, −0.90, 67.2), DESI + CMB + Union3
1056.16 → 1040.32 (15.84, 3.6σ; −0.68, −1.02, 66.0; DESI: 17.4, 3.8σ). Coalition values for the DES-SN5YR game with the
full CMB as the player: v(DESI) 4.66, v(CMB) 3.17, v(SN) 5.88, v(DESI+CMB) 10.09, v(CMB+SN) 8.19, v(DESI+SN) 13.57,
v(all) 19.87. Shapley values: **DESI 7.88 (40 %), CMB 4.45 (22 %), DES-SN5YR 7.54 (38 %)**; compressed game 7.99 / 2.17 /
7.85. Per-dataset gains at the full best fits: CMB + DES-SN5YR: plik_lite +1.2, TT +0.5, lensing +0.4, DESY5 +6.1; CMB
alone: plik_lite +1.7, TT +0.9, EE +0.4, lensing +0.2. Running at the time of writing: CMB + Union3 (both models), for the
Union3 game; DESI + CMB + Pantheon+ (queue C).

### Blind finding check, pass 2 (2026-11-01, on the revised section, context-free Opus referee)

VERDICT: novel but not sound: the Shapley values give weight 1/3 to singleton coalitions whose profile Δχ² is set by the prior-box edges (most visibly the full-CMB-alone value behind "the CMB doubled"), the error bars assume the fitted w0wa model is true, and several stated conclusions (CMB smallest in every game, no zero player, agreement with Karim et al., "not a systematic") go beyond or contradict the note's own intervals.

Decided: the revision pass is spent and the leading objection is a method point (the prior-box dependence of the
singleton coalitions' profile Δχ², which the Shapley formula weights by 1/3 in a three-player game); recorded as a miss,
parked, the path to pass under the section in `finding.md`. The four sentences the referee named as going beyond the
tables are corrected there with dated brackets; the verdict stands on the section as refereed. Headline of this work:
not passed.

### Row 7: the full-CMB three-player game for Union3 (after the gate; not refereed)

CMB + Union3 with the full Planck 2018 likelihood: 1037.97 → 1031.63 (Δχ² 6.35, 2.0σ; w0 = −0.69, wa = −1.13, H0 = 66.7).
Coalition values: v(DESI) 4.66, v(CMB) 3.17, v(Union3) 3.40, v(DESI+CMB) 10.09, v(CMB+Union3) 6.35, v(DESI+Union3) 10.07,
v(all) 15.84. Shapley values: **DESI 6.98 (44 %), CMB 4.37 (28 %), Union3 4.48 (28 %)**; compressed game 7.03 / 2.05 / 4.75.
Same pattern as the DES-SN5YR game: the full likelihood doubles the CMB's value (1.1 of it the prior-edge CMB-alone
singleton the referee named), DESI's value unchanged, the supernova value down 0.3, the grand total up 2.0. The Pantheon+
pair (queue C) was stopped at the budget. All full-likelihood fits: `data/r123/out/fullcmb/*.json`.
