# Results (R-124, cosmo-deep-4, 2026-11-01)

NANOGrav 15-year continuous-wave limits on the periodic-quasar binary candidates, as a population. This file is the
worker's own. Script: `periodic_cw_limits.py` (usage in its docstring); data under `data/ng15/` and `data/periodic/`
(gitignored). Reachability from the cloud (contrary to the brief's expectation): zenodo.org answers 200 and serves the
NANOGrav 15-year releases (timing data 16051178, the CW analysis products 8067506 with the pickled `enterprise` pulsars and
the sky-and-frequency 95 % limit map, the noise spectra 8092346); data.nanograv.org answers but its 15yr page is a shell;
every VizieR and CDS host (cds.unistra.fr, cfa.harvard.edu, nao.ac.jp, idia.ac.za, the TAP endpoint, cdsarc ftp) is cut by the
proxy mid-handshake ("tunnel closed"); OUP and IOP table hosts answer 403/404; arXiv e-prints answer 200, so the candidate
tables were taken from the papers' LaTeX sources.

## Trajectory check (before any compute; `blind-review.md`, context-free Opus, 2026-11-01)

| proposal | verdict (verbatim) |
|---|---|
| Implied strain per candidate (period, z, mass; mass ratio and mass error marginalised) against the 15-year sky-dependent 95 % limits (cheap pass), then fixed-position fixed-frequency `enterprise` upper limits on the 20-50 strongest with scrambled-sky nulls (targeted pass); the claim the table and the population fraction excluded with its interval | VERDICT: partly novel: NANOGrav 15-yr targeted CW search of 114 periodic-AGN SMBHB candidates (mostly CRTS) with per-target strain/chirp-mass 95% limits, Bayes factors and random-sky nulls, arXiv:2508.16534 (and implied-strain vs PTA sensitivity for 149 CRTS/PTF/PS1 candidates, arXiv:2009.11865); what remains open is an explicit population excluded-fraction from EM-mass-vs-limit comparison with mass/mass-ratio uncertainties propagated, coverage of the PTF/PS1-MDS/DES/ZTF catalogues if absent from the 114, and the period-doubling (f_gw = 1/P) case |

Read in full afterwards (the referee could not open the PDF): Agarwal et al. 2025 (arXiv:2508.16534, ApJL 998 L11, v2
2026-01-13) is NANOGrav's own targeted search on the 15-year data for 114 targets: the 111 CRTS candidates of Graham et al.
2015, the two OVRO blazars (PKS 2131−021, PKS J0805−0111) and 3C 66B; fixed sky position, distance and f_gw = 2/P_obs;
95 % strain and chirp-mass limits per target with the EM total mass beside them (Tables 6-8); Bayes factors (CURN and HD
backgrounds; mean 0.73 ± 0.32 under HD, all consistent with noise after trials); sky-shuffle and dropout tests on the two
mild outliers (J1536+0441, J0729+4008); the only exclusion a part of 3C 66B's parameter space. It does not cover PTF, PS1
MDS, DES+SDSS or ZTF, does not state a population fraction, and does not treat the period-doubling case. Step 4 of the
brief (the `enterprise` targeted pass) is therefore published for the catalogue that matters; the open part is the cheap
pass on the rest and the tally. Decision (agent, before computing): the open part is computed because it is cheap (no
PTA likelihood, the released limit map and the published targeted limits), and whether it stands alone as a result is put
to the finding check rather than decided here.

## The cheap pass

| # | test | script | statistic | significance | what it means |
|---|---|---|---|---|---|
| P0 | Convention check against Agarwal et al.'s own pairs: the chirp mass that follows from each target's h95, f_gw and D_L (flat LCDM, H0 = 70, Om = 0.3) against their quoted Mc95 | `periodic_cw_limits.py --check` | Mc(h95)/Mc95 = 0.659 median, 0.51-0.88 (16-84 %); times (1+z): 1.244, 1.240-1.250 | the ratio is constant to 1 % once multiplied by (1+z) | their Mc95 is the redshifted (observer-frame) chirp mass up to a fixed 24 % convention factor (their D_L or h0 normalisation); the strain limits h95 are used directly and the EM masses are redshifted on our side, so the comparison does not depend on it. The released sky map's limit at the 114 positions and frequencies is 1.09 × their targeted h95 (median; 0.49-7.96): the targeted limits are about equal to the sky-dependent all-sky map, and their quoted 2.2 × improvement is against the sky-averaged curve |
| P1 | CRTS (Graham et al. 2015) through the published targeted limits: implied h0 per target from the EM total mass (0.4 dex), q log-uniform 0.1-1, f_gw = 2/P_obs, against Agarwal et al.'s h95 at that target; 101 of the 113 appendix rows have an EM mass (3C 66B is tabulated separately in the paper and left out) | `periodic_cw_limits.py` | excluded at 95 %: **0 of 101** (Clopper-Pearson 95 %: 0-3.6 %); expected number excluded with the nuisances marginalised 0.93; h95 / median h0: min 5.8, median 210; Mc95 / Mc_EM(q = 1): min 4.5, median 29, none below 1; 8 targets with P_excl > 0.05 | nothing excluded; the closest targets (SDSS J131706.19+271416.7, J134855.27−0321, J140704.43+273556.6, HS 1630+2355, J114857.33+160023.1, J160730.33+144904.3, all logM 9.8-9.95 at z 0.8-2.7) have P_excl 0.06-0.13 | the tally Agarwal et al. did not state: their own numbers exclude no CRTS candidate at its claimed parameters; the limits sit a factor 6-15 above the strongest candidates' median implied strain and 4.5-30 above in chirp mass |
| P2 | The same through the released all-sky sky-and-frequency map (22 bins to 29.7 nHz, 192 pixels), f_gw = 2/P_obs and f_gw = 1/P_obs (period doubling) | `periodic_cw_limits.py` | 2/P: 0 of 101 (0-3.6 %), expected 0.82, UL/h0 min 6.0, median 226; 1/P: 0 of 101, expected 0.55, min 10.6, median 297 | the period-doubling case is further from exclusion (h0 falls by 2^(2/3) and the map's limits do not improve as fast between 14 and 7 nHz) | neither frequency convention changes the answer |
| P3 | PTF (Charisi et al. 2016, 50 candidates with masses and observed periods, 33 surviving their extended light curves), PS1 MDS (Liu et al. 2019, 26), DES+SDSS (Chen et al. 2020, 5), ZTF (Chen et al. 2024, the 10 of the 11 printed rows of Table 1 with a mass; the full 184-row table is behind OUP and CDS, both unreachable) against the sky map; above 29.7 nHz the last bin's map scaled by f / 27.5 nHz (white-noise regime), flagged | `periodic_cw_limits.py` | PTF: 0 of 50 (0-7.1 %), expected 0.09, UL/h0 min 4.3 (SDSS J145859.07+153144.7, logM 10.03 at z 2.9, f 73 nHz, extrapolated limit), among the 33 survivors min 22.4, expected 0.01; PS1: 0 of 26 (0-13 %), expected 0.02, min 27; DES: 0 of 5, min 387; ZTF (10): 0 of 10, min 12 | the non-CRTS catalogues are higher-frequency (23-177 nHz at 2/P) and lower-mass; every candidate's median implied strain is 4-1000 × below the limit | the catalogues Agarwal et al. did not cover add no exclusion and could not: at their frequencies the 15-year limits are 1-3 × 10⁻¹⁴ and their masses 10⁸-10⁹ M⊙ |
| P4 | Pooled: all catalogues, both conventions | `periodic_cw_limits.py` | 0 of 192 candidates excluded (0-1.9 %), expected 1.1 (2/P, with the CRTS targeted limits) with the nuisances marginalised, 0.7 at 1/P | a population statement of "none" with an upper bound of 2 % | with the strongest candidates' limits a factor ≥ 4-6 above the median implied strain, the targeted `enterprise` pass of the brief (a further factor ~2 at best over the sky map, per Agarwal et al.) cannot reach an exclusion; not run |

The 192 pooled is 101 + 50 + 26 + 5 + 10. Per-candidate table: `data/periodic/candidates_cw.csv` (gitignored; regenerated by
the script in under a minute).
| P5 | The tally note (the new last section of `finding.md`) through the blind gate, the section alone in a file | Agent tool, Opus | referee's verdict (verbatim): "VERDICT: sound but not novel: Xin, Mingarelli & Hazboun 2021, arXiv:2009.11865 (with NANOGrav targeted search Agarwal et al. 2025, arXiv:2508.16534)" | MISS (not novel); the row is closed | the referee's reading: the per-target null is published (arXiv:2508.16534) and the population comparison was published for the 11-year data (arXiv:2009.11865; Sesana et al. 2018 already say each CRTS and PTF candidate is individually compatible with the PTA limits), so the 15-year update on the other catalogues is "a mechanical update"; presentation faults named, not fixed (the section stays as reviewed): the Clopper-Pearson interval is the wrong statistic for a deterministic tally with per-candidate probabilities, the "two to three data sets away" forecast misreads the targeted paper's factor 2.2 (which is targeted-vs-sky-averaged, not per generation), and the constant 1.24 in the Mc95 convention check is 2^(3/10), i.e. a √2 strain normalisation that would move the closest margin from 5.8 to 4.1 without an exclusion |

Outcome: no gate pass. The question of the brief is answered by the literature (the targeted pass) and by this tally (the open
part): no published periodic-quasar candidate is within reach of the 15-year limits at its published mass; the population
test of these candidates is the background (Sesana et al. 2018; Casey-Clyde et al. 2025), not continuous waves. The
`enterprise` runs of step 4 were not started (row P4). Budget used: about 2.5 h of 6.
