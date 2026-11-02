# Results (R-151, cosmo-deep-20, 2026-11-02)

Is thin-stream density power at each stream's own epicyclic scale, or broadband? A cross-stream test of the Ibata et al.
(arXiv:2002.01488, GD-1 power is epicyclic) and Banik et al. (arXiv:1911.02662, GD-1 + Pal 5 power needs subhalos,
n_sub/n_CDM = 0.4 +0.3 -0.2) readings (`open-problems/round2/p8-streams.md`; trajectory check "partly novel": only the
cross-stream epicyclic-vs-broadband test is open). This file is the worker's own. Step 0 is the feasibility gate: forecast
the separation; below 2.5 sigma, stop and report.

Scripts: `r151_epicycle.py` (predicted epicyclic scale per stream), `r151_pepper.py` (galpy streampepperdf mocks,
subhalo-perturbed density), `r151_forecast.py` (the separation forecast). Data and mocks under `data/r151/` (gitignored),
logs `logs/r151_*.log`.

## F0. Access (fetch times; `data/` started empty)

| what | where | time |
|---|---|---|
| arXiv PDFs 2002.01488, 1911.02662, 2603.09305, 2206.04241, 2510.07511, 2110.06950 | arxiv.org/pdf | 1-2 s each |
| galpy 1.12 | PyPI | ~1 min |
| `streampepperdf` (not in any released galpy: never merged into galpy's main line) | github.com/jobovy/streamtools (`galpy.util.bovy_*` imports patched to `galpy.util.conversion/coords`) | 2 s clone |
| Banik et al.'s own GD-1 analysis code and density files | github.com/nbanik/Analyze_GD-1 | 2 s clone |

No survey photometry was downloaded in step 0: the member counts and noise come from the published density profiles
(Patrick et al. 2206.04241, Boone et al. 2510.07511, Ibata et al. 2002.01488), which is all a forecast needs.

## F1. Predicted epicyclic spacing per stream (`r151_epicycle.py`, `r151_epicycle.json`)

Kuepper et al. (2008, 2012) spacing, as Ibata et al. use it (their eqs 12-13): y_C = (4 pi Omega/kappa)
|1 - 4 Omega^2/kappa^2| x_E, x_E = 2.88 (G M / (4 Omega^2 - kappa^2))^(1/3). Omega and kappa in MWPotential2014 (Banik's
potential) at pericentre and at the guiding radius (r_peri + r_apo)/2; orbits from Li et al. 2022 (S5, arXiv:2110.06950
Table 2), GD-1 from Ibata et al. Progenitor mass: GC streams 1e4-1e5 Msun; dwarf streams their total (dark) mass, 1e7-1e9
(Indus 3e8 in Yang et al. 2603.09305). Window = the stream's observed extent (Patrick et al. 2206.04241 lengths).

| stream | type | r_peri-r_apo (kpc) | y_C (kpc) | y_C (deg along the window) | window | cycles in window |
|---|---|---|---|---|---|---|
| ATLAS-Aliqa Uma | GC | 12.1-38.2 | 1.39-5.11 | 3.3-12.1 | 37 deg / 15.6 kpc | 3.1-11.2 |
| Phoenix | GC | 12.7-19.1 | 1.44-3.65 | 2.7-6.7 | 14 deg / 7.4 kpc | 2.0-5.1 |
| GD-1 | GC | 14.0-26.0 | 1.55-4.30 | 10.3-28.6 | 80 deg / ~12 kpc | 2.8-7.7 |
| Jhelum | dwarf | 9.1-27.8 | 11-87 | 61-481 | 29 deg / 5.3 kpc | 0.1-0.5 |
| Indus | dwarf | 12.3-17.9 | 30-76 | 110-275 | 20 deg / 5.6 kpc | 0.1-0.2 |
| Turranburra | dwarf | 8.3-35.2 | 10-98 | 21-203 | 17 deg / 8.2 kpc | 0.1-0.8 |

Check: GD-1's measured 2.64 +- 0.18 kpc (Ibata et al.) sits inside its predicted 1.55-4.30 kpc, at M ~ 3-5e4 Msun.

Two facts here shape the whole test:
1. **The dwarf streams have no testable epicyclic scale in their DES windows.** With a dark progenitor the formula puts
   the spacing at tens of kpc, longer than the whole window. Indus's N-body peaks (Yang et al.) are ~20 deg apart, so
   its 20 deg DES window holds about one, and that spacing comes from a simulation of the late, mostly stripped dwarf,
   not from the orbit-and-mass formula the method fixes in advance. Jhelum and Turranburra are the same, and Jhelum is
   also perturbed by Sagittarius. For all three, "the epicyclic band" is the fundamental mode of the window, which is
   exactly where broadband (subhalo) power is largest, so they cannot separate the readings.
2. **The GC streams' bands overlap in kpc.** y_C goes as (G M)^(1/3) (r/v_c)^(2/3): the stream-to-stream differences
   from orbit radius (factor ~1.2-1.7) are smaller than each stream's own spread from the unknown progenitor mass
   (factor 2.15 per decade of mass). ATLAS-AU, Phoenix and GD-1 all predict ~1.4-5 kpc. Different orbits predicting
   different scales is the lever the proposal relied on, and this sample barely has it. In degrees the bands differ,
   but only through distance, and that is a projection, not dynamics.

## F2. Member counts and 1D density noise (from published profiles; no fetch needed)

| stream | excess members in the window | background in the window | per-bin fractional noise (40 bins) | source |
|---|---|---|---|---|
| GD-1 | 1500 | 1500 | 0.23 | Ibata et al. 3-sigma Gaia+PS1 samples (868-1171); checked against the Gaia+PS1 profile in Banik et al.'s repo (`fromDenis/GD1_stream_track_dens.txt`): 17-18 % error per 2 deg bin over phi1 = -60..-10, i.e. N_eff ~ 800 background-free, the same noise as (1500, 1500) |
| ATLAS-Aliqa Uma | 300 | 900 | 0.73 | stellar mass 1.14e4 Msun (Patrick et al.), 3.3x Phoenix at ~1.2x the distance |
| Phoenix | 100 | 300 | 1.26 | DES r < 24.2-24.5: 5-7 stars/deg over 14 deg (Boone et al. 2510.07511 Fig 10; Tavangar et al. 2022) |
| Jhelum, Indus, Turranburra | n/a | n/a | n/a | no epicyclic band inside the window (F1) |

Gaia DR3 proper motions reach only G < 20.5-21, the brightest few percent of these DES members, so they barely cut the
background at DES depth. The background-to-member ratio of 3 for the DES streams is the matched-filter contrast in Boone
et al.'s synthetic test (a 25 % excess) and the Phoenix maps.

## F3. Forecast (`r151_pepper.py`, `r151_forecast.py`; logs `r151_forecast_v2.log`, `r151_forecast_counts.log`)

**Mocks.** galpy streampepperdf (jobovy/streamtools), Bovy et al. 2017 / Banik et al. conventions: MWPotential2014,
trailing arm, 5 Gyr, CDM subhalo density per mass decade 0.3 (10^6.5/M) per (8 kpc)^3, Hernquist r_s = 1.05 kpc
(M/1e8)^0.5, 1e5-1e9 Msun, impacts to 5 r_s, sigma = 120 km/s, 8 impact times. GD-1 (peri 14, apo 26 kpc, arm 15.4 kpc):
20.1 impacts expected at CDM; 60 realisations each at n/n_CDM = 0.2 and 1.0; delta_rms over the arm 0.08 at CDM.
ATLAS-AU's eccentric orbit fails galpy's stream-track setup (`ValueError: Found time value not in the integration time
domain`, twice, at 26 and 11 track chunks), and the Phoenix setup did not finish in ~25 min, so both use random segments of
the GD-1 mock cut to their arm lengths (scales kept in kpc). Their orbits are similar (peri 12-14 kpc), and only the shape
of the subhalo spectrum enters a band fraction.

**Estimator (as the method would fix it).** Members drawn from (1 + delta) plus a uniform background, background
subtracted, normalised by a cubic (Banik et al.), Hann-tapered periodogram in path length, modes with wavelength >= 1 kpc,
shot noise subtracted. Statistic f = excess power in the stream's predicted epicyclic band / excess power in all modes;
combined f = the mean over streams. Separation z = sqrt(2) Phi^-1(AUC), AUC = P(f_epi > f_sub). The draws are 2000 per
model per stream. The ratio has heavy tails when the total excess nears zero, so a mean/std z understates the
separation; the AUC does not.

**Epicycle model, chosen to favour separation.** All the excess power is put in one sinusoid at y_C (log-uniform in the
predicted band, random phase), with the same total variance as the subhalo field. Harmonics, envelopes and washing out
away from the progenitor spread a real epicyclic pattern's power. Bar, GMC and Sagittarius perturbations add broadband
power to both families. Both effects can only lower z, so these numbers are upper bounds.

| amplitude (total delta_rms) | subhalo shape | GD-1 | ATLAS-AU | Phoenix | **combined** |
|---|---|---|---|---|---|
| CDM-native (0.08) | n/n_CDM = 1 | 0.27 | 0.02 | 0.00 | **0.18** |
| 0.1 | 1 | 0.48 | 0.03 | 0.05 | **0.34** |
| 0.2 | 1 | 1.26 | 0.07 | 0.07 | **0.66** |
| 0.3 (Ibata's high-contrast peaks) | 1 | 1.50 | 0.35 | 0.15 | **0.81** |
| 0.3 | 0.2 | 1.45 | 0.37 | 0.17 | **0.86** |
| 0.3, no noise at all | 1 | 2.27 | 2.09 | 1.35 | **3.22** |
| 0.3, no noise at all | 0.2 | 2.35 | 2.39 | 1.55 | **3.42** |
| 0.3, members x4, background = members | 1 | 1.68 | 1.49 | 0.92 | **1.64** |
| 0.3, members x10, no background | 1 | 2.08 | 1.67 | 0.98 | **1.89** |

**Shuffled-band null (what it could test).** Scored in another stream's band, each stream keeps most of its own band modes:
GD-1 in ATLAS-AU's band 5/5, in Phoenix's 4/5; ATLAS-AU in GD-1's 7/8, in Phoenix's 6/8; Phoenix in GD-1's 2/3, in
ATLAS-AU's 3/3. The bands are nearly the same set of modes, so swapping scales cannot remove a signal. The cross-stream
lever does not exist in this sample.

## F4. Gate verdict: stop (forecast separation <= 0.9 sigma, below the 2.5 sigma line)

With the data that exist, the joint test cannot separate "power at each stream's own epicyclic scale" from "broadband
power". The combined separation is at most 0.8-0.9 sigma, and that is with the epicycle model built to favour
separation. Three things set this, in order:
1. **Only the GC streams have a testable band, and their bands overlap** (F1). The dwarf streams' predicted spacing is
   longer than their windows. The GC streams all predict 1.4-5 kpc because y_C ~ M^(1/3) (r/v_c)^(2/3): the unknown
   progenitor mass (x2.15 per decade) outweighs the orbit differences (x1.2-1.7). So the band covers 3-8 of each
   stream's 7-15 usable modes, a broadband spectrum already puts 0.2-0.6 of its power there, and the shuffle null has
   nothing to remove.
2. **The DES streams are shot-noise dominated** at any plausible amplitude: ATLAS-AU and Phoenix add <= 0.4 sigma each.
   Nearly all the information is GD-1's, the stream the dispute is already about.
3. **The mode count caps it even without noise**: 3.2-3.4 sigma combined, and with 10x the members and no background,
   still 1.9 sigma.

What would change the answer (for a future proposal, not this job): GC streams with a **surviving progenitor of measured
bound mass** (Pal 5, NGC 5466, M68/Fjorm, NGC 3201/Gjoll), where y_C is predicted to about +-20 % from the orbit alone,
and whose bands differ by more than that. But Banik et al. already concede epicyclic structure near a living progenitor,
and the dispute is about the old, mixed parts of a dissolved stream. A test that is sharp only near progenitors does not
reach the disputed regime. That tension is the real obstacle, and it is a finding of this gate.

Decided here (loose ends): no survey photometry was fetched and no stream power spectrum was measured (step 1 would
have been the first). No `finding.md` section or finding check was run, because there is nothing to claim, as with
R-148. The Jhelum, Indus and Turranburra mocks were not run because their bands fall outside their windows whatever the
mocks show. A units slip in the first GD-1 run (Bovy's `nsubhalo` multiplied by ro^3, giving 10,286 impacts) was caught
from the impact count and fixed before any forecast used it.
