# Round 8, harvester C: unretested claims in parameter tensions, their claimed sources, 21 cm / PTA, and catalogue oddities (2026-11-02)

Brief: claims at >= 3 sigma (2018-2026) that nobody has re-tested on independent data, where public data now allow a
re-test that separates "real" from "chance" at >= 3 sigma (the Giant Arc shape, R-155). Two gates checked before listing:
enough independent objects in the claimed volume, and the new data measure the same quantity. Search and reasoning only;
no code run, no data downloaded. Skipped as tried: everything in rounds 1-7 and the "Tried by us" list in
`../../open-problems.md` (SN host age E, DES5Y low-z offset, local void, bulk flows, dipoles, quasar Hubble diagram, MPTA
amplitude and hotspot, PTA x LSS, EDGES, P1D running, BAO caps).

## Survivors: none

The claims in these areas fall into three groups, none of the Giant Arc shape:

1. **Re-tested within ~2 years by the survey teams**, usually with overlapping data and a settled answer (DES5Y/Pantheon+
   offset; Dainotti's H0(z) on DES; ST polarisation in NG15; CPTA's HD evidence by the other PTAs as a whole).
2. **Model-mediated claims with no frozen template**: the "anomaly" is a fit parameter (H0(z) slope, M_B transition,
   Omega_m running), so new data move the parameter but cannot separate "real" from "chance" without the same model
   choices.
3. **No independent data at the claim's scale**: Cepheid-calibrator transitions need new calibrators; 21 cm claims need
   raw-visibility reprocessing; a per-bin PTA claim at 14 nHz lies where the other public PTAs have little signal.

## Screened and dropped

| Claim (arXiv, significance) | Why dropped | Queries |
|---|---|---|
| CPTA DR1 Hellings-Downs at 14 nHz, 4.6 sigma (2306.16216), from one frequency (1.5/T, T = 3.4 yr) | **Fails 5.** No bin-by-bin re-test exists, but in NG15/EPTA DR2/PPTA DR3 (public) 14 nHz is bin ~7. There the GWB power (falling as f^-13/3) is below the white- and red-noise floor, and the published free spectra show the signal only in the lowest ~5 bins. At CPTA's amplitude (log10 A ~ -14.4, near NG15's), the HD S/N in that one bin from public data is ~1 sigma (estimate), so it cannot separate the readings. The other PTAs' broadband HD evidence (3-4 sigma) already answers the physical question. | WebSearch (ext) "CPTA 14 nHz ... independent test"; review 2603.13643 |
| Running H0(z) in SN Ia, alpha ~ 0.01 at 2-5.8 sigma (Dainotti 2103.02117, 2411.05744; master sample 2501.11772) | **Fails 2 and 4.** Already re-run on DES, JLA and Pantheon+ (2501.11772: 1 sigma to > 3 sigma, depending on sample). The slope is degenerate with M_B, Omega_m and w, so it isn't a frozen template; "chance" isn't separable from "model choice". | WebSearch "Dainotti ... H0(z) ... DES"; 2602.22840, 2510.14390 |
| SH0ES calibrator transition (M_B / Cepheid R_W at ~20-50 Mpc), Perivolaropoulos & Skara 2109.04406 (3 sigma in Cepheid params), 2208.11169 (Delta chi2 ~ -15), 2504.19807 (withdrawn) | **Fails 3.** The transition sits in the calibrator sample (~40 hosts). Independent objects in the claimed range need new Cepheid/TRGB hosts. ZTF/ATLAS Hubble-flow SNe measure M_B only through H0, not the calibrator step. | arXiv API `abs:transition AND abs:SnIa AND (absolute OR magnitude) AND abs:Mpc` (3 hits, all SH0ES-only) |
| SN Ia Hubble residuals vs large-scale environment (voids/filaments/groups) | **Fails 1.** No >= 3 sigma claim exists. The 2026 TITAN "cousins" (2609.13525) is a scatter result at 2.3-2.8 sigma; the cluster vs field result (2306.01088) shows minor differences; local-environment steps (LsSFR, local colour, age) are re-tested many times (ZTF DR2 2605.06799, HST/MUSE 2606.22173: no trend > 3 sigma). | WebSearch (ext) "Hubble residuals environment ... filaments voids"; arXiv API `abs:supernova AND abs:"Hubble residuals" AND (environment OR density OR filament OR void)` (21 entries) |
| DES5Y vs Pantheon+ 0.04 mag low/high-z offset; Union3 w0wa 3.8 sigma | **Fails 2.** Traced to selection bias and host-mass choices (2501.06664, 2601.19424), combined in Unite (2609.05053). Same SNe, not independent objects. | WebSearch (ext) "SN Ia subsample anomaly ... Pantheon+ DES-SN5YR Union3" |
| DESI DR1 LRG z = 0.51 BAO bin (~3 sigma off BOSS/Planck) | **Fails 2.** DR2 moved it back toward BOSS. There is no independent spectroscopic volume at that z beyond BOSS, which is the side of the disagreement already measured. | reasoning (known DR2 result) |
| Non-tensorial (ST) correlation in NANOGrav 12.5 yr (Chen et al. 2101.06869, Bayes factor ~100) | **Fails 2.** Re-tested by NANOGrav 15 yr (no evidence for non-Einsteinian polarisations). | reasoning (known literature) |
| MeerKAT interferometric 21 cm auto-power detection at z = 0.32/0.44 (Paul et al. 2301.11943, ~8-11 sigma); CHIME 21 cm x eBOSS | **Fails 3.** The independent data are raw visibilities that need a full foreground pipeline, beyond a small team. CHIME x DESI already extends CHIME x eBOSS. These are detections, not anomalies, so "chance" isn't the live reading. | WebSearch (ext) "tentative detection 21 cm ... cross-correlation" |
| ZTF SN Ia DR2 low-z (0.02-0.04) excess at 5-7 sigma (Lordet 2604.12714) | **No dispute.** The claim paper attributes it to known local structure. A re-test on post-2020 BTS SNe would confirm known LSS, not separate a contested reading. | (01-h0 map row 1) |
| Lya mean-transmission dip at z ~ 3.2 (He II reionisation) | **Outside the window**: claims date from 2003-2013, and later results disagree. DESI DR1 could re-test it, but it isn't a 2018-2026 claim. | reasoning |
| Generic "unexpected/anomalous ... sigma" in DESI/eROSITA/ZTF/Euclid survey papers | **No new claim found.** The 2026 local-void and bulk-flow papers (2602.03928, 2602.21975) are model fits to structures we already tried. | WebSearch (ext) "anomalous large structure ... DESI OR eROSITA OR Euclid"; arXiv API `(DESI OR eROSITA OR ZTF OR Euclid) AND (unexpected OR puzzling OR anomalous) AND sigma` (the query parsed badly: 1 off-target hit) |

Search budget: 6 WebSearch (5 extended), 3 arXiv API queries, prior-round maps read for overlap.

## Note for the next round

Tension "sources" are fit parameters, not positional templates, so they don't fit the Giant Arc shape. In this area, the
independent data that would reopen it are new calibrator hosts (JWST Cycle 4-5 SH0ES/CCHP releases), IPTA DR3 (for the
CPTA 14 nHz bin, if its sensitivity there improves enough), and Gaia DR4 (December 2026: the quasar proper-motion
quadrupole and anisotropic expansion).
