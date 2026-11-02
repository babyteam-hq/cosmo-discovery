# Results (R-153, cosmo-deep-21, 2026-11-02): ALP decay line in SPHEREx QR2 images of dwarf spheroidals

**Stopped at the Step 0 feasibility gate.** On real QR2 images, the expected 95% limit at 1.5-2.5 eV is
1.7-2.1e-11 GeV^-1. That is the eight-dwarf combination, with an optimistic bound for six of them. It does not reach
below the weaker disputed JWST limit (arXiv:2609.10674: 2.2e-11 at 1.5 eV, 1.3e-11 at 2.0 eV, 1.2e-11 at 2.5 eV)
except at 1.5 eV, where the two are equal. It is 30-50x above the strong claim (arXiv:2503.11753, ~5e-13). So the
search cannot arbitrate between them. No limit was computed, the method check (1b) was not run, and nothing went in
`finding.md`.

This file is the worker's own. Scripts:
- `alp_list.py`: the image lists.
- `alp_reduce.py`: one L2 image to (sky cell x ln-lambda bin) sums.
- `alp_pattern.py`: the fixed detector pattern.
- `alp_run.sh`, `alp_prun.sh`: launchers.
- `alp_model.py`: line physics and D-factors.
- `alp_cells.py`: OFF subtraction, per-cell continuum, cleaning.
- `alp_search.py`: templates, stellar nuisance, sign-flip errors, profile likelihood.
- `alp_feas.py`, `alp_gate.py`: the Step 0 numbers.
- `alp_diag_*.py`: the diagnostics quoted below.

Data under `data/alp/` (gitignored).

## Trajectory check (run in the dispatching session, 2026-11-01; not rerun)

"VERDICT: partly novel: whether SPHEREx dSph data can reach g_agamma of a few 1e-12 over 0.5-3 eV is already forecast
(arXiv:2412.12286), and the JWST order-of-magnitude disagreement is already set out and traced to instrument modelling
(arXiv:2609.10674, vs 2503.11753/2310.15395/2503.14582); what remains open is an actual decay-line search in SPHEREx Quick Release
data with same-exposure ON-OFF, which no paper has done, as a JWST-independent arbiter at about 1-3 eV"

## What the papers give (read in full by a sub-agent; numbers marked "fig." read off figures, +-10%)

- **2412.12286** (Regis, Taoso & Terol Calvo).
  - Targets and D-factors: the eight dSphs of Song, Hiroshima & Murase (arXiv:2401.15606, Table 2: J^dec within
    θ_Δ ~ 0.2-2.3 deg). Boötes I is not among them.
  - Noise: 27 nW m^-2 sr^-1 per pixel per FWHM channel at 0.75 µm (2-yr all-sky, photon + read).
  - Statistic: χ² = Σ Φ²/σ² = 3.84. A forecast, with no systematics: "background subtraction" and "masking" are not
    discussed, and flat field and starlight are not mentioned.
  - dSph forecast (fig.): 3.7e-11 (1 eV), 2.7e-11 (1.5), 2.1e-11 (2), 1.7e-11 (2.5), 1.4e-11 (3 eV). Draco and UMi
    each come within ~1.4x of the combination.
- **2609.10674** (Gemmell, Dessert, Caputo & Foster).
  - NIRSpec final limits (fig. 9 right): 2.3e-11 (1.0 eV), 2.2e-11 (1.5), 1.3e-11 (2.0), 1.2e-11 (2.5). Before
    their spurious-signal correction: 1.2e-11, 1.2e-11, 7.2e-12, 6.4e-12. NIRSpec ends at ~2.65 eV.
  - Pinetti 2503.11753 as drawn there: ~6.7e-13, 5.6e-13, 6.7e-13, 4.2e-13.
  - Reference ids confirmed (HTML numbering): [49] 2310.15395, [63] 2503.11753, [68] 2311.04987, [69] 2503.14582.
  - The disagreement is traced to "ten datasets with anomalously small reported uncertainties". Sensitivity saturates
    at N_obs ≳ 25.
- **MUSE-Faint** (2307.07403: 2.7-5.3 eV; 2310.01996: 2.65-5.27 eV): no numerical g in the abstracts.
- **D-factors.** Pace & Strigari (arXiv:1802.06811, Table 4) give log10 D(0.5°) per dwarf; the values are in
  `alp_model.py`. Draco is 18.54 (+0.11 -0.14) and Ursa Minor is 18.20 (+0.14 -0.08).

## Step 0: feasibility (2026-11-02, 06:10-07:55Z)

**Fetches.**
- Image lists (IRSA TAP `spherex.obscore` + `spherex.artifact`, one survey week at a time) took 32 min. The full
  artifact table drops the connection through the proxy (`SSLEOFError`), so it has to go week by week.
- QR2 L2 images within 2.3° of each dwarf, all six detectors:

  | dwarf | images |
  |---|---|
  | Draco | 486 |
  | Ursa Minor | 519 |
  | Sculptor | 542 |
  | Fornax | 492 |
  | Sextans | 318 |
  | Carina | 583 |
  | Leo I | 307 |
  | Leo II | 407 |
  | Boötes I | 495 (not in Song+2024; not used) |

- Downloads run at ~3 s per 72 MB image per stream. Reduction is CPU-bound at ~1 image/s on 4 cores. Draco's 486
  took 8 min. The pattern pass over 1,228 images (every third image of each list) took 15 min.

**Line physics.**
- Γ = g² m³/64π. I_line = Γ D_Ω/4π (energy intensity; two photons of m/2).
- Each pixel's passband is a Gaussian with FWHM equal to the pipeline bandwidth. The pipeline gives λ/Δλ = 41.6 in
  D2, from the per-pixel 'W' WCS.
- Velocity dispersion is negligible; the systemic velocity shifts the line by v/c.
- νI_ν at the peak is 0.94 R I_line. At m = 2 eV and g = 1e-11, Draco's centre has 0.4 nW m^-2 sr^-1.

**The ON-OFF and its noise, as measured.**

1. **First pass: rings.**
   - Method: ON (r < 0.6°) minus OFF (1.0-1.75°), in the same exposure and the same ln λ bin (0.002 wide). The
     stack was fitted with a local quadratic plus the line.
   - The per-exposure ON-OFF scatter is 11x the pipeline VARIANCE (χ²/dof = 130 about the stack mean).
   - The line amplitude's bootstrap error over exposures is 7-15x the white one.
   - Draco alone: g95 ≈ 2.5-3e-11 at 1.5-2.5 eV.
   - The excess is the same in a starlight-free ring (0.6-1.0°), so it is not the dwarf's light.
   - This was a superseded first reduction; that version of the reducer was not kept.
2. **Per-sky-cell continuum (0.05° cells).**
   - Each cell is seen at ~80 wavelengths per detector, as the LVF moves it. A Legendre polynomial of degree 4 in
     ln λ per detector removes fixed sky structure.
   - It did not lower the noise: χ²/dof is 15 before and after the polynomial.
   - The excess is coherent across the ON cells of one (exposure, bin): var(mean)/(var/n) = 13.5.
   - It is not removed by a per-exposure plane or quadratic background fitted to the OFF cells.
3. **Tails** (`alp_diag_noise.py`).
   - 1.2% of entries carry 72% of the χ², with symmetric tails. The robust core is 1.8x the pipeline σ.
   - A 5-robust-σ clip about each cell's continuum drops 0.3-3% of entries per detector. Each detector's variance is
     then rescaled by its robust χ² (1.7-6x).
   - The clip barely changes the line error: the remaining excess is coherent within exposures.
4. **A fixed detector-frame pattern** (`alp_diag_pattern.py`, `alp_diag_pattern_onoff.py`).
   - Method: remove each image's sky spectrum (the median of each ln λ bin), then average in detector coordinates.
   - Two half-stacks of 40 Draco D1 images correlate at 0.98, with 35 nW rms in 8x8-pixel blocks. That is nearly all
     of a single image's 37 nW block rms.
   - The pattern repeats across sky groups: north/south/ecliptic correlation is 0.96-0.98 in D1 and 0.72-0.94 in
     D2-D6. Its amplitude partly scales with zodi.
   - Iso-wavelength lines are curved, by about 1% in λ across a row, so ON and OFF at one wavelength sit on
     different detector rows. A pattern that differs between them enters every exposure's ON-OFF.
   - Subtracting the pattern, built from the other half of the images, cuts the per-image ON-OFF scatter from 36 to
     13 nW in D1.
   - The reduction subtracts this pattern scaled per image (`alp_reduce.py --pattern`). The pattern is built from
     every third image of the seven *other* dwarfs, so a dwarf's line cannot leak into its own correction.
5. **Final Step 0 chain** (`alp_search.py`, `alp_feas.py`).
   - The model is the per-cell continuum, plus the line template, plus a stellar nuisance. The line template is
     (D_Ω(cell) − D_Ω(OFF)) × the passband. The stellar nuisance is each cell's own fitted ON-OFF continuum × the
     same passband at λ0, with a free amplitude per mass.
   - Errors come from 1,000 random sign flips of whole exposures' contributions.

**Draco and Ursa Minor, expected g95 (A = 0) per mass range, pattern subtracted.**

| dwarf | m_a (eV) | from pipeline variance (photon + read) | white, robust-rescaled | sign-flip (real noise) | flip / pipeline σ |
|---|---|---|---|---|---|
| Draco (486 images) | 1.0-1.5 | 7.5e-12 | 9.5e-12 | 2.6e-11 | 11 |
| | 1.5-2.0 | 6.3e-12 | 9.1e-12 | 2.1e-11 | 11 |
| | 2.0-2.5 | 5.3e-12 | 7.8e-12 | 1.8e-11 | 15 |
| | 2.5-3.3 | 3.5e-12 | 5.5e-12 | 1.7e-11 | 21 |
| Ursa Minor (475 used) | 1.0-1.5 | 9.7e-12 | 1.6e-11 | 4.6e-11 | 20 |
| | 1.5-2.0 | 2.9e-11 | 3.5e-11 | 5.1e-11 | 2.4 (D2 has only 54 UMi images) |
| | 2.0-2.5 | 9.6e-12 | 1.2e-11 | 3.2e-11 | 13 |

The photon-noise forecast on the real coverage is ~5-6e-12 at 1.5-2.5 eV for Draco alone. That is close to the
proposal's 4e-12, but the data are 3.3-3.9x worse in g, 11-15x in line amplitude. Two things set the floor:
- QR2 gives only ~20-25 exposures with a given wavelength on a dwarf.
- Exposure-coherent residuals remain after the pattern.

The measured numbers happen to sit near 2412.12286's dSph forecast (2.1e-11 at 2 eV).

**Starlight.**
- The ON-OFF continuum of Draco's central 0.1° is 21, 10 and 5 nW m^-2 sr^-1 at 1.0, 1.6 and 2.2 µm, after the
  pipeline's source mask. At 0.1-0.2° it is ≤ 1 nW. Ursa Minor's centre is 8, 5 and 0.8 nW at the same wavelengths.
- The narrow stellar-feature amplitude β is consistent with zero at every mass. β/σ has an rms of 0.6-1.1 for Draco
  and 0.5-1.6 for Ursa Minor, with |β| typically 3-7% of the local continuum.
- Adding the nuisance moves Draco's expected limit by ≤ 8%.
- So the stellar residual is below the noise. It is not what limits the reach.

**Joint expected limit** (`alp_gate.py`; profile likelihood, D-factor log-error profiled, Asimov A = 0):
- Draco + UMi use their measured sign-flip σ.
- The six other dwarfs get a flattering bound: Draco's σ scaled by D(0.5°) and √(exposures), with no zodi or
  coverage penalty, although Sextans and the Leos sit near the ecliptic.

| m_a (eV) | Draco + UMi | 8 dwarfs (flattering) | 2609.10674 final (weaker JWST) | 2503.11753 (strong) | 8-dwarf / weaker |
|---|---|---|---|---|---|
| 1.00 | 2.8e-11 | 2.6e-11 | 2.3e-11 | 6.7e-13 | 1.12 |
| 1.25 | 2.8e-11 | 2.5e-11 | 2.2e-11 | 6.2e-13 | 1.10 |
| 1.50 | 2.3e-11 | 2.1e-11 | 2.2e-11 | 5.6e-13 | 0.95 |
| 1.75 | 2.3e-11 | 2.1e-11 | 1.8e-11 | 6.2e-13 | 1.18 |
| 2.00 | 1.9e-11 | 1.7e-11 | 1.3e-11 | 6.7e-13 | 1.32 |
| 2.25 | 2.4e-11 | 2.1e-11 | 1.3e-11 | 5.5e-13 | 1.72 |
| 2.50 | 1.9e-11 | 1.7e-11 | 1.2e-11 | 4.2e-13 | 1.41 |

**Gate: fails.** Over 1.5-2.5 eV the expected limit is 0.95-1.7x the weaker disputed JWST limit, equal only at
1.5 eV, and 30-50x above the strong one. A null would not test 2503.11753's claim. Its interval would also not
improve on 2609.10674's limit, except marginally at 1.5 eV. The only possible outcome is an independent limit of
~2e-11 that agrees with the weaker side, which arbitrates nothing.

**What would change it.**
- The line amplitude's error falls as N_exp^-1/2, so g falls as N_exp^-1/4. Reaching 1.2e-11 at 2-2.5 eV takes ~4x
  QR2's exposures on Draco, beyond the 2-year mission (~2x).
- Reaching 2503.11753's level would take ~10^6 x. No SPHEREx dataset can do it.
- A deeper exposure-level calibration could close part of the 11-15x gap to photon noise. The candidates are a
  time-resolved or per-survey-week pattern, and stray light and persistence models. The residual's origin is not
  identified: it is coherent within an exposure, varies on 1-6 ln-λ bins (10-60 rows), and survives a plane fit.
  Even photon-limited data (~5e-12, Draco + UMi ~4e-12) cannot test 2503.11753. They would only arbitrate the
  middle ground (5e-12 to 1.2e-11) that neither side's limit covers.

## Loose ends

- The fake-dwarf image lists (`alp_list.py fake`; four blank fields at matched ecliptic latitude, ~7° away) were
  started and stopped at the gate: only fk_draco was listed (673 images).
- Sculptor, Fornax, Sextans, Carina, Leo I and Leo II were not reduced. Their contribution above is the flattering
  scaling.
- The pattern uses every third image and a per-image scale of one template. A per-pixel offset + zodi model (the
  `lin` option in `alp_reduce.load_pattern`) was written but not used.
