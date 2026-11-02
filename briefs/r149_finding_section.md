# Photometric-redshift isolation leaks unevenly by galaxy type, but not in lensing amplitude: on DESI DR1 BGS spectroscopic truth, KiDS-bright-style isolation leaves red galaxies with a real > 10 % neighbour within 0.3 Mpc 1.7 times as often as blue ones (32 % against 19 %); a ΛCDM mock that reproduces this leakage predicts that it shifts the early/late lensing-RAR offset by −0.001 dex (68 %: −0.005 to +0.002), at most 0.09 of the published 0.14-0.19 dex in any variant

## Claim

Brouwer et al. 2021 (arXiv:2106.11677) found isolated KiDS-1000 early- and late-type galaxies ≥ 6σ apart in the weak-lensing
radial acceleration relation (RAR). Their measure is log10⟨g_obs,E/g_obs,L⟩: 0.14 dex (Sérsic split) and 0.19 dex (u−r split) at
g_bar > 10⁻¹³ m s⁻² (R < 0.3 Mpc/h70), and 0.17 and 0.27 dex over all 15 g_bar bins. Their isolation uses photometric redshifts: no
neighbour with > 10 % of the lens stellar mass within 3 Mpc/h70 in 3D. They assumed photo-z errors degrade the isolation of both
types equally (their Sec. 5.3). Mistele et al. 2024 (arXiv:2310.15248) find one relation with 4 Mpc/h70 isolation and their own
early-type masses.

Two results, in this order.

**1. The leakage is type-dependent (measured).** The KiDS-bright isolation is emulated on the DESI DR1 BGS group catalogue: the
published KiDS-bright photo-z error distribution, 5 seeds, about 780,000 candidate lenses per seed. Each photo-z-isolated lens is
then checked against the spectroscopic redshifts. Types are mass-matched as Brouwer did. Interval: sky jackknife (50 regions) plus
seed scatter, statistical only.

| split | fraction with a real > 10 % neighbour (500 km/s) | red / early | blue / late | difference (95 %) | ratio (95 %) |
|---|---|---|---|---|---|
| rest-frame g−r | within 0.3 Mpc | 0.319 | 0.192 | 0.128 (0.124-0.132) | 1.67 (1.64-1.69) |
| rest-frame g−r | within 0.5 Mpc | 0.453 | 0.303 | 0.149 (0.145-0.153) | 1.49 (1.48-1.51) |
| rest-frame g−r | within 1 Mpc | 0.651 | 0.517 | 0.134 (0.129-0.139) | 1.26 (1.25-1.27) |
| rest-frame g−r | lens is a hidden satellite of a group with such a neighbour | 0.290 | 0.113 | 0.177 (0.174-0.181) | 2.57 (2.53-2.62) |
| spectral QUIESCENT | within 0.3 / 0.5 Mpc | 0.327 / 0.461 | 0.200 / 0.314 | 0.127 / 0.147 | 1.63 / 1.47 |
| spectral QUIESCENT | hidden satellite | 0.310 | 0.118 | 0.192 (0.188-0.196) | 2.62 (2.57-2.68) |

Brouwer's equal-contamination assumption fails in fraction. The gap is set by red lenses being 2.6 times as often unrecognised
satellites. A 4 Mpc criterion lowers the fractions by ≤ 0.04. A Legacy Survey morphology split (DEV against EXP/REX) narrows the
gap (0.27 / 0.39 against 0.20 / 0.32 at 0.3 / 0.5 Mpc, unmatched).

**2. Its lensing effect is not type-dependent (calibrated on a mock).** The same isolation is run on the UniverseMachine DR1
catalogue of the Bolshoi-Planck simulation at z = 0.21, with a BGS- or KiDS-bright-like magnitude limit. That catalogue has
halos, stellar masses and star-formation rates, and its quenched fractions are fit to SDSS by halo and central/satellite status.
The mock reproduces the measured leakage (quenched vs star-forming, BGS-like): 0.37 / 0.48 against 0.22 / 0.32 within 0.3 /
0.5 Mpc, hidden satellites 0.30 against 0.14 (DESI: 0.33 / 0.46 against 0.20 / 0.31; 0.31 against 0.12).

The excess surface density (ESD) is computed from the mock's own halos, with no assumed amplitude. Per type, A_T is the stacked
ESD of photo-z-isolated lenses over that of truly isolated lenses of the same type. Both samples are binned in Brouwer's g_bar
with photometric masses and reweighted to one stellar-mass distribution. The leakage's contribution to the published statistic is
then δ = log10 mean_k (A_E,k / A_L,k). That is Brouwer's log⟨g_E/g_L⟩ with photo-z isolation minus the same with true isolation.
It is evaluated with the KiDS-bright-like limit (r < 20, z_v = 0.15 / 0.25 / 0.35).

| amplitude | criterion | δ at g_bar > 10⁻¹³ (68 %; 95 %) | of 0.19 / 0.14 | δ over all bins (68 %; 95 %) | of 0.27 / 0.17 |
|---|---|---|---|---|---|
| mock | 3 Mpc | −0.001 (−0.005 to 0.002; −0.008 to 0.006) | −0.01 / −0.01 | −0.008 (−0.020 to 0.003; −0.031 to 0.014) | −0.03 / −0.05 |
| rescaled to DESI's fractions (0.5 Mpc) | 3 Mpc | −0.002 (−0.005 to 0.001; −0.008 to 0.005) | −0.01 / −0.01 | −0.010 (−0.021 to 0.001; −0.032 to 0.011) | −0.04 / −0.06 |
| rescaled to DESI's hidden-satellite fractions | 3 Mpc | 0.004 (0.001 to 0.007; −0.003 to 0.011) | 0.02 / 0.03 | 0.010 (−0.001 to 0.021; −0.013 to 0.032) | 0.04 / 0.06 |
| mock | 4 Mpc | 0.006 (−0.002 to 0.013; −0.009 to 0.020) | 0.03 / 0.04 | 0.000 (−0.014 to 0.015) | 0.00 / 0.00 |
| rescaled, hidden satellites | 4 Mpc | 0.013 (0.005 to 0.021; −0.003 to 0.029) | 0.07 / 0.09 | 0.018 (0.002 to 0.033) | 0.07 / 0.11 |

- **Radius dependence (3 Mpc).** log10(A_E/A_L) stays within ±0.005 dex from 0.04 to 0.7 Mpc (errors 0.002-0.02). It turns
  negative beyond 1 Mpc: −0.019, −0.028 and −0.076 (± 0.04-0.06) at 1.3, 1.7 and 2.2 Mpc. Each type's own amplitude rises from
  1.01 at 0.04 Mpc to 1.12 at 0.23 Mpc, 1.30 at 0.54 Mpc and 2.1-2.5 at 2.2 Mpc, so the leakage is large, but it is nearly the
  same fraction for both types.
- **Why it cancels.** Red lenses collect more contaminating mass, but at fixed stellar mass their own (truly isolated) ESD is
  larger in about the same proportion: in the mock, log10(ESD_E/ESD_L) of truly isolated lenses is 0.09-0.21 at g_bar > 10⁻¹³ and
  0.25-0.42 below.

The claim is the two tables. Type-dependent photo-z leakage is real and measured, but in a ΛCDM mock that reproduces it, it
changes the published early/late offset by at most 0.07 (u−r) or 0.09 (Sérsic) of its value at g_bar > 10⁻¹³ in any
variant (95 % upper ends ≤ 0.15 and ≤ 0.21). Over all bins it changes it by −0.06 to +0.11. So the published split is not a leakage artefact. This does not say what the split is.

## Data and method

- **Data.** DESI DR1 BGS group VAC (bgs-groups-dr1-v1.0; r < 19.5, 2.23 M galaxies): fastspecfit stellar masses, rest-frame
  g−r, QUIESCENT, group membership and central/satellite. Lenses are spectroscopic, 0.05 < z < 0.3. In the photo-z view Brouwer's
  cuts apply (0.1 < z_ph < 0.5, 8.5 < log M*_ph < 11). All galaxies with a redshift are neighbours. H0 = 70, Ωm = 0.3, physical
  distances.
- **Photo-z emulation (data and mock alike).** Bilicki et al. 2021's (arXiv:2101.06010) fit to KiDS-bright minus GAMA, a
  generalized Lorentzian (a = 2.613, s = 0.0149) whose wings are the outliers, drawn per galaxy. Photo masses rescale by D_L² and
  get 0.12 dex scatter, as in Brouwer's MICE emulation. Isolation is in 3D, from photo-z positions, at 3 and 4 Mpc and 10 %.
- **Truth.** For the contamination statistic: a real neighbour with M* > 0.1 M*_lens within projected R and |Δv| < 500 km/s. For
  the mock reference: the same 3D criterion with true masses in true redshift space (Brouwer's App. A reference).
- **Mock.** UniverseMachine DR1 (Behroozi et al. 2019, arXiv:1806.07893) on Bolshoi-Planck, 250 Mpc/h box, snapshot z = 0.208.
  Stellar mass is UM's observed-mass column, and the type is quenched if sSFR < 10⁻¹¹ yr⁻¹.
  - Apparent magnitudes come from BGS itself: at each virtual redshift, m_r − DM is fit as linear in log M* per spectral type,
    with Gaussian scatter, by maximum likelihood with the r = 19.5 truncation.
  - Each box axis serves as the line of sight, observed at virtual redshifts z_v (BGS-like: 0.10 / 0.20 / 0.28 weighted by the
    DESI lens n(z); KiDS-like: 0.15 / 0.25 / 0.35, weights 1/4, 1/2, 1/4), with no evolution between them.
  - Not a light cone. MICE (Brouwer's mock) was not reachable without an account.
- **ESD from the mock.** Every distinct halo above 10^11 M☉ within projected 8 Mpc and ±30 Mpc of the lens enters, as an NFW
  truncated at its catalogue virial radius (Dutton & Macciò 2014 concentration), at its true projected offset. A satellite lens
  also gets its own subhalo (current mass and radius), and every lens gets its baryons as a point mass. The lens's own halo is
  therefore whatever the mock gives it, and the reference is per type.
- **Rescaling.** Each type's excess (A_T − 1) is multiplied by f_DESI / f_mock for that type, at 0.5 Mpc, 0.3 Mpc or for hidden
  satellites, measured with the BGS-like limit. The factors are 0.95 / 0.98, 0.89 / 0.91 and 1.05 / 0.84 (quenched / star-forming).
- **Intervals.** Mock: 16 sky-plane jackknife regions per realisation (3 lines of sight × 2 photo-z seeds); variance = mean
  jackknife variance + scatter over realisations. DESI: sky jackknife + seed scatter.
- **Fixed before any number was computed** (results-R147.md, timestamped 04:43 UTC): mock, magnitude model, virtual redshifts,
  photo-z and true isolation, ESD construction, amplitude and offset statistics, rescaling rule, intervals, variants, and what
  each outcome would mean. One implementation choice was made after: the four ESD samples are reweighted to one common
  stellar-mass shape, so that photo-z and truly isolated stacks of a type compare like with like.

## Checks

| check (KiDS-like unless stated) | δ at g_bar > 10⁻¹³ / all bins |
|---|---|
| main, 3 Mpc | −0.001 / −0.008 |
| BGS-like selection (r < 19.5, DESI n(z)) | 0.001 / −0.004 |
| untruncated NFW (z_v = 0.20, one realisation; baseline 0.003 / −0.011) | 0.003 / −0.005 |
| halos down to 10^10.5 M☉ (same) | 0.003 / −0.007 |
| common reference instead of per-type (the mock's excess over the truly isolated star-forming profile) | 0.010 / 0.082 |
| same, 4 Mpc | 0.023 / 0.090 |

The last two rows are not the claim. They show what a type-independent reference does: it gives back a spurious 0.01-0.02 dex at
small radius and 0.08-0.09 dex over all bins. That is the size of the error a common reference would make, and it is a quarter
of what an uncalibrated halo-model amplitude gave at small radius.

**Amplitude benchmark.** Over all lenses, with no type split, the mock's photo-z-isolated ESD is 1.22 times the truly isolated
ESD at 0.3 Mpc, 1.55 at 0.75 Mpc and 3.0 at 3 Mpc (mean 1.81 ± 0.07 over 0.3-3 Mpc). Brouwer's MICE test reports about 1.3
beyond 0.3 Mpc/h70. The mock therefore has more leakage mass at large radius than MICE, or Brouwer's figure refers to a narrower
range. Since the claim is a ratio of two types' amplitudes that agree at every radius, a common rescaling of the excess leaves
it near zero.

## What it means

- At fixed stellar mass, photometric isolation leaves red, quiescent galaxies with a real massive neighbour 1.5-1.7 times as
  often as blue ones, mostly because they are more often unrecognised satellites. Brouwer's equal-contamination assumption is
  wrong in fraction.
- It is nearly harmless for the early/late comparison. In ΛCDM as UniverseMachine populates it, the fractional ESD a type gains
  from leakage is the same for both types to ≤ 0.005 dex inside 0.7 Mpc. Red lenses' larger contamination is matched by their
  heavier own halos.
- The published 0.14-0.19 dex at g_bar > 10⁻¹³ is therefore not a leakage artefact, at the ≤ 0.09 level (≤ 0.21 at 95 %). The
  same holds for the full-range 0.17-0.27 dex at the ≤ 0.11 level, except beyond 1 Mpc, where the mock's leakage slightly
  favours late types.
- An observation outside the claim: Brouwer's statistic applied to the mock's truly isolated lenses gives 0.14 dex
  (g_bar > 10⁻¹³) and 0.26 dex (all bins) between quenched and star-forming. This is ΛCDM's halo-mass gap at fixed stellar
  mass, as the mock sets it. It is close to the published values, but it was not pre-registered and is not tested here.

## Weakest points

- **One mock, one galaxy model.** The cancellation rests on UniverseMachine's halo-mass gap between quenched and star-forming
  galaxies at fixed stellar mass, and on its satellite quenching. The mock matches DESI's leakage fractions to 0.02-0.04, which
  tests the satellite side but not the own-halo side.
  - If real isolated early and late types had the same halos (the strict "one RAR" null), the relevant number would be the
    common-reference row: 0.010-0.023 dex, 0.06-0.17 of the published value at small radius.
  - So the bound is ≤ 0.09 in the mock's universe and ≤ 0.17 under the strict null. Neither reaches the published offset.
- **Not a light cone.** A single z = 0.21 snapshot is observed at virtual redshifts; projection depth is ±30 Mpc of halos, with
  no matter outside halos above 10^11 M☉. The box (305 physical Mpc) wraps the photo-z line-of-sight error (~65 Mpc) about
  8 % of the time.
- **The benchmark is not reproduced.** The mock's all-lens amplitude at 0.3-3 Mpc is 1.8, against MICE's 1.3.
- **Sample stand-ins.** sSFR quenching in the mock and rest-frame g−r or QUIESCENT in DESI stand in for observer-frame u−r. DEV
  against EXP/REX stands in for Sérsic n ≷ 2, and gets no mock counterpart. The DESI interval is statistical only: the
  emulation's systematics (BGS depth and fibre incompleteness in dense regions, the error model) move the fractions by ≤ 0.03 in
  the variants run.

## To settle it

- The same pipeline on a second galaxy model with a different halo-mass gap by colour (e.g. a semi-analytic model, or MICE
  itself). The cancellation should then track that gap.
- A KiDS-bright × GAMA (or DESI DR2 BGS in KiDS-N) spectroscopic subsample, measuring the isolated lensing RAR by type with both
  isolation selections, would measure δ directly.
