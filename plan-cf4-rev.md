# Analysis plan: does the Cosmicflows-4 velocity field predict supernova Hubble residuals as well as a correct CF4-like reconstruction would?

Fixed before any calibrated amplitude is computed.

## Question

Type Ia supernova Hubble residuals at 0.01 < z < 0.1 carry the radial peculiar velocity of each host. A reconstructed velocity
field predicts that velocity. Fit the residuals against the field's prediction with a free scale β (β = 1 if the field is
exactly the true velocity). A real reconstruction is not the truth: it is a shrunk posterior mean from noisy, incomplete data,
and the scale a *correct* reconstruction returns through this estimator is not known a priori. The question is whether the
Cosmicflows-4 fields (the CF4++ posterior-mean grid of Courtois et al. 2025, arXiv:2502.01308) return the β a correct
reconstruction with Cosmicflows-4's data quality would, and the same for the 2M++ field (Carrick et al. 2015, arXiv:1504.04627)
against a correct 2M++-like reconstruction, which is the check that the estimator can be trusted at all.

## Samples

1. **ZTF SN Ia DR2** (Rigault et al. 2025, public "lite" release): 1,122 SNe Ia at 0.015 < z < 0.1 (and the 551 at z < 0.06),
   host spectroscopic redshifts, SALT2 x0, x1, c; cuts as released for cosmology (snia-cosmo, coverage and fit-quality flags,
   |x1| < 3, |c| < 0.3). Not used by Cosmicflows-4 (whose SN distances predate DR2). Residual y = −2.5 log x0 − μ(z_CMB), with
   α and β_SALT free (nuisance columns), diagonal errors (SALT covariance, 300 km/s velocity noise, redshift error, σ_int set
   for χ²/dof = 1) and a host-mass step nuisance (global host mass > 10, DR2 host table). 113 ZTF hosts are CF4 galaxies (TF/FP
   distances, matched within 30" and 1,000 km/s): their small-scale velocity is shared with the CF4 datum in the mocks, and a
   variant drops them. DR2 photometry is not cosmology-grade (the release says so); only relative magnitudes over the sky matter here.
   The only sample independent of Cosmicflows-4. **The headline is ZTF at z < 0.06** (551 SNe; see the power forecast below:
   the z < 0.1 sample's null is too broad to separate anything at 2σ).
2. **Pantheon+** (Scolnic et al. 2022): 618 SNe at 0.01 < z_CMB < 0.1, calibrators removed. Residual y = m_b,corr − μ(z_CMB)
   (no 2M++ correction in y; but m_b,corr comes from the release's SALT2mu/BBC fit run on the 2M++-corrected z_HD, so a little
   2M++ information can leak into y: for this reason the stop rule R1 rests on ZTF, and Pantheon+'s 2M++ check is reported, not
   used to stop). Covariance: STAT+SYS minus the VPEC systematic (the 2M++-based peculiar-velocity systematic, from
   the release's systematic-grouping files: G_VPEC minus the mean of seven groupings whose systematic is negligible; 38 % of the
   low-z SYS). The 250 km/s diagonal velocity-noise term in STAT is kept (field-agnostic). 602 of 618 are matched to a
   Cosmicflows-4 galaxy within 120" and 1,000 km/s (CF4 used their SN distances: the samples share inputs).
3. **Pantheon+ minus CF4 inputs**: empty. Every Pantheon+ SN at 0.01 < z < 0.1 is a CF4 SN-distance input (`cf4_notin.py`,
   VizieR J/ApJ/944/94; the 16 that fail the match above are hosts at larger separation, not independent objects). Not run;
   instead the shared inputs are modelled in the null (below).
4. **DES-SN5YR Dovekie** low-z (191 SNe; DES distances, DES covariance as released): a pipeline check only, not an independent
   sample (same objects as Pantheon+ largely; its covariance keeps DES's own VPEC term; mock noise drawn independently).

## Fields under test

- CF4++ posterior-mean velocity grid (128³, 1000 Mpc/h, supergalactic, CMB frame), read by trilinear interpolation at the SN's
  redshift-space position (comoving distance of z_CMB, Ω_m = 0.334).
- 2M++ (Carrick et al. 2015 public grid, 257³, ±200 Mpc/h, Galactic, β* = 0.43 with V_ext added), same reading; outside the
  ±200 Mpc/h cube the field is V_ext alone. The same field on every sample, ZTF included (axis order checked: correlation with
  Pantheon+'s own VPEC column 0.92 in the documented order vs 0.46 transposed).

**Shared small-scale velocity.** Every object's 250 km/s unresolved velocity is split into a 150 km/s part shared by all
members of its CF4 group (1PGC; SNe inherit their matched host's group) and the rest per object. SN pairs in the same group
get the shared term in their covariance (Pantheon+: added to the matrix; ZTF: a full matrix replaces the diagonal), in data
and mocks alike.

## Statistic

The field's radial velocity at the SN positions is split by least squares into a monopole a·r (dropped), a dipole B·r̂ and the
rest. Each part is turned into a distance-modulus correction c = μ(z_HD) − μ(z_CMB), (1+z_HD) = (1+z_CMB)/(1+v/c). Joint GLS:
y = nuisance + β_dipole c_dipole + β_rest c_rest (both at once, GLS). Also the whole-field β, and a predictive statistic: the whitened correlation between the residual and
c_rest after the nuisance and dipole columns are projected out (ρ_rest), with its own null. β measures a scale; ρ measures
whether the field's non-dipole structure predicts the residuals at all.

## Null (route 2: a CF4-like mock chain built here; no public CF4 posterior realisations exist)

Per realisation (300 for the baseline, 60 for each variant):
1. Truth: ΛCDM Gaussian density on a 256³ grid (1000 Mpc/h, 3.9 Mpc/h cells, so the truth has power below the
   reconstruction's 7.8 Mpc/h cells), observer at the centre, Eisenstein-Hu no-wiggle P(k), Ω_m 0.31, σ8 0.81, linear
   velocities, f = Ω_m^0.55, plus 250 km/s unresolved velocity per object.
2. Redshift-space positions: every object (CF4 galaxy, SN) sits at its real observed position r_obs = cz/H0; its true
   position solves r_true + (v_r(r_true) + v_nl)/H0 = r_obs; the truth is read at r_true, the reconstructions at r_obs, as the
   data are. The real CF4 positions carry the survey's own selection by construction.
3. CF4-like catalogue: the 55,694 real CF4 galaxies with a distance error, each observed as a radial velocity with its own
   catalogue error ((ln10/5) e_DM V_cmb) plus the 250 km/s unresolved velocity. For the Pantheon+ hosts, the CF4 observation
   carries the same mock light-curve noise and small-scale velocity as the mock SN residual (shared inputs).
4. CF4-like field: the Wiener filter (posterior mean under the prior and the stated errors) on the CF4++ grid (128³).
5. 2M++-like field (Carrick et al.'s steps): galaxies Poisson-sampled from a log-normal field n(r) exp(bδ − b²σ²/2) (b = 1.2,
   2MASS K-band Schechter n(r), K < 11.5 or 12.5 by region, nothing beyond 200 Mpc/h), moved to redshift space; the zone of
   avoidance (|b| < 5°, < 10° within 30° of the Galactic centre) filled by cloning the mirrored strips; δ_g smoothed by a
   4 Mpc/h Gaussian; three iterations of the real-space reconstruction; a scale and an external dipole fitted to the mock CF4
   velocities within 100 Mpc/h **and the Pantheon+ SNe at z < 0.067 with their shared mock noise** (Carrick et al. fitted
   SFI++ and the A1 SNe, which are Pantheon+ objects). Not emulated: luminosity weighting (number weighting instead).
6. Mock residuals at each sample's real positions (Pantheon+: its covariance minus the velocity block, plus the shared
   250 km/s velocity; DES: its covariance; ZTF: its diagonal errors, then the ZTF magnitude selection (logistic in peak m_B,
   m50 = 19.0, width 0.15, the observed m_B ceiling at z > 0.07) applied by redrawing the noise of rejected SNe, and σ_int
   refitted per mock as for the data).

Variants (one change each), so the null's dependence on each mock choice is measured, not assumed:
- CF4 reconstruction: WF prior σ8 0.70 and 0.90 (truth 0.81); WF prior shape n_s = 0.90; the WF field smoothed by an 8 Mpc/h
  Gaussian; grouped input (CF4 groups by 1PGC, one datum per group, shared small-scale velocity); log-normal distance errors
  with no Malmquist correction; true errors 1.3× the stated ones the filter uses; an FP zero-point step (FP-only distances
  +0.05 mag north of the equator, −0.05 south: the 0.1 mag step measured in this repo, R-119); no shared inputs; the first
  chain (128³ truth, real-space positions).
- More CF4-chain variants (method check pass 2): unresolved velocity 150 and 400 km/s (in the data model and the filter); the
  truth's σ8 at 0.74 and 0.86 with the prior matched; "clustered" positions (each tracer, host and SN moved to a cell within
  ±7.8 Mpc/h of its real position with weight exp(1.2 δ), so they sit in the mock's overdensities while the survey selection is
  kept on larger scales; the reconstruction rebuilt on the moved positions).
- 2M++-like: bias 1.0 and 1.4.
- ZTF (fit stage): per-ZTF-field zero points of 0.01 and 0.02 mag injected; selection m50 = 18.8 and width 0.30; an E(B−V)
  nuisance column; the 113 CF4-host SNe dropped. The selection redraws the noise only, not the velocity term (a small effect
  at z < 0.1).
The null for a field is the distribution over realisations of β through the same estimator when the field is the CF4-like
(respectively 2M++-like) reconstruction of that realisation's truth. Also recorded: the truth injected as the field.

- Method check pass 3: a 2LPT truth ("lpt2": particles displaced and moved at second order from the same Gaussian field; every
  CF4 galaxy, host and SN is a particle drawn from the 5×5×5 cells around its real position, weighted by the particle count,
  so tracers sit in collapsing structure with 2LPT, not linear, velocities; the 2M++-like galaxies are the particles thinned
  to n(r), bias 1) through both chains; and "wfest" (log-normal distance moduli through the real-data Watkins-Feldman estimator
  and H0 rescale), which is the null for WF(real CF4). The ZTF selection's acceptance uses the full mock magnitude (true
  distance modulus, so the velocity term, plus noise); only the noise is redrawn because each position's velocity is fixed.

**The same Wiener filter on the real CF4 catalogue** (method check pass 2, point 1): the filter is applied to the real
CF4 velocities (Watkins-Feldman estimator u = V ln(V/H0 d), a global H0 rescale removed), giving a field WF(real CF4) whose
null is the "wfest" variant (same estimator, H0 rescale and errors on mock distance moduli), with the baseline alongside. It is put through the estimator on every sample as a third field, and compared
with CF4++ at field level (rms radial velocity in shells, correlation at the SN positions). The χ² per datum of the real data
about it is compared with the mocks' (stated errors, and errors 1.3× stated). This separates "the CF4 data, filtered
correctly, fail" from "CF4++'s algorithm departs from a correct filter".

What the chain does not reproduce, stated as the result's scope: CF4++'s own algorithm (an HMC with its own prior, Malmquist
treatment and zero-point fit); the null is the β of a correct linear Gaussian posterior mean with CF4's data and the listed
departures from it. N-body nonlinear structure (infall around clusters) beyond the sub-grid Gaussian power, the varied
dispersion and the clustered-positions variant. The 2M++-like chain is number-weighted, not luminosity-weighted, so R1's
"estimator validated" is conditional on that choice.

Power forecast for every sample (54 baseline realisations of the final chain, null only, before any data fit): CF4-like
β_rest null on ZTF z < 0.1 0.86 ± 0.42 (too broad: no observed value in [0, 0.86] could be 2σ away), on ZTF z < 0.06
1.01 ± 0.20, on Pantheon+ 1.39 ± 0.15, on DES 0.93 ± 0.33; ρ_rest null: ZTF z < 0.06 0.21 ± 0.05, Pantheon+ 0.45 ± 0.04 (high:
the reconstruction carries each SN's own noise). 2M++-like β_rest null: ZTF z < 0.06 0.95 ± 0.32, Pantheon+ 1.04 ± 0.23; the
2M++-like dipole nulls have sd ≈ 0.8-1.3 and constrain nothing. Hence the headline moved from ZTF z < 0.1 to ZTF z < 0.06,
the tightest sample independent of CF4, and R1 applies to it and to Pantheon+.

Feasibility forecast (40 realisations of the first chain, Pantheon+ only): CF4-like null β_rest 1.35 ± 0.17, β_dipole
0.98 ± 0.31; truth as the field 1.01 ± 0.11 / 1.03 ± 0.39 in this joint fit (the 0.88 quoted earlier for the truth was a
separate rest-only fit with the dipole left in the residual: a different statistic). The dipole cannot separate anything.

## What is reported

Per sample and field: β_obs (joint fit), the null mean and spread, the ratio β_obs / ⟨β_null⟩ with a 68 % interval from the
null spread (which holds the SN noise and cosmic variance), and the fraction of null realisations at or below β_obs.

**Decision rules (fixed now).**
- R1, the 2M++ check first: 2M++'s observed β_rest and β_dipole on ZTF z < 0.06 and on Pantheon+, each against the baseline
  2M++-like null (Gaussian fit to 300 realisations, checked against the empirical rank). If any is outside the central 99 %
  (two-sided p < 0.01), the estimator is not trusted on that sample and the analysis stops there; that is the result. (A ZTF-only
  rule, as pass 2 suggested, would be uninformative: the forecast's ZTF z < 0.1 2M++-like null is 0.82 ± 0.40. A Pantheon+
  failure can come from the estimator or from the BBC z_HD leak; either way it stops.) The dipole parts are weak gates (null sd
  ≈ 1) and are reported as such. A pass validates the estimator and the residuals only, not the
  CF4-like stand-in, and is conditional on the 2M++-like chain's number weighting.
- R2, the CF4++ claim, as a statement about scale: "the non-dipole part of CF4++ carries less of the SN-measured flow per unit
  of its own amplitude than a correct CF4-like reconstruction does (over-amplified or decorrelated)" only if, on ZTF (z < 0.06),
  β_rest,obs is below the null mean by more than 2 null-sd (one-sided p < 0.025) under **every** CF4 and ZTF variant (the worst
  case is the variant closest to the data). The predictive version ("decorrelated") additionally needs ρ_rest,obs below its
  null by more than 2 null-sd under every variant. Two-sided: β_rest,obs above the null mean by more than 2 null-sd under every
  variant reads "CF4++'s non-dipole part is shrunk more than a correct reconstruction". "Consistent with a correct
  reconstruction" if within 2 null-sd under the baseline. Anything else is "inconclusive". Pantheon+ is reported the same way,
  second. The same rules are applied to WF(real CF4), whose null is the same.
- R3, the paired difference: significant only if Δ_obs is outside the central 95 % of the baseline null Δ distribution and the
  bootstrap 95 % interval excludes the null mean.
- The ratio β_obs / ⟨β_null⟩ is quoted with the baseline null (68 % interval from the null sd) and the range of ratios across
  variants.

**Paired CF4++ − 2M++ difference** per sample, ZTF first: Δβ_rest = β_rest(CF4++) − β_rest(2M++) from separate fits on the
same SNe, its error from (a) a pairs bootstrap (1,000 draws) over the rows of the whitened system (≈ SNe), refitting both
fields on each draw, and (b) a joint fit with both fields' columns at once (the GLS covariance of the two coefficients);
and the null distribution of Δ from the mocks (both reconstructions of the same truth on the same SNe).

**Bulk-flow direction**: the free bulk-flow vector fitted to the residuals (3 parameters, GLS) and the CF4++ mean vector at
the same SNe; the angle between them, its error from the same bootstrap, and its null distribution in the mocks (SN-fitted
vector vs the CF4-like field's vector).

## Trials

2 fields × 3 samples (ZTF, ZTF z < 0.06, Pantheon+; DES as pipeline check) × 2 parts (dipole, rest). The headline is fixed
now: ZTF (z < 0.1) β_rest ratio for CF4++, with Pantheon+ second; anything else is supporting.

## What each reading predicts

- CF4++ is a correct CF4-like reconstruction: its ratio ≈ 1 in every sample (rest and dipole), as 2M++'s is.
- CF4++'s non-dipole part is miscalibrated (structure uncorrelated with the true flow, or amplified): its β_rest ratio ≪ 1 in
  ZTF and Pantheon+, while 2M++'s ratio ≈ 1; the paired Δ below its null.
- The estimator or the residuals are at fault: 2M++'s ratio far from 1 too (the stop condition).
- Shared inputs matter: the CF4-like null is higher on Pantheon+ than on ZTF at comparable noise.
