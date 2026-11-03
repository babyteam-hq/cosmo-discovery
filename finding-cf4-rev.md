# Supernova Hubble residuals against reconstructed velocity fields, with a calibrated null: the best-validated field (2M++) fails its own null, and so does a correct filter of the Cosmicflows-4 data, so the test cannot grade the Cosmicflows-4 field

2026-11-03 (R-168). Scripts: `cf4rev_common.py`, `cf4rev_chain.py`, `cf4rev_fit.py`, `cf4rev_realwf.py`, `cf4rev_table.py`;
plan fixed before any data fit in `plan-cf4-rev.md`; rows in `results-cf4-rev.md`.

## Claim

Type Ia supernova Hubble residuals at 0.01 < z < 0.1 (ZTF SN Ia DR2 at z < 0.06, 551 SNe, the sample independent of
Cosmicflows-4; Pantheon+, 618 SNe) are regressed, by generalised least squares, on a reconstructed velocity field's
prediction at each SN, split into the field's bulk-flow dipole at the SN positions and the non-dipole rest, each with a free
scale β (joint fit). A real reconstruction is a shrunk posterior mean, so β = 1 is not what a correct field returns. We
calibrate the null with mock chains built on the real data's geometry: a ΛCDM truth (256³ in 1000 Mpc/h, sub-grid power,
redshift-space positions, a 250 km/s small-scale velocity shared in part within CF4 groups); a CF4-like reconstruction (the
Wiener-filter posterior mean from the 55,694 real Cosmicflows-4 galaxies observed with their own catalogue errors, the
Pantheon+ hosts carrying the SNe's own mock noise because they are CF4 inputs); and a 2M++-like reconstruction (log-normal
galaxies from a 2MASS K-band n(r), redshift space, zone-of-avoidance cloning, 4 Mpc/h smoothing, iterative real-space
reconstruction, scale and external dipole fitted to mock velocities and to the Pantheon+ SNe below z = 0.067, as Carrick et
al. 2015 did). Residuals are mocked with each sample's own error model (ZTF: host-mass step, magnitude selection, σ_int refit
per mock; Pantheon+: its covariance with the 2M++ peculiar-velocity systematic removed).

The pre-registered gate is the 2M++ check: 2M++'s observed β must lie within the central 99 % of its own 2M++-like null on
ZTF z < 0.06 and on Pantheon+ (Gaussian fit to 300 realisations, checked against the empirical rank), or the estimator is
not trusted. **It fails on both samples**, on the bulk-flow part:

| sample | 2M++ part | observed β | 2M++-like null (300) | z | two-sided p (Gaussian / rank) | verdict |
|---|---|---|---|---|---|---|
| ZTF z<0.06 | bulk flow | 3.95 ± 0.85 | 0.47 ± 0.74 | +4.7 | < 10⁻⁵ / 0 of 300 | fail |
| ZTF z<0.06 | non-dipole | 0.31 ± 0.19 | 0.96 ± 0.35 | −1.9 | 0.063 / 0.027 | pass |
| Pantheon+ | bulk flow | 2.49 ± 0.35 | 0.67 ± 0.67 | +2.7 | 0.007 / 0.007 | fail |
| Pantheon+ | non-dipole | 0.43 ± 0.13 | 1.03 ± 0.31 | −1.9 | 0.055 / 0.007 | pass (Gaussian); the rank test would fail it |

The failure does not depend on the 2M++-like mock choices: with a 2LPT particle truth (tracers, hosts and SNe drawn from the
particles near their real positions, 40 realisations) the 2M++ bulk-flow part is +3.6 null-sd (ZTF z < 0.06) and +3.5
(Pantheon+) and its non-dipole part −3.0 and −2.4; with bias 1.0, +3.1 and +3.9 (non-dipole −2.4, −3.2); with bias 1.4 the
null has outliers (sd 2.3), and on its median and robust sd the bulk-flow parts are +6.3 and +4.1.

The predictive statistic says the same of the non-dipole part: the whitened correlation between the residuals and 2M++'s
non-dipole prediction is 0.07 (ZTF z < 0.06) and 0.14 (Pantheon+), against 0.21 ± 0.05 and 0.31 ± 0.06 for a correct
2M++-like field (−2.7 and −3.2 null-sd).

So the analysis stops where the plan says it stops: with the null calibrated, the estimator still does not return for the
best-validated field what a correct reconstruction of that kind would, and no β can be read as a deficiency of Cosmicflows-4.

The rest of the table, recorded but not interpreted (CF4-like null for CF4++ and for WF(real CF4), the chain's own Wiener
filter applied to the real Cosmicflows-4 catalogue):

| sample | field | β_rest observed | null | ratio | z |
|---|---|---|---|---|---|
| ZTF z<0.06 | CF4++ | 0.19 ± 0.08 | 0.98 ± 0.21 | 0.19 ± 0.22 | −3.8 |
| ZTF z<0.06 | WF(real CF4) | 0.23 ± 0.23 | 0.98 ± 0.21 | 0.23 ± 0.22 | −3.6 |
| ZTF z<0.06 | 2M++ | 0.31 ± 0.19 | 0.96 ± 0.35 | 0.33 ± 0.36 | −1.9 |
| Pantheon+ | CF4++ | 0.23 ± 0.06 | 1.38 ± 0.14 | 0.17 ± 0.10 | −8.4 |
| Pantheon+ | WF(real CF4) | 0.72 ± 0.14 | 1.38 ± 0.14 | 0.52 ± 0.10 | −4.8 |
| Pantheon+ | 2M++ | 0.43 ± 0.13 | 1.03 ± 0.31 | 0.42 ± 0.30 | −1.9 |
| ZTF z<0.1 | CF4++ | 0.27 ± 0.06 | 0.96 ± 0.42 | 0.28 ± 0.44 | −1.6 |
| ZTF z<0.1 | WF(real CF4) | −0.19 ± 0.14 | 0.96 ± 0.42 | −0.19 ± 0.44 | −2.7 |

(β_dipole for CF4++: 0.40 ± 0.12 on ZTF z < 0.06 and 0.89 ± 0.13 on Pantheon+, against nulls of 1.03 ± 0.82 and
0.96 ± 0.41: the dipole part constrains nothing.)

Every field falls short of its calibrated null on its non-dipole part: CF4++, the 2M++ field, and the null's own linear
Wiener filter applied to the real Cosmicflows-4 catalogue. A shortfall shared by all three reconstructions, including
one built by the null's own method from the real data, is a property of the data or the residuals relative to the mocks
(their ΛCDM linear truth, their error models), not of CF4++. The paired CF4++ − 2M++ difference on the same SNe
(bootstrap over the whitened rows) is −0.16 ± 0.23 (ZTF z < 0.06) and −0.21 ± 0.14 (Pantheon+) for the non-dipole part: no
difference between the fields (null −0.6 to +0.6 and −0.3 to +0.8 at 95 %). The bulk-flow parts do differ (Δβ_dipole −3.6 ± 0.7
and −1.6 ± 0.3, outside their nulls), but that is the 2M++ dipole failing R1, not a CF4 result. The direction offset the
earlier note reported without an error: the SNe's own bulk flow (free vector, GLS) is 453 ± 79 km/s toward (l, b) =
(310, 5) ± 11° on ZTF z < 0.06 and 281 ± 41 km/s toward (312, 20) ± 12° on Pantheon+; CF4++'s vector at the same SNe points
68° (bootstrap 61-74°) and 38° (30-46°) away, but a correct CF4-like field at those SNe is that far from the SN fit in a
quarter of realisations (null median 39° and 23°; ranks 0.76, 0.75). The offset is not significant.

## What it means

The earlier version of this note (first section above) read β_rest = 0.10 to 0.26 for CF4++ against a null of 0.88 that
treated the field as the truth. With a null built from CF4-like and 2M++-like reconstructions, a correct field returns
β_rest ≈ 1.0 on ZTF and 1.4 on Pantheon+ (above 1 because the reconstruction carries each Pantheon+ SN's own noise), and
CF4++ falls 4 to 8 null-sd below; but 2M++ falls 1.9 null-sd below its own null (with a rank p of 0.007 on Pantheon+), its
bulk-flow part lies 4.7 null-sd above it on ZTF and 2.7 on Pantheon+ (the SNe's own bulk flow, 450 ± 80 km/s toward (l, b) = (310, 5) at z < 0.06, is several times what a
2M++-like reconstruction of a ΛCDM sky carries there), and the null's own Wiener filter run on the real Cosmicflows-4
catalogue falls 3.6 to 4.8 null-sd below. The mocks therefore do not describe the real residuals well enough to grade any field.
The candidates: the real flow is larger and less linear than a ΛCDM Gaussian truth on these scales (the bulk-flow excess), or
the SN residuals carry sky-correlated structure the error models omit (the ZTF field jackknife gives β_rest errors about twice
the GLS ones). Either way, the β shortfall of CF4++ cannot be attributed to CF4++.

Two field-level facts stand on their own, without the estimator. CF4++ carries 2 to 4 times the radial velocity amplitude
that the linear Wiener filter extracts from the same catalogue (rms 463 vs 204 km/s at the Pantheon+ SNe, 693 vs 239 at ZTF),
and beyond 150 Mpc/h the two are anti-correlated (−0.34, −0.41 in shells); and the real catalogue's scatter about the
filter's prediction is smaller than its stated errors imply (χ²/datum 0.87 against the mocks' 1.00 ± 0.003), so the
catalogue errors are not underestimated. These say CF4++ differs from a linear filter of its data, not that it is wrong.

## Weakest points

- The mock truth is a Gaussian (and, as a variant, 2LPT) ΛCDM field; if the real Universe's flow on 100-200 Mpc/h scales is
  larger than ΛCDM allows, no ΛCDM null calibrates a field fit to it, and the 2M++ dipole failure may be that physics rather
  than an estimator fault. The rule treats both the same.
- The 2M++-like chain is number-weighted (2M++ is luminosity-weighted) and its n(r) is a Schechter model, not the catalogue's;
  the fitted scale (0.4 to 1.3 across realisations) shows the shot-noise dilution Carrick et al.'s fitted β* also absorbs.
  The R1 failure survives bias 1.0 and 1.4 and a 2LPT truth.
- CF4++'s own HMC (prior, Malmquist treatment, zero-point fit) is not emulated; the null is a correct linear posterior mean.
- Pantheon+'s m_b,corr comes from a BBC fit on the 2M++-corrected z_HD; the 2M++ systematic is removed from the covariance
  by differencing grouping files, an approximation (38 % of SYS at low z).
- The CF4-only mock variants (prior, smoothing, grouping, error model, zero-point step, v_nl, σ8, clustered positions, the
  Watkins-Feldman estimator null for WF(real CF4)) were built but not run once R1 failed: the rule stops before they matter.
  The 2M++-chain variants that bear on R1 were run (below).

## To settle it

A null the real 2M++ field passes: a mock truth whose bulk flow on 100-200 Mpc/h matches what the SNe see (constrained
realisations conditioned on the SN or CF4 bulk flow), or an N-body light-cone with 2M++ and CF4 selection; and SN residual
mocks with the sky-correlated scatter the ZTF field jackknife shows. Until the best-validated field passes, the test grades
nothing.

## Blind review (2026-11-03, finding check, this file alone)

VERDICT: novel but not sound: the gate that "fails" uses a 2M++-like null that does not copy how the real Carrick field and its external dipole were built, and it cannot separate an estimator fault from a real bulk flow larger than ΛCDM allows, so the stated diagnosis does not follow.

Strongest objection: the real 2M++ field is luminosity-weighted with β* and V_ext fitted to SFI++ and the First Amendment SNe, while the 2M++-like mocks are number-weighted and fit V_ext to mock Pantheon+ SNe; the bulk-flow part that fails is mostly that fitted V_ext, and the SNe's own 450 km/s flow restates the known bulk-flow excess over ΛCDM, so the gate mixes estimator and cosmology. Also named: dipole leakage into β_rest in the joint fit (a free SN dipole would test it), Hollinger & Hudson 2024 (arXiv:2312.03904) as prior N-body calibration of a 2M++-like reconstruction, the ZTF field jackknife errors (about twice the GLS ones) that would shrink every z, and the Pantheon+ BBC-on-z_HD leak.
