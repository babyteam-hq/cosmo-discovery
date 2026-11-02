"""R-160 Step 0: feasibility. Pair counts and the null scatter of the pairwise amplitude C, from positions and redshifts
with *random* angles only (no real angle enters an alignment statistic here), plus the SPICE-RACS -> JVAS angle
transfer factor rho = <cos 2(chi0_SPICE - PA_JVAS)> on the same sources (a calibration, disclosed).

Power: the published amplitude dC_JVAS (r160_jvas.py: 0.0146 +- 0.0049 at n_v = 10, ~3.3 deg) carried to SPICE-RACS
angles is dC_JVAS * rho^2 (each source's angle correlates with its JVAS angle by rho). Reading 2D: pairs within
3.3 deg. Reading 3D: pairs within R comoving; amplitude the same dC_JVAS * rho^2 (the smallest 3D reading) or
dC_JVAS * rho^2 / f(R), f = the fraction of 3.3-deg pairs of the same sources that lie within R in 3D (the reading
in which the JVAS 2D excess is made by 3D-close pairs alone).

Usage: python -I r160_feas.py data/r160
"""
import sys

import numpy as np
from astropy.table import Table

import r160_common as c

DC_JVAS, E_DC_JVAS, THETA_J = 0.0146, 0.0049, 3.3
SCALES = (150.0, 300.0, 500.0)
MIN_ANG = 1 / 60.0  # deg: drop pairs closer than 1' (same object, lensed images)


def null_sigma(x, i, j, nreal=300, seed=1):
    pg = c.PairGeom(x, i, j)
    rng = np.random.default_rng(seed)
    vals = pg.C_many(rng.uniform(0, np.pi, (nreal, len(x))))
    return vals.mean(), vals.std()


FP_BINS = np.array([0, 0.02, 0.05, 1.0])


def rho_by_fp(j):
    """rho per fractional-polarisation bin (SPICE fracpol), from the JVAS same-source set."""
    fp = np.asarray(j["fracpol"], float)
    d = np.radians(np.asarray(j["derot_polangle"], float) - np.asarray(j["jvas_PA"], float))
    out = []
    for lo, hi in zip(FP_BINS[:-1], FP_BINS[1:]):
        m = (fp >= lo) & (fp < hi)
        out.append(np.mean(np.cos(2 * d[m])))
        print(f"  rho(fracpol {lo:.2f}-{hi:.2f}) = {out[-1]:.3f} (N = {m.sum()})")
    return np.array(out)


def weighted_snr(s, pairs, rho_bins):
    """S/N at the published amplitude for the optimal pair weight w = rho_i rho_j (pairs treated as independent,
    var(cos) = 1/2): S/N = dC sqrt(2 sum (rho_i rho_j)^2)."""
    if len(pairs) == 0:
        return 0.0
    r = rho_bins[np.clip(np.digitize(np.asarray(s["fracpol"], float), FP_BINS) - 1, 0, len(rho_bins) - 1)]
    w = r[pairs[:, 0]] * r[pairs[:, 1]]
    return DC_JVAS * np.sqrt(2 * np.sum(w ** 2))


def report(name, s, rho, rho_bins=None):
    x = c.unit(s["ra"], s["dec"])
    dc = c.comoving_distance(np.asarray(s["z"], float))
    P = x * dc[:, None]
    print(f"\n[{name}] N = {len(s)}; z median {np.median(s['z']):.2f}")
    out = {}
    p2 = c.pairs_within(x, c.chord(THETA_J), x_unit=x, min_ang_deg=MIN_ANG)
    m, sd = null_sigma(x, p2[:, 0], p2[:, 1])
    exp = DC_JVAS * rho ** 2
    print(f"  2D theta <= {THETA_J} deg: pairs {len(p2)}; sigma(C) {sd:.4f}; expected {exp:.4f} -> {exp/sd:.1f} sigma;"
          f" rho = 1 (no transfer loss) {DC_JVAS/sd:.1f} sigma; rho_i rho_j-weighted {weighted_snr(s, p2, rho_bins):.1f} sigma")
    if "z" not in s.colnames or not np.all(np.isfinite(s["z"])):
        return out
    for R in SCALES:
        p = c.pairs_within(P, R, x_unit=x, min_ang_deg=MIN_ANG)
        if len(p) < 2:
            print(f"  3D R <= {R:.0f}: pairs {len(p)}")
            continue
        m, sd = null_sigma(x, p[:, 0], p[:, 1])
        # f(R): of the 3.3-deg pairs, the fraction within R in 3D
        d3 = np.linalg.norm(P[p2[:, 0]] - P[p2[:, 1]], axis=1)
        f = np.mean(d3 <= R)
        print(f"  3D R <= {R:.0f} Mpc: pairs {len(p)}; sources in a pair {len(np.unique(p))}; sigma(C) {sd:.4f};"
              f" expected (amplitude as published) {exp:.4f} -> {exp/sd:.1f} sigma; f(R) = {f:.4f} -> 3D-only reading"
              f" {min(exp / max(f, 1e-9), 1):.3f} -> {min(exp / max(f, 1e-9), 1)/sd:.1f} sigma; rho = 1: as published"
              f" {DC_JVAS/sd:.1f}, 3D-only {min(DC_JVAS / max(f, 1e-9), 1)/sd:.1f} sigma; weighted as published"
              f" {weighted_snr(s, p, rho_bins):.1f} sigma")
        out[R] = (len(p), sd)
    return out


def main(d):
    s = Table.read(f"{d}/sample.fits")
    # transfer factor on JVAS-matched sources (all, not only quasars)
    j = s[s["jvas_idx"] >= 0]
    dchi = np.radians(np.asarray(j["derot_polangle"], float) - np.asarray(j["jvas_PA"], float))
    rho = np.mean(np.cos(2 * dchi))
    rho_e = np.std(np.cos(2 * dchi)) / np.sqrt(len(j))
    dchi_obs = np.radians(np.asarray(j["polangle"], float) - np.asarray(j["jvas_PA"], float))
    print(f"JVAS same-source transfer: N = {len(j)}; rho = <cos 2(chi0 - PA_JVAS)> = {rho:.3f} +- {rho_e:.3f};"
          f" (observed-frame angle, not derotated: {np.mean(np.cos(2*dchi_obs)):.3f})")
    rb = rho_by_fp(j)
    q = s[np.isfinite(s["z"])]
    report("all matched quasars", q, rho, rb)
    report("flat-spectrum quasars (headline)", q[q["flat"]], rho, rb)
    report("flat & compact quasars", q[q["flat"] & q["compact"]], rho, rb)
    # what would have the power: the JVAS claim itself is 2D (no redshifts); all flat-spectrum SPICE sources
    fl = s[s["flat"]].copy()
    fl["z"] = np.nan
    report("all flat-spectrum SPICE-RACS sources, 2D only (no redshift; outside the brief's sample)", fl, rho, rb)


if __name__ == "__main__":
    main(sys.argv[1])
