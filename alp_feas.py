"""R-153 Step 0: feasibility numbers per dwarf from the reduced images (pattern-subtracted).

Prints, per dwarf: exposures, the expected 95% limit from (a) the pipeline variance alone (photon + read noise as the
pipeline reports it), (b) white noise rescaled by the robust chi2 of the continuum residuals, (c) the sign-flip
ensemble (the real noise, correlated within exposures); the dwarf's starlight (ON-OFF continuum) in the inner cells and
the fitted narrow stellar-feature amplitude beta (as a fraction of the local continuum) with its flip sigma.
Saves data/alp/search/feas_<name>.npz. Usage: python -u alp_feas.py <name> [<red-dir>]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import alp_cells as C  # noqa: E402
import alp_search as S  # noqa: E402


def main():
    name = sys.argv[1]
    red = sys.argv[2] if len(sys.argv) > 2 else S.REDP
    ent, masses, T, N, summ = S.prepare(name, red=red)
    dec = S.decompose(ent, T, N)
    rng = np.random.default_rng(153)
    sg = rng.choice([-1.0, 1.0], size=(S.NFLIP, len(dec["exps"])))
    A = S.solve(dec)
    Af = S.solve(dec, sg)
    sF = np.nanstd(Af, 0)
    sW = S.white_sigma(dec)
    scale = np.ones(len(masses))
    # pipeline-only sigma: undo the per-detector variance rescaling at the detector holding each mass's line
    sc = {k: v[1] for k, v in summ.items()}
    lam0 = 2 * 1.23984198 / masses
    for k, (lo, hi) in C.DET_EDGES.items():
        scale[(lam0 >= lo) & (lam0 < hi)] = sc.get(k, 1.0)
    sP = sW / np.sqrt(scale)
    # stellar feature amplitude beta: solve for the nuisance coefficient too
    TT, TN, NN = dec["TT"], dec["TN"], dec["NN"]
    det = TT * NN - TN ** 2

    def beta(Ty, Ny):
        return (TT * Ny - TN * Ty) / det
    b = beta(dec["Ty_e"].sum(0), dec["Ny_e"].sum(0))
    bf = np.nanstd(beta(sg @ dec["Ty_e"], sg @ dec["Ny_e"]), 0)
    cont = ent["y"] - C.continuum_residuals(ent, S.DEG)
    star = {}
    for lo, hi in ((0.9, 1.1), (1.5, 1.7), (2.1, 2.3), (3.3, 3.6), (4.5, 4.8)):
        m = (ent["lam"] >= lo) & (ent["lam"] < hi)
        star[f"{lo}-{hi}"] = [float(np.nanmedian(cont[m & (ent["theta"] < a2) & (ent["theta"] >= a1)]))
                              for a1, a2 in ((0, 0.1), (0.1, 0.2), (0.3, 0.6))]
    os.makedirs(S.OUT, exist_ok=True)
    np.savez(os.path.join(S.OUT, f"feas_{name}.npz"), masses=masses, A=A, sF=sF, sW=sW, sP=sP, beta=b, beta_sF=bf,
             nexp=len(dec["exps"]))
    g = lambda s: 1e-11 * np.sqrt(1.645 * s)
    print(f"{name}: exposures {len(dec['exps'])}, entries {len(ent['y'])}, clean {summ}")
    print("  starlight ON-OFF continuum nW (theta 0-0.1, 0.1-0.2, 0.3-0.6 deg):", {k: [round(x, 2) for x in v] for k, v in star.items()})
    for lo, hi in ((0.5, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 2.5), (2.5, 3.3)):
        s = (masses >= lo) & (masses < hi)
        print(f"  m {lo}-{hi} eV: g95 expected  pipeline-var {np.nanmedian(g(sP[s])):.2e}  white(rescaled) "
              f"{np.nanmedian(g(sW[s])):.2e}  sign-flip {np.nanmedian(g(sF[s])):.2e}   flip/pipeline sigma "
              f"{np.nanmedian(sF[s]/sP[s]):.1f};  stellar beta/sigma rms {np.nanstd(b[s]/bf[s]):.2f}, "
              f"|beta| median {np.nanmedian(np.abs(b[s])):.4f}")


if __name__ == "__main__":
    main()
