"""R-153 Step 0: the expected joint 95% limit (Asimov: A_hat = 0) for Draco + Ursa Minor from their sign-flip errors,
the D-factor log-uncertainty (Pace & Strigari D(0.5 deg), mean of the two-sided errors) profiled; and an upper bound on
what the other six Song+2024 dwarfs could add, scaling Draco's measured sigma by D(0.5)^2 and exposure count
(sigma_d ~ sigma_Draco x (D_Draco / D_d) x sqrt(N_Draco / N_d)), which flatters them (no zodi or coverage penalty).
Compared with arXiv:2609.10674 (final, Fig. 9 right, read off) and arXiv:2503.11753 (as drawn there).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import alp_model as M  # noqa: E402
import alp_search as S  # noqa: E402

JWST_WEAK = {1.0: 2.3e-11, 1.5: 2.2e-11, 2.0: 1.3e-11, 2.5: 1.2e-11}    # 2609.10674 final
JWST_STRONG = {1.0: 6.7e-13, 1.5: 5.6e-13, 2.0: 6.7e-13, 2.5: 4.2e-13}  # 2503.11753 as drawn in 2609.10674
NEXP = {"draco": 486, "umi": 519, "sculptor": 542, "fornax": 492, "sextans": 318, "carina": 583, "leo1": 307, "leo2": 407}


def main():
    f = {n: np.load(os.path.join(S.OUT, f"feas_{n}.npz")) for n in ("draco", "umi")}
    masses = f["draco"]["masses"]
    lsig = {n: np.mean(M.DWARFS[n]["ps_err"]) for n in M.DWARFS}
    rows = []
    for i, m in enumerate(masses):
        sig = [f["draco"]["sF"][i], f["umi"]["sF"][i]]
        g2, _, _ = S.profile_limit([0.0, 0.0], sig, [lsig["draco"], lsig["umi"]])
        sd = f["draco"]["sF"][i]
        extra = [sd * 10 ** (M.DWARFS["draco"]["ps"][2] - M.DWARFS[n]["ps"][2]) * np.sqrt(NEXP["draco"] / NEXP[n])
                 for n in NEXP if n not in ("draco", "umi")]
        g8, _, _ = S.profile_limit([0.0] * 8, sig + extra, [lsig["draco"], lsig["umi"]] +
                                   [lsig[n] for n in NEXP if n not in ("draco", "umi")])
        rows.append((m, 1e-11 * np.sqrt(g2), 1e-11 * np.sqrt(g8)))
    rows = np.array(rows)
    np.save(os.path.join(S.OUT, "gate_expected.npy"), rows)
    print("m_a eV | expected g95 Draco+UMi | + six others (flattering bound) | 2609.10674 final | 2503.11753")
    for m0 in (1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5):
        s = np.abs(rows[:, 0] - m0) < 0.06
        a, b = np.nanmedian(rows[s, 1]), np.nanmedian(rows[s, 2])
        w = np.interp(m0, list(JWST_WEAK), list(JWST_WEAK.values()))
        st = np.interp(m0, list(JWST_STRONG), list(JWST_STRONG.values()))
        print(f"{m0:.2f} | {a:.2e} | {b:.2e} | {w:.1e} | {st:.1e}  -> 8-dwarf / weak JWST = {b/w:.2f}")


if __name__ == "__main__":
    main()
