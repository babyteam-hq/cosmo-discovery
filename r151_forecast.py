"""R-151 step 0c: forecast separation of "power at each stream's own epicyclic scale" from "broadband power".

Inputs: r151_epicycle.json (predicted epicyclic band per stream, r151_epicycle.py) and the streampepperdf mocks
(data/r151/pepper_<stream>.npz, r151_pepper.py). For a stream without its own mock, the GD-1 mock's
fractional density field is cut to a random segment of the stream's arm length (its scales kept in kpc).

Estimator (fixed here, as the method would use it): counts in n_bins equal bins of path length s over the arm;
background subtracted with its expectation; normalised by a cubic fit (Banik et al.); periodogram of delta(s)
with a Hann taper; modes n = 1..n_max with wavelength L/n >= lam_min (1 kpc). Shot noise (N_s + N_b)/N_s^2 per
mode (scaled by the taper) is subtracted.

Statistic: f = sum_{band} (P_n - N) / sum_{n<=n_max} (P_n - N); band = modes whose wavelength lies in the stream's
predicted epicyclic range. Combined: mean of f over streams.

Models: (sub) streampepperdf density at CDM-relative abundance 0.2 and 1.0; (epi) BEST CASE for separation: a
single sinusoid at y_C (drawn log-uniform in the predicted range, random phase) holding the same total variance as
the subhalo field at abundance 1. A real epicyclic pattern (harmonics, envelope, washing out away from the
progenitor; Banik et al.) spreads its power, and bar / GMC / Sagittarius perturbations add broadband power to both
families, so both effects only lower the separation computed here.

Separation: z = sqrt(2) Phi^-1(AUC), AUC = P(f_epi > f_sub) (distribution-free; the ratio has heavy tails).
"""
import json
import sys
import numpy as np

rng = np.random.default_rng(151)
EPI = json.load(open("research/cosmo-discovery/r151_epicycle.json"))

# member counts in the arm (excess stars) and background in the same window; sources in results-R151.md F2
COUNTS = {
    "Phoenix": (100, 300),        # DES Y3/Y6 r<24.2-24.5: ~5-7 stars/deg over 14 deg (Boone+ 2510.07511 Fig 10)
    "ATLAS-AliqaUma": (300, 900),  # stellar mass 3.3x Phoenix at ~1.2x the distance (Patrick+ 2206.04241 Table)
    "GD-1": (1500, 1500),         # Gaia DR2+PS1 3-sigma sample 868-1171 (Ibata+ 2020); deeper with Legacy Surveys
}
NBINS = 40
LAM_MIN = 1.0
NMC = 2000


def periodogram(d):
    w = np.hanning(len(d))
    F = np.fft.rfft(d * w)
    return np.abs(F) ** 2 / np.sum(w**2)  # per-mode power, white noise of variance v -> v


def measure(delta_true, L, Ns, Nb):
    """Draw counts from (1+delta) members + uniform background, return the band-fraction ingredients."""
    n = len(delta_true)
    lam_m = Ns * (1 + delta_true) / delta_true.size / np.mean(1 + delta_true)
    c = rng.poisson(np.clip(lam_m, 0, None)) + rng.poisson(Nb / n, size=n)
    m = c - Nb / n
    x = np.linspace(-1, 1, n)
    smooth = np.polyval(np.polyfit(x, m, 3), x)
    smooth = np.where(np.abs(smooth) < 1e-3 * Ns / n, 1e-3 * Ns / n, smooth)
    d = m / smooth - 1
    P = periodogram(d)
    noise = (Ns + Nb) / Ns**2 * n  # per-bin variance of d ~ (Ns+Nb)/n / (Ns/n)^2 ; periodogram of white -> variance
    return P, noise


def band_modes(name, L):
    lo, hi = EPI[name]["yC_kpc"]
    nmax = int(np.floor(L / LAM_MIN))
    modes = np.arange(1, nmax + 1)
    lam = L / modes
    inband = (lam >= lo) & (lam <= hi)
    return modes, inband


def fraction(P, noise, modes, inband):
    ex = P[modes] - noise
    tot = ex.sum()
    return ex[inband].sum() / tot if tot > 0 else np.nan


def load_mock(name):
    try:
        z = np.load(f"research/cosmo-discovery/data/r151/pepper_{name}.npz")
        src = name
    except FileNotFoundError:
        z = np.load("research/cosmo-discovery/data/r151/pepper_GD-1.npz")
        src = "GD-1 (segment)"
    return z, src


def resample(z, key, n, L=None):
    """Mock fractional density on n equal path-length bins: over the whole arm, or (template from another
    stream's mock) over a random segment of length L kpc of it, keeping the field's scales in kpc."""
    s = z["s"]
    d0 = z["dens0"]
    out = []
    for D in z[key]:
        frac = np.where(d0 > 0, D / np.where(d0 > 0, d0, 1), 1) - 1
        if L is not None and L < s[-1] - s[0]:
            a = rng.uniform(s[0], s[-1] - L)
            sg = np.linspace(a, a + L, n + 1)
        else:
            sg = np.linspace(s[0], s[-1], n + 1)
        sc = 0.5 * (sg[1:] + sg[:-1])
        out.append(np.interp(sc, s, frac))
    return np.array(out)


def sep(e, s_):
    """Distribution-free separation: AUC = P(f_epi > f_sub) over all pairs, z = sqrt(2) Phi^-1(AUC)
    (equal to d' for two equal-width Gaussians). Undefined draws (total excess power <= 0) count as f = 0,
    i.e. 'no epicyclic concentration seen'."""
    from scipy.stats import norm, mannwhitneyu
    e = np.where(np.isfinite(e), e, 0.0)
    s_ = np.where(np.isfinite(s_), s_, 0.0)
    auc = mannwhitneyu(e, s_, alternative="greater").statistic / (e.size * s_.size)
    auc = min(auc, 1 - 0.5 / e.size)
    return np.sqrt(2) * norm.ppf(auc), np.median(e), np.median(s_)


def run(names, amp, noiseless=False, shape="1.0"):
    """Per stream and combined separation with both models scaled to total delta rms `amp` (None: CDM native)."""
    fe_all, fs_all, per = [], [], {}
    for name in names:
        L = EPI[name]["L_kpc"]
        Ns, Nb = COUNTS[name]
        if noiseless:
            Ns, Nb = 10**9, 0
        modes, inband = band_modes(name, L)
        z, src = load_mock(name)
        fields = resample(z, f"dens_{shape}", NBINS, None if src == name else L)
        rms = np.sqrt(np.mean(fields.var(axis=1)))
        A = rms if amp is None else amp
        fields = fields * (A / rms)
        fs = np.array([fraction(*measure(fields[k % len(fields)], L, Ns, Nb), modes, inband) for k in range(NMC)])
        lo, hi = EPI[name]["yC_kpc"]
        sc = (np.arange(NBINS) + 0.5) / NBINS * L
        fe = []
        for k in range(NMC):
            y = np.exp(rng.uniform(np.log(lo), np.log(hi)))
            d = np.sqrt(2) * A * np.cos(2 * np.pi * sc / y + rng.uniform(0, 2 * np.pi))
            fe.append(fraction(*measure(d, L, Ns, Nb), modes, inband))
        fe = np.array(fe)
        per[name] = sep(fe, fs) + (src, int(inband.sum()), len(modes))
        fe_all.append(fe)
        fs_all.append(fs)
    ce = np.mean(np.nan_to_num(np.array(fe_all)), axis=0)
    cs = np.mean(np.nan_to_num(np.array(fs_all)), axis=0)
    return per, sep(ce, cs)


def main():
    names = sys.argv[1:] or list(COUNTS)
    print("separation (sigma) of best-case epicycle vs subhalo (streampepperdf shape), median f_epi / f_sub")
    for shape in ("1.0", "0.2"):
        for amp, nl in ((None, False), (0.1, False), (0.2, False), (0.3, False), (0.3, True)):
            per, comb = run(names, amp, nl, shape)
            lab = f"shape n/n_CDM={shape} amp={'CDM-native' if amp is None else amp}{' noiseless' if nl else ''}"
            txt = "; ".join(f"{n} {v[0]:.2f} ({v[1]:.2f}/{v[2]:.2f}, band {v[4]}/{v[5]} modes)" for n, v in per.items())
            print(f"{lab:42s} combined {comb[0]:.2f} ({comb[1]:.2f}/{comb[2]:.2f}) | {txt}", flush=True)
    # shuffle check: overlap of predicted bands (mode sets) between streams
    for i, n in enumerate(names):
        for m in names:
            if m == n:
                continue
            L = EPI[n]["L_kpc"]
            lo, hi = EPI[m]["yC_kpc"]
            lam = L / np.arange(1, int(L / LAM_MIN) + 1)
            own = band_modes(n, L)[1]
            print(f"shuffle: {n} scored in {m}'s band: {np.sum(own & (lam >= lo) & (lam <= hi))}/{own.sum()} of its own band modes shared")


if __name__ == "__main__":
    main()
