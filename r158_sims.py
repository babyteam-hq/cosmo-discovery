"""R-158 null 4: simulated achiral galaxies convolved with each dataset's empirical PSF, through prep and both classifiers.
Model: exponential disc + n=4 bulge + elliptical ring sharing the disc's inclination and PA, plus a bar at its own PA.
Two centred elliptical components are mirror-symmetric about their bisector, so every model is achiral; any classifier
asymmetry on them comes from the PSF (or pixelisation). Also run with the mirrored PSF (expect the opposite sign).
Sizes and S/N are drawn from the dataset's own prep distribution. Usage: r158_sims.py cosmos|euclid N"""
import os, sys, json, numpy as np, torch
from scipy import ndimage as ndi
from scipy.signal import fftconvolve
from r158_common import DATA
import r158_classify as C
ds, N = sys.argv[1], int(sys.argv[2])
rng = np.random.default_rng(11)
from r158_psf import psf_from, stars
st, crop = stars(ds); psf, nst = psf_from(st, crop)
from r158_common import prep_subset
pr = prep_subset(ds)
def spiral(n, r50, q, pa, hand, amp, pitch, bt, nclump, rng):
    """Chiral model: exponential disc whose light is modulated by an m=2 logarithmic spiral of handedness `hand` (+1/-1),
    plus a bulge and random clumps; inclined by q at PA. Mirror of the model = the model with -hand."""
    yy, xx = np.mgrid[0:n, 0:n] - (n - 1) / 2.0
    u = xx * np.cos(pa) + yy * np.sin(pa); v = (-xx * np.sin(pa) + yy * np.cos(pa)) / q
    r = np.hypot(u, v) + 1e-3; th = np.arctan2(v, u)
    Rd = r50 / 1.678
    disc = np.exp(-r / Rd) * (1 + amp * np.cos(2 * (th - hand * np.log(r / Rd) / np.tan(pitch))))
    Rb = 0.25 * Rd; img = (1 - bt) * disc / disc.sum() + bt * np.exp(-7.67 * ((r / Rb) ** 0.25)) / np.exp(-7.67 * ((r / Rb) ** 0.25)).sum()
    for _ in range(nclump):
        rc = rng.uniform(0.3, 2.5) * Rd; tc = rng.uniform(0, 2 * np.pi)
        uc, vc = rc * np.cos(tc), rc * np.sin(tc)
        img += rng.uniform(0.01, 0.05) * np.exp(-0.5 * ((u - uc) ** 2 + (v - vc) ** 2) / (0.15 * Rd) ** 2) / (2 * np.pi * (0.15 * Rd) ** 2 * q)
    return img
def model(n, r50, q, pa, pab, bt, fbar, fring, rr):
    yy, xx = np.mgrid[0:n, 0:n] - (n - 1) / 2.0
    def ell(pa_, q_):
        u = xx * np.cos(pa_) + yy * np.sin(pa_); v = -xx * np.sin(pa_) + yy * np.cos(pa_); return np.hypot(u, v / q_)
    Rd = r50 / 1.678
    re = ell(pa, q)
    img = (1 - bt - fbar - fring) * np.exp(-re / Rd) / (2 * np.pi * Rd ** 2 * q)
    Rb = 0.25 * Rd; img += bt * np.exp(-7.67 * ((ell(pa, max(q, 0.7)) / Rb) ** 0.25 - 1)) / (2 * np.pi * Rb ** 2 * 200)
    rb = ell(pab, 0.25 * q + 0.05); img += fbar * np.exp(-0.5 * (rb / Rd) ** 2) / (2 * np.pi * Rd ** 2 * 0.3)
    img += fring * np.exp(-0.5 * ((re - rr * Rd) / (0.25 * Rd)) ** 2) / (2 * np.pi * rr * Rd * 0.6 * Rd * q)
    return img
net = C.make_net(); net.load_state_dict(torch.load(os.path.join(DATA, "net.pt"))); net.eval()
INDIV = len(sys.argv) > 3 and sys.argv[3] == "indiv"
if INDIV:
    kern = []
    for s_ in st:
        k_, n_ = psf_from([s_], crop)
        if n_: kern.append(k_)
    print("individual star kernels", len(kern))
res = {}
for tag, P in (("psf", psf), ("mirror_psf", psf[:, ::-1])):
    X = []
    rng = np.random.default_rng(11)
    for i in range(N):
        j = rng.integers(len(pr["r50"])); r50 = pr["r50"][j]; snr = pr["snr"][j]
        n = int(pr["px"][j])
        q = rng.uniform(0.4, 1); pa = rng.uniform(0, np.pi); pab = rng.uniform(0, np.pi)
        m = model(n, r50, q, pa, pab, rng.uniform(0, 0.4), rng.uniform(0, 0.2) * (rng.random() < 0.5),
                  rng.uniform(0, 0.25) * (rng.random() < 0.5), rng.uniform(1.5, 3))
        if INDIV:
            K = kern[rng.integers(len(kern))]; m = fftconvolve(m, K if tag == "psf" else K[:, ::-1], mode="same")
        else:
            m = fftconvolve(m, P, mode="same")
        x0, i0 = C.prep(m + rng.normal(0, 1e-6 * m.max(), m.shape))
        if "snr" not in i0: X.append(None); continue
        sig = 1e-6 * m.max() * i0["snr"] / snr
        x, _ = C.prep(m + rng.normal(0, sig, m.shape))
        X.append(x)
    ok = np.array([x is not None for x in X]); X = np.array([x if x is not None else np.zeros((64, 64), np.float32) for x in X])
    with torch.no_grad():
        f = torch.cat([net(torch.tensor(X[b:b + 256])[:, None]) for b in range(0, len(X), 256)]).numpy()
    S = np.array([C.fft_handedness(x) for x in X])
    res[tag] = {"f": f[ok].tolist(), "S": S[ok].tolist()}
# chiral injections: true handedness drawn with P(+1) = 1/2 (eps, armless rate) at the sample's R50 and S/N
rng = np.random.default_rng(12)
X, H, R, SN = [], [], [], []
for i in range(N):
    j = rng.integers(len(pr["r50"])); r50 = pr["r50"][j]; snr = pr["snr"][j]; n = int(pr["px"][j])
    hand = 1 if rng.random() < 0.5 else -1
    m = spiral(n, r50, rng.uniform(0.4, 1), rng.uniform(0, np.pi), hand, rng.uniform(0.15, 0.8),
               np.deg2rad(rng.uniform(10, 35)), rng.uniform(0, 0.4), rng.integers(0, 6), rng)
    m = fftconvolve(m, psf, mode="same")
    x0, i0 = C.prep(m + rng.normal(0, 1e-6 * m.max(), m.shape))
    if "snr" not in i0: continue
    x, _ = C.prep(m + rng.normal(0, 1e-6 * m.max() * i0["snr"] / snr, m.shape))
    if x is None: continue
    X.append(x); H.append(hand); R.append(r50); SN.append(snr)
X = np.array(X)
with torch.no_grad():
    f = torch.cat([net(torch.tensor(X[b:b + 256])[:, None]) for b in range(0, len(X), 256)]).numpy()
S = np.array([C.fft_handedness(x) for x in X])
res["chiral"] = {"f": f.tolist(), "S": S.tolist(), "hand": H, "r50": R, "snr": SN}
json.dump({"n_stars": nst, "psf_kernel": psf.tolist(), **res}, open(os.path.join(DATA, f"sims_{ds}{'_indiv' if INDIV else ''}.json"), "w"))
print(ds, "stars in PSF", nst, "sims", len(res["psf"]["f"]))
