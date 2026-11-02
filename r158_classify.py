"""R-158 classifiers. Everything here is exactly mirror-equivariant by construction:

prep(stamp): background, segmentation, R50 and resampling on a grid symmetric about the stamp centre, so
    prep(fliplr(stamp)) == fliplr(prep(stamp)) up to float rounding.
fft_handedness(img): logarithmic-spiral Fourier power at positive vs negative pitch (m = 1..3), exactly odd under a flip.
Net: a CNN g antisymmetrised over reflection, f(x) = [g(x) - g(fliplr x)]/2, so f(fliplr x) = -f(x) exactly
    (rotation invariance is learned from rotation augmentation; it is not needed for exactness).
Sign convention: f > 0 and S > 0 mean the label the network was trained to call +1 (GZ1 'CW', as seen on the sky).
Stamps must be in sky parity (det CD < 0) with array row index increasing to the north-ish (FITS, origin lower)."""
import numpy as np
from scipy import ndimage as ndi

G = 64            # output grid
HALF = 3.0        # output half-width in units of R50

def sky_parity(a, header):
    cd = np.array([[header.get("CD1_1", header.get("CDELT1", -1) * header.get("PC1_1", 1)), header.get("CD1_2", 0.0)],
                   [header.get("CD2_1", 0.0), header.get("CD2_2", header.get("CDELT2", 1) * header.get("PC2_2", 1))]])
    return a[:, ::-1] if np.linalg.det(cd) > 0 else a

def prep(a, min_r50=1.5):
    """a: 2D stamp centred on the target (square, any size). Returns (64x64 float32 image, info dict) or (None, info)."""
    a = np.nan_to_num(np.asarray(a, dtype=np.float64))
    n = a.shape[0]
    c = (n - 1) / 2.0
    yy, xx = np.mgrid[0:n, 0:n]
    rr = np.hypot(xx - c, yy - c)
    border = rr > 0.42 * n
    b = a[border]
    med = np.median(b); mad = 1.4826 * np.median(np.abs(b - med))
    for _ in range(3):
        k = np.abs(b - med) < 3 * mad + 1e-30
        med = np.median(b[k]); mad = 1.4826 * np.median(np.abs(b[k] - med))
    a = a - med
    s = ndi.gaussian_filter(a, 1.0, mode="nearest")
    sb = s[border]; ss = 1.4826 * np.median(np.abs(sb - np.median(sb))) + 1e-30
    lab, nl = ndi.label(s > 2.5 * ss)
    ic = int(round(c)); core = lab[int(np.floor(c)):int(np.ceil(c)) + 1, int(np.floor(c)):int(np.ceil(c)) + 1]
    ids = [v for v in np.unique(core) if v > 0]
    if not ids:
        win = lab[max(0, ic - 3):ic + 4, max(0, ic - 3):ic + 4]
        ids = [v for v in np.unique(win) if v > 0]
    if not ids:
        return None, {"why": "no detection"}
    main = np.isin(lab, ids)
    other = (lab > 0) & ~main
    other = ndi.binary_dilation(other, iterations=2) & ~main
    a[other] = 0.0
    f = np.clip(a[main], 0, None); r = rr[main]
    o = np.argsort(r); cf = np.cumsum(f[o])
    if cf[-1] <= 0:
        return None, {"why": "no flux"}
    r50 = float(r[o][np.searchsorted(cf, 0.5 * cf[-1])])
    r90 = float(r[o][min(len(o) - 1, np.searchsorted(cf, 0.9 * cf[-1]))])
    snr = float(cf[-1] / (mad * np.sqrt(main.sum()) + 1e-30))
    info = {"r50_pix": r50, "r90_pix": r90, "snr": snr, "npix": int(main.sum()), "edge": bool(main[border].any())}
    if r50 < min_r50:
        info["why"] = "small"; return None, info
    W = HALF * r50
    g = c + W * (2 * np.arange(G) / (G - 1) - 1)
    gy, gx = np.meshgrid(g, g, indexing="ij")
    out = ndi.map_coordinates(a, [gy, gx], order=1, mode="constant", cval=0.0)
    m = ndi.map_coordinates(main.astype(float), [gy, gx], order=1, mode="constant", cval=0.0) > 0.5
    scale = np.percentile(out[m], 99) if m.sum() > 10 else out.max()
    if not scale > 0:
        return None, {"why": "scale"}
    out = np.arcsinh(10 * out / scale) / np.arcsinh(10)
    return out.astype(np.float32), info

_gy, _gx = np.mgrid[0:G, 0:G]
_c = (G - 1) / 2.0
_R = np.hypot(_gx - _c, _gy - _c); _T = np.arctan2(_gy - _c, _gx - _c)
_R50 = (G - 1) / 2.0 / HALF
_ANN = (_R > 0.5 * _R50) & (_R < 2.8 * _R50)
_P = np.arange(1, 21, dtype=float)

def fft_handedness(img):
    """S = (P+ - P-)/(P+ + P-) with P+- the log-spiral power at pitch +-p (p = 1..20), m = 1, 2, 3.
    Uses the (ln r, theta) measure: pixel weight 1/r^2. A mirror (theta -> pi - theta) maps p to -p, so S is odd."""
    I = np.clip(img[_ANN], 0, None) / _R[_ANN] ** 2
    lr = np.log(_R[_ANN]); th = _T[_ANN]
    tot = I.sum() + 1e-30
    Pp = Pm = 0.0
    for m in (1, 2, 3):
        e = np.exp(-1j * m * th)
        ph = np.exp(-1j * np.outer(_P, lr))
        Ap = (ph * (I * e)).sum(1) / tot        # pitch +p
        Am = (np.conj(ph) * (I * e)).sum(1) / tot  # pitch -p
        Pp += (np.abs(Ap) ** 2).sum(); Pm += (np.abs(Am) ** 2).sum()
    return float((Pp - Pm) / (Pp + Pm + 1e-30))

def make_net():
    import torch, torch.nn as nn
    class G_(nn.Module):
        def __init__(s):
            super().__init__()
            def blk(i, o): return [nn.Conv2d(i, o, 3, padding=1), nn.BatchNorm2d(o), nn.ReLU()]
            s.f = nn.Sequential(*blk(1, 12), nn.MaxPool2d(2), *blk(12, 24), *blk(24, 24), nn.MaxPool2d(2),
                                *blk(24, 48), *blk(48, 48), nn.MaxPool2d(2), *blk(48, 48), nn.AdaptiveAvgPool2d(1))
            s.h = nn.Linear(48, 1)
        def forward(s, x): return s.h(s.f(x).flatten(1))[:, 0]
    class Net(nn.Module):
        def __init__(s):
            super().__init__(); s.g = G_()
        def forward(s, x):
            B = x.shape[0]
            o = s.g(torch.cat([x, torch.flip(x, (3,))], 0)).view(2, B)
            return (o[0] - o[1]) / 2.0
    return Net()

def fft_handedness_batch(X, chunk=4096):
    """Vectorised fft_handedness over an array of 64x64 images (identical formula)."""
    lr = np.log(_R[_ANN]); th = _T[_ANN]; w = 1.0 / _R[_ANN] ** 2
    ph = np.exp(-1j * np.outer(_P, lr))                  # (20, Nann)
    E = [np.exp(-1j * m * th) for m in (1, 2, 3)]
    out = np.empty(len(X))
    for b in range(0, len(X), chunk):
        I = np.clip(np.asarray(X[b:b + chunk], dtype=np.float64).reshape(-1, G * G)[:, _ANN.ravel()], 0, None) * w
        tot = I.sum(1, keepdims=True) + 1e-30
        Pp = Pm = 0.0
        for e in E:
            Ie = I * e
            Ap = Ie @ ph.T / tot; Am = Ie @ np.conj(ph).T / tot
            Pp = Pp + (np.abs(Ap) ** 2).sum(1); Pm = Pm + (np.abs(Am) ** 2).sum(1)
        out[b:b + chunk] = (Pp - Pm) / (Pp + Pm + 1e-30)
    return out
