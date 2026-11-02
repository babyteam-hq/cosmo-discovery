"""R-158: train the D4-symmetrised, reflection-antisymmetric CNN on GZ1 handedness labels (Legacy DR10 r stamps).
Augmentation: rotation, scale, blur and noise (to reach space-telescope z~1 resolution and S/N); flip+label flip.
Saves the model and a validation table (accuracy vs |f|, clean and degraded)."""
import os, sys, numpy as np, torch, torch.nn.functional as F, time
from r158_common import DATA
import r158_classify as C
torch.set_num_threads(4); torch.manual_seed(0); np.random.seed(0)
d = np.load(os.path.join(DATA, "gz1", "train.npz"))
pf = os.path.join(DATA, "gz1", "prepped.npz")
if not os.path.exists(pf):
    X, Y, keep = [], [], []
    for i, (a, y) in enumerate(zip(d["x"], d["y"])):
        x, info = C.prep(a)
        if x is None: continue
        X.append(x); Y.append(y); keep.append(i)
    np.savez_compressed(pf, x=np.array(X), y=np.array(Y), keep=np.array(keep))
p = np.load(pf); X, Y = torch.tensor(p["x"])[:, None], torch.tensor((p["y"] > 0).astype(np.float32))
n = len(X); idx = torch.randperm(n); nv = n // 7
vi, ti = idx[:nv], idx[nv:]
print("train", len(ti), "val", nv, "frac CW", Y.mean().item(), flush=True)

def blur(x, s):
    if s < 0.3: return x
    k = int(3 * s) * 2 + 1; t = torch.arange(k) - k // 2; g = torch.exp(-t ** 2 / (2 * s * s)); g = g / g.sum()
    x = F.conv2d(x, g.view(1, 1, 1, k), padding=(0, k // 2)); return F.conv2d(x, g.view(1, 1, k, 1), padding=(k // 2, 0))

def augment(x, y, strong=1.0):
    B = x.shape[0]
    th = torch.rand(B) * 2 * np.pi; sc = 1 + (torch.rand(B) - 0.5) * 0.3
    A = torch.zeros(B, 2, 3); A[:, 0, 0] = torch.cos(th) * sc; A[:, 0, 1] = -torch.sin(th) * sc
    A[:, 1, 0] = torch.sin(th) * sc; A[:, 1, 1] = torch.cos(th) * sc
    x = F.grid_sample(x, F.affine_grid(A, x.shape, align_corners=True), align_corners=True)
    out = x.clone()
    lev = np.random.randint(0, 4, B)   # blur levels 0, 1, 2, 3 grid px (resolution)
    for L in (1, 2, 3):
        k = torch.tensor(lev == L)
        if k.any(): out[k] = blur(x[k], L * strong)
    out = out + torch.randn_like(out) * torch.rand(B, 1, 1, 1) * 0.12 * strong
    fl = torch.rand(B) < 0.5
    out[fl] = torch.flip(out[fl], (3,)); y = y.clone(); y[fl] = 1 - y[fl]
    return out, y

net = C.make_net(); opt = torch.optim.AdamW(net.parameters(), 2e-3, weight_decay=1e-4)
E = int(sys.argv[1]) if len(sys.argv) > 1 else 12
sched = torch.optim.lr_scheduler.OneCycleLR(opt, 3e-3, total_steps=E * (len(ti) // 64 + 1))
t0 = time.time()
for ep in range(E):
    net.train(); perm = ti[torch.randperm(len(ti))]; L = 0
    for b in range(0, len(perm), 64):
        j = perm[b:b + 64]; xb, yb = augment(X[j], Y[j])
        loss = F.binary_cross_entropy_with_logits(net(xb), yb)
        opt.zero_grad(); loss.backward(); opt.step(); sched.step(); L += loss.item() * len(j)
    net.eval()
    with torch.no_grad():
        fv = torch.cat([net(X[vi[b:b + 256]]) for b in range(0, nv, 256)])
    acc = ((fv > 0).float() == Y[vi]).float().mean().item()
    print(f"ep {ep} loss {L / len(ti):.4f} val acc {acc:.4f} t {time.time() - t0:.0f}s", flush=True)
torch.save(net.state_dict(), os.path.join(DATA, "net.pt"))
np.savez(os.path.join(DATA, "val_idx.npz"), vi=vi.numpy())
