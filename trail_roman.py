"""Is the DF2-DF4 line special in an independent LSB-galaxy catalogue? Roman, Castilla & Pascual-Granado 2021 (A&A 656, A44;
arXiv:2110.09527) detected 42 LSBGs (mu_g(0) > 24, r_e > 5") over 6x6 deg around NGC 1052 with one uniform procedure.
 (1) count the Roman objects on the strip of half-width w about the line through DF2 (RCP 29) and DF4 (RCP 12), over the
     trail's length (RCP 32 to DF7 ends, padded), excluding the two that define it;
 (2) null A (van Dokkum+2022's): rotate every object about NGC 1052 by a random angle (keeps the radial distribution),
     redraw the line through the rotated DF2 and DF4, count;
 (3) null B (look-elsewhere): the same count for the line through every pair of catalogue objects within the group region;
     where does the DF2-DF4 pair rank;
 (4) the same with the large LSBGs only (r_e > 15"), which the paper calls over-represented near NGC 1052.
usage: python trail_roman.py
"""
import numpy as np, itertools
rng = np.random.default_rng(1)
# id, RA, Dec, r_e("), mu0, note  (Roman+2021 Tables 1-2, transcribed via ar5iv)
CAT = [(1,37.2722,-10.6033,6.0,25.9,''),(2,38.2332,-9.0575,5.4,26.5,''),(3,38.5704,-8.9949,8.3,25.0,''),(4,39.0332,-9.9845,6.0,24.2,''),
(5,39.0430,-8.3067,5.2,24.5,''),(6,39.2677,-5.2661,9.2,26.3,''),(7,39.4813,-6.2566,10.8,24.1,''),(8,39.5983,-8.2217,5.6,26.2,''),
(9,39.6241,-7.9257,12.8,25.9,'DF7'),(10,39.6952,-6.2366,8.1,24.5,''),(11,39.8028,-8.1408,10.0,25.9,'DF5'),(12,39.8128,-8.1160,16.8,24.3,'DF4'),
(13,39.8270,-7.5369,11.0,25.4,''),(14,39.8617,-7.3707,12.3,26.9,''),(15,39.9104,-7.4737,8.1,25.6,''),(16,39.9139,-8.2285,6.1,24.3,'Ta21-12000'),
(17,39.9696,-8.2121,6.4,26.6,''),(18,40.0194,-8.4461,24.7,26.3,'DF1'),(19,40.0341,-7.9473,8.6,24.6,''),(20,40.0820,-7.9847,6.8,24.1,''),
(21,40.1200,-8.2434,8.9,26.0,''),(22,40.1502,-9.4977,5.2,25.5,''),(23,40.1719,-11.0829,9.9,25.5,''),(24,40.1894,-7.6470,7.1,25.0,'DF8'),
(25,40.2277,-5.3300,6.4,24.8,''),(26,40.2897,-8.2968,9.1,24.7,'RCP26'),(27,40.3126,-7.4934,6.8,25.5,''),(28,40.4215,-8.3475,7.1,26.1,''),
(29,40.4451,-8.4028,21.3,24.7,'DF2'),(30,40.4475,-8.7854,7.6,24.4,''),(31,40.6033,-9.4483,6.5,24.7,''),(32,40.6202,-8.3768,23.0,27.8,'RCP32'),
(33,40.6504,-8.0426,9.9,25.3,''),(34,40.6583,-7.3381,17.8,24.1,''),(35,40.6963,-7.7721,8.6,24.3,''),(36,40.7636,-8.0140,6.7,25.5,''),
(37,40.8702,-7.8732,7.0,25.9,''),(38,40.9303,-6.9219,8.4,25.5,''),(39,41.1604,-7.1764,11.2,24.1,''),(40,41.4845,-7.6496,6.4,24.1,''),
(41,41.8189,-8.2923,5.5,24.5,''),(42,43.2093,-7.8612,5.9,25.3,'')]
N1052 = (40.2700, -8.2558)
def tangent(ra, dec, ra0=N1052[0], dec0=N1052[1]):
    ra, dec, ra0, dec0 = map(np.radians, (ra, dec, ra0, dec0))
    cosc = np.sin(dec0)*np.sin(dec) + np.cos(dec0)*np.cos(dec)*np.cos(ra-ra0)
    xi = np.cos(dec)*np.sin(ra-ra0)/cosc; eta = (np.cos(dec0)*np.sin(dec) - np.sin(dec0)*np.cos(dec)*np.cos(ra-ra0))/cosc
    return np.degrees(xi)*60, np.degrees(eta)*60   # arcmin, xi toward +RA
ids = np.array([c[0] for c in CAT]); xi, eta = tangent(np.array([c[1] for c in CAT]), np.array([c[2] for c in CAT]))
re_ = np.array([c[3] for c in CAT]); note = [c[5] for c in CAT]
P = np.column_stack([xi, eta])
i2, i4, i32, i7 = [int(np.where(ids == k)[0][0]) for k in (29, 12, 32, 9)]

def strip_count(P, a, b, w, pad, exclude):
    """objects within perpendicular distance w of the line through P[a], P[b], and along it between the trail ends (+pad)."""
    d = P[b] - P[a]; u = d/np.linalg.norm(d); n = np.array([-u[1], u[0]])
    t = (P - P[a]) @ u; s = (P - P[a]) @ n
    # trail length: from RCP32 to DF7 projected, padded
    t_lo, t_hi = sorted([(P[i32]-P[a])@u, (P[i7]-P[a])@u]); t_lo -= pad; t_hi += pad
    m = (np.abs(s) <= w) & (t >= t_lo) & (t <= t_hi); m[list(exclude)] = False
    return m

W, PAD = 5.2, 5.0
m = strip_count(P, i2, i4, W, PAD, (i2, i4))
print(f"Roman+2021: {len(ids)} LSBGs. On the DF2-DF4 strip (half-width {W}', trail length RCP32..DF7 +{PAD}' each end), excluding DF2 and DF4:")
print("  ", [(int(ids[k]), note[k] or '-', round(float(re_[k]), 1)) for k in np.where(m)[0]])
n_obs = int(m.sum())
# null A: rotate about NGC 1052
def rotate(P, th):
    c, s = np.cos(th), np.sin(th); return P @ np.array([[c, -s], [s, c]]).T
cnt = []
for _ in range(20000):
    th = rng.uniform(0, 2*np.pi, len(P)); Q = np.column_stack([P[:, 0]*np.cos(th) - P[:, 1]*np.sin(th), P[:, 0]*np.sin(th) + P[:, 1]*np.cos(th)])
    cnt.append(int(strip_count(Q, i2, i4, W, PAD, (i2, i4)).sum()))
cnt = np.array(cnt)
print(f"null A (random position angles about NGC 1052, radii kept): mean {cnt.mean():.2f}, P(>= {n_obs}) = {np.mean(cnt >= n_obs):.4f}")
# null A2: rotate all but DF2/DF4 (line fixed), i.e. are the OTHER objects on the line by chance
cnt2 = []
for _ in range(20000):
    th = rng.uniform(0, 2*np.pi, len(P)); th[[i2, i4]] = 0
    Q = np.column_stack([P[:, 0]*np.cos(th) - P[:, 1]*np.sin(th), P[:, 0]*np.sin(th) + P[:, 1]*np.cos(th)])
    cnt2.append(int(strip_count(Q, i2, i4, W, PAD, (i2, i4)).sum()))
cnt2 = np.array(cnt2)
print(f"null A2 (DF2-DF4 line fixed, others rotated): mean {cnt2.mean():.2f}, P(>= {n_obs}) = {np.mean(cnt2 >= n_obs):.4f}")
# null B: every pair within 1.5 deg of NGC 1052 defines a line; count others on a strip of the same length centred on the pair's midpoint
r = np.hypot(xi, eta); core = np.where(r < 90)[0]
L = abs((P[i32]-P[i2]) @ ((P[i4]-P[i2])/np.linalg.norm(P[i4]-P[i2]))) + abs((P[i7]-P[i2]) @ ((P[i4]-P[i2])/np.linalg.norm(P[i4]-P[i2]))) + 2*PAD
def strip_count_len(P, a, b, w, L, exclude):
    d = P[b] - P[a]; u = d/np.linalg.norm(d); n = np.array([-u[1], u[0]]); mid = (P[a]+P[b])/2
    t = (P - mid) @ u; s = (P - mid) @ n
    m = (np.abs(s) <= w) & (np.abs(t) <= L/2); m[list(exclude)] = False; return m
pair_counts = []
for a, b in itertools.combinations(core, 2):
    if np.linalg.norm(P[a]-P[b]) < 5: continue
    pair_counts.append((int(strip_count_len(P, a, b, W, L, (a, b)).sum()), int(ids[a]), int(ids[b])))
pc = np.array([c[0] for c in pair_counts]); obs_pair = strip_count_len(P, i2, i4, W, L, (i2, i4)).sum()
print(f"null B (lines through every pair of the {len(core)} LSBGs within 1.5 deg of NGC 1052, {len(pc)} pairs, same strip length {L:.0f}'): DF2-DF4 strip holds {obs_pair}; "
      f"pairs with >= {obs_pair}: {np.mean(pc >= obs_pair):.3f}; max {pc.max()}; top pairs: {sorted(pair_counts, reverse=True)[:6]}")
# large LSBGs only
big = re_ > 15
print(f"large LSBGs (r_e > 15\"): {int(big.sum())} objects: ", [(int(ids[k]), note[k] or '-') for k in np.where(big)[0]], "; on the strip:", [(int(ids[k]), note[k] or '-') for k in np.where(m & big)[0]])
cntb = []
for _ in range(20000):
    th = rng.uniform(0, 2*np.pi, len(P)); th[[i2, i4]] = 0
    Q = np.column_stack([P[:, 0]*np.cos(th) - P[:, 1]*np.sin(th), P[:, 0]*np.sin(th) + P[:, 1]*np.cos(th)])
    cntb.append(int((strip_count(Q, i2, i4, W, PAD, (i2, i4)) & big).sum()))
cntb = np.array(cntb); nb = int((m & big).sum())
print(f"  null A2 for large ones: mean {cntb.mean():.2f}, P(>= {nb}) = {np.mean(cntb >= nb):.4f}")
# strip-width sensitivity
for w in (3.0, 5.2, 8.0, 12.0):
    mm = strip_count(P, i2, i4, w, PAD, (i2, i4)); c = []
    for _ in range(5000):
        th = rng.uniform(0, 2*np.pi, len(P)); th[[i2, i4]] = 0
        Q = np.column_stack([P[:, 0]*np.cos(th) - P[:, 1]*np.sin(th), P[:, 0]*np.sin(th) + P[:, 1]*np.cos(th)])
        c.append(int(strip_count(Q, i2, i4, w, PAD, (i2, i4)).sum()))
    c = np.array(c); print(f"  half-width {w:4.1f}': on strip {int(mm.sum())}, null mean {c.mean():.2f}, P = {np.mean(c >= mm.sum()):.4f}")
# null C (proper look-elsewhere): in each rotated catalogue, the richest strip over ALL pair-lines among the core objects,
# compared with the observed richest (which is the trail line, whichever pair defines it)
def max_over_pairs(Q, core, w, L):
    best = 0
    for a, b in itertools.combinations(core, 2):
        if np.linalg.norm(Q[a]-Q[b]) < 5: continue
        c = int(strip_count_len(Q, a, b, w, L, (a, b)).sum()); best = max(best, c)
    return best
for w in (3.0, 5.2, 8.0):
    obs_max = max_over_pairs(P, core, w, L)
    mx = []
    for _ in range(2500):
        th = rng.uniform(0, 2*np.pi, len(P))
        Q = np.column_stack([P[:, 0]*np.cos(th) - P[:, 1]*np.sin(th), P[:, 0]*np.sin(th) + P[:, 1]*np.cos(th)])
        mx.append(max_over_pairs(Q, core, w, L))
    mx = np.array(mx)
    print(f"null C (richest strip over all pair-lines, half-width {w}'): observed max {obs_max}; null max mean {mx.mean():.2f}, P(>= {obs_max}) = {np.mean(mx >= obs_max):.3f}  (2500 rotations)")
