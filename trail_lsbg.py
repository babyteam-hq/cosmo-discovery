"""NGC 1052 trail, sweep-2.md row 34 (R-117 did this and lost it): do two further independent LSBG catalogues hold more
objects on the extended DF2-DF4 line, beyond RCP32 (east) and DF7 (west)? Catalogues (VizieR J/ApJS/252/18, Tanoglidis+2021,
the DES Y3 LSBG catalogue, 20,977 objects, r_e > 2.5", mu_eff,g > 24.3; and J/ApJS/267/27 table2, SMUDGes, Zaritsky+2023,
r_e > 5.3", mu_0,g > 24), fetched over a 4 x 3 deg box about (40.1, -8.2) from vizier.idia.ac.za into data/lsbg/.
 (1) which of the 12 trail objects each catalogue holds (match within 0.5');
 (2) objects within +-5.2' of the DF2-DF4 line, by position t along it (t = 0 at the DF2-DF4 midpoint, +toward DF2/RCP32):
     the trail's span (DF7 at t = -36', RCP32 at +30') and the extensions to +-120';
 (3) the count on each extension against the radius-preserving rotation null (van Dokkum+2022's: rotate the catalogue about
     NGC 1052, keep the line), and the DF-class subset (r_e > 15", mu_0 or mu_eff > 24).
usage: python trail_lsbg.py data/lsbg
"""
import sys, numpy as np

N1052 = (40.2700, -8.2558)
TRAIL = {'RCP32': (40.6202, -8.3768), 'DF2': (40.4451, -8.4028), 'RCP28': (40.4215, -8.3475), 'RCP26': (40.2897, -8.2968), 'RCP21': (40.1200, -8.2434),
         'DF9': (40.0292, -8.2290), 'RCP17': (39.9696, -8.2121), 'Ta21-12000': (39.9139, -8.2285), 'DF4': (39.8128, -8.1160), 'DF5': (39.8028, -8.1408),
         'LEDA4014647': (39.7021, -8.0494), 'DF7': (39.6241, -7.9257)}
W = 5.2; EXT = 120.0

def tangent(ra, dec, ra0=N1052[0], dec0=N1052[1]):
    ra, dec, ra0, dec0 = map(np.radians, (ra, dec, ra0, dec0))
    cosc = np.sin(dec0)*np.sin(dec) + np.cos(dec0)*np.cos(dec)*np.cos(ra - ra0)
    xi = np.cos(dec)*np.sin(ra - ra0)/cosc; eta = (np.cos(dec0)*np.sin(dec) - np.sin(dec0)*np.cos(dec)*np.cos(ra - ra0))/cosc
    return np.degrees(xi)*60, np.degrees(eta)*60

def read_tsv(path, want):
    """VizieR asu-tsv: several blocks (tables); return the block whose header holds every name in `want`, as dict of columns."""
    blocks, cur = [], []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            if cur: blocks.append(cur); cur = []
            continue
        cur.append(line.rstrip('\n').split('\t'))
    if cur: blocks.append(cur)
    for b in blocks:
        hdr = [h.strip() for h in b[0]]
        if all(w in hdr for w in want):
            rows = [r for r in b[3:] if len(r) == len(hdr)]
            return {h: np.array([r[i].strip() for r in rows]) for i, h in enumerate(hdr)}
    raise SystemExit(f'{path}: no table with {want}')

def num(a):
    out = np.full(len(a), np.nan)
    for i, v in enumerate(a):
        try: out[i] = float(v)
        except ValueError: pass
    return out

def main(ddir):
    cats = {}
    d = read_tsv(f'{ddir}/J_ApJS_252_18.tsv', ['RAJ2000', 'DEJ2000', 'gReff'])
    cats['DES Y3 LSBG (Tanoglidis+2021)'] = dict(ra=num(d['RAJ2000']), dec=num(d['DEJ2000']), re=num(d['gReff']), mu=num(d['gmu0']), name=d['ID'])  # gReff arcsec (Galfit), gmu0 central SB
    s = read_tsv(f'{ddir}/J_ApJS_267_27_table2.tsv', ['RAJ2000', 'DEJ2000', 'Re'])
    cats['SMUDGes (Zaritsky+2023)'] = dict(ra=num(s['RAJ2000']), dec=num(s['DEJ2000']), re=num(s['Re']), mu=num(s['mu0g']), name=s['SMDG'])
    tx, ty = tangent(np.array([v[0] for v in TRAIL.values()]), np.array([v[1] for v in TRAIL.values()])); tn = list(TRAIL)
    a, b = np.array([tx[tn.index('DF2')], ty[tn.index('DF2')]]), np.array([tx[tn.index('DF4')], ty[tn.index('DF4')]])
    mid = (a + b)/2; u = (a - b)/np.linalg.norm(a - b); n = np.array([-u[1], u[0]])
    tt, ss = (np.column_stack([tx, ty]) - mid)@u, (np.column_stack([tx, ty]) - mid)@n
    print('trail objects along the DF2-DF4 line (t along, s across, arcmin): ' + ', '.join(f'{k} ({t:+.1f}, {sv:+.1f})' for k, t, sv in zip(tn, tt, ss)))
    tmin, tmax = tt.min(), tt.max(); rng = np.random.default_rng(34)
    for cname, c in cats.items():
        ok = np.isfinite(c['ra']) & np.isfinite(c['dec']); x, y = tangent(c['ra'][ok], c['dec'][ok]); re_, mu, nm = c['re'][ok], c['mu'][ok], c['name'][ok]
        P = np.column_stack([x, y]); t, sp = (P - mid)@u, (P - mid)@n; r = np.hypot(x, y)
        # (1) which trail objects it holds
        held = []
        for k, (xk, yk) in zip(tn, zip(tx, ty)):
            dd = np.hypot(x - xk, y - yk); j = np.argmin(dd)
            if dd[j] < 0.5: held.append(f'{k} (r_e {re_[j]:.1f}", mu {mu[j]:.1f})')
        print(f'\n{cname}: {ok.sum()} objects in the box (r_e 5-95 %: {np.nanpercentile(re_, 5):.1f}-{np.nanpercentile(re_, 95):.1f}"); holds {len(held)} of the 12 trail objects: ' + '; '.join(held))
        strip = np.abs(sp) < W; big = (re_ > 15) & (mu > 24)
        # (2) objects on the strip, listed by t
        for j in np.where(strip & (np.abs(t) < EXT))[0]:
            where = 'trail span' if tmin - 1 < t[j] < tmax + 1 else ('east extension (beyond RCP32)' if t[j] > 0 else 'west extension (beyond DF7)')
            known = [k for k, xk, yk in zip(tn, tx, ty) if np.hypot(x[j] - xk, y[j] - yk) < 0.5]
            print(f'  t={t[j]:+6.1f}\' s={sp[j]:+5.1f}\' RA={c["ra"][ok][j]:.4f} Dec={c["dec"][ok][j]:.4f} r_e={re_[j]:5.1f}" mu={mu[j]:5.1f} {nm[j]} [{where}]' + (f' = {known[0]}' if known else (' NEW, DF-class' if big[j] else ' new, small')))
        # (3) counts on the extensions vs the rotation null, all and DF-class
        for label, sel in (('all', np.ones_like(strip)), ('DF-class (r_e > 15", mu > 24)', big)):
            east = strip & sel & (t > tmax + 1) & (t < EXT); west = strip & sel & (t < tmin - 1) & (t > -EXT)
            ne, nw = east.sum(), west.sum(); nulls = np.empty((2000, 2))
            for i in range(2000):
                th = rng.uniform(0, 2*np.pi, len(r)); Q = np.column_stack([r*np.cos(th), r*np.sin(th)]); tq, sq = (Q - mid)@u, (Q - mid)@n
                st = (np.abs(sq) < W) & sel; nulls[i] = [(st & (tq > tmax + 1) & (tq < EXT)).sum(), (st & (tq < tmin - 1) & (tq > -EXT)).sum()]
            print(f'  {label}: east extension {ne} (null {nulls[:, 0].mean():.2f}, P(>=) {np.mean(nulls[:, 0] >= ne):.3f}); west extension {nw} (null {nulls[:, 1].mean():.2f}, P(>=) {np.mean(nulls[:, 1] >= nw):.3f}); both {ne + nw} (null {nulls.sum(1).mean():.2f}, P(>=) {np.mean(nulls.sum(1) >= ne + nw):.3f})')
            # and the trail span itself, excluding the known trail objects, as the same catalogue's view of the alignment
            span = strip & sel & (t > tmin - 1) & (t < tmax + 1)
            known = np.array([any(np.hypot(x[j] - xk, y[j] - yk) < 0.5 for xk, yk in zip(tx, ty)) for j in range(len(x))])
            print(f'  {label}: on the trail span {span.sum()}, of which known trail objects {np.sum(span & known)}, other {np.sum(span & ~known)}')

main(sys.argv[1])
