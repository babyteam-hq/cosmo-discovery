"""R-153 Step 0 diagnostic: the noise of the cell entries before the pattern (data/alp/red/draco, reduced without
--pattern): chi2 by pixel count and detector, the coherence of ON cells within one (exposure, bin), and the tails
(results-R153.md: 1.2% of entries carry 72% of chi2).
"""
import sys; sys.path.insert(0,'.')
import numpy as np, alp_cells as C
d=C.load('draco'); ent=C.offsub(d); sel=ent['theta']<0.6
ent={k:(v[sel] if isinstance(v,np.ndarray) else v) for k,v in ent.items()}
deg=4; w=1/ent['var']; res=np.full(len(w),np.nan)
order=np.argsort(ent['cell'],kind='stable'); cells,st=np.unique(ent['cell'][order],return_index=True); en=np.r_[st[1:],len(order)]
for c,a,b in zip(cells,st,en):
    idx=order[a:b]
    if len(idx)<20: continue
    B=C.basis(ent['det'][idx],ent['lam'][idx],deg); ww=w[idx]
    F=np.linalg.pinv(B.T@(B*ww[:,None])); y=ent['y'][idx]; res[idx]=y-B@(F@(B.T@(ww*y)))
ok=np.isfinite(res); z=res*np.sqrt(w)
print('chi2 per entry overall', np.mean(z[ok]**2))
for lo,hi in [(4,10),(10,30),(30,100),(100,300),(300,3000)]:
    m=ok&(ent['n']>=lo)&(ent['n']<hi); print('n %d-%d: entries %d, <z^2> %.1f, rms res nW %.3f, rms pred %.3f'%(lo,hi,m.sum(),np.mean(z[m]**2),np.std(res[m]),np.sqrt(np.mean(ent['var'][m]))))
for k in range(1,7):
    m=ok&(ent['det']==k); print('det',k,'<z^2> %.1f'%np.mean(z[m]**2))
# coherence: per (e,lb) mean residual across cells vs expectation
key=ent['e'].astype(np.int64)*958+ent['lb']; uk,inv=np.unique(key[ok],return_inverse=True)
cnt=np.bincount(inv); sm=np.bincount(inv,weights=res[ok]); s2=np.bincount(inv,weights=res[ok]**2)
mm=cnt>=5
# ratio of variance of the group mean to (mean of variance)/cnt  -> 1 if independent
vm=(sm[mm]/cnt[mm])**2; vi=s2[mm]/cnt[mm]/cnt[mm]
print('coherence across cells in same (e,lb): var(mean)/(var/n) = %.2f'%(vm.mean()/vi.mean()))
# radial structure: mean residual per exposure by theta bins
for lo,hi in [(0,0.15),(0.15,0.3),(0.3,0.45),(0.45,0.6)]:
    m=ok&(ent['theta']>=lo)&(ent['theta']<hi); print('theta %.2f-%.2f mean res %.3f nW, rms %.3f'%(lo,hi,np.mean(res[m]),np.std(res[m])))
np.save('data/alp/res_draco.npy',res)
