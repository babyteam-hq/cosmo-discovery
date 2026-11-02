"""R-153 Step 0 diagnostic: per-image ON-OFF scatter before and after subtracting the pattern built from the other half
of the same 40 Draco D1 images (results-R153.md: 36 -> 13 nW). Same inputs as alp_diag_pattern.py.
"""
import glob, sys, numpy as np, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from astropy.io import fits
import alp_reduce as R
files=sorted(glob.glob('./data/alp/raw/d1/*.fits'))
BB=[0,1,2,4,5,6,7,9,10,11,15,17,19,21]
ims=[]; geo=[]
for f in files:
    h=fits.open(f); im=h['IMAGE'].data.astype(float); zo=h['ZODI'].data.astype(float); fl=h['FLAGS'].data
    bad=np.zeros(fl.shape,bool)
    for b in BB: bad|=((fl>>b)&1).astype(bool)
    th,cell,lam,bw=R.geometry(h,260.0517,57.9153)
    conv=1e-20*2.998e14/lam*1e9
    r=np.where(bad|~np.isfinite(im),np.nan,(im-zo)*conv)
    lb=np.floor((np.log(lam)-R.L0)/R.DL).astype(int)
    o=np.full(r.shape,np.nan)
    for b in np.unique(lb):
        m=(lb==b)&np.isfinite(r)
        if m.sum()<500: continue
        o[m]=r[m]-np.median(r[m])
    s=1.4826*np.nanmedian(np.abs(o)); clip=np.abs(o)>5*s
    r[clip]=np.nan; o[clip]=np.nan
    ims.append((r,o,th,lb))
n=len(ims)
PA=np.nanmean([ims[i][1] for i in range(0,n,2)],0); PB=np.nanmean([ims[i][1] for i in range(1,n,2)],0)
def onoff(r,th,lb):
    on=(th<0.6)&np.isfinite(r); off=(th>=1.0)&(th<1.75)&np.isfinite(r)
    bins=np.unique(lb[on]); out={}
    for b in bins:
        a=on&(lb==b); c=off&(lb==b)
        if a.sum()>200 and c.sum()>200: out[b]=r[a].mean()-r[c].mean()
    return out
for corr in (False,True):
    allv={}
    for i,(r,o,th,lb) in enumerate(ims):
        rr=r-(PB if i%2==0 else PA) if corr else r
        for b,v in onoff(rr,th,lb).items(): allv.setdefault(b,[]).append(v)
    dev=[np.array(v)-np.mean(v) for v in allv.values() if len(v)>=3]
    dev=np.concatenate(dev)
    print('pattern-corrected' if corr else 'uncorrected', 'per-image ON-OFF per bin scatter about the mean: rms %.2f nW (N=%d)'%(dev.std()*np.sqrt(1), len(dev)))
