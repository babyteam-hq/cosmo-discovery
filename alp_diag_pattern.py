"""R-153 Step 0 diagnostic: the fixed detector pattern (half-split of 40 Draco D1 images, detector coordinates).

Fetch first: awk -F, '$2==1{print $4}' data/alp/list_draco.csv | head -40 | (cd data/alp/raw/d1 && xargs -P 8 -n 1 curl -sS -O)
Prints the half-stack pattern rms and correlation (results-R153.md: 35 nW, 0.98). Run from research/cosmo-discovery.
"""
import glob, numpy as np, warnings
warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
import logging; logging.getLogger('astropy').setLevel(logging.ERROR)
files=sorted(glob.glob('./data/alp/raw/d1/*.fits'))
BB=[0,1,2,4,5,6,7,9,10,11,15,17,19,21]
stk=[]; 
for f in files:
    h=fits.open(f); im=h['IMAGE'].data.astype(float); zo=h['ZODI'].data.astype(float); fl=h['FLAGS'].data
    bad=np.zeros(fl.shape,bool)
    for b in BB: bad|=((fl>>b)&1).astype(bool)
    ww=WCS(h['IMAGE'].header,fobj=h,key='W'); g=np.arange(0,2048,8).clip(0,2039)
    yc,xc=np.meshgrid(g,g,indexing='ij'); lam,_=ww.pixel_to_world_values(xc.ravel(),yc.ravel()); lam=lam.reshape(xc.shape)
    lamf=np.repeat(np.repeat(lam,8,0),8,1)[:2040,:2040]
    r=np.where(bad|~np.isfinite(im),np.nan,(im-zo))
    lb=np.floor(np.log(lamf)/0.002).astype(int)
    # subtract the per-lambda-bin median (the exposure's sky spectrum incl. offset): leaves detector pattern + sky structure
    out=np.full(r.shape,np.nan)
    for b in np.unique(lb):
        m=(lb==b)&np.isfinite(r)
        if m.sum()<500: continue
        out[m]=r[m]-np.median(r[m])
    conv=1e-20*2.998e14/lamf*1e9
    # clip sources residuals
    o=out*conv; s=1.4826*np.nanmedian(np.abs(o)); o[np.abs(o)>5*s]=np.nan
    stk.append(o.reshape(255,8,255,8))
    print(f.split('/')[-1][:30], 'rms %.1f'%np.nanstd(o), flush=True) if len(stk)<3 else None
S=np.array(stk)  # (n,255,8,255,8)
bm=np.nanmean(S,axis=(2,4))  # per-image 8x8 block means (n,255,255)
A=np.nanmean(bm[0::2],0); B=np.nanmean(bm[1::2],0)
m=np.isfinite(A)&np.isfinite(B)
print('half-stack block(8px) pattern: rms A %.2f B %.2f, corr %.2f, rms(A-B)/2 %.2f nW'%(A[m].std(),B[m].std(),np.corrcoef(A[m],B[m])[0,1],(A-B)[m].std()/2))
print('single image 8px-block rms %.2f nW'%np.nanstd(bm[0]))
# along-row profile: rows structure
ra=np.nanmean(A,1); rb=np.nanmean(B,1); mm=np.isfinite(ra)&np.isfinite(rb); print('row-profile corr %.2f, rms %.2f'%(np.corrcoef(ra[mm],rb[mm])[0,1], ra[mm].std()))

