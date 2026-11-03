import re, sys
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import built_on  # the opening "Built on" section (built_on.json)
src=open('site/research-cf4-supernovae.html',encoding='utf-8').read()
head=src[:src.index('<main class="paper">')]
T="ACT DR6 CMB lensing measures the eRASS1 cluster mass scale at 25 to 30 per cent: the pre-registered samples give b<sub>M</sub> = 0.65 to 0.70 &plusmn; 0.19 &plusmn; 0.12, excluding neither the eROSITA weak-lensing calibration nor a halving of it"
Tplain="ACT DR6 CMB lensing measures the eRASS1 cluster mass scale at 25 to 30 per cent: the pre-registered samples give b_M = 0.65 to 0.70 ± 0.19 ± 0.12, excluding neither the eROSITA weak-lensing calibration nor a halving of it"
desc="ACT DR6 CMB-lensing convergence stacked on 2,444 eRASS1 galaxy clusters detects them at 5σ and measures the mass scale of the eROSITA weak-lensing calibration at 25 to 30 per cent: b_M = 0.65 to 0.70 ± 0.19 (statistical) ± 0.12 (model), consistent with the eRASS1 scale and with a 20 to 35 per cent lower one."
head=head.replace("Type Ia supernovae do not see the non-dipole velocity structure of the Cosmicflows-4 reconstructions",Tplain)
head=re.sub(r'<meta name="description" content="[^"]*">',f'<meta name="description" content="{desc}">',head)
head=re.sub(r'<meta property="og:description" content="[^"]*">',f'<meta property="og:description" content="{desc}">',head)
head=head.replace("research-cf4-supernovae.html","research-erass1-act-lensing.html")
# The paper's body HTML, the one argument (its first build's copy is gone; the page has since been edited in place).
body=open(sys.argv[1],encoding='utf-8').read().replace('@@TITLE@@',T)
open('site/research-erass1-act-lensing.html','w',encoding='utf-8').write(built_on.insert(head+body, 'erass1-act-lensing'))
s=open('site/sitemap.xml',encoding='utf-8').read()
if 'research-erass1-act-lensing' not in s:
    s=s.replace('  <url><loc>https://hermesai.dev/research-cf4-supernovae.html</loc></url>\n','  <url><loc>https://hermesai.dev/research-cf4-supernovae.html</loc></url>\n  <url><loc>https://hermesai.dev/research-erass1-act-lensing.html</loc></url>\n')
    open('site/sitemap.xml','w',encoding='utf-8').write(s)
print('written',len(head+body))
