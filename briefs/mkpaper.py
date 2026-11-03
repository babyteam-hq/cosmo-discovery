import re, sys
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import built_on  # the opening "Built on" section (built_on.json)
# The head and styles are the ACT page's own (they were the CF4 paper's, which is now a withdrawal notice, R-169).
src=open('site/research-erass1-act-lensing.html',encoding='utf-8').read()
head=src[:src.index('<main class="paper">')]
T="ACT DR6 CMB lensing detects eRASS1 clusters at 5&sigma;; their mass rescaling, b<sub>M</sub> = 0.65 to 0.70 &plusmn; 0.22 for the SZ-faint clusters, carries an unmodelled selection bias and does not yet test the eROSITA weak-lensing calibration"
Tplain="ACT DR6 CMB lensing detects eRASS1 clusters at 5σ; their mass rescaling, b_M = 0.65 to 0.70 ± 0.22 for the SZ-faint clusters, carries an unmodelled selection bias and does not yet test the eROSITA weak-lensing calibration"
desc="ACT DR6 CMB-lensing convergence stacked on 2,444 eRASS1 galaxy clusters outside the SZ-cluster holes detects them at 5σ, with a mass rescaling b_M = 0.65 to 0.70, ± 0.17 to 0.19 (statistical) ± 0.12 (model). Removing the SZ-detected clusters biases b_M by an amount the model does not include, and the estimator's response inside the holes is uncalibrated, so this is a detection, not yet a test of the eRASS1 mass scale."
head=re.sub(r'<meta name="description" content="[^"]*">',f'<meta name="description" content="{desc}">',head)
head=re.sub(r'<meta property="og:description" content="[^"]*">',f'<meta property="og:description" content="{desc}">',head)
# The paper's body HTML, the one argument (its first build's copy is gone; the page has since been edited in place).
body=open(sys.argv[1],encoding='utf-8').read().replace('@@TITLE@@',T)
open('site/research-erass1-act-lensing.html','w',encoding='utf-8').write(built_on.insert(head+body, 'erass1-act-lensing'))
s=open('site/sitemap.xml',encoding='utf-8').read()
if 'research-erass1-act-lensing' not in s:
    s=s.replace('</urlset>','  <url><loc>https://hermesai.dev/research-erass1-act-lensing.html</loc></url>\n</urlset>')
    open('site/sitemap.xml','w',encoding='utf-8').write(s)
print('written',len(head+body))
