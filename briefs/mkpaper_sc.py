"""Builds site/research-erass1-superclusters.html (R-142) from the numbers in results-R140.md (D4, D5) and results-R137.md (D4).

The head and styles are the last paper's (site/research-erass1-act-lensing.html); the body is below; the figure is a static
SVG drawn here from the landed numbers. Run from the repo root: python research/cosmo-discovery/briefs/mkpaper_sc.py
"""
import re
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import built_on  # the opening "Built on" section (built_on.json)

SITE = 'site/research-erass1-superclusters.html'
T = ("ACT DR6 CMB lensing over 365 eRASS1 superclusters, against the &Lambda;CDM expectation from the Websky light-cone: no "
     "lensing excess beyond &Lambda;CDM outside the members' one- and two-halo model, with 95 per cent upper limits of 0.9 and "
     "0.3 of the member model")
TPLAIN = ("ACT DR6 CMB lensing over 365 eRASS1 superclusters, against the ΛCDM expectation from the Websky light-cone: no "
          "lensing excess beyond ΛCDM outside the members' one- and two-halo model, with 95 per cent upper limits of 0.9 and "
          "0.3 of the member model")
DESC = ("The ACT DR6 convergence over 365 eRASS1 superclusters, beyond their member clusters' one- and two-halo model and "
        "against matched isolated clusters, agrees with the ΛCDM expectation from the Websky light-cone within 1.1σ in every "
        "pre-registered row: relative to ΛCDM the excess is +0.18 ± 0.37 of the member model in the members' hull and "
        "−0.43 ± 0.38 in a 20-Mpc disc.")

# Per sample and aperture: label, n, observed E - E_control (R-137 D4), LCDM mean, one-sample sd, LCDM 95 % range (R-140 D4),
# Delta and its sigma, 95 % interval on Delta, percentile (R-140 D5).
ROWS = [
    ('&ge; 2 members, hull + 5 cMpc', '&ge; 2 members, hull', 365, +0.29, +0.107, 0.363, (-0.59, +0.83), +0.18, 0.37, (-0.55, +0.91), 69),
    ('&ge; 2 members, disc 20 cMpc', '&ge; 2 members, disc', 365, -0.31, +0.120, 0.368, (-0.56, +0.87), -0.43, 0.38, (-1.17, +0.31), 12),
    ('&ge; 3 members, hull + 5 cMpc', '&ge; 3 members, hull', 119, +0.27, +0.129, 0.495, (-0.81, +1.12), +0.14, 0.51, (-0.86, +1.14), 61),
    ('&ge; 3 members, disc 20 cMpc', '&ge; 3 members, disc', 119, -0.05, +0.204, 0.603, (-0.96, +1.35), -0.25, 0.63, (-1.48, +0.97), 35),
]


def num(v, places=2):
    """A signed number with a typographic minus."""
    return '0' if round(v, places) == 0 else f'{v:+.{places}f}'.replace('-', '&#8722;')


def figure():
    x0, x1, vmin, vmax, top, rh = 250, 740, -1.25, 1.5, 12, 58
    X = lambda v: x0 + (v - vmin) / (vmax - vmin) * (x1 - x0)
    H = top + len(ROWS) * rh
    out = []
    for v in (-1.0, -0.5, 0.0, 0.5, 1.0, 1.5):
        out.append(f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{top}" y2="{H}" style="stroke:var(--line)" stroke-width="1"/>')
        out.append(f'<text x="{X(v):.1f}" y="{H + 18}" text-anchor="middle" font-size="12" style="fill:var(--ink-2)">{num(v, 1)}</text>')
    out.append(f'<line x1="{X(0):.1f}" x2="{X(0):.1f}" y1="{top}" y2="{H}" style="stroke:var(--ink-2)" stroke-width="1.5" '
               'stroke-dasharray="3 4"/>')
    out.append(f'<text x="{(x0 + x1) / 2:.0f}" y="{H + 38}" text-anchor="middle" font-size="12" style="fill:var(--ink-2)">'
               'E &#8722; E_control, in units of the member model</text>')
    for i, (_, short, n, obs, mu, sd, r95, d, ds, d95, pct) in enumerate(ROWS):
        y = top + i * rh + rh / 2
        if i:
            out.append(f'<line x1="8" x2="{x1}" y1="{top + i * rh}" y2="{top + i * rh}" style="stroke:var(--line)" stroke-width="1"/>')
        out.append(f'<text x="{x0 - 14}" y="{y - 3:.0f}" text-anchor="end" font-size="13" style="fill:var(--ink)">{short}</text>')
        out.append(f'<text x="{x0 - 14}" y="{y + 13:.0f}" text-anchor="end" font-size="12" style="fill:var(--ink-2)">'
                   f'n = {n}; &#916; = {num(d)} &#177; {ds:.2f}</text>')
        tip = (f'{short}: &#923;CDM {num(mu)}, one-sample sd {sd:.2f}, 95 % [{num(r95[0])}, {num(r95[1])}]; '
               f'observed {num(obs)} (percentile {pct})')
        out.append(f'<g><title>{tip}</title>')
        out.append(f'<line x1="{X(r95[0]):.1f}" x2="{X(r95[1]):.1f}" y1="{y:.1f}" y2="{y:.1f}" style="stroke:var(--muted)" '
                   'stroke-width="1.5" stroke-linecap="round"/>')
        out.append(f'<line x1="{X(mu - sd):.1f}" x2="{X(mu + sd):.1f}" y1="{y:.1f}" y2="{y:.1f}" style="stroke:var(--muted)" '
                   'stroke-width="7" stroke-linecap="round" opacity="0.55"/>')
        out.append(f'<line x1="{X(mu):.1f}" x2="{X(mu):.1f}" y1="{y - 8:.1f}" y2="{y + 8:.1f}" style="stroke:var(--ink-2)" '
                   'stroke-width="2"/>')
        out.append(f'<circle cx="{X(obs):.1f}" cy="{y:.1f}" r="6" style="fill:var(--series-1);stroke:var(--surface)" '
                   'stroke-width="2"/>')
        out.append('</g>')
    return (f'<svg viewBox="0 0 760 {H + 46}" role="img" aria-label="Observed E minus E_control per sample and aperture '
            'against the LCDM distribution for a data-sized sample">\n' + '\n'.join(out) + '\n</svg>')


def table():
    rows = []
    for label, _, n, obs, mu, sd, r95, d, ds, d95, pct in ROWS:
        rows.append(f'<tr><th scope="row">{label}</th><td>{n}</td><td>{num(obs)}</td><td>{num(mu)}; {sd:.2f}</td>'
                    f'<td><b>{num(d)} &plusmn; {ds:.2f}</b></td><td>[{num(d95[0])}, {num(d95[1])}]</td><td>{pct}</td></tr>')
    return '\n'.join(rows)


src = open('site/research-erass1-act-lensing.html', encoding='utf-8').read()
head = src[:src.index('<main class="paper">')]
old_title = re.search(r'<title>([^<]*)</title>', head).group(1)
head = head.replace(old_title, TPLAIN)
head = re.sub(r'<meta name="description" content="[^"]*">', f'<meta name="description" content="{DESC}">', head)
head = re.sub(r'<meta property="og:description" content="[^"]*">', f'<meta property="og:description" content="{DESC}">', head)
head = head.replace('research-erass1-act-lensing.html', 'research-erass1-superclusters.html')
head = head.replace('.legend span::before{', '.legend .band::before{width:22px;height:7px;border-radius:4px;opacity:.55}\n.legend span::before{')
body = open('research/cosmo-discovery/briefs/paper_sc_body.html', encoding='utf-8').read()
body = body.replace('@@TITLE@@', T).replace('@@FIGURE@@', figure()).replace('@@TABLE@@', table())
open(SITE, 'w', encoding='utf-8', newline='\n').write(built_on.insert(head + body, 'erass1-superclusters'))
s = open('site/sitemap.xml', encoding='utf-8').read()
if 'research-erass1-superclusters' not in s:
    a = '  <url><loc>https://hermesai.dev/research-erass1-act-lensing.html</loc></url>\n'
    s = s.replace(a, a + '  <url><loc>https://hermesai.dev/research-erass1-superclusters.html</loc></url>\n')
    open('site/sitemap.xml', 'w', encoding='utf-8', newline='\n').write(s)
print('written', SITE, len(head + body))
