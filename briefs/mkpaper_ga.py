"""Builds site/research-giant-arc.html (R-157) from the numbers in results-R155.md (F5, F6) and finding.md's R-155 section.

The head and styles are the last paper's (site/research-erass1-superclusters.html); the body is paper_ga_body.html; the figure
is a static SVG drawn here from the landed numbers. Run from the repo root: python research/cosmo-discovery/briefs/mkpaper_ga.py
"""
import re
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import built_on  # the opening "Built on" section (built_on.json)

SITE = 'site/research-giant-arc.html'
LAST = 'site/research-erass1-superclusters.html'
T = ("The Giant Arc does not reappear on independent sightlines: on 3,822 DESI DR1 quasar sightlines that SDSS never observed, "
     "its fitted Mg II overdensity is A = &#8722;0.23 (95 per cent: &#8722;0.90 to +0.63), excluding the published "
     "&delta;&rho;/&rho; = 1.3 at 3.2&sigma;")
TPLAIN = ("The Giant Arc does not reappear on independent sightlines: on 3,822 DESI DR1 quasar sightlines that SDSS never "
          "observed, its fitted Mg II overdensity is A = −0.23 (95 per cent: −0.90 to +0.63), excluding the published "
          "δρ/ρ = 1.3 at 3.2σ")
DESC = ("Counting DESI DR1 Mg II absorbers on 3,822 quasar sightlines through the Giant Arc's frozen hull that are in neither "
        "SDSS DR16Q nor the Zhu & Ménard searched lists: 8 observed against 8.86 expected. The fitted overdensity is "
        "A = −0.23 (68 per cent −0.71 to +0.26), excluding the published δρ/ρ = 1.3 at 3.2σ; a weaker real structure "
        "(A ≲ 0.6) is not excluded.")

# Per cell: label, short label, sightlines, N_obs, N_exp, delta, shift mean (sd), Delta [95 %], A, A 68 %, A 95 %, sigma vs 1.3.
# results-R155.md F6; A's 68 % for the non-primary cells is Delta's 68 % over f_H (F5: 0.976, 0.878, 0.730, 0.515).
ROWS = [
    ('<b>CHMS (primary)</b>', 'CHMS (primary)', 3822, 8, 8.86, -0.097, (-0.100, 0.341), (-0.88, +0.62), -0.23, (-0.71, +0.26), (-0.90, +0.63), 3.2),
    ('spheres R = 47.5 Mpc', 'R = 47.5 Mpc', 3533, 4, 4.13, -0.032, (-0.085, 0.467), (-0.97, +0.98), -0.02, (-0.59, +0.57), (-1.10, +1.11), 2.3),
    ('spheres R = 70 Mpc', 'R = 70 Mpc', 4806, 6, 9.64, -0.378, (-0.113, 0.336), (-0.91, +0.41), -0.45, (-1.00, +0.07), (-1.25, +0.56), 3.3),
    ('spheres R = 95 Mpc', 'R = 95 Mpc', 6281, 15, 16.39, -0.085, (-0.111, 0.265), (-0.43, +0.63), +0.05, (-0.54, +0.60), (-0.84, +1.23), 2.2),
    ('CHMS, W &gt; 0.6 &Aring;', 'CHMS, W &gt; 0.6 &#197;', 3822, 8, 8.40, -0.048, (-0.068, 0.337), (-0.77, +0.69), -0.10, (-0.60, +0.37), (-0.79, +0.70), 2.9),
]
PUBLISHED = 1.3


def num(v, places=2):
    """A signed number with a typographic minus."""
    return '0' if round(v, places) == 0 else f'{v:+.{places}f}'.replace('-', '&#8722;')


def figure():
    x0, x1, vmin, vmax, top, rh = 250, 740, -1.4, 1.6, 26, 58
    X = lambda v: x0 + (v - vmin) / (vmax - vmin) * (x1 - x0)
    H = top + len(ROWS) * rh
    out = []
    for v in (-1.0, -0.5, 0.0, 0.5, 1.0, 1.5):
        out.append(f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{top}" y2="{H}" style="stroke:var(--line)" stroke-width="1"/>')
        out.append(f'<text x="{X(v):.1f}" y="{H + 18}" text-anchor="middle" font-size="12" style="fill:var(--ink-2)">{num(v, 1)}</text>')
    out.append(f'<line x1="{X(0):.1f}" x2="{X(0):.1f}" y1="{top}" y2="{H}" style="stroke:var(--ink-2)" stroke-width="1.5" '
               'stroke-dasharray="2 4"/>')
    out.append(f'<text x="{X(0):.1f}" y="{top - 10}" text-anchor="middle" font-size="12" style="fill:var(--ink-2)">chance, A = 0</text>')
    out.append(f'<line x1="{X(PUBLISHED):.1f}" x2="{X(PUBLISHED):.1f}" y1="{top}" y2="{H}" style="stroke:var(--series-2)" '
               'stroke-width="2" stroke-dasharray="7 4"/>')
    out.append(f'<text x="{X(PUBLISHED):.1f}" y="{top - 10}" text-anchor="middle" font-size="12" style="fill:var(--series-2)">'
               'published, 1.3</text>')
    out.append(f'<text x="{(x0 + x1) / 2:.0f}" y="{H + 38}" text-anchor="middle" font-size="12" style="fill:var(--ink-2)">'
               'fitted overdensity A in the hull (&#948;&#961;/&#961;)</text>')
    for i, (_, short, n, nobs, nexp, d, sh, d95, a, a68, a95, sig) in enumerate(ROWS):
        y = top + i * rh + rh / 2
        if i:
            out.append(f'<line x1="8" x2="{x1}" y1="{top + i * rh}" y2="{top + i * rh}" style="stroke:var(--line)" stroke-width="1"/>')
        out.append(f'<text x="{x0 - 14}" y="{y - 3:.0f}" text-anchor="end" font-size="13" style="fill:var(--ink)">{short}</text>')
        out.append(f'<text x="{x0 - 14}" y="{y + 13:.0f}" text-anchor="end" font-size="12" style="fill:var(--ink-2)">'
                   f'{nobs} seen, {nexp:.2f} expected; 1.3 at {sig:.1f}&#963;</text>')
        tip = (f'{short}: {n:,} sightlines, {nobs} absorbers against {nexp:.2f} expected; A = {num(a)}, '
               f'68 % [{num(a68[0])}, {num(a68[1])}], 95 % [{num(a95[0])}, {num(a95[1])}]; published 1.3 excluded at {sig:.1f}&#963;')
        out.append(f'<g><title>{tip}</title>')
        out.append(f'<line x1="{X(a95[0]):.1f}" x2="{X(a95[1]):.1f}" y1="{y:.1f}" y2="{y:.1f}" style="stroke:var(--series-1)" '
                   'stroke-width="1.5" stroke-linecap="round"/>')
        out.append(f'<line x1="{X(a68[0]):.1f}" x2="{X(a68[1]):.1f}" y1="{y:.1f}" y2="{y:.1f}" style="stroke:var(--series-1)" '
                   'stroke-width="7" stroke-linecap="round" opacity="0.55"/>')
        out.append(f'<circle cx="{X(a):.1f}" cy="{y:.1f}" r="6" style="fill:var(--series-1);stroke:var(--surface)" '
                   'stroke-width="2"/>')
        out.append('</g>')
    return (f'<svg viewBox="0 0 760 {H + 46}" role="img" aria-label="The fitted Mg II overdensity A of the Giant Arc per cell, '
            'with 68 and 95 per cent intervals, against zero and the published 1.3">\n' + '\n'.join(out) + '\n</svg>')


def table():
    rows = []
    for label, _, n, nobs, nexp, d, sh, d95, a, a68, a95, sig in ROWS:
        rows.append(f'<tr><th scope="row">{label}</th><td>{n:,}</td><td>{nobs}</td><td>{nexp:.2f}</td><td>{num(d)}</td>'
                    f'<td>{num(sh[0])} ({sh[1]:.2f})</td><td>[{num(d95[0])}, {num(d95[1])}]</td>'
                    f'<td><b>{num(a)}</b> [{num(a95[0])}, {num(a95[1])}]</td><td>{sig:.1f}&sigma;</td></tr>')
    return '\n'.join(rows)


src = open(LAST, encoding='utf-8').read()
head = src[:src.index('<main class="paper">')]
old_title = re.search(r'<title>([^<]*)</title>', head).group(1)
head = head.replace(old_title, TPLAIN)
head = re.sub(r'<meta name="description" content="[^"]*">', lambda m: f'<meta name="description" content="{DESC.replace("&", "&amp;")}">', head)
head = re.sub(r'<meta property="og:description" content="[^"]*">',
              lambda m: f'<meta property="og:description" content="{DESC.replace("&", "&amp;")}">', head)
head = head.replace('research-erass1-superclusters.html', 'research-giant-arc.html')
body = open('research/cosmo-discovery/briefs/paper_ga_body.html', encoding='utf-8').read()
body = body.replace('@@TITLE@@', T).replace('@@FIGURE@@', figure()).replace('@@TABLE@@', table())
open(SITE, 'w', encoding='utf-8', newline='\n').write(built_on.insert(head + body, 'giant-arc'))
s = open('site/sitemap.xml', encoding='utf-8').read()
if 'research-giant-arc' not in s:
    a = '  <url><loc>https://hermesai.dev/research-erass1-superclusters.html</loc></url>\n'
    assert a in s
    s = s.replace(a, a + '  <url><loc>https://hermesai.dev/research-giant-arc.html</loc></url>\n')
    open('site/sitemap.xml', 'w', encoding='utf-8', newline='\n').write(s)
print('written', SITE, len(head + body))
