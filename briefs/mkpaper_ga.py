"""Builds site/research-giant-arc.html: the revised Giant Arc paper (R-172, from finding-ga-rev.md and results-R155.md F8).

The first version (R-157, "not seen at its published strength") was found on re-check to compare against a
selection-inflated amplitude; its one revision (R-171) measured that shrinkage by injection and passed the finding check
with the claim that the test cannot yet tell a real Arc from chance. This builds that revision as the paper, written as its
interval. The head and styles are the superclusters paper's (site/research-erass1-superclusters.html); the body is
paper_ga_body.html; the figure is a static SVG drawn here from the landed numbers (garev_infer.py's output, results-R155.md
F8). Links to the public repository point at its tag v1.3, the first export holding the revision's files.
Run from the repo root: python research/cosmo-discovery/briefs/mkpaper_ga.py
"""
import re
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import built_on  # the opening "Built on" section (built_on.json)

SITE = 'site/research-giant-arc.html'
LAST = 'site/research-erass1-superclusters.html'
TAG = 'v1.3'
REPO = f'https://github.com/babyteam-hq/cosmo-discovery/blob/{TAG}'
TREE = f'https://github.com/babyteam-hq/cosmo-discovery/tree/{TAG}'
T = ("On new DESI sightlines the Giant Arc cannot yet be told from chance: a real Arc recovered by the discovery recipe "
     "predicts &delta; = +0.62 (68 per cent: +0.25 to +1.05) in its hull, and DESI DR1 sees 8 Mg II absorbers against 8.86 "
     "expected (&delta; = &#8722;0.10), 1.3&sigma; below a real Arc and 0.5&sigma; below chance")
TPLAIN = ("On new DESI sightlines the Giant Arc cannot yet be told from chance: a real Arc recovered by the discovery recipe "
          "predicts δ = +0.62 (68 per cent: +0.25 to +1.05) in its hull, and DESI DR1 sees 8 Mg II absorbers against 8.86 "
          "expected (δ = −0.10), 1.3σ below a real Arc and 0.5σ below chance")
DESC = ("Injecting real overdensities into the SDSS catalogue the Giant Arc was found in, and keeping those the discovery "
        "recipe recovers as Arc-like, a real Arc predicts δ = +0.62 (68 per cent +0.25 to +1.05) on 3,822 DESI DR1 quasar "
        "sightlines SDSS never observed. DESI sees 8 absorbers against 8.86 expected: 1.3σ below a real Arc, 0.5σ below "
        "chance. The data cannot yet tell; DESI DR2 would bring the two readings to about 2σ apart.")

# results-R155.md F8 / finding-ga-rev.md, from garev_infer.py: label, short label, GA-like recoveries, A_true 16/50/84,
# predicted delta 16/50/84, P(N <= 8 | real), its sigma, chance delta (bg under chance), P(N <= 8 | chance), LR real/chance.
ROWS = [
    ('<b>H: uniform in the Arc&rsquo;s hull (primary)</b>', 'uniform in the hull (primary)', 280, (0.5, 1.3, 2.0),
     (0.25, 0.62, 1.05), 0.10, 1.3, 0.15, 0.32, 0.37),
    ('S: member spheres only (clumpy)', 'clumpy: member spheres', 170, (1.6, 3.0, 4.0), (0.15, 0.26, 0.42), 0.22, 0.8, 0.24,
     0.24, 0.91),
    ('L: hull of 60.8 Mpc spheres (diffuse)', 'diffuse: 60.8 Mpc spheres', 200, (0.5, 1.0, 1.3), (0.22, 0.56, 1.05), 0.12, 1.2,
     0.15, 0.32, 0.42),
    ('H, SDSS placement weighted by quasar g', 'hull, g-weighted SDSS', 244, (0.5, 1.0, 1.6), (0.19, 0.48, 0.86), 0.13, 1.1,
     0.14, 0.33, 0.46),
]
OBSERVED = -0.10
PUBLISHED = 1.3


def num(v, places=2):
    """A signed number with a typographic minus."""
    return '0' if round(v, places) == 0 else f'{v:+.{places}f}'.replace('-', '&#8722;')


def figure():
    x0, x1, vmin, vmax, top, rh = 250, 740, -0.4, 1.4, 26, 58
    X = lambda v: x0 + (v - vmin) / (vmax - vmin) * (x1 - x0)
    H = top + len(ROWS) * rh
    out = []
    for v in (-0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25):
        out.append(f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{top}" y2="{H}" style="stroke:var(--line)" stroke-width="1"/>')
        out.append(f'<text x="{X(v):.1f}" y="{H + 18}" text-anchor="middle" font-size="12" style="fill:var(--ink-2)">{num(v)}</text>')
    out.append(f'<line x1="{X(OBSERVED):.1f}" x2="{X(OBSERVED):.1f}" y1="{top}" y2="{H}" style="stroke:var(--ink)" '
               'stroke-width="2"/>')
    out.append(f'<text x="{X(OBSERVED):.1f}" y="{top - 10}" text-anchor="middle" font-size="12" style="fill:var(--ink)">'
               'DESI sees &#8722;0.10</text>')
    out.append(f'<line x1="{X(PUBLISHED):.1f}" x2="{X(PUBLISHED):.1f}" y1="{top}" y2="{H}" style="stroke:var(--ink-2)" '
               'stroke-width="1.5" stroke-dasharray="2 4"/>')
    out.append(f'<text x="{X(PUBLISHED):.1f}" y="{top - 10}" text-anchor="middle" font-size="12" style="fill:var(--ink-2)">'
               'published, 1.3</text>')
    out.append(f'<text x="{(x0 + x1) / 2:.0f}" y="{H + 38}" text-anchor="middle" font-size="12" style="fill:var(--ink-2)">'
               'absorber excess &#948; in the hull on new DESI sightlines</text>')
    for i, (_, short, n, at, dp, p, sig, ch, pch, lr) in enumerate(ROWS):
        y = top + i * rh + rh / 2
        if i:
            out.append(f'<line x1="8" x2="{x1}" y1="{top + i * rh}" y2="{top + i * rh}" style="stroke:var(--line)" stroke-width="1"/>')
        out.append(f'<text x="{x0 - 14}" y="{y - 3:.0f}" text-anchor="end" font-size="13" style="fill:var(--ink)">{short}</text>')
        out.append(f'<text x="{x0 - 14}" y="{y + 13:.0f}" text-anchor="end" font-size="12" style="fill:var(--ink-2)">'
                   f'seen at {sig:.1f}&#963; below a real Arc</text>')
        tip = (f'{short}: {n} Arc-like recoveries; a real Arc predicts &#948; = {num(dp[1])} (68 % {num(dp[0])} to '
               f'{num(dp[2])}); chance {num(ch)}; P(N &#8804; 8 | real) = {p:.2f}, | chance = {pch:.2f}; '
               f'likelihood ratio real/chance {lr:.2f}')
        out.append(f'<g><title>{tip}</title>')
        out.append(f'<line x1="{X(dp[0]):.1f}" x2="{X(dp[2]):.1f}" y1="{y:.1f}" y2="{y:.1f}" style="stroke:var(--series-1)" '
                   'stroke-width="7" stroke-linecap="round" opacity="0.55"/>')
        out.append(f'<circle cx="{X(dp[1]):.1f}" cy="{y:.1f}" r="6" style="fill:var(--series-1);stroke:var(--surface)" '
                   'stroke-width="2"/>')
        out.append(f'<circle cx="{X(ch):.1f}" cy="{y:.1f}" r="5" style="fill:var(--surface);stroke:var(--series-2)" '
                   'stroke-width="2"/>')
        out.append('</g>')
    return (f'<svg viewBox="0 0 760 {H + 46}" role="img" aria-label="For each injected shape, the absorber excess a real '
            'Giant Arc predicts on new DESI sightlines (68 per cent band), the excess chance predicts, and the observed '
            '&#8722;0.10">\n' + '\n'.join(out) + '\n</svg>')


def table():
    rows = []
    for label, _, n, at, dp, p, sig, ch, pch, lr in ROWS:
        rows.append(f'<tr><th scope="row">{label}</th><td>{n}</td><td>{at[0]:.1f} / {at[1]:.1f} / {at[2]:.1f}</td>'
                    f'<td><b>{num(dp[0])} / {num(dp[1])} / {num(dp[2])}</b></td><td>{p:.2f} ({sig:.1f}&sigma;)</td>'
                    f'<td>{num(ch)}</td><td>{pch:.2f}</td><td>{lr:.2f}</td></tr>')
    return '\n'.join(rows)


src = open(LAST, encoding='utf-8').read()
head = src[:src.index('<main class="paper">')]
old_title = re.search(r'<title>([^<]*)</title>', head).group(1)
head = head.replace(old_title, TPLAIN)
head = re.sub(r'<meta name="description" content="[^"]*">', lambda m: f'<meta name="description" content="{DESC.replace("&", "&amp;")}">', head)
head = re.sub(r'<meta property="og:description" content="[^"]*">',
              lambda m: f'<meta property="og:description" content="{DESC.replace("&", "&amp;")}">', head)
head = re.sub(r'<meta property="og:title" content="[^"]*">', lambda m: f'<meta property="og:title" content="{TPLAIN}">', head)
head = head.replace('research-erass1-superclusters.html', 'research-giant-arc.html')
assert 'noindex' not in head
body = open('research/cosmo-discovery/briefs/paper_ga_body.html', encoding='utf-8').read()
body = (body.replace('@@TITLE@@', T).replace('@@FIGURE@@', figure()).replace('@@TABLE@@', table())
        .replace('@@REPO@@', REPO).replace('@@TREE@@', TREE))
assert '@@' not in body
open(SITE, 'w', encoding='utf-8', newline='\n').write(built_on.insert(head + body, 'giant-arc'))
s = open('site/sitemap.xml', encoding='utf-8').read()
if 'research-giant-arc' not in s:
    a = '  <url><loc>https://hermesai.dev/research-erass1-superclusters.html</loc></url>\n'
    assert a in s
    s = s.replace(a, a + '  <url><loc>https://hermesai.dev/research-giant-arc.html</loc></url>\n')
    open('site/sitemap.xml', 'w', encoding='utf-8', newline='\n').write(s)
print('written', SITE, len(head + body))
