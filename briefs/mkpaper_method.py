"""Builds site/research-methodology.html (R-162): how the cosmology papers are made, with the referee prompts verbatim.

The head and styles are the Giant Arc paper's (site/research-giant-arc.html); the body is paper_method_body.html; the three
prompts are read from blind-review.md here, so the page never drifts from the rules. Links to the public repository point at
the tag v1. Run from the repo root: python research/cosmo-discovery/briefs/mkpaper_method.py
"""
import html
import re

SITE = 'site/research-methodology.html'
LAST = 'site/research-giant-arc.html'
REPO = 'https://github.com/babyteam-hq/cosmo-discovery/tree/v1'
T = 'How the HermesAI cosmology papers are made: the blind gate, the open-problem map and the record of passes and misses'
DESC = ('Methods behind the HermesAI cosmology papers: AI agents choose, run and write up tests on public data, and a '
        'context-free referee judges each one three times (trajectory, method, finding) with fixed verdict lines. The prompts '
        'verbatim, how candidates are found and dispatched, and the record of passes, misses and gate stops.')
PROMPTS = {'PROMPT_TRAJECTORY': 'trajectory check', 'PROMPT_METHOD': 'method check', 'PROMPT_FINDING': 'finding check'}


def prompt(rules, name):
    """The indented block under '## Prompt (<name>...)' in blind-review.md, dedented by its four spaces."""
    m = re.search(r'^## Prompt \(' + re.escape(name) + r'[^)]*\)\n\n((?:    .*\n|\n)+)', rules, re.M)
    assert m, name
    return '\n'.join(line[4:] for line in m.group(1).rstrip('\n').split('\n'))


rules = open('research/cosmo-discovery/blind-review.md', encoding='utf-8').read()
src = open(LAST, encoding='utf-8').read()
head = src[:src.index('<main class="paper">')]
old_title = re.search(r'<title>([^<]*)</title>', head).group(1)
head = head.replace(old_title, T)
head = re.sub(r'<meta name="description" content="[^"]*">', lambda m: f'<meta name="description" content="{DESC}">', head)
head = re.sub(r'<meta property="og:description" content="[^"]*">', lambda m: f'<meta property="og:description" content="{DESC}">', head)
head = head.replace('research-giant-arc.html', 'research-methodology.html')
head = head.replace('</style>', 'pre.prompt{font-family:ui-monospace,Consolas,Menlo,monospace;font-size:13px;line-height:1.5;'
                    'white-space:pre-wrap;background:var(--surface);border-top:1px solid var(--line);'
                    'border-bottom:1px solid var(--line);padding:12px;margin:12px 0 20px}\n</style>')
body = open('research/cosmo-discovery/briefs/paper_method_body.html', encoding='utf-8').read()
body = body.replace('@@TITLE@@', T).replace('@@REPO@@', REPO)
for key, name in PROMPTS.items():
    body = body.replace(f'@@{key}@@', html.escape(prompt(rules, name), quote=False))
assert '@@' not in body
open(SITE, 'w', encoding='utf-8', newline='\n').write(head + body)
s = open('site/sitemap.xml', encoding='utf-8').read()
if 'research-methodology' not in s:
    a = '  <url><loc>https://hermesai.dev/research-giant-arc.html</loc></url>\n'
    assert a in s
    s = s.replace(a, a + '  <url><loc>https://hermesai.dev/research-methodology.html</loc></url>\n')
    open('site/sitemap.xml', 'w', encoding='utf-8', newline='\n').write(s)
print('written', SITE, len(head + body))
