"""The "Built on" section every research page opens with (Will, 2026-11-02): the people, teams and surveys whose data and
methods the work uses, first and foremost, and the work offered back to them.

The credits are in built_on.json, keyed by page (site/research-<key>.html). The generators (mkpaper_*.py) call insert() on
the page they build, so a new paper carries the section once it has an entry; a page with no generator is updated in
place by running this file:

    python3 research/cosmo-discovery/briefs/built_on.py            # every page with an entry
    python3 research/cosmo-discovery/briefs/built_on.py giant-arc  # one page

insert() is idempotent: it replaces the section between its markers, or adds it after the abstract (after the date line
on a withdrawal notice), and adds the section's styles to the head once.
"""
import html
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
START, END = '<!-- built-on -->', '<!-- /built-on -->'
CSS_MARK = '/* built-on */'
CSS = (CSS_MARK + '\n.built-on{margin:28px auto 0;max-width:640px;font-size:15px;line-height:1.55;border-left:3px solid var(--line,#ccc);'
       'padding:2px 0 2px 16px}\n.built-on h2{margin:0 0 8px}\n.built-on ul{margin:6px 0 10px;padding-left:20px}'
       '\n.built-on li{margin:4px 0}\n.built-on p{margin:6px 0}\n')
ASK = ('This work is offered back to them. If it uses your data or method, we would welcome your review, and would gladly '
       'have you co-author it, direct what is done next, or take it over; or tell us what to correct or take down. '
       'Write to <a href="mailto:team@hermesai.dev">team@hermesai.dev</a>.')


def credits():
    return json.load(open(os.path.join(HERE, 'built_on.json'), encoding='utf-8'))


def section(key, lead=None):
    c = credits()[key]
    items = []
    for e in c['credits']:
        cite = re.sub(r'arXiv:(\d{4}\.\d{4,5})', r'<a href="https://arxiv.org/abs/\1">arXiv:\1</a>', html.escape(e['cite'], quote=False))
        if e.get('arxiv') and e['arxiv'] not in cite:
            cite += f' (<a href="https://arxiv.org/abs/{e["arxiv"]}">arXiv:{e["arxiv"]}</a>)'
        items.append(f'  <li><b>{html.escape(e["who"], quote=False)}</b>: {html.escape(e["what"], quote=False)} ({cite}).</li>')
    lead = lead or c.get('lead') or ('This work stands on data and methods built by human researchers:')
    return (f'{START}\n<section class="built-on" id="built-on">\n<h2>Built on</h2>\n<p>{lead}</p>\n<ul>\n'
            + '\n'.join(items) + f'\n</ul>\n<p>{ASK}</p>\n</section>\n{END}')


def insert(page, key):
    """The page with its Built-on section (and the section's styles) in place."""
    sec = section(key)
    if START in page:
        page = page[:page.index(START)] + sec + page[page.index(END) + len(END):]
    elif '</header>' in page:
        i = page.index('</header>') + len('</header>')
        page = page[:i] + '\n' + sec + page[i:]
    else:  # a withdrawal notice: after its date line
        m = re.search(r'<p class="date">.*?</p>', page, re.S)
        page = page[:m.end()] + '\n' + sec + page[m.end():]
    if CSS_MARK not in page:
        page = page.replace('</style>', CSS + '</style>', 1)
    return page


def apply(key, root='.'):
    path = os.path.join(root, 'site', f'research-{key}.html')
    page = open(path, encoding='utf-8').read()
    out = insert(page, key)
    if out != page:
        open(path, 'w', encoding='utf-8', newline='\n').write(out)
    return out != page


if __name__ == '__main__':
    for k in (sys.argv[1:] or list(credits())):
        print(k, 'updated' if apply(k) else 'unchanged')
