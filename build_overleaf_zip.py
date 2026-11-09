"""Zip the report and every figure it includes, ready to upload to Overleaf.

Usage (from the repository root):  python build_overleaf_zip.py
Writes scudem_problem_a_overleaf.zip with scudem_problem_a_report.tex and figures/<each included PNG>.
"""
import os
import re
import sys
import zipfile

ROOT = os.path.dirname(os.path.abspath(__file__))
TEX = 'scudem_problem_a_report.tex'
OUT = 'scudem_problem_a_overleaf.zip'

src = open(os.path.join(ROOT, TEX), encoding='utf-8').read()
body = '\n'.join(line for line in src.split('\n') if not line.lstrip().startswith('%'))
figs = sorted(set(re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}', body)))
missing = [f for f in figs if not os.path.exists(os.path.join(ROOT, 'figures', f))]
if missing:
    sys.exit('missing figures: %s' % missing)
with zipfile.ZipFile(os.path.join(ROOT, OUT), 'w', zipfile.ZIP_DEFLATED) as z:
    z.write(os.path.join(ROOT, TEX), TEX)
    for f in figs:
        z.write(os.path.join(ROOT, 'figures', f), 'figures/' + f)
print('%s: %s + %d figures' % (OUT, TEX, len(figs)))
