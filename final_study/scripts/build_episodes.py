"""Rebuild socket episodes from the raw Ghidini casefile workbook, keeping child IDs and dates.

Rules (same as the original data pack):
  * an episode starts at every visit where a new prosthesis is provided;
  * event = first later visit whose clinical evaluation is Growth with subcategory tight and/or short;
  * censored = next new prosthesis for a non-growth reason, or the child's last visit;
  * zero-length episodes are dropped.
Adds: child_id, start_date, end_date, years_since_first_visit (of the episode start), episode_index,
age_at_amputation, n_episodes_child.
Usage: python -I build_episodes.py <xlsx> <out_csv>
"""
import sys, csv, re, datetime as dt
import openpyxl
sys.stdout.reconfigure(encoding='utf-8')
src, out = sys.argv[1], sys.argv[2]
wb = openpyxl.load_workbook(src, read_only=True, data_only=True)

ages = {}
for r in list(wb['Demographics'].iter_rows(values_only=True))[1:]:
    if r[0] and str(r[0]).startswith('CF'):
        a = r[4]
        ages[r[0]] = float(a) if isinstance(a, (int, float)) else None

LEVEL = {'TT': 'transtibial', 'KD': 'knee_disarticulation', 'TF': 'transfemoral'}

def norm(x):
    return '' if x is None else str(x).strip().lower()

def growth_cause(what, sub):
    """Return 'tight', 'short', 'tight+short' or '' for a visit."""
    # Containment, not positional pairing: the subcategory order does not reliably follow the evaluation
    # order in the casefiles (e.g. "Broken, Growth" / "Short, unspecified what broken").
    w, s = norm(what), norm(sub)
    if 'growth' not in w:
        return ''
    tight = 'tight' in s
    short = 'short' in s
    return 'tight+short' if (tight and short) else 'tight' if tight else 'short' if short else ''

def is_new(solution, what):
    return 'new prosthesis' in norm(solution) or 'new prosthesis' in norm(what)

rows_out = []
for sheet, level in LEVEL.items():
    rows = list(wb[sheet].iter_rows(values_only=True))[2:]
    visits = {}
    for r in rows:
        if not r[0] or not isinstance(r[3], dt.datetime):
            continue
        visits.setdefault(r[0], []).append(r)
    for cid, vs in visits.items():
        vs.sort(key=lambda r: r[3])
        first = vs[0][3]
        last = vs[-1][3]
        starts = [i for i, v in enumerate(vs) if is_new(v[8], v[6])]
        for k, i in enumerate(starts):
            t0 = vs[i][3]
            end, event, cause, action = None, 0, '', ''
            for j in range(i + 1, len(vs)):
                c = growth_cause(vs[j][6], vs[j][7])
                if c:
                    end, event, cause = vs[j][3], 1, c
                    action = 'new_prosthesis' if is_new(vs[j][8], vs[j][6]) else 'adjustment'
                    break
                if is_new(vs[j][8], vs[j][6]):
                    end = vs[j][3]
                    break
            if end is None:
                end = last
            days = (end - t0).days
            if days <= 0:
                continue
            rows_out.append(dict(child_id=cid, level=level, start_date=t0.date().isoformat(), end_date=end.date().isoformat(),
                                 days=days, months=round(days / 30.4375, 2), event=event, cause=cause, action=action,
                                 years_since_first_visit=round((t0 - first).days / 365.25, 3),
                                 age_at_amputation=ages.get(cid), episode_index=k + 1))
n_by_child = {}
for r in rows_out:
    n_by_child[r['child_id']] = n_by_child.get(r['child_id'], 0) + 1
for r in rows_out:
    r['n_episodes_child'] = n_by_child[r['child_id']]
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()))
    w.writeheader(); w.writerows(rows_out)

ev = [r for r in rows_out if r['event'] == 1]
from collections import Counter
print('episodes', len(rows_out), Counter(r['level'] for r in rows_out))
print('events', len(ev), Counter(r['cause'] for r in ev), Counter(r['action'] for r in ev))
print('children', len(n_by_child), '| episodes with known age at amputation', sum(r['age_at_amputation'] is not None for r in rows_out))
