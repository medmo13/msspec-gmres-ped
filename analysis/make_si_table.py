"""Table of the ten extra clusters of the Supporting Information (theta scans from runs/p100).

Writes manuscript/si_table.tex (tabular body) and analysis/out/si_scans.csv. Numbers come from the saved
result.npz files of runs/p100/<case>/{lu_0,gmres_0,series1_0,series2_0,series3_0}.
"""
import csv
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from evaluate import metrics  # noqa: E402
from extra_metrics import measures  # noqa: E402

CASES = ['L36', 'L25', 'L26', 'L01', 'L03', 'H30', 'H46', 'H26', 'H41', 'H14']
mat = {r['mat_id']: r for r in csv.DictReader(open(ROOT / 'platform/materials_v2.csv'))}
R = ROOT / 'runs/p100'
rows = []
for c in CASES:
    ref = np.load(R / c / 'lu_0/result.npz')
    th, y = ref['theta'], ref['cross_section']
    out = dict(case=c, **{k: mat[c][k] for k in ('name', 'prototype', 'level', 'natoms', 'ke', 'rho')})
    for m in ('gmres', 'series1', 'series2', 'series3'):
        a = np.load(R / c / (m + '_0/result.npz'))['cross_section']
        out[m + '_eps2'] = 100 * metrics(a, y)['eps_2']
        out[m + '_pearson'] = measures(th, y, a)[3]
    rows.append(out)
with open(ROOT / 'analysis/out/si_scans.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def g(x):
    return '%.0f' % x if x >= 10 else '%.1f' % x if x >= 1 else '%.2f' % x


lines = []
for r in rows:
    lines.append('%s & %s & %s & %s & %.0f & %.2f & %s & %s / %s / %s & %s' % (
        r['name'], r['prototype'], r['level'].replace('_', '\\_'), r['natoms'], float(r['ke']), float(r['rho']), g(r['gmres_eps2']),
        g(r['series1_eps2']), g(r['series2_eps2']), g(r['series3_eps2']), '%.4f' % r['gmres_pearson']))
(ROOT / 'manuscript/si_table.tex').write_text(' \\\\\n'.join(lines) + '\n')
print(list(rows[0]))
print('\n'.join(lines))
