"""Four more accuracy measures against dense LU, computed from saved result.npz only.

For every case with a dense-LU run (runs/<set>/<case>/lu_0) and every method run that gave
a result, with ref the LU cross section and a the test cross section on the same theta points:
    eps_rel   100 ||a - ref|| / ||ref||                     (percent; equals 100 eps_2 of evaluate.py)
    R         theta-adapted Pendry R-factor, 0 is a perfect match:
                  L = (dY/dtheta) / (Y^2 + V^2),  V = 0.05 max|ref|,  derivative along theta in degrees
                  R = sum (L_a - L_ref)^2 / sum (L_a^2 + L_ref^2)
              Pendry's original R is defined along the energy; this is the same form along theta.
    R2        sqrt( sum (a - ref)^2 / sum (ref - mean ref)^2 ), not normalised, can exceed 1
    pearson   Pearson correlation of a and ref
Same definitions as compare() in msspec_SE/msspec_standard/material_platform/platform_common.py
(eps_rel, R); R2 and pearson as asked by the author.

Label gmres_e2c67 is GMRES at tolerance 1e-2 with a cap of 12 steps where that converged, and the cap-67 run
(set p100_cap67) for the clusters that did not. A converged run stops at the same step for either cap.
Label gmres_e6c67 is the same for tolerance 1e-6 (cap 30 where it converged, set p100_tol6_cap67 for the others).
Writes analysis/out/p100_extra_metrics.csv (one row per run) and p100_extra_summary.csv (median and quartiles
per label). Usage: python analysis/extra_metrics.py [case ...]
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
# (set, method run name, label)
RUNS = [('p100', 'gmres_0', 'gmres_e2'), ('p100_tol3', 'gmres_0', 'gmres_e3'),
        ('p100_tol6', 'gmres_0', 'gmres_e6'), ('p100', 'series1_0', 'series1'),
        ('p100', 'series2_0', 'series2'), ('p100', 'series3_0', 'series3'),
        ('p100_cap67', 'gmres_0', 'gmres_e2_cap67'),
        ('p100_tol6_cap67', 'gmres_0', 'gmres_e6_cap67')]
FIELDS = ['set', 'case', 'label', 'eps_rel_pct', 'R', 'R2', 'pearson']


def measures(theta, ref, a):
    eps_rel = 100 * np.linalg.norm(a - ref) / np.linalg.norm(ref)
    v = 0.05 * np.max(np.abs(ref))
    l_ref = np.gradient(ref, theta) / (ref ** 2 + v ** 2)
    l_a = np.gradient(a, theta) / (a ** 2 + v ** 2)
    r = np.sum((l_a - l_ref) ** 2) / np.sum(l_a ** 2 + l_ref ** 2)
    r2 = np.sqrt(np.sum((a - ref) ** 2) / np.sum((ref - ref.mean()) ** 2))
    pear = np.corrcoef(a, ref)[0, 1]
    return float(eps_rel), float(r), float(r2), float(pear)


def ok_run(d):
    rj, npz = d / 'row.json', d / 'result.npz'
    if not (rj.exists() and npz.exists()):
        return None
    if str(json.loads(rj.read_text()).get('error', 0)) not in ('0', 'False', ''):
        return None
    return np.load(npz)


def main(only):
    rows = []
    for s, name, label in RUNS:
        sd = ROOT / 'runs' / s
        if not sd.exists():
            continue
        for case in sorted(p for p in sd.iterdir() if p.is_dir()):
            if only and case.name not in only:
                continue
            ref, a = ok_run(case / 'lu_0'), ok_run(case / name)
            if ref is None or a is None or not np.all(np.isfinite(a['cross_section'])):
                continue
            e, r, r2, p = measures(ref['theta'], ref['cross_section'], a['cross_section'])
            rows.append(dict(set=s, case=case.name, label=label, eps_rel_pct=e, R=r, R2=r2, pearson=p))
    have = {r['case'] for r in rows if r['label'] == 'gmres_e2'}
    for r in [r for r in rows if r['label'] == 'gmres_e2_cap67']:
        if r['case'] not in have:
            rows.append(dict(r, label='gmres_e2c67'))
    rows += [dict(r, label='gmres_e2c67') for r in rows if r['label'] == 'gmres_e2']
    have = {r['case'] for r in rows if r['label'] == 'gmres_e6'}
    rows += [dict(r, label='gmres_e6c67') for r in rows if r['label'] == 'gmres_e6']
    rows += [dict(r, label='gmres_e6c67') for r in rows if r['label'] == 'gmres_e6_cap67' and r['case'] not in have]
    out = ROOT / 'analysis/out/p100_extra_metrics.csv'
    if not only:
        summary = []
        for lab in dict.fromkeys(r['label'] for r in rows):
            s = [r for r in rows if r['label'] == lab]
            d = dict(label=lab, n=len(s))
            for m in ('eps_rel_pct', 'R', 'R2', 'pearson'):
                v = np.array([r[m] for r in s], dtype=float)
                d[m + '_med'], d[m + '_q25'], d[m + '_q75'] = np.median(v), np.percentile(v, 25), np.percentile(v, 75)
            summary.append(d)
        with open(ROOT / 'analysis/out/p100_extra_summary.csv', 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)
        with open(out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS)
            w.writeheader(); w.writerows(rows)
        print('wrote', out, len(rows), 'rows; columns:', FIELDS)
    else:
        for r in rows:
            print(r)


if __name__ == '__main__':
    main(set(sys.argv[1:]))
