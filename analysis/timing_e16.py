"""Timing and memory of the E16 cases: warm-up run (rep 'w') and repeats 0 to 2 (A2.8).

Reads runs/e16_timing/<case>/<method>_<rep>/row.json (written by code/runners/run_set.py).
Writes analysis/out/e16_timing.csv (one row per case and method) and analysis/out/e16_timing.npz.
"""
import csv
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / 'runs' / 'e16_timing'
CASES = ['H11', 'H10', 'H08', 'L45']
METHODS = ['lu', 'gmres', 'series1', 'series2', 'series3']


def row(case, method, rep):
    p = RUN / case / ('%s_%s' % (method, rep)) / 'row.json'
    return json.loads(p.read_text())


def main():
    out, arrays = [], {}
    for c in CASES:
        for m in METHODS:
            w = row(c, m, 'w')
            reps = [row(c, m, r) for r in range(3)]
            wall = np.array([float(r['wall_s']) for r in reps])
            rss = np.array([float(r['peak_rss_mib']) for r in reps])
            load_b = np.array([float(r['load_before']) for r in reps])
            load_a = np.array([float(r['load_after']) for r in reps])
            rec = dict(case=c, method=m, warmup_wall_s=float(w['wall_s']),
                       wall_median_s=float(np.median(wall)), wall_min_s=float(wall.min()),
                       wall_max_s=float(wall.max()), peak_rss_mib=float(rss.max()),
                       load_before_max=float(load_b.max()), load_after_max=float(load_a.max()),
                       natoms=reps[0]['natoms'], rc=max(int(r['rc']) for r in reps + [w]))
            out.append(rec)
            arrays['%s_%s_wall' % (c, m)] = wall
            arrays['%s_%s_rss' % (c, m)] = rss
    (ROOT / 'analysis' / 'out').mkdir(parents=True, exist_ok=True)
    with open(ROOT / 'analysis' / 'out' / 'e16_timing.csv', 'w', newline='') as fh:
        w_ = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w_.writeheader(); w_.writerows(out)
    np.savez(ROOT / 'analysis' / 'out' / 'e16_timing.npz', **arrays)
    print('wrote analysis/out/e16_timing.csv', len(out), 'rows and e16_timing.npz')


if __name__ == '__main__':
    main()
