"""A4.8: compare each P100 dense-LU theta scan with the saved prod3000 inversion scan.

Usage: python analysis/check_lu_prod3000.py <runs/set> <out.csv>
For each case folder with lu_0/result.npz, the prod3000 scan_inv.npz of the same run_tag
(from lu_0/row.json) is loaded. Columns: case, run_tag, status, theta_equal, rel_l2,
max_rel, flag (rel_l2 > 1e-4). A flag is listed, not a failure.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

PROD = Path('/home/mohammed/Desktop/REAL-WORK/msspec_ML/prod3000')
LIMIT = 1e-4
FIELDS = ['case', 'run_tag', 'status', 'theta_equal', 'rel_l2', 'max_rel', 'flag']


def compare(a, ref):
    a = np.asarray(a, float)
    ref = np.asarray(ref, float)
    rel_l2 = float(np.linalg.norm(a - ref) / np.linalg.norm(ref))
    with np.errstate(divide='ignore', invalid='ignore'):
        max_rel = float(np.max(np.abs(a - ref) / np.abs(ref)))
    return rel_l2, max_rel


def check_set(set_dir, prod=PROD):
    out = []
    for case in sorted(p for p in Path(set_dir).iterdir() if p.is_dir()):
        lu = case / 'lu_0'
        if not (lu / 'row.json').exists():
            continue
        row = json.loads((lu / 'row.json').read_text())
        rec = {f: '' for f in FIELDS}
        rec.update(case=case.name, run_tag=row.get('run_tag', ''))
        ref_file = Path(prod) / 'runs' / rec['run_tag'] / 'scan_inv.npz'
        if not (lu / 'result.npz').exists():
            rec['status'] = 'no_lu_result'
        elif not ref_file.exists():
            rec['status'] = 'no_prod3000_file'
        else:
            mine, ref = np.load(lu / 'result.npz'), np.load(ref_file)
            if mine['cross_section'].shape != ref['cross_section'].shape:
                rec['status'] = 'shape_differs'
            else:
                rec['theta_equal'] = bool(np.allclose(mine['theta'], ref['theta'], atol=1e-9))
                rec['rel_l2'], rec['max_rel'] = compare(mine['cross_section'], ref['cross_section'])
                rec['status'] = 'compared'
                rec['flag'] = bool(rec['rel_l2'] > LIMIT or not rec['theta_equal'])
        out.append(rec)
    return out


def main(set_dir, out_csv):
    rows = check_set(set_dir)
    with open(out_csv, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    flagged = [r['case'] for r in rows if r['flag'] is True]
    print('wrote', out_csv, len(rows), 'rows; columns:', FIELDS)
    print('compared:', sum(r['status'] == 'compared' for r in rows), 'flagged:', flagged)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
