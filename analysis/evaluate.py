"""Accuracy metrics of PROTOCOL.md section 3, computed from saved arrays only.

Run folder layout (written by code/runners/run_set.py):
    runs/<set>/<case>/<method>_<rep>/row.json   one result row (rc, error, error_reason, ...)
    runs/<set>/<case>/<method>_<rep>/result.npz arrays cross_section, direct_signal, theta, phi

The reference of every case is its dense-LU run (method 'lu', repetition 0). Any other
run is compared with it.

Definitions (arrays a and ref are intensities on the same directions):
    eps_2           ||a - ref||_2 / ||ref||_2
    eps_inf         max |a - ref| / max(|ref|, 0.01 * max|ref|)       (primary, floor 1 %)
    eps_inf_nofloor max |a - ref| / |ref|                             (supporting information)
    R_chi           sum (chi - chi_ref)^2 / sum (chi^2 + chi_ref^2),  chi = (I - I0) / I0,
                    I0 the direct signal; points with I0 <= 0 are left out.
    pass            eps_2 < 0.10 and eps_inf < 0.10, all values finite, and the run did not fail.
diverged        series methods only (A4.3): non-finite result, or eps_2 >= 1. Empty for other methods.
Warm-up runs (rep label 'w') are ignored.
A failed run (non-zero return code, error flag, missing or non-finite result) counts as a
fail with its reason; it is never replaced.

Usage: python analysis/evaluate.py <runs/set> <out.csv>
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

FLOOR = 0.01
LIMIT = 0.10
DIVERGED = 1.0


def metrics(a, ref, i0=None, i0_ref=None, floor=FLOOR):
    a = np.asarray(a, dtype=float)
    ref = np.asarray(ref, dtype=float)
    if a.shape != ref.shape:
        raise ValueError('shape %s against reference %s' % (a.shape, ref.shape))
    n_bad = int(np.sum(~np.isfinite(a)))
    out = {'n_dir': int(a.size), 'n_nonfinite': n_bad}
    if n_bad:
        out.update(eps_2=np.inf, eps_inf=np.inf, eps_inf_nofloor=np.inf, r_chi=np.inf)
        return out
    diff = np.abs(a - ref)
    out['eps_2'] = float(np.linalg.norm(a - ref) / np.linalg.norm(ref))
    out['eps_inf'] = float(np.max(diff / np.maximum(np.abs(ref), floor * np.max(np.abs(ref)))))
    with np.errstate(divide='ignore', invalid='ignore'):
        rel = np.where(ref == 0, np.where(diff == 0, 0.0, np.inf), diff / np.abs(ref))
    out['eps_inf_nofloor'] = float(np.max(rel))
    out['r_chi'] = r_chi(a, ref, i0, i0_ref)
    return out


def r_chi(a, ref, i0, i0_ref):
    if i0 is None or i0_ref is None:
        return float('nan')
    i0 = np.asarray(i0, dtype=float)
    i0_ref = np.asarray(i0_ref, dtype=float)
    ok = (i0 > 0) & (i0_ref > 0)
    if not ok.any():
        return float('nan')
    chi = (a[ok] - i0[ok]) / i0[ok]
    chi_ref = (ref[ok] - i0_ref[ok]) / i0_ref[ok]
    den = np.sum(chi ** 2 + chi_ref ** 2)
    return float(np.sum((chi - chi_ref) ** 2) / den) if den > 0 else float('nan')


def passes(m):
    return bool(m['n_nonfinite'] == 0 and m['eps_2'] < LIMIT and m['eps_inf'] < LIMIT)


def diverged(method, rec):
    """True/False for series methods, '' for all others (A4.3)."""
    if not method.startswith('series'):
        return ''
    if rec.get('error_reason') == 'nonfinite_result' or rec.get('n_nonfinite', 0) > 0:
        return True
    return bool(rec.get('eps_2', 0.0) >= DIVERGED)


def load_run(d):
    row = json.loads((d / 'row.json').read_text())
    npz = d / 'result.npz'
    arrays = dict(np.load(npz)) if npz.exists() else None
    return row, arrays


def evaluate_set(set_dir):
    set_dir = Path(set_dir)
    out = []
    for case in sorted(p for p in set_dir.iterdir() if p.is_dir()):
        runs = {d.name: load_run(d) for d in sorted(case.iterdir())
                if d.is_dir() and (d / 'row.json').exists() and not d.name.endswith('_w')}
        ref_row, ref_arr = runs.get('lu_0', (None, None))
        ref_ok = ref_arr is not None and ref_row.get('error') in (0, '0', False, None, '')
        for name, (row, arr) in runs.items():
            method, rep = name.rsplit('_', 1)
            rec = {'case': case.name, 'method': method, 'rep': rep, 'rc': row.get('rc', ''),
                   'error_reason': row.get('error_reason', '')}
            failed = str(row.get('error', 0)) not in ('0', 'False', '') or arr is None
            if method == 'lu':
                st = 'reference' if not failed else (
                    'lu_not_run' if row.get('run_status') == 'lu_not_run' else 'lu_failed')
                rec.update(status=st, passed='')
            elif not ref_ok:
                rec.update(status='no_reference', passed='')
            elif failed:
                rec.update(status='method_failure', passed=False)
            else:
                m = metrics(arr['cross_section'], ref_arr['cross_section'],
                            arr.get('direct_signal'), ref_arr.get('direct_signal'))
                rec.update(m)
                rec.update(status='ok', passed=passes(m))
            rec['diverged'] = diverged(method, rec)
            out.append(rec)
    return out


FIELDS = ['case', 'method', 'rep', 'rc', 'status', 'passed', 'error_reason', 'n_dir',
          'n_nonfinite', 'diverged', 'eps_2', 'eps_inf', 'eps_inf_nofloor', 'r_chi']


def main(set_dir, out_csv):
    rows = evaluate_set(set_dir)
    with open(out_csv, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    print('wrote', out_csv, len(rows), 'rows; columns:', FIELDS)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
