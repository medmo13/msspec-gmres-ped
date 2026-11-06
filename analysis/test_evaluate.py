"""Tests of analysis/evaluate.py on synthetic arrays with known errors.

Run: python analysis/test_evaluate.py   (work files go to platform/work/)
"""
import json
import tempfile
from pathlib import Path

import numpy as np

import evaluate as ev

ROOT = Path(__file__).resolve().parents[1]


def close(x, y, tol=1e-12):
    assert abs(x - y) <= tol * max(1.0, abs(y)), (x, y)


def main():
    x = np.linspace(0, 1, 200)
    ref = 1.0 + 0.5 * np.sin(6 * x)
    # 1. uniform 5 % error: eps_2 = eps_inf = 0.05, floor not engaged
    m = ev.metrics(1.05 * ref, ref)
    close(m['eps_2'], 0.05); close(m['eps_inf'], 0.05); close(m['eps_inf_nofloor'], 0.05)
    assert ev.passes(m)
    # 2. 12 % error fails the gate
    assert not ev.passes(ev.metrics(1.12 * ref, ref))
    # 3. floor: one tiny reference value, absolute error 1e-3, max|ref| = 1
    r2 = np.ones(50); r2[7] = 1e-4
    a2 = r2.copy(); a2[7] += 1e-3
    m = ev.metrics(a2, r2)
    close(m['eps_inf'], 1e-3 / 0.01)         # divided by the floor 0.01 * max
    close(m['eps_inf_nofloor'], 10.0)
    # 4. non-finite value: fail, counted
    a3 = ref.copy(); a3[3] = np.nan
    m = ev.metrics(a3, ref)
    assert m['n_nonfinite'] == 1 and not ev.passes(m) and np.isinf(m['eps_2'])
    # 5. R_chi with I0 = 1: chi = I - 1. ref chi = 0.2 everywhere, a chi = 0.1 everywhere
    i0 = np.ones(10)
    m = ev.metrics(1.1 * i0, 1.2 * i0, i0, i0)
    close(m['r_chi'], (0.1 ** 2) / (0.1 ** 2 + 0.2 ** 2))
    # 6. identical arrays: zero error
    m = ev.metrics(ref, ref, i0=np.ones_like(ref), i0_ref=np.ones_like(ref))
    assert m['eps_2'] == 0 and m['eps_inf'] == 0 and m['r_chi'] == 0
    # 7. evaluate_set: reference, passing run, failing accuracy, error row, non-finite run, no result
    work = ROOT / 'platform' / 'work'
    work.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='eval_test_', dir=work) as td:
        case = Path(td) / 'caseA'

        def put(name, arr, **row):
            d = case / name
            d.mkdir(parents=True)
            row.setdefault('rc', 0); row.setdefault('error', 0)
            (d / 'row.json').write_text(json.dumps(row))
            if arr is not None:
                np.savez(d / 'result.npz', cross_section=arr, direct_signal=np.ones_like(arr))
        put('lu_0', ref)
        put('gmres_0', 1.03 * ref)
        put('k0_0', 1.5 * ref)
        put('gmres_1', None, rc=99, error=1, error_reason='none_available')
        bad = ref.copy(); bad[0] = np.inf
        put('series1_0', bad)
        put('series2_0', 1.5 * ref + 2.0)                    # eps_2 well above 1 -> diverged
        put('series3_0', 1.02 * ref)                         # fine -> not diverged
        put('series1_w', 1.02 * ref)                         # warm-up folder: must be ignored
        put('gmres_w', 1.02 * ref)
        rows = {(r['method'], r['rep']): r for r in ev.evaluate_set(td)}
        assert rows[('lu', '0')]['status'] == 'reference'
        r = rows[('gmres', '0')]; assert r['status'] == 'ok' and r['passed'] is True
        close(r['eps_2'], 0.03)
        r = rows[('k0', '0')]; assert r['status'] == 'ok' and r['passed'] is False
        r = rows[('gmres', '1')]
        assert r['status'] == 'method_failure' and r['passed'] is False \
            and r['error_reason'] == 'none_available'
        r = rows[('series1', '0')]; assert r['passed'] is False and r['n_nonfinite'] == 1
        assert ('series1', 'w') not in rows and ('gmres', 'w') not in rows
        assert rows[('series1', '0')]['diverged'] is True       # NaN/inf value
        assert rows[('series2', '0')]['eps_2'] >= 1 and rows[('series2', '0')]['diverged'] is True
        assert rows[('series3', '0')]['diverged'] is False
        assert rows[('gmres', '0')]['diverged'] == '' and rows[('k0', '0')]['diverged'] == ''
        # a series run that failed with nonfinite_result is diverged, though it has no array
        put('series2_1', None, rc=0, error=1, error_reason='nonfinite_result')
        rows = {(r['method'], r['rep']): r for r in ev.evaluate_set(td)}
        assert rows[('series2', '1')]['status'] == 'method_failure' and rows[('series2', '1')]['diverged'] is True
        # a case without a reference is reported, not scored
        (case / 'lu_0' / 'result.npz').unlink()
        rows = {(r['method'], r['rep']): r for r in ev.evaluate_set(td)}
        assert rows[('gmres', '0')]['status'] == 'no_reference'
        assert rows[('lu', '0')]['status'] == 'lu_failed'
        # a skipped LU (memory rule) is 'lu_not_run', not a failure of any method
        row = json.loads((case / 'lu_0' / 'row.json').read_text())
        row.update(run_status='lu_not_run', error=0)
        (case / 'lu_0' / 'row.json').write_text(json.dumps(row))
        rows = {(r['method'], r['rep']): r for r in ev.evaluate_set(td)}
        assert rows[('lu', '0')]['status'] == 'lu_not_run'
        assert rows[('gmres', '0')]['status'] == 'no_reference'
    print('PASS: evaluate.py tests (8 groups)')


if __name__ == '__main__':
    main()
