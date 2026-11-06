"""Tests of analysis/stats.py and make_numbers.py on made-up data with known counts.

Run: python analysis/test_stats.py   (work files go to platform/work/)
"""
import csv
import math
import tempfile
from pathlib import Path

import numpy as np

import make_numbers as mn
import stats as st

ROOT = Path(__file__).resolve().parents[1]


def row(case, method, status, passed='', reason='', eps_2='', eps_inf='', diverged='', rep='0'):
    return dict(case=case, method=method, rep=rep, status=status, passed=passed, error_reason=reason,
                eps_2=eps_2, eps_inf=eps_inf, diverged=diverged)


def main():
    # 1. Wilson, 8 of 10, from the textbook formula written out here
    z, n, k = 1.959964, 10, 8
    p = k / n
    d = 1 + z * z / n
    lo = (p + z * z / (2 * n) - z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / d
    hi = (p + z * z / (2 * n) + z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / d
    r, l, h = st.wilson(8, 10)
    assert abs(l - lo) < 1e-12 and abs(h - hi) < 1e-12 and r == 0.8
    assert round(l, 3) == 0.490 and round(h, 3) == 0.943, (l, h)
    assert st.wilson(0, 10)[1] == 0 or abs(st.wilson(0, 10)[1]) < 1e-12
    assert math.isnan(st.wilson(0, 0)[0])
    # 2. bootstrap: same seed, same answer; resamples materials
    vals = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06]
    grp = ['a', 'a', 'b', 'b', 'c', 'c']
    r1 = st.boot_median_ci(vals, grp, n_boot=500)
    r2 = st.boot_median_ci(vals, grp, n_boot=500)
    assert r1 == r2 and r1[0] == np.median(vals) and 0.01 <= r1[1] <= r1[2] <= 0.06
    # 3. Wilcoxon on log ratios: a constant ratio 0.5 gives median ratio 0.5 and a small p
    t_ref = np.arange(1.0, 13.0)
    med, pval = st.wilcoxon_log_ratio(0.5 * t_ref * np.exp(0.01 * np.sin(np.arange(12))), t_ref)
    assert abs(med - 0.5) < 0.01 and pval < 0.01
    # 4. P1 to P3 on a made-up table with known counts
    cases = {}
    for i in range(10):
        cases['C%02d' % i] = dict(mat_id='C%02d' % i, stratum='high' if i < 6 else 'low', name='M%d' % (i // 2))
    rows = []
    for i in range(10):
        c = 'C%02d' % i
        rows.append(row(c, 'lu', 'reference'))
        if i == 9:                                              # no reference: excluded
            rows.append(row(c, 'gmres', 'no_reference'))
            continue
        if i == 8:                                              # method failure (low)
            rows.append(row(c, 'gmres', 'method_failure', False, 'none_available'))
        elif i in (0, 7):                                       # accuracy fail
            rows.append(row(c, 'gmres', 'ok', False, '', 0.2, 0.3))
        else:
            rows.append(row(c, 'gmres', 'ok', True, '', 0.01 * (i + 1), 0.02 * (i + 1)))
        rows.append(row(c, 'series1', 'ok', i >= 5, '', 0.05, 0.05, False))
        rows.append(row(c, 'series2', 'ok', False, '', 3.0, 5.0, True))
        rows.append(row(c, 'series3', 'method_failure', False, 'nonfinite_result', '', '', True))
    rows.append(row('C00', 'gmres', 'ok', True, '', 0.0, 0.0, rep='1'))        # rep 1 is not counted
    got = {(r['endpoint'], r['method'], r['stratum'], r['what']): r for r in st.endpoints(rows, cases)}
    g = got[('P1', 'gmres', 'all', 'pass')]
    assert (g['k'], g['n'], g['excluded']) == (6, 9, 1), g          # 9 scored: 6 pass, 2 acc fail, 1 failure
    g = got[('P1', 'gmres', 'high', 'pass')]; assert (g['k'], g['n']) == (5, 6), g
    g = got[('P1', 'gmres', 'low', 'pass')]; assert (g['k'], g['n'], g['excluded']) == (1, 3, 1), g
    g = got[('P2', 'gmres', 'all', 'method_failure')]
    assert (g['k'], g['n']) == (1, 9) and g['note'] == 'none_available=1', g
    g = got[('P3', 'series1', 'all', 'pass')]; assert (g['k'], g['n']) == (4, 9), g
    g = got[('P3', 'series2', 'all', 'diverged')]; assert (g['k'], g['n']) == (9, 9), g
    g = got[('P3', 'series3', 'low', 'diverged')]; assert (g['k'], g['n']) == (3, 3), g
    s1 = got[('S1', 'gmres', 'all', 'median_eps_2')]
    assert s1['n'] == 8 and s1['lo'] <= s1['rate'] <= s1['hi']
    # 5. CSV round trip and numbers.tex (digits spelled out, values traced to the CSV)
    work = ROOT / 'platform' / 'work'
    work.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='stats_test_', dir=work) as td:
        out = Path(td) / 'out'
        out.mkdir()
        with open(out / 'toy_endpoints.csv', 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=st.FIELDS)
            w.writeheader()
            w.writerows(st.endpoints(rows, cases))
        tex = Path(td) / 'numbers.tex'
        import sys
        old = sys.argv
        sys.argv = ['make_numbers', '--out-dir', str(out), '--tex', str(tex)]
        mn.main()
        sys.argv = old
        text = tex.read_text()
        assert '\\newcommand{\\ToyPOneGmresAllPassK}{6}' in text, text[:600]
        assert '\\newcommand{\\ToyPOneGmresAllPassN}{9}' in text
        assert '\\newcommand{\\ToyPOneGmresAllPassRate}{66.7}' in text
        assert not any(c.isdigit() for line in text.splitlines()[1:] for c in line.split('}{')[0])
    print('PASS: stats.py and make_numbers.py tests (5 groups)')


if __name__ == '__main__':
    main()
