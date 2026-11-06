"""Endpoints P1 to P3 and summary S1 of PROTOCOL.md (A2.6, A4.3) from an evaluate.py CSV.

Usage: python analysis/stats.py <evaluate.csv> [--cases platform/materials_v2.csv] [--out analysis/out]
Writes <out>/<set>_endpoints.csv, <set> = the evaluate CSV name without '_eval'.
Rows with rep 0 only. Cases without a usable LU reference (no_reference) are excluded and
counted in column 'excluded'. A method failure counts as a fail. Rates have Wilson 95 %
intervals. S1: median eps_2 and eps_inf over GMRES runs that gave a result, with a
percentile bootstrap that resamples materials (10,000 resamples, seed 20261030).
Wilcoxon on log time ratios (used for the timing subset, step 5) is a function only.
"""
import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

SEED = 20261030
N_BOOT = 10000
Z = 1.959964
ROOT = Path(__file__).resolve().parents[1]
FIELDS = ['endpoint', 'method', 'stratum', 'what', 'k', 'n', 'rate', 'lo', 'hi', 'excluded', 'note']


def wilson(k, n, z=Z):
    if n == 0:
        return float('nan'), float('nan'), float('nan')
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, centre - half, centre + half


def boot_median_ci(values, groups, n_boot=N_BOOT, seed=SEED):
    """Median and percentile 95 % interval, resampling whole groups (materials)."""
    by = defaultdict(list)
    for v, g in zip(values, groups):
        by[g].append(v)
    keys = sorted(by)
    arrs = [np.asarray(by[k], float) for k in keys]
    rng = np.random.default_rng(seed)
    meds = np.empty(n_boot)
    for i in range(n_boot):
        pick = rng.integers(0, len(arrs), len(arrs))
        meds[i] = np.median(np.concatenate([arrs[j] for j in pick]))
    return float(np.median(values)), float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


def wilcoxon_log_ratio(t_method, t_ref):
    """Paired Wilcoxon signed-rank test on log(t_method / t_ref). Returns median ratio, p."""
    from scipy.stats import wilcoxon
    d = np.log(np.asarray(t_method, float)) - np.log(np.asarray(t_ref, float))
    res = wilcoxon(d)
    return float(np.exp(np.median(d))), float(res.pvalue)


def flag(x):
    return str(x) == 'True'


def load_cases(path):
    return {r['mat_id']: r for r in csv.DictReader(open(path))}


def endpoints(rows, cases):
    rows = [r for r in rows if r['rep'] == '0' and r['method'] != 'lu']
    out = []
    strata = [('all', None), ('high', 'high'), ('low', 'low')]

    def pick(method, st):
        return [r for r in rows if r['method'] == method
                and (st is None or cases[r['case']]['stratum'] == st)]

    def add(endpoint, method, stname, what, k, n, excl, note=''):
        p, lo, hi = wilson(k, n)
        out.append(dict(endpoint=endpoint, method=method, stratum=stname, what=what, k=k, n=n,
                        rate=p, lo=lo, hi=hi, excluded=excl, note=note))

    for stname, st in strata:
        rs = pick('gmres', st)
        scored = [r for r in rs if r['status'] in ('ok', 'method_failure')]
        excl = len(rs) - len(scored)
        add('P1', 'gmres', stname, 'pass', sum(flag(r['passed']) for r in scored), len(scored), excl)
        fails = [r for r in scored if r['status'] == 'method_failure']
        reasons = Counter(r['error_reason'] for r in fails)
        add('P2', 'gmres', stname, 'method_failure', len(fails), len(scored), excl,
            '; '.join('%s=%d' % kv for kv in sorted(reasons.items())))
        for order in (1, 2, 3):
            m = 'series%d' % order
            rs = pick(m, st)
            scored = [r for r in rs if r['status'] in ('ok', 'method_failure')]
            excl = len(rs) - len(scored)
            add('P3', m, stname, 'pass', sum(flag(r['passed']) for r in scored), len(scored), excl)
            add('P3', m, stname, 'diverged', sum(flag(r['diverged']) for r in scored), len(scored), excl)
    # S1: accuracy of GMRES runs that gave a result
    ok = [r for r in rows if r['method'] == 'gmres' and r['status'] == 'ok']
    for stname, st in strata:
        sel = [r for r in ok if st is None or cases[r['case']]['stratum'] == st]
        for metric in ('eps_2', 'eps_inf'):
            if not sel:
                continue
            vals = [float(r[metric]) for r in sel]
            groups = [cases[r['case']]['name'] for r in sel]
            med, lo, hi = boot_median_ci(vals, groups)
            out.append(dict(endpoint='S1', method='gmres', stratum=stname, what='median_' + metric,
                            k='', n=len(vals), rate=med, lo=lo, hi=hi, excluded='',
                            note='bootstrap by material, %d resamples, seed %d' % (N_BOOT, SEED)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('eval_csv')
    ap.add_argument('--cases', default=str(ROOT / 'platform' / 'materials_v2.csv'))
    ap.add_argument('--out', default=str(ROOT / 'analysis' / 'out'))
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.eval_csv)))
    out = endpoints(rows, load_cases(a.cases))
    name = Path(a.eval_csv).stem
    name = name[:-5] if name.endswith('_eval') else name
    Path(a.out).mkdir(parents=True, exist_ok=True)
    dest = Path(a.out) / (name + '_endpoints.csv')
    with open(dest, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(out)
    print('wrote', dest, len(out), 'rows; columns:', FIELDS)


if __name__ == '__main__':
    main()
