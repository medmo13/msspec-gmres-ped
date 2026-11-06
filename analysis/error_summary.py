"""Error of each method against dense LU, over the 97 cases with a reference.

Reads analysis/out/{p100,p100_tol3,p100_tol6}_eval.csv (evaluate.py); gmres_e2c67 adds the cap-67 rerun p100_cap67_eval.csv. For GMRES at tolerance 1e-2, 1e-3,
1e-6 and the series of order 1 to 3, and for all, strong (rho > 1) and weak cases, writes
analysis/out/p100_error_summary.csv: number of cases, cases with a result, median and quartiles of eps_2 and
eps_inf over the cases with a result, and the number of cases with eps_2 above 1 (error larger than the signal).
Usage: python analysis/error_summary.py
"""
import csv
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'analysis' / 'out'
mat = {r['mat_id']: r for r in csv.DictReader(open(ROOT / 'platform/materials_v2.csv'))}
SETS = [('gmres_e2', 'p100_eval.csv', 'gmres'), ('gmres_e3', 'p100_tol3_eval.csv', 'gmres'),
        ('gmres_e6', 'p100_tol6_eval.csv', 'gmres'), ('series1', 'p100_eval.csv', 'series1'),
        ('series2', 'p100_eval.csv', 'series2'), ('series3', 'p100_eval.csv', 'series3'),
        ('gmres_e2c67', 'p100_eval.csv', 'gmres'), ('gmres_e6c67', 'p100_tol6_eval.csv', 'gmres')]


def num(x):
    return float(x) if x not in ('', 'inf', 'nan') else np.inf


def main():
    out = []
    for name, f, m in SETS:
        rs = [r for r in csv.DictReader(open(OUT / f))
              if r['method'] == m and r['status'] in ('ok', 'method_failure') and r.get('rep', '0') == '0']
        if name in ('gmres_e2c67', 'gmres_e6c67'):      # first cap where it converged, the cap-67 rerun for the others
            rf = 'p100_cap67_eval.csv' if name == 'gmres_e2c67' else 'p100_tol6_cap67_eval.csv'
            rerun = {r['case']: r for r in csv.DictReader(open(OUT / rf)) if r['method'] == 'gmres'}
            rs = [rerun[r['case']] if r['status'] == 'method_failure' and r['case'] in rerun else r for r in rs]
        for st in ('all', 'high', 'low'):
            s = [r for r in rs if st == 'all' or mat[r['case']]['stratum'] == st]
            ok = [r for r in s if r['status'] == 'ok']
            e2 = np.array([num(r['eps_2']) for r in ok])
            ei = np.array([num(r['eps_inf']) for r in ok])
            out.append(dict(method=name, stratum=st, n=len(s), with_result=len(ok),
                            e2_median=np.median(e2), e2_q25=np.percentile(e2, 25), e2_q75=np.percentile(e2, 75),
                            e2_max=e2.max(), einf_median=np.median(ei), einf_q25=np.percentile(ei, 25),
                            einf_q75=np.percentile(ei, 75), e2_above_1=int((e2 >= 1).sum())))
    with open(OUT / 'p100_error_summary.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
    print('columns:', list(out[0].keys()))
    for r in out:
        if r['stratum'] == 'all':
            print(r['method'], r['n'], r['with_result'], '%.3g' % r['e2_median'], '%.3g' % r['einf_median'])


if __name__ == '__main__':
    main()
