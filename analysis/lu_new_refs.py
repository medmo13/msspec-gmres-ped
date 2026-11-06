"""Error of GMRES and the series against the new dense-LU references on the clusters that earlier
runs had marked too large for LU (4 large LiF-type clusters; 3 strong scatterers). Writes analysis/out/lu_new_refs.csv."""
import csv, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from evaluate import metrics
R = Path(__file__).resolve().parents[1] / 'runs'
rows = []
def cs(p):
    return np.load(p / 'result.npz')['cross_section']
def rowinfo(setname, case, method):
    for x in csv.DictReader(open(R / setname / 'rows.csv')):
        if x['case'] == case and x['method'] == method:
            return x
for case, lu_set, srcs in [(c, 'big_lu', {'gmres': 'big_ref', 'series1': 'big', 'series2': 'big', 'series3': 'big'}) for c in ['B01', 'B02', 'B03', 'B04']] + \
                          [(c, 'p100_lu_excluded', {m: 'p100' for m in ['gmres', 'series1', 'series2', 'series3']}) for c in ['H19', 'H37', 'H49']]:
    ref = cs(R / lu_set / case / 'lu_0')
    li = rowinfo(lu_set, case, 'lu')
    for m, s in srcs.items():
        a = cs(R / s / case / (m + '_0'))
        met = metrics(a, ref)
        gi = rowinfo(s, case, m)
        rows.append(dict(case=case, natoms=li['natoms'], ke=li['ke'], rho=li['rho'], method=m, src=s,
                         tol=gi['tol'], eps_2=met['eps_2'], eps_inf=met['eps_inf'], nonfinite=met['n_nonfinite'],
                         lu_wall_s=li['wall_s'], lu_peak_mib=li['peak_rss_mib'], gm_wall_s=gi['wall_s']))
out = Path(__file__).parent / 'out' / 'lu_new_refs.csv'
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows:
    print(r['case'], r['natoms'], 'rho=%.2f' % float(r['rho']), r['method'], 'tol', r['tol'], 'eps2=%.2e epsinf=%.2e' % (r['eps_2'], r['eps_inf']), 'nonfin', r['nonfinite'])
print('wrote', out)
