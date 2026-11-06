"""Series of order 4 on a subset (D43): error against dense LU and wall time of the series of order 1 to 4 and of GMRES (tol 1e-2).

Order 4 was run on the first 13 clusters of the test set only (runs/series4); it was stopped because it is much slower than orders 1 to 3.
Orders 1 to 3 and GMRES come from runs/p100 (GMRES: runs/p100_cap67 where it replaces a failed first run). Wall times are from different
sessions, not a paired timing. Writes analysis/out/series4_subset.csv and figures/series4_subset.png/.pdf.
"""
import csv
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from evaluate import metrics  # noqa: E402
R = ROOT / 'runs'


def last(s, m):
    d = {}
    for x in csv.DictReader(open(R / s / 'rows.csv')):
        if x['method'] == m:
            d[x['case']] = (s, x)
    return d


def cs(s, c, m):
    return np.load(R / s / c / (m + '_0') / 'result.npz')['cross_section']


s4 = last('series4', 's4') if False else last('series4', 'series4')
ser = {k: last('p100', 'series%d' % k) for k in (1, 2, 3)}
gm = last('p100', 'gmres')
gm.update({c: v for c, v in last('p100_cap67', 'gmres').items()})
rows = []
for c in sorted(s4):
    lu = cs('p100', c, 'lu')
    rec = dict(case=c, natoms=int(s4[c][1]['natoms']), rho=float(s4[c][1]['rho']))
    for k in (1, 2, 3):
        s, x = ser[k][c]
        rec['eps2_s%d' % k] = 100 * metrics(cs(s, c, 'series%d' % k), lu)['eps_2']
        rec['wall_s%d' % k] = float(x['wall_s'])
    rec['eps2_s4'] = 100 * metrics(cs('series4', c, 'series4'), lu)['eps_2']
    rec['wall_s4'] = float(s4[c][1]['wall_s'])
    s, x = gm[c]
    rec['eps2_g'] = 100 * metrics(cs(s, c, 'gmres'), lu)['eps_2']
    rec['wall_g'] = float(x['wall_s'])
    rows.append(rec)
out = ROOT / 'analysis/out/series4_subset.csv'
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print('columns:', list(rows[0]))
for r in rows:
    print(r['case'], r['natoms'], '%.2f' % r['rho'], ' '.join('%.3g' % r[k] for k in r if k.startswith(('eps2', 'wall'))))

fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
cm = plt.get_cmap('viridis')
rmin, rmax = min(r['rho'] for r in rows), max(r['rho'] for r in rows)
for r in rows:
    col = cm((r['rho'] - rmin) / (rmax - rmin))
    ax[0].plot([1, 2, 3, 4], [r['eps2_s%d' % k] for k in (1, 2, 3, 4)], 'o-', color=col, lw=1, ms=3)
    ax[1].plot([1, 2, 3, 4], [r['wall_s%d' % k] for k in (1, 2, 3, 4)], 'o-', color=col, lw=1, ms=3)
    ax[0].plot([5.2], [r['eps2_g']], 's', color=col, ms=4)
    ax[1].plot([5.2], [r['wall_g']], 's', color=col, ms=4)
for a in ax:
    a.set_yscale('log')
    a.set_xticks([1, 2, 3, 4, 5.2])
    a.set_xticklabels(['1', '2', '3', '4', 'GMRES'])
    a.set_xlabel('series order (squares: GMRES, tolerance $10^{-2}$)')
    for sp in ('top', 'right'):
        a.spines[sp].set_visible(False)
ax[0].axhline(10, color='k', ls=':', lw=1)
ax[0].set_ylabel(r'$\varepsilon_2$ against dense LU (%)')
ax[1].set_ylabel('wall time per scan (s)')
ax[0].set_title('(a) error')
ax[1].set_title('(b) time')
sm = plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(rmin, rmax))
fig.colorbar(sm, ax=ax, label=r'$\rho$', fraction=0.03, pad=0.02)
fig.savefig(ROOT / 'figures/series4_subset.png', dpi=200, bbox_inches='tight')
fig.savefig(ROOT / 'figures/series4_subset.pdf', bbox_inches='tight')
