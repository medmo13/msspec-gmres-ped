"""Worst clusters (lowest Pearson r against dense LU) for GMRES stopped after 1 to 4 steps.

Reads runs/p100_e3/<case>/{lu_0,k1_0..k4_0}. For each k the cluster with the lowest r is taken (duplicates merged).
Writes platform/work/worst_clusters.png and prints the plotted r values.
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / 'runs' / 'p100_e3'
KS = ['k1_0', 'k2_0', 'k3_0', 'k4_0']
COL = ['#CFE0F2', '#8DB5DE', '#3775BA', '#06294F']


def load(c, m):
    return np.load(R / c / m / 'result.npz'), json.loads((R / c / m / 'row.json').read_text())


r = {}
for c in sorted(p.name for p in R.iterdir()):
    if not (R / c / 'lu_0' / 'result.npz').exists():
        continue
    a = load(c, 'lu_0')[0]['cross_section']
    for k in KS:
        r[(c, k)] = np.corrcoef(a, load(c, k)[0]['cross_section'])[0, 1]
worst = []
for k in KS:
    c = min((c for c, kk in r if kk == k), key=lambda c: r[(c, k)])
    if c not in worst:
        worst.append(c)
fig, axes = plt.subplots(len(worst), 2, figsize=(12, 3.3 * len(worst)), squeeze=False)
for i, c in enumerate(worst):
    lu, row = load(c, 'lu_0')
    mat = row['run_tag'].split('_ke_')[1].split('_o_')[0].split('_', 1)[1].replace('_', ' ')
    th = lu['theta']
    for ax, log in zip(axes[i], (True, False)):
        ax.plot(th, lu['cross_section'], color='#B64342', lw=2.4, label='dense LU')
        for k, col in zip(KS, COL):
            ax.plot(th, load(c, k)[0]['cross_section'], color=col, lw=1.3, ls='--',
                    label='%s step%s, r=%.2f' % (k[1], '' if k[1] == '1' else 's', r[(c, k)]))
        ax.set_xlabel('Theta (degrees)'); ax.set_ylabel('Intensity (arb. units)')
        if log:
            ax.set_yscale('log'); ax.set_ylim(lu['cross_section'].max() * 1e-4, None)
    axes[i][0].set_title('%s, %s atoms, %s eV, rho %.2f' % (mat, row['natoms'], row['ke'], float(row['rho'])), loc='left', fontsize=10)
    axes[i][0].legend(frameon=False, fontsize=7, loc='lower left')
    print(c, mat, row['natoms'], row['ke'], row['rho'], {k: round(r[(c, k)], 3) for k in KS})
for ax in axes.ravel():
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
fig.tight_layout()
fig.savefig(ROOT / 'platform/work/worst_clusters.png', dpi=130)
