"""Time and peak-memory ratio of dense LU to GMRES (tol 1e-2) over the clusters with a dense-LU reference.

Reads runs/p100/rows.csv (LU and GMRES, one run each) and runs/p100_cap67/rows.csv (GMRES rerun with cap 67 for
the clusters that hit the first cap of 12). Writes analysis/out/p100_cost.csv, figures/p100_cost.png and .pdf.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
TIMED = {'H11', 'H10', 'H08', 'L45'}
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 11, 'axes.linewidth': 1.5,
                     'axes.spines.top': False, 'axes.spines.right': False})


def rd(p):
    return list(csv.DictReader(open(ROOT / p)))


def main():
    p = {(r['case'], r['method']): r for r in rd('runs/p100/rows.csv')}
    cap = {r['case']: r for r in rd('runs/p100_cap67/rows.csv') if r['method'] == 'gmres'}
    out = []
    for (c, m), r in sorted(p.items()):
        if m != 'gmres':
            continue
        lu = p.get((c, 'lu'))
        g = r if r['run_status'] == 'ok' else cap.get(c)
        if lu is None or lu['run_status'] != 'ok' or g is None or g['run_status'] != 'ok':
            continue
        out.append({'case': c, 'natoms': int(g['natoms']), 'rho': float(g['rho']),
                    'lu_wall_s': float(lu['wall_s']), 'gmres_wall_s': float(g['wall_s']),
                    'lu_rss_mib': float(lu['peak_rss_mib']), 'gmres_rss_mib': float(g['peak_rss_mib']),
                    'time_ratio': float(lu['wall_s']) / float(g['wall_s']),
                    'mem_ratio': float(lu['peak_rss_mib']) / float(g['peak_rss_mib'])})
    with open(ROOT / 'analysis/out/p100_cost.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    n = np.array([o['natoms'] for o in out])
    strong = np.array([o['rho'] > 1 for o in out])
    timed = np.array([o['case'] in TIMED for o in out])
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, key, lab, ttl in ((axs[0], 'time_ratio', 'Time of dense LU / time of GMRES', '(a) Time'),
                              (axs[1], 'mem_ratio', 'Peak memory of dense LU / of GMRES', '(b) Peak memory')):
        y = np.array([o[key] for o in out])
        ax.scatter(n[~strong], y[~strong], s=26, color='#0F4D92', alpha=0.75, label='weak ($\\rho<1$)')
        ax.scatter(n[strong], y[strong], s=26, color='#B64342', alpha=0.75, label='strong ($\\rho>1$)')
        ax.scatter(n[timed], y[timed], s=90, facecolors='none', edgecolors='black', linewidths=1.3, label='timed emitters')
        ax.axhline(1, color='#555555', linewidth=1, linestyle='--')
        ax.set_yscale('log')
        ax.set_xlabel('Atoms in the cluster')
        ax.set_ylabel(lab)
        ax.set_title(ttl, loc='left')
        ax.grid(axis='y', which='major', color='#CFCECE', linewidth=0.6)
        ax.set_axisbelow(True)
    axs[0].legend(frameon=False, fontsize=9, loc='upper left')
    fig.tight_layout()
    fig.savefig(ROOT / 'figures/p100_cost.png', dpi=300)
    fig.savefig(ROOT / 'figures/p100_cost.pdf')
    print('wrote figures/p100_cost.png, .pdf and analysis/out/p100_cost.csv', len(out), 'clusters')


if __name__ == '__main__':
    main()
