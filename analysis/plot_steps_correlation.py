"""Pearson correlation of GMRES after 1 to 3 steps and of the series of order 1 to 3 with the dense-LU theta scan.

Reads runs/p100_steps/<case>/k{1,2,3}_0, runs/p100/<case>/series{1,2,3}_0 and the LU run runs/p100/<case>/lu_0
(result.npz and row.json), for every cluster with a finite result. Clusters are split at rho = 1.
Writes figures/steps_correlation.png/.pdf and analysis/out/steps_correlation.csv (case, rho, method, pearson).
"""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / 'runs'
METHODS = [('G1', 'p100_steps', 'k1_0'), ('G2', 'p100_steps', 'k2_0'), ('G3', 'p100_steps', 'k3_0'),
           ('S1', 'p100', 'series1_0'), ('S2', 'p100', 'series2_0'), ('S3', 'p100', 'series3_0')]
COLOR = {'G': '#0F4D92', 'S': '#9A4D8E'}

plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 10, 'axes.linewidth': 1.2,
                     'axes.spines.top': False, 'axes.spines.right': False})


def main():
    rows = []
    for sd in sorted((RUNS / 'p100_steps').iterdir()):
        ref = RUNS / 'p100' / sd.name / 'lu_0'
        if not (ref / 'result.npz').exists():
            continue
        rho = float(json.loads((ref / 'row.json').read_text())['rho'])
        y = np.load(ref / 'result.npz')['cross_section']
        for lab, s, name in METHODS:
            f = RUNS / s / sd.name / name / 'result.npz'
            if not f.exists():
                continue
            a = np.load(f)['cross_section']
            if np.all(np.isfinite(a)):
                rows.append(dict(case=sd.name, rho=rho, method=lab, pearson=float(np.corrcoef(a, y)[0, 1])))
    out = ROOT / 'analysis' / 'out' / 'steps_correlation.csv'
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['case', 'rho', 'method', 'pearson']); w.writeheader(); w.writerows(rows)
    print('columns: case, rho, method, pearson;', len(rows), 'rows')
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.3), sharey=True)
    rng = np.random.default_rng(1)
    for ax, (title, sel) in zip(axes, [(r'weak scatterers, $\rho<1$', lambda r: r < 1),
                                       (r'strong scatterers, $\rho\geq1$', lambda r: r >= 1)]):
        for i, (lab, _, _) in enumerate(METHODS):
            v = np.array([r['pearson'] for r in rows if r['method'] == lab and sel(r['rho'])])
            ax.scatter(i + rng.uniform(-0.25, 0.25, len(v)), v, s=10, color=COLOR[lab[0]], alpha=0.45, lw=0)
            med = np.median(v)
            ax.hlines(med, i - 0.38, i + 0.38, color='k', lw=2)
            ax.text(i, 1.035, '%.2f' % med, ha='center', fontsize=7)
            print(title, lab, len(v), 'median %.4f' % med)
        ax.set_xticks(range(6)); ax.set_xticklabels([m[0] for m in METHODS])
        ax.set_title(title, fontsize=10)
        ax.axhline(1, color='0.6', lw=0.8, ls=':')
        ax.set_ylim(-0.5, 1.08)
    axes[0].set_ylabel('Pearson correlation with dense LU')
    fig.text(0.5, -0.02, 'G1 to G3: GMRES after 1 to 3 steps. S1 to S3: series of order 1 to 3. One dot per cluster; black bar and number: median.',
             ha='center', fontsize=7)
    fig.tight_layout()
    fig.savefig(ROOT / 'figures' / 'steps_correlation.png', dpi=200, bbox_inches='tight')
    fig.savefig(ROOT / 'figures' / 'steps_correlation.pdf', bbox_inches='tight')


if __name__ == '__main__':
    main()
