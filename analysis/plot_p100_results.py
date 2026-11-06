"""P100 results: error of every case against dense LU for GMRES (three tolerances) and the series, and error against rho.

Reads analysis/out/p100_eval.csv (evaluate.py), runs/p100/rows.csv (rho, atoms).
Writes figures/p100_results.png and .pdf.
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

COLOR = {'gmres': '#0F4D92', 'series1': '#42949E', 'series2': '#8BCF8B', 'series3': '#9A4D8E'}
LABEL = {'gmres': 'GMRES, tol 1e-2', 'series1': 'Series 1', 'series2': 'Series 2', 'series3': 'Series 3'}
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 12, 'axes.linewidth': 1.2,
                     'axes.spines.top': False, 'axes.spines.right': False})


def load(f, m):
    return [r for r in csv.DictReader(open(ROOT / 'analysis/out' / f))
            if r['method'] == m and r['status'] == 'ok' and r.get('rep', '0') == '0'
            and r['eps_2'] not in ('', 'inf', 'nan')]


EXTRA = {'p100_eval.csv': 'p100_cap67_eval.csv', 'p100_tol6_eval.csv': 'p100_tol6_cap67_eval.csv'}


def load_all(f, m):
    # GMRES at 1e-2 and 1e-6: add the clusters that converged in the cap-67 rerun
    rs = load(f, m)
    if m == 'gmres' and f in EXTRA:
        have = {r['case'] for r in rs}
        rs += [r for r in load(EXTRA[f], m) if r['case'] not in have]
    return rs


def main():
    mat = {r['mat_id']: r for r in csv.DictReader(open(ROOT / 'platform/materials_v2.csv'))}
    rho_any = {r['mat_id']: float(r['rho']) for r in csv.DictReader(open(ROOT / 'platform/materials_v2.csv'))}
    sets = [('GMRES\ntol 1e-2', 'p100_eval.csv', 'gmres', '#0F4D92'), ('GMRES\ntol 1e-3', 'p100_tol3_eval.csv', 'gmres', '#3775BA'),
            ('GMRES\ntol 1e-6', 'p100_tol6_eval.csv', 'gmres', '#06294F'), ('Series\norder 1', 'p100_eval.csv', 'series1', '#42949E'),
            ('Series\norder 2', 'p100_eval.csv', 'series2', '#8BCF8B'), ('Series\norder 3', 'p100_eval.csv', 'series3', '#9A4D8E')]
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.8), gridspec_kw={'width_ratios': [1.15, 1, 1.0]})

    # Left: error of every case against dense LU, per method (box: quartiles; points: cases).
    rng = np.random.default_rng(1)
    for k, (lab, f, m, col) in enumerate(sets):
        rs = load_all(f, m)
        e = np.array([max(float(r['eps_2']), 1e-8) for r in rs])
        ax1.boxplot([e], positions=[k], widths=0.55, showfliers=False, patch_artist=True,
                    boxprops={'facecolor': col, 'alpha': 0.35, 'edgecolor': col}, medianprops={'color': 'black', 'lw': 1.6},
                    whiskerprops={'color': col}, capprops={'color': col})
        ax1.scatter(k + rng.uniform(-0.18, 0.18, len(e)), e, s=8, color=col, alpha=0.7, lw=0)
    ax1.axhline(1.0, color='grey', ls=':', lw=1); ax1.text(-0.4, 1.25, 'error = signal', ha='left', fontsize=8, color='grey')
    ax1.set_yscale('log'); ax1.set_ylim(1e-8 * 1.5, 1e3)
    ax1.set_xticks(range(len(sets))); ax1.set_xticklabels([x[0] for x in sets], fontsize=10)
    ax1.set_ylabel('relative L2 error, against dense LU')
    ax1.set_title('(a) Error of every case, per method', loc='left')

    # Right: error against rho, one point per case.
    for lab, f, m, col in sets[:1] + sets[3:]:
        pts = [(rho_any[r['case']], float(r['eps_2'])) for r in load_all(f, m)]
        p = np.array(pts)
        ax2.scatter(p[:, 0], np.maximum(p[:, 1], 1e-6), s=16, color=col, label=lab.replace('\n', ' '), alpha=0.8)
    ax2.axvline(1.0, color='grey', ls=':', lw=1)
    ax2.set_xscale('log'); ax2.set_yscale('log')
    ax2.set_xlabel('rho (spectral radius of G0 T); dotted line: rho = 1')
    ax2.set_ylabel('relative L2 error against dense LU')
    ax2.set_title('(b) Error against rho', loc='left')
    ax2.legend(frameon=False, fontsize=8, loc='lower right')

    # Third panel: 1 - r (Pearson correlation of the scan with the dense-LU scan) of every case, same methods
    pe = {}
    for r in csv.DictReader(open(ROOT / 'analysis/out/p100_extra_metrics.csv')):
        pe.setdefault(r['label'], {})[r['case']] = 1 - float(r['pearson'])
    labs = ['gmres_e2c67', 'gmres_e3', 'gmres_e6c67', 'series1', 'series2', 'series3']
    for k, ((lab, f, m, col), key) in enumerate(zip(sets, labs)):
        v = np.array([max(x, 1e-14) for x in pe[key].values()])
        ax3.boxplot([v], positions=[k], widths=0.55, showfliers=False, patch_artist=True,
                    boxprops={'facecolor': col, 'alpha': 0.35, 'edgecolor': col}, medianprops={'color': 'black', 'lw': 1.6},
                    whiskerprops={'color': col}, capprops={'color': col})
        ax3.scatter(k + rng.uniform(-0.18, 0.18, len(v)), v, s=8, color=col, alpha=0.7, lw=0)
    ax3.set_yscale('log')
    ax3.set_xticks(range(len(sets))); ax3.set_xticklabels([x[0] for x in sets], fontsize=10)
    ax3.set_ylabel('1 - r (Pearson), against dense LU')
    ax3.set_title('(c) Pearson correlation', loc='left')
    fig.tight_layout()
    fig.savefig(ROOT / 'figures/p100_results.png', dpi=200)
    fig.savefig(ROOT / 'figures/p100_results.pdf')
    print('wrote figures/p100_results.png and .pdf')


if __name__ == '__main__':
    main()
