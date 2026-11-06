"""Large low-rho clusters (about 270 atoms, above the dense-LU memory limit): GMRES at tolerance 1e-3 and the series of
order 1 to 3 against dense LU (runs/big_lu, D39).

Reads runs/big_lu/<case>/lu_0 (reference; run_set.py D39) and runs/big_ref (row info) and runs/big/<case>/{gmres_0,series1_0..series3_0}.
Left: polar-angle scans. Right: eps_2 (percent) of each method against the reference (analysis/extra_metrics.measures).
Writes figures/big_clusters.png/.pdf and analysis/out/big_clusters.csv (case, material, atoms, rho, method, eps_2, 1-r).
"""
import csv
import json
import sys
from pathlib import Path
import numpy as np
import matplotlib
import matplotlib.ticker
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from extra_metrics import measures  # noqa: E402
from cluster_ase import draw  # noqa: E402

CASES = ['B01', 'B02', 'B03', 'B04']
SHOWN = ['B01', 'B03', 'B04']          # B02 (LiF at 787 eV) differs from B01 only in energy; it stays in the CSV, not in the figure
M = [('GMRES, tol 1e-3', 'big', 'gmres_0', '#0F4D92'), ('Series, order 1', 'big', 'series1_0', '#42949E'),
     ('Series, order 2', 'big', 'series2_0', '#8BCF8B'), ('Series, order 3', 'big', 'series3_0', '#9A4D8E')]
REF = '#B64342'
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 10, 'axes.linewidth': 1.2,
                     'axes.spines.top': False, 'axes.spines.right': False})


def load(s, c, m):
    d = ROOT / 'runs' / s / c / m
    return np.load(d / 'result.npz'), json.loads((d / 'row.json').read_text())


def main():
    fig = plt.figure(figsize=(12, 3.9 * len(SHOWN)))
    outer = fig.add_gridspec(len(SHOWN), 1, hspace=0.35, top=0.99, bottom=0.075)
    table = []
    for c in CASES:                      # the CSV has every case
        ref, _ = load('big_lu', c, 'lu_0'); _, row = load('big_ref', c, 'gmres_0')
        mat = row['run_tag'].split('_ke_')[1].split('_o_')[0].split('_', 1)[1].replace('_', ' ')
        for lab, ss, m, col in M:
            z, _ = load(ss, c, m)
            e, r, r2, p = measures(ref['theta'], ref['cross_section'], z['cross_section'])
            table.append(dict(case=c, material=mat, natoms=row['natoms'], rho=row['rho'], method=lab, eps2_pct=e, one_minus_r=1 - p))
    for i, c in enumerate(SHOWN):        # the figure has the cases in SHOWN
        ref, _ = load('big_lu', c, 'lu_0'); _, row = load('big_ref', c, 'gmres_0')
        theta, ref_cs = ref['theta'], ref['cross_section']
        mat = row['run_tag'].split('_ke_')[1].split('_o_')[0].split('_', 1)[1].replace('_', ' ')
        inner = outer[i].subgridspec(1, 2, width_ratios=[1.3, 1.0], wspace=0.22)
        right = inner[1].subgridspec(2, 1, height_ratios=[1.0, 2.3], hspace=0.2)
        draw([fig.add_subplot(right[0].subgridspec(1, 3, wspace=0.05)[0, j]) for j in range(3)], c, row['run_tag'])
        a1, a2 = fig.add_subplot(inner[0]), fig.add_subplot(right[1])
        a1.plot(theta, ref_cs, color=REF, lw=2.4, label='dense LU (reference)')
        eps, om = [], []
        for lab, ss, m, col in M:
            z, _ = load(ss, c, m)
            a1.plot(theta, z['cross_section'], color=col, lw=1.3, ls='--', label=lab)
            e, r, r2, p = measures(theta, ref_cs, z['cross_section'])
            eps.append(e); om.append(max(1 - p, 1e-16))
        a1.set_xlabel('Theta (degrees)'); a1.set_ylabel('Intensity (arb. units)'); a1.set_yscale('log'); a1.set_ylim(bottom=ref_cs.max() * 1e-5)
        a1.set_title('%s, %s atoms, rho %.2f, %s eV' % (mat, row['natoms'], float(row['rho']), row['ke']), loc='left', fontsize=10)
        if i == 0:
            handles = a1.get_legend_handles_labels()
        a2.bar(np.arange(4), eps, color=[m_[-1] for m_ in M], width=0.75)
        a2.set_yscale('log'); a2.set_ylabel(r'$\varepsilon_2$ against reference (%)')
        a2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        a2.set_xticks(np.arange(4)); a2.set_xticklabels(['G', 'S1', 'S2', 'S3'], fontsize=9)
        a2.set_title('error of each method', loc='left', fontsize=10)
        a2r = a2.twinx(); a2r.spines['right'].set_visible(True)          # Pearson 1 - r on the right axis
        a2r.plot(np.arange(4), om, 'D', color='k', ms=4); a2r.set_yscale('log'); a2r.set_ylabel(r'$1-r$ (diamonds)')
        a2r.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    fig.legend(*handles, frameon=False, fontsize=8, ncol=5, loc='lower center', bbox_to_anchor=(0.5, 0.02))
    fig.text(0.5, 0.0, 'G: GMRES, tolerance 1e-3.  S1 to S3: series of order 1 to 3.  Reference: dense LU.  Bars: $\\varepsilon_2$; diamonds: $1-r$ (Pearson).',
             ha='center', fontsize=9)
    fig.savefig(ROOT / 'figures/big_clusters.png', dpi=200, bbox_inches='tight'); fig.savefig(ROOT / 'figures/big_clusters.pdf', bbox_inches='tight')
    with open(ROOT / 'analysis/out/big_clusters.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0].keys())); w.writeheader(); w.writerows(table)
    print('columns:', list(table[0].keys()))
    for t in table:
        print(t['case'], t['material'], t['natoms'], 'rho %.3f' % float(t['rho']), t['method'], 'eps2 %.3g %%' % t['eps2_pct'])


if __name__ == '__main__':
    main()
