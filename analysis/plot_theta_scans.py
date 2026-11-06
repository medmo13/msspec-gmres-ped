"""Theta scans of the E16 cases: intensity against theta for LU, GMRES and series 1 to 3.

Reads runs/e16_timing/<case>/<method>_0/{result.npz,row.json} (repetition 0 only).
Under each scan one panel per measure against the same case's LU run: eps_2 and eps_inf (percent), R
(theta-adapted Pendry), R2 and 1 - Pearson r (analysis/extra_metrics.py measures()).
Writes figures/theta_scans.png and figures/theta_scans.pdf.
"""
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
from evaluate import metrics  # noqa: E402
from extra_metrics import measures  # noqa: E402
from cluster_ase import draw  # noqa: E402

RUN = ROOT / 'runs' / 'e16_timing'
CASES = ['H11', 'H10', 'H08', 'L45']
METHODS = ['lu', 'gmres', 'series1', 'series2', 'series3']
LABEL = {'lu': 'LU (reference)', 'gmres': 'GMRES, tol 1e-2', 'series1': 'Series 1',
         'series2': 'Series 2', 'series3': 'Series 3'}
COLOR = {'lu': '#B64342', 'gmres': '#0F4D92', 'series1': '#42949E', 'series2': '#8BCF8B', 'series3': '#9A4D8E'}
STYLE = {'lu': '-', 'gmres': '--', 'series1': ':', 'series2': ':', 'series3': ':'}

PANELS = {'eps2': r'$\varepsilon_2$ (%)', 'epsinf': r'$\varepsilon_\infty$ (%)', 'R': '$R$',
          'R2': '$R_2$', '1-r': '$1-r$'}

plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 10, 'axes.linewidth': 1.2,
                     'axes.spines.top': False, 'axes.spines.right': False})


def load(case, method):
    d = RUN / case / (method + '_0')
    z = np.load(d / 'result.npz')
    row = json.loads((d / 'row.json').read_text())
    return z, row


def describe(row):
    # Material name from run_tag: r_<id>_s_<id>_ke_<energy>_<material text>_o_<point group>_phi_<phi>
    tag = row['run_tag']
    mat = tag.split('_ke_')[1].split('_o_')[0].split('_', 1)[1].replace('_', ' ')
    return r'%s, %s atoms, %s eV, $\rho$ = %.2f' % (mat, row['natoms'], row['ke'], float(row['rho']))


def main(cases=CASES, run=None, stem='theta_scans'):
    global RUN
    if run:
        RUN = ROOT / 'runs' / run
    CASES_ = cases
    fig = plt.figure(figsize=(12, 3.9 * len(CASES_)))
    outer = fig.add_gridspec(len(CASES_), 1, hspace=0.35, top=0.99, bottom=0.045)
    first = None
    for i, case in enumerate(CASES_):
        inner = outer[i].subgridspec(1, 2, width_ratios=[1.3, 1.0], wspace=0.22)
        right = inner[1].subgridspec(2, 1, height_ratios=[1.0, 2.3], hspace=0.2)
        cl = right[0].subgridspec(1, 3, wspace=0.05)
        bars = right[1].subgridspec(2, 3, hspace=0.9, wspace=0.75)
        axc = [fig.add_subplot(cl[0, j]) for j in range(3)]
        draw(axc, case)
        z_ref, row_ref = load(case, 'lu')
        ax = fig.add_subplot(inner[0])
        theta = z_ref['theta']
        ax.plot(theta, z_ref['cross_section'], color=COLOR['lu'], lw=2.2, label=LABEL['lu'])
        vals = {k: [] for k in PANELS}
        for m in METHODS[1:]:
            z, _ = load(case, m)
            ax.plot(theta, z['cross_section'], color=COLOR[m], ls=STYLE[m], lw=1.4, label=LABEL[m])
            e = metrics(z['cross_section'], z_ref['cross_section'])
            _, r_, r2_, p_ = measures(theta, z_ref['cross_section'], z['cross_section'])
            vals['eps2'].append(100 * e['eps_2']); vals['epsinf'].append(100 * e['eps_inf'])
            vals['R'].append(r_); vals['R2'].append(r2_); vals['1-r'].append(max(1 - p_, 1e-16))
        ax.set_ylabel('Intensity (arb. units)'); ax.set_yscale('log'); ax.set_ylim(bottom=ax.get_ylim()[1] * 1e-5)
        ax.set_title(describe(row_ref), loc='left', fontsize=10)
        ax.set_xlabel('Theta (degrees)')
        axl = fig.add_subplot(bars[1, 2]); axl.axis('off')
        axl.legend(*ax.get_legend_handles_labels(), frameon=False, fontsize=7, loc='center left')
        for j, (k, title) in enumerate(PANELS.items()):
            axm = fig.add_subplot(bars[j // 3, j % 3])
            axm.bar(np.arange(4), vals[k], color=[COLOR[m] for m in METHODS[1:]], width=0.75)
            axm.set_yscale('log')
            axm.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            axm.set_xticks(np.arange(4))
            axm.set_xticklabels(['G', 'S1', 'S2', 'S3'], fontsize=8)
            axm.set_title(title, fontsize=9)
            if k in ('eps2', 'epsinf'):
                axm.set_ylabel('%', fontsize=8)
            axm.tick_params(axis='y', labelsize=7)
    fig.text(0.5, 0.0, 'G: GMRES, tol 1e-2.  S1 to S3: series of order 1 to 3.', ha='center', fontsize=9)
    out = ROOT / 'figures'
    fig.savefig(out / (stem + '.png'), dpi=200, bbox_inches='tight')
    fig.savefig(out / (stem + '.pdf'), bbox_inches='tight')
    print('wrote figures/%s.png and .pdf' % stem)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        main(sys.argv[3].split(','), sys.argv[1], sys.argv[2])
    else:
        main()
