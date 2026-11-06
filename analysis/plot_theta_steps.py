"""Theta scans for GMRES stopped after 1 to 4 steps, and for the series of order 1 to 3, against dense LU.

Reads runs/p100_e3/<case>/{lu_0,k1_0..k4_0} (GMRES after k steps) and runs/p100/<case>/{lu_0,series1_0..series3_0}.
Same four emitters as plot_theta_scans.py. Per emitter: GMRES patterns, series patterns, and eps_2 and 1 - r
bars for all seven methods (eps_2 as in extra_metrics.measures, percent).
Writes figures/theta_scans_steps.png and .pdf and analysis/out/theta_scans_steps.csv (the plotted numbers).
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

CASES = ['H11', 'H10', 'H08', 'L45']
G = [('GMRES, 1 step', 'p100_e3', 'k1_0', '#CFE0F2'), ('GMRES, 2 steps', 'p100_e3', 'k2_0', '#8DB5DE'),
     ('GMRES, 3 steps', 'p100_e3', 'k3_0', '#3775BA'), ('GMRES, 4 steps', 'p100_e3', 'k4_0', '#06294F')]
S = [('Series, order 1', 'p100', 'series1_0', '#42949E'), ('Series, order 2', 'p100', 'series2_0', '#8BCF8B'),
     ('Series, order 3', 'p100', 'series3_0', '#9A4D8E')]
REF = '#B64342'
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 10, 'axes.linewidth': 1.2,
                     'axes.spines.top': False, 'axes.spines.right': False})


def load(s, case, m):
    d = ROOT / 'runs' / s / case / m
    return np.load(d / 'result.npz'), json.loads((d / 'row.json').read_text())


def describe(row):
    mat = row['run_tag'].split('_ke_')[1].split('_o_')[0].split('_', 1)[1].replace('_', ' ')
    return '%s, %s atoms, %s eV' % (mat, row['natoms'], row['ke'])


def main():
    fig = plt.figure(figsize=(14, 3.9 * len(CASES)))
    outer = fig.add_gridspec(len(CASES), 1, hspace=0.35, top=0.99, bottom=0.04)
    table = []
    for i, case in enumerate(CASES):
        inner = outer[i].subgridspec(1, 3, width_ratios=[1, 1, 0.9], wspace=0.3)
        right = inner[2].subgridspec(2, 1, height_ratios=[1.0, 2.3], hspace=0.25)
        draw([fig.add_subplot(right[0].subgridspec(1, 3, wspace=0.05)[0, j]) for j in range(3)], case)
        a1, a2, a3 = fig.add_subplot(inner[0]), fig.add_subplot(inner[1]), fig.add_subplot(right[1])
        ref_s = {s: load(s, case, 'lu_0')[0] for s in ('p100', 'p100_e3')}
        theta = ref_s['p100']['theta']
        title = describe(load('p100', case, 'lu_0')[1])
        for ax, group in ((a1, G), (a2, S)):
            ax.plot(theta, ref_s['p100']['cross_section'], color=REF, lw=2.4, label='LU (reference)')
            for lab, s, m, col in group:
                ax.plot(theta, load(s, case, m)[0]['cross_section'], color=col, lw=1.3, ls='--', label=lab)
            ax.set_xlabel('Theta (degrees)'); ax.set_ylabel('Intensity (arb. units)'); ax.set_yscale('log')
        lo, hi = ref_s['p100']['cross_section'].max() * 1e-4, max(a1.get_ylim()[1], a2.get_ylim()[1])
        a1.set_ylim(lo, hi); a2.set_ylim(lo, hi)
        a1.set_title(title, loc='left', fontsize=10)
        eps, om = [], []
        for lab, s, m, col in G + S:
            z = load(s, case, m)[0]
            e, r, r2, p = measures(theta, ref_s[s]['cross_section'], z['cross_section'])
            eps.append(e); om.append(max(1 - p, 1e-16))
            table.append(dict(case=case, method=lab, eps2_pct=e, one_minus_r=1 - p))
        x = np.arange(len(eps))
        a3.bar(x, eps, color=[c for *_, c in G + S], width=0.75)
        a3.set_yscale('log'); a3.set_ylabel(r'$\varepsilon_2$ against LU (%)')
        a3.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        a3.set_xticks(x); a3.set_xticklabels(['G1', 'G2', 'G3', 'G4', 'S1', 'S2', 'S3'], fontsize=8)
        a3.set_title('error of each method', loc='left', fontsize=10)
        a3r = a3.twinx(); a3r.spines['right'].set_visible(True)         # Pearson 1 - r on the right axis
        a3r.plot(x, om, 'D', color='k', ms=4); a3r.set_yscale('log'); a3r.set_ylabel(r'$1-r$ (diamonds)')
        a3r.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        if i == 0:
            a1.legend(frameon=False, fontsize=7, loc='lower left'); a2.legend(frameon=False, fontsize=7, loc='lower left')
    fig.text(0.5, 0.0, 'G1 to G4: GMRES stopped after 1 to 4 steps.  S1 to S3: series of order 1 to 3.  Bars: $\\varepsilon_2$; diamonds: $1-r$ (Pearson r against dense LU).', ha='center', fontsize=9)
    fig.savefig(ROOT / 'figures/theta_scans_steps.png', dpi=200, bbox_inches='tight'); fig.savefig(ROOT / 'figures/theta_scans_steps.pdf', bbox_inches='tight')
    with open(ROOT / 'analysis/out/theta_scans_steps.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0].keys())); w.writeheader(); w.writerows(table)
    print('columns:', list(table[0].keys()))
    for t in table:
        print(t['case'], t['method'], '%.3g' % t['eps2_pct'])


if __name__ == '__main__':
    main()
