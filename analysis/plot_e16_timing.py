"""Figure of the E16 timing and memory study (A2.8): median wall time and peak RSS.

Reads analysis/out/e16_timing.csv and analysis/out/e16_timing.npz (from timing_e16.py).
Writes figures/e16_timing.png and figures/e16_timing.pdf.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
CASES = ['H11', 'H10', 'H08', 'L45']
CASE_LABEL = {'H11': 'BaTiO$_3$\n118 atoms', 'H10': 'NaS\n96 atoms', 'H08': 'ReMg\n96 atoms', 'L45': 'InTe\n51 atoms'}
METHODS = ['lu', 'gmres', 'series1', 'series2', 'series3']
METHOD_LABEL = {'lu': 'LU (reference)', 'gmres': 'GMRES, tol 1e-2', 'series1': 'Series 1',
                'series2': 'Series 2', 'series3': 'Series 3'}
COLOR = {'lu': '#B64342', 'gmres': '#0F4D92', 'series1': '#42949E', 'series2': '#8BCF8B', 'series3': '#9A4D8E'}

plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 11, 'axes.linewidth': 1.5,
                     'axes.spines.top': False, 'axes.spines.right': False})


def main():
    rows = {(r['case'], r['method']): r for r in csv.DictReader(open(ROOT / 'analysis/out/e16_timing.csv'))}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))
    width = 0.16
    x = np.arange(len(CASES))
    for k, m in enumerate(METHODS):
        off = (k - 2) * width
        med = [float(rows[(c, m)]['wall_median_s']) for c in CASES]
        lo = [float(rows[(c, m)]['wall_median_s']) - float(rows[(c, m)]['wall_min_s']) for c in CASES]
        hi = [float(rows[(c, m)]['wall_max_s']) - float(rows[(c, m)]['wall_median_s']) for c in CASES]
        ax1.bar(x + off, med, width, yerr=[lo, hi], color=COLOR[m], label=METHOD_LABEL[m],
                capsize=3, error_kw={'linewidth': 1})
        rss = [float(rows[(c, m)]['peak_rss_mib']) for c in CASES]
        ax2.bar(x + off, rss, width, color=COLOR[m], label=METHOD_LABEL[m])
    ax1.set_yscale('log')
    ax1.set_ylabel('Wall time per calculation (s, median of 3)')
    ax1.set_title('(a) Time per calculation', loc='left')
    ax2.set_yscale('log')
    ax2.set_ylabel('Peak memory (MiB, largest of 3)')
    ax2.set_title('(b) Peak memory', loc='left')
    for ax in (ax1, ax2):
        ax.set_xticks(x)
        ax.set_xticklabels([CASE_LABEL[c] for c in CASES])
        ax.grid(axis='y', which='major', color='#CFCECE', linewidth=0.6)
        ax.set_axisbelow(True)
    ax1.legend(frameon=False, fontsize=9, loc='upper right')
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    out = ROOT / 'figures'
    fig.savefig(out / 'e16_timing.png', dpi=300)
    fig.savefig(out / 'e16_timing.pdf')
    print('wrote figures/e16_timing.png and .pdf')


if __name__ == '__main__':
    main()
