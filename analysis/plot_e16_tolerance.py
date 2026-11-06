"""Figure and summary npz for the E16 step and tolerance study (D31).

Reads the three evaluate CSVs (platform/work/<set>_eval.csv, from analysis/evaluate.py) and
writes analysis/out/e16_tolerance.npz and figures/e16_tolerance.{png,pdf}.
Colours follow the figures4papers PALETTE (blue for GMRES, red for baselines).
"""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
PAL = {'blue_main': '#0F4D92', 'blue_secondary': '#3775BA', 'red_strong': '#B64342',
       'red_2': '#E9A6A1', 'neutral': '#CFCECE', 'teal': '#42949E'}
CASES = [('H11', 'barium titanate'), ('H10', 'sodium sulfide'), ('H08', 'rhenium-magnesium'), ('L45', 'indium telluride')]
STEPS = ['k0', 'k1', 'k2', 'k3']
GM = [('gmres_t2', 'e16_tol2', 'gmres', '1e-2'), ('gmres_t3', 'e16_tol3', 'gmres', '1e-3'),
      ('gmres_t6', 'e16_tol6', 'gmres', '1e-6')]


def load(name):
    return {(r['case'], r['method']): r for r in csv.DictReader(open(ROOT / 'platform' / 'work' / (name + '_eval.csv')))}


def eps(rec):
    return float(rec['eps_2']) if rec and rec['status'] == 'ok' else np.nan


def main():
    t2, t3, t6 = load('e16_tol2'), load('e16_tol3'), load('e16_tol6')
    series = ['series1', 'series2', 'series3']
    # table of values (full precision), saved as npz
    arrays = {}
    for c, _ in CASES:
        arrays[c + '_steps_eps2'] = np.array([eps(t2[(c, s)]) for s in STEPS])
        arrays[c + '_gmres_eps2'] = np.array([eps(t2[(c, 'gmres')]), eps(t3[(c, 'gmres')]), eps(t6[(c, 'gmres')])])
        arrays[c + '_series_eps2'] = np.array([eps(t2[(c, s)]) for s in series])
    arrays['gmres_tolerances'] = np.array([1e-2, 1e-3, 1e-6])
    arrays['steps'] = np.array([0, 1, 2, 3])
    (ROOT / 'analysis' / 'out').mkdir(parents=True, exist_ok=True)
    np.savez(ROOT / 'analysis' / 'out' / 'e16_tolerance.npz', **arrays)

    labels = ['0', '1', '2', '3', 'GMRES\n1e-2', 'GMRES\n1e-3', 'GMRES\n1e-6', 'series\n1', 'series\n2', 'series\n3']
    x = np.arange(len(labels))
    fig, ax = plt.subplots(1, 4, figsize=(18, 4.9), sharey=True)
    for a, (c, name) in zip(ax, CASES):
        steps = arrays[c + '_steps_eps2']
        g = arrays[c + '_gmres_eps2']
        s = arrays[c + '_series_eps2']
        vals = np.concatenate([steps, g, s])
        cols = [PAL['red_2']] * 4 + [PAL['blue_main']] * 3 + [PAL['teal']] * 3
        for i, v in enumerate(vals):
            if np.isfinite(v):
                a.bar(i, v, color=cols[i])
            else:
                a.text(i, 3e-7, 'failed', rotation=90, ha='center', va='bottom', fontsize=8, color=PAL['red_strong'])
        a.axhline(0.1, color='k', ls=':', lw=1)
        a.set_yscale('log'); a.set_ylim(1e-7, 20)
        a.set_xticks(x); a.set_xticklabels(labels, fontsize=7)
        a.set_title(name, fontsize=11)
    ax[0].set_ylabel('error against LU, eps_2 (log)')
    fig.suptitle('E16 study, cap 30. Light red: GMRES stopped after k steps (k = 0 to 3). Blue: GMRES to a tolerance. '
                 'Teal: series expansion. Dotted: 10 % pass limit.', fontsize=10)
    fig.tight_layout()
    out = ROOT / 'figures'
    out.mkdir(exist_ok=True)
    fig.savefig(out / 'e16_tolerance.png', dpi=200)
    fig.savefig(out / 'e16_tolerance.pdf')
    print('wrote figures/e16_tolerance.png and .pdf; analysis/out/e16_tolerance.npz')


if __name__ == '__main__':
    main()
