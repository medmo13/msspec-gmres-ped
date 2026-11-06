"""E3 figure from analysis/out/<tag>_curves.csv (made by e3_convergence.py).

A: relative preconditioned residual against step (all cases, median line), B: cross-section error
eps_2 against step, C: eps_2 against residual with the y = x line,
D: ratio eps_2 / residual per case (spread of the constant in Appendix B).
Also writes analysis/out/<tag>_summary.csv with the numbers the paper quotes.
"""
import csv, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
tag = sys.argv[1] if len(sys.argv) > 1 else 'p100_e3'
rows = list(csv.DictReader(open(ROOT / 'analysis/out' / (tag + '_curves.csv'))))
mat = {r['mat_id']: r for r in csv.DictReader(open(ROOT / 'platform/materials_v2.csv'))}
BLUE, RED, TEAL, GREY = '#0F4D92', '#B64342', '#42949E', '#CFCECE'
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 10, 'axes.linewidth': 1.2,
                     'axes.spines.top': False, 'axes.spines.right': False})
cases = sorted({r['case'] for r in rows})
strong = {c: float(mat[c]['rho']) >= 1 for c in cases}
steps = list(range(7))


def get(c, key):
    d = {r['step']: float(r[key]) for r in rows if r['case'] == c}
    return d


fig = plt.figure(figsize=(14, 8.2)); gs = fig.add_gridspec(2, 6, hspace=0.38, wspace=1.2)
ax = np.array([[fig.add_subplot(gs[0, 0:3]), fig.add_subplot(gs[0, 3:6])], [fig.add_subplot(gs[1, 0:2]), fig.add_subplot(gs[1, 2:4])]])
axE = fig.add_subplot(gs[1, 4:6])
for (a, key, lab, ttl) in ((ax[0, 0], 'resid_max', 'relative preconditioned residual', 'A  Residual falls with each step'),
                           (ax[0, 1], 'eps_2', 'cross-section error eps_2', 'B  The error falls with it')):
    for c in cases:
        d = get(c, key)
        y = [d[str(k)] for k in steps if str(k) in d]
        a.plot(range(len(y)), y, color=RED if strong[c] else BLUE, alpha=0.25, lw=0.8)
    for sel, col, nm in ((True, RED, 'rho >= 1 (strong)'), (False, BLUE, 'rho < 1 (weak)')):
        med = [np.median([get(c, key)[str(k)] for c in cases if strong[c] == sel and str(k) in get(c, key)]) for k in steps]
        a.plot(steps, med, color=col, lw=2.2, marker='o', label='median, ' + nm)
    a.set_yscale('log'); a.set_xlabel('GMRES steps (0 = preconditioner alone)'); a.set_ylabel(lab); a.set_title(ttl, loc='left', fontsize=10)
axB = ax[0, 1]; axB.legend(frameon=False, fontsize=8.5)
axC = ax[1, 0]
xs, ys = [], []
for c in cases:
    for r in rows:
        if r['case'] == c and r['step'] != '0':
            x, y = float(r['resid_max']), float(r['eps_2'])
            axC.scatter(x, y, s=10, color=RED if strong[c] else BLUE, alpha=0.5, lw=0)
            xs.append(x); ys.append(y)
axC.plot([1e-4, 3], [1e-4, 3], color='k', lw=0.8, ls=':')
axC.set_xscale('log'); axC.set_yscale('log'); axC.set_xlabel('relative preconditioned residual'); axC.set_ylabel('eps_2')
axC.set_title('C  Error against residual, all steps', loc='left', fontsize=10)
axD = ax[1, 1]
fin = [(c, float(r['eps_2']) / float(r['resid_max'])) for c in cases for r in rows if r['case'] == c and r['step'] == 'final' and float(r['resid_max']) > 0]
q = sorted(v for _, v in fin)
axD.bar(range(len(q)), q, color=[GREY] * len(q)); axD.set_yscale('log'); axD.axhline(1, color='k', lw=1, ls=':')
axD.set_xlabel('cases, sorted'); axD.set_ylabel('eps_2 / residual at the stop'); axD.set_title('D  Error per unit residual at the stopping point', loc='left', fontsize=10)
# E: Pearson 1 - r against dense LU after each step (from the saved scans of runs/<tag>/<case>/k<j>_0)
sdir = ROOT / 'runs' / tag
omr = {}
for c in cases:
    ref = np.load(sdir / c / 'lu_0' / 'result.npz')['cross_section']
    omr[c] = {k: 1 - np.corrcoef(np.load(sdir / c / ('k%d_0' % k) / 'result.npz')['cross_section'], ref)[0, 1]
              for k in steps if (sdir / c / ('k%d_0' % k) / 'result.npz').exists()}
    axE.plot(list(omr[c]), list(omr[c].values()), color=RED if strong[c] else BLUE, alpha=0.25, lw=0.8)
for sel, col in ((True, RED), (False, BLUE)):
    axE.plot(steps, [np.median([omr[c][k] for c in cases if strong[c] == sel and k in omr[c]]) for k in steps], color=col, lw=2.2, marker='o')
axE.set_yscale('log'); axE.set_xlabel('GMRES steps (0 = preconditioner alone)'); axE.set_ylabel('1 - r (Pearson)')
axE.set_title('E  Pearson correlation falls with it', loc='left', fontsize=10)
with open(ROOT / 'analysis/out' / (tag + '_pearson.csv'), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['case', 'step', 'one_minus_r'])
    [w.writerow([c, k, v]) for c in cases for k, v in omr[c].items()]
fig.savefig(ROOT / 'figures' / (tag + '.png'), dpi=200, bbox_inches='tight'); fig.savefig(ROOT / 'figures' / (tag + '.pdf'), bbox_inches='tight')
xs, ys = np.array(xs), np.array(ys)
from scipy.stats import spearmanr
summ = {'n_cases': len(cases), 'spearman_resid_eps2': float(spearmanr(xs, ys)[0]),
        'median_ratio_final': float(np.median(q)) if q else float('nan'), 'max_ratio_final': float(max(q)) if q else float('nan'),
        'min_ratio_final': float(min(q)) if q else float('nan')}
for k in (3, 6):
    summ['median_resid_step%d' % k] = float(np.median([get(c, 'resid_max')[str(k)] for c in cases]))
    summ['median_eps2_step%d' % k] = float(np.median([get(c, 'eps_2')[str(k)] for c in cases]))
with open(ROOT / 'analysis/out' / (tag + '_summary.csv'), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['key', 'value']); [w.writerow([k, v]) for k, v in summ.items()]
print(summ)
