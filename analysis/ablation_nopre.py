"""Ablation: GMRES with the tree preconditioner against GMRES with M = I (D40), same restart 6, cap 67.

Tree rows: the cap-100 rerun (runs/p100_tol6_cap100) over the cap-67 rerun (runs/p100_cap67, runs/p100_tol6_cap67) over the first run (runs/p100, runs/p100_tol6).
No-tree rows: runs/nopre_tol2, runs/nopre_tol6. Reference: dense LU of runs/p100 (97 cases). Steps are deterministic; wall times
come from different sessions and are not a paired timing. Writes analysis/out/ablation_nopre.csv and figures/ablation_nopre.png/.pdf.
"""
import csv
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
R = ROOT / 'runs'


def last_rows(s, m):
    d = {}
    p = R / s / 'rows.csv'
    if p.exists():
        for x in csv.DictReader(open(p)):
            if x['method'] == m:
                d[x['case']] = (s, x)
    return d


def tree_rows(base, rerun, rerun2=None):
    d = last_rows(base, 'gmres')
    d.update(last_rows(rerun, 'gmres'))
    if rerun2:
        d.update(last_rows(rerun2, 'gmres'))     # cap-100 rerun of the 1e-6 cases that did not converge at 67 (D41)
    return d


rho = {x['mat_id']: float(x['rho']) for x in csv.DictReader(open(ROOT / 'platform/materials_v2.csv'))} if 'mat_id' in next(csv.DictReader(open(ROOT / 'platform/materials_v2.csv'))) else {}
lu = {}
for c, (s, x) in last_rows('p100', 'lu').items():
    if x['run_status'] == 'ok':
        lu[c] = np.load(R / s / c / 'lu_0' / 'result.npz')['cross_section']
for c, (s, x) in last_rows('p100_lu_excluded', 'lu').items():      # LU run later (D39)
    lu[c] = np.load(R / s / c / 'lu_0' / 'result.npz')['cross_section']
rows = []
for tol, tb, tr, tr2, ns, *ns2 in [('1e-2', 'p100', 'p100_cap67', None, 'nopre_tol2'),
                                 ('1e-6', 'p100_tol6', 'p100_tol6_cap67', 'p100_tol6_cap100', 'nopre_tol6', 'nopre_tol6_cap100', 'nopre_tol6_cap200')]:
    T, N = tree_rows(tb, tr, tr2), last_rows(ns, 'nopre')
    for n2 in ns2:
        if n2:
            N.update(last_rows(n2, 'nopre'))
    for c in sorted(N):
        s_n, xn = N[c]
        s_t, xt = T[c]
        rec = dict(tol=tol, case=c, natoms=xn['natoms'], rho=xn['rho'], ke=xn['ke'])
        for lab, s, x, m in [('tree', s_t, xt, 'gmres'), ('nopre', s_n, xn, 'nopre')]:
            ok = x['run_status'] == 'ok'
            rec[lab + '_ok'] = int(ok)
            rec[lab + '_steps'] = int(x['iters_max']) if ok and x['iters_max'] else ''
            rec[lab + '_wall_s'] = x['wall_s']
            rec[lab + '_load'] = x['load_before']
            if ok and c in lu:
                met = metrics(np.load(R / s / c / (m + '_0') / 'result.npz')['cross_section'], lu[c])
                rec[lab + '_eps2'], rec[lab + '_epsinf'] = met['eps_2'], met['eps_inf']
            else:
                rec[lab + '_eps2'] = rec[lab + '_epsinf'] = ''
        rows.append(rec)
out = ROOT / 'analysis/out/ablation_nopre.csv'
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print('columns:', list(rows[0]))

fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.8))
col = {'1e-2': '#0F4D92', '1e-6': '#B64342'}
for tol in ('1e-2', '1e-6'):
    R_ = [r for r in rows if r['tol'] == tol and r['tree_ok'] and r['nopre_ok']]
    st, sn = np.array([r['tree_steps'] for r in R_]), np.array([r['nopre_steps'] for r in R_])
    ax[0].scatter(st, sn, s=14, color=col[tol], alpha=0.8, label='tolerance %s' % tol)
    rr = np.array([float(r['rho']) for r in R_])
    ax[1].scatter(rr, sn / st, s=14, color=col[tol], alpha=0.8)
    e = [r for r in R_ if r['tree_eps2'] != '' and r['nopre_eps2'] != '']
    ax[2].scatter([float(r['tree_eps2']) * 100 for r in e], [float(r['nopre_eps2']) * 100 for r in e], s=14, color=col[tol], alpha=0.8)
    f = [r for r in rows if r['tol'] == tol]
    print('tol', tol, 'cases', len(f), 'both ok', len(R_), 'tree fail', sum(not r['tree_ok'] for r in f), 'nopre fail', sum(not r['nopre_ok'] for r in f),
          '| steps ratio nopre/tree median %.2f min %.2f max %.2f' % (np.median(sn / st), (sn / st).min(), (sn / st).max()),
          '| median steps tree %.0f nopre %.0f' % (np.median(st), np.median(sn)),
          '| median eps2 tree %.3g nopre %.3g' % (np.median([float(r['tree_eps2']) for r in e]) * 100, np.median([float(r['nopre_eps2']) for r in e]) * 100),
          '| nopre steps<=tree in %d' % int(np.sum(sn <= st)))
    for r in f:
        if not r['tree_ok'] or not r['nopre_ok']:
            print('   not ok:', r['case'], r['natoms'], 'rho %.2f' % float(r['rho']), 'tree_ok', r['tree_ok'], 'nopre_ok', r['nopre_ok'])
for a_ in ax[:1]:
    lim = [1, 70]; a_.plot(lim, lim, 'k:', lw=1); a_.set_xscale('log'); a_.set_yscale('log')
    a_.set_xlabel('steps with tree'); a_.set_ylabel('steps without tree (M = I)'); a_.legend(frameon=False, fontsize=8)
ax[1].axhline(1, color='k', ls=':', lw=1); ax[1].set_xscale('log'); ax[1].set_xticks([0.4, 0.7, 1, 2, 3]); ax[1].set_xticklabels(['0.4', '0.7', '1', '2', '3']); ax[1].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax[1].set_xlabel(r'$\rho$'); ax[1].set_ylabel('steps without tree / steps with tree')
lim = [1e-6, 30]; ax[2].plot(lim, lim, 'k:', lw=1); ax[2].set_xscale('log'); ax[2].set_yscale('log')
ax[2].set_xlabel(r'$\varepsilon_2$ with tree (%)'); ax[2].set_ylabel(r'$\varepsilon_2$ without tree (%)')
for a_, t in zip(ax, ['(a) steps', '(b) step ratio against rho', '(c) error against dense LU']):
    a_.set_title(t, loc='left', fontsize=10); a_.spines[['top', 'right']].set_visible(False)
fig.tight_layout(); fig.savefig(ROOT / 'figures/ablation_nopre.png', dpi=200); fig.savefig(ROOT / 'figures/ablation_nopre.pdf')
