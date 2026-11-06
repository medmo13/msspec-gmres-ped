"""Why 13 of 97 P100 LU scans differ from the prod3000 scan_inv.npz (relative L2 > 1e-4).

Reads analysis/out/p100_lu_check.csv, runs/p100/<case>/lu_0/result.npz, the prod3000 scan_inv.npz and
the sensitivity test platform/work/perturb/<case>/{base,pert}.npz (all atoms shifted by +-1e-6 A,
made by platform/work/perturb_test.py). Writes analysis/out/lu_prod3000_causes.csv and
figures/lu_vs_prod3000.png/.pdf.
"""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
PROD = Path('/home/mohammed/Desktop/REAL-WORK/msspec_ML/prod3000/runs')
mat = {r['mat_id']: r for r in csv.DictReader(open(ROOT / 'platform/materials_v2.csv'))}
chk = [r for r in csv.DictReader(open(ROOT / 'analysis/out/p100_lu_check.csv')) if r['rel_l2']]
BLUE, RED, TEAL, GREY = '#0F4D92', '#B64342', '#42949E', '#CFCECE'
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 10, 'axes.linewidth': 1.2,
                     'axes.spines.top': False, 'axes.spines.right': False})


def rel(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


rows = []
for r in chk:
    c = r['case']
    p = np.load(PROD / mat[c]['run_tag'] / 'scan_inv.npz')
    o = np.load(ROOT / 'runs/p100' / c / 'lu_0/result.npz')
    rows.append(dict(case=c, name=mat[c]['name'], natoms=mat[c]['natoms'], flag=r['flag'] == 'True',
                     rel_cs=rel(o['cross_section'], p['cross_section']),
                     rel_direct=rel(o['direct_signal'], p['direct_signal'])))
pert = {}
for c in ('H17', 'L12', 'L34', 'H42', 'L43', 'L45'):
    d = ROOT / 'platform/work/perturb' / c
    if (d / 'pert.npz').exists():
        b, q = np.load(d / 'base.npz'), np.load(d / 'pert.npz')
        pert[c] = rel(q['cross_section'], b['cross_section'])
for x in rows:
    x['rel_cs_after_1e-6A_shift'] = pert.get(x['case'], '')
with open(ROOT / 'analysis/out/lu_prod3000_causes.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

fig = plt.figure(figsize=(11, 8.2))
gs = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.28)
# A: two example scans
axA = fig.add_subplot(gs[0, 0])
for c, off, col in (('H17', 0, RED), ('L12', 1.4, BLUE)):
    p = np.load(PROD / mat[c]['run_tag'] / 'scan_inv.npz'); o = np.load(ROOT / 'runs/p100' / c / 'lu_0/result.npz')
    s = 1 / p['cross_section'].max()
    axA.plot(p['theta'], p['cross_section'] * s + off, color='k', lw=1.6, label='saved earlier scan' if off == 0 else None)
    axA.plot(o['theta'], o['cross_section'] * s + off, color=col, lw=1.2, ls='--', label='our LU, same cluster' if off == 0 else None)
    axA.text(-88, off + 1.02, '%s, %s atoms: L2 diff %.1f %%' % (mat[c]['name'], mat[c]['natoms'],
             100 * [x['rel_cs'] for x in rows if x['case'] == c][0]), fontsize=8.5)
axA.set_xlabel('theta (deg)'); axA.set_ylabel('intensity / max (curves shifted)'); axA.set_ylim(0, 2.9)
axA.legend(frameon=False, fontsize=8.5, loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2); axA.set_title('A  The curves look alike; the sizes differ', loc='left', fontsize=10)
# B: all 97, ranked
axB = fig.add_subplot(gs[0, 1])
allv = sorted([(float(r['rel_l2']), r['case'], r['flag'] == 'True') for r in chk])
xs = np.arange(len(allv))
axB.bar(xs, [v[0] for v in allv], color=[RED if v[2] else GREY for v in allv], width=0.9)
axB.set_yscale('log'); axB.axhline(1e-4, color='k', lw=1, ls=':'); axB.text(1, 1.3e-4, 'flag level 1e-4', fontsize=8.5)
axB.set_xlabel('97 cases, sorted'); axB.set_ylabel('L2 difference, our LU vs earlier scan')
axB.set_title('B  84 cases agree to ~1e-6; 13 do not', loc='left', fontsize=10)
# C: direct signal vs cross section
axC = fig.add_subplot(gs[1, 0])
fl = [x for x in rows if x['flag']]
axC.scatter([max(x['rel_direct'], 1e-8) for x in fl], [x['rel_cs'] for x in fl], color=RED, s=36, zorder=3)
for x in fl:
    pass
axC.plot([1e-8, 1e-1], [1e-8, 1e-1], color='k', lw=0.8, ls=':')
axC.set_xscale('log'); axC.set_yscale('log'); axC.set_xlim(1e-8, 1e-1); axC.set_ylim(5e-5, 1e-1)
axC.set_xlabel('difference of the direct signal (no multiple scattering)'); axC.set_ylabel('difference of the full scan')
axC.text(1.3e-8, 6e-5, '"1e-8" = identical', fontsize=8)
axC.set_title('C  The flagged cases already differ before any solve', loc='left', fontsize=10)
# D: sensitivity to a 1e-6 A shift
axD = fig.add_subplot(gs[1, 1])
names = list(pert); vals = [pert[c] for c in names]
cols = [RED if any(x['case'] == c and x['flag'] for x in rows) else TEAL for c in names]
axD.bar(range(len(names)), vals, color=cols)
axD.set_yscale('log'); axD.set_xticks(range(len(names))); axD.set_xticklabels([mat[c]['name'] for c in names], fontsize=8)
axD.set_ylabel('change of our own LU scan'); axD.set_xlabel('after moving every atom by +-0.000001 A')
axD.set_title('D  A rounding-size shift gives the same size of change', loc='left', fontsize=10)
fig.savefig(ROOT / 'figures/lu_vs_prod3000.png', dpi=200, bbox_inches='tight')
fig.savefig(ROOT / 'figures/lu_vs_prod3000.pdf', bbox_inches='tight')
print('flagged', len(fl), 'of', len(chk), '; perturbation', {k: '%.1e' % v for k, v in pert.items()})
