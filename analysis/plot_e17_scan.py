"""Energy scans of the strong scatterers done so far (four 10 to 100 eV scans and the 3-energy test of a 277-atom cluster), 10 to 100 eV in 5 eV steps (set e17_scan; DEVIATIONS D37).

Reads runs/e17_scan/<case>/{lu,gmres,series3}_0 and runs/e17_scan/rows.csv; writes analysis/out/e17_scan.csv and
figures/e17_scan.png/.pdf. Columns: (a) cross section at theta 0 against energy, (b) error against dense LU at each
energy, (c) GMRES steps (largest and smallest of the 10 right-hand sides per energy).
"""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
B = ROOT / 'runs' / 'e17_scan'
CASES = [('H42', 'Ru hcp (0001), polar, 26 atoms', B), ('H17', 'PbB rocksalt (011), low-index, 48 atoms', B),
         ('H28', 'TiO$_2$ rutile, high-index, 85 atoms', B), ('H48', 'Hf hcp, semipolar, 106 atoms', B),
         ('B01', 'LiF rocksalt (m201), 277 atoms', ROOT / 'runs' / 'e17_big')]
plt.rcParams.update({'font.size': 10, 'axes.linewidth': 1.2, 'axes.spines.top': False, 'axes.spines.right': False})
fig, axs = plt.subplots(len(CASES), 2, figsize=(9.5, 3.1 * len(CASES)), sharex=True)
out = []
rout = []
smax = 0
for i, (c, title, B) in enumerate(CASES):
    ld = lambda m: np.load(B / c / (m + '_0') / 'result.npz') if (B / c / (m + '_0') / 'result.npz').exists() else None
    lu, g, s = ld('lu'), ld('gmres'), ld('series3')
    gdir = B / c / 'gmres_0'
    if g is None:                       # cap 12 failed: use the cap-67 rerun (PROTOCOL practice, DEVIATIONS D37)
        gdir = ROOT / 'runs' / 'e17_scan_cap67' / c / 'gmres_0'
        g = np.load(gdir / 'result.npz')
    ax = axs[i]
    ref = lu if lu is not None else g
    E = ref['energy']
    if lu is not None:
        ax[0].plot(E, lu['cross_section'], '-', lw=4.5, color='#E8A9A8', label='dense LU')
    if g is not None:
        ax[0].plot(E, g['cross_section'], 's', ms=5, mfc='none', mew=1.3, color='#0F4D92', label='GMRES, tol $10^{-2}$')
    if s is not None:
        ax[0].plot(E, s['cross_section'], '^:', ms=3, color='#42949D', label='series, order 3')
    ax[0].set_yscale('log'); ax[0].set_ylabel('Cross section'); ax[0].set_title(title, loc='left', fontsize=10)
    rtxt = []
    if lu is not None:
        for arr, mk, col, lab in ((g, 's-', '#0F4D92', 'GMRES'), (s, '^-', '#42949D', 'series, order 3')):
            if arr is None:
                continue
            err = 100 * np.abs(arr['cross_section'] - lu['cross_section']) / np.abs(lu['cross_section'])
            out += [dict(case=c, method=lab, energy=e, err_pct=x) for e, x in zip(E, err)]
            rr = np.corrcoef(arr['cross_section'], lu['cross_section'])[0, 1]     # Pearson r over the energy scan
            rout.append(dict(case=c, method=lab, one_minus_r=1 - rr))
            rtxt.append('%s %.1e' % ('GMRES' if lab == 'GMRES' else 'series 3', 1 - rr))
        ax[0].text(0.98, 0.04, '$1-r$: ' + ', '.join(rtxt), transform=ax[0].transAxes, ha='right', va='bottom', fontsize=8)
    if g is not None:
        row = json.load(open(gdir / 'row.json'))
        it = np.array([int(x) for x in row['iters_list'].split()]).reshape(len(E), -1)
        ax[1].plot(E, it.max(axis=1), 'o-', ms=3, color='#0F4D92', label='largest of 10 right-hand sides')
        ax[1].plot(E, it.min(axis=1), 'o--', ms=3, color='#9CB8D6', label='smallest')
        smax = max(smax, it.max())
        out += [dict(case=c, method='gmres_steps_max', energy=e, err_pct=x) for e, x in zip(E, it.max(axis=1))]
    ax[1].set_ylabel('GMRES steps')
    if i == 0:
        for a in ax:
            a.legend(frameon=False, fontsize=8)
for a in axs[:, 1]:
    a.set_ylim(0, smax + 1)   # one common scale for the step counts
for a in axs[-1]:
    a.set_xlabel('Kinetic energy (eV)')
fig.tight_layout()
fig.savefig(ROOT / 'figures' / 'e17_scan.png', dpi=200); fig.savefig(ROOT / 'figures' / 'e17_scan.pdf')
with open(ROOT / 'analysis' / 'out' / 'e17_scan.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['case', 'method', 'energy', 'err_pct']); w.writeheader(); w.writerows(out)
with open(ROOT / 'analysis' / 'out' / 'e17_scan_r.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['case', 'method', 'one_minus_r']); w.writeheader(); w.writerows(rout)
print('wrote figures/e17_scan.png, analysis/out/e17_scan.csv', len(out), 'rows')
