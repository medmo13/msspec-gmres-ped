"""Plot of the E17 energy-scan small test (Ru hcp 0001, 26 atoms, normal emission, 5 energies). Test only, not a result.

Reads runs/e17_scan_test/H42/{lu,gmres,series3}_0 and writes platform/work/e17_scan_test.png.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
B = ROOT / 'runs' / 'e17_scan_test' / 'H42'
plt.rcParams.update({'font.size': 11, 'axes.linewidth': 1.3, 'axes.spines.top': False, 'axes.spines.right': False})
R = {m: np.load(B / (m + '_0') / 'result.npz') for m in ('lu', 'gmres', 'series3')}
E = R['lu']['energy']
iters = np.array([int(x) for x in json.load(open(B / 'gmres_0' / 'row.json'))['iters_list'].split()]).reshape(len(E), -1)
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
ax[0].plot(E, R['lu']['cross_section'], 'o-', color='#B64342', label='dense LU')
ax[0].plot(E, R['gmres']['cross_section'], 's--', color='#0F4D92', label='GMRES, tol $10^{-2}$')
ax[0].plot(E, R['series3']['cross_section'], '^:', color='#42949D', label='series, order 3')
ax[0].set_yscale('log'); ax[0].set_xlabel('Kinetic energy (eV)'); ax[0].set_ylabel('Cross section at $\\theta=0$')
ax[0].set_title('(a) Energy scan', loc='left'); ax[0].legend(frameon=False, fontsize=9)
for m, c, mk in (('gmres', '#0F4D92', 's'), ('series3', '#42949D', '^')):
    err = 100 * np.abs(R[m]['cross_section'] - R['lu']['cross_section']) / np.abs(R['lu']['cross_section'])
    ax[1].plot(E, err, mk + '-', color=c, label={'gmres': 'GMRES', 'series3': 'series, order 3'}[m])
ax[1].axhline(10, color='#555555', ls='--', lw=1); ax[1].set_yscale('log')
ax[1].set_xlabel('Kinetic energy (eV)'); ax[1].set_ylabel('Error at each energy (%)')
ax[1].set_title('(b) Error against dense LU', loc='left'); ax[1].legend(frameon=False, fontsize=9)
ax[2].plot(E, iters.max(axis=1), 'o-', color='#0F4D92', label='largest of 10 right-hand sides')
ax[2].plot(E, iters.min(axis=1), 'o--', color='#9CB8D6', label='smallest')
ax[2].axhline(12, color='#555555', ls='--', lw=1); ax[2].set_ylim(0, 13)
ax[2].set_xlabel('Kinetic energy (eV)'); ax[2].set_ylabel('GMRES steps'); ax[2].set_title('(c) Steps per energy (cap 12)', loc='left')
ax[2].legend(frameon=False, fontsize=9)
fig.tight_layout(); out = ROOT / 'platform' / 'work' / 'e17_scan_test.png'; fig.savefig(out, dpi=150); print('wrote', out)
print('energies', E.tolist(), 'max steps', iters.max(axis=1).tolist())
