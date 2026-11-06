"""Table-of-contents graphic (8.25 x 4.45 cm, JCTC): the spanning tree on the BaTiO3 cluster (case H11) and the median
error of the series of order 1 to 3 and of GMRES against dense LU. Errors come from analysis/out/p100_error_summary.csv
(fractions; shown in percent); time and memory ratios from manuscript/numbers_results.tex (CostTimeMed, CostMemMed).
Writes figures/toc_graphic.png/.pdf."""
import csv
import json
import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import cdist
from ase import Atoms
from ase.data.colors import jmol_colors
from ase.visualize.plot import plot_atoms

ROOT = Path(__file__).resolve().parents[1]
PROD = Path('/home/mohammed/Desktop/REAL-WORK/msspec_ML/prod3000/runs')
BLUE, RED, GREY = '#0F4D92', '#B64342', '#CFCECE'
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 7, 'axes.linewidth': 0.8,
                     'axes.spines.top': False, 'axes.spines.right': False})

mat = {r['mat_id']: r for r in csv.DictReader(open(ROOT / 'platform/materials_v2.csv'))}['H11']
d = PROD / mat['run_tag']
cfg = json.loads((d / 'config.json').read_text())
lines = (d / 'cluster.xyz').read_text().splitlines()
n = int(lines[0])
sym = [lines[2 + i].split()[0] for i in range(n)]
pos = np.array([[float(x) for x in lines[2 + i].split()[1:4]] for i in range(n)])
em = int(cfg['absorber_idx'])
T = minimum_spanning_tree(cdist(pos, pos)).toarray()
edges = list(zip(*np.nonzero(T)))
assert len(edges) == n - 1

fig = plt.figure(figsize=(8.25 / 2.54, 4.45 / 2.54))
ax = fig.add_axes([0.0, 0.05, 0.5, 0.88])
atoms = Atoms(sym, positions=pos)
atoms.rotate(-90, 'x', rotate_cell=False)
xy = atoms.positions[:, :2]
cols = [tuple(jmol_colors[z]) for z in atoms.numbers]
cols[em] = (1.0, 0.0, 1.0)
rad = [0.32] * n
rad[em] = 0.5
plot_atoms(atoms, ax, radii=rad, colors=cols, show_unit_cell=0)
centers = np.array([pt.center for pt in ax.patches])
off = centers.mean(axis=0) - xy.mean(axis=0)
for i, j in edges:
    ax.plot(xy[[i, j], 0] + off[0], xy[[i, j], 1] + off[1], color='#222222', lw=0.5, zorder=10, alpha=0.85)
ax.set_aspect('equal'); ax.axis('off')
ax.set_title('spanning-tree preconditioner', fontsize=7, loc='center')

rows = {r['method']: r for r in csv.DictReader(open(ROOT / 'analysis/out/p100_error_summary.csv')) if r['stratum'] == 'all'}
labs = ['series 1', 'series 2', 'series 3', 'GMRES']
vals = [100 * float(rows[k]['e2_median']) for k in ('series1', 'series2', 'series3', 'gmres_e2c67')]
ax = fig.add_axes([0.62, 0.30, 0.36, 0.55])
ax.bar(range(4), vals, color=[GREY, GREY, GREY, RED], width=0.7)
ax.set_yscale('log'); ax.set_xticks(range(4)); ax.set_xticklabels(labs, rotation=45, ha='right', fontsize=6)
ax.set_ylabel('median error (%)', fontsize=6.5); ax.tick_params(labelsize=6)
for x, v in enumerate(vals):
    ax.text(x, v * 1.15, '%.3g' % v, ha='center', va='bottom', fontsize=6)
ax.set_ylim(0.2, max(vals) * 4)
ax.set_title('%s clusters, vs dense LU' % rows['gmres_e2c67']['n'], fontsize=7, loc='center')
tex = (ROOT / 'manuscript/numbers_results.tex').read_text()
tm, mm = (re.search(r'\\newcommand\{\\%s\}\{([^}]*)\}' % k, tex).group(1) for k in ('CostTimeMed', 'CostMemMed'))
fig.text(0.74, 0.03, '%s$\\times$ faster, %s$\\times$ less memory' % (tm, mm), ha='center', fontsize=6, color=BLUE)
fig.savefig(ROOT / 'figures/toc_graphic.png', dpi=300)
fig.savefig(ROOT / 'figures/toc_graphic.pdf')
print('values', dict(zip(labs, vals)), 'time', tm, 'mem', mm)
