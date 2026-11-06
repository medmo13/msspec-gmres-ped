"""Fig. 1 of the paper: spanning tree on one test cluster and the block pattern of M.

Case H11 (barium titanate perovskite, 1s, 118 atoms). The cluster is read from the prod3000 run folder
(cluster.xyz and config.json, the same files platform/rho_probe.py read_cluster uses, parsed here
without importing MsSpec). The tree is the minimum-distance spanning tree (Prim, start at atom 0, as in
MS_TREE_SERIES_BATCH). It is built with scipy and checked against a plain Prim loop; edge count and total
length are printed. Writes figures/tree_cluster.png/.pdf and analysis/out/tree_cluster_numbers.csv.
"""
import csv
import json
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
CASE = 'H11'
LMAX = 8
BLUE, RED, TEAL, GREY = '#0F4D92', '#B64342', '#42949E', '#CFCECE'
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 10, 'axes.linewidth': 1.2,
                     'axes.spines.top': False, 'axes.spines.right': False})

mat = {r['mat_id']: r for r in csv.DictReader(open(ROOT / 'platform/materials_v2.csv'))}[CASE]
d = PROD / mat['run_tag']
cfg = json.loads((d / 'config.json').read_text())
lines = (d / 'cluster.xyz').read_text().splitlines()
n = int(lines[0])
sym = [lines[2 + i].split()[0] for i in range(n)]
pos = np.array([[float(x) for x in lines[2 + i].split()[1:4]] for i in range(n)])
em = int(cfg['absorber_idx'])
assert n == int(mat['natoms']) and sym[em] == cfg['absorber_element']

D = cdist(pos, pos)
T = minimum_spanning_tree(D).toarray()
edges_sp = sorted(tuple(sorted(e)) for e in zip(*np.nonzero(T)))


def prim(D):
    n = len(D)
    seen = [0]
    dist = D[0].copy()
    par = np.zeros(n, int)
    order, edges = [0], []
    done = np.zeros(n, bool)
    done[0] = True
    for _ in range(n - 1):
        dd = np.where(done, np.inf, dist)
        j = int(np.argmin(dd))
        edges.append(tuple(sorted((int(par[j]), j))))
        done[j] = True
        order.append(j)
        upd = D[j] < dist
        dist = np.where(upd, D[j], dist)
        par = np.where(upd, j, par)
    return order, sorted(edges)


order, edges_pr = prim(D)
len_sp = sum(D[i, j] for i, j in edges_sp)
len_pr = sum(D[i, j] for i, j in edges_pr)
print('atoms', n, 'edges scipy', len(edges_sp), 'edges Prim', len(edges_pr),
      'total length scipy %.4f Prim %.4f A' % (len_sp, len_pr))
assert len(edges_sp) == n - 1 == len(edges_pr) and abs(len_sp - len_pr) < 1e-6

el = sorted(set(sym))
fig = plt.figure(figsize=(11, 4.4))
ax = fig.add_axes([0.03, 0.17, 0.36, 0.73])
# (a) ASE drawing (jmol colours, emitter magenta) of the side view; the cluster is turned by -90 degrees about x so that z is up
atoms = Atoms(sym, positions=pos)
atoms.rotate(-90, 'x', rotate_cell=False)
xy = atoms.positions[:, :2]
cols = [tuple(jmol_colors[z]) for z in atoms.numbers]
cols[em] = (1.0, 0.0, 1.0)
rad = [0.32] * n
rad[em] = 0.5          # equal small spheres so that the tree edges stay visible; the emitter is larger
plot_atoms(atoms, ax, radii=rad, colors=cols, show_unit_cell=0)
centers = np.array([pt.center for pt in ax.patches])
off = centers.mean(axis=0) - xy.mean(axis=0)          # plot_atoms shifts the picture by a constant; same shift for the edges
assert np.allclose(np.sort(centers[:, 0]), np.sort(xy[:, 0] + off[0]), atol=1e-6)
for i, j in edges_pr:
    ax.plot(xy[[i, j], 0] + off[0], xy[[i, j], 1] + off[1], color='#222222', lw=0.8, zorder=10, alpha=0.85, solid_capstyle='round')
ax.set_aspect('equal'); ax.axis('off')
for e, z in sorted({(a, Atoms(a).numbers[0]) for a in el}):
    ax.plot([], [], 'o', color=jmol_colors[z], ms=7, mec='#555555', mew=0.4, label=e)
ax.plot([], [], 'o', color='magenta', ms=7, mec='#555555', mew=0.4, label='emitter')
ax.plot([], [], '-', color='#333333', lw=1.0, label='tree edge')
ax.legend(frameon=False, fontsize=8, ncol=3, loc='upper center', bbox_to_anchor=(0.5, 0.02), columnspacing=0.9, handletextpad=0.3)
ax.set_title('(a) Tree edges, side view', loc='left', fontsize=10)

rank = np.empty(n, int); rank[order] = np.arange(n)
ax = fig.add_axes([0.45, 0.20, 0.25, 0.70])
full = np.ones((n, n))
ax.imshow(full, cmap='Greys', vmin=0, vmax=3, origin='upper', extent=(0, n, n, 0))
for i, j in edges_pr:
    for a, b in ((i, j), (j, i)):
        ax.add_patch(plt.Rectangle((rank[b], rank[a]), 1, 1, color=RED, lw=0))
for k in range(n):
    ax.add_patch(plt.Rectangle((k, k), 1, 1, color=BLUE, lw=0))
ax.set_xlim(0, n); ax.set_ylim(n, 0); ax.set_aspect('equal')
ax.set_xlabel('atom j (visit order)'); ax.set_ylabel('atom i (visit order)')
ax.set_title('(b) Atom blocks of M', loc='left', fontsize=10)
ax.plot([], [], 's', color=BLUE, label='diagonal'); ax.plot([], [], 's', color=RED, label='tree edge')
ax.plot([], [], 's', color=GREY, label='dropped')
ax.legend(frameon=False, fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=3)

ax = fig.add_axes([0.78, 0.20, 0.20, 0.70])
q = (LMAX + 1) ** 2
mem_full = 16.0 * (n * q) ** 2 / 2 ** 20
mem_m = 16.0 * q * q * (n + 2 * (n - 1)) / 2 ** 20
ax.bar([0, 1], [mem_full, mem_m], color=[GREY, RED], width=0.6)
ax.set_xticks([0, 1]); ax.set_xticklabels(['full $A$', 'blocks of $M$'])
ax.set_ylabel('Largest storage (MiB)')
for x, v in ((0, mem_full), (1, mem_m)):
    ax.text(x, v, '%.0f' % v, ha='center', va='bottom', fontsize=9)
ax.set_title('(c) Storage at lmax 8', loc='left', fontsize=10)
fig.savefig(ROOT / 'figures/tree_cluster.png', dpi=200)
fig.savefig(ROOT / 'figures/tree_cluster.pdf')
nn = [np.sort(D[i])[1] for i in range(n)]
with open(ROOT / 'analysis/out/tree_cluster_numbers.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['case', 'natoms', 'lmax', 'block_size', 'edges', 'tree_length_A', 'mem_full_mib', 'mem_blocks_mib',
                'mem_ratio', 'edge_len_max_A', 'edge_len_median_A'])
    el_ = [D[i, j] for i, j in edges_pr]
    w.writerow([CASE, n, LMAX, q, len(edges_pr), '%.4f' % len_pr, '%.2f' % mem_full, '%.2f' % mem_m,
                '%.2f' % (mem_full / mem_m), '%.3f' % max(el_), '%.3f' % np.median(el_)])
print('mem full %.1f MiB, M blocks %.1f MiB' % (mem_full, mem_m))
