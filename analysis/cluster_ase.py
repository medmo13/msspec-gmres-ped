"""Draw a test cluster with ASE (three views), emitter marked.

The cluster is read from the prod3000 run folder of the case (cluster.xyz and config.json, as in
analysis/plot_tree_cluster.py). Atoms are drawn with ase.visualize.plot.plot_atoms in the style of the
v1 msspec_gpu figures (see draw()).
"""
import csv
import json
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.visualize.plot import plot_atoms

ROOT = Path(__file__).resolve().parents[1]
PROD = Path('/home/mohammed/Desktop/REAL-WORK/msspec_ML/prod3000/runs')


def load(case, tag=None):
    # tag: run_tag of a cluster that is not in materials_v2.csv (the large clusters of plot_big.py)
    mat = {r['mat_id']: r for r in csv.DictReader(open(ROOT / 'platform/materials_v2.csv'))}.get(case, {'run_tag': tag})
    d = PROD / mat['run_tag']
    cfg = json.loads((d / 'config.json').read_text())
    lines = (d / 'cluster.xyz').read_text().splitlines()
    n = int(lines[0])
    sym = [lines[2 + i].split()[0] for i in range(n)]
    pos = [[float(x) for x in lines[2 + i].split()[1:4]] for i in range(n)]
    em = int(cfg['absorber_idx'])
    assert sym[em] == cfg['absorber_element'] and ('natoms' not in mat or n == int(mat['natoms']))
    return Atoms(sym, positions=pos), em, mat


def draw(axs, case, tag=None):
    """Top, side and perspective views with ASE's renderer, as draw_cluster() of
    msspec_gpu/archive/v1_msspec_gpu_inversion/harness/plots_physics.py: jmol colours, spheres of 0.55 covalent
    radius, emitter in magenta."""
    from ase.data import covalent_radii
    from ase.data.colors import jmol_colors
    atoms, em, mat = load(case, tag)
    cols = [tuple(jmol_colors[z]) for z in atoms.numbers]
    cols[em] = (1.0, 0.0, 1.0)
    rad = [0.55 * covalent_radii[z] for z in atoms.numbers]
    for ax, rot, t in zip(axs, ('', '-90x', '-65x,-25y'), ('top (along z)', 'side', 'perspective')):
        plot_atoms(atoms, ax, rotation=rot, radii=rad, colors=cols, show_unit_cell=0)
        ax.set_title(t, fontsize=6.5)
        ax.axis('off')
    return mat, atoms.get_chemical_symbols()[em]
