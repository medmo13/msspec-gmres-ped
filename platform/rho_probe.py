"""Spectral radius rho(G0T) of one prod3000 run, via MsSpec EIG with algorithm='power'.

Prints one JSON line: run_tag, ke, lmax, rho, wall_s.
Usage: python platform/rho_probe.py <run_dir> <ke_eV> <lmax>
Scratch output goes to platform/work/<run_tag>/ (never into the read-only prod3000 tree).
"""
import sys, time, json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'code' / 'msspec_gmres' / 'src'))
from ase import Atoms
from msspec.calculator import MSSPEC


def read_cluster(run_dir):
    """Read a prod3000 cluster.xyz (count line, comment line, then species x y z ...).

    The absorber comes from config.json ('absorber_idx'), not from a flag column:
    the extra columns in cluster.xyz differ between runs. The atom at that index is
    checked against config.json's absorber_element and absorber_xyz.
    """
    cfg = json.loads((Path(run_dir) / 'config.json').read_text())
    lines = (Path(run_dir) / 'cluster.xyz').read_text().splitlines()
    n = int(lines[0])
    symbols, pos = [], []
    for i in range(n):
        f = lines[2 + i].split()
        symbols.append(f[0])
        pos.append([float(x) for x in f[1:4]])
    idx = int(cfg['absorber_idx'])
    if symbols[idx] != cfg['absorber_element'] or \
            np.linalg.norm(np.array(pos[idx]) - np.array(cfg['absorber_xyz'])) > 1e-3:
        raise ValueError('absorber mismatch in %s' % run_dir)
    atoms = Atoms(symbols=symbols, positions=pos)
    atoms.absorber = idx
    return atoms


def main(run_dir, ke, lmax):
    run_dir = Path(run_dir)
    run_tag = run_dir.name
    atoms = read_cluster(run_dir)
    folder = ROOT / 'platform' / 'work' / run_tag
    calc = MSSPEC(spectroscopy='EIG', algorithm='power', folder=str(folder), txt=None)
    calc.tmatrix_parameters.set_parameter('lmax_mode', 'imposed')
    calc.tmatrix_parameters.set_parameter('lmaxt', int(lmax))
    calc.set_atoms(atoms)
    t0 = time.time()
    calc.get_eigen_values(kinetic_energy=float(ke))
    wall = time.time() - t0
    res = Path(calc.tmp_folder) / 'output' / 'results.dat'
    rho = float(res.read_text().split()[1])
    return {'run_tag': run_tag, 'ke': float(ke), 'lmax': int(lmax),
            'natoms': len(atoms), 'rho': rho, 'wall_s': round(wall, 3)}


if __name__ == '__main__':
    print(json.dumps(main(sys.argv[1], sys.argv[2], sys.argv[3])))
