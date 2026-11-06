"""E16 small test (A3.6): about 100 directions, one small Cu cluster.

Runs LU, the preconditioner alone (k=0), GMRES after k=1..4 steps, and the
converged GMRES. Each run is one fresh process. Writes platform/work/e16_small/.
Not a reported result and not a timing benchmark.
"""
import csv, json, os, re, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT.parent / 'platform' / 'work' / 'e16_small'
SRC = ROOT / 'msspec_gmres' / 'src'
NTHETA = 101

WORKER = r'''
import json, sys
import numpy as np
cfg = json.loads(sys.argv[1])
sys.path.insert(0, cfg['src'])
from ase.build import bulk
from msspec.calculator import MSSPEC
from msspec.utils import center_cluster, cut_sphere, cut_plane, get_atom_index
a0 = 3.61
cluster = bulk('Cu', a=a0, cubic=True).repeat((2, 2, 2))
center_cluster(cluster)
cluster = cut_sphere(cluster, radius=a0 + 0.01)
cluster = cut_plane(cluster, z=0.1)
cluster.absorber = get_atom_index(cluster, 0, 0, -a0)
print('NATOMS', len(cluster))
if cfg['algorithm'] == 'gmres':
    calc = MSSPEC(spectroscopy='PED', algorithm='gmres', folder=cfg['folder'], txt=None)
    calc.calculation_parameters.set_parameter('linear_solver_tolerance', 1e-2)
    calc.calculation_parameters.set_parameter('linear_solver_max_iterations', 12)
else:
    calc = MSSPEC(spectroscopy='PED', algorithm='inversion', folder=cfg['folder'], txt=None)
calc.set_atoms(cluster)
data = calc.get_theta_scan(level='2p3/2', theta=np.linspace(-70, 70, cfg['nt']),
                           phi=0, kinetic_energy=300.0)
np.save(cfg['out'], np.asarray(data[-1].cross_section))
'''

CASES = [('lu', 'inversion', None)] + [('k%d' % k, 'gmres', k) for k in range(5)] + \
        [('conv', 'gmres', None)]


def run(name, alg, k):
    out = OUT / (name + '.npy')
    cfg = {'src': str(SRC), 'algorithm': alg, 'folder': str(OUT / ('calc_' + name)),
           'out': str(out), 'nt': NTHETA}
    env = dict(os.environ, TMPDIR=str(OUT), OMP_NUM_THREADS='8',
               MSSPEC_TREE_ACTION_THREADS='8', OPENBLAS_NUM_THREADS='1' if alg == 'gmres' else '8')
    for v in ('MSSPEC_TREE_STOP_AFTER', 'MSSPEC_TREE_ITER_DUMP', 'MSSPEC_FALLBACK_DENSE'):
        env.pop(v, None)
    if k is not None:
        env['MSSPEC_TREE_STOP_AFTER'] = str(k)
    t0 = time.time()
    p = subprocess.run([sys.executable, '-c', WORKER, json.dumps(cfg)], env=env,
                       capture_output=True, text=True)
    wall = time.time() - t0
    txt = p.stdout + p.stderr
    (OUT / (name + '.log')).write_text(txt)
    its = [int(m.group(1)) for m in re.finditer(r'TREE_GMRES iterations,rhs=\s+(\d+)', txt)]
    return p.returncode, wall, its, ('stop_after=' in txt), out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print('work dir', OUT, 'load', os.getloadavg()[0], flush=True)
    res = {}
    rows = []
    for name, alg, k in CASES:
        rc, wall, its, announced, out = run(name, alg, k)
        res[name] = np.load(out) if rc == 0 and out.exists() else None
        rows.append(dict(case=name, k='' if k is None else k, rc=rc, wall_s=round(wall, 1),
                         n_iter_lines=len(its), iters_min=min(its) if its else '',
                         iters_max=max(its) if its else '', stop_announced=announced))
        print(rows[-1], flush=True)
    ref = res['lu']
    for r in rows:
        a = res[r['case']]
        if a is None or ref is None:
            r.update(eps2_pct='', epsinf_pct='')
            continue
        r['eps2_pct'] = 100 * np.linalg.norm(a - ref) / np.linalg.norm(ref)
        fl = 0.01 * np.max(np.abs(ref))
        r['epsinf_pct'] = 100 * np.max(np.abs(a - ref) / np.maximum(np.abs(ref), fl))
    with open(OUT / 'e16_small.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) + ['eps2_pct', 'epsinf_pct'])
        w.writeheader(); w.writerows(rows)
    print('done')


if __name__ == '__main__':
    main()
