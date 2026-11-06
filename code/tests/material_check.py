"""Material check: GMRES build against the original build.

For each material, runs in fresh subprocesses:
  orig direct     original build, algorithm='inversion'        (dense LU)
  new  direct     code/msspec_gmres build, algorithm='inversion'
  orig tree       original build, linear_solver='tree_gmres'
  new  tree       code/msspec_gmres build, algorithm='gmres'

Checks:
  direct: new == orig, bit for bit (gate)
  tree:   new and orig each within 1e-4 relative L2 of the LU reference (gate)
  tree:   new == orig to 1e-12 is reported but is NOT a gate. The C1 in-cycle
          early stop changes the iteration path on purpose (the original build
          runs all MAXIT iterations; the new build stops once the estimate
          meets TOL), so the two tree results differ at the solver-tolerance
          level by design.

Records wall time and peak RSS per run. Writes a JSON report next to this file.

usage: MSSPEC_ORIG_SRC=<pre-GMRES msspec src folder> python tests/material_check.py [Ag Ni Cu]
"""
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CODE = HERE.parent
NEW_SRC = CODE / 'msspec_gmres' / 'src'
# The pre-GMRES build is not in this workspace (no git history, no clean copy),
# so its path comes from the environment. main() stops if it is missing.
ORIG_SRC = Path(os.environ['MSSPEC_ORIG_SRC']) if os.environ.get('MSSPEC_ORIG_SRC') else None
WORK = CODE.parent / 'platform' / 'work'

TOL = 1e-5
CAP = 12
RESTART = 6
A0 = {'Ag': 4.09, 'Ni': 3.52, 'Cu': 3.61}

WORKER = r'''
import json, resource, sys, time
import numpy as np
cfg = json.loads(sys.argv[1])
sys.path.insert(0, cfg['src'])
from ase.build import bulk
from msspec.calculator import MSSPEC
from msspec.utils import center_cluster, cut_sphere, cut_plane, get_atom_index

a0 = cfg['a0']
cluster = bulk(cfg['elem'], a=a0, cubic=True).repeat((2, 2, 2))
center_cluster(cluster)
cluster = cut_sphere(cluster, radius=a0 + 0.01)
cluster = cut_plane(cluster, z=0.1)
cluster.absorber = get_atom_index(cluster, 0, 0, -a0)

if cfg['mode'] == 'orig_tree':
    calc = MSSPEC(spectroscopy='PED', algorithm='inversion', folder=cfg['folder'], txt=None)
    calc.calculation_parameters.set_parameter('linear_solver', 'tree_gmres')
elif cfg['mode'] == 'new_tree':
    calc = MSSPEC(spectroscopy='PED', algorithm='gmres', folder=cfg['folder'], txt=None)
else:
    calc = MSSPEC(spectroscopy='PED', algorithm='inversion', folder=cfg['folder'], txt=None)
if cfg['mode'] in ('orig_tree', 'new_tree'):
    calc.calculation_parameters.set_parameter('linear_solver_tolerance', cfg['tol'])
    calc.calculation_parameters.set_parameter('linear_solver_max_iterations', cfg['cap'])
    calc.calculation_parameters.set_parameter('linear_solver_restart', cfg['restart'])
calc.set_atoms(cluster)
t0 = time.time()
data = calc.get_theta_scan(level='2p3/2', theta=np.linspace(-70, 70, 71),
                           phi=0, kinetic_energy=300.0)
wall = time.time() - t0
np.save(cfg['out'], np.asarray(data[-1].cross_section))
rss_mib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
print(json.dumps({'wall_s': wall, 'rss_mib': rss_mib}))
'''


def run(mode, src, elem, workdir):
    out = os.path.join(workdir, '%s_%s.npy' % (elem, mode))
    cfg = {'src': str(src), 'mode': mode, 'elem': elem, 'a0': A0[elem],
           'folder': os.path.join(workdir, 'calc_%s_%s' % (elem, mode)),
           'out': out, 'tol': TOL, 'cap': CAP, 'restart': RESTART}
    env = dict(os.environ)
    if mode in ('orig_tree', 'new_tree'):
        env['MSSPEC_FALLBACK_DENSE'] = '0'
    proc = subprocess.run([sys.executable, '-c', WORKER, json.dumps(cfg)],
                          capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError('%s %s failed:\n%s' % (mode, elem, proc.stderr[-2000:]))
    stats = json.loads(proc.stdout.strip().splitlines()[-1])
    stats['iters'] = [int(m) for m in re.findall(r'TREE_GMRES iterations,rhs=\s+(\d+)', proc.stdout)]
    return np.load(out), stats


def main(elems):
    report = {'settings': {'tol': TOL, 'cap': CAP, 'restart': RESTART,
                           'theta_points': 71, 'kinetic_energy_eV': 300.0},
              'materials': {}}
    if ORIG_SRC is None or not ORIG_SRC.is_dir():
        sys.exit('set MSSPEC_ORIG_SRC to the pre-GMRES msspec src folder')
    WORK.mkdir(parents=True, exist_ok=True)
    ok = True
    with tempfile.TemporaryDirectory(prefix='material_check_', dir=WORK) as workdir:
        print('temp dir:', workdir)
        for elem in elems:
            res = {}
            for mode, src in (('orig_direct', ORIG_SRC), ('new_direct', NEW_SRC),
                              ('orig_tree', ORIG_SRC), ('new_tree', NEW_SRC)):
                try:
                    arr, stats = run(mode, src, elem, workdir)
                except RuntimeError as exc:
                    if not mode.endswith('_tree'):
                        raise
                    # Method failure (fallback is off): recorded, never replaced by LU.
                    res[mode] = {'error': str(exc)}
                    continue
                res[mode] = {'arr': arr, **stats}
            failed = {m: v['error'] for m, v in res.items() if 'error' in v}
            if failed:
                ok = False
                report['materials'][elem] = {'method_failures': failed, 'pass': False}
                print('%-3s method failure: %s' % (elem, ', '.join(sorted(failed))))
                continue
            d_ident = np.array_equal(res['new_direct']['arr'], res['orig_direct']['arr'])
            ref = res['orig_tree']['arr']
            t_rel = float(np.max(np.abs(res['new_tree']['arr'] - ref) / np.maximum(np.abs(ref), 1e-300)))
            lu = res['new_direct']['arr']
            rel_lu = {m: float(np.linalg.norm(res[m]['arr'] - lu) / np.linalg.norm(lu))
                      for m in ('orig_tree', 'new_tree')}
            entry = {
                'direct_bit_identical': bool(d_ident),
                'tree_max_rel_diff': t_rel,
                'tree_pass': bool(t_rel <= 1e-12),
                'tree_rel_l2_vs_lu': rel_lu,
                'tree_iterations': {m: res[m]['iters'] for m in ('orig_tree', 'new_tree')},
                'runs': {m: {'wall_s': round(res[m]['wall_s'], 3),
                             'rss_mib': round(res[m]['rss_mib'], 1)}
                         for m in ('orig_direct', 'new_direct', 'orig_tree', 'new_tree')},
            }
            acc_ok = all(v <= 1e-4 for v in rel_lu.values())
            entry['tree_accuracy_pass'] = bool(acc_ok)
            ok = ok and d_ident and acc_ok
            report['materials'][elem] = entry
            print('%-3s direct bit-identical=%s  tree vs orig max rel diff=%.3e (not a gate)' % (
                elem, d_ident, t_rel))
            print('    rel_L2 vs LU: orig=%.3e new=%.3e  accuracy gate=%s' % (
                rel_lu['orig_tree'], rel_lu['new_tree'], acc_ok))
            print('    iterations orig=%s new=%s' % (entry['tree_iterations']['orig_tree'],
                                                 entry['tree_iterations']['new_tree']))
    out_path = HERE / 'material_check_report.json'
    out_path.write_text(json.dumps(report, indent=2) + '\n')
    print('report:', out_path)
    print('PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:] or ['Ag', 'Ni', 'Cu']))
