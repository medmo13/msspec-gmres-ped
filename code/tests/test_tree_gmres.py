"""Regression tests for the tree-preconditioned GMRES solver (PED only).

Each case runs in its own subprocess, because a Fortran STOP ends the process.
Run with:  python tests/test_tree_gmres.py   (or pytest tests/)

Cases:
  lu            dense LU reference, algorithm='inversion'
  tree          tree-GMRES, default settings
  tree_dump     tree-GMRES with MSSPEC_TREE_ITER_DUMP=1 (must equal 'tree' exactly)
  forced_fail   cap 1 iteration, tol 1e-12: tree fails. No dense fallback exists,
                so the process stops with status none_available (method failure).
"""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'msspec_gmres' / 'src'
# Work files stay in the project (project rule: never /tmp or the harness scratchpad).
WORK = ROOT.parent / 'platform' / 'work'

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

if cfg['algorithm'] == 'gmres':
    calc = MSSPEC(spectroscopy='PED', algorithm='gmres', folder=cfg['folder'], txt=None)
    calc.calculation_parameters.set_parameter('linear_solver_tolerance', cfg['tol'])
    calc.calculation_parameters.set_parameter('linear_solver_max_iterations', cfg['cap'])
    calc.calculation_parameters.set_parameter('linear_solver_restart', cfg['restart'])
    if cfg['fixed'] is not None:
        calc.calculation_parameters.set_parameter('linear_solver_fixed_steps', cfg['fixed'])
    if cfg['axial'] is False:
        calc.calculation_parameters.set_parameter('linear_solver_axial_translation', False)
    if cfg['dump']:
        calc.calculation_parameters.set_parameter('linear_solver_iterate_dump', True)
        calc.calculation_parameters.set_parameter('linear_solver_verify', True)
else:
    calc = MSSPEC(spectroscopy='PED', algorithm='inversion', folder=cfg['folder'], txt=None)
calc.set_atoms(cluster)
data = calc.get_theta_scan(level='2p3/2', theta=np.linspace(-70, 70, 15),
                           phi=0, kinetic_energy=300.0)
np.save(cfg['out'], np.asarray(data[-1].cross_section))
'''

TOL = 1e-12
CAP = 1


def run_case(name, algorithm, workdir, env_extra=None, tol=1e-2, cap=12, restart=6, fixed=None, axial=True, dump=False):
    out = os.path.join(workdir, name + '.npy')
    folder = os.path.join(workdir, 'calc_' + name)
    cfg = {'src': str(SRC), 'algorithm': algorithm, 'folder': folder,
           'out': out, 'tol': tol, 'cap': cap, 'restart': restart,
           'fixed': fixed, 'axial': axial, 'dump': dump}
    env = dict(os.environ)
    env.pop('MSSPEC_TREE_ITER_DUMP', None)
    env.pop('MSSPEC_TREE_STOP_AFTER', None)
    env.pop('MSSPEC_FALLBACK_DENSE', None)
    env.update(env_extra or {})
    proc = subprocess.run([sys.executable, '-c', WORKER, json.dumps(cfg)],
                          env=env, capture_output=True, text=True)
    stdout = proc.stdout + proc.stderr
    arr = np.load(out) if proc.returncode == 0 else None
    return proc.returncode, stdout, arr, folder


def rel_l2(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def _read_vec(lines, start, n):
    """Read n complex values written as 're im' pairs from lines[start:]."""
    vals = [complex(*map(float, lines[start + k].split())) for k in range(n)]
    return np.array(vals)


def parse_dump(path):
    """Return (iterates, final) from msspec_tree_iter.dat.

    iterates[rhs] = [(iter, resid, x), ...] and final[rhs] = x, where x is the
    complex vector of one block.
    """
    with open(path) as fh:
        lines = fh.read().splitlines()
    iterates, final = {}, {}
    i = 0
    while i < len(lines):
        head = lines[i].split()
        if head[0] == 'ITERATE':
            _, rhs, it, _step, n, resid = head
            x = _read_vec(lines, i + 1, int(n))
            iterates.setdefault(int(rhs), []).append((int(it), float(resid), x))
            i += 1 + int(n)
        elif head[0] == 'FINAL':
            _, rhs, n = head
            final[int(rhs)] = _read_vec(lines, i + 1, int(n))
            i += 1 + int(n)
        else:
            i += 1
    return iterates, final


def main():
    results = {}
    WORK.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='tree_gmres_test_', dir=WORK) as workdir:
        rc, out, lu, _ = run_case('lu', 'inversion', workdir)
        assert rc == 0, out
        rc, out, tree, _ = run_case('tree', 'gmres', workdir)
        assert rc == 0, out
        err = rel_l2(tree, lu)
        assert err < 1e-2, 'tree vs LU rel_L2=%g' % err
        results['tree_vs_lu'] = err

        rc, out, tree_dump, folder = run_case(
            'tree_dump', 'gmres', workdir, env_extra={'MSSPEC_TREE_ITER_DUMP': '1'})
        assert rc == 0, out
        assert np.array_equal(tree_dump, tree), 'dump changes the solution'
        results['dump_bit_identical'] = True

        rc, out, py_dump, py_folder = run_case('py_dump', 'gmres', workdir, dump=True)
        assert rc == 0, out
        assert np.array_equal(py_dump, tree), 'Python dump/verify parameters change the solution'
        assert os.path.isfile(os.path.join(py_folder, 'output', 'msspec_tree_iter.dat')), 'Python dump parameter wrote no file'
        assert 'TREE_SERIES' in out, 'Python verify parameter printed no residual'
        results['python_dump_verify'] = True

        dump = os.path.join(folder, 'output', 'msspec_tree_iter.dat')
        assert os.path.isfile(dump), 'dump file missing: ' + dump
        iterates, final = parse_dump(dump)
        assert iterates, 'dump has no ITERATE blocks'
        assert set(final) == set(iterates), 'FINAL blocks do not match RHS list'
        worst = 0.0
        for rhs, seq in iterates.items():
            assert seq[0][0] == 0, 'rhs %d first block is iteration %d' % (rhs, seq[0][0])
            # GMRES with x0 = 0 minimises over alpha*M^-1 b at step 1, so the
            # residual cannot grow from iteration 0 on (same absolute unit).
            res = [r for _, r, _ in seq]
            assert all(b <= a * (1 + 1e-6) for a, b in zip(res, res[1:])), \
                'rhs %d residual grows: %s' % (rhs, res)
            last_resid = seq[-1][1]
            assert last_resid < 1e-2, 'rhs %d last resid %g' % (rhs, last_resid)
            # The last iterate of each RHS must be the solution the solver returns.
            scale = max(np.max(np.abs(final[rhs])), 1e-300)
            worst = max(worst, np.max(np.abs(seq[-1][2] - final[rhs])) / scale)
        assert worst <= 1e-14, 'last iterate differs from returned X: %g' % worst
        results['dump_rhs'] = len(iterates)
        results['last_iterate_vs_returned'] = worst

        # E16 switch: STOP_AFTER=k must return exactly the dumped iterate k,
        # and the observable must come from that iterate.
        for k in (0, 2):
            rc, out, arr_k, folder_k = run_case(
                'stop%d' % k, 'gmres', workdir,
                env_extra={'MSSPEC_TREE_ITER_DUMP': '1',
                           'MSSPEC_TREE_STOP_AFTER': str(k)})
            assert rc == 0, out
            assert 'stop_after=' in out, 'stop mode not announced'
            _, fin_k = parse_dump(os.path.join(folder_k, 'output', 'msspec_tree_iter.dat'))
            assert set(fin_k) == set(iterates), 'RHS list differs in stop mode'
            worst_k = 0.0
            for rhs, x in fin_k.items():
                ref = [v for (i, _, v) in iterates[rhs] if i == k]
                assert ref, 'rhs %d has no iterate %d in the normal run' % (rhs, k)
                scale = max(np.max(np.abs(ref[0])), 1e-300)
                worst_k = max(worst_k, np.max(np.abs(x - ref[0])) / scale)
            assert worst_k <= 1e-12, 'stop %d differs from dumped iterate: %g' % (k, worst_k)
            assert not np.array_equal(arr_k, tree), 'stop %d equals the converged result' % k
            results['stop%d_vs_dumped_iterate' % k] = worst_k
            results['stop%d_vs_lu' % k] = rel_l2(arr_k, lu)
        for bad in ('99', 'abc'):
            rc, out, _, _ = run_case('stopbad', 'gmres', workdir,
                                     env_extra={'MSSPEC_TREE_STOP_AFTER': bad})
            assert rc != 0 and 'none_available' in out, 'bad value %r not refused' % bad
        results['stop_bad_values_refused'] = True

        # The same switch through the Python parameter linear_solver_fixed_steps.
        for k in (0, 1, 2, 3):
            rc, out, arr_p, _ = run_case('fixed%d' % k, 'gmres', workdir, fixed=k)
            assert rc == 0, out
            rc, out, arr_e, _ = run_case('fixedenv%d' % k, 'gmres', workdir,
                                         env_extra={'MSSPEC_TREE_STOP_AFTER': str(k)})
            assert rc == 0, out
            assert np.array_equal(arr_p, arr_e), 'fixed_steps=%d differs from the environment switch' % k
            results['fixed_steps%d_vs_lu' % k] = rel_l2(arr_p, lu)
        rc, out, arr_g, _ = run_case('gaunt', 'gmres', workdir, axial=False)
        assert rc == 0, out
        results['direct_gaunt_vs_axial'] = rel_l2(arr_g, tree)
        assert results['direct_gaunt_vs_axial'] < 1e-6, results['direct_gaunt_vs_axial']

        rc, out, _, _ = run_case('forced_fail', 'gmres', workdir,
                                 tol=TOL, cap=CAP)
        assert rc != 0, 'expected a method failure (no fallback exists)'
        assert 'none_available' in out, out
        assert 'dense_lu' not in out, out
        results['forced_fail_exit'] = rc
    return results


def test_tree_gmres():
    main()


# One case per name, with the arguments main() uses for it.
SINGLE_CASES = {
    'lu': dict(algorithm='inversion'),
    'tree': dict(algorithm='gmres'),
    'tree_dump': dict(algorithm='gmres', env_extra={'MSSPEC_TREE_ITER_DUMP': '1'}),
    'forced_fail': dict(algorithm='gmres', tol=TOL, cap=CAP),
}


def run_one(name):
    """Run one case alone and check only what that case must show."""
    spec = SINGLE_CASES[name]
    kwargs = dict(env_extra=spec.get('env_extra'), tol=spec.get('tol', 1e-2),
                  cap=spec.get('cap', 12))
    WORK.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='tree_gmres_one_', dir=WORK) as workdir:
        rc, out, arr, _ = run_case(name, spec['algorithm'], workdir, **kwargs)
        if name == 'forced_fail':
            assert rc != 0, 'expected a method failure (no fallback exists)'
            assert 'none_available' in out, out
            assert 'dense_lu' not in out, out
        else:
            assert rc == 0, out
        if name in ('tree', 'tree_dump'):
            _, lu_out, lu, _ = run_case('lu', 'inversion', workdir)
            assert lu is not None, lu_out
            err = rel_l2(arr, lu)
            print('rel_L2 vs LU = %g' % err)
    print(out.strip().splitlines()[-1] if out.strip() else '')
    print('PASS: %s' % name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tree-GMRES regression tests.')
    parser.add_argument('--case', choices=sorted(SINGLE_CASES),
                        help='run only this case (default: all cases)')
    args = parser.parse_args()
    if args.case:
        run_one(args.case)
    else:
        res = main()
        for key, val in res.items():
            print('%-20s %s' % (key, val))
        print('PASS: tree-GMRES regression tests')
