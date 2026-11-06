"""E3 runner: one GMRES run per case with the iterate dump on (PROTOCOL.md E3, PLAN.md C2).

Writes <out>/<case>/gmres_dump/{result.npz, msspec_tree_iter.dat, log.txt}. The dump gives the
absolute preconditioned residual after every step; the cross section after step k comes from
the k-step runs (method k<digit> in run_set.py). Diagnostic only, never timed.

Usage: python code/runners/e3_run.py --out DIR --mat-ids H42 --scan theta:-70:70:100
"""
import argparse, csv, json, os, subprocess, sys
from pathlib import Path

CODE = Path(__file__).resolve().parents[1]
ROOT = CODE.parent
sys.path.insert(0, str(CODE / 'runners'))
import run_set as rs  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('--out', required=True)
ap.add_argument('--mat-ids', required=True)
ap.add_argument('--scan', default='theta:-90:89:180')
ap.add_argument('--tol', type=float, default=1e-2)
ap.add_argument('--cap', type=int, default=12)
ap.add_argument('--restart', type=int, default=6)
a = ap.parse_args()
scan = rs.parse_scan(a.scan)
cases = {c['mat_id']: c for c in csv.DictReader(open(ROOT / 'platform/materials_v2.csv'))}
for mid in a.mat_ids.split(','):
    case = cases[mid]
    d = Path(a.out) / mid / 'gmres_dump'
    d.mkdir(parents=True, exist_ok=True)
    cfg = {'src': str(CODE / 'msspec_gmres/src'), 'platform': str(ROOT / 'platform'), 'case': case,
           'method': 'gmres', 'scan': scan, 'run_dir': str(rs.PROD / 'runs' / case['run_tag']),
           'folder': str(d / 'calc'), 'out': str(d / 'result.npz'),
           'tol': a.tol, 'cap': a.cap, 'restart': a.restart}
    env = rs.run_env('gmres', None)
    env['MSSPEC_DEBYE_QD'] = case['debye_qd']
    env['MSSPEC_TREE_ITER_DUMP'] = '1'
    env['TMPDIR'] = str(rs.WORK)
    rc, wall, rss, _ = rs.execute([sys.executable, '-W', 'ignore', str(CODE / 'runners/worker.py'),
                                   json.dumps(cfg)], env, d / 'log.txt', False)
    print(mid, 'rc', rc, 'wall %.1f s' % wall, 'dump files:',
          [str(p.relative_to(d)) for p in d.rglob('msspec_tree_iter.dat')], flush=True)
