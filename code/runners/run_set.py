"""Resumable runner for the benchmark sets (PROTOCOL.md A2.11, PLAN.md section 10).

One run = one (case, method, repetition) in a fresh process (worker.py). Every run writes
at once, whatever the outcome:
    <out>/<set>/<case>/<method>_<rep>/row.json    the result row
    <out>/<set>/<case>/<method>_<rep>/result.npz  the arrays (only if the run gave a result)
    <out>/<set>/<case>/<method>_<rep>/log.txt     full stdout and stderr (git-ignored)
and appends the row to <out>/<set>/rows.csv. A run whose row.json exists is skipped.

Methods: lu | gmres | nopre (GMRES with M = I, D40) | k0..k4 (GMRES stopped after exactly k steps, E16) | series1..series8.
Gates, all recorded in the row:
  - the 1-minute load average must be <= --max-load at the start of each run (refuse, or wait
    up to --wait-load minutes);
  - the source and binary hashes must equal code/BUILD_MANIFEST.json and start with --expect-hash;
  - dense LU: predicted memory 0.667 MiB * atoms^2 (A2.9, amended by D28: swap no longer
    blocks it). It runs when prediction + 3 GiB < MemAvailable, and it is killed when
    MemAvailable falls below 2 GiB. Otherwise the row is 'lu_not_run' (not a method failure).
A failed GMRES solve (non-zero return code, 'none_available') is an error row, never replaced.
--warmup runs every method once per case first, in the same rotated order, with rep label 'w'
(folder <method>_w, ignored by the evaluator). Without it, times are not benchmark times (A2.8).

Usage: python code/runners/run_set.py --set NAME --methods lu,gmres --scan theta:-90:89:180 \
           [--mat-ids H01,L02 | --stratum high] [--reps 1] [--warmup] [--dry-run]
"""
import argparse
import csv
import json
import os
import re
import signal
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

CODE = Path(__file__).resolve().parents[1]
ROOT = CODE.parent
sys.path.insert(0, str(CODE / 'tests'))
PROD = Path('/home/mohammed/Desktop/REAL-WORK/msspec_ML/prod3000')
WORK = ROOT / 'platform' / 'work'
LU_MIB_PER_ATOM2 = 0.667
LU_FREE_LIMIT_MIB = 8192
LU_MARGIN_MIB = 3072
LU_STOP_MIB = 2048
THREADS = 8

FIELDS = ['series_src', 'set', 'case', 'run_tag', 'stratum', 'method', 'rep', 'order_pos', 'run_status', 'error',
          'error_reason', 'rc', 'status_line', 'wall_s', 'worker_calc_s', 'ms_tree_setup_s',
          'ms_solve_s', 'peak_rss_mib', 'iters_min', 'iters_max', 'iters_list', 'iters_equal_k', 'n_nonfinite',
          'natoms', 'ke', 'rho', 'lmax', 'use_debye_model', 'scan', 'tol', 'cap', 'restart',
          'pred_lu_mib', 'mem_avail_before_mib', 'swap_on', 'load_before', 'load_after',
          'governor', 'git_head', 'source_hash', 'start_utc']


def meminfo():
    d = {}
    for line in open('/proc/meminfo'):
        k, v = line.split(':')
        d[k] = int(v.split()[0]) / 1024.0
    return d


def build_hashes():
    import write_build_manifest as wbm
    _, src = wbm.source_hashes()
    return src, wbm.binary_hashes()


SERIES_SRC = '/home/mohammed/Desktop/REAL-WORK/msspec_ML/FINAL/07_msspec_dev/msspec_python3/src'
SERIES_SO_SHA = 'e009c24cbdbaec58'     # first 16 hex of sha256 of _phd_se_noso_nosp_nosym.so (D25)


def check_series_src(src):
    import hashlib
    so = Path(src) / 'msspec/spec/fortran/_phd_se_noso_nosp_nosym.so'
    h = hashlib.sha256(so.read_bytes()).hexdigest()
    if h[:16] != SERIES_SO_SHA:
        sys.exit('REFUSED: series binary %s differs from D25 (%s)' % (h[:16], SERIES_SO_SHA))


def check_build(expect):
    manifest = json.loads((CODE / 'BUILD_MANIFEST.json').read_text())
    src, bins = build_hashes()
    if src != manifest['source_tree_sha256']:
        sys.exit('REFUSED: source hash %s differs from BUILD_MANIFEST.json %s'
                 % (src[:12], manifest['source_tree_sha256'][:12]))
    if bins != manifest['binary_files_sha256']:
        sys.exit('REFUSED: a binary differs from BUILD_MANIFEST.json (rebuild, then write the manifest)')
    if not src.startswith(expect):
        sys.exit('REFUSED: source hash %s does not start with the expected %s' % (src[:12], expect))
    return src


def parse_scan(text):
    p = text.split(':')
    if p[0] == 'theta' and len(p) == 4:
        return {'type': 'theta', 'theta': [float(p[1]), float(p[2]), int(p[3])]}
    if p[0] == 'thetaphi' and len(p) == 7:
        if float(p[1]) < 0:
            sys.exit('bad --scan %r: a map starts at theta 0 (negative theta repeats directions)' % text)
        return {'type': 'thetaphi', 'theta': [float(p[1]), float(p[2]), int(p[3])],
                'phi': [float(p[4]), float(p[5]), int(p[6])]}
    if p[0] == 'energy' and len(p) == 5:
        return {'type': 'energy', 'energy': [float(p[1]), float(p[2]), int(p[3])],
                'thetas': [float(t) for t in p[4].split(',')]}
    sys.exit('bad --scan %r: use theta:start:stop:n, thetaphi:t0:t1:nt:p0:p1:np or energy:e0:e1:n:theta1,theta2' % text)


def valid_method(m):
    return m in ('lu', 'gmres', 'nopre') or re.fullmatch(r'k[0-9]+', m) or re.fullmatch(r'series[1-8]', m)


def lu_decision(pred, mi):
    """('run' | 'not_run', reason) for dense LU under A2.9."""
    if pred < LU_FREE_LIMIT_MIB:
        return 'run', ''
    if mi['MemAvailable'] < pred + LU_MARGIN_MIB:
        return 'not_run', 'low_ram (%.0f MiB free, need %.0f)' % (mi['MemAvailable'], pred + LU_MARGIN_MIB)
    return 'run', ''


def wait_for_load(max_load, wait_min):
    t_end = time.time() + 60 * wait_min
    while os.getloadavg()[0] > max_load:
        if time.time() >= t_end:
            sys.exit('REFUSED: load average %.2f above %.2f at start (A2.11)' % (os.getloadavg()[0], max_load))
        time.sleep(30)


def run_env(method, k):
    env = dict(os.environ, TMPDIR=str(WORK), PYTHONPATH='')
    for v in ('MSSPEC_TREE_STOP_AFTER', 'MSSPEC_TREE_ITER_DUMP', 'MSSPEC_FALLBACK_DENSE',
              'MSSPEC_DEBYE_QD', 'MSSPEC_TREE_NOPRECOND'):
        env.pop(v, None)
    if method == 'nopre':
        env['MSSPEC_TREE_NOPRECOND'] = '1'      # M = I: GMRES without the tree preconditioner (D40)
    if method in ('gmres', 'nopre') or k is not None:
        env.update(OMP_NUM_THREADS=str(THREADS), MSSPEC_TREE_ACTION_THREADS=str(THREADS),
                   OPENBLAS_NUM_THREADS='1')
    else:
        env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS=str(THREADS))
    if k is not None:
        env['MSSPEC_TREE_STOP_AFTER'] = str(k)
    return env


def execute(cmd, env, log_path, watch_ram):
    """Run cmd; return (rc, wall_s, peak_rss_mib, killed_reason). Peak RSS from wait4."""
    t0 = time.time()
    with open(log_path, 'w') as log:
        proc = subprocess.Popen(cmd, env=env, stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True)
        killed = ''
        last_mem = 0.0
        while True:
            pid, status, ru = os.wait4(proc.pid, os.WNOHANG)
            if pid:
                break
            if watch_ram and time.time() - last_mem > 1.0:
                last_mem = time.time()
                if meminfo()['MemAvailable'] < LU_STOP_MIB and not killed:
                    os.killpg(proc.pid, signal.SIGKILL)
                    killed = 'stopped_low_ram'
            time.sleep(0.05)
    rc = os.waitstatus_to_exitcode(status)
    return rc, time.time() - t0, ru.ru_maxrss / 1024.0, killed


def parse_log(text):
    iters = [int(m.group(1)) for m in re.finditer(r'TREE_GMRES iterations,rhs=\s+(\d+)', text)]
    status = [l.strip() for l in text.splitlines() if 'MSSPEC_SOLVER_STATUS' in l]
    setup = solve = 0.0
    for m in re.finditer(r'MSSPEC_TIMING je,tree_setup_s,solve_s,n,na=\s+\d+\s+(\S+)\s+(\S+)', text):
        setup += float(m.group(1)); solve += float(m.group(2))
    calc_s = re.findall(r'WORKER_CALC_S (\S+)', text)
    return iters, status, setup, solve, (float(calc_s[-1]) if calc_s else '')


def one_run(a, case, method, rep, order_pos, set_dir, scan, src_hash, git_head):
    run_dir = set_dir / case['mat_id'] / ('%s_%s' % (method, rep))
    if (run_dir / 'row.json').exists():
        return None
    run_dir.mkdir(parents=True, exist_ok=True)
    wait_for_load(a.max_load, a.wait_load)
    mi = meminfo()
    k = int(method[1:]) if re.fullmatch(r'k[0-9]+', method) else None
    pred = a.lu_mib_per_atom2 * int(case['natoms']) ** 2
    row = {f: '' for f in FIELDS}
    row.update({'set': a.set, 'case': case['mat_id'], 'run_tag': case['run_tag'],
                'stratum': case['stratum'], 'method': method, 'rep': rep, 'order_pos': order_pos,
                'error': 0, 'natoms': case['natoms'], 'ke': case['ke'], 'rho': case['rho'],
                'lmax': case['lmax'], 'use_debye_model': case['use_debye_model'],
                'scan': a.scan, 'pred_lu_mib': round(pred, 1),
                'mem_avail_before_mib': round(mi['MemAvailable']), 'swap_on': int(mi['SwapTotal'] > 0),
                'load_before': round(os.getloadavg()[0], 2), 'git_head': git_head,
                'source_hash': src_hash[:12], 'start_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                'governor': Path('/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor').read_text().strip()})
    if method in ('gmres', 'nopre') or k is not None:
        row.update(tol=a.tol, cap=a.cap, restart=a.restart)
    watch = False
    if method == 'lu':
        decision, why = lu_decision(pred, mi)
        if decision == 'not_run':
            row.update(run_status='lu_not_run', error_reason=why)
            return finish(row, run_dir, set_dir)
        watch = pred >= LU_FREE_LIMIT_MIB
    use_alt = a.series_src and method.startswith('series')
    row['series_src'] = a.series_src if use_alt else ''
    cfg = {'src': a.series_src if use_alt else str(CODE / 'msspec_gmres' / 'src'), 'platform': str(ROOT / 'platform'),
           'case': case, 'method': method, 'scan': scan, 'run_dir': str(PROD / 'runs' / case['run_tag']),
           'folder': str(run_dir / 'calc'), 'out': str(run_dir / 'result.npz'),
           'tol': a.tol, 'cap': a.cap, 'restart': a.restart}
    env = run_env(method, k)
    env['MSSPEC_DEBYE_QD'] = case['debye_qd']
    rc, wall, rss, killed = execute([sys.executable, '-W', 'ignore', str(CODE / 'runners' / 'worker.py'),
                                     json.dumps(cfg)], env, run_dir / 'log.txt', watch)
    text = (run_dir / 'log.txt').read_text(errors='replace')
    iters, status, setup, solve, calc_s = parse_log(text)
    row.update(rc=rc, wall_s=round(wall, 2), worker_calc_s=calc_s, peak_rss_mib=round(rss, 1),
               ms_tree_setup_s=round(setup, 3), ms_solve_s=round(solve, 3),
               status_line=status[-1] if status else '', load_after=round(os.getloadavg()[0], 2))
    if iters:
        row.update(iters_min=min(iters), iters_max=max(iters), iters_list=' '.join(str(i) for i in iters))
    result = run_dir / 'result.npz'
    reason = killed
    if not reason and rc != 0:
        reason = 'none_available' if any('none_available' in s for s in status) else 'rc=%d' % rc
    if not reason and not result.exists():
        reason = 'no_result'
    if result.exists():
        row['n_nonfinite'] = int(np.sum(~np.isfinite(np.load(result)['cross_section'])))
    if not reason and row.get('n_nonfinite'):
        reason = 'nonfinite_result'
    if k is not None:
        row['iters_equal_k'] = int(bool(iters) and all(i == k for i in iters)) if k > 0 else int(
            bool(iters) and all(i == 0 for i in iters))
        if not reason and not row['iters_equal_k']:
            reason = 'iterations_not_k'
    row.update(run_status='ok' if not reason else ('lu_stopped' if killed else 'method_failure'),
               error=int(bool(reason)), error_reason=reason)
    if method == 'lu' and reason:
        row['run_status'] = 'lu_not_run' if killed else 'lu_failed'
    if not a.keep_calc:
        subprocess.run(['rm', '-rf', str(run_dir / 'calc')])
    return finish(row, run_dir, set_dir)


def finish(row, run_dir, set_dir):
    (run_dir / 'row.json').write_text(json.dumps(row, indent=1))
    rows_csv = set_dir / 'rows.csv'
    new = not rows_csv.exists()
    with open(rows_csv, 'a', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if new:
            w.writeheader()
        w.writerow(row)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--set', required=True)
    ap.add_argument('--cases', default=str(ROOT / 'platform' / 'materials_v2.csv'))
    ap.add_argument('--mat-ids', default='')
    ap.add_argument('--stratum', default='')
    ap.add_argument('--methods', default='lu,gmres')
    ap.add_argument('--reps', type=int, default=1)
    ap.add_argument('--warmup', action='store_true', help="one untimed run per method first, rep label 'w' (A2.8)")
    ap.add_argument('--scan', default='theta:-90:89:180')
    ap.add_argument('--out-root', default=str(ROOT / 'runs'))
    ap.add_argument('--tol', type=float, default=1e-2)
    ap.add_argument('--cap', type=int, default=67)
    ap.add_argument('--restart', type=int, default=6)
    ap.add_argument('--max-load', type=float, default=1.5)
    ap.add_argument('--wait-load', type=float, default=0.0, help='minutes to wait for the load gate')
    ap.add_argument('--expect-hash', default='b64d4c64')
    ap.add_argument('--series-src', default=SERIES_SRC, help='read-only MsSpec src used only for series methods (D25)')
    ap.add_argument('--lu-mib-per-atom2', type=float, default=LU_MIB_PER_ATOM2, help='LU memory model; energy scans below 100 eV measured about 0.12')
    ap.add_argument('--keep-calc', action='store_true')
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()

    methods = a.methods.split(',')
    bad = [m for m in methods if not valid_method(m)]
    if bad:
        sys.exit('unknown methods: %s' % bad)
    scan = parse_scan(a.scan)
    cases = list(csv.DictReader(open(a.cases)))
    if a.mat_ids:
        want = a.mat_ids.split(',')
        missing = [w for w in want if w not in {c['mat_id'] for c in cases}]
        if missing:
            sys.exit('unknown mat ids: %s' % missing)
        cases = [c for c in cases if c['mat_id'] in want]
    if a.stratum:
        cases = [c for c in cases if c['stratum'] == a.stratum]
    set_dir = Path(a.out_root) / a.set
    src_hash = check_build(a.expect_hash)
    if any(m.startswith('series') for m in a.methods.split(',')) and a.series_src:
        check_series_src(a.series_src)
    git_head = subprocess.run(['git', '-C', str(ROOT), 'rev-parse', '--short', 'HEAD'],
                              capture_output=True, text=True).stdout.strip()
    WORK.mkdir(parents=True, exist_ok=True)
    print('set %s: %d cases x %d methods x %d reps; build %s; load %.2f; out %s; tmp %s'
          % (a.set, len(cases), len(methods), a.reps, src_hash[:12], os.getloadavg()[0], set_dir, WORK),
          flush=True)
    n = 0
    for ci, case in enumerate(cases):
        for rep in (['w'] if a.warmup else []) + list(range(a.reps)):
            shift = (ci + (0 if rep == 'w' else rep)) % len(methods)          # alternate the method order (A2.8)
            for pos, method in enumerate(methods[shift:] + methods[:shift]):
                done = (set_dir / case['mat_id'] / ('%s_%s' % (method, rep)) / 'row.json').exists()
                if a.dry_run:
                    extra = ''
                    if method == 'lu':
                        extra = ' LU:' + str(lu_decision(a.lu_mib_per_atom2 * int(case['natoms']) ** 2, meminfo()))
                    print('%s %s_%s pos %d %s%s' % (case['mat_id'], method, rep, pos,
                                                    'done' if done else 'todo', extra))
                    continue
                row = one_run(a, case, method, rep, pos, set_dir, scan, src_hash, git_head)
                if row:
                    n += 1
                    print('%s %s_%s rc=%s %s %s wall=%ss rss=%s MiB it=%s-%s %s'
                          % (row['case'], method, rep, row['rc'], row['run_status'], row['error_reason'],
                             row['wall_s'], row['peak_rss_mib'], row['iters_min'], row['iters_max'],
                             row['status_line'][-40:]), flush=True)
    print('finished; %d runs executed' % n)


if __name__ == '__main__':
    main()
