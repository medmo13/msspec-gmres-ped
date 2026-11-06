"""Resumable scan of rho(G0T) over candidate prod3000 runs.

One run per distinct (name, prototype, surface_orientation) with natoms <= 150,
in a fixed seeded order. Each run is computed in its own process (rho_probe.py),
so a crash or a Fortran STOP does not stop the scan. Rows already in rho_scan.csv
are skipped.

Usage: python platform/rho_scan.py [--limit N]
"""
import argparse, csv, json, os, random, resource, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROD = Path('/home/mohammed/Desktop/REAL-WORK/msspec_ML/prod3000')
PY = ROOT / 'code' / '_venv' / 'bin' / 'python'
OUT = ROOT / 'platform' / 'rho_scan.csv'
SEED = 20261024
NATOMS_MAX = 150
FIELDS = ['run_tag', 'name', 'prototype', 'surface_orientation', 'orientation_class',
          'level', 'absorber_idx', 'phi', 'ke', 'radius', 'natoms', 'lmax',
          'rho', 'wall_s', 'peak_rss_mib', 'rc']


def candidates(order='random'):
    """One run per distinct material. order='random': seeded shuffle.
    order='ke': lowest kinetic energy first (rho falls as ke rises in the first scan),
    so the most likely rho >> 1 runs are computed first.
    """
    rows = list(csv.DictReader(open(PROD / 'data.csv')))
    rows = [r for r in rows if int(r['natoms']) <= NATOMS_MAX]
    if order == 'ke':
        rows.sort(key=lambda r: (float(r['ke']), -int(r['natoms']), r['run_tag']))
    else:
        rng = random.Random(SEED)
        rng.shuffle(rows)
    seen, pool = set(), []
    for r in rows:
        key = (r['name'], r['prototype'], r['surface_orientation'])
        if key not in seen:
            seen.add(key)
            pool.append(r)
    return pool


def done_tags():
    """Runs with a finished row (rc 0). Failed rows stay in the CSV and are retried."""
    if not OUT.exists():
        return set()
    return {r['run_tag'] for r in csv.DictReader(open(OUT)) if r['rc'] == '0'}


def run_one(r):
    run_dir = PROD / 'runs' / r['run_tag']
    proc = subprocess.run([str(PY), str(ROOT / 'platform' / 'rho_probe.py'),
                           str(run_dir), r['ke'], r['lmax']],
                          capture_output=True, text=True)
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024.0
    out = {k: r.get(k, '') for k in FIELDS}
    out['rc'] = proc.returncode
    out['peak_rss_mib'] = round(rss, 1)
    lines = [l for l in proc.stdout.splitlines() if l.startswith('{')]
    if proc.returncode == 0 and lines:
        j = json.loads(lines[-1])
        out['rho'] = j['rho']
        out['wall_s'] = j['wall_s']
    else:
        out['rho'] = ''
        out['wall_s'] = ''
        (ROOT / 'platform' / 'work').mkdir(parents=True, exist_ok=True)
        with open(ROOT / 'platform' / 'work' / 'rho_errors.log', 'a') as fh:
            fh.write('=== %s rc=%d\n%s\n' % (r['run_tag'], proc.returncode,
                                             (proc.stdout + proc.stderr)[-2000:]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--order', choices=['random', 'ke'], default='random')
    args = ap.parse_args()
    pool = candidates(args.order)
    if args.limit:
        pool = pool[:args.limit]
    skip = done_tags()
    new_file = not OUT.exists()
    with open(OUT, 'a', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if new_file:
            w.writeheader()
        for i, r in enumerate(pool):
            if r['run_tag'] in skip:
                continue
            row = run_one(r)
            w.writerow(row)
            fh.flush()
            print('[%d/%d] %s rho=%s rc=%s' % (i + 1, len(pool), r['run_tag'],
                                              row['rho'], row['rc']), flush=True)


if __name__ == '__main__':
    main()
