"""Check platform/materials_v2.csv (PROTOCOL.md A2.4, A2.5). Prints every count it asserts.

Checks: 50 rows per stratum; one row per (name, prototype) and none in both strata;
every run_tag is an rc 0 row of rho_scan.csv with the same rho, natoms, ke; high rho >= 1.2,
low rho < 1, natoms <= 150; per run folder: config.json agrees on lmax, level, theta,
use_debye_model and natoms, cluster.xyz has natoms atoms, and the absorber matches config.json.
With --recompute-rho N, the spectral radius of the first N rows is computed again with
rho_probe.py (power method) and compared with the scan value.

Usage: python platform/check_materials.py [--recompute-rho N]
"""
import argparse
import csv
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'platform'))
CSV_IN = ROOT / 'platform' / 'materials_v2.csv'
SCAN = ROOT / 'platform' / 'rho_scan.csv'
PROD = Path('/home/mohammed/Desktop/REAL-WORK/msspec_ML/prod3000')


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)
    print('ok  ', msg)


def main(n_rho):
    rows = list(csv.DictReader(open(CSV_IN)))
    print('columns:', list(rows[0]))
    scan = {r['run_tag']: r for r in csv.DictReader(open(SCAN)) if r['rc'] == '0'}
    strata = {s: [r for r in rows if r['stratum'] == s] for s in ('high', 'low')}
    check(len(strata['high']) == 50 and len(strata['low']) == 50,
          'rows per stratum: high %d, low %d' % (len(strata['high']), len(strata['low'])))
    mats = {s: {(r['name'], r['prototype']) for r in v} for s, v in strata.items()}
    check(all(len(mats[s]) == 50 for s in mats), 'distinct (name, prototype): 50 in each stratum')
    check(not (mats['high'] & mats['low']), 'materials in both strata: 0')
    check(len({r['run_tag'] for r in rows}) == 100, 'distinct run tags: 100')
    check(all(r['run_tag'] in scan for r in rows), 'every run tag is an rc 0 scan row')
    for r in rows:
        s = scan[r['run_tag']]
        assert float(r['rho']) == float(s['rho']) and r['natoms'] == s['natoms'] \
            and float(r['ke']) == float(s['ke']), r['run_tag']
    print('ok   rho, natoms, ke equal the scan values in all 100 rows')
    check(all(float(r['rho']) >= 1.2 for r in strata['high']), 'high stratum: all rho >= 1.2 (min %.3f)'
          % min(float(r['rho']) for r in strata['high']))
    check(all(float(r['rho']) < 1.0 for r in strata['low']), 'low stratum: all rho < 1 (max %.3f)'
          % max(float(r['rho']) for r in strata['low']))
    check(all(int(r['natoms']) <= 150 for r in rows), 'natoms <= 150 in all rows (max %d)'
          % max(int(r['natoms']) for r in rows))

    import rho_probe
    seen = Counter()
    for r in rows:
        d = PROD / 'runs' / r['run_tag']
        cfg = json.loads((d / 'config.json').read_text())
        assert cfg['natoms'] == int(r['natoms']), r['run_tag']
        atoms = rho_probe.read_cluster(d)            # checks absorber element and position
        assert len(atoms) == int(r['natoms']) and atoms.absorber == int(r['absorber_idx']), r['run_tag']
        assert str(cfg['use_debye_model']) == r['use_debye_model'], r['run_tag']
        assert float(cfg['temperature']) == float(r['temperature']), r['run_tag']
        seen[(cfg['lmax'], cfg['level'] == r['level'], cfg['use_debye_model'],
              tuple(cfg['theta']))] += 1
    print('ok   config.json, cluster.xyz and absorber agree for all 100 rows')
    print('settings seen (lmax, level matches row, use_debye_model, theta):', dict(seen))
    check(len(seen) == 1, 'one settings set across the 100 cases (A2.4)')
    over = [r['mat_id'] for r in rows if float(r['pred_lu_mib']) >= 8192]
    print('info %d cases predict LU >= 8 GiB (A2.9 rule needs swap off and free RAM): high %d, low %d'
          % (len(over), sum(1 for m in over if m[0] == 'H'), sum(1 for m in over if m[0] == 'L')))
    for r in rows[:n_rho]:
        out = subprocess.run([sys.executable, str(ROOT / 'platform' / 'rho_probe.py'),
                              str(PROD / 'runs' / r['run_tag']), r['ke'], r['lmax']],
                             capture_output=True, text=True)
        j = json.loads([l for l in out.stdout.splitlines() if l.startswith('{')][-1])
        check(abs(j['rho'] - float(r['rho'])) < 1e-6,
              'rho recomputed %s: %.9f vs scan %.9f' % (r['mat_id'], j['rho'], float(r['rho'])))
    print('PASS')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--recompute-rho', type=int, default=0)
    main(ap.parse_args().recompute_rho)
