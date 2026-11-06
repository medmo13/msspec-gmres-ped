"""Build platform/materials_v2.csv: the P100 selection of PROTOCOL.md A2.5.

Source: the rc 0 rows of platform/rho_scan.csv. Each case is one prod3000 run folder
(read-only, config.json and sample.json inside).

High stratum: rho >= 1.2, at most 150 atoms. Low stratum: rho < 1, at most 150 atoms.
One row per distinct material (name, prototype) in each stratum, and no material in
both strata (the high stratum is chosen first; its materials leave the low pool).
Within a material the row is drawn at random (seeded). The 50 materials of a stratum are
then chosen round-robin over prototypes, and inside a prototype round-robin over
orientation classes, each list shuffled with the seed. So every prototype and orientation
class that exists in the pool is represented before any is repeated.

Also writes what the runner needs, so it never imports the sibling code: temperature,
Debye temperature, Debye wave vector (from the sibling materials.py, read-only), the
Debye flag, level, phi, and the predicted dense-LU memory.

Usage: python platform/select2.py
"""
import csv
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCAN = ROOT / 'platform' / 'rho_scan.csv'
OUT = ROOT / 'platform' / 'materials_v2.csv'
PROD = Path('/home/mohammed/Desktop/REAL-WORK/msspec_ML/prod3000')
SEED = 20261030
NPER = 50
RHO_HIGH = 1.2
RHO_LOW = 1.0
MAX_ATOMS = 150
LU_MIB_PER_ATOM2 = 0.667     # PROTOCOL.md A2.9

FIELDS = ['mat_id', 'stratum', 'run_tag', 'name', 'prototype', 'surface_orientation',
          'orientation_class', 'level', 'absorber_idx', 'phi', 'ke', 'natoms', 'lmax', 'rho',
          'temperature', 'debye_temperature', 'use_debye_model', 'debye_qd',
          'pred_lu_mib', 'seed']


def load_debye_qd():
    sys.path.insert(0, str(PROD))
    import materials
    return materials.debye_qd


def pool(rows, test):
    """Rows passing `test`, one random row per (name, prototype)."""
    rng = random.Random(SEED)
    by_mat = defaultdict(list)
    for r in rows:
        if test(r):
            by_mat[(r['name'], r['prototype'])].append(r)
    return {m: rng.choice(rs) for m, rs in sorted(by_mat.items())}


def choose(mats, n, rng):
    """n materials, round-robin over prototypes, then over orientation classes."""
    tree = defaultdict(lambda: defaultdict(list))
    for m, r in sorted(mats.items()):
        tree[r['prototype']][r['orientation_class']].append(m)
    for proto in tree.values():
        for lst in proto.values():
            rng.shuffle(lst)
    protos = sorted(tree)
    rng.shuffle(protos)
    cls_order = {p: sorted(tree[p]) for p in protos}
    for p in protos:
        rng.shuffle(cls_order[p])
    picked, turn = [], defaultdict(int)
    while len(picked) < n and any(tree[p][c] for p in protos for c in tree[p]):
        for p in protos:
            if len(picked) >= n:
                break
            order = cls_order[p]
            for k in range(len(order)):
                c = order[(turn[p] + k) % len(order)]
                if tree[p][c]:
                    picked.append(tree[p][c].pop())
                    turn[p] = (turn[p] + k + 1) % len(order)
                    break
    return picked


def main():
    rows = [r for r in csv.DictReader(open(SCAN)) if r['rc'] == '0']
    rejected = [r for r in csv.DictReader(open(SCAN)) if r['rc'] != '0']
    rng = random.Random(SEED)
    high_pool = pool(rows, lambda r: float(r['rho']) >= RHO_HIGH and int(r['natoms']) <= MAX_ATOMS)
    high = choose(high_pool, NPER, rng)
    low_pool = pool(rows, lambda r: float(r['rho']) < RHO_LOW and int(r['natoms']) <= MAX_ATOMS
                    and (r['name'], r['prototype']) not in set(high))
    low = choose(low_pool, NPER, rng)
    print('scan rc 0 rows %d, rc != 0 rows %d (excluded, listed): %s'
          % (len(rows), len(rejected), [(r['name'], r['natoms']) for r in rejected]))
    print('high pool %d materials -> chosen %d; low pool %d materials -> chosen %d'
          % (len(high_pool), len(high), len(low_pool), len(low)))
    debye_qd = load_debye_qd()
    out = []
    for stratum, picks, src in (('high', high, high_pool), ('low', low, low_pool)):
        for m in picks:
            r = src[m]
            d = PROD / 'runs' / r['run_tag']
            cfg = json.loads((d / 'config.json').read_text())
            spec = json.loads((d / 'sample.json').read_text())['struct']
            n = int(r['natoms'])
            out.append({
                'mat_id': '%s%02d' % (stratum[0].upper(), len(out) % NPER + 1),
                'stratum': stratum, 'run_tag': r['run_tag'], 'name': r['name'],
                'prototype': r['prototype'], 'surface_orientation': r['surface_orientation'],
                'orientation_class': r['orientation_class'], 'level': r['level'],
                'absorber_idx': r['absorber_idx'], 'phi': r['phi'], 'ke': r['ke'],
                'natoms': n, 'lmax': r['lmax'], 'rho': r['rho'],
                'temperature': cfg['temperature'], 'debye_temperature': cfg['debye_temperature'],
                'use_debye_model': cfg['use_debye_model'], 'debye_qd': '%.6f' % debye_qd(spec),
                'pred_lu_mib': round(LU_MIB_PER_ATOM2 * n * n, 1), 'seed': SEED})
    with open(OUT, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(out)
    print('wrote', OUT, len(out), 'rows')


if __name__ == '__main__':
    main()
