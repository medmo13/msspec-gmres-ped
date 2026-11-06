"""Cases CSV for the large low-rho clusters (above the dense-LU memory limit, series converges).

Reads platform/work/big_probe_list.csv (prod3000 runs of 250+ atoms, light elements) and
platform/work/big_rho.jsonl (rho of each, from rho_probe.py). Keeps the clusters with rho below
--max-rho, lowest rho first, and writes platform/materials_big.csv in the materials_v2.csv format
(the runner's --cases). Prints the table it wrote.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROD = Path('/home/mohammed/Desktop/REAL-WORK/msspec_ML/prod3000')
sys.path.insert(0, str(ROOT / 'platform'))
import select2  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('--max-rho', type=float, default=0.7)
ap.add_argument('--n', type=int, default=5)
a = ap.parse_args()
rho = {j['run_tag']: j for j in (json.loads(l) for l in open(ROOT / 'platform/work/big_rho.jsonl') if l.startswith('{"run_tag"'))}
probe = [r for r in csv.DictReader(open(ROOT / 'platform/work/big_probe_list.csv')) if r['run_tag'] in rho]
probe = sorted([r for r in probe if rho[r['run_tag']]['rho'] < a.max_rho], key=lambda r: rho[r['run_tag']]['rho'])[:a.n]
debye_qd = select2.load_debye_qd()
out = []
for k, r in enumerate(probe):
    d = PROD / 'runs' / r['run_tag']
    cfg = json.loads((d / 'config.json').read_text())
    spec = json.loads((d / 'sample.json').read_text())['struct']
    n = int(r['natoms'])
    out.append({'mat_id': 'B%02d' % (k + 1), 'stratum': 'low', 'run_tag': r['run_tag'], 'name': r['name'],
                'prototype': r['prototype'], 'surface_orientation': cfg['surface_orientation'],
                'orientation_class': '', 'level': r['level'], 'absorber_idx': cfg['absorber_idx'],
                'phi': cfg.get('phi', r['run_tag'].split('_phi_')[1]), 'ke': r['ke'], 'natoms': n, 'lmax': r['lmax'],
                'rho': rho[r['run_tag']]['rho'], 'temperature': cfg['temperature'],
                'debye_temperature': cfg['debye_temperature'], 'use_debye_model': cfg['use_debye_model'],
                'debye_qd': '%.6f' % debye_qd(spec), 'pred_lu_mib': round(select2.LU_MIB_PER_ATOM2 * n * n, 1),
                'seed': select2.SEED})
with open(ROOT / 'platform/materials_big.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=select2.FIELDS); w.writeheader(); w.writerows(out)
print('wrote platform/materials_big.csv; columns:', select2.FIELDS)
for o in out:
    print(o['mat_id'], o['name'], o['prototype'], o['level'], o['natoms'], 'atoms rho', round(float(o['rho']), 3), 'ke', o['ke'], 'pred LU MiB', o['pred_lu_mib'])
