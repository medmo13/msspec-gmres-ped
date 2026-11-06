"""Choose the four clusters of the energy series (E17): list size, emitter depth proxies and orientation class.

Coordination = atoms within 1.2 x the nearest-neighbour distance of the emitter; radial = distance of the emitter
from the cluster centre / largest atom distance from the centre (same proxies as PROTOCOL A3.2b).
Prints the table of all 100 cases; writes platform/work/e17_candidates.csv.
"""
import csv, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'platform'))
from rho_probe import read_cluster
PROD = Path('/home/mohammed/Desktop/REAL-WORK/msspec_ML/prod3000/runs')

rows = list(csv.DictReader(open(ROOT / 'platform' / 'materials_v2.csv')))
out = []
for r in rows:
    at = read_cluster(PROD / r['run_tag'])
    p = at.get_positions(); e = p[at.absorber]
    d = np.linalg.norm(p - e, axis=1); d[at.absorber] = np.inf
    coord = int(np.sum(d < 1.2 * d.min()))
    c = p.mean(axis=0)
    rad = float(np.linalg.norm(e - c) / np.max(np.linalg.norm(p - c, axis=1)))
    out.append(dict(mat_id=r['mat_id'], name=r['name'], prototype=r['prototype'], orient=r['surface_orientation'],
                    oclass=r['orientation_class'], level=r['level'], natoms=int(r['natoms']), ke=float(r['ke']),
                    rho=float(r['rho']), coord=coord, radial=round(rad, 2)))
with open(ROOT / 'platform/work/e17_candidates.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
for o in sorted(out, key=lambda o: o['natoms']):
    print(o['mat_id'], o['name'], o['prototype'], o['orient'], o['oclass'], o['level'], o['natoms'], o['ke'], round(o['rho'], 2), o['coord'], o['radial'])
