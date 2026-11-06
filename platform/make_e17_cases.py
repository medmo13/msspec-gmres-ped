"""Cases of the energy series (E17): four clusters of materials_v2.csv x the energy grid of PROTOCOL A3.3.

Writes platform/e17_cases.csv (same columns as materials_v2.csv; mat_id = <id>_e<ke>, ke = grid energy, rho left empty)
and platform/work/e17_probe_cmds.txt (one 'run_tag ke lmax mat_id' per line for the rho probe).
Choice (sizes 26 to 106 atoms, four orientation classes, emitter depth proxies from platform/select_e17.py;
nothing above 110 atoms, so dense LU always fits in memory):
  H42 Ru hcp 0001 polar, 26 atoms, coordination 10, radial 0.38
  H17 PbB rocksalt 011 low_index, 48 atoms, coordination 4, radial 0.82
  H28 TiO2 rutile m122 high_index, 85 atoms, coordination 2, radial 0.74
  H48 Hf hcp 10m1m1 semipolar, 106 atoms, coordination 12, radial 0.33
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHOSEN = ['H42', 'H17', 'H28', 'H48']
GRID = [10, 20, 50, 100, 150, 200]
rows = {r['mat_id']: r for r in csv.DictReader(open(ROOT / 'platform' / 'materials_v2.csv'))}
fields = list(next(iter(rows.values())))
out, cmds = [], []
for m in CHOSEN:
    for ke in GRID:
        r = dict(rows[m]); r['mat_id'] = '%s_e%03d' % (m, ke); r['ke'] = str(float(ke)); r['rho'] = ''
        out.append(r); cmds.append('%s %d %s %s' % (r['run_tag'], ke, r['lmax'], r['mat_id']))
with open(ROOT / 'platform' / 'e17_cases.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(out)
(ROOT / 'platform' / 'work' / 'e17_probe_cmds.txt').write_text('\n'.join(cmds) + '\n')
print(len(out), 'cases')
