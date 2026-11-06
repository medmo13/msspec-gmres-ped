"""Build platform/materials_100.csv: 50 high-rho cases and 50 low-rho cases.

High group: the 50 highest-rho cases of msspec_SE high_rho/ (platform/high_group_top50.csv,
read-only source, see DEVIATIONS.md D11). No atom-count cap for this group.
Low group: prod3000 runs from platform/rho_scan.csv with rc 0 and rho < 1, at most
150 atoms, one run per distinct (name, prototype), taken in scan order (lowest ke first).

If the low group has fewer than 50 distinct materials, the table is written with the
shortfall and the count is printed; the shortfall is recorded in DEVIATIONS.md.

Usage: python platform/select.py
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HIGH = ROOT / 'platform' / 'high_group_top50.csv'
SCAN = ROOT / 'platform' / 'rho_scan.csv'
OUT = ROOT / 'platform' / 'materials_100.csv'
FIELDS = ['mat_id', 'run_tag', 'formula', 'name', 'prototype', 'surface_orientation',
          'orientation_class', 'level', 'absorber_idx', 'ke', 'phi', 'radius', 'natoms',
          'lmax', 'rho_G0T', 'group', 'source']
NLOW = 50


def high_rows():
    rows = []
    for r in csv.DictReader(open(HIGH)):
        rows.append({'run_tag': r['run_dir'], 'formula': '', 'name': r['material'],
                     'prototype': '', 'surface_orientation': '', 'orientation_class': '',
                     'level': '', 'absorber_idx': '', 'ke': r['ke'], 'phi': '',
                     'radius': '', 'natoms': r['natoms'], 'lmax': '', 'rho_G0T': r['rho'],
                     'group': 'high', 'source': 'msspec_SE high_rho (D11)'})
    return rows


def low_rows():
    seen, rows = set(), []
    for r in csv.DictReader(open(SCAN)):
        if r['rc'] != '0' or float(r['rho']) >= 1.0 or int(r['natoms']) > 150:
            continue
        key = (r['name'], r['prototype'])
        if key in seen:
            continue
        seen.add(key)
        rows.append({'run_tag': r['run_tag'], 'formula': '', 'name': r['name'],
                     'prototype': r['prototype'],
                     'surface_orientation': r['surface_orientation'],
                     'orientation_class': r['orientation_class'], 'level': r['level'],
                     'absorber_idx': r['absorber_idx'], 'ke': r['ke'], 'phi': r['phi'],
                     'radius': r['radius'], 'natoms': r['natoms'], 'lmax': r['lmax'],
                     'rho_G0T': r['rho'], 'group': 'low', 'source': 'prod3000 rho_scan'})
        if len(rows) == NLOW:
            break
    return rows


def main():
    high = high_rows()
    low = low_rows()
    rows = high + low
    for i, r in enumerate(rows, 1):
        r['mat_id'] = 'M%03d' % i
    with open(OUT, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print('high rows: %d' % len(high))
    print('low rows:  %d (target %d)' % (len(low), NLOW))
    print('total:     %d' % len(rows))
    if len(low) < NLOW:
        print('SHORTFALL: %d low-rho materials missing' % (NLOW - len(low)))


if __name__ == '__main__':
    main()
