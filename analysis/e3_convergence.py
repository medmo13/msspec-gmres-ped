"""E3: observable error against the preconditioned residual, step by step.

Per case it reads
  <set>/<case>/lu_0/result.npz            dense LU scan (reference)
  <set>/<case>/k<j>_0/result.npz          GMRES stopped after exactly j steps (j = 0..6)
  <set>/<case>/gmres_0/result.npz         GMRES as in the paper (tol, cap, restart from PROTOCOL.md)
  <e3>/<case>/gmres_dump/calc/output/msspec_tree_iter.dat   absolute preconditioned residual per step
and writes analysis/out/<tag>_curves.csv and figures/<tag>.png/.pdf.

The relative residual of step j for one right-hand side is RESID_j / ||M^-1 b||, where
||M^-1 b|| is the norm of the iteration-0 vector in the dump (it is the norm used by the solver
stopping rule). Several right-hand sides (emitter channels) exist per energy; the curve uses the
worst one (largest relative residual). A right-hand side that stopped early keeps its last value.
Usage: python analysis/e3_convergence.py --set-dir RUNS/e3 --dump-dir RUNS/e3_dump --tag e3
"""
import argparse, csv, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from evaluate import metrics, LIMIT  # noqa: E402

KMAX = 6


def read_dump(path):
    """{(irhs, iter): resid}, {irhs: ||iterate0||}, final iteration count per rhs."""
    lines = Path(path).read_text().splitlines()
    resid, norm0, last = {}, {}, {}
    i = 0
    while i < len(lines):
        t = lines[i].split()
        if t[0] == 'ITERATE':
            irhs, it, j, n, r = int(t[1]), int(t[2]), int(t[3]), int(t[4]), float(t[5])
            if (irhs, it) not in resid:
                resid[(irhs, it)] = r
                if it == 0:
                    v = np.array([l.split() for l in lines[i + 1:i + 1 + n]], dtype=float)
                    norm0[irhs] = float(np.sqrt(np.sum(v ** 2)))
                last[irhs] = max(last.get(irhs, 0), it)
            i += n + 1
        elif t[0] == 'FINAL':
            i += int(t[2]) + 1
        else:
            i += 1
    return resid, norm0, last


def curve(resid, norm0, last, k):
    vals = []
    for irhs, n0 in norm0.items():
        it = min(k, last[irhs])          # stopped early: keep the last value
        vals.append(resid[(irhs, it)] / n0)
    return max(vals), float(np.median(vals))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--set-dir', required=True)
    ap.add_argument('--dump-dir', required=True)
    ap.add_argument('--tag', default='e3')
    ap.add_argument('--names', default='')
    a = ap.parse_args()
    sd, dd = Path(a.set_dir), Path(a.dump_dir)
    rows = []
    for cdir in sorted(p for p in sd.iterdir() if (p / 'lu_0' / 'result.npz').exists()):
        c = cdir.name
        dump = dd / c / 'gmres_dump/calc/output/msspec_tree_iter.dat'
        if not dump.exists():
            continue
        ref = np.load(cdir / 'lu_0/result.npz')['cross_section']
        resid, norm0, last = read_dump(dump)
        print(c, 'dump rhs', len(norm0), 'steps per rhs', sorted(set(last.values())))
        for k in range(KMAX + 1):
            f = cdir / ('k%d_0' % k) / 'result.npz'
            if not f.exists():
                continue
            m = metrics(np.load(f)['cross_section'], ref)
            rmax, rmed = curve(resid, norm0, last, k)
            rows.append(dict(case=c, step=k, resid_max=rmax, resid_med=rmed, eps_2=m['eps_2'], eps_inf=m['eps_inf']))
        f = cdir / 'gmres_0' / 'result.npz'
        if f.exists():
            m = metrics(np.load(f)['cross_section'], ref)
            kk = max(last.values())
            rmax, rmed = curve(resid, norm0, last, kk)
            rows.append(dict(case=c, step='final', resid_max=rmax, resid_med=rmed, eps_2=m['eps_2'], eps_inf=m['eps_inf']))
    out = ROOT / 'analysis/out' / (a.tag + '_curves.csv')
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print('wrote', out, len(rows), 'rows')


if __name__ == '__main__':
    main()
