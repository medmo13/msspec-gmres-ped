"""P100 GMRES variants: steps 0 to 3, tolerance 1e-2, 1e-3, 1e-6 (A5).

Reads analysis/out/{p100,p100_tol3,p100_tol6,p100_steps}_eval.csv (evaluate.py) and
runs/<set>/rows.csv (wall time, peak memory). Settings with no LU reference are left out.
Writes analysis/out/p100_variants.csv, analysis/out/p100_variants.npz, figures/p100_variants.png/.pdf.
"""
import csv
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from stats import wilson  # noqa: E402

# (label, set, method, colour)
SETTINGS = [('0 steps', 'p100_steps', 'k0', '#CFCECE'), ('1 step', 'p100_steps', 'k1', '#E9A6A1'),
            ('2 steps', 'p100_steps', 'k2', '#B64342'), ('3 steps', 'p100_steps', 'k3', '#8B2A2A'),
            ('tol 1e-2', 'p100', 'gmres', '#3775BA'), ('tol 1e-3', 'p100_tol3', 'gmres', '#0F4D92'),
            ('tol 1e-6', 'p100_tol6', 'gmres', '#06294F')]
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 10, 'axes.linewidth': 1.2,
                     'axes.spines.top': False, 'axes.spines.right': False})


def main():
    out, arrs = [], {}
    for label, s, m, col in SETTINGS:
        ev = [r for r in csv.DictReader(open(ROOT / 'analysis/out' / ('%s_eval.csv' % s)))
              if r['method'] == m and r['status'] in ('ok', 'method_failure')]
        rows = {r['case']: r for r in csv.DictReader(open(ROOT / 'runs' / s / 'rows.csv')) if r['method'] == m}
        n = len(ev)
        k = sum(1 for r in ev if r['passed'] == 'True')
        fail = sum(1 for r in ev if r['status'] == 'method_failure')
        e2 = np.array([float(r['eps_2']) for r in ev if r['status'] == 'ok'])
        extra = {'p100': 'p100_cap67', 'p100_tol6': 'p100_tol6_cap67'}.get(s)
        if extra and m == 'gmres':  # add the clusters that converged in the cap-67 rerun (error only)
            have = {r['case'] for r in ev if r['status'] == 'ok'}
            e2 = np.append(e2, [float(r['eps_2']) for r in csv.DictReader(open(ROOT / 'analysis/out' / ('%s_eval.csv' % extra)))
                                if r['method'] == m and r['status'] == 'ok' and r['case'] not in have])
        om = []                                   # 1 - Pearson r against dense LU, from the saved scans
        for r in ev:
            if r['status'] != 'ok':
                continue
            d = ROOT / 'runs' / s / r['case']
            om.append(1 - np.corrcoef(np.load(d / (m + '_0') / 'result.npz')['cross_section'],
                                      np.load(d / 'lu_0' / 'result.npz')['cross_section'])[0, 1])
        if extra and m == 'gmres':
            for r in csv.DictReader(open(ROOT / 'analysis/out' / ('%s_eval.csv' % extra))):
                if r['method'] == m and r['status'] == 'ok' and r['case'] not in have:
                    d = ROOT / 'runs' / extra / r['case']
                    om.append(1 - np.corrcoef(np.load(d / 'gmres_0' / 'result.npz')['cross_section'],
                                              np.load(d.parent.parent / s / r['case'] / 'lu_0' / 'result.npz')['cross_section'])[0, 1])
        om = np.array(om)
        wall = np.array([float(rows[r['case']]['wall_s']) for r in ev])
        rss = np.array([float(rows[r['case']]['peak_rss_mib']) for r in ev])
        rate, lo, hi = wilson(k, n)
        out.append(dict(setting=label, n=n, passed=k, method_failures=fail, rate=rate, lo=lo, hi=hi,
                        median_eps_2=float(np.median(e2)), median_one_minus_r=float(np.median(om)), median_wall_s=float(np.median(wall)),
                        median_rss_mib=float(np.median(rss)), max_rss_mib=float(rss.max())))
        key = label.replace(' ', '_')
        arrs[key + '_eps_2'] = e2; arrs[key + '_one_minus_r'] = om; arrs[key + '_wall_s'] = wall; arrs[key + '_rss_mib'] = rss
    with open(ROOT / 'analysis/out/p100_variants.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
    np.savez(ROOT / 'analysis/out/p100_variants.npz', **arrs)

    fig, (b, c, d) = plt.subplots(1, 3, figsize=(13, 4.2))
    x = np.arange(len(SETTINGS))
    for i, (o, (label, _, _, col)) in enumerate(zip(out, SETTINGS)):
        d.bar(i, o['median_rss_mib'], color=col, width=0.7)
    data = [arrs[l.replace(' ', '_') + '_eps_2'] for l, _, _, _ in SETTINGS]
    bp = b.boxplot(data, positions=x, widths=0.6, patch_artist=True, showfliers=True,
                   flierprops={'markersize': 2.5}, medianprops={'color': 'black'})
    for p, (_, _, _, col) in zip(bp['boxes'], SETTINGS):
        p.set_facecolor(col)
    b.set_yscale('log')
    bq = c.boxplot([np.maximum(arrs[l.replace(' ', '_') + '_one_minus_r'], 1e-12) for l, _, _, _ in SETTINGS], positions=x, widths=0.6,
                   patch_artist=True, showfliers=True, flierprops={'markersize': 2.5}, medianprops={'color': 'black'})
    for p_, (_, _, _, col) in zip(bq['boxes'], SETTINGS):
        p_.set_facecolor(col)
    c.set_yscale('log'); c.set_ylabel(r'$1-r$ (Pearson), against dense LU'); c.set_title('(b) Pearson correlation per case', loc='left')
    b.set_ylabel(r'Relative error $\varepsilon_2$ against dense LU'); b.set_title('(a) L2 error per case, against dense LU', loc='left')
    d.set_ylabel('Median peak memory (MiB)'); d.set_title('(c) Peak memory', loc='left')
    for ax in (b, c, d):
        ax.set_xticks(x); ax.set_xticklabels([s[0] for s in SETTINGS], rotation=25)
    fig.tight_layout()
    fig.savefig(ROOT / 'figures/p100_variants.png', dpi=200); fig.savefig(ROOT / 'figures/p100_variants.pdf')
    print('wrote p100_variants.csv/.npz and figures/p100_variants.png/.pdf')


if __name__ == '__main__':
    main()
