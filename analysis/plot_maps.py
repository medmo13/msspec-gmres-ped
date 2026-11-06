"""Stereographic maps of dense LU, GMRES after 1 to 3 steps, converged GMRES and the series of order 1 to 3.

Usage: python analysis/plot_maps.py <set dir> <case> <out stem> [gmres tolerance label]
Reads <set dir>/<case>/{lu_0,k1_0,k2_0,k3_0,gmres_0,series1_0..series3_0}/result.npz (cross_section, theta, phi).
Style after msspec_gpu/archive/v1_msspec_gpu_inversion/harness/plots_physics.py (draw_cluster, stereo_panel): ASE cluster
views with the emitter in magenta; stereographic projection r = tan(theta/2), rings every 15 degrees, grey intensity, RdBu_r
differences; the maps are interpolated for display only, eps_2 uses the native points (analysis/evaluate.metrics).
Row 1: intensity normalised to the LU maximum. Row 2: difference to dense LU in percent of the LU maximum, saturated at
+-10 percent. Writes figures/<out stem>.png/.pdf.
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from evaluate import metrics  # noqa: E402
from cluster_ase import draw  # noqa: E402

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
                     'font.size': 7, 'axes.linewidth': 0.8, 'mathtext.default': 'regular'})


def stereo_panel(ax, theta_deg, phi_deg, vals, vmin, vmax, cmap, levels=None):
    """Stereographic projection r = tan(theta/2); scattered samples are triangulated and resampled on a pixel grid
    for display only (copied from plots_physics.stereo_panel, v1 msspec_gpu)."""
    th, ph = np.radians(theta_deg), np.radians(phi_deg)
    r = np.tan(th / 2)
    xy = np.column_stack([r * np.cos(ph), r * np.sin(ph)])
    _, idx = np.unique(np.round(xy, 10), axis=0, return_index=True)
    tri = mtri.Triangulation(xy[idx, 0], xy[idx, 1])
    rmax = np.tan(np.radians(theta_deg.max()) / 2)
    xx, yy = np.meshgrid(np.linspace(-rmax, rmax, 220), np.linspace(-rmax, rmax, 220))
    rr = np.hypot(xx, yy)
    z = mtri.LinearTriInterpolator(tri, vals[idx])(xx, yy)
    z = np.ma.masked_where((rr > rmax) | (rr < np.tan(np.radians(2) / 2)), z)
    art = ax.pcolormesh(xx, yy, z, cmap=cmap, vmin=vmin, vmax=vmax, shading='auto', rasterized=True)
    ax.set_aspect('equal')
    lim = 1.12 * rmax
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.axis('off')
    for t in (30, 60, 90):
        rt = np.tan(np.radians(t) / 2)
        ax.add_patch(plt.Circle((0, 0), rt, fill=False, ec='#767676', alpha=.4, lw=.4))
    for p in range(0, 360, 90):
        a = np.radians(p)
        ax.plot([0, rmax * np.cos(a)], [0, rmax * np.sin(a)], color='#767676', alpha=.25, lw=.4)
    return art


def main(set_dir, case, stem, tol_label='tol 10$^{-2}$'):
    names = [('lu_0', 'Dense LU'), ('k1_0', 'GMRES\n1 step'), ('k2_0', 'GMRES\n2 steps'), ('k3_0', 'GMRES\n3 steps'),
             ('gmres_0', 'GMRES\n' + tol_label), ('series1_0', 'Series\norder 1'), ('series2_0', 'Series\norder 2'),
             ('series3_0', 'Series\norder 3')]
    z = {n: np.load(Path(set_dir) / case / n / 'result.npz') for n, _ in names}
    ref = z['lu_0']['cross_section']
    th, ph = z['lu_0']['theta'], z['lu_0']['phi']
    fig = plt.figure(figsize=(7.2, 4.9))
    gs = fig.add_gridspec(4, 8, height_ratios=[1.3, 1, 1, 0.07], hspace=0.22, wspace=0.05,
                          left=0.06, right=0.99, top=0.95, bottom=0.07)
    axs = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[0, 4:6])]
    mat, esym = draw(axs, case)
    axi = fig.add_subplot(gs[0, 6:8]); axi.axis('off')
    axi.text(0, 0.85, '%s (%s)' % (mat['name'], mat['prototype']), fontsize=8, fontweight='bold', va='top')
    axi.text(0, 0.6, '%s atoms\n%.0f eV\nemitter %s %s\n(magenta)' % (mat['natoms'], float(mat['ke']), esym, mat['level']),
             fontsize=7, va='top', linespacing=1.5)
    fig.text(0.01, 0.975, '(a)', fontsize=9, fontweight='bold', va='top')
    fig.text(0.01, 0.60, '(b)', fontsize=9, fontweight='bold', va='top')
    for j, (n, title) in enumerate(names):
        a = z[n]['cross_section']
        ax = fig.add_subplot(gs[1, j])
        ci = stereo_panel(ax, th, ph, a / ref.max(), 0, 1, 'gray')
        e = 'reference' if n == 'lu_0' else (r'$\varepsilon_2$ = %.3g %%' % (100 * metrics(a, ref)['eps_2'])
                                             + '\n$1-r$ = %.2g' % (1 - np.corrcoef(a, ref)[0, 1]))
        ax.set_title(title + '\n' + e, fontsize=6.5, pad=3, linespacing=1.3)
        if n != 'lu_0':
            axd = fig.add_subplot(gs[2, j])
            cd = stereo_panel(axd, th, ph, 100 * (a - ref) / ref.max(), -10, 10, 'RdBu_r')
    cax = fig.add_subplot(gs[3, 0:3]); cb = fig.colorbar(ci, cax=cax, orientation='horizontal', ticks=[0, 0.5, 1])
    cb.set_label('intensity / LU maximum', fontsize=6.5, labelpad=1); cb.ax.tick_params(labelsize=6, length=2)
    cax2 = fig.add_subplot(gs[3, 4:8]); cb2 = fig.colorbar(cd, cax=cax2, orientation='horizontal', ticks=[-10, -5, 0, 5, 10])
    cb2.set_label('difference to LU (% of LU maximum, saturated at $\\pm$10)', fontsize=6.5, labelpad=1)
    cb2.ax.tick_params(labelsize=6, length=2)
    fig.text(0.035, 0.38, 'Intensity', rotation=90, va='center', fontsize=7)
    fig.text(0.035, 0.20, 'Difference', rotation=90, va='center', fontsize=7)
    out = ROOT / 'figures' / stem
    fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
    fig.savefig(str(out) + '.pdf', bbox_inches='tight')
    print('wrote', out)


if __name__ == '__main__':
    main(*sys.argv[1:])
