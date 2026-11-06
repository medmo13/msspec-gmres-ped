"""Write manuscript/numbers_results.tex: macros for the Results section, read only from analysis/out/*.csv.

Sources: p100_endpoints.csv, p100_eval.csv, p100_tol3_eval.csv, p100_tol6_eval.csv, p100_steps_eval.csv,
p100_variants.csv, e16_timing.csv, p100_lu_check.csv; stratum and rho from platform/materials_v2.csv.
Names use letters only. Suffixes: K count, N size, Pct percent, Lo/Hi interval in percent.
Usage: python analysis/make_results_numbers.py
"""
import csv
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'analysis' / 'out'
mat = {r['mat_id']: r for r in csv.DictReader(open(ROOT / 'platform/materials_v2.csv'))}
M = {}
DIG = {'0': 'Zero', '1': 'One', '2': 'Two', '3': 'Three', '4': 'Four', '5': 'Five', '6': 'Six', '7': 'Seven', '8': 'Eight', '9': 'Nine'}


def word(t):
    return ''.join(DIG.get(c, c) for c in t)
NUM = {1: 'One', 2: 'Two', 3: 'Three'}
ST = {'all': 'All', 'high': 'High', 'low': 'Low'}


def put(name, val):
    assert name not in M, name
    M[name] = val


def pct(x, d=1):
    t = ('%.' + str(d) + 'f') % (100 * float(x))
    return t[1:] if t.strip('-0.') == '' and t[0] == '-' else t


def rows(f):
    return [r for r in csv.DictReader(open(OUT / f)) if r.get('rep', '0') == '0']


# endpoints
for r in csv.DictReader(open(OUT / 'p100_endpoints.csv')):
    ep, m, st, what = r['endpoint'], r['method'], ST[r['stratum']], r['what']
    if ep == 'P1':
        b = 'ResGmres' + st
        put(b + 'PassK', r['k']); put(b + 'N', r['n']); put(b + 'PassPct', pct(r['rate']))
        put(b + 'PassLo', pct(r['lo'])); put(b + 'PassHi', pct(r['hi'])); put(b + 'Excl', r['excluded'])
    elif ep == 'P2':
        put('ResGmres' + st + 'FailK', r['k']); put('ResGmres' + st + 'FailPct', pct(r['rate']))
    elif ep == 'P3':
        b = 'ResSeries' + NUM[int(m[-1])] + st + ('Pass' if what == 'pass' else 'Div')
        put(b + 'K', r['k']); put(b + 'Pct', pct(r['rate'])); put(b + 'Lo', pct(r['lo'])); put(b + 'Hi', pct(r['hi']))
    elif ep == 'S1':
        key = 'EpsTwo' if what.endswith('eps_2') else 'EpsInf'
        d = 2 if key == 'EpsTwo' else 1
        b = 'Res' + key + 'Median' + st
        put(b, pct(r['rate'], d)); put(b + 'Lo', pct(r['lo'], d)); put(b + 'Hi', pct(r['hi'], d))
        put(b + 'N', r['n'])

# GMRES details at tolerance 1e-2
ev = rows('p100_eval.csv')
g = [r for r in ev if r['method'] == 'gmres']
fails = [r for r in g if r['status'] == 'method_failure']
okno = [r for r in g if r['status'] == 'ok' and r['passed'] != 'True']
put('ResGmresFailRhoMin', '%.2f' % min(float(mat[r['case']]['rho']) for r in fails))
put('ResGmresFailRhoMax', '%.2f' % max(float(mat[r['case']]['rho']) for r in fails))
put('ResGmresFailAtomsMax', str(max(int(mat[r['case']]['natoms']) for r in fails)))
put('ResGmresOkNotPassK', str(len(okno)))
put('ResGmresOkNotPassEpsTwoMax', pct(max(float(r['eps_2']) for r in okno)))
put('ResGmresOkNotPassEpsInfMin', pct(min(float(r['eps_inf']) for r in okno)))
put('ResGmresOkNotPassEpsInfMax', pct(max(float(r['eps_inf']) for r in okno)))
put('ResGmresPassEpsTwoMax', pct(max(float(r['eps_2']) for r in g if r['passed'] == 'True')))
put('ResLuNotRunK', str(sum(r['status'] == 'lu_not_run' for r in ev if r['method'] == 'lu')))
sp = [r for r in ev if r['method'].startswith('series') and r['passed'] == 'True']
put('ResSeriesPassCases', str(len({r['case'] for r in sp})))
put('ResSeriesPassRho', '%.2f' % float(mat[sp[0]['case']]['rho']))
put('ResSeriesPassAtoms', mat[sp[0]['case']]['natoms'])
lowser = [r for r in ev if r['method'] == 'series1' and r['status'] == 'ok' and mat[r['case']]['stratum'] == 'low']
rho1 = [r for r in lowser if float(mat[r['case']]['rho']) < 1]
put('ResSeriesOneLowRhoBelowOneN', str(len(rho1)))
put('ResSeriesOneLowRhoBelowOneDivK', str(sum(r['diverged'] == 'True' for r in rho1)))

# variants: all and per stratum
sets = {'ZeroSteps': ('p100_steps_eval.csv', 'k0'), 'OneStep': ('p100_steps_eval.csv', 'k1'),
        'TwoSteps': ('p100_steps_eval.csv', 'k2'), 'ThreeSteps': ('p100_steps_eval.csv', 'k3'),
        'TolEThree': ('p100_tol3_eval.csv', 'gmres'), 'TolESix': ('p100_tol6_eval.csv', 'gmres'),
        'TolETwo': ('p100_eval.csv', 'gmres')}
var = {r['setting']: r for r in csv.DictReader(open(OUT / 'p100_variants.csv'))}
vkey = {'ZeroSteps': '0 steps', 'OneStep': '1 step', 'TwoSteps': '2 steps', 'ThreeSteps': '3 steps',
        'TolETwo': 'tol 1e-2', 'TolEThree': 'tol 1e-3', 'TolESix': 'tol 1e-6'}
for name, (f, meth) in sets.items():
    rs = [r for r in rows(f) if r['method'] == meth and r['status'] in ('ok', 'method_failure')]
    for st in ('high', 'low'):
        s = [r for r in rs if mat[r['case']]['stratum'] == st]
        put('Var%s%sPassK' % (name, ST[st]), str(sum(r['passed'] == 'True' for r in s)))
        put('Var%s%sN' % (name, ST[st]), str(len(s)))
        put('Var%s%sFailK' % (name, ST[st]), str(sum(r['status'] == 'method_failure' for r in s)))
    v = var[vkey[name]]
    put('Var%sPassK' % name, v['passed']); put('Var%sN' % name, v['n']); put('Var%sPassPct' % name, pct(v['rate']))
    put('Var%sPassLo' % name, pct(v['lo'])); put('Var%sPassHi' % name, pct(v['hi']))
    put('Var%sFailK' % name, v['method_failures'])
    put('Var%sEpsTwoMedian' % name, '%.3g' % (100 * float(v['median_eps_2'])))
    put('Var%sWall' % name, '%.1f' % float(v['median_wall_s']))
    put('Var%sRss' % name, '%.0f' % float(v['median_rss_mib']))

# timing and memory (E16, median of 3 repeats after warm-up)
tm = list(csv.DictReader(open(OUT / 'e16_timing.csv')))
by = {(r['case'], r['method']): r for r in tm}
for c in sorted({r['case'] for r in tm}):
    cw = word(c)
    put('Tim%sAtoms' % cw, by[(c, 'lu')]['natoms'])
    for m, nm in (('lu', 'Lu'), ('gmres', 'Gmres'), ('series1', 'SeriesOne'), ('series3', 'SeriesThree')):
        put('Tim%s%sWall' % (cw, nm), '%.1f' % float(by[(c, m)]['wall_median_s']))
        put('Tim%s%sRss' % (cw, nm), '%.0f' % float(by[(c, m)]['peak_rss_mib']))
    put('Tim%sSpeedup' % cw, '%.1f' % (float(by[(c, 'lu')]['wall_median_s']) / float(by[(c, 'gmres')]['wall_median_s'])))
    put('Tim%sMemRatio' % cw, '%.1f' % (float(by[(c, 'lu')]['peak_rss_mib']) / float(by[(c, 'gmres')]['peak_rss_mib'])))
put('TimLoadMax', '%.2f' % max(max(float(r['load_before_max']), float(r['load_after_max'])) for r in tm))
put('TimLuLoadMax', '%.2f' % max(float(r['load_after_max']) for r in tm if r['method'] == 'lu'))

# LU cross-check against the saved prod3000 scans
lc = [r for r in csv.DictReader(open(OUT / 'p100_lu_check.csv')) if r['rel_l2']]
fl = [r for r in lc if r['flag'] == 'True']
put('LuCheckN', str(len(lc))); put('LuCheckFlagK', str(len(fl))); put('LuCheckOkK', str(len(lc) - len(fl)))
put('LuCheckFlagMaxPct', pct(max(float(r['rel_l2']) for r in fl)))
put('LuCheckOkMax', '%.1e' % max(float(r['rel_l2']) for r in lc if r['flag'] != 'True'))

# direct signal (no multiple scattering) of the flagged cases (plot_lu_prod3000.py output)
cz = list(csv.DictReader(open(OUT / 'lu_prod3000_causes.csv')))
fz = [r for r in cz if r['flag'] == 'True']
same = [r for r in fz if float(r['rel_direct']) == 0.0]
put('LuCheckDirectDiffK', str(sum(float(r['rel_direct']) > 0 for r in fz)))
put('LuCheckDirectSameK', str(len(same)))
put('LuCheckDirectSameMin', '%.1f\\times10^{-4}' % (1e4 * min(float(r['rel_cs']) for r in same)))
put('LuCheckDirectSameMax', '%.1f\\times10^{-4}' % (1e4 * max(float(r['rel_cs']) for r in same)))

# sensitivity to a 1e-6 Angstrom shift of all atoms (platform/work/perturb_test.py)
import numpy as np
for c in ('H17', 'L12', 'L34', 'H42', 'L43', 'L45'):
    d = ROOT / 'platform/work/perturb' / c
    if not (d / 'base.npz').exists():   # scratch data of a private check, not part of this record
        continue
    b, q = np.load(d / 'base.npz')['cross_section'], np.load(d / 'pert.npz')['cross_section']
    put('Pert%s' % word(c), '%.1e' % (np.linalg.norm(q - b) / np.linalg.norm(b)))

# the four emitters of the theta-scan figure
matv = {r['mat_id']: r for r in csv.DictReader(open(ROOT / 'platform/materials_v2.csv'))}
ev = {(r['case'], r['method']): r for r in rows('p100_eval.csv')}
for c in ('H11', 'H10', 'H08', 'L45'):
    put('Scan%sRho' % word(c), '%.2f' % float(matv[c]['rho']))
    for m, nm in (('gmres', 'Gmres'), ('series1', 'SeriesOne'), ('series2', 'SeriesTwo'), ('series3', 'SeriesThree')):
        put('Scan%s%sEpsTwo' % (word(c), nm), pct(ev[(c, m)]['eps_2']))
        put('Scan%s%sEpsInf' % (word(c), nm), pct(ev[(c, m)]['eps_inf']))

# spanning tree of the example cluster (analysis/plot_tree_cluster.py)
tc = list(csv.DictReader(open(OUT / 'tree_cluster_numbers.csv')))[0]
put('TreeAtoms', tc['natoms']); put('TreeEdges', tc['edges']); put('TreeBlock', tc['block_size'])
put('TreeLmax', tc['lmax']); put('TreeMemFull', '%.0f' % float(tc['mem_full_mib']))
put('TreeMemBlocks', '%.0f' % float(tc['mem_blocks_mib'])); put('TreeMemRatio', '%.0f' % float(tc['mem_ratio']))
put('TreeEdgeMax', '%.1f' % float(tc['edge_len_max_A'])); put('TreeEdgeMedian', '%.1f' % float(tc['edge_len_median_A']))

# E3: error against the preconditioned residual (analysis/e3_convergence.py, plot_e3.py)
e3 = list(csv.DictReader(open(OUT / 'p100_e3_curves.csv')))
e3s = {r['key']: float(r['value']) for r in csv.DictReader(open(OUT / 'p100_e3_summary.csv'))}
fin = [r for r in e3 if r['step'] == 'final']
steps = [r for r in e3 if r['step'] != '0']          # steps 1 to 6 and the converged points
res_f = [float(r['resid_max']) for r in fin]
e2_f = [float(r['eps_2']) for r in fin]
ei_f = [float(r['eps_inf']) for r in fin]
put('EThreeCases', str(len({r['case'] for r in e3})))
put('EThreeConvK', str(len(fin)))
put('EThreeResidMax', '%.1f\\times10^{-2}' % (100 * max(res_f)))
put('EThreeEpsTwoMaxPct', pct(max(e2_f)))
put('EThreeEpsTwoPassK', str(sum(x < 0.10 for x in e2_f)))
put('EThreeEpsInfPassK', str(sum(x < 0.10 for x in ei_f)))
put('EThreeEpsInfFailK', str(sum(x >= 0.10 for x in ei_f)))
put('EThreeEpsTwoBelowResidK', str(sum(e < r for e, r in zip(e2_f, res_f))))
put('EThreeRatioMedian', '%.2f' % e3s['median_ratio_final'])
put('EThreeRatioMax', '%.1f' % e3s['max_ratio_final'])
put('EThreeSpearman', '%.2f' % e3s['spearman_resid_eps2'])
small = [(float(r['resid_max']), float(r['eps_2'])) for r in steps if float(r['resid_max']) <= 1e-2]
put('EThreeSmallResidN', str(len(small)))
put('EThreeSmallResidEpsTwoMaxPct', pct(max(e for _, e in small)))
mid = [(float(r['resid_max']), float(r['eps_2'])) for r in steps if float(r['resid_max']) <= 0.1]
put('EThreeMidResidN', str(len(mid)))
put('EThreeMidResidPassPct', pct(sum(e < 0.10 for _, e in mid) / len(mid), 1))
put('EThreeMidResidEpsTwoMaxPct', pct(max(e for _, e in mid), 0))
put('EThreeStepThreeResidMedPct', pct(e3s['median_resid_step3'], 1))
put('EThreeStepThreeEpsTwoMedPct', pct(e3s['median_eps2_step3'], 1))
put('EThreeStepSixResidMedPct', pct(e3s['median_resid_step6'], 1))
put('EThreeStepSixEpsTwoMedPct', pct(e3s['median_eps2_step6'], 1))


# error of every method against dense LU, in percent (analysis/error_summary.py)
import math
def fpct(x):
    v = 100 * float(x)
    if v < 0.01:
        e = int(math.floor(math.log10(v)))
        return '\\ensuremath{%.1f\\times10^{%d}}' % (v / 10 ** e, e)
    return '%.2g' % v if v < 10 else '%.0f' % v
EM = {'gmres_e2': 'GmresETwo', 'gmres_e2c67': 'GmresETwoCap', 'gmres_e6c67': 'GmresESixCap', 'gmres_e3': 'GmresEThree', 'gmres_e6': 'GmresESix',
      'series1': 'SeriesOne', 'series2': 'SeriesTwo', 'series3': 'SeriesThree'}
for r in csv.DictReader(open(OUT / 'p100_error_summary.csv')):
    b = 'Err%s%s' % (EM[r['method']], ST[r['stratum']])
    put(b + 'N', r['n']); put(b + 'Got', r['with_result']); put(b + 'AboveOne', r['e2_above_1'])
    put(b + 'ETwoMed', fpct(r['e2_median'])); put(b + 'ETwoQOne', fpct(r['e2_q25'])); put(b + 'ETwoQThree', fpct(r['e2_q75']))
    put(b + 'ETwoMax', fpct(r['e2_max']))
    put(b + 'EInfMed', fpct(r['einf_median'])); put(b + 'EInfQOne', fpct(r['einf_q25'])); put(b + 'EInfQThree', fpct(r['einf_q75']))

# the four extra measures (analysis/extra_metrics.py): median over the clusters with a result
def f3(x):
    return '%.3g' % x if x >= 1e-2 else '\\ensuremath{%.1f\\times10^{%d}}' % (x / 10 ** math.floor(math.log10(x)), math.floor(math.log10(x)))
for r in csv.DictReader(open(OUT / 'p100_extra_summary.csv')):
    if r['label'] not in EM:
        continue
    b = 'Ext%s' % EM[r['label']]
    put(b + 'N', r['n'])
    put(b + 'ERelMed', fpct(float(r['eps_rel_pct_med']) / 100)); put(b + 'RMed', f3(float(r['R_med'])))
    put(b + 'RTwoMed', f3(float(r['R2_med']))); put(b + 'OneMinusRMed', f3(1 - float(r['pearson_med'])))
    put(b + 'PearsonMed', '%.3f' % float(r['pearson_med']))

# rerun of the clusters without a result at cap 12, with cap 67 (set p100_cap67)
import json
_rr = [json.loads(f.read_text()) for f in sorted((ROOT / 'runs/p100_cap67').glob('*/gmres_0/row.json'))]
put('CapRerunN', str(len(_rr)))
put('CapRerunCap', str(max(int(r['cap']) for r in _rr)))
put('CapRerunOk', str(sum(str(r['error']) in ('0', 'False') for r in _rr)))
put('CapRerunItersMin', str(min(int(r['iters_min']) for r in _rr)))
put('CapRerunItersMax', str(max(int(r['iters_max']) for r in _rr)))
_re = [float(r['eps_2']) for r in csv.DictReader(open(OUT / 'p100_cap67_eval.csv')) if r['method'] == 'gmres']
put('CapRerunEpsTwoMax', fpct(max(_re)))
_r6 = [json.loads(f.read_text()) for f in sorted((ROOT / 'runs/p100_tol6_cap67').glob('*/gmres_0/row.json'))]
put('CapTwoFail', str(sum(str(r['error']) not in ('0', 'False') for r in _rr)))
put('CapSixN', str(len(_r6)))
put('ErrGmresESixCapAllNMissing', str(sum(str(r['error']) not in ('0', 'False') for r in _r6)))
put('CapSixOk', str(sum(str(r['error']) in ('0', 'False') for r in _r6)))
put('CapSixFail', str(sum(str(r['error']) not in ('0', 'False') for r in _r6)))
_ok6 = [r for r in _r6 if str(r['error']) in ('0', 'False')]
put('CapSixItersMin', str(min(int(r['iters_min']) for r in _ok6)))
put('CapSixItersMax', str(max(int(r['iters_max']) for r in _ok6)))
_f6 = [r for r in _r6 if str(r['error']) not in ('0', 'False')]
put('CapSixFailRho', ', '.join('%.2f' % float(r['rho']) for r in _f6))
put('CapSixFailAtoms', ', '.join(r['natoms'] for r in _f6))

# wall time of the clusters that needed more than the first cap of 30 steps at 1e-6 (set p100_tol6_cap67),
# against dense LU and the series of the same cluster in set p100 (one run each; wall_s is the whole process)
import statistics as _stt
def _w(case, m):
    return float(json.loads((ROOT / 'runs/p100' / case / (m + '_0') / 'row.json').read_text())['wall_s'])
_tt = [(r['case'], int(r['iters_max']), int(r['natoms']), float(r['wall_s'])) for r in _ok6]
_rl = [g / _w(c, 'lu') for c, _, _, g in _tt]
_rs = [g / _w(c, 'series3') for c, _, _, g in _tt]
put('CapSixTimeLuMed', '%.2f' % _stt.median(_rl)); put('CapSixTimeLuMin', '%.2f' % min(_rl)); put('CapSixTimeLuMax', '%.2f' % max(_rl))
put('CapSixTimeSerMed', '%.1f' % _stt.median(_rs)); put('CapSixTimeSerMin', '%.1f' % min(_rs)); put('CapSixTimeSerMax', '%.1f' % max(_rs))
_wc, _wi, _wa, _wg = max(_tt, key=lambda t: t[1])
put('CapSixWorstAtoms', str(_wa)); put('CapSixWorstGmres', '%.1f' % _wg)
put('CapSixWorstLu', '%.1f' % _w(_wc, 'lu')); put('CapSixWorstSerOne', '%.1f' % _w(_wc, 'series1')); put('CapSixWorstSerThree', '%.1f' % _w(_wc, 'series3'))

# time and peak memory of dense LU against GMRES, tol 1e-2, over the clusters with a dense-LU reference
# (analysis/plot_p100_cost.py; one run each in the campaign of the test set)
_pc = list(csv.DictReader(open(OUT / 'p100_cost.csv')))
_tr = [float(r['time_ratio']) for r in _pc]
_mr = [float(r['mem_ratio']) for r in _pc]
put('CostN', str(len(_pc)))
put('CostTimeMed', '%.1f' % _stt.median(_tr)); put('CostTimeMin', '%.2f' % min(_tr)); put('CostTimeMax', '%.1f' % max(_tr))
put('CostMemMed', '%.1f' % _stt.median(_mr)); put('CostMemMin', '%.1f' % min(_mr)); put('CostMemMax', '%.1f' % max(_mr))
_small = [r for r in _pc if int(r['natoms']) < 40]
_large = [r for r in _pc if int(r['natoms']) >= 100]
put('CostSmallN', str(len(_small))); put('CostLargeN', str(len(_large)))
put('CostSmallTimeMed', '%.1f' % _stt.median(float(r['time_ratio']) for r in _small))
put('CostSmallMemMed', '%.1f' % _stt.median(float(r['mem_ratio']) for r in _small))
put('CostLargeTimeMed', '%.1f' % _stt.median(float(r['time_ratio']) for r in _large))
put('CostLargeMemMed', '%.1f' % _stt.median(float(r['mem_ratio']) for r in _large))

# clusters whose GMRES time is the rerun with the larger cap: ratio if the failed first run (runs/p100) is added
_p1 = {(r['case'], r['method']): r for r in csv.DictReader(open(ROOT / 'runs/p100/rows.csv'))}
_rr = [r for r in _pc if _p1[(r['case'], 'gmres')]['run_status'] != 'ok']
_rw = [float(r['lu_wall_s']) / (float(r['gmres_wall_s']) + float(_p1[(r['case'], 'gmres')]['wall_s'])) for r in _rr]
put('CostRerunN', str(len(_rr))); put('CostRerunWithFirstMin', '%.2f' % min(_rw)); put('CostRerunWithFirstMax', '%.2f' % max(_rw))

# Pearson correlation of the same comparisons, as 1 - r (analysis/extra_metrics.py, p100_extra_metrics.csv)
_pe = {(r['case'], r['label']): r for r in csv.DictReader(open(OUT / 'p100_extra_metrics.csv')) if r['set'] == 'p100'}
for c in ('H11', 'H10', 'H08', 'L45'):
    for lab, nm in (('gmres_e2', 'Gmres'), ('series1', 'SeriesOne'), ('series2', 'SeriesTwo'), ('series3', 'SeriesThree')):
        assert abs(float(_pe[(c, lab)]['eps_rel_pct']) - 100 * float(ev[(c, lab.replace('gmres_e2', 'gmres'))]['eps_2'])) < 1e-6, (c, lab)
        put('Scan%s%sOneMinusR' % (word(c), nm), f3(1 - float(_pe[(c, lab)]['pearson'])))

# theta scans after 1 to 4 GMRES steps and series order 1 to 3 (analysis/plot_theta_steps.py)
_cn = {'H11': 'HOneOne', 'H10': 'HOneZero', 'H08': 'HZeroEight', 'L45': 'LFourFive'}
_mn = {'GMRES, 1 step': 'GOne', 'GMRES, 2 steps': 'GTwo', 'GMRES, 3 steps': 'GThree', 'GMRES, 4 steps': 'GFour',
       'Series, order 1': 'SOne', 'Series, order 2': 'STwo', 'Series, order 3': 'SThree'}
for _r in csv.DictReader(open(OUT / 'theta_scans_steps.csv')):
    put('Step%s%sEpsTwo' % (_cn[_r['case']], _mn[_r['method']]), f3(float(_r['eps2_pct'])))
    put('Step%s%sOneMinusR' % (_cn[_r['case']], _mn[_r['method']]), f3(float(_r['one_minus_r'])))

# large clusters, series against GMRES tol 1e-3 (analysis/plot_big.py)
_bn = {'B01': 'BOne', 'B02': 'BTwo', 'B03': 'BThree', 'B04': 'BFour'}
_bm = {'GMRES, tol 1e-3': 'G', 'Series, order 1': 'SOne', 'Series, order 2': 'STwo', 'Series, order 3': 'SThree'}
for _r in csv.DictReader(open(OUT / 'big_clusters.csv')):
    put('Big%s%sEpsTwo' % (_bn[_r['case']], _bm[_r['method']]), f3(float(_r['eps2_pct'])))
    put('Big%s%sOneMinusR' % (_bn[_r['case']], _bm[_r['method']]), f3(float(_r['one_minus_r'])))
    if _r['method'] == 'GMRES, tol 1e-3':
        put('Big%sAtoms' % _bn[_r['case']], _r['natoms'])
        put('Big%sRho' % _bn[_r['case']], '%.2f' % float(_r['rho']))

# Pearson correlation with dense LU after 1 to 3 steps and series 1 to 3 (analysis/plot_steps_correlation.py)
_cr = {}
for _r in csv.DictReader(open(OUT / 'steps_correlation.csv')):
    _cr.setdefault((_r['method'], float(_r['rho']) >= 1), []).append(float(_r['pearson']))
_cm = {'G1': 'GOne', 'G2': 'GTwo', 'G3': 'GThree', 'S1': 'SOne', 'S2': 'STwo', 'S3': 'SThree'}
for (_m, _strong), _v in _cr.items():
    put('Corr%s%sMed' % (_cm[_m], 'Strong' if _strong else 'Weak'), '%.2f' % float(np.median(_v)))
_ca = {}
for (_m, _strong), _v in _cr.items():
    _ca.setdefault(_m, []).extend(_v)
for _m, _v in _ca.items():
    put('Corr%sAllMed' % _cm[_m], '%.2f' % float(np.median(_v)))
put('CorrAllN', str(len(_ca['G1'])))
put('CorrWeakN', str(len(_cr[('G1', False)])))
put('CorrStrongN', str(len(_cr[('G1', True)])))

# Stereographic maps of the four E16 emitters (analysis/out/maps_e16.csv, written by the map evaluation)
_cn = {'H11': 'HOneOne', 'H10': 'HOneZero', 'H08': 'HZeroEight', 'L45': 'LFourFive'}
_mm = {'G1': 'GOne', 'G2': 'GTwo', 'G3': 'GThree', 'Gtol2': 'GTolTwo', 'Gtol3': 'GTolThree', 'S1': 'SOne', 'S2': 'STwo', 'S3': 'SThree'}
for _r in csv.DictReader(open(OUT / 'maps_e16.csv')):
    put('Map%s%sEpsTwo' % (_cn[_r['case']], _mm[_r['method']]), f3(float(_r['eps2_pct'])))
    put('Map%s%sPearson' % (_cn[_r['case']], _mm[_r['method']]), '%.3f' % float(_r['pearson']))
    if _r['method'] == 'G1':
        put('Map%sDirs' % _cn[_r['case']], _r['n_dir'])

# largest direction-wise difference of the GMRES tol 1e-2 map from the LU map, in percent of the LU maximum
_mx = 0.0
for _c in _cn:
    _ref = np.load(ROOT / 'runs/maps_e16' / _c / 'lu_0/result.npz')['cross_section']
    _g = np.load(ROOT / 'runs/maps_e16' / _c / 'gmres_0/result.npz')['cross_section']
    _mx = max(_mx, float(np.max(np.abs(_g - _ref)) / _ref.max() * 100))
put('MapGTolTwoMaxDiff', '%.1f' % _mx)

# ten extra clusters of the Supporting Information (analysis/out/si_scans.csv, from analysis/make_si_table.py)
_si = list(csv.DictReader(open(OUT / 'si_scans.csv')))
put('SIRhoMin', '%.2f' % min(float(r['rho']) for r in _si)); put('SIRhoMax', '%.2f' % max(float(r['rho']) for r in _si))
put('SIAtomsMin', str(min(int(r['natoms']) for r in _si))); put('SIAtomsMax', str(max(int(r['natoms']) for r in _si)))
put('SIKeMin', '%.0f' % min(float(r['ke']) for r in _si)); put('SIKeMax', '%.0f' % max(float(r['ke']) for r in _si))
put('SIMaterialsN', str(len({r['name'] for r in _si}))); put('SIStructuresN', str(len({r['prototype'] for r in _si})))

# numbers that the text used to type by hand (large clusters, SI selection pool, map ranges)
import glob as _glob
_bg = [json.loads(open(f).read()) for f in sorted(_glob.glob(str(ROOT / 'runs/big_ref/B0*/gmres_0/row.json')))]
put('BigGOneMinusRMax', f3(max(float(r['one_minus_r']) for r in csv.DictReader(open(OUT / 'big_clusters.csv')) if r['method'].startswith('GMRES'))))
put('BigAtomsMin', str(min(int(r['natoms']) for r in _bg))); put('BigAtomsMax', str(max(int(r['natoms']) for r in _bg)))
_lu = list(csv.DictReader(open(ROOT / 'runs/big_lu/rows.csv')))      # measured dense LU (D39)
put('BigLuGiBMin', '%.1f' % min(float(r['peak_rss_mib']) / 1024 for r in _lu)); put('BigLuGiBMax', '%.1f' % max(float(r['peak_rss_mib']) / 1024 for r in _lu))
put('BigLuWallMin', '%.0f' % min(float(r['wall_s']) for r in _lu)); put('BigLuWallMax', '%.0f' % max(float(r['wall_s']) for r in _lu))
# the three strong scatterers that the LU memory estimate had excluded (analysis/lu_new_refs.py)
_nr = list(csv.DictReader(open(OUT / 'lu_new_refs.csv')))
_hx = [r for r in _nr if r['case'].startswith('H')]
put('LuNewAtomsMin', str(min(int(r['natoms']) for r in _hx))); put('LuNewAtomsMax', str(max(int(r['natoms']) for r in _hx)))
put('LuNewRhoMin', '%.2f' % min(float(r['rho']) for r in _hx)); put('LuNewRhoMax', '%.2f' % max(float(r['rho']) for r in _hx))
put('LuNewPeakMiB', '%.0f' % max(float(r['lu_peak_mib']) for r in _hx)); put('LuNewWallS', '%.0f' % max(float(r['lu_wall_s']) for r in _hx))
_hg = [r for r in _hx if r['method'] == 'gmres']
put('LuNewGEpsTwoMin', f3(100 * min(float(r['eps_2']) for r in _hg))); put('LuNewGEpsTwoMax', f3(100 * max(float(r['eps_2']) for r in _hg)))
put('LuNewGEpsInfMin', f3(100 * min(float(r['eps_inf']) for r in _hg))); put('LuNewGEpsInfMax', f3(100 * max(float(r['eps_inf']) for r in _hg)))
_hs = [r for r in _hx if r['method'] == 'series3']
put('LuNewSThreeEpsTwoMin', f3(100 * min(float(r['eps_2']) for r in _hs))); put('LuNewSThreeEpsTwoMax', f3(100 * max(float(r['eps_2']) for r in _hs)))
_pool = [d.name for d in (ROOT / 'runs/p100').iterdir() if d.is_dir() and all((d / (m + '_0/result.npz')).exists()
         for m in ('lu', 'gmres', 'series1', 'series2', 'series3'))]
put('SIPoolN', str(len([c for c in _pool if c not in ('H11', 'H10', 'H08', 'L45')])))
put('SIRhoWeakMax', '%.2f' % max(float(r['rho']) for r in _si if float(r['rho']) < 1))
put('SIRhoStrongMin', '%.2f' % min(float(r['rho']) for r in _si if float(r['rho']) >= 1))
_g3 = [float(_r['eps2_pct']) for _r in csv.DictReader(open(OUT / 'maps_e16.csv')) if _r['method'] == 'G3' and _r['case'] != 'L45']
put('MapStrongGThreeMin', f3(min(_g3))); put('MapStrongGThreeMax', f3(max(_g3)))

# energy scans of five clusters (analysis/plot_e17_scan.py, plot_e17_timing.py)
_en = [r for r in csv.DictReader(open(OUT / 'e17_scan.csv'))]
_et = {(r['case'], r['method']): r for r in csv.DictReader(open(OUT / 'e17_timing.csv'))}
_EK = {'H42': 'Ru', 'H17': 'PbB', 'H28': 'TiO', 'H48': 'Hf', 'B01': 'LiF'}
_er = {(r['case'], r['method']): float(r['one_minus_r']) for r in csv.DictReader(open(OUT / 'e17_scan_r.csv'))}
_eg_all, _es_all = [], []
for _c, _k in _EK.items():
    _g = [float(r['err_pct']) for r in _en if r['case'] == _c and r['method'] == 'GMRES']
    _s = [float(r['err_pct']) for r in _en if r['case'] == _c and r['method'] == 'series, order 3']
    _st = [float(r['err_pct']) for r in _en if r['case'] == _c and r['method'] == 'gmres_steps_max']
    _eg_all += _g; _es_all += _s
    put('EnScan%sAtoms' % _k, _et[(_c, 'lu')]['natoms'])
    put('EnScan%sGMed' % _k, '%.1f' % _stt.median(_g)); put('EnScan%sGMax' % _k, '%.1f' % max(_g))
    put('EnScan%sSMed' % _k, '%.0f' % _stt.median(_s)); put('EnScan%sSMin' % _k, '%.1f' % min(_s)); put('EnScan%sSMax' % _k, '%.0f' % max(_s))
    put('EnScan%sStepsMax' % _k, '%.0f' % max(_st)); put('EnScan%sStepsMin' % _k, '%.0f' % min(_st))
    put('EnScan%sGOneMinusR' % _k, f3(_er[(_c, 'GMRES')])); put('EnScan%sSOneMinusR' % _k, f3(_er[(_c, 'series, order 3')]))
    put('EnScan%sStepsTen' % _k, '%.0f' % _st[0])
    put('EnScan%sNEnergies' % _k, str(len(_g)))
    for _m, _n in (('lu', 'Lu'), ('gmres', 'Gmres'), ('series3', 'Ser')):
        put('EnScan%s%sWall' % (_k, _n), '%.0f' % float(_et[(_c, _m)]['wall_median_s']))
        put('EnScan%s%sRss' % (_k, _n), '%.0f' % float(_et[(_c, _m)]['peak_rss_mib']))
    put('EnScan%sSpeedup' % _k, '%.1f' % (float(_et[(_c, 'lu')]['wall_median_s']) / float(_et[(_c, 'gmres')]['wall_median_s'])))
    put('EnScan%sMemRatio' % _k, '%.1f' % (float(_et[(_c, 'lu')]['peak_rss_mib']) / float(_et[(_c, 'gmres')]['peak_rss_mib'])))
    put('EnScan%sGmresMin' % _k, '%.1f' % (float(_et[(_c, 'gmres')]['wall_median_s']) / 60))
    put('EnScan%sLuMin' % _k, '%.0f' % (float(_et[(_c, 'lu')]['wall_median_s']) / 60))
put('EnScanCapExceedN', str(sum(int(M['EnScan%sStepsMax' % k]) > 12 for k in _EK.values())))
put('EnScanGOneMinusRMin', f3(min(v for (c, m), v in _er.items() if m == 'GMRES'))); put('EnScanGOneMinusRMax', f3(max(v for (c, m), v in _er.items() if m == 'GMRES')))
put('EnScanSOneMinusRMin', f3(min(v for (c, m), v in _er.items() if m != 'GMRES'))); put('EnScanSOneMinusRMax', f3(max(v for (c, m), v in _er.items() if m != 'GMRES')))
put('EnScanEnergiesN', str(len([r for r in _en if r['case'] == 'H42' and r['method'] == 'GMRES'])))
put('EnScanPointsN', str(len(_eg_all)))
put('EnScanGAboveTenN', str(sum(x > 10 for x in _eg_all))); put('EnScanGMedAll', '%.1f' % _stt.median(_eg_all))
put('EnScanGBelowFiveN', str(sum(x < 5 for x in _eg_all)))
put('EnScanSAboveTenN', str(sum(x > 10 for x in _es_all)))
put('EnScanSBelowTenN', str(sum(x <= 10 for x in _es_all)))

# series of order 4 on a subset (analysis/series4_subset.py, D43)
_s4 = list(csv.DictReader(open(OUT / 'series4_subset.csv')))
put('SFourN', str(len(_s4)))
put('SFourRhoMin', '%.2f' % min(float(r['rho']) for r in _s4)); put('SFourRhoMax', '%.2f' % max(float(r['rho']) for r in _s4))
put('SFourAtomsMin', str(min(int(r['natoms']) for r in _s4))); put('SFourAtomsMax', str(max(int(r['natoms']) for r in _s4)))
put('SFourEpsTwoMin', '%.0f' % min(float(r['eps2_s4']) for r in _s4)); put('SFourEpsTwoMax', '%.0f' % max(float(r['eps2_s4']) for r in _s4))
put('SFourBelowTenN', str(sum(float(r['eps2_s4']) <= 10 for r in _s4)))
put('SFourGEpsTwoMax', '%.1f' % max(float(r['eps2_g']) for r in _s4))
_rt = [float(r['wall_s4']) / float(r['wall_s3']) for r in _s4]
put('SFourTimeRatioMin', '%.0f' % min(_rt)); put('SFourTimeRatioMax', '%.0f' % max(_rt))
put('SFourWallMin', '%.0f' % min(float(r['wall_s4']) for r in _s4)); put('SFourWallMax', '%.0f' % max(float(r['wall_s4']) for r in _s4))
put('SThreeWallMax', '%.0f' % max(float(r['wall_s3']) for r in _s4))
put('SFourGWallMax', '%.0f' % max(float(r['wall_g']) for r in _s4))
_b3 = sum(float(r['eps2_s4']) < float(r['eps2_s3']) for r in _s4); put('SFourBetterThanThreeN', str(_b3))
_s5 = list(csv.DictReader(open(OUT / 'series5_stopped.csv')))[0]
put('SFiveStopAtoms', _s5['natoms']); put('SFiveStopMin', _s5['elapsed_min_when_stopped'])

# ablation: GMRES with and without the tree preconditioner (analysis/ablation_nopre.py, D40 to D42)
_ab = list(csv.DictReader(open(OUT / 'ablation_nopre.csv')))
for _tl, _k in (('1e-2', 'Two'), ('1e-6', 'Six')):
    _a = [r for r in _ab if r['tol'] == _tl and r['tree_ok'] == '1' and r['nopre_ok'] == '1']
    _q = [int(r['nopre_steps']) / int(r['tree_steps']) for r in _a]
    put('Abl%sN' % _k, str(len(_a)))
    put('Abl%sRatioMed' % _k, '%.2f' % _stt.median(_q)); put('Abl%sRatioMin' % _k, '%.2f' % min(_q)); put('Abl%sRatioMax' % _k, '%.2f' % max(_q))
    put('Abl%sStepsTreeMed' % _k, '%.0f' % _stt.median(int(r['tree_steps']) for r in _a)); put('Abl%sStepsNoMed' % _k, '%.0f' % _stt.median(int(r['nopre_steps']) for r in _a))
    put('Abl%sStepsTreeMax' % _k, str(max(int(r['tree_steps']) for r in _a))); put('Abl%sStepsNoMax' % _k, str(max(int(r['nopre_steps']) for r in _a)))
    put('Abl%sEpsTreeMed' % _k, f3(100 * _stt.median(float(r['tree_eps2']) for r in _a))); put('Abl%sEpsNoMed' % _k, f3(100 * _stt.median(float(r['nopre_eps2']) for r in _a)))
    put('Abl%sNoLowerN' % _k, str(sum(int(r['nopre_steps']) <= int(r['tree_steps']) for r in _a)))
_h = [r for r in _ab if r['tol'] == '1e-6' and r['case'] == 'H24'][0]
put('AblHBigAtoms', _h['natoms']); put('AblHBigRho', '%.2f' % float(_h['rho']))
put('AblHBigStepsTree', _h['tree_steps']); put('AblHBigStepsNo', _h['nopre_steps'])
put('AblHBigEpsNo', f3(100 * float(_h['nopre_eps2']))); put('AblHBigEpsTree', f3(100 * float(_h['tree_eps2'])))

# rerun with a larger cap of the two H24 runs that did not converge at the cap of the main tables (sets p100_tol3_cap100, p100_tol6_cap100)
import sys as _sys; _sys.path.insert(0, str(ROOT / "analysis"))
from evaluate import metrics as _metrics
_lu24 = np.load(ROOT / 'runs/p100/H24/lu_0/result.npz')['cross_section']
for _s, _k in [('p100_tol3_cap100', 'CapHThree'), ('p100_tol6_cap100', 'CapHSix')]:
    _r = json.loads((ROOT / 'runs' / _s / 'H24/gmres_0/row.json').read_text())
    assert str(_r['error']) in ('0', 'False')
    _m = _metrics(np.load(ROOT / 'runs' / _s / 'H24/gmres_0/result.npz')['cross_section'], _lu24)
    put(_k + 'Cap', str(int(_r['cap']))); put(_k + 'ItersMin', str(int(_r['iters_min']))); put(_k + 'ItersMax', str(int(_r['iters_max'])))
    put(_k + 'EpsTwo', fpct(_m['eps_2']))
put('CapHAtoms', str(int(_r['natoms']))); put('CapHRho', '%.2f' % float(_r['rho']))

lines = ['% Generated by analysis/make_results_numbers.py. Do not edit.']
lines += ['\\newcommand{\\%s}{%s}' % (k, v) for k, v in M.items()]
(ROOT / 'manuscript/numbers_results.tex').write_text('\n'.join(lines) + '\n')
print('wrote manuscript/numbers_results.tex', len(M), 'macros')
