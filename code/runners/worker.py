"""One PED calculation in a fresh process (called by run_set.py, never by hand).

Builds the calculator the way the prod3000 pipeline did (lmax 8 imposed, Debye model,
'debye_waller' damping, temperature and Debye temperature from the case), runs one scan
and saves result.npz. A Fortran STOP ends this process; the runner reads the return code.

argv[1] is a JSON object: src, platform, case (one row of materials_v2.csv), method,
scan, folder, out, tol, cap, restart.
"""
import json
import resource
import sys
import time

import numpy as np

# The Fortran expansion code (correlated Debye model) overflows the default 8 MB stack.
_soft, _hard = resource.getrlimit(resource.RLIMIT_STACK)
if _soft != _hard:
    resource.setrlimit(resource.RLIMIT_STACK, (_hard, _hard))

cfg = json.loads(sys.argv[1])
sys.path.insert(0, cfg['src'])
sys.path.insert(0, cfg['platform'])
from msspec.calculator import MSSPEC          # noqa: E402
from rho_probe import read_cluster            # noqa: E402

case, method, scan = cfg['case'], cfg['method'], cfg['scan']
atoms = read_cluster(cfg['run_dir'])
assert len(atoms) == int(case['natoms']) and atoms.absorber == int(case['absorber_idx'])

t0 = time.time()
if method == 'lu':
    calc = MSSPEC(spectroscopy='PED', algorithm='inversion', folder=cfg['folder'], txt=None)
elif method in ('gmres', 'nopre') or (method[0] == 'k' and method[1:].isdigit()):
    calc = MSSPEC(spectroscopy='PED', algorithm='gmres', folder=cfg['folder'], txt=None)
    calc.calculation_parameters.set_parameter('linear_solver_tolerance', cfg['tol'])
    calc.calculation_parameters.set_parameter('linear_solver_max_iterations', cfg['cap'])
    calc.calculation_parameters.set_parameter('linear_solver_restart', cfg['restart'])
elif method.startswith('series'):
    calc = MSSPEC(spectroscopy='PED', algorithm='expansion', folder=cfg['folder'], txt=None)
    calc.calculation_parameters.scattering_order = int(method[6:])
    calc.spec_parameters.calc_npath = 0
else:
    raise SystemExit('unknown method ' + method)
calc.tmatrix_parameters.lmax_mode = 'imposed'
calc.tmatrix_parameters.lmaxt = int(case['lmax'])
calc.calculation_parameters.use_debye_model = case['use_debye_model'] == 'True'
calc.calculation_parameters.vibrational_damping = 'debye_waller'
calc.calculation_parameters.temperature = float(case['temperature'])
calc.calculation_parameters.debye_temperature = float(case['debye_temperature'])
calc.set_atoms(atoms)

ke, level, phi = float(case['ke']), case['level'], float(case['phi'])
if scan['type'] == 'energy':
    # One energy scan: all energies in one calculation (E17), at the case's phi and the listed polar angles.
    data = calc.get_energy_scan(phi=phi, theta=scan['thetas'] if len(scan['thetas']) > 1 else scan['thetas'][0],
                                level=level, kinetic_energy=tuple(scan['energy']))
elif scan['type'] == 'theta':
    data = calc.get_theta_scan(phi=phi, theta=np.linspace(*scan['theta']), level=level,
                               kinetic_energy=ke)
else:
    # Maps need room for the stereographic angle grid; the default NPH_M (8000) is too
    # small. Same run-time override as the stereo examples in msspec_gpu (D29).
    data = calc.get_theta_phi_scan(phi=np.linspace(*scan['phi']), theta=np.linspace(*scan['theta']),
                                   level=level, kinetic_energy=ke, malloc={'NPH_M': 200000, 'NDIM_M': 400000})
ds = data.get_last_dset()
np.savez(cfg['out'], cross_section=np.asarray(ds.cross_section, dtype=float),
         direct_signal=np.asarray(ds.direct_signal, dtype=float),
         theta=np.asarray(ds.theta, dtype=float), phi=np.asarray(ds.phi, dtype=float),
         energy=np.asarray(ds.energy, dtype=float))
print('WORKER_CALC_S %.3f' % (time.time() - t0))
