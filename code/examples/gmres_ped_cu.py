"""Photoelectron diffraction of a Cu cluster with the tree-preconditioned GMRES solver.

Run from a work folder:  PYTHONPATH=<repo>/code/msspec_gmres/src python <repo>/code/examples/gmres_ped_cu.py (writes calc/Cu there)
The same script with algorithm='inversion' gives the dense LU reference.
"""

import numpy as np
from ase.build import bulk
from msspec.calculator import MSSPEC
from msspec.utils import center_cluster, cut_sphere, cut_plane, get_atom_index


# Cu cluster of 4 planes with the emitter 3.61 A below the surface
a0 = 3.61
cluster = bulk('Cu', a=a0, cubic=True).repeat((2, 2, 2))
center_cluster(cluster)
cluster = cut_sphere(cluster, radius=a0 + 0.01)
cluster = cut_plane(cluster, z=0.1)
cluster.absorber = get_atom_index(cluster, 0, 0, -a0)

# The solver is chosen by the algorithm; its settings are calculation parameters
calc = MSSPEC(spectroscopy='PED', algorithm='gmres', folder='calc/Cu', txt=None)
params = calc.calculation_parameters

params.linear_solver_tolerance = 1e-2       # stop when the residual is below this
params.linear_solver_max_iterations = 67    # give up after this many steps
params.linear_solver_restart = 6            # steps before GMRES restarts

# Overrides the tolerance: stop after exactly 2 steps (0 = preconditioner only, -1 = off)
params.linear_solver_fixed_steps = 2

calc.set_atoms(cluster)
data = calc.get_theta_scan(level='2p3/2', theta=np.linspace(-70, 70, 15),
                           phi=0, kinetic_energy=300.0)
print(np.asarray(data[-1].cross_section))
