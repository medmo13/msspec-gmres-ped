# coding: utf8

from msspec.calculator import MSSPEC
from ase import Atoms

# Create an atomic chain O-Rh
cluster = Atoms(['O', 'Rh'], positions = [(0,0,0), (0,0,4.)])

# Create the calculator
calc = MSSPEC(spectroscopy = 'PED')
calc.set_atoms(cluster)
cluster.absorber = 0

# compute
data = calc.get_scattering_factors(level='1s', kinetic_energy=723)

data.view()
