# coding: utf8

from msspec.calculator import MSSPEC
from ase.build import fcc111, add_adsorbate
import numpy as np

data = None
all_z = np.arange(1.10, 1.65, 0.05)
for zi, z0 in enumerate(all_z):
    # construct the cluster
    cluster = fcc111('Rh', size = (2,2,1))
    cluster.pop(3)
    add_adsorbate(cluster, 'O', z0, position = 'fcc')
    cluster.absorber = len(cluster) - 1

    # Define a calculator
    calc = MSSPEC(spectroscopy='PED', algorithm='inversion')
    calc.set_atoms(cluster)

    # Compute
    data = calc.get_theta_phi_scan(level='1s', kinetic_energy=723, data=data,
                                   malloc={'NPH_M': 8000})
    dset = data[-1]
    dset.title = "{:d}) z = {:.2f} angstroms".format(zi, z0)

data.view()
