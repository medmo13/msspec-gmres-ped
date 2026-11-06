# coding: utf8

from msspec.calculator import MSSPEC
from msspec.utils import hemispherical_cluster, get_atom_index

from ase.build  import bulk
from ase.visualize import view

a0 = 3.6 # The lattice parameter in angstroms

# Create the copper cubic cell
copper = bulk('Cu', a=a0, cubic=True)
cluster = hemispherical_cluster(copper, planes=3, emitter_plane=2)

# Set the absorber (the deepest atom centered in the xy-plane)
cluster.absorber = get_atom_index(cluster, 0, 0, 0)

# Create a calculator for the PhotoElectron Diffration
calc = MSSPEC(spectroscopy='PED')

# Set the cluster to use for the calculation
calc.set_atoms(cluster)

# Run the calculation
data = calc.get_theta_scan(level='2p3/2')

# Show the results
data.view()

# Clean temporary files
calc.shutdown()
