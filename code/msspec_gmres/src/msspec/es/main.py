# coding: utf-8

from ase import Atoms
from ase.io import write,read
from ase.visualize import view
from scipy.spatial import ConvexHull
from es_mod import empty_spheres as esph


"""=============Generate empty spheres in copper cluster
Structure = read('cluster_examples/copper.xyz')
struct = np.ndarray.tolist(Structure.positions)
set = esph.Delaunay_Intersphere(struct)
Set=Atoms(positions=set)
view(Structure+Set)
view(Set)
#"""#====================================================

from msspec.calculator import MSSPEC
from msspec.utils import *

#"""=============Use Python MsSpec
cluster = read('cluster_examples/GeCl4.xyz')

# Set the absorber (the deepest atom centered in the xy-plane)
cluster.absorber = 0
# Create a calculator for the PhotoElectron Diffration
calc = MSSPEC(spectroscopy='PED')
# Set the cluster to use for the calculation
calc.set_atoms(cluster)

# Run the calculation
data = calc.get_theta_scan(level='2p3/2', kinetic_energy=[320,325,5])

# Show the results
data.view()
#"""#===============================

#"""===================MsSpec on ClusterC Test :
cluster = read('ClusterFinal.xyz')
# Set the absorber (the deepest atom centered in the xy-plane)
cluster.absorber = 0
# Create a calculator for the PhotoElectron Diffration
calc = MSSPEC(spectroscopy='PED')
# Set the cluster to use for the calculation
calc.set_atoms(cluster)

# Run the calculation
data = calc.get_theta_scan(level='2p3/2', kinetic_energy=[320,325,5])

# Show the results
data.view()
#"""#===============================