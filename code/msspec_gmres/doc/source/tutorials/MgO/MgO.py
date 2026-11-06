# coding: utf-8

from msspec.calculators import PED
from msspec.utils import *

from ase.build  import bulk
from ase.visualize import view

from matplotlib import pyplot as plt
import sys

a0 = 4.21 # The lattice parameter in angstroms
n  = 2    # An integer useful to tweak the size of
          # the cluster

#all_n = (1., 1.5, 2.)
all_n = (1., 1.5)
level = '2p3/2'

if False:
    for n_idx, n in enumerate(all_n):
        # Create the copper cubic cell
        cluster = bulk('MgO', crystalstructure = 'rocksalt', a = a0, cubic = True)
        # Repeat the cell many times along x, y and z
        x = int(4*n)
        cluster = cluster.repeat((x, x, x))
        # Put the center of the structure at the origin
        center_cluster(cluster)
        # Keep atoms inside a given radius
        cluster = cut_sphere(cluster, radius = n*a0 + .01)
        # Keep only atoms below the plane z <= 0
        cluster = cut_plane(cluster, z = 0.1)

        # Set the absorber (for example the deeper atom in z,
        # centered in the x-y plane)
        cluster.absorber = get_atom_index(cluster, 0,0,-a0)
        #view(cluster)
        #continue
        # Create a calculator for the PhotoElectron Diffration
        calc = PED(folder = 'calc' + str(n_idx))
        # Set the cluster to use for the calculation
        calc.set_atoms(cluster)

        # Run the calculation
        calc.get_theta_scan(level = level)

        calc.save('data%s.msp' % str(n_idx))


#import sys
#sys.exit()

calc = PED()
ax = plt.subplot(111)

for i, n in enumerate(all_n):
    calc.load('data{:d}.msp'.format(i))
    x, y, _ = calc.get_results()
    #y = y/max(y)
    ax.plot(x, y, label = '{:d} atoms (n = {:.1f})'.format(len(calc.get_atoms()), n) )

ax.legend()
ax.grid(True)
ax.set_xlabel(r'$\theta$ ($\degree$)')
ax.set_ylabel(r'Normalized intensity')
ax.set_title(r'Polar scan of Mg({}) @ {:.2f} eV'.format(level, calc.parameters['kinetic_energy'][0]))
plt.show()

