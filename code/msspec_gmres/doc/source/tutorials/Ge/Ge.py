# -*- encoding: utf-8 -*-
# vim: set fdm=indent ts=4 sw=4 sts=4 et ai tw=80 cc=+0 mouse=a nu : #

from msspec.calculators import PED, XRaySource
from msspec.utils import *

from ase.build  import bulk
from ase.build  import fcc111, add_adsorbate
from ase.visualize import view
from ase.data import reference_states, atomic_numbers

from matplotlib import pyplot as plt
import numpy as np
import sys


n  = 1
symbol = 'Ge'
cluster = bulk(symbol, cubic = True)
cluster = cluster.repeat((10, 10, 10))
center_cluster(cluster)
a0 = reference_states[atomic_numbers[symbol]]['a']
cluster = cut_sphere(cluster, radius = a0 * (n+ .01))
cluster = cut_plane(cluster, z = 0.1)
view(cluster)
#sys.exit(0)
#cluster.absorber = 12
#cluster.absorber = get_atom_index(cluster, 0,0,-1*n*a0)


calc = PED(folder = './Ge')
calc.set_atoms(cluster)

calc.set_calculation_parameters(scattering_order = 3,RA_cut_off = 1, 
        mean_square_vibration = 0.0)
#calc.set_source_parameters(theta = -55., phi = 0.)
    
calc.set_source_parameters(energy = 100.)
calc.set_muffintin_parameters(interstitial_potential = 0)


if False:
    level = '3d'
    all_be    = np.linspace(27., 32., 5)
    all_theta = 76.
    data = []

    for be in all_be:
        theta, Xsec, dirsig = calc.get_theta_scan(theta = 76., level = level, 
                binding_energy = be)
        print '<'*80
        print Xsec
        data.append(Xsec - dirsig)

    data = np.array(data).flatten()
    print data
    ax = plt.subplot(111)
    ax.plot(all_be, data)
    plt.show()

if True:
    ax = plt.subplot(111)
    for V in [0., 15.]:
        level = '3d'
        all_theta = np.arange(0., 80., 1.)
        phi = 0.

        calc.set_muffintin_parameters(interstitial_potential = V)

        cluster.absorber = 18
        surface_data = calc.get_theta_scan(theta = all_theta, level = level, 
                    phi = phi)#, plane = -1)

        cluster.absorber = 12
        bulk_data = calc.get_theta_scan(theta = all_theta, level = level, 
                    phi = phi)#, plane = -1)

        theta0, Xsec0, dirsig0 = surface_data
        theta1, Xsec1, dirsig1 = bulk_data
        #plt.plot(theta0, Xsec0, label = 'surface')
        #plt.plot(theta1, Xsec1, label = 'bulk')
        plt.plot(theta1, Xsec0/Xsec1, label = 'Vint = {:.2f} eV'.format(V))

    plt.ylim(None, 8)
    plt.legend()
    plt.show()


if False:
    # compute
    level = '2s'
    ke    = 156.5

    all_theta = np.arange(0.,91.)
    theta, phi, Xsec, dirsig = calc.get_stereo_scan(theta = all_theta, level = level,
                                        kinetic_energy = ke)



    # plot
    X, Y = np.meshgrid(np.radians(phi), theta)
    ax = plt.subplot(111, projection = 'polar')
    ax.pcolormesh(X, Y, Xsec, cmap='gray')
    #plt.title(r'z_0 = %.3f $\AA$' % z0)

    plt.show()

