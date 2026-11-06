# coding: utf8

from ase.build import bulk
from msspec.utils import *

a = 4.21                                                           # The lattice parameter
epsilon = 0.01                                                     # a small usefull value
MgO = bulk('MgO', a=a, crystalstructure='rocksalt', cubic=True)    # ASE magic here

# shaping the cluster in half sphere
MgO = MgO.repeat((10, 10, 10))
center_cluster(MgO)
MgO = cut_plane(MgO, z=epsilon)
MgO = cut_sphere(MgO, radius=3.*a + epsilon)

print len(MgO)
raise

from ase.io import write
MgO.rotate(-65, 'x')
for a in range(0, 360, 5):
	write('img_{:03d}.pov'.format(a), MgO, camera_type='perspective', display=False, rotation='{:d}y'.format(a), 
	      transparent=True)

