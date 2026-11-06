#!/usr/bin/env python
# coding: utf-8
#
# Copyright © 2016-2020 - Rennes Physics Institute
#
# This file is part of msspec.
#
# msspec is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# msspec is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this msspec.  If not, see <http://www.gnu.org/licenses/>.
#
# Source file  : src/msspec/utils.py
# Last modified: Wed, 26 Feb 2025 11:15:03 +0100
# Committed by : Sylvain Tricot <sylvain.tricot@univ-rennes.fr>


"""
Module utils
============


"""


import os
import re

import numpy as np
import ase.atom
from ase import Atom
from ase import Atoms

try:
    from msspec.iodata import Data
except ImportError as err:
    print(err)

from msspec.misc import LOGGER


class ForeignPotential(object):
    def __init__(self):
        self.data = Data(title='Foreign Potential')
        # Load exported potentials
        # phagen_data is a dictionary with:
        # self.phagen_data = {
        #     'VintTotal'  : <float>,
        #     'VintCoulomb': <float>,
        #     'RHOint'     : <float>,
        #     'types': [
        #        {
        #            'atom_type': <int>,
        #            'Z'        : <int>,
        #            'RWS'      : <float>,
        #            'data'     : <np.array(..., 4, dtype=float)>
        #        },
        #        ...
        #        ...
        #        ...
        #        ]
        #     }
        self.phagen_data = {'types': []}

    def write(self, filename, prototypical_atoms):
        LOGGER.debug("Writing Phagen input potential file: {}".format(filename))

        def DEPRECATEDappend_atom_potential(atom):
            Z = atom.number
            # Find the right type (Z) in the phagen_data
            itypes = []
            for i, t in enumerate(self.phagen_data['types']):
                if t['Z'] == Z:
                    itypes.append(i)
            # Check now that we have only one type in the list
            # otherwise we do not know yet how to deal with this.
            assert len(itypes) > 0, "Cannot find the data for atom with Z={}".format(Z)
            assert len(itypes) == 1, "Too many datasets for atom with Z={}".format(Z)
            # So far so good, let's write the block
            t = self.phagen_data['types'][itypes[0]]
            s = "{:<7d}{:<10d}{:1.4f}\n".format(
                t['Z'], len(t['data']), t['RWS'])
            line_fmt = "{:+1.8e} " * 4 + "\n"
            for row in t['data']:
                s += line_fmt.format(*row)
            return s

        def append_atom_potential(atom):
            line_fmt = "{:+1.8e} " * 4 + "\n"
            atom_type = atom.get('atom_type')
            assert atom_type != None, "Unable get the atom type!"
            for t in self.phagen_data['types']:
                if t['atom_type'] == atom_type:
                    s = "{:<7d}{:<10d}{:1.4f}\n".format(
                        t['Z'], len(t['data']), t['RWS'])
                    for row in t['data']:
                        s += line_fmt.format(*row)
                    return s

        content = ""
        # Start by writing the header line
        content += "{:.2f}  {:.4f}   {:.2f}\n".format(
            self.phagen_data['VintCoulomb'],
            self.phagen_data['RHOint'],
            self.phagen_data['VintTotal'])
        # Then for each atom in the given prototypical cluster,
        # write the data block
        for atom in prototypical_atoms:
            content += append_atom_potential(atom)

        # Write the content to filename
        try:
            with open(filename, 'r') as fd:
                old_content = fd.read()
        except IOError:
            old_content = ''

        modified = False
        if content != old_content:
            with open(filename, 'w') as fd:
                fd.write(content)
            modified = True

        return modified


class SPRKKRPotential(ForeignPotential):
    def __init__(self, atoms, potfile, *exported_files):
        super().__init__()
        self.atoms = atoms
        self.potfile = potfile
        self.load_sprkkr_atom_types()
        for f in exported_files:
            LOGGER.info("Loading file {}...".format(f))
            # get the IT from the filename
            m=re.match('SPRKKR-IT_(?P<IT>\d+)-PHAGEN.*', os.path.basename(f))
            it = int(m.group('IT'))

            # load the content of the header (2 lines)
            with open(f, 'r') as fd:
                first_line, second_line = [fd.readline().strip()
                                           for _ in range(2)]

            # load Coulomb and Total interstitial potential
            pattern = (r'#\s*VMTZ_TOTAL\s*=\s*(?P<TOTAL>.*?)\s+'
                       r'VMTZ_Coulomb\s*=\s*(?P<COULOMB>.*?)\s+.*$')
            m = re.match(pattern, first_line)
            self.phagen_data.update(VintCoulomb=float(m.group('COULOMB')),
                                    VintTotal=float(m.group('TOTAL')),
                                    RHOint=0.)

            # load Z, Wigner-Seitz radius from 2nd line of header
            type_data = {}
            _ = re.split(r'\s+', second_line.strip("# "))
            type_data.update(atom_type=int(it), Z=int(_[0]), RWS=float(_[3]))

            # load the data
            data = np.loadtxt(f, comments='#')
            type_data.update(data=data)

            self.phagen_data['types'].append(type_data)

        # store the sprkkr type number in the sprkkr_info property of each
        # atom in the provided self.atoms

    def load_sprkkr_atom_types(self):
        def read(content, pattern, types):
            # compile the pattern for regex matching
            pat = re.compile(pattern, re.MULTILINE)
            # get the keys and data as list of strings
            keys = re.split(r'\s+',
                            pat.search(content).group('KEYS').strip(' \n'))
            txt = pat.search(content).group('DATA').strip('\n').split('\n')
            data = []
            for line in txt:
                # Unpacking values
                values_txt = re.split(r'\s+', line.strip())
                data_dict = {}
                for i, _ in enumerate(values_txt):
                    # Type casting
                    value = types[i].__call__(_)
                    data_dict[keys[i]] = value
                # push to the list
                data.append(data_dict)
            return data

        # load info in *.pot file
        LOGGER.info("Loading SPRKKR *.pot file {}...".format(self.potfile))
        with open(self.potfile, 'r') as fd:
            content = fd.read()

        self.sites_data = read(content,
                               (r'^\s*SITES\s*\n((.*\n)+?\s*(?P<KEYS>IQ.*)\n'
                                r'(?P<DATA>(.*\n)+?))\*+'),
                               [int] + [float] * 3)
        self.types_data = read(content,
                               (r'^\s*TYPES\s*\n(\s*(?P<KEYS>IT.*)\n'
                                r'(?P<DATA>(.*\n)+?))\*+'),
                               [int, str] + [int] * 4)
        self.occ_data = read(content,
                             (r'^\s*OCCUPATION\s*\n(\s*(?P<KEYS>IQ.*)\n'
                              r'(?P<DATA>(.*\n)+?))\*+'),
                             [int] * 5 + [float])

        LOGGER.debug("SITES:")
        for _ in self.sites_data:
            LOGGER.debug(_)

        LOGGER.debug("TYPES:")
        for _ in self.types_data:
            LOGGER.debug(_)

        LOGGER.debug("OCCUPATION:")
        for _ in self.occ_data:
            LOGGER.debug(_)

        for site in self.sites_data:
            IQ = site['IQ']
            # Get the Atom at x, y, z
            x = site['QBAS(X)']
            y = site['QBAS(Y)']
            z = site['QBAS(Z)']
            i = get_atom_index(self.atoms, x, y, z, scaled=True)
            for occupation in self.occ_data:
                if occupation['IQ'] == IQ:
                    IT = occupation['ITOQ']
                    atom = self.atoms[i]
                    atom.set('atom_type', IT)
            LOGGER.debug("Site #{} is type #{}, atom {}".format(IQ, IT, atom))



class EmptySphere(Atom):
    def __init__(self, *args, **kwargs):
        Atom.__init__(self, *args, **kwargs)
        self.symbol = 'X'


def get_atom_index(atoms, x, y, z, scaled=False):
    """ Return the index of the atom that is the closest to the coordiantes
    given as parameters.

    :param ase.Atoms atoms: an ASE Atoms object
    :param float x: the x position in angstroms
    :param float y: the y position in angstroms
    :param float z: the z position in angstroms
    :param bool scaled: whether the x,y,z coordinates are scaled

    :return: the index of the atom as an integer
    :rtype: int
    """
    # get all distances
    if scaled:
        positions = atoms.get_scaled_positions()
    else:
        positions = atoms.get_positions()
    d = np.linalg.norm(positions - np.array([x, y, z]), axis=1)
    # get the index of the min distance
    i = np.argmin(d)
    # return the index
    return i


def center_cluster(atoms, invert=False):
    """ Centers an Atoms object by translating it so the origin is roughly
    at the center of the cluster.
    The function supposes that the cluster is wrapped inside the unit cell,
    with the origin being at the corner of the cell.
    It is used in combination with the cut functions, which work only if the
    origin is at the center of the cluster
    :param ase.Atoms atoms: an ASE Atoms object
    :param bool invert: if True, performs the opposite translation
                        (uncentering the cluster)
    """
    for i, cell_vector in enumerate(atoms.get_cell()):
        if invert:
            atoms.translate(0.5*cell_vector)
        else:
            atoms.translate(-0.5*cell_vector)


def cut_sphere(atoms, radius, center=(0, 0, 0)):
    """ Removes all the atoms of an Atoms object outside a sphere with a
    given radius

    :param ase.Atoms atoms: an ASE Atoms object
    :param float radius: the radius of the sphere

    :return: The modified atom cluster
    :rtype: ase.Atoms
    """
    assert radius >= 0, "Please give a positive radius value"
    radii = np.linalg.norm(atoms.positions - center, axis=1)
    indices = np.where(radii <= radius)[0]
    return atoms[indices]


def cut_cylinder(atoms, axis="z", radius=None):
    """ Removes all the atoms of an Atoms object outside a cylinder with a
    given axis and radius

    :param ase.Atoms atoms: an ASE Atoms object
    :param str axis: string "x", "y", or "z". The axis of the cylinder,
                     "z" by default
    :param float radius: the radius of the cylinder

    :return: The modified atom cluster
    :rtype: ase.Atoms
    """
    if axis not in ('z',):
        raise ValueError("axis value != 'z' is not supported yet.")
    X, Y, Z = atoms.positions.T
    R = np.sqrt(X**2 + Y **2)
    T = np.arctan2(Y, X)
    i = np.where(R <= radius)[0]
    return atoms[i]


def cut_cone(atoms, radius, z=0):
    """Shapes the cluster as a cone.

    Keeps all the atoms of the input Atoms object inside a cone of based
    radius *radius* and of height *z*.

    :param atoms: The cluster to modify.
    :type atoms: :py:class:`ase.Atoms`
    :param radius: The base cone radius in :math:`\mathring{A}`. # noqa: W605
    :type radius: float
    :param z: The height of the cone in :math:`\mathring{A}`.    # noqa: W605
    :type z: float
    :return: A new cluster.
    :rtype: :py:class:`ase.Atoms`
    """
    new_atoms = atoms.copy()
    max_theta = np.arctan(radius/(-z))

    u = np.array((0, 0, -z))
    normu = np.linalg.norm(u)
    new_atoms.translate(u)
    indices = []
    for i in range(len(new_atoms)):
        v = new_atoms[i].position
        normv = np.linalg.norm(v)

        _ = np.dot(u, v)/normu/normv
        if _ == 0:
            print(v)
        theta = np.arccos(_)
        if theta <= max_theta:
            indices.append(i)

    new_atoms = new_atoms[indices]
    new_atoms.translate(-u)

    return new_atoms


def cut_plane(atoms, x=None, y=None, z=None):
    """ Removes the atoms whose coordinates are higher (or lower, for a
     negative cutoff value) than the coordinates given for every dimension.

     For example,

     .. code-block:: python

        cut_plane(atoms, x=[-5,5], y=3.6, z=0)
        # every atom whose x-coordinate is higher than 5 or lower than -5,
        # and/or y-coordinate is higher than 3.6, and/or z-coordinate is higher
        # than 0 is deleted.

     :param ase.Atoms atoms: an ASE Atoms object
     :param int x: x cutoff value
     :param int y: y cutoff value
     :param int z: z cutoff value

     :return: The modified atom cluster
     :rtype: ase.Atoms
     """
    dim_names = ('x', 'y', 'z')
    dim_values = [x, y, z]
    for i, (name, value) in enumerate(zip(dim_names, dim_values)):
        assert isinstance(value, (int, float, list, tuple, type(None))), \
            "Wrong type"
        if isinstance(value, (tuple, list)):
            assert len(value) == 2 and np.all(
                [isinstance(el, (int, float, type(None))) for el in value]), \
                "Wrong length"
        else:
            try:
                if value >= 0:
                    dim_values[i] = [-np.inf, value]
                else:
                    dim_values[i] = [value, np.inf]
            except Exception:
                dim_values[i] = [value, value]

        if dim_values[i][0] is None:
            dim_values[i][0] = -np.inf
        if dim_values[i][1] is None:
            dim_values[i][1] = np.inf

    dim_values = np.array(dim_values)

    X, Y, Z = atoms.positions.T
    i0 = np.where(X                     >= dim_values[0, 0])[0]
    i1 = np.where(X[i0]                 <= dim_values[0, 1])[0]
    i2 = np.where(Y[i0][i1]             >= dim_values[1, 0])[0]
    i3 = np.where(Y[i0][i1][i2]         <= dim_values[1, 1])[0]
    i4 = np.where(Z[i0][i1][i2][i3]     >= dim_values[2, 0])[0]
    i5 = np.where(Z[i0][i1][i2][i3][i4] <= dim_values[2, 1])[0]
    indices = np.arange(len(atoms))[i0][i1][i2][i3][i4][i5]

    return atoms[indices]


def hemispherical_cluster(cluster, emitter_tag=0, emitter_plane=0, diameter=0,
                          planes=0, shape='spherical'):

    """Creates and returns a cluster based on an Atoms object and some
    parameters.

    :param cluster: the Atoms object used to create the cluster
    :type  cluster: Atoms object
    :param emitter_tag: the tag of your emitter
    :type  emitter_tag: integer
    :param diameter: the diameter of your cluster in Angströms
    :type  diameter: float
    :param planes: the number of planes of your cluster
    :type  planes: integer
    :param emitter_plane: the plane where your emitter will be starting by 0
                          for the first plane
    :type  emitter_plane: integer

    See :ref:`hemispherical_cluster_faq` for more informations.
    """

    def get_xypos(cluster, ze, symbol=None):
        nmin = np.inf

        for atom in cluster:
            if ze - eps < atom.z < ze + eps and (atom.symbol == symbol or
                                                 symbol is None):
                n = np.sqrt(atom.x**2 + atom.y**2)
                if (n < nmin):
                    nmin = n
                    iatom = atom.index

        pos = cluster.get_positions()[iatom]
        tx, ty = pos[0], pos[1]
        return tx, ty

    cell = cluster.get_cell()

    eps = 0.01              # a useful small value
    c = cell[:, 2].max()    # a lattice parameter
    a = cell[:, 0].max()    # a lattice parameter

    # the number of planes in the cluster
    p = len(np.unique(np.round(cluster.get_positions()[:, 2], 4)))
    # the symbol of your emitter
    symbol = cluster[np.where(cluster.get_tags() == emitter_tag)[0][0]].symbol

    if shape.lower() in ('spherical'):
        assert (diameter != 0 or planes != 0), \
            "At least one of diameter or planes parameter must be use."
    elif shape.lower() in ('cylindrical'):
        assert (diameter != 0 and planes != 0), \
            "Diameter and planes parameters must be defined for cylindrical shape."

    if diameter == 0:
        # calculate the minimal diameter according to the number of planes
        min_diameter = 1+2*(planes*c/p+1)
    else:
        min_diameter = diameter

    # number of repetition in each direction
    rep = int(3*min_diameter/min(a, c))
    #print("rep = ", rep)

    # repeat the cluster
    cluster = cluster.repeat((rep, rep, rep))

    # center the cluster
    center_cluster(cluster)

    # reset the cell
    cluster.set_cell(cell)

    # cut the cluster so that we have a centered surface
    cluster = cut_plane(cluster, z=eps)

    # positions where atoms have the tag of the  emitter_tag
    i = np.where(cluster.get_tags() == emitter_tag)

    # an array of all unique z corresponding to where we have the right
    # atom's tag
    all_ze = np.sort(np.unique(np.round(cluster.get_positions()[:, 2][i], 4)))

    # an array of all unique z
    all_z = np.sort(np.unique(np.round(cluster.get_positions()[:, 2], 4)))

    # calculate the number of planes above the emitter's plane
    n = np.where(all_z == all_z.max())[0][0] - np.where(
        all_z == all_ze.max())[0][0]

    # the height of the emitter's plane
    ze = all_ze.max()

    # if the number of planes above the emitter's plane is smaller than it must
    # be, recalculate n and ze
    while n < emitter_plane:
        all_ze = all_ze[:-1]
        n = np.where(all_z == all_z.max())[0][0] - np.where(
            all_z == all_ze.max())[0][0]
        ze = all_ze.max()

    # values of x and y of the emitter
    tx, ty = get_xypos(cluster, ze, symbol)

    # center the cluster on the emitter
    Atoms.translate(cluster, [-tx, -ty, 0])

    # calculate where to cut to get the right number of planes above the
    # emitter
    z_cut = all_z[np.where(all_z == all_ze.max())[0][0] + emitter_plane]

    # translate the surface at z=0
    Atoms.translate(cluster, [0, 0, -z_cut])

    # cut the planes above those we want to keep
    cluster = cut_plane(cluster, z=eps)

    radius = diameter/2
    if planes != 0:
        # an array of all unique remaining z
        all_z = np.sort(np.unique(np.round(cluster.get_positions()[:, 2], 4)))

        zplan = all_z[-planes]
        xplan, yplan = get_xypos(cluster, zplan)
        radius = np.sqrt(xplan**2 + yplan**2 + zplan**2)

        if diameter != 0 and shape in ('spherical'):
            assert (radius <= diameter/2), ("The number of planes is too high "
                                            "compared to the diameter.")
        radius = max(radius, diameter/2)

    if shape in ('spherical'):
        # cut a sphere in our cluster with the diameter which is indicate in
        # the parameters
        cluster = cut_sphere(cluster, radius=radius + eps)
    elif shape in ('cylindrical'):
        # cut a sphere in our cluster with the diameter which is indicate in
        # the parameters
        cluster = cut_cylinder(cluster, radius=radius + eps)
    else:
        raise NameError('Unkknown shape specifier: \"{}\"'.format(shape))

    if planes != 0:
        # calculate where to cut to get the right number of planes
        positions = np.unique(np.round(cluster.get_positions()[:, 2], 4))
        zcut = np.sort(positions)[::-1][planes-1] - eps

        # cut the right number of planes
        cluster = cut_plane(cluster, z=zcut)

    # an array of all unique remaining z
    all_z = np.sort(np.unique(np.round(cluster.get_positions()[:, 2], 4)))

    assert emitter_plane < len(all_z), ("There are not enough existing "
                                            "plans.")
    ze = all_z[- emitter_plane - 1]         # the z-coordinate of the emitter
    Atoms.translate(cluster, [0, 0, -ze])   # put the emitter in (0,0,0)

    return cluster

def shape_cluster(primitive, emitter_tag=0, emitter_plane=0, diameter=0,
                          planes=0, shape='spherical'):

    """Creates and returns a cluster based on an Atoms object and some
    parameters.

    :param cluster: the Atoms object used to create the cluster
    :type  cluster: Atoms object
    :param emitter_tag: the tag of your emitter
    :type  emitter_tag: integer
    :param diameter: the diameter of your cluster in Angströms
    :type  diameter: float
    :param planes: the number of planes of your cluster
    :type  planes: integer
    :param emitter_plane: the plane where your emitter will be starting by 0
                          for the first plane
    :type  emitter_plane: integer

    See :ref:`hemispherical_cluster_faq` for more informations.
    """
    # We need the radius of the cluster and the number of planes
    if shape.lower() in ('ispherical', 'cylindrical'):
        assert (nplanes != 0 and diameter != 0), "nplanes and diameter cannot be zero for '{}' shape".format(shape)
    elif shape.lower() in ('spherical'):
        if diameter <= 0:
            # find the diameter based on the number of planes
            assert planes != 0, "planes should be > 0"


    n = 3
    natoms = 0
    while True:
        n += 2
        cluster = primitive.copy()
        # Repeat the primitive cell
        cluster = cluster.repeat((n, n, n))
        center_cluster(cluster)

        # Find the emitter closest to the origin
        all_tags = cluster.get_tags()
        are_emitters = all_tags == emitter_tag
        _ie = np.linalg.norm(cluster[are_emitters].positions, axis=1).argmin()
        ie = np.nonzero(are_emitters)[0][_ie]
        # Translate the cluster to this emitter position
        cluster.translate(-cluster[ie].position)
        # cut plane at surface and at bottom
        all_z = np.unique(cluster.positions[:,2])
        try:
            zsurf = all_z[all_z >= 0][emitter_plane]
        except IndexError:
            # There are not enough planes above the emitter
            zsurf = all_z.max()
        try:
            zbottom = all_z[all_z <= 0][::-1][planes - (emitter_plane+1)]
        except IndexError:
            # There are not enough planes below the emitter
            zbottom = all_z.min()
        cluster = cut_plane(cluster, z=[zbottom,zsurf])
        # spherical shape
        if shape.lower() in ('spherical'):
            cluster = cut_sphere(cluster, radius=diameter/2, center=(0,0,zsurf))
        if shape.lower() in ('ispherical'):
            cluster = cut_sphere(cluster, radius=diameter/2, center=(0,0,0))
        elif shape.lower() in ('cylindrical'):
            cluster = cut_cylinder(cluster, radius=diameter/2)
        else:
            raise NameError("Unknown shape")
        cluster.set_cell(primitive.cell)
        if len(cluster) <= natoms:
            break
        else:
            natoms = len(cluster)


    return cluster











