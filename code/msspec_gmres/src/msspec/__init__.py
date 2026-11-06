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
# msspec is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this msspec.  If not, see <http://www.gnu.org/licenses/>.
#
# Source file  : src/msspec/__init__.py
# Last modified: Mon, 27 Sep 2021 17:49:48 +0200
# Committed by : sylvain tricot <sylvain.tricot@univ-rennes1.fr>


import ase

from msspec.misc import LOGGER
from msspec.version import __version__


def init_msspec():
    LOGGER.debug('Initialization of the msspec module')
    ase.atom.names['mt_radius'] = ('mt_radii', 0.)
    ase.atom.names['mt_radius_scale'] = ('mt_radii_scale', 1.)
    ase.atom.names['proto_index'] = ('proto_indices', 1)
    ase.atom.names['mean_square_vibration'] = ('mean_square_vibrations', 0.)
    ase.atom.names['forward_angle'] = ('forward_angles', 20.)
    ase.atom.names['backward_angle'] = ('backward_angles', 20.)
    ase.atom.names['RA_cut_off'] = ('RA_cuts_off', 1)
    ase.atom.names['atom_type'] = ('atom_types', None)
    ase.atoms.Atoms.absorber = None

init_msspec()
