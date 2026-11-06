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
# Source file  : src/msspec/parameters.py
# Last modified: Tue, 15 Feb 2022 15:37:28 +0100
# Committed by : Sylvain Tricot <sylvain.tricot@univ-rennes1.fr>


"""
Module parameters
=================

"""


import re
import textwrap

import ase
import numpy as np
from terminaltables import AsciiTable

from msspec.misc import get_level_from_electron_configuration
from msspec.misc import LOGGER
from msspec.misc import UREG
from msspec.misc import XRaySource
from msspec.utils import ForeignPotential
from msspec.utils import SPRKKRPotential

class Parameter(object):
    def __init__(self, name, types=None, limits=(None, None),
                 unit=None, allowed_values=None, default=None,
                 pattern=None, fmt='s', binding=None, private=False, group="",
                 doc='Sorry, no more help for this parameter.'):
        self.name = name
        self.group = group
        if isinstance(types, (tuple, list)):
            self.allowed_types = types
        else:
            self.allowed_types = (types,)
        self.low_limit, self.high_limit = limits
        self.unit = unit
        self.allowed_values = allowed_values
        self.default = default
        self.pattern = pattern
        self.fmt = fmt
        self.binding = binding
        self.private = private
        self.docstring = textwrap.dedent(doc)
        self._value = None

        self.value = default

    def serialize(self):
        data = {
            'name': self.name,
            'types': self.allowed_types,
            'limits': (self.low_limit, self.high_limit),
            'unit': self.unit,
            'allowed_values': self.allowed_values,
            'default': self.default,
            'pattern': self.pattern,
            'fmt': self.fmt,
            # 'binding': self.binding,
            'private': self.private,
            'group': self.group,
            'doc': self.docstring
        }
        return data

    def convert(self, value):
        if hasattr(value, 'magnitude'):
            if self.unit != None:
                return value.to(self.unit).magnitude
        return value

    def assert_message(self, msg, *args):
        s = '\'{}\': {}'.format(self.name, msg)
        s = s.format(*args)

        try:
            allowed_types = '\n'.join(
                ['- ' + str(_) for _ in self.allowed_types])
        except TypeError:
            allowed_types = self.allowed_types

        try:
            allowed_values = '\n'.join(
                ['- ' + str(_) for _ in self.allowed_values])
        except TypeError:
            allowed_values = self.allowed_values

        data = [
            ['PROPERTY', 'VALUE'],
            ['Name', '\'{}\''.format(self.name)],
            ['Allowed types', '{}'.format(allowed_types)],
            ['Limits', '{} <= value <= {}'.format(self.low_limit,
                                                    self.high_limit)],
            ['Unit', '{}'.format(self.unit)],
            ['Allowed values', '{}'.format(allowed_values)],
            ['Default', '{}'.format(self.default)]
        ]
        t = AsciiTable(data)
        table = '\n'.join(['\t' * 2 + _ for _ in t.table.split('\n')])

        s = "{}\n\n{}\n{}".format(s, table, self.docstring)

        return s

    def check(self, value):

        # convert if a unit was given
        _value = self.convert(value)

        if self.allowed_types != None:
            # val = value
            # if hasattr(value, 'magnitude'):
            #    val = value.magnitude
            try:
                if isinstance(_value, bool):
                    assert bool in self.allowed_types
                assert isinstance(_value, self.allowed_types)
            except AssertionError:
                raise AssertionError(self.assert_message(
                    'Type {} is not an allowed type for this option '
                    '({} allowed)',
                    str(type(_value)), str(self.allowed_types)))
        # if not isinstance(_value, (list, tuple, np.ndarray)):
        #    _values = [_value,]

        # if isinstance(_value, tuple):
        #    _values = list(_value)


        # for i, val in enumerate(_values):
        for val in np.array(_value).flatten():
            if self.low_limit != None:
                assert val >= self.low_limit, self.assert_message(
                    'Value must be >= {:s}',
                    str(self.low_limit))
            if self.high_limit != None:
                assert val <= self.high_limit, self.assert_message(
                    'Value must be <= {:s}',
                    str(self.high_limit))
            if self.allowed_values != None: #  and isinstance(val, tuple(type(x) for x in self.allowed_values)):
                assert val in self.allowed_values, self.assert_message(
                    'This value is not allowed. Please choose between '
                    'one of {:s}',
                    str(self.allowed_values))
            if self.pattern != None:
                p = re.compile(self.pattern)
                m = p.match(val)
                assert m != None, self.assert_message(
                    'Wrong string format')
                # _value[i] = val

        return _value

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value):
        v = self.check(value)
        self._value = v
        # if hasattr(value, 'magnitude'):
        #    self._value = value.magnitude
        # else:
        #    self._value = value
        if self.binding:
            new_value = None
            try:
                new_value = self.binding(self)
            except AttributeError as err:
                LOGGER.warning(err)
            if new_value is not None:
                self._value = new_value
                # LOGGER.debug('{:s} = {:s}'.format(self.name, str(self.value)))

    def get(self, *args):
        return self.value

    def set(self, *args):
        value = args[-1]
        self.value = value

    def reset(self):
        self.value = self.default

    def __str__(self):
        fmt = '{:' + self.fmt + '}'
        try:
            return fmt.format(self._value)
        except ValueError:
            return str(self._value)

    def __len__(self):
        try:
            return len(self._value)
        except TypeError:
            return 1

class BaseParameters(object):
    __isfrozen = False

    def __init__(self):
        self._parameters = []

    def add_parameters(self, *parameters):
        cls = self.__class__

        for p in parameters:
            p.group = self.__class__.__name__
            try:
                p.binding = getattr(self, 'bind_' + p.name)
            except AttributeError:
                pass

            if not p.private:
                setattr(self, 'set_' + p.name, p.set)
                setattr(self, 'get_' + p.name, p.get)
                setattr(cls, p.name, property(fset=p.set, fget=p.get,
                                              doc=p.docstring))
            self._parameters.append(p)

    def get_parameter(self, name):
        retval = None
        for p in self._parameters:
            if p.name == name:
                retval = p
        if retval is None:
            raise KeyError('No such parameter: {}'.format(name))
        return retval

    def set_parameter(self, name, value, force=False):
        found = False
        for p in self._parameters:
            if p.name == name:
                if p.private:
                    if force:
                        LOGGER.debug("forcing the value of the private "
                                     "parameter '%s' to the value = %s",
                                     p.name, value)
                    else:
                        err_msg = ("Cannot change the value of the private "
                                   "parameter '{}'").format(p.name)
                        LOGGER.error("Cannot change the value of a private "
                                     "parameter")
                        raise NameError(err_msg)
                p.value = value
                found = True
        if found is False:
            err_msg = 'No such parameter: {}'.format(name)
            LOGGER.error("Unknwon parameter's name!")
            raise ValueError(err_msg)

    def freeze(self, frozen=True):
        self.__isfrozen = frozen

    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            data = '\n'.join(
                ["\t\t- '{}'".format(p.name) for p in self._parameters])
            err_msg = """
            '{}' is not an allowed attribute of {} class.
            Please use one of:\n{}""".format(key, self.__class__.__name__,
                                             data)
            LOGGER.error('Unknown attribute!')
            raise AttributeError(err_msg)
        object.__setattr__(self, key, value)

    def __iter__(self):
        for x in range(len(self._parameters)):
            yield self._parameters[x]

    def __getitem__(self, index):
        return self._parameters[index]



class PhagenParameters(BaseParameters):
    def __init__(self):
        parameters = (
            Parameter('version', types=(str,), fmt='>3s', default='1.1'),
            Parameter('calctype', allowed_values=('xpd', 'xas', 'aed', 'led',
                                                  'rex', 'els', 'e2e'#, 'dos'
                                                 ),
                      types=(str,), default='xpd'),
            Parameter('expmode', allowed_values=('cis', 'cfs', 'cel'),
                      types=(str,), default='cis'),
            Parameter('coor', allowed_values=('angs', 'au'), types=(str,),
                      default='angs'),
            Parameter('enunit', allowed_values=('ryd', 'eV'), fmt='>3s',
                      types=(str,), default='ryd'),
            Parameter('einc', types=(int, float), limits=(0, None), fmt='.1f',
                      default=700.),
            Parameter('esct', types=(int, float), limits=(0, None), fmt='.1f',
                      default=580.),
            Parameter('scangl', types=(int, float), limits=(0, 360), fmt='.2f',
                      default=0.),
            Parameter('lambda', types=(int, float), fmt='.1f', default=0.),
            Parameter('emin', types=(int, float), limits=(0, None), fmt='.4f',
                      default=13.5236),
            Parameter('emax', types=(int, float), limits=(0, None), fmt='.4f',
                      default=13.5236),
            Parameter('delta', types=(int, float), limits=(0, None), fmt='.4f',
                      default=0.3),
            Parameter('cip', types=(int, float), limits=(0, None), fmt='.4f',
                      default=0.),
            Parameter('potgen', allowed_values=('in', 'ex'), types=(str,),
                      default='in'),
            Parameter('potype', allowed_values=('hdrel', 'hedin', 'xalph',
                                                'dhrel', 'dhcmp',
                                                'lmto', 'msf', 'spkkr'
                                                ),
                      types=(str,), fmt='>5s', default='hedin'),
            Parameter('relc', allowed_values=('nr', 'sr', 'so'), types=(str,),
                      default='nr'),
            Parameter('norman', allowed_values=('stdcrm', 'scaled', 'extrad'),
                      types=(str,), default='stdcrm'),
            Parameter('ovlpfac', types=(int, float), limits=(0, None),
                      fmt='.4f', default=0.),
            Parameter('gamma', types=(int, float), limits=(None, None),
                      fmt='.2f', default=0.03),
            Parameter('charelx', allowed_values=('ex', 'gs'), types=(str,),
                      default='gs'),
            Parameter('ionzst', allowed_values=('neutral', 'ionic'),
                      types=(str,), fmt='>7s', default='neutral'),
            Parameter('eikappr', allowed_values=('yes', 'no'), types=(str,),
                      fmt='>3s', default='no'),
            Parameter('db', types=(int, float), fmt='.2f', default=0.01),
            Parameter('optrsh', allowed_values=('y', 'n'), types=(str,),
                      default='n'),
            Parameter('rsh', types=(float,), limits=(0., None), fmt='.1f',
                      default=0.9),
            Parameter('lmax_mode', types=(int,), limits=(0, 2), fmt='d',
                      default=2),
            Parameter('lmaxt', types=(int,), limits=(0, None), fmt='d',
                      default=20),
            Parameter('edge', types=(str,), fmt='<2s', default='k'),
            Parameter('edge1', types=(str,), fmt='<2s',
                      default='l1'),
            Parameter('edge2', types=(str,), fmt='<2s', default='l1'),
            Parameter('l2h', types=(int,), fmt='d', default=4),

            Parameter('ionicity', types=dict, default={}),
            Parameter('noproto', allowed_values=('.true.', '.false.'),
                      types=(str,), fmt='>7s', default='.true.'),
            #Parameter('absorber', types=(int,), limits=(1, None), fmt='d', default=1),
            #Parameter('nosym', types=(str,), allowed_values=('.true.', '.false.'), fmt='s', default='.true.'),
            #Parameter('outersph', types=(str,), allowed_values=('.true.', '.false.'), fmt='s', default='.false.'),
            Parameter('atoms', types=ase.atoms.Atoms, default=ase.atoms.Atoms()))
        BaseParameters.__init__(self)
        self.add_parameters(*parameters)
        self.freeze()

class PhagenMallocParameters(BaseParameters):
    def __init__(self):
        parameters = (
            Parameter('nat_', types=int, limits=(1, None), default=1550),
            Parameter('ua_', types=int, limits=(1, None), default=1550),
            Parameter('neq_', types=int, limits=(1, None), default=48),
            Parameter('rdx_', types=int, limits=(1, None), default=1500),
            Parameter('lmax_', types=int, limits=(1, None), default=60),
            Parameter('nef_', types=int, limits=(1, None), default=200),
            Parameter('lexp_', types=int, limits=(1, None), default=10),
            Parameter('nep_', types=int, limits=(1, None), default=1000),
        )
        BaseParameters.__init__(self)
        self.add_parameters(*parameters)
        self.freeze()

class SpecParameters(BaseParameters):
    def __init__(self):
        parameters = (
            Parameter('calctype_spectro', types=str,
                      allowed_values=('PHD', 'AED', 'XAS', 'LED', 'ACS', 'EIG'),
                      default='PHD'),
            Parameter('calctype_ispin', types=int, limits=(0, 1), default=0,
                      fmt='d'),
            Parameter('calctype_idichr', types=int, limits=(0, 2), default=0,
                      fmt='d'),
            Parameter('calctype_ipol', types=int, limits=(-1, 2), default=0,
                      fmt='d'),
            Parameter('calctype_iamp', types=int, limits=(0, 1), default=1,
                      fmt='d'),

            Parameter('ped_li', types=str, default='1s'),
            Parameter('ped_so', types=str, default='1/2'),
            Parameter('ped_initl', types=int, limits=(-1, 2), default=2,
                      fmt='d'),
            Parameter('ped_iso', types=int, limits=(-1, 2), default=0,
                      fmt='d'),
            Parameter('ped_iphi', types=int, limits=(-1, 3), default=3,
                      fmt='d'),
            Parameter('ped_itheta', types=int, limits=(-1, 3), default=0,
                      fmt='d'),
            Parameter('ped_ie', types=int, limits=(0, 4), default=0,
                      fmt='d'),
            Parameter('ped_ifthet', types=int, limits=(0, 1), default=0,
                      fmt='d'),
            Parameter('ped_nphi', types=int, limits=(1, None), default=360,
                      fmt='d'),
            Parameter('ped_ntheta', types=int, limits=(1, None), default=45,
                      fmt='d'),
            Parameter('ped_ne', types=int, limits=(1, None), default=1,
                      fmt='d'),
            Parameter('ped_nfthet', types=int, limits=(1, None), default=1,
                      fmt='d'),
            Parameter('ped_phi0', types=float, limits=(0., 360.), default=0.,
                      fmt='.2f'),
            Parameter('ped_theta0', types=float, limits=(-360., 360.),
                      default=-70., fmt='.2f'),
            Parameter('ped_e0', types=float, limits=(0., None), default=316.4,
                      fmt='.2f'),
            Parameter('ped_r0', types=float, default=0.5,
                      fmt='.3f'),
            Parameter('ped_phi1', types=float, limits=(0., 360.), default=0.,
                      fmt='.2f'),
            Parameter('ped_theta1', types=float, limits=(-360., 360.),
                      default=69., fmt='.2f'),
            Parameter('ped_e1', types=float, limits=[0., None], default=316.4,
                      fmt='.2f'),
            Parameter('ped_r1', types=float, default=-1.0,
                      fmt='.3f'),
            Parameter('ped_thlum', types=float, default=-55.0,
                      fmt='.2f'),
            Parameter('ped_philum', types=float, default=0.,
                      fmt='.2f'),
            Parameter('ped_elum', types=float, default=1253.6,
                      fmt='.2f'),
            Parameter('ped_imod', types=int, default=1,
                      fmt='d'),
            Parameter('ped_imoy', types=int, default=0,
                      fmt='d'),
            Parameter('ped_accept', types=float, default=0.,
                      fmt='.2f'),
            Parameter('ped_ichkdir', types=int, default=0,
                      fmt='d'),

            Parameter('leed_iphi', types=int, limits=(-1, 3), default=-1,
                      fmt='d'),
            Parameter('leed_itheta', types=int, limits=(-1, 3), default=0,
                      fmt='d'),
            Parameter('leed_ie', types=int, limits=(0, 4), default=0,
                      fmt='d'),
            Parameter('leed_ifthet', types=int, limits=(0, 1), default=0,
                      fmt='d'),
            Parameter('leed_nphi', types=int, limits=(1, None), default=1,
                      fmt='d'),
            Parameter('leed_ntheta', types=int, limits=(1, None), default=140,
                      fmt='d'),
            Parameter('leed_ne', types=int, limits=(1, None), default=1,
                      fmt='d'),
            Parameter('leed_nfthet', types=int, limits=(1, None), default=1,
                      fmt='d'),
            Parameter('leed_phi0', types=float, limits=(0., 360.), default=0.,
                      fmt='.2f'),
            Parameter('leed_theta0', types=float, limits=(-360., 360.),
                      default=-70., fmt='.2f'),
            Parameter('leed_e0', types=float, limits=(0., None), default=100.,
                      fmt='.2f'),
            Parameter('leed_r0', types=float, default=0.5,
                      fmt='.3f'),
            Parameter('leed_phi1', types=float, limits=(0., 360.), default=0.,
                      fmt='.2f'),
            Parameter('leed_theta1', types=float, limits=(-360., 360.),
                      default=69., fmt='.2f'),
            Parameter('leed_e1', types=float, limits=[0., None], default=316.4,
                      fmt='.2f'),
            Parameter('leed_r1', types=float, default=-1.0,
                      fmt='.3f'),
            Parameter('leed_thini', types=float, default=-55.0,
                      fmt='.2f'),
            Parameter('leed_phiini', types=float, default=0.,
                      fmt='.2f'),
            Parameter('leed_imod', types=int, default=1,
                      fmt='d'),
            Parameter('leed_imoy', types=int, default=0,
                      fmt='d'),
            Parameter('leed_accept', types=float, default=0.,
                      fmt='.2f'),
            Parameter('leed_ichkdir', types=int, default=0,
                      fmt='d'),

            Parameter('exafs_edge', types=str, default='L1',
                      pattern=r'((K|K1)|[KLMNO][1-9])', fmt='<2s'),
            Parameter('exafs_initl', types=int, limits=[-1, 2],
                      default=1, fmt='d'),
            Parameter('exafs_thlum', types=float, default=-55.0,
                      fmt='.2f'),
            Parameter('exafs_philum', types=float, default=0.,
                      fmt='.2f'),
            Parameter('exafs_ne', types=int, limits=[1, None],
                      default=1, fmt='d'),
            Parameter('exafs_ekini', types=float, limits=[0, None],
                      default=200., fmt='.2f'),
            Parameter('exafs_ekfin', types=float, limits=[0, None],
                      default=600., fmt='.2f'),
            Parameter('exafs_ephini', types=float, limits=[0, None],
                      default=1486.7, fmt='.2f'),

            Parameter('aed_edgec', types=str, default='L2',
                      pattern=r'((K|K1)|[KLMNO][1-9])', fmt='<2s'),
            Parameter('aed_edgei', types=str, default='M2',
                      pattern=r'((K|K1)|[KLMNO][1-9])', fmt='<2s'),
            Parameter('aed_edgea', types=str, default='M2',
                      pattern=r'((K|K1)|[KLMNO][1-9])', fmt='<2s'),
            Parameter('aed_imult', types=int, limits=[0, 1], default=1,
                      fmt='d'),
            Parameter('aed_mult', types=str, default='1D2',
                      pattern=r'\d[SPDFG]\d', fmt='s'),
            Parameter('aed_iphi', types=int, limits=(-1, 3), default=0,
                      fmt='d'),
            Parameter('aed_itheta', types=int, limits=(-1, 3), default=1,
                      fmt='d'),
            Parameter('aed_ifthet', types=int, limits=(0, 1), default=0,
                      fmt='d'),
            Parameter('aed_iint', types=int, limits=(0, 3), default=0,
                      fmt='d'),
            Parameter('aed_nphi', types=int, limits=(1, None), default=1,
                      fmt='d'),
            Parameter('aed_ntheta', types=int, limits=(1, None), default=1,
                      fmt='d'),
            Parameter('aed_nfthet', types=int, limits=(1, None), default=1,
                      fmt='d'),
            Parameter('aed_phi0', types=float, limits=(0., 360.), default=0.,
                      fmt='.2f'),
            Parameter('aed_theta0', types=float, limits=(-360., 360.),
                      default=45., fmt='.2f'),
            Parameter('aed_r0', types=float, default=0.5,
                      fmt='.3f'),
            Parameter('aed_phi1', types=float, limits=(0., 360.), default=0.,
                      fmt='.2f'),
            Parameter('aed_theta1', types=float, limits=(-360., 360.),
                      default=70., fmt='.2f'),
            Parameter('aed_r1', types=float, default=-1.0,
                      fmt='.3f'),
            Parameter('aed_imod', types=int, default=1,
                      fmt='d'),
            Parameter('aed_imoy', types=int, default=0,
                      fmt='d'),
            Parameter('aed_accept', types=float, default=1.,
                      fmt='.2f'),
            Parameter('aed_ichkdir', types=int, default=0,
                      fmt='d'),

            Parameter('eigval_ne', types=int, limits=[1, None], default=1,
                      fmt='d'),
            Parameter('eigval_ekini', types=float, limits=[0., None],
                      default=10., fmt='.2f'),
            Parameter('eigval_ekfin', types=float, limits=[0., None],
                      default=10., fmt='.2f'),
            Parameter('eigval_idamp', types=int, limits=[0, 3],
                      default=3, fmt='d'),
            Parameter('eigval_ispectrum_ne', types=int, limits=[0, 1],
                      default=1, fmt='d'),
            Parameter('eigval_ipwm', types=int, limits=[-4, 4], default=4,
                      fmt='d'),
            Parameter('eigval_method', types=str, default='EPSI',
                      allowed_values=['AITK', 'RICH', 'SALZ', 'EPSI', 'EPSG',
                                      'RHOA', 'THET', 'LEGE', 'CHEB', 'OVER',
                                      'DURB', 'DLEV', 'TLEV', 'ULEV', 'VLEV',
                                      'ELEV', 'EULE', 'GBWT', 'VARI', 'ITHE',
                                      'EALG']),
            Parameter('eigval_acc', types=float, limits=[0., None],
                      default=0.001, fmt='.5f'),
            Parameter('eigval_expo', types=float, default=1., fmt='.3f'),
            Parameter('eigval_nmax', types=int, limits=[1, None], default=200,
                      fmt='d'),
            Parameter('eigval_niter', types=int, limits=[1, None], default=10,
                      fmt='d'),
            Parameter('eigval_ntable', types=int, limits=[1, None], default=3,
                      fmt='d'),
            Parameter('eigval_shift', types=float, default=0.,
                      fmt='.3f'),
            Parameter('eigval_ixn', types=int, limits=[1, 5], default=1,
                      fmt='d'),
            Parameter('eigval_iva', types=int, limits=[1, 5], default=1,
                      fmt='d'),
            Parameter('eigval_ign', types=int, limits=[1, 7], default=1,
                      fmt='d'),
            Parameter('eigval_iwn', types=int, limits=[1, 6], default=1,
                      fmt='d'),
            Parameter('eigval_l', types=int, limits=[0, 1], default=0,
                      fmt='d'),
            Parameter('eigval_alpha', types=float, default=1., fmt='.2f'),
            Parameter('eigval_beta', types=float, default=1., fmt='.2f'),

            Parameter('calc_no', types=int, limits=[0, 8], default=1, fmt='d'),
            Parameter('calc_ndif', types=int, limits=[1, 10], default=3,
                      fmt='d'),
            Parameter('calc_ispher', types=int, limits=[0, 1], default=1,
                      fmt='d'),
            Parameter('calc_igr', types=int, limits=[0, 2], default=0, fmt='d'),

            Parameter('calc_iren', types=int, limits=[0, 8], default=1, fmt='d'),
            Parameter('calc_nren', types=int, limits=[1, None], default=1, fmt='d'),
            Parameter('calc_renr', types=float, limits=[None, None], default=1., fmt='.3f'),
            Parameter('calc_reni', types=float, limits=[None, None], default=0., fmt='.3f'),

            Parameter('calc_isflip', types=int, limits=[0, 1], default=0,
                      fmt='d'),
            Parameter('calc_irdia', types=int, limits=[0, 1], default=0,
                      fmt='d'),
            Parameter('calc_itrtl', types=int, limits=[1, 9], default=7,
                      fmt='d'),
            Parameter('calc_itest', types=int, limits=[0, 2], default=0,
                      fmt='d'),
            Parameter('calc_isom', types=int, limits=[0, 2], default=0,
                      fmt='d'),
            Parameter('calc_nonvol', types=int, limits=[0, 1], default=1,
                      fmt='d'),
            Parameter('calc_npath', types=int, limits=[0, None], default=100,
                      fmt='d'),
            Parameter('calc_vint', types=(float,int), default=0., fmt='.2f'),
            Parameter('calc_ifwd', types=int, limits=[0, 1], default=0,
                      fmt='d'),
            Parameter('calc_nthout', types=int, limits=[0, None], default=1,
                      fmt='d'),
            Parameter('calc_ino', types=int, limits=[0, None], default=0,
                      fmt='d'),
            Parameter('calc_ira', types=int, limits=[0, 1], default=0, fmt='d'),
            Parameter('calc_ipw', types=int, limits=[0, 1], default=0, fmt='d'),
            Parameter('calc_ncut', types=int, limits=[0, 10], default=2,
                      fmt='d'),
            Parameter('calc_pctint', types=float, limits=[1e-4, 999.9999],
                      default=0.01, fmt='.4f'),
            Parameter('calc_ipp', types=int, limits=[1, 2], default=1, fmt='d'),
            Parameter('calc_ilength', types=int, limits=[0, 1], default=0,
                      fmt='d'),
            Parameter('calc_rlength', types=float, limits=[0, None],
                      default=10., fmt='.2f'),
            Parameter('calc_unlength', types=str, allowed_values=['ANG', ],
                      default='ANG'),
            Parameter('calc_idwsph', types=int, limits=[0, 1], default=0,
                      fmt='d'),
            Parameter('calc_ispeed', types=int, limits=[None, 1], default=1,
                      fmt='d'),
            Parameter('calc_iatt', types=int, limits=[0, 1], default=1,
                      fmt='d'),
            Parameter('calc_iprint', types=int, limits=[0, 3], default=1,
                      fmt='d'),
            Parameter('calc_idcm', types=int, limits=[0, 2], default=0, fmt='d'),
            Parameter('calc_td', types=float, limits=[0., None], default=420.,
                      fmt='.2f'),
            Parameter('calc_t', types=float, limits=[0., None], default=293.,
                      fmt='.2f'),
            Parameter('calc_rsj', types=float, limits=[0., None], default=1.2,
                      fmt='.2f'),

            Parameter('calc_ilpm', types=int, limits=[-1, 2], default=2,
                      fmt='d'),
            Parameter('calc_xlpm0', types=float, limits=[0., None], default=15.,
                      fmt='.2f'),

            Parameter('input_data', types=str, default='../input/spec.dat'),
            Parameter('input_unit00', types=int, default=5, fmt='d'),
            Parameter('input_tl', types=str, default='../output/tmatrix.tl'),
            Parameter('input_unit01', types=int, default=1, fmt='d'),
            Parameter('input_rad', types=str, default='../output/tmatrix.rad'),
            Parameter('input_unit02', types=int, default=3, fmt='d'),
            Parameter('input_cluster', types=str,
                      default='../output/cluster.clu'),
            Parameter('input_unit03', types=int, default=4, fmt='d'),
            Parameter('input_adsorbate', types=str,
                      default='../input/adsorbate.pos'),
            Parameter('input_unit04', types=int, default=2, fmt='d'),
            Parameter('input_kdirs', types=str, default='../input/kdirs.dat'),
            Parameter('input_unit05', types=int, default=11, fmt='d'),

            Parameter('input2_tl', types=str, default='../output/tmatrix.tl'),
            Parameter('input2_unit06', types=int, default=12, fmt='d'),
            Parameter('input2_rad', types=str, default='../output/tmatrix.rad'),
            Parameter('input2_unit07', types=int, default=13, fmt='d'),
            Parameter('input2_kdirs', types=str, default='../input/kdirs.dat'),
            Parameter('input2_unit08', types=int, default=14, fmt='d'),

            Parameter('output_log', types=str, default='../output/spec.log'),
            Parameter('output_unit09', types=int, default=6, fmt='d'),
            Parameter('output_res', types=str, default='../output/results.dat'),
            Parameter('output_unit10', types=int, default=9, fmt='d'),
            Parameter('output_sf', types=str, default='../output/facdif1.dat'),
            Parameter('output_unit11', types=int, default=8, fmt='d'),
            Parameter('output_augclus', types=str,
                      default='../output/augclus.clu'),
            Parameter('output_unit12', types=int, default=10, fmt='d'),

            Parameter('extra_atoms', types=ase.atoms.Atoms,
                      default=ase.atoms.Atoms()),
            Parameter('extra_nat', types=int, limits=[1, None], default=1,
                      fmt='d'),
            Parameter('extra_energies', types=(list, tuple), default=[0., ]),
            Parameter('extra_nlmax', types=int, default=10),
            Parameter('extra_level', types=str, default='1s'),
            Parameter('extra_parameters', types=dict, default={})
        )
        BaseParameters.__init__(self)
        self.add_parameters(*parameters)
        self.freeze()

class SpecMallocParameters(BaseParameters):
    def __init__(self):
        parameters = (
            Parameter('NATP_M', types=int, limits=(1, None), default=20),
            Parameter('NATCLU_M', types=int, limits=(1, None), default=300),
            Parameter('NAT_EQ_M', types=int, limits=(1, None), default=16),
            Parameter('N_CL_L_M', types=int, limits=(1, None), default=1),
            Parameter('NE_M', types=int, limits=(1, None), default=100),
            Parameter('NL_M', types=int, limits=(1, None), default=50),
            Parameter('LI_M', types=int, limits=(1, None), default=7),
            Parameter('NEMET_M', types=int, limits=(1, None), default=2),
            Parameter('NO_ST_M', types=int, limits=(1, None), default=2),
            Parameter('NDIF_M', types=int, limits=(1, None), default=10),
            Parameter('NSO_M', types=int, limits=(1, None), default=2),
            Parameter('NTEMP_M', types=int, limits=(1, None), default=1),
            Parameter('NODES_EX_M', types=int, limits=(1, None), default=3),
            Parameter('NSPIN_M', types=int, limits=(1, None), default=1),
            Parameter('NTH_M', types=int, limits=(1, None), default=2000),
            Parameter('NPH_M', types=int, limits=(1, None), default=2000),
            Parameter('NDIM_M', types=int, limits=(1, None), default=100000),
            Parameter('N_TILT_M', types=int, limits=(1, None), default=11),
            Parameter('N_ORD_M', types=int, limits=(1, None), default=200),
            Parameter('NPATH_M', types=int, limits=(1, None), default=500),
            Parameter('NGR_M', types=int, limits=(1, None), default=10),
        )
        BaseParameters.__init__(self)
        self.add_parameters(*parameters)
        self.freeze()



class GlobalParameters(BaseParameters):
    def __init__(self, phagen_parameters=None, spec_parameters=None):
        parameters = (
            Parameter('spectroscopy', types=str, allowed_values=(
                'PED', 'AED', 'LEED', 'EXAFS', 'EIG'), default='PED',
                doc="""
                There are 4 choices for the spectroscopy option:

                    - '**PED**', for Photo Electron Diffraction
                    - '**AED**', for Auger Electron Diffraction
                    - '**LEED**', for Low Energy Electron Diffraction
                    - '**EXAFS**', for the Extended X-ray Absorption Fine Structure

                Additionally, a 5th keyword **EIG** is used to get deeper information
                about the convergence of the eigen values of multiple scattering
                matrix.

                The value is case insensitive.
                """),
            Parameter('algorithm', types=str, allowed_values=('expansion',
                                                              'inversion',
                                                              'gmres',
                                                              'correlation',
                                                              'power'),
                      default='expansion', doc="""
                      You can choose the algorithm used for the computation of the scattering path operator.

                        - '**inversion**', for the classical matrix inversion method
                        - '**gmres**', for the restarted GMRES solver preconditioned by the exact tree factorization, with dense LU as the fallback (only for spectroscopy='PED')
                        - '**expansion**', for the Rehr-Albers series expansion
                        - '**correlation**', for the correlation-expansion algorithm
                        - '**power**', for the power method approximation scheme (only for spectroscopy='EIG')

                      The series expansion algorithm is well suited for high energy since the number of terms
                      required decreases as the energy increases. The exact solution is obtained by the matrix inversion
                      method but should be preferably used for lower energy.
                      """),
            Parameter('polarization', types=(type(None), str),
                      allowed_values=(None, 'linear_qOz', 'linear_xOy',
                                      'circular'), default=None, doc="""
                                      Specify the polarization of the incident light.

                                        - **None**, for unpolarized light
                                        - '**linear_qOz**' for a polarization vector in the :math:`(\\vec{q}0z)` plane
                                        - '**linear_xOy**' for a polarization vector in the :math:`(x0y)` plane
                                        - '**circular**' for circular dichroism

                                      """),
            Parameter('dichroism', types=(type(None), str),
                      allowed_values=(None, 'natural', 'sum_over_spin',
                                      'spin_resolved'), default=None, doc="""
                                      Used to perform dichroic calculations. The default (None) is to disable this.
                                      """),
            Parameter('spinpol', types=bool, default=False, doc="""
            Enable or disbale spin-resolved calculations.
            """),
            Parameter('folder', types=str, default='./calc', doc="""
            This parameter is the path to the temporary folder used for the calculations. If you do not change this
            parameter between calculations, all the content will be overridden. This is usually not a problem, since the
            whole bunch of data created during a computation is not meant to be saved. But you may want to anyway by
            changing it to another path.

            This folder is not automatically removed after a computation. It is removed when calling the :meth:`shutdown`
            calculator method:

            .. code-block:: python

                calc = MSSPEC() # the './calc' folder is created
                # do your calculation here
                calc.shutdown() # the folder is removed


            """)
        )
        BaseParameters.__init__(self)
        self.add_parameters(*parameters)
        self.phagen_parameters = phagen_parameters
        self.spec_parameters = spec_parameters
        self.spec_parameters.extra_parameters['global'] = self
        self.freeze()

    def bind_spectroscopy(self, p):
        mapping = {
            'PED': ('xpd', 'PHD'),
            'AED': ('aed', 'AED'),
            'LEED': ('led', 'LED'),
            'EXAFS': ('xas', 'XAS'),
            'EIG': ('xpd', 'EIG'),
        }
        phagen_calctype, spec_calctype = mapping[p.value]
        self.phagen_parameters.calctype = phagen_calctype
        self.spec_parameters.calctype_spectro = spec_calctype

    def bind_polarization(self, p):
        if p.value is None:
            ipol = 0
        elif p.value == 'linear_qOz':
            ipol = 1
        elif p.value == 'linear_xOy':
            ipol = -1
        elif p.value == 'circular':
            ipol = 2
        self.spec_parameters.calctype_ipol = ipol

    def bind_spinpol(self, p):
        if p.value == True:
            LOGGER.error('Spin polarization is not yet enabled in the Python version.')
            raise NotImplementedError

    def bind_dichroism(self, p):
        if p.value is not None:
            LOGGER.error('Dichroism is not yet enabled in the Python version.')
            raise NotImplemented


class MuffintinParameters(BaseParameters):
    def __init__(self, phagen_parameters, spec_parameters):
        parameters = (
            Parameter('charge_relaxation', types=bool, default=True, doc="""
            Used to specify whether the charge density of the photoabsorbing atom, which is used
            to construct the potential, is allowed to relax around the core hole or not.
            """),
            Parameter('ionicity', types=dict, default={}, doc="""
            A dictionary to specify the ionicity of each kind of atoms. If empty, the atoms are considered to be
            neutral. Otherwise, each key must be a chemical symbol and the corresponding value should be the fraction
            of electrons added (negative) or substracted (positive) with respect to neutrality.
            As an example for a LaFeO\ :sub:`3` cluster::

            >>> calc.muffintin_parameters.ionicity = {'La': 0.15, 'Fe': 0.15, 'O': -0.1}

            means that 0.15 electrons have been substracted from La, likewise from Fe. Neutrality implies that 0.3
            electrons have to be added to oxygen atoms.

            """),
            Parameter('relativistic_mode', types=str,
                      allowed_values=('non_relativistic', 'scalar_relativistic',
                                      'spin_orbit_resolved'),
                      default='non_relativistic', doc="""
                      To tell whether to use the scalar relativistic approximation or not.
                      """),
            Parameter('radius_overlapping', types=float, limits=(0., 1.),
                      default=0., doc="""
                      to allow atomic spheres to overlapp with each other. The value is a percentage, 1. means 100%.
                      """),
            Parameter('interstitial_potential', types=(int, float),
                      unit=UREG.eV, default=0., doc="""
                      The average interstitial potential (or inner potential) expressed in eV. It is used to compute
                      the refraction at the surface.
                      """),
            Parameter('hydrogen_radius', types=(int, float), default=0.9,
                      limits=(0., None), unit=UREG.angstroms, doc="""
                      The program can have difficulties to find the radius of the hydrogen atom (small atom). You can
                      specify here a value for the radius. If you set it to 'None', the calculation of the muffin-tin
                      radius of H atoms will be left to the program.
                      """)
        )
        BaseParameters.__init__(self)
        self.add_parameters(*parameters)
        self.phagen_parameters = phagen_parameters
        self.spec_parameters = spec_parameters
        self.freeze()

    def bind_charge_relaxation(self, p):
        mapping = {True: 'ex', False: 'gs'}
        self.phagen_parameters.charelx = mapping[p.value]

    def bind_ionicity(self, p):
        self.phagen_parameters.ionicity = p.value

    def bind_relativistic_mode(self, p):
        mapping = {'non_relativistic': 'nr', 'scalar_relativistic': 'sr',
                   'spin_orbit_resolved': 'so'}
        self.phagen_parameters.relc = mapping[p.value]

    def bind_radius_overlapping(self, p):
        self.phagen_parameters.ovlpfac = p.value

    def bind_interstitial_potential(self, p):
        self.spec_parameters.calc_vint = p.value

    def bind_hydrogen_radius(self, p):
        self.phagen_parameters.rsh = p.value

class TMatrixParameters(BaseParameters):
    def __init__(self, phagen_parameters):
        parameters = (
            Parameter('potential', types=(str, ForeignPotential),
                      default='muffin_tin',
                      doc="""
                      This option allows to choose which kind of potential is
                      used to compute the T-Matrix. For now, only the internal
                      *Muffin-Tin* potential is supported.""",
                      ),
            Parameter('exchange_correlation', types=str,
                      allowed_values=('hedin_lundqvist_real',
                                      'hedin_lundqvist_complex',
                                      'x_alpha_real',
                                      'dirac_hara_real', 'dirac_hara_complex'),
                      default='hedin_lundqvist_complex', doc="""
                      Set the type of exchange and correlation potential that will be used.
                      """),
            Parameter('imaginery_part', types=(int, float), default=0., doc="""
            This value is added to complex potentials to account for core-hole lifetime and experimental resolution
            broadening effects.
            """),
            Parameter('lmax_mode', types=str,
                      allowed_values=('max_ke', 'true_ke', 'imposed'),
                      default='true_ke', doc="""
                      This allows to control the number of basis functions used to expand the wave function on each
                      atom. It can be:

                      - '**imposed**', and will be equal to the *lmaxt* parameter (see below). It is therefore independent
                        on both the energy and atom type.
                      - '**max_ke**', in this case :math:`l_{max}` is computed according to the formula
                        :math:`l_{max} = kr_{at} + 1` where :math:`k=E^{1/2}_{max}` with :math:`E_{max}` being the
                        maximum kinetic energy. In this case :math:`l_{max}` is independent of the energy but
                        depends on the atom number.
                      - '**true_ke**', in this case :math:`l_{max}` is computed according to the formula
                        :math:`l_{max} = kr_{at} + 1` where :math:`k=E^{1/2}_k` with :math:`E_k` being the kinetic
                        energy. In this case :math:`l_{max}` depends both on the energy and the atom number.

                      """),
            Parameter('lmaxt', types=int, limits=(1, None), default=19, doc="""
            The value of :math:`l_{max}` if *lmax_mode = 'imposed'*
            """),
            Parameter('tl_threshold', types=(type(None), float), default=None, doc="""
            This option allows to control the number of basis function by defining a threshold value for the *tl*.
            For example::

            >>> calc.tmatrix_parameters.tl_threshold = 1e-6

            will remove all *tl* with a value :math:`< 1.10^{-6}`

            .. note::
                This option is compatible with any modes of the *lmax_mode* option.

                """),
            Parameter('max_tl', types=(type(None), dict), default=None, doc="""
            *max_tl* is used to sepcify a maximum number of basis functions to use for each kind of atoms. For example,
            in the case of an MgO cluster, you could write::

            >>> calc.muffintin_parameters.max_tl = {'Mg': 20, 'O', 15}

            to tell the program to use at most 20 *tl* for Mg and 15 for O. It acts like a filter, meaning that if you
            use this option, you are not required to sepcif a value for each kind of atoms in your cluster. You can
            restrict the number of *tl* only for one type of atom for example.

            .. note::
                This option is compatible with any modes of the *lmax_mode* option.

            """)
        )
        BaseParameters.__init__(self)
        self.add_parameters(*parameters)
        self.phagen_parameters = phagen_parameters
        self.freeze()

    def bind_potential(self, p):
        mapping = {SPRKKRPotential: 'spkkr'}
        if isinstance(p.value, str):
            assert p.value.lower() in ('muffin_tin'), p.assert_message(
                "\"{}\" is not an known potential type", p.value)
            self.phagen_parameters.potgen = 'in'
        elif isinstance(p.value, ForeignPotential):
            self.phagen_parameters.potgen = 'ex'
            potype = mapping[p.value.__class__]
            self.phagen_parameters.potype = potype



    def bind_exchange_correlation(self, p):
        potential = self.get_parameter('potential').value
        if potential == 'muffin_tin':
            mapping = {
                'hedin_lundqvist_real': 'hdrel',
                'hedin_lundqvist_complex': 'hedin',
                'x_alpha_real': 'xalph',
                'dirac_hara_real': 'dhrel',
                'dirac_hara_complex': 'dhcmp'
            }
            self.phagen_parameters.potype = mapping[p.value]
        else:
            mapping = {SPRKKRPotential: 'spkkr'}
            self.phagen_parameters.potype = mapping[potential.value.__class__]

    def bind_potential_file(self, p):
        pass

    def bind_imaginery_part(self, p):
        self.phagen_parameters.gamma = p.value

    def bind_lmax_mode(self, p):
        mapping = {
            'imposed': 0,
            'max_ke': 1,
            'true_ke': 2
        }
        self.phagen_parameters.lmax_mode = mapping[p.value]

    def bind_lmaxt(self, p):
        self.phagen_parameters.lmaxt = p.value

    def bind_max_tl(self, p):
        cluster = self.phagen_parameters.get_parameter('atoms').value
        # issue a warning if a chemical symbol entered in max_tl is not
        # a chemical symbol in the cluster
        if isinstance(p.value, dict):
            for symbol in list(p.value.keys()):
                if symbol not in np.unique(cluster.get_chemical_symbols()):
                    LOGGER.warning('You provided a maximum tl value for '
                                   '\'%s\' atoms, but there is no such '
                                   'chemical symbol in your cluster.', symbol)

class SourceParameters(BaseParameters):
    def __init__(self, global_parameters=None, phagen_parameters=None, spec_parameters=None):
        parameters = (
            Parameter('energy', types=(list, tuple, int, float),
                      limits=(0, None), unit=UREG.eV, doc="""
                      The photon energy in eV.

                      Common  laboratories X-ray source Mg |kalpha| and Al |kalpha| lines are
                      defined in the :py:class:`msspec.misc.XRaySource` class. For example:

                      .. highlight:: python

                      ::

                          >>> from msspec.calculator import MSSPEC
                          >>> calc = MSSPEC()
                          >>> calc.source_parameters.energy = XRaySource.MG_KALPHA
                          >>> print calc.source_parameters.energy
                          1253.6


                      """,
                      default=XRaySource.MG_KALPHA),
            Parameter('theta', types=(int, float), limits=(-180., 180.),
                      unit=UREG.degree, default=-55., doc="""
                      The polar angle of the photon incidence (in degrees). Please refer to
                      :ref:`this figure <ped_full_picture>` for questions regarding the proper
                      orientation.
                      """),
            Parameter('phi', types=(int, float), limits=(-180., 180.),
                      unit=UREG.degree, default=0., doc="""
                      The azimuthal angle of the photon incidence (in degrees). Please refer to
                      :ref:`this figure <ped_full_picture>` for questions regarding the proper
                      orientation.
                      """),
        )
        BaseParameters.__init__(self)
        self.add_parameters(*parameters)
        self.global_parameters = global_parameters
        self.phagen_parameters = phagen_parameters
        self.spec_parameters = spec_parameters
        self.freeze()

    def bind_energy(self, p):
        spectro = self.global_parameters.spectroscopy
        if spectro in ('PED',):
            assert isinstance(p.value, (int, float)), (
                'Only a single value for the light is allowed in '
                'PhotoElectron Diffraction spectroscopy')
            self.spec_parameters.ped_elum = float(p.value)
            LOGGER.info('Incomming photon energy set to %s', p.value * p.unit)
        else:
            LOGGER.warning('Setting the source energy is pointless in %s '
                           'spectroscopy. Statement ignored.', spectro)

    def bind_theta(self, p):
        spectro = self.global_parameters.spectroscopy
        if spectro in ('PED', 'EXAFS', 'LEED'):
            self.spec_parameters.ped_thlum = float(p.value)
            self.spec_parameters.leed_thini = float(p.value)
            self.spec_parameters.exafs_thlum = float(p.value)
            LOGGER.info('Incomming photon beam theta angle set to %s',
                        p.value * p.unit)
        else:
            LOGGER.warning('Setting the source beam polar (theta) angle is '
                           'pointless in %s spectroscopy. Statement ignored.',
                           spectro)

    def bind_phi(self, p):
        spectro = self.global_parameters.spectroscopy
        if spectro in ('PED', 'EXAFS', 'LEED'):
            self.spec_parameters.ped_philum = float(p.value)
            self.spec_parameters.leed_phiini = float(p.value)
            self.spec_parameters.exafs_philum = float(p.value)
            LOGGER.info('Incomming photon beam theta angle set to %s',
                        p.value * p.unit)
        else:
            LOGGER.warning('Setting the source beam azimutal (phi) angle is '
                    'pointless in %s spectroscopy. Statement ignored.',
                           spectro)

class DetectorParameters(BaseParameters):
    def __init__(self, global_parameters, phagen_parameters, spec_parameters):
        parameters = (
            Parameter('angular_acceptance', types=(int, float),
                      unit=UREG.degree, limits=(0., None), default=0.,
                      doc="""
                      The angular acceptance of the detector in degrees used
                      when the *average_sampling* option is enabled below.
                      """),
            Parameter('average_sampling', types=(type(None), str),
                      allowed_values=(None, 'low', 'medium', 'high'),
                      default=None, doc="""
                      Used to averaged the signal over directions lying in the
                      cone of half-angle *angular_acceptance*. The number of
                      directions to take into account depends on the choosen
                      value:

                      - **None**, for no averaging at all
                      - '**low**', to average over 5 directions
                      - '**medium**', to average over 13 directions
                      - '**high**', to average over 49 directions
                      """),
            Parameter('rotate', types=bool, default=False, doc="""
                      When False, the sample is rotated when doing a scan (the
                      usual way). Otherwise, the detector is rotated.
                      """)
        )
        BaseParameters.__init__(self)
        self.add_parameters(*parameters)
        self.global_parameters = global_parameters
        self.phagen_parameters = phagen_parameters
        self.spec_parameters = spec_parameters
        self.freeze()

    def bind_angular_acceptance(self, p):
        spectro = self.global_parameters.spectroscopy
        if spectro in ('PED', 'LEED', 'AED'):
            self.spec_parameters.ped_accept = float(p.value)
            self.spec_parameters.leed_accept = float(p.value)
            self.spec_parameters.aed_accept = float(p.value)
            LOGGER.info('Detector angular acceptance set to %s', p.value * p.unit)
        else:
            LOGGER.warning('Setting the detector angular acceptance is ignored '
                           'in %s spectroscopy.', spectro)

    def bind_average_sampling(self, p):
        spectro = self.global_parameters.spectroscopy
        if spectro in ('PED', 'LEED', 'AED'):
            imoy = p.allowed_values.index(p.value)
            self.spec_parameters.ped_imoy = imoy
            self.spec_parameters.leed_imoy = imoy
            self.spec_parameters.aed_imoy = imoy
            LOGGER.info('Detector average sampling set to %s', p.value.upper())
        else:
            LOGGER.warning('Setting the average sampling is ignored in %s '
                           'spectroscopy.', spectro)

    def bind_rotate(self, p):
        spectro = self.global_parameters.spectroscopy
        if spectro in ('PED', 'LEED', 'AED'):
            imod = int(not p.value)
            self.spec_parameters.ped_imod = imod
            self.spec_parameters.leed_imod = imod
            self.spec_parameters.aed_imod = imod
            if imod == 0:
                LOGGER.info('The DETECTOR will rotate')
            else:
                LOGGER.info('The SAMPLE will rotate')
        else:
            LOGGER.warning('Setting the detector ratation is ignored in %s '
                           'spectroscopy.', spectro)

class ScanParameters(BaseParameters):
    def __init__(self, global_parameters, phagen_parameters, spec_parameters):
        parameters = (
            Parameter('type', allowed_values=('theta', 'phi', 'theta_phi',
                                              'energy', 'scatf'),
                      types=str, default='theta'),
            Parameter('theta', types=(np.ndarray, list, tuple, int, float),
                      unit=UREG.degree, limits=(-90., 90.),
                      default=np.linspace(0., 88., 45)),
            Parameter('phi', types=(np.ndarray, list, tuple, int, float),
                      unit=UREG.degree, limits=(0., 360.),
                      default=np.array([0.])),
            Parameter('kinetic_energy', types=(list, tuple, int, float),
                      unit=UREG.eV, limits=(0., None),
                      default=200., doc="""
                      if given as a list or tuple:
                        * with 2 elements, 10 points of energy will be generated
                          with the first element as the starting energy and the
                          second element as the stopping energy.
                        * with 3 elements, the first element is the starting energy
                          the second one is the stopping energy and the last one is
                          the number of points.

                      if given as a float or integer, there will be one point
                      for the kinetic energy.
                      """),
            Parameter('ke_array', types=np.ndarray, unit=UREG.eV,
                      default=np.array([200., ]), private=True)
        )
        BaseParameters.__init__(self)
        self.add_parameters(*parameters)
        self.phagen_parameters = phagen_parameters
        self.spec_parameters = spec_parameters
        self.global_parameters = global_parameters

        self.spec_parameters.extra_parameters['scan'] = self
        self.freeze()

    def bind_type(self, p):
        spectro = self.global_parameters.spectroscopy
        calculation_parameters = self.spec_parameters.extra_parameters[
            'calculation']

        scantype = p.value
        if spectro in ('EXAFS',):
            LOGGER.warning('No scan type are allowed in %s spectroscopy. This '
                           'will be ignored.', spectro)
            return
        if spectro in ('AED',) and scantype == 'energy':
            LOGGER.error('Energy scan is not compatible with %s '
                         'spectroscopy.', spectro)
            raise ValueError
        params = {
            'phi': (2, 0, 0, 0),
            'theta': (0, 2, 0, 0),
            'energy': (0, 0, 1, 0),
            'scatf': (0, 1, 0, 1),
            'theta_phi': (-1, 0, 0, 0)}
        self.spec_parameters.ped_iphi = params[scantype][0]
        self.spec_parameters.ped_itheta = params[scantype][1]
        self.spec_parameters.ped_ie = params[scantype][2]
        self.spec_parameters.ped_ifthet = params[scantype][3]

        self.spec_parameters.leed_iphi = params[scantype][0]
        self.spec_parameters.leed_itheta = params[scantype][1]
        self.spec_parameters.leed_ie = params[scantype][2]
        self.spec_parameters.leed_ifthet = params[scantype][3]

        self.spec_parameters.aed_iphi = params[scantype][0]
        self.spec_parameters.aed_itheta = params[scantype][1]
        self.spec_parameters.aed_ifthet = params[scantype][3]

        if scantype == 'scatf':
            self.spec_parameters.ped_theta0 = -360.
            self.spec_parameters.ped_theta1 = 0.
            self.spec_parameters.ped_ntheta = 577
            self.spec_parameters.ped_nfthet = 577

            self.spec_parameters.leed_theta0 = -360.
            self.spec_parameters.leed_theta1 = 0.
            self.spec_parameters.leed_ntheta = 577
            self.spec_parameters.leed_nfthet = 577

            self.spec_parameters.aed_theta0 = -360.
            self.spec_parameters.aed_theta1 = 0.
            self.spec_parameters.aed_ntheta = 577
            self.spec_parameters.aed_nfthet = 577

        if scantype == "scatf":
            calculation_parameters.set_parameter('basis_functions',
                                                 'plane_wave', force=True)
        else:
            calculation_parameters.set_parameter('basis_functions',
                                                 'spherical', force=True)

        LOGGER.info('\'%s\' scan type choosen.', p.value)

    def bind_theta(self, p):
        spectro = self.global_parameters.spectroscopy
        scantype = self.get_parameter('type').value

        #if scantype == 'scatf':
        #    comment = 'with the scattering factor scan type.'
        #if spectro == 'EXAFS':
        #    comment = 'in EXAFS spetroscopy.'
        #if spectro in ('EXAFS',) or scantype in ('scatf',):
        #    msg = 'Setting the theta angle is not possible %s' % comment
        #    LOGGER.error('Incompatible options!')
        #    raise ValueError(msg)

        # p._value = np.array(p.value, dtype=float).flatten()
        arr = np.array(p.value, dtype=float).flatten()

        theta0 = arr[0]
        theta1 = arr[-1]
        ntheta = len(arr)

        self.spec_parameters.ped_theta0 = theta0
        self.spec_parameters.ped_theta1 = theta1
        self.spec_parameters.ped_ntheta = ntheta

        self.spec_parameters.leed_theta0 = theta0
        self.spec_parameters.leed_theta1 = theta1
        self.spec_parameters.leed_ntheta = ntheta

        self.spec_parameters.aed_theta0 = theta0
        self.spec_parameters.aed_theta1 = theta1
        self.spec_parameters.aed_ntheta = ntheta

        LOGGER.info('theta = %s (%d points)', arr * p.unit, len(arr))
        return arr

    def bind_phi(self, p):
        spectro = self.global_parameters.spectroscopy
        scantype = self.get_parameter('type').value

       #if scantype == 'scatf':
       #    comment = 'with scattering factor scan type.'
       #if spectro == 'EXAFS':
       #    comment = 'in EXAFS spetroscopy.'
       #if spectro in ('EXAFS',) or scantype in ('scatf',):
       #    msg = 'Setting the phi angle is not possible %s' % comment
       #    LOGGER.error('Incompatible options')
       #    raise ValueError(msg)

        arr = np.array(p.value, dtype=float).flatten()

        phi0 = arr[0]
        phi1 = arr[-1]
        nphi = len(arr)

        self.spec_parameters.ped_phi0 = phi0
        self.spec_parameters.ped_phi1 = phi1
        self.spec_parameters.ped_nphi = nphi

        self.spec_parameters.leed_phi0 = phi0
        self.spec_parameters.leed_phi1 = phi1
        self.spec_parameters.leed_nphi = nphi

        self.spec_parameters.aed_phi0 = phi0
        self.spec_parameters.aed_phi1 = phi1
        self.spec_parameters.aed_nphi = nphi

        LOGGER.info('phi = %s (%d points)', arr * p.unit, len(arr))
        return arr

    def bind_kinetic_energy(self, p):
        npts = 10
        energies = [0, 0, npts]
        if isinstance(p.value, (list, tuple, np.ndarray)):
            assert len(p.value) <= 3, (
                'kinetic_energy -- if given as a list -- must have at most '
                '3 elements, in which case:\n'
                '* the first element is the starting energy,\n'
                '* the second one is the stopping energy,\n'
                '* and the last one is the number of points.\n'
                'If the last element is omitted, %s points are '
                'assumed' % str(npts))
            energies[0] = float(p.value[0])
            energies[1] = float(p.value[0])
            if len(p.value) > 1:
                energies[1] = float(p.value[1])
            if len(p.value) > 2:
                energies[-1] = int(p.value[-1])
        else:
            energies[0] = float(p.value)
            energies[1] = float(p.value)
            energies[2] = 1

        ke_array_eV = np.linspace(*energies)
        ke_array_ry = (ke_array_eV * p.unit).to('rydberg')

        emin = np.min(ke_array_ry)
        emax = np.max(ke_array_ry)
        npts = len(ke_array_ry)

        assert npts > 0, 'The number of energy points must be >= 1'
        if npts == 1:
            delta = 0.3 * UREG.rydberg
        else:
            delta = (emax - emin) / (npts - 1)

        self.phagen_parameters.emin = emin.magnitude
        self.phagen_parameters.emax = emax.magnitude
        self.phagen_parameters.delta = float(delta.magnitude)

        self.spec_parameters.ped_e0 = energies[0]
        self.spec_parameters.ped_e1 = energies[1]
        self.spec_parameters.ped_ne = energies[2]

        self.spec_parameters.eigval_ekini = energies[0]
        self.spec_parameters.eigval_ekfin = energies[1]
        self.spec_parameters.eigval_ne = energies[2]

        self.set_parameter('ke_array', ke_array_eV, force=True)

        LOGGER.info('Kinetic energy = %s (%d points)', ke_array_eV, npts)

class CalculationParameters(BaseParameters):
    def __init__(self, global_parameters, phagen_parameters, spec_parameters):
        parameters = (
            Parameter('RA_cutoff', types=int, limits=(0, 8), default=1,
                      doc="""
            The Rehr-Albers cut-off parameter which controls the degree of
            sphericity introduced in the description of the basis functions
            used to expand the wave function around each atomic center.
            It is only meaningful for the series expansion algorithm.
            Its value is limited to 8 but it is rarely necessary to go beyond
            2 or 3."""),
            Parameter('scattering_order', types=int, limits=(1, 10), default=3,
                      doc="""
            The scattering order. Only meaningful for the 'expansion' algorithm.
            Its value is limited to 10."""),
            Parameter('renormalization_mode', allowed_values=(None, 'G_n', 'Sigma_n',
                                                              'Z_n', 'Pi_1', 'Lowdin',
                                                              'SigmaPi_2',
                                                              'SigmaPi_2_a25',
                                                              'SigmaPi_2_a75'),
		              types=(type(None), str), default=None,
                      doc="""
            Enable the calculation of the coefficients for the renormalization of
            the multiple scattering series.
            You can choose :math:`G_n`, :math:`\\Sigma_n`, :math:`Z_n`,
            :math:`\\Pi_1`, Lowdin :math:`K^2`, or the experimental fixed
            :math:`n=2`, :math:`\\alpha=1/2` averaged-product blend
            :math:`\\mathrm{SigmaPi}_2`. The latter is available only for PED
            expansion and sums transformed Neumann powers from 0 through
            ``calc_nren`` inclusive. Experimental variants
            ``SigmaPi_2_a25`` and ``SigmaPi_2_a75`` set :math:`\\alpha=0.25`
            and :math:`0.75` respectively."""),
            Parameter('renormalization_omega', types=(int,float,complex),
		              default=1.+0j,
                      doc="""
            The :math:`\\omega` coefficient used to initialize the
            renormalization alogorithm."""),
            Parameter('linear_solver_tolerance', types=(int, float),
                      limits=(1.e-12, 1.e-1), default=1.e-2, fmt='.3e',
                      doc="Relative residual tolerance for the GMRES solver of ``algorithm='gmres'``."),
            Parameter('linear_solver_max_iterations', types=int,
                      limits=(1, 200), default=67, fmt='d',
                      doc="Maximum iteration cap per RHS for ``algorithm='gmres'``; GMRES stops at tolerance (default 67, the largest value used in the test set; at most 200)."),
            Parameter('linear_solver_restart', types=int,
                      limits=(2, 12), default=6, fmt='d',
                      doc="Restart length of the GMRES solver used by ``algorithm='gmres'``."),
            Parameter('linear_solver_fixed_steps', types=int,
                      limits=(-1, 12), default=-1, fmt='d',
                      doc="""Study option of ``algorithm='gmres'``: return the iterate after exactly
            this many GMRES steps instead of solving to the tolerance. 0 returns the preconditioner
            alone, 1 one step, and so on up to the restart length. -1 (default) solves to the tolerance."""),
            Parameter('linear_solver_axial_translation', types=bool, default=True,
                      doc="""Matrix-free action of ``algorithm='gmres'``. True (default): rotate each
            pair onto the z axis, translate along z and rotate back (Rehr and Albers). False: direct
            Gaunt-coefficient evaluation."""),
            Parameter('linear_solver_threads', types=int, limits=(0, None), default=0, fmt='d',
                      doc="""Number of OpenMP threads of the matrix-free action of ``algorithm='gmres'``.
            0 (default) leaves the environment variable MSSPEC_TREE_ACTION_THREADS unchanged."""),
            Parameter('linear_solver_verify', types=bool, default=False,
                      doc="""Compute the true residual ||A x - b|| / ||b|| of the answer of ``algorithm='gmres'``
            with one extra full matrix action, and print it. It does not change the answer."""),
            Parameter('linear_solver_diagnostics', types=bool, default=False,
                      doc="""Extra measurements of ``algorithm='gmres'`` (norms of the kept and omitted tree blocks).
            Implies ``linear_solver_verify``. Slow; for study only."""),
            Parameter('linear_solver_iterate_dump', types=bool, default=False,
                      doc="""Write the whole solution vector after every GMRES step (and the final one) to the file
            msspec_tree_iter.dat of ``algorithm='gmres'``. The files are large."""),
            Parameter('RA_cutoff_damping', types=int, limits=(0, None),
                      default=0, doc="""
            The Rehr-Albers truncation order. If > 0, the *RA_cutoff* is
            decreased by 1 every *i*\ :sup:`th` scatterer until 0, where
            *i* = *RA_cutoff_damping*."""),
            Parameter('spin_flip', types=bool, default=False,
            doc="""
            This parameter tells if spin-flip is authorized or not in the
            scattering process.

            :Note:

            This option works only if the spin polarization is
            enabled in your calculator object (see spinpol_)."""),
            Parameter('integrals', types=str, allowed_values=('all',
                                                              'diagonal'),
                      default='all', doc="""
            This option allows to take into account all four radial integrals
            (R++, R+-, R-+ and R--) in the calculation or only the diagonal
            radial integrals (R++ and R--) which are generally much larger.

            .. note::

                This option works only if the spin polarization is
                enabled in your calculator object.

            """),
            Parameter('path_filtering', types=(type(None), str, tuple, list),
                      allowed_values=(None,
                                      'forward_scattering',
                                      'backward_scattering',
                                      'distance_cutoff',
                                      'plane_wave_normal',
                                      'plane_wave_spin_averaged'),
                      default=None, doc="""
            Used to activate some filters. It is possible to specify several
            of them by grouping them in a tuple or a list. For example::

            >>> my_filters = ('forward_scattering', 'backward_scattering')
            >>> calc.calculation_parameters.path_filtering = my_filters

            """),
            Parameter('off_cone_events', types=int, limits=(0, None), default=1,
                      doc="""
            Used in conjunction with the '*forward_scattering*' filter.
            If the number of scattering processes outside the forward (or
            backward) scattering cone is greater than this number, then the
            path is rejected and its contribution to the scattering path
            operator won’t be computed.
            """),
            Parameter('scattering_order_cutoff', types=int, limits=(0, 10),
                      default=2, doc="""
            Used in conjunction with the ‘*plane_wave_normal*’ filter. It states
            to activate the plane wave approximation (which is fast but
            less accurate) to compute the contribution when the scattering order
            is greater than this value."""),

            Parameter('distance', types=(int, float), limits=(0, None),
                      unit=UREG.angstroms, default=10., doc="""
            Used with the '*distance_cut_off*' filter. Paths whose length is
            larger than this value are simply rejected."""),
            Parameter('vibrational_damping', types=(type(None), str),
                      allowed_values=(None, 'debye_waller', 'averaged_tl'),
                      default='debye_waller', doc="""
            Tells how to compute the effect of atomic vibrations. It can be:

                - '**debye_waller**' for using the Debye Waller model.
                - '**averaged_tl**' to use the more correct averaging over T-matrix elements."""),
            Parameter('temperature', types=(int, float), limits=(0, None),
                      unit=UREG.degK, default=293., doc="""
            The temperature of the cluster. Used when *use_debye_model* = True
            """),
            Parameter('debye_temperature', types=(int, float), limits=(0, None),
                      unit=UREG.degK, default=420., doc="""
            The Debye temperature used for the calculation of the mean square
            displacements if *use_debye_model* = True"""),
            Parameter('use_debye_model', types=bool, default=False,
            doc="""
            No matter the way you compute the effect of atomic vibrations,
            you need the mean square displacements of atoms. It can be computed
            internally following the Debye model if you set this option to True.
            """),
            Parameter('vibration_scaling', types=(int, float),
                      limits=(0., None), default=1.2, doc="""
            Used to simulate the fact that surface atoms vibrate more than
            bulk ones. It is a factor applied to the mean square displacements.
            """),
            Parameter('basis_functions', types=str, allowed_values=(
                'plane_wave', 'spherical'), default='spherical', private=True),
            Parameter('cutoff_factor', types=(int, float),
                      limits=(1e-4, 999.9999), default=0.01, private=False),
            Parameter('mean_free_path', types=(int, float, str),
                      default='SeahDench', #allowed_values=('mono', 'SeahDench'),
                      doc="""
            The electron mean free path value. You can either:
            - Enter a value (in Angströms), in this case any value <=0 will disable the damping
            - Enter the keyword 'mono' to use a formula valid only for monoelemental samples
            - Enter the keyword 'SeahDench' to use the Seah and Dench formula.

            .. note::

                The mean free path is only taken into account when the input T-matrix corresponds
                to a real potential as, when the potential is complex, this damping is taken care
                of by the imaginery part othe potential.

            """),

        )
        BaseParameters.__init__(self)
        self.add_parameters(*parameters)
        self.global_parameters = global_parameters
        self.phagen_parameters = phagen_parameters
        self.spec_parameters = spec_parameters
        self.spec_parameters.extra_parameters['calculation'] = self
        self.freeze()

    def bind_RA_cutoff(self, p):
        self.spec_parameters.calc_no = p.value
        LOGGER.info('Rehr-Albers cutoff parameter set to %s', p.value)

    def bind_scattering_order(self, p):
        self.spec_parameters.calc_ndif = p.value
        LOGGER.info('Scattering order set to %s', p.value)

    def bind_renormalization_mode(self, p):
        iren_map = {None     : 0,
                    'G_n'    : 1,
                    'Sigma_n': 2,
                    'Z_n'    : 3,
                    'Pi_1'   : 4,
                    'Lowdin' : 5,
                    'SigmaPi_2': 6,
                    'SigmaPi_2_a25': 7,
                    'SigmaPi_2_a75': 8}
        # Check that the method is neither 'Z_n' nor 'K^2' for other
        # 'spetroscopy' than EIG
        try:
            if (self.global_parameters.spectroscopy == 'PED'
                and self.global_parameters.algorithm == 'expansion'):
                    # Parameter validation already enforces allowed_values.
                    pass
            elif (self.global_parameters.spectroscopy == 'EIG'
                and self.global_parameters.algorithm == 'inversion'):
                    assert( p.value in p.allowed_values
                            and p.value not in ('SigmaPi_2',
                                                'SigmaPi_2_a25',
                                                'SigmaPi_2_a75') )
            else:
                assert( p.value is None )
            self.spec_parameters.calc_iren = iren_map[p.value]
            LOGGER.info("Renormalization activated with \'{}\' method".format(p.value))
        except AssertionError:
            LOGGER.error("Renormalization method \'{}\' is not compatible with this "
                         "spectroscopy/algorithm ({!r}/{!r})!".format(
                             p.value, self.global_parameters.spectroscopy,
                             self.global_parameters.algorithm))
            exit(1)


    def bind_renormalization_omega(self, p):
        omega = complex(p.value)
        self.spec_parameters.calc_renr = omega.real
        self.spec_parameters.calc_reni = omega.imag
        LOGGER.info("Renormalization omega set to \'{}\'".format(p.value))

    def bind_RA_cutoff_damping(self, p):
        self.spec_parameters.calc_ino = p.value
        LOGGER.info('Rehr-Albers cutoff damping  set to %s', p.value)

    def bind_basis_functions(self, p):
        if p.value == 'plane_wave':
            ispher = 0
        elif p.value == 'spherical':
            ispher = 1
        self.spec_parameters.calc_ispher = ispher
        LOGGER.info('Type of basis functions: \'%s\'', p.value)

    def bind_spin_flip(self, p):
        isflip = int(p.value)
        if self.global_parameters.spinpol is False and isflip:
            err_msg = (
                "'{}' is ignored since the 'spinpol' global parameter is set "
                "to False. Enable spin polarization in the constructor of "
                "your Calculator if you want to use this option."
            ).format(p.name)
            LOGGER.error("Incompatible options!")
            raise ValueError(err_msg)
        self.spec_parameters.calc_isflip = isflip
        LOGGER.info('Spin-flip set to %s', p.value)

    def bind_integrals(self, p):
        irdia = 0 if p.value == 'all' else 1
        self.spec_parameters.calc_irdia = irdia
        LOGGER.info('Radial integrals taken into account: %s', p.value)
        if self.global_parameters.spinpol is False:
            LOGGER.warning(
                "'{}' is ignored since the 'spinpol' global parameter is set "
                "to False. Enable spin polarization in the constructor of "
                "your Calculator if you want to use this option.".format(p.name))

    def bind_path_filtering(self, p):
        ifwd = ipw = ilength = 0
        ipp = 1

        if p.value != None:
            if ('plane_wave_spin_averaged' in p.value and 'plane_wave_normal' in
                p.value):
                err_msg = (
                    "Only one plane wave filter (either 'plane_wave_normal' or "
                    "'plane_wave_spin_averaged') can be used at a time (along "
                    "with other filters if needed).")
                LOGGER.error('Incompatible options!')
                raise ValueError(err_msg)
            if 'forward_scattering' in p.value:
                ifwd = 1
            if 'backward_scattering' in p.value:
                ifwd = 1
            if 'distance_cutoff' in p.value:
                ilength = 1
            if 'plane_wave_normal' in p.value:
                ipw = 1
                ipp = 1
            if 'plane_wave_spin_averaged' in p.value:
                ipw = 1
                ipp = 2

        self.spec_parameters.calc_ifwd = ifwd
        self.spec_parameters.calc_ipw = ipw
        self.spec_parameters.calc_ilength = ilength
        self.spec_parameters.calc_ipp = ipp
        LOGGER.info('Filters activated: %s', p.value)

    def bind_off_cone_events(self, p):
        self.spec_parameters.calc_nthout = p.value
        LOGGER.info('Off cone events set to %s', p.value)
        path_filtering = self.get_parameter('path_filtering').value
        f = 'forward_scattering'
        if path_filtering is not None and f not in path_filtering:
            LOGGER.warning("'%s' option set but ignored since the "
                           "'%s' filter' is not activated",
                           p.name, f)

    def bind_scattering_order_cutoff(self, p):
        self.spec_parameters.calc_ncut = p.value
        LOGGER.info('Scattering order cutoff set to %s', p.value)
        path_filtering = self.get_parameter('path_filtering').value
        f = 'plane_wave_normal'
        if path_filtering is not None and f not in path_filtering:
            LOGGER.warning("'%s' option set but ignored since the "
                           "'%s' filter' is not activated",
                           p.name, f)

    def bind_cutoff_factor(self, p):
        self.spec_parameters.calc_pctint = float(p.value)
        LOGGER.info('Cutoff factor set tp %s', p.value)

    def bind_distance(self, p):
        self.spec_parameters.calc_rlength = float(p.value)
        LOGGER.info('Distance cutoff set to %s', p.value * p.unit)
        path_filtering = self.get_parameter('path_filtering').value
        f = 'distance_cutoff'
        if path_filtering is not None and f not in path_filtering:
            LOGGER.warning("'%s' option set but ignored since the "
                           "'%s' filter' is not activated",
                           p.name, f)

    def bind_vibrational_damping(self, p):
        #self.spec_parameters.calc_ispeed = 1
        if p.value is None:
            LOGGER.info('Vibrational damping disabled')
            # no spec parameters updated here since it is done
            # when writing the data input file by setting the
            # sample temperature to 0K and using the debye_model
            return
        if p.value == 'debye_waller':
            self.spec_parameters.calc_idwsph = 0
        elif p.value == 'averaged_tl':
            self.spec_parameters.calc_idwsph = 1
        LOGGER.info('Vibrational damping activated with \'%s\' model', p.value)

    def bind_temperature(self, p):
        self.spec_parameters.calc_t = float(p.value)
        LOGGER.info('Sample temperature set to %s', p.value * p.unit)
        if not self.get_parameter('use_debye_model').value:
            LOGGER.warning("The sample temperature was set, but will be "
                           "ignored since 'use_debye_model' parameter is "
                           "False.")

    def bind_debye_temperature(self, p):
        self.spec_parameters.calc_td = float(p.value)
        LOGGER.info('Sample Debye temperature set to %s', p.value * p.unit)
        if not self.get_parameter('use_debye_model').value:
            LOGGER.warning("The sample Debye temperature was set, but will be "
                           "ignored since 'use_debye_model' parameter is "
                           "False.")

    def bind_use_debye_model(self, p):
        if p.value:
            self.spec_parameters.calc_idcm = 1
        else:
            self.spec_parameters.calc_idcm = 0
        LOGGER.info('use of the Debye model for mean square displacements: %s',
                    p.value)
        if self.get_parameter('vibrational_damping').value is None:
            LOGGER.warning("'use_debye_model' parameter was set but will be "
                           "ignored as long as vibrational_damping is "
                           "disbaled.")

    def bind_vibration_scaling(self, p):
        self.spec_parameters.calc_rsj = float(p.value)
        LOGGER.info('Vibration scaling set to: %s', p.value)

    def bind_mean_free_path(self, p):
        if isinstance(p.value, str):
            if p.value == 'mono':
                self.spec_parameters.calc_ilpm = 1
            elif p.value == 'SeahDench':
                self.spec_parameters.calc_ilpm = 2
        else:
            if p.value <= 0:
                self.spec_parameters.calc_ilpm = -1
            else:
                self.spec_parameters.calc_ilpm = 0
                self.spec_parameters.calc_xlpm0 = float(p.value)
        LOGGER.info('Mean free path set to: %s', str(p.value))


class PEDParameters(BaseParameters):
    def __init__(self, phagen_parameters, spec_parameters):
        parameters = (
            Parameter('level', types=str, pattern=r'\d+[spdfgSPDFG](\d/2)?$',
                      default='1s', doc="""
            The level is the electronic level where the electron comes from.
            It is written: *nlJ*
            where:

            - *n* is the principal quantum number
            - *l* is the orbital quantum number
            - *J* is the spin-orbit component

            Example::

            >>> calc.spectroscopy_parameters.level = '2p3/2'
            >>> calc.spectroscopy_parameters.level = '2p' # is equivalent to '2p1/2'

            """),
            Parameter('final_state', types=int, limits=(-1, 2), default=2),
            Parameter('spin_orbit', types=(type(None), str),
                      allowed_values=(None, 'single', 'both'), default=None),
        )
        BaseParameters.__init__(self)
        self.add_parameters(*parameters)
        self.phagen_parameters = phagen_parameters
        self.spec_parameters = spec_parameters

    def bind_level(self, p):
        edge = get_level_from_electron_configuration(p.value)
        self.phagen_parameters.edge = edge

        li, so = re.match(r'(^\d+[spdfg])(.*$)', p.value).groups()
        if so == '':
            so = '1/2'

        self.spec_parameters.ped_li = li
        self.spec_parameters.ped_so = so
        self.spec_parameters.extra_level = p.value

    def bind_final_state(self, p):
        self.spec_parameters.ped_initl = p.value

    def bind_spin_orbit(self, p):
        somap = {
            None: 0,
            'single': 1,
            'both': 2}
        self.spec_parameters.ped_iso = somap[p.value]


class EIGParameters(BaseParameters):
    def __init__(self, phagen_parameters, spec_parameters):
        parameters = (
            Parameter('level', types=str, pattern=r'\d+[spdfgSPDFG](\d/2)?$',
                      default='1s', doc="""
            The level is the electronic level where the electron comes from.
            It is written: *nlJ*
            where:

            - *n* is the principal quantum number
            - *l* is the orbital quantum number
            - *J* is the spin-orbit component

            Example::

            >>> calc.spectroscopy_parameters.level = '2p3/2'
            >>> calc.spectroscopy_parameters.level = '2p' # is equivalent to '2p1/2'

            """),
            Parameter('final_state', types=int, limits=(-1, 2), default=2),
            Parameter('spin_orbit', types=(type(None), str),
                      allowed_values=(None, 'single', 'both'), default=None),
            Parameter('kernel_matrix_spectrum', types=(bool,), default=False, doc="""
            Whether to output the kernel matrix spectrum for each energy point.
            """),
            Parameter('method', types=(str,), default='EPSI',
                      allowed_values=['AITK', 'RICH', 'SALZ', 'EPSI', 'EPSG',
                                      'RHOA', 'THET', 'LEGE', 'CHEB', 'OVER',
                                      'DURB', 'DLEV', 'TLEV', 'ULEV', 'VLEV',
                                      'ELEV', 'EULE', 'GBWT', 'VARI', 'ITHE',
                                      'EALG'],
                      doc="""The convergence acceleration scheme to be used.
            """),
        )
        BaseParameters.__init__(self)
        self.add_parameters(*parameters)
        self.phagen_parameters = phagen_parameters
        self.spec_parameters = spec_parameters

    def bind_level(self, p):
        edge = get_level_from_electron_configuration(p.value)
        self.phagen_parameters.edge = edge

        li, so = re.match(r'(^\d+[spdfg])(.*$)', p.value).groups()
        if so == '':
            so = '1/2'

        self.spec_parameters.ped_li = li
        self.spec_parameters.ped_so = so
        self.spec_parameters.extra_level = p.value

    def bind_method(self, p):
        self.spec_parameters.eigval_method = p.value

    def bind_kernel_matrix_spectrum(self, p):
        value = int(p.value)
        self.spec_parameters.eigval_ispectrum_ne = value


class CompCurveParameters(BaseParameters):
    def __init__(self):
        parameters = (
            Parameter('general_npar', types=int, default=1, limits=[1, 2],
                      fmt='d'),
            Parameter('general_norm', types=int, default=0, limits=[0, 5],
                      fmt='d'),
            Parameter('general_iscale', types=int, default=0, limits=[0, 1],
                      fmt='d'),
            Parameter('general_inorm', types=int, default=1, limits=[-2, 2],
                      fmt='d'),
            Parameter('general_isym', types=int, default=0, limits=[0, 2],
                      fmt='d'),
            Parameter('general_sym', types=(int, float), default=180.,
                      limits=[0., 360.], fmt='.2f'),
            Parameter('general_iposi', types=int, default=0, limits=[0, 1],
                      fmt='d'),
            Parameter('general_idist', types=int, default=0, limits=[0, 4],
                      fmt='d'),
            Parameter('general_icur', types=int, default=0, limits=[0, 4],
                      fmt='d'),
            Parameter('general_isa', types=int, default=0, limits=[0, 4],
                      fmt='d'),
            Parameter('general_iprint', types=int, default=1, limits=[0, 1],
                      fmt='d'),

            Parameter('weights_iweight', types=int, default=0, limits=[0, 8],
                      fmt='d'),
            Parameter('weights_alpha', types=(int, float), default=1.,
                      fmt='.2f'),
            Parameter('weights_beta', types=(int, float), default=0.33,
                      fmt='.2f'),
            Parameter('weights_sigma', types=(int, float), default=0.5,
                      fmt='.2f'),
            Parameter('weights_ishift', types=int, default=0, limits=[0, 1],
                      fmt='d'),
            Parameter('weights_maxw', types=(int, float), default=30,
                      fmt='.2f'),

            Parameter('rfc_vi', types=(int, float), default=12, fmt='.2f'),

            Parameter('sim_alphas', types=(int, float), default=1., fmt='.2f'),
            Parameter('sim_betas', types=(int, float), default=1., fmt='.2f'),
            Parameter('sim_nbins', types=int, default=30, limits=[1, None],
                      fmt='d'),

            Parameter('dist_alphad', types=(int, float), default=.5,
                      fmt='.2f'),
            Parameter('dist_ibeta', types=int, default=2, fmt='d'),
            Parameter('dist_l', types=int, default=2, fmt='d'),
            Parameter('dist_sigmad', types=(int, float), default=1, fmt='.2f'),

            Parameter('gof_nbing', types=int, default=30, limits=[1, None],
                      fmt='d'),
            Parameter('gof_alphag', types=(int, float), default=1.5, fmt='.2f'),

            Parameter('kdist_alphak', types=(int, float), default=.5,
                      fmt='.2f'),
            Parameter('kdist_l', types=int, default=2, fmt='d'),
            Parameter('kdist_sigmak', types=(int, float), default=5.5,
                      fmt='.2f'),

            Parameter('mom_ngrid', types=int, default=75, fmt='d'),
            Parameter('mom_nmom', types=int, default=75, fmt='d'),
            Parameter('mom_basis', types=str,
                      allowed_values=['GEOM', 'LEGE', 'CHEB', 'KRAW', 'HAHN',
                                      'MEIX', 'CHAR', 'SHMA'],
                      default='KRAW', fmt='s'),
            Parameter('mom_ialg', types=int, limits=[1, 3], default=1, fmt='d'),
            Parameter('mom_mu', types=(int, float), default=.5, fmt='.2f'),
            Parameter('mom_nu', types=(int, float), default=.5, fmt='.2f'),

            Parameter('chords_ichord', types=int, default=3, limits=[1, 3],
                      fmt='d'),
            Parameter('chords_method', types=str,
                      allowed_values=['SIN', 'HIS', 'SUM'], default='SUM',
                      fmt='s'),
            Parameter('chords_value', types=int, limits=[1, 3], default=1,
                      fmt='d'),
            Parameter('chords_nbinc', types=int, default=30, limits=[1, None],
                      fmt='d'),

            Parameter('codes_nconnect', types=int, allowed_values=[3, 5, 9],
                      default=9, fmt='d'),
            Parameter('codes_scalec', types=(int, float), default=1,
                      fmt='.2f'),

            Parameter('cont_nbin', types=int, limits=[1, None], default=66,
                      fmt='d'),
            Parameter('cont_nlen', types=int, limits=[1, None], default=4,
                      fmt='d'),
            Parameter('cont_shan', types=str,
                      allowed_values=['CDIS', 'TANG', 'CURV', 'TRAR', 'BEAS',
                                      '8CCH', 'CLEN', 'CANG', 'ACDI', 'FOUR'],
                      default='TANG',
                      fmt='s'),
            Parameter('cont_ifou', types=int, limits=[1, 4], default=1,
                      fmt='d'),
            Parameter('cont_inorm', types=int, limits=[1, 4], default=2,
                      fmt='d'),

            Parameter('exp_filename', types=str, default='exp/experiment.txt',
                      fmt='s'),
            Parameter('calc_filename', types=(list,),
                      default=['calc/calculation.txt',],
                      fmt='s'),
            Parameter('calc_param1', types=(float, int), default=0.,
                      fmt='.2f'),
            Parameter('calc_param2', types=(float, int), default=0.,
                      fmt='.2f'),
        )
        BaseParameters.__init__(self)
        self.add_parameters(*parameters)
        self.freeze()

class CompCurveGeneralParameters(BaseParameters):
    def __init__(self, compcurve_parameters):
        parameters = (
            Parameter('normalization', types=(type(None), str),
                      allowed_values=(None, "variance", "area", "max",
                                      "decimal_scaling", "zero_one"),
                      default="max"),
            Parameter('rescale', types=bool, default=True),
            Parameter('function', types=str, allowed_values=(
                "coordinates", "chi", "cdf", "curvature", "signature"),
                default="chi"),
        )
        BaseParameters.__init__(self)
        self.add_parameters(*parameters)
        self.compcurve_parameters = compcurve_parameters
        self.freeze()

    def bind_normalization(self, p):
        value = p.allowed_values.index(p.value)
        self.compcurve_parameters.general_norm = value
        LOGGER.info("Curve Comparison: Normalization mode set to "
                    "\"{}\"".format(p.value))

    def bind_rescale(self, p):
        self.compcurve_parameters.general_iscale = int(p.value)
        state = "deactivated"
        if p.value:
            state = "activated"
        LOGGER.info("Curve Comparison: Rescaling of data {}".format(state))

    def bind_function(self, p):
        value = p.allowed_values.index(p.value)
        self.compcurve_parameters.general_icur = value
        LOGGER.info("Curve Comparison: Type of data used for comparison "
                    "set to \"{}\"".format(p.value))
