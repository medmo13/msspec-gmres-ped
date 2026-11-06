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
# Source file  : src/msspec/misc.py
# Last modified: Mon, 27 Sep 2021 17:49:48 +0200
# Committed by : sylvain tricot <sylvain.tricot@univ-rennes1.fr>


"""
Module misc
===========

"""


import logging
from pint import UnitRegistry
import numpy as np
import inspect
import re
import os


class XRaySource(object):
    MG_KALPHA = 1253.6
    AL_KALPHA = 1486.6
    def __init__(self):
        pass

UREG = UnitRegistry()
#UREG.define('rydberg = c * h * rydberg_constant = Ry')
#UREG.define('bohr_radius = 4 * pi * epsilon_0 * hbar**2 / electron_mass / e**2 = a0')

#logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger('msspec')

np.set_printoptions(formatter={'float': lambda x:'%.2f' % x}, threshold=5)

def set_log_level(level):
    lvl = getattr(logging, level.upper())
    LOGGER.setLevel(lvl)

def set_log_output(stream):
    LOGGER.parent.handlers[0].stream = stream

def get_call_info(frame):
    args, varargs, varkw, loc = inspect.getargvalues(frame)
    _, _, function, _, _ = inspect.getframeinfo(frame)
    s = '%s called with:\n' % function

    for kws in (args, varargs, varkw):
        if kws != None:
            if isinstance(kws, (tuple, list)):
                for kw in kws:
                    s += '\t\t%s = %r\n' % (kw, loc[kw])
            else:
                s += '\t\t%s = %r\n' % (kws, loc[kws])
    return s.strip()


def log_process_output(process, logger=None, severity=logging.INFO):
    if logger == None:
        logger = logging
    else:
        logger = logging.getLogger(logger)
        logger.setLevel(LOGGER.getEffectiveLevel())

    for line in iter(process.stdout.readline, b''):
        logger.log(severity, line.rstrip('\n'))

    process.stdout.close()
    process.wait()


def get_level_from_electron_configuration(notation):
    l_letters = 'spdfghiklmnoqrtuv'
    n_letters = 'klmnopqrstuv'
    pattern = re.compile(r'(?P<n>\d+)(?P<l>[%s%s])'
                         r'((?P<j>\d)/2)?' % (l_letters, l_letters.upper()))
    m = pattern.match(notation)
    assert m, 'Not a valid notation'
    n, l, _, j = m.groups()
    n = int(n)
    l = l_letters.index(l.lower())
    assert (l < n), 'Invalid l'
    j1 = abs(l - 0.5)
    j2 = abs(l + 0.5)
    if j:
        j = int(j) / 2.
    else:
        j = j1
    assert j in (j1, j2), 'Invalid j'
    letter = n_letters[n - 1]
    subscript = str(int((2 * (l + j) + 1) / 2))
    if letter == 'k':
        subscript = ''
    return '{}{}'.format(letter, subscript)
