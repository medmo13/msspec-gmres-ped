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
# Source file  : src/msspec/tests.py
# Last modified: Mon, 27 Sep 2021 17:49:48 +0200
# Committed by : sylvain tricot <sylvain.tricot@univ-rennes1.fr>


import os
import sys
import unittest
from hashlib import md5
from pickle import dumps

from ase.build  import bulk

from msspec.calculator import MSSPEC
from msspec.misc import set_log_level
from msspec.utils import *

set_log_level('error')
RESULTS_FILENAME = os.path.join(os.path.dirname(__file__), 'results.txt')


def perform_test(obj, funcname, filename=RESULTS_FILENAME):
    f = getattr(obj, '_' + funcname)
    output = md5(dumps(f())).hexdigest()
    results = {}
    with open(filename, 'r') as fd:
        for line in fd:
            k, v = line.split()
            results[k] = v
    k = '{}.{}'.format(obj.__class__.__name__, funcname)
    obj.assertEqual(output, results[k])


class UtilsTestCase(unittest.TestCase):
    def setUp(self):
        self.a0 = 3.6
        self.cluster = bulk('Cu', a=self.a0, cubic=True)
        self.cluster = self.cluster.repeat((4, 4, 4))
        center_cluster(self.cluster)

    def _test_base_cluster(self):
        output = self.cluster
        output = output.get_positions()
        return output

    def _test_cut_sphere(self):
        output = cut_sphere(self.cluster, radius=self.a0 + .01)
        output = output.get_positions()
        return output

    def _test_cut_plane(self):
        output = cut_plane(self.cluster, z=0.1)
        output = output.get_positions()
        return output


    def test_base_cluster(self):
        return perform_test(self, 'test_base_cluster')

    def test_cut_sphere(self):
        return perform_test(self, 'test_cut_sphere')

    def test_cut_plane(self):
        return perform_test(self, 'test_cut_plane')

    def runTest(self):
        pass

class PEDCopperTestCase(unittest.TestCase):
    def setUp(self):
        a0 = 3.6 # The lattice parameter in angstroms

        # Create the copper cubic cell
        cluster = bulk('Cu', a=a0, cubic=True)
        # Repeat the cell many times along x, y and z
        cluster = cluster.repeat((4, 4, 4))
        # Put the center of the structure at the origin
        center_cluster(cluster)
        # Keep atoms inside a given radius
        cluster = cut_sphere(cluster, radius=a0 + .01)
        # Keep only atoms below the plane z <= 0
        cluster = cut_plane(cluster, z=0.1)

        # Set the absorber (the deepest atom centered in the xy-plane)
        cluster.absorber = get_atom_index(cluster, 0, 0, -a0)
        # Create a calculator for the PhotoElectron Diffration
        self.calc = MSSPEC(spectroscopy='PED', txt=None)
        # Set the cluster to use for the calculation
        self.calc.set_atoms(cluster)

    def _test_theta_scan(self):
        # Run the calculation
        data = self.calc.get_theta_scan(level='2p3/2')
        return data[-1].cross_section
        #data.save('theta.hdf5')
        #return data[-1]

    def _test_phi_scan(self):
        # Run the calculation
        data = self.calc.get_phi_scan(level='2p3/2', theta=35.)
        return data[-1].cross_section

    def test_theta_scan(self):
        perform_test(self, 'test_theta_scan')

    def test_phi_scan(self):
        perform_test(self, 'test_phi_scan')

    def runTest(self):
        pass

    def tearDown(self):
        self.calc.shutdown()


def run_tests():
    suite = unittest.TestSuite()

    suite.addTest(UtilsTestCase('test_base_cluster'))
    suite.addTest(UtilsTestCase('test_cut_sphere'))
    suite.addTest(UtilsTestCase('test_cut_plane'))

    suite.addTest(PEDCopperTestCase('test_phi_scan'))
    suite.addTest(PEDCopperTestCase('test_theta_scan'))
    runner = unittest.TextTestRunner(verbosity=2)
    rc = runner.run(suite)
    #print("rc = ", rc)
    exit(not(rc.wasSuccessful()))

def create_tests_results():
    with open(RESULTS_FILENAME, 'w') as fd:
        pass

    def create_results(obj, *funcs, **kwargs):
        obj.setUp()
        with open(RESULTS_FILENAME, 'a') as fd:
            for funcname in funcs:
                clsname = obj.__class__.__name__
                print(('Generating results for test: {}.{}'.format(clsname, funcname)))
                o = getattr(obj, '_' + funcname)()
                k = '{}.{}'.format(obj.__class__.__name__, funcname)
                fd.write('{} {}\n'.format(k, md5(dumps(o)).hexdigest()))
        obj.tearDown()

    tc = UtilsTestCase()
    create_results(tc, 'test_base_cluster')
    create_results(tc, 'test_cut_sphere')
    create_results(tc, 'test_cut_plane')

    tc = PEDCopperTestCase()
    create_results(tc, 'test_theta_scan')
    tc = PEDCopperTestCase()
    create_results(tc, 'test_phi_scan')

def delete_results_file():
    os.remove(RESULTS_FILENAME)


if __name__ == "__main__":
    run_tests()
