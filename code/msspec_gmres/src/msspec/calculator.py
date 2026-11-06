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
# Source file  : src/msspec/calculator.py
# Last modified: Tue, 25 Oct 2022 16:21:38 +0200
# Committed by : Sylvain Tricot <sylvain.tricot@univ-rennes1.fr> 1666707698 +0200


"""
Module calculator
=================

This module contains different classes used to define a new calculator for
specific spectroscopies understood by MsSpec.

These spectroscopies are listed :ref:`here <globalparameters-spectroscopy>`.

There is one *calculator* class for each spectroscopy. The class name is based
on the spectroscopy name. For instance, the class for PhotoElectron Diffraction
is called :py:class:`_PED`.

The helper function :py:func:`calculator.MSSPEC` is used to create objects from
these classes by passing the kind of spectroscopy as a keyword argument.

For more information on MsSpec, follow this
`link <https://ipr.univ-rennes1.fr/msspec>`__

"""

import inspect
import logging
import os
import re
import shutil
import sys
import time
from collections import OrderedDict
from datetime import datetime
from io import StringIO
from shutil import copyfile
from shutil import rmtree
from subprocess import PIPE
from subprocess import Popen

import ase.atom
import ase.atoms
import ase.data
from ase.io.extxyz import write_xyz
import numpy as np
from ase.calculators.calculator import Calculator
from terminaltables.ascii_table import AsciiTable

from msspec import iodata
from msspec.calcio import CompCurveIO
from msspec.calcio import PhagenIO
from msspec.calcio import SpecIO
from msspec.data import electron_be
from msspec.misc import get_call_info
from msspec.misc import get_level_from_electron_configuration
from msspec.misc import log_process_output
from msspec.misc import LOGGER
from msspec.misc import set_log_output
from msspec.misc import UREG
from msspec.misc import XRaySource
from msspec.parameters import CalculationParameters
from msspec.parameters import CompCurveParameters
from msspec.parameters import CompCurveGeneralParameters
from msspec.parameters import DetectorParameters
from msspec.parameters import EIGParameters
from msspec.parameters import GlobalParameters
from msspec.parameters import MuffintinParameters
from msspec.parameters import PEDParameters
from msspec.parameters import PhagenMallocParameters
from msspec.parameters import PhagenParameters
from msspec.parameters import ScanParameters
from msspec.parameters import SourceParameters
from msspec.parameters import SpecMallocParameters
from msspec.parameters import SpecParameters
from msspec.parameters import TMatrixParameters
from msspec.phagen.fortran.libphagen import main as do_phagen
from msspec.spec.fortran import _eig_mi
from msspec.spec.fortran import _eig_pw
from msspec.spec.fortran import _phd_mi_noso_nosp_nosym
from msspec.spec.fortran import _phd_gm_noso_nosp_nosym
from msspec.spec.fortran import _phd_se_noso_nosp_nosym
from msspec.spec.fortran import _phd_ce_noso_nosp_nosym
from msspec.spec.fortran import _comp_curves
from msspec.utils import get_atom_index


def _call_with_environment(function, args, values):
    """Call an in-process Fortran entry point with scoped environment flags."""
    previous = {key: os.environ.get(key) for key in values}
    os.environ.update(values)
    try:
        return function(*args)
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


class _MSCALCULATOR(Calculator):
    """
    This class defines an ASE calculator for doing Multiple scattering
    calculations.
    """
    implemented_properties = ['', ]
    __data = {}

    def __init__(self, spectroscopy='PED', algorithm='expansion',
                 polarization=None, dichroism=None, spinpol=False,
                 folder='./calc', txt='-', **kwargs):
        stdout = sys.stdout
        if isinstance(txt, str) and txt != '-':
                stdout = open(txt, 'w')
        #elif isinstance(txt, buffer):
        #    stdout = txt
        elif txt == None:
            stdout = open('/dev/null', 'a')
        #set_log_output(stdout)
        ########################################################################
        LOGGER.debug('Initialization of %s', self.__class__.__name__)
        LOGGER.debug(get_call_info(inspect.currentframe()))
        ########################################################################
        # Init the upper class
        Calculator.__init__(self, **kwargs)

        ########################################################################
        LOGGER.debug('    create low level parameters')
        ########################################################################
        self.phagen_parameters = PhagenParameters()
        self.phagen_malloc_parameters = PhagenMallocParameters()
        self.spec_parameters = SpecParameters()
        self.spec_malloc_parameters = SpecMallocParameters()

        ########################################################################
        LOGGER.debug('    create higher level parameters')
        ########################################################################
        self.tmatrix_parameters = TMatrixParameters(self.phagen_parameters)
        self.muffintin_parameters = MuffintinParameters(self.phagen_parameters,
                                                        self.spec_parameters)


        self.global_parameters = GlobalParameters(self.phagen_parameters,
                                                  self.spec_parameters)

        if spectroscopy == 'PED':
            self.spectroscopy_parameters = PEDParameters(self.phagen_parameters,
                                                         self.spec_parameters)
        elif spectroscopy == 'EIG':
            self.spectroscopy_parameters = EIGParameters(self.phagen_parameters,
                                                         self.spec_parameters)
        else:
            raise NameError('No such spectroscopy')

        self.source_parameters = SourceParameters(self.global_parameters,
                                                  self.phagen_parameters,
                                                  self.spec_parameters)

        self.detector_parameters = DetectorParameters(self.global_parameters,
                                                      self.phagen_parameters,
                                                      self.spec_parameters)

        self.scan_parameters = ScanParameters(self.global_parameters,
                                              self.phagen_parameters,
                                              self.spec_parameters)

        self.calculation_parameters = CalculationParameters(
            self.global_parameters, self.phagen_parameters, self.spec_parameters)

        # initialize all parameters with defaults values
        LOGGER.info("Set default values =========================================")
        for p in (list(self.global_parameters) +
                  list(self.muffintin_parameters) +
                  list(self.tmatrix_parameters) +
                  list(self.source_parameters) +
                  list(self.detector_parameters) +
                  list(self.scan_parameters) +
                  list(self.calculation_parameters) +
                  list(self.spectroscopy_parameters)):
            p.set(p.default)
        LOGGER.info("End of default values ======================================")

        # updated global parameters with provided keywords
        self.global_parameters.spectroscopy = spectroscopy
        self.global_parameters.algorithm = algorithm
        if algorithm == 'gmres' and spectroscopy != 'PED':
            raise ValueError("algorithm='gmres' is only available for "
                             "spectroscopy='PED'")
        self.global_parameters.polarization = polarization
        self.global_parameters.dichroism = dichroism
        self.global_parameters.spinpol = spinpol
        self.global_parameters.folder = folder



        self.phagenio = PhagenIO(self.phagen_parameters,
                                 self.phagen_malloc_parameters)
        self.specio = SpecIO(self.spec_parameters,
                             self.spec_malloc_parameters,
                             self.phagenio)

        ########################################################################
        LOGGER.debug('    create a space dedicated to the calculation')
        ########################################################################
        self.init_folder = os.getcwd()
        #self.msspec_folder = os.path.join(MSSPEC_ROOT)
        self.tmp_folder = os.path.abspath(folder)
        LOGGER.debug('    folder: \'%s\'', self.tmp_folder)
        if not os.path.exists(self.tmp_folder):
            os.makedirs(self.tmp_folder)
            os.makedirs(os.path.join(self.tmp_folder, 'input'))
            os.makedirs(os.path.join(self.tmp_folder, 'output'))
        #copyfile(os.path.join(self.msspec_folder, 'ase', 'Makefile'),
        #         os.path.join(self.tmp_folder, 'Makefile'))

        #os.chdir(self.tmp_folder)

        #inv = cor = 'NO'
        #if algorithm == 'expansion':
        #    pass
        #elif algorithm == 'inversion':
        #    inv = 'YES'
        #elif algorithm == 'correlation':
        #    cor = 'YES'

        # spin orbit resolved (not yet)
        #sorb = 'NO'

        # spin resolved
        #dichro_spinpol = False
        #if dichroism in ('sum_over_spin', 'spin_resolved'):
            dichro_spinpol = True

        #spin = 'NO'
        #if spinpol or dichro_spinpol:
        #    spin = 'YES'

        #if spin == 'YES':
        #    LOGGER.error('Option not implemented!')
        #    raise NotImplementedError(
        #        'Spin polarization is not implemeted yet!')


        calctype_spectro = self.spec_parameters.get_parameter('calctype_spectro')
        calctype_spectro = calctype_spectro.value
        #self._make_opts = (MSSPEC_ROOT, calctype_spectro, inv, cor,
        #                   spin, sorb, self.tmp_folder)

        # Initialize the working environment
        #self._make('init')

        self.modified = False

        self.resources = {}
        ########################################################################
        LOGGER.debug('    initialization done.\n')
        ########################################################################

    def _make(self, target):
        LOGGER.debug(get_call_info(inspect.currentframe()))
        os.chdir(self.tmp_folder)
        cmd = ("make__SPACE__ROOT_FOLDER=\"{}\"__SPACE__SPEC=\"{}\"__SPACE__INV=\"{}\"__SPACE__COR=\"{"
               "}\"__SPACE__"
               "SPIN=\"{}\"__SPACE__SO=\"{}\"__SPACE__CALC_FOLDER=\"{}\"__SPACE__{}").format(*(self._make_opts + (
            target,))).split('__SPACE__')
        #cmd = cmd.replace(' ', '\ ')
        #cmd = cmd.split('__SPACE__')

        LOGGER.debug('    the full command is: %s', cmd)

        child = Popen(cmd,stdout=PIPE, stderr=PIPE)
        logger_name = LOGGER.name
        if target == 'tmatrix':
            logger_name = 'Phagen'
        elif target == 'compute':
            logger_name = 'Spec'

        log_process_output(child, logger=logger_name)
        os.chdir(self.init_folder)

        if child.returncode != 0:
            LOGGER.error("Unable to successfully run the target: {}".format(target))
            sys.exit(1)


    def _guess_ke(self, level):
        """ Try to guess the kinetic energy based on the level and
            the source energy. If the kinetic energy cannot be infered
            because the level is not reported in the database, the returned
            value is None.
        """
        try:
            state = get_level_from_electron_configuration(level)
            absorber_atomic_number = self.atoms[self.atoms.absorber].number
            lines = electron_be[absorber_atomic_number]
            binding_energy = lines[state]
        except KeyError:
            # unable to find a binding energy in the database
            return None

        # let's assume work function energy (in eV)
        wf = 4.5
        source_energy = self.source_parameters.get_parameter('energy').value
        ke = source_energy - binding_energy - wf
        #return np.array(ke, dtype=float).flatten()
        return ke


    def run_phagen(self):
        input_fname = os.path.join(self.tmp_folder, 'input/input.ms')
        modified = self.phagenio.write_input_file(filename=input_fname)

        # Write the external potential file
        pot = self.tmatrix_parameters.potential
        if pot != 'muffin_tin':
            inpot_fname = os.path.join(self.tmp_folder, 'output/fort.2')
            mflag = pot.write(inpot_fname, self.phagenio.prototypical_atoms)
            modified = modified or mflag
            # Make a copy of the input potential file
            dest = os.path.join(self.tmp_folder, 'input/extpot.dat')
            shutil.copy(inpot_fname, dest)

        # Run if input files changed
        os.chdir(os.path.join(self.tmp_folder, 'output'))
        modified = True
        if modified:
            # Run phagen
            do_phagen()
            # Copy some output files to be more explicit
            shutil.copy('clus/clus.out', 'cluster.clu')
            shutil.copy('fort.35', 'tmatrix.tl')
            shutil.copy('fort.55', 'tmatrix.rad')

        # Load data
        self.phagenio.load_tl_file('tmatrix.tl')
        self.phagenio.load_cluster_file('cluster.clu')
        self.phagenio.load_prototypical_atoms('div/molinpot3.out')
        # Change back to the init_folder
        os.chdir(self.init_folder)


    def run_spec(self, malloc={}):
        def get_li(level):
            orbitals = 'spdfghi'
            m = re.match(r'\d(?P<l>[%s])(\d/2)?' % orbitals, level)
            return orbitals.index(m.group('l'))

        #include_fname = os.path.join(self.tmp_folder, 'src/spec.inc')
        input_fname = os.path.join(self.tmp_folder, 'input/spec.dat')
        kdirs_fname = os.path.join(self.tmp_folder, 'input/kdirs.dat')

        mod0 = self.specio.write_input_file(filename=input_fname)
        #mod1 = self.specio.write_include_file(filename=include_fname)
        mod2 = self.specio.write_kdirs_file(filename=kdirs_fname)

        #self.modified = self.modified or mod0 or mod1 or mod2
        self.modified = self.modified or mod0 or mod2

        #self._make('tmatrix')
        #self._make('bin/spec')
        #t0 = time.time()
        #self._make('compute')
        #t1 = time.time()
        #self.resources['spec_time'] = t1 - t0
        self.modified = True
        if self.modified:
            #self.get_tmatrix()
            t0 = time.time()
            os.chdir(os.path.join(self.tmp_folder, 'output'))
            # set/get the dimension values
            requirements = OrderedDict({
               'NATP_M'     : self.phagenio.nat,
               'NATCLU_M'   : len(self.atoms),
               'NAT_EQ_M'   : self.phagenio.nateqm,
               'N_CL_L_M'   : 1,
               'NE_M'       : self.phagenio.ne,
               'NL_M'       : self.phagenio.nlmax + 1,
               'LI_M'       : get_li(self.spec_parameters.extra_level) + 1,
               'NEMET_M'    : 1,
               'NO_ST_M'    : self.spec_parameters.calc_no,
               'NDIF_M'     : 10,
               'NSO_M'      : 2,
               'NTEMP_M'    : 1,
               'NODES_EX_M' : 3,
               'NSPIN_M'    : 1, # to change for spin dependent
               'NTH_M'      : 2000,
               'NPH_M'      : 8000,
               'NDIM_M'     : 100000,
               'N_TILT_M'   : 11, # to change see extdir.f
               'N_ORD_M'    : 250,
               'NPATH_M'    : self.spec_malloc_parameters.NPATH_M,
               'NGR_M'      : 10,})

            # update with provided values
            for key, value in malloc.items():
                requirements[key] = value

            # set some automatic values for memory allocation
            for key, value in requirements.items():
                setattr(self.spec_malloc_parameters, key, value)

            # Get the spec function to launch
            if self.global_parameters.spectroscopy == 'PED':
                if self.global_parameters.algorithm == 'expansion':
                    do_spec = _phd_se_noso_nosp_nosym.run
                elif self.global_parameters.algorithm == 'inversion':
                    do_spec = _phd_mi_noso_nosp_nosym.run
                elif self.global_parameters.algorithm == 'gmres':
                    do_spec = _phd_gm_noso_nosp_nosym.run
                elif self.global_parameters.algorithm == 'correlation':
                    do_spec = _phd_ce_noso_nosp_nosym.run
                else:
                    LOGGER.error("\'{}\' spectroscopy with \'{}\' algorithm is not "
                                 "an allowed combination.".format(self.global_parameters.spectroscopy,
                                 self.global_parameters.algorithm))
                    raise ValueError
            elif self.global_parameters.spectroscopy == 'EIG':
                if self.global_parameters.algorithm == 'inversion':
                    do_spec = _eig_mi.run
                elif self.global_parameters.algorithm == 'power':
                    do_spec = _eig_pw.run
                else:
                    LOGGER.error("\'{}\' spectroscopy with \'{}\' algorithm is not "
                                 "an allowed combination.".format(self.global_parameters.spectroscopy,
                                 self.global_parameters.algorithm))
                    raise ValueError
            else:
                LOGGER.error("\'{}\' spectroscopy is not supported yet.".format(self.global_parameters.spectroscopy))
                raise NotImplementedError

            # cannot use this, unfortunately
            #do_spec(*requirements.values())
            spec_args = (
               requirements['NATP_M'],
               requirements['NATCLU_M'],
               requirements['NAT_EQ_M'],
               requirements['N_CL_L_M'],
               requirements['NE_M'],
               requirements['NL_M'],
               requirements['LI_M'],
               requirements['NEMET_M'],
               requirements['NO_ST_M'],
               requirements['NDIF_M'],
               requirements['NSO_M'],
               requirements['NTEMP_M'],
               requirements['NODES_EX_M'],
               requirements['NSPIN_M'],
               requirements['NTH_M'],
               requirements['NPH_M'],
               requirements['NDIM_M'],
               requirements['N_TILT_M'],
               requirements['N_ORD_M'],
               requirements['NPATH_M'],
               requirements['NGR_M'])

            if self.global_parameters.algorithm == 'gmres':
                solver_env = {
                    'MSSPEC_TEST_TREE_SERIES': '1',
                    'MSSPEC_TREE_SERIES_ORDER': str(
                        self.calculation_parameters.linear_solver_max_iterations),
                    'MSSPEC_TREE_TOL': str(
                        self.calculation_parameters.linear_solver_tolerance),
                    'MSSPEC_TREE_GMRES': '1',
                    'MSSPEC_TREE_GMRES_RESTART': str(
                        self.calculation_parameters.linear_solver_restart),
                }
                # Optional settings are written only when they differ from the default,
                # so a value already present in the environment is kept.
                cp = self.calculation_parameters
                if cp.linear_solver_fixed_steps >= 0:
                    solver_env['MSSPEC_TREE_STOP_AFTER'] = str(cp.linear_solver_fixed_steps)
                if not cp.linear_solver_axial_translation:
                    solver_env['MSSPEC_TREE_AXIAL_TRANSLATION'] = '0'
                if cp.linear_solver_threads > 0:
                    solver_env['MSSPEC_TREE_ACTION_THREADS'] = str(cp.linear_solver_threads)
                if cp.linear_solver_verify:
                    solver_env['MSSPEC_TREE_VERIFY'] = '1'
                if cp.linear_solver_diagnostics:
                    solver_env['MSSPEC_TREE_DIAGNOSTICS'] = '1'
                if cp.linear_solver_iterate_dump:
                    solver_env['MSSPEC_TREE_ITER_DUMP'] = '1'
                _call_with_environment(do_spec, spec_args, solver_env)
            else:
                do_spec(*spec_args)

            t1 = time.time()
            self.resources['spec_time'] = t1 - t0
            os.chdir(self.init_folder)

    def get_tmatrix(self):
        LOGGER.info("Getting the TMatrix...")
        LOGGER.debug(get_call_info(inspect.currentframe()))

        # If using an external potential, Phagen should be run twice,
        # the first time with the internal MT potential.
        pot = self.tmatrix_parameters.potential
        if pot != 'muffin_tin':
            LOGGER.info("Phagen is running first with MuffinTin potential...")
            self.tmatrix_parameters.potential = "muffin_tin"
            self.tmatrix_parameters.exchange_correlation = "hedin_lundqvist_complex"
            self.run_phagen()
            self.tmatrix_parameters.potential = pot

        self.run_phagen()

        tl = self.phagenio.tl
        tl_threshold = self.tmatrix_parameters.get_parameter('tl_threshold')
        if tl_threshold.value != None:
            LOGGER.debug("    applying tl_threshold to %s...",
                         tl_threshold.value)
            go_on = True
            while go_on:
                go_on = False
                for ia in range(self.phagenio.nat):
                    for ie in range(self.phagenio.ne):
                        last_tl = tl[ia][ie][-1, -2:]
                        # convert to complex
                        last_tl = last_tl[0] + last_tl[1] * 1j
                        if np.abs(last_tl) < tl_threshold.value:
                            # remove last line of tl
                            tl[ia][ie] = tl[ia][ie][:-1, :]
                            go_on = True

        max_tl = self.tmatrix_parameters.get_parameter('max_tl').value
        cluster = self.phagen_parameters.get_parameter('atoms').value
        #proto_indices = cluster.get_array('proto_indices')
        proto_indices = np.array(
            [atom.get('proto_index') for atom in self.atoms])

        if max_tl != None:
            LOGGER.debug("    applying max_tl: %s", max_tl)
            for ia in range(self.phagenio.nat):
                for ie in range(self.phagenio.ne):
                    try:
                        # for each set of tl:
                        # 1. get the symbol of the prototipical atom
                        j = np.where(proto_indices == ia+1)[0]
                        symbol = cluster[j][0].symbol
                        # 2. get the number of max tl allowed
                        ntl = max_tl[symbol]
                        # 3. reshape the tl set accordingly
                        tl[ia][ie] = tl[ia][ie][:ntl, :]
                    except KeyError:
                        pass

        self.phagenio.write_tl_file(
            os.path.join(self.tmp_folder, 'output/tmatrix.tl'))

        # update spec extra parameters here
        self.spec_parameters.set_parameter('extra_nat', self.phagenio.nat)
        self.spec_parameters.set_parameter('extra_nlmax', self.phagenio.nlmax)



    def set_atoms(self, atoms):
        """Defines the cluster on which the calculator will work.

        :param atoms: The cluster to attach the calculator to.
        :type atoms: :py:class:`ase.Atoms`
        """
        if atoms.absorber == None:
            LOGGER.error("You must define the absorber before setting the atoms to the"
                         "calculator.")
        self.atoms = atoms
        self.phagen_parameters.set_parameter('atoms', atoms)
        self.spec_parameters.set_parameter('extra_atoms', atoms)


    def get_parameters(self):
        """Get all the defined parameters in the calculator.

        :return: A list of all parameters objects.
        :rtype: List of :py:class:`parameters.Parameter`

        """
        _ = []
        for section in ('global', 'muffintin', 'tmatrix', 'spectroscopy',
                        'source', 'detector', 'scan', 'calculation'):
            parameters = getattr(self, section + '_parameters')
            for p in parameters:
                _.append(p)
        return _

    def add_cluster_to_dset(self, dset):
        clusbuf = StringIO()
        self.atoms.info['absorber'] = self.atoms.absorber
        #self.atoms.write(clusbuf, format='xyz')
        write_xyz(clusbuf, self.atoms)
        dset.add_parameter(group='Cluster', name='cluster', value=clusbuf.getvalue(), hidden="True")

    def shutdown(self):
        """Removes the temporary folder and all its content.

        The user may whish to keep the calculation folder (see :ref:`globalparameters-folder`) so it is not removed
        at the end of the calculation. The calculation folder contains raw results from *Phagen* and *Spec* programs as
        well as their input files and configuration. It allows the program to save some time by not repeating some
        tasks (such as the Fortran code generation, building the binaries, computing things that did not changed
        between runs...).
        Calling this function at the end of the calculation will erase this calculation folder.

        .. warning::

            Calling this function will erase the produced data without prompting you for confirmation,
            so take care of explicitly saving your results in your script, by using the
            :py:func:`iodata.Data.save` method for example.

        """
        LOGGER.info('Deleting temporary files...')
        rmtree(self.tmp_folder)

class _PED(_MSCALCULATOR):
    """This class creates a calculator object for PhotoElectron DIffraction
    spectroscopy.

    :param algorithm: The algorithm to use for the computation. See
        :ref:`globalparameters-algorithm` for more details about the allowed
        values and the type.

    :param polarization: The incoming light polarization (see
        :ref:`globalparameters-polarization`)

    :param dichroism: Wether to enable or not the dichroism (see
        :ref:`globalparameters-dichroism`)

    :param spinpol: Enable or disable the spin polarization in the calculation
        (see :ref:`globalparameters-spinpol`)

    :param folder: The path to the temporary folder for the calculations. See
        :ref:`globalparameters-folder`

    :param txt: The name of a file where to redirect standard output. The string
        '-' will redirect the standard output to the screen (default).
    :type txt: str

    .. note::

        This class constructor is not meant to be called directly by the user.
        Use the :py:func:`MSSPEC` to instanciate any calculator.


    """
    def __init__(self, algorithm='expansion', polarization=None, dichroism=None,
                 spinpol=False, folder='./calc', txt='-'):
        _MSCALCULATOR.__init__(self, spectroscopy='PED', algorithm=algorithm,
                              polarization=polarization, dichroism=dichroism,
                              spinpol=spinpol, folder=folder, txt=txt)

        self.iodata = iodata.Data('PED Simulation')

    def _get_scan(self, scan_type='theta', phi=0,
                  theta=np.linspace(-70, 70, 141), level=None,
                  kinetic_energy=None, data=None,
                  malloc={}, other_parameters={}):
        LOGGER.info("Computting the %s scan...", scan_type)
        if data:
            self.iodata = data

        if kinetic_energy is None:
            # try to guess the kinetic energy
            kinetic_energy = self._guess_ke(level)

        # if still None...
        if kinetic_energy is None:
            LOGGER.error('Unable to guess the kinetic energy!')
            raise ValueError('You must define a kinetic_energy value.')

        # update the parameters
        self.scan_parameters.set_parameter('kinetic_energy', kinetic_energy)
        all_ke = self.scan_parameters.get_parameter('ke_array')
        if np.any(all_ke.value < 0):
            LOGGER.error('Source energy is not high enough or level too deep!')
            raise ValueError('Kinetic energy is < 0! ({})'.format(
                kinetic_energy))
        self.scan_parameters.set_parameter('type', scan_type)

        # make sure there is only one energy point in scatf scan
        if scan_type == 'scatf':
            assert len(all_ke) == 1, ('kinetic_energy should not be an array '
                                      'in scatf scan')


        if scan_type != 'scatf':
            self.scan_parameters.set_parameter('phi', phi)
            self.scan_parameters.set_parameter('theta', theta)

        self.spectroscopy_parameters.set_parameter('level', level)

        # It is still possible to modify any option right before runing phagen
        # and spec
        for k, v in other_parameters.items():
            grp_str, param_str = k.split('.')
            grp = getattr(self, grp_str)
            grp.set_parameter(param_str, v, force=True)

        self.get_tmatrix()
        self.run_spec(malloc)

        # Now load the data
        ndset = len(self.iodata)
        dset = self.iodata.add_dset('{} scan [{:d}]'.format(scan_type, ndset))
        for p in self.get_parameters():
            bundle = {'group': str(p.group),
                      'name': str(p.name),
                      'value': str(p.value),
                      'unit': '' if p.unit is None else str(p.unit)}
            dset.add_parameter(**bundle)
        if scan_type in ('theta', 'phi', 'energy'):
            results_fname = os.path.join(self.tmp_folder, 'output/results.dat')
            data = self.specio.load_results(results_fname)
            for _plane, _theta, _phi, _energy, _dirsig, _cs in data.T:
                if _plane == -1:
                    dset.add_row(theta=_theta, phi=_phi, energy=_energy, cross_section=_cs, direct_signal=_dirsig)
        elif scan_type in ('scatf',):
            results_fname = os.path.join(self.tmp_folder, 'output/facdif1.dat')
            data = self.specio.load_facdif(results_fname)
            data = data[:, [1, 4, 5, 6, 8]].T
            _proto, _sf_real, _sf_imag, _theta, _energy = data
            _sf = _sf_real + _sf_imag * 1j
            dset.add_columns(proto_index=_proto, sf_real=np.real(_sf),
                             sf_imag=np.imag(_sf), sf_module=np.abs(_sf),
                             theta=_theta, energy=_energy)
        elif scan_type in ('theta_phi',):
            results_fname = os.path.join(self.tmp_folder, 'output/results.dat')
            data = self.specio.load_results(results_fname)
            #theta_c, phi_c = data[[2, 3], :]
            #xsec_c = data[-1, :]
            #dirsig_c = data[-2, :]

            #dset.add_columns(theta=theta_c)
            #dset.add_columns(phi=phi_c)
            #dset.add_columns(cross_section=xsec_c)
            #dset.add_columns(direct_signal=dirsig_c)
            for _plane, _theta, _phi, _energy, _dirsig, _cs in data.T:
                if _plane == -1:
                    dset.add_row(theta=_theta, phi=_phi, energy=_energy, cross_section=_cs,
                                 direct_signal=_dirsig)

        # create a view
        title = ''
        for ke in all_ke.value:
            if scan_type == 'theta':
                absorber_symbol = self.atoms[self.atoms.absorber].symbol
                title = 'Polar scan of {}({}) at {:.2f} eV'.format(
                    absorber_symbol, level, ke)
                xlabel = r'Angle $\theta$($\degree$)'
                ylabel = r'Signal (a. u.)'

                view = dset.add_view("E = {:.2f} eV".format(ke), title=title,
                                     xlabel=xlabel, ylabel=ylabel)
                for angle_phi in self.scan_parameters.get_parameter(
                    'phi').value:
                    where = ("energy=={:.2f} and phi=={:.2f}"
                             "").format(ke, angle_phi)
                    legend = r'$\phi$ = {:.1f} $\degree$'.format(angle_phi)
                    view.select('theta', 'cross_section', where=where,
                                legend=legend)
            if scan_type == 'phi':
                absorber_symbol = self.atoms[self.atoms.absorber].symbol
                title = 'Azimuthal scan of {}({}) at {:.2f} eV'.format(
                    absorber_symbol, level, ke)
                xlabel = r'Angle $\phi$($\degree$)'
                ylabel = r'Signal (a. u.)'

                view = dset.add_view("E = {:.2f} eV".format(ke), title=title,
                                     xlabel=xlabel, ylabel=ylabel)
                for angle_theta in self.scan_parameters.get_parameter(
                    'theta').value:
                    where = ("energy=={:.2f} and theta=={:.2f}"
                             "").format(ke, angle_theta)
                    legend = r'$\theta$ = {:.1f} $\degree$'.format(angle_theta)
                    view.select('phi', 'cross_section', where=where,
                                legend=legend)

            if scan_type == 'theta_phi':
                absorber_symbol = self.atoms[self.atoms.absorber].symbol
                title = ('Stereographic projection of {}({}) at {:.2f} eV'
                         '').format(absorber_symbol, level, ke)
                xlabel = r'Angle $\phi$($\degree$)'
                ylabel = r'Signal (a. u.)'

                view = dset.add_view("E = {:.2f} eV".format(ke), title=title,
                                     xlabel=xlabel, ylabel=ylabel,
                                     projection='stereo', colorbar=True, autoscale=False)
                view.select('theta', 'phi', 'cross_section')


        if scan_type == 'scatf':
            for i in range(self.phagenio.nat):
                proto_index = i+1
                title = 'Scattering factor at {:.3f} eV'.format(kinetic_energy)

                mini = min(map(np.min, [dset.sf_real, dset.sf_imag, dset.sf_module]))
                maxi = max(map(np.max, [dset.sf_real, dset.sf_imag, dset.sf_module]))
                view = dset.add_view("Proto. atom #{:d}".format(proto_index),
                                     title=title, projection='polar',
                                     ylim=[mini, maxi])
                where = "proto_index=={:d}".format(proto_index)
                view.select('theta', 'sf_module', where=where,
                            legend=r'$|f(\theta)|$')
                view.select('theta', 'sf_real', where=where,
                            legend=r'$\Re(f(\theta))$')
                view.select('theta', 'sf_imag', where=where,
                            legend=r'$\Im(f(\theta))$')

        if scan_type == 'energy':
            absorber_symbol = self.atoms[self.atoms.absorber].symbol
            title = (r'Energy scan of {}({}) at $\theta$={:.2f}$\degree$ and '
                      '$\phi$={:.2f}$\degree$').format(
                     absorber_symbol, level, theta, phi)
            xlabel = r'Photoelectron kinetic energy (eV)'
            ylabel = r'Signal (a. u.)'

            view = dset.add_view("EnergyScan".format(ke), title=title,
                                 xlabel=xlabel, ylabel=ylabel)
            view.select('energy', 'cross_section')

        # save the cluster
        #clusbuf = StringIO()
        #self.atoms.info['absorber'] = self.atoms.absorber
        #self.atoms.write(clusbuf, format='xyz')
        #dset.add_parameter(group='Cluster', name='cluster', value=clusbuf.getvalue(), hidden="True")
        self.add_cluster_to_dset(dset)

        LOGGER.info('%s scan computing done!', scan_type)

        return self.iodata

    def get_potential(self, atom_index=None, data=None, units={'energy': 'eV', 'space': 'angstrom'}):
        """Computes the coulombic part of the atomic potential.

        :param atom_index: The atom indices to get the potential of, either as a list or as a single integer
        :param data: The data object to store the results to
        :param units: The units to be used. A dictionary with the keys 'energy' and 'space'
        :return: A Data object
        """
        LOGGER.info("Getting the Potential...")
        LOGGER.debug(get_call_info(inspect.currentframe()))

        _units = {'energy': 'eV', 'space': 'angstrom'}
        _units.update(units)

        if data:
            self.iodata = data

        self.run_phagen()

        filename = os.path.join(self.tmp_folder, 'output/tmatrix.tl')
        tl = self.phagenio.load_tl_file(filename)

        filename = os.path.join(self.tmp_folder, 'output/cluster.clu')
        self.phagenio.load_cluster_file(filename)

        if self.phagen_parameters.potgen in ('in'):
            filename = os.path.join(self.tmp_folder, 'output/plot/plot_vc.dat')
        else:
            filename = os.path.join(self.tmp_folder, 'output/plot/plot_v.dat')
        pot_data = self.phagenio.load_potential_file(filename)

        cluster = self.phagen_parameters.get_parameter('atoms').value

        dset = self.iodata.add_dset('Potential [{:d}]'.format(len(self.iodata)))
        r = []
        v = []
        index = np.empty((0,1), dtype=int)

        absorber_position = cluster[cluster.absorber].position
        for _pot_data in pot_data:
            # find the proto index of these data
            at_position = (_pot_data['coord'] * UREG.bohr_radius).to('angstrom').magnitude + absorber_position
            at_index = get_atom_index(cluster, *at_position)
            at_proto_index = cluster[at_index].get('proto_index')
            #values = np.asarray(_pot_data['values'])
            values = _pot_data['values']
            index = np.append(index, np.ones(values.shape[0], dtype=int) * at_proto_index)
            r = np.append(r, (values[:, 0] * UREG.bohr_radius).to(_units['space']).magnitude)
            v = np.append(v, (values[:, 1] * UREG.rydberg).to(_units['energy']).magnitude)

        dset.add_columns(distance=r, potential=v, index=index)
        view = dset.add_view('potential data', title='Potential energy of atoms',
                             xlabel='distance from atomic center [{:s}]'.format(_units['space']),
                             ylabel='energy [{:s}]'.format(_units['energy']), scale='linear',
                             autoscale=True)

        if atom_index == None:
            for i in range(pot_data[len(pot_data) - 1]['index']):
                view.select('distance', 'potential', where="index=={:d}".format(i),
                            legend="Atom index #{:d}".format(i + 1))
        else:
            for i in atom_index:
                view.select('distance', 'potential', where="index=={:d}".format(cluster[i].get('proto_index') - 1),
                            legend="Atom index #{:d}".format(i))

        return self.iodata

    def get_scattering_factors(self, level='1s', kinetic_energy=None,
                               data=None, **kwargs):
        """Computes the scattering factors of all prototypical atoms in the
        cluster.

        This function computes the real and imaginery parts of the scattering
        factor as well as its modulus for each non symetrically equivalent atom
        in the cluster. The results are stored in the *data* object if provided
        as a parameter.

        :param level: The electronic level. See :ref:`pedparameters-level`.
        :param kinetic_energy: see :ref:`scanparameters-kinetic_energy`.
        :param data: a :py:class:`iodata.Data` object to append the results to
            or None.

        :returns: The modified :py:class:`iodata.Data` object passed as an
            argument or a new :py:class:`iodata.Data` object.

        """
        data = self._get_scan(scan_type='scatf', level=level, data=data,
                              kinetic_energy=kinetic_energy, **kwargs)
        return data

    def get_theta_scan(self, phi=0, theta=np.linspace(-70, 70, 141),
                       level=None, kinetic_energy=None, data=None, **kwargs):
        """Computes a polar scan of the emitted photoelectrons.

        :param phi: The azimuthal angle in degrees. See
            :ref:`scanparameters-phi`.
        :param theta: All the values of the polar angle to be computed. See
            :ref:`scanparameters-theta`.
        :param level: The electronic level. See :ref:`pedparameters-level`.
        :param kinetic_energy: see :ref:`scanparameters-kinetic_energy`.
        :param data: a :py:class:`iodata.Data` object to append the results to
            or None.

        :returns: The modified :py:class:`iodata.Data` object passed as an
            argument or a new :py:class:`iodata.Data` object.

        """
        data = self._get_scan(scan_type='theta', level=level, theta=theta,
                              phi=phi, kinetic_energy=kinetic_energy,
                              data=data, **kwargs)
        return data

    def get_phi_scan(self, phi=np.linspace(0, 359, 359), theta=0,
                     level=None, kinetic_energy=None, data=None, **kwargs):
        """Computes an azimuthal scan of the emitted photoelectrons.

        :param phi: All the values of the azimuthal angle to be computed. See
            :ref:`scanparameters-phi`.
        :param theta: The polar angle in degrees. See
            :ref:`scanparameters-theta`.
        :param level: The electronic level. See :ref:`pedparameters-level`.
        :param kinetic_energy: see :ref:`scanparameters-kinetic_energy`.
        :param data: a :py:class:`iodata.Data` object to append the results to
            or None.

        :returns: The modified :py:class:`iodata.Data` object passed as an
            argument or a new :py:class:`iodata.Data` object.

        """
        data = self._get_scan(scan_type='phi', level=level, theta=theta,
                              phi=phi, kinetic_energy=kinetic_energy,
                              data=data, **kwargs)
        return data

    def get_theta_phi_scan(self, phi=np.linspace(0, 360),
                           theta=np.linspace(0, 90, 45), level=None,
                           kinetic_energy=None, data=None, **kwargs):
        """Computes a stereographic scan of the emitted photoelectrons.

        The azimuth ranges from 0 to 360° and the polar angle ranges from 0 to
        90°.

        :param level: The electronic level. See :ref:`pedparameters-level`.
        :param kinetic_energy: see :ref:`scanparameters-kinetic_energy`.
        :param data: a :py:class:`iodata.Data` object to append the results to
            or None.

        :returns: The modified :py:class:`iodata.Data` object passed as an
            argument or a new :py:class:`iodata.Data` object.

        """
        data = self._get_scan(scan_type='theta_phi', level=level, theta=theta,
                              phi=phi, kinetic_energy=kinetic_energy, data=data,
                              **kwargs)
        return data

    def get_energy_scan(self, phi=0, theta=0,
                        level=None, kinetic_energy=None, data=None, **kwargs):
        """Computes an energy scan of the emitted photoelectrons.

        :param phi: All the values of the azimuthal angle to be computed. See
            :ref:`scanparameters-phi`.
        :param theta: The polar angle in degrees. See
            :ref:`scanparameters-theta`.
        :param level: The electronic level. See :ref:`pedparameters-level`.
        :param kinetic_energy: see :ref:`scanparameters-kinetic_energy`.
        :param data: a :py:class:`iodata.Data` object to append the results to
            or None.

        :returns: The modified :py:class:`iodata.Data` object passed as an
            argument or a new :py:class:`iodata.Data` object.

        """
        data = self._get_scan(scan_type='energy', level=level, theta=theta,
                              phi=phi, kinetic_energy=kinetic_energy,
                              data=data, **kwargs)
        return data


class _EIG(_MSCALCULATOR):
    """
    .. note::

        This class constructor is not meant to be called directly by the user.
        Use the :py:func:`MSSPEC` to instanciate any calculator.

    """
    def __init__(self, algorithm='inversion', polarization=None, dichroism=None,
                 spinpol=False, folder='./calc', txt='-'):
        _MSCALCULATOR.__init__(self, spectroscopy='EIG', algorithm=algorithm,
                               polarization=polarization, dichroism=dichroism,
                               spinpol=spinpol, folder=folder, txt=txt)
        if algorithm not in ('inversion', 'power'):
            LOGGER.error("Only the 'inversion' or the 'power' algorithms "
                         "are supported in EIG spectroscopy mode")
            exit(1)
        self.iodata = iodata.Data('EIG Simulation')
        # ensure no core-hole is created by using the "led" (LEED) mode
        # of Phagen
        self.phagen_parameters.calctype = "led"


    def get_eigen_values(self, kinetic_energy, data=None):
        LOGGER.info("Computting the eigen values...")
        if data:
            self.iodata = data

        # update the parameters
        self.scan_parameters.set_parameter('kinetic_energy', kinetic_energy)
        all_ke = self.scan_parameters.get_parameter('ke_array').value
        if np.any(all_ke < 0):
            LOGGER.error('Source energy is not high enough or level too deep!')
            raise ValueError('Kinetic energy is < 0! ({})'.format(
                kinetic_energy))

        # ensure that if the alogorithm is 'inversion' the IPWM must be set to 0
        if self.global_parameters.algorithm.lower() == 'inversion':
            self.spec_parameters.eigval_ipwm = 0

        self.get_tmatrix()
        self.run_spec()

        # Now create a dataset
        ndset = len(self.iodata)
        dset = self.iodata.add_dset('Eigen values calculation [{:d}]'.format(ndset))
        for p in self.get_parameters():
            bundle = {'group': str(p.group),
                      'name': str(p.name),
                      'value': str(p.value),
                      'unit': '' if p.unit is None else str(p.unit)}
            dset.add_parameter(**bundle)

        results_fname = os.path.join(self.tmp_folder, 'output/results.dat')
        data = self.specio.load_results(results_fname)

        energy, spec_rad = data
        dset.add_columns(energy=energy, spectral_radius=spec_rad)

        # create a view
        title = 'Spectral radii versus kinetic energy'
        xlabel = r'Kinetic energy (eV)'
        ylabel = r'Spectral radius'
        view = dset.add_view("Spectral radius", title=title, xlabel=xlabel, ylabel=ylabel,
                             marker='o')
        view.select('energy', 'spectral_radius')

        # add the cluster to the dataset
        self.add_cluster_to_dset(dset)


        return self.iodata


def MSSPEC(spectroscopy='PED', **kwargs):
    """ The MsSpec calculator constructor.

    Instanciates a calculator for the given spectroscopy.

    :param spectroscopy: See :ref:`globalparameters-spectroscopy`.
    :param kwargs: Other keywords are passed to the spectroscopy-specific
        constructor.
    :returns: A :py:class:`calculator._MSCALCULATOR` object.
    """
    module = sys.modules[__name__]
    cls = getattr(module, '_' + spectroscopy)
    return cls(**kwargs)


class RFACTOR(object):
    def __init__(self, folder='./rfc'):
        self.iodata = iodata.Data('R-Factor analysis')
        self.folder = folder
        self._params = CompCurveParameters()
        self.general_parameters = CompCurveGeneralParameters(
            self._params)
        self.io = CompCurveIO(self._params)

        # prepare the working area
        if not os.path.exists(self.folder):
            os.makedirs(self.folder)
            os.makedirs(os.path.join(self.folder, 'rfc'))
            os.makedirs(os.path.join(self.folder, 'exp/div'))
            os.makedirs(os.path.join(self.folder, 'calc/div'))
            os.makedirs(os.path.join(self.folder, 'plot'))

        # store the r-factor analysis results
        # self.variables = {'variable0_name': [value0, value1, ...],
        #                   'variable1_name': [value0, value1, ...],
        #                   ...}
        self.variables = {}
        # self.rfc = [[rf0_0, rf0_1, ... , rf0_11],    <-- run 0
        #             [rf1_0, rf1_1, ... , rf1_11],    <-- run 1
        #             ............................,
        #             [rfn_0, rfn_1, ... , rfn_11]]    <-- run n
        self.rfc = []

        # The x and y array to compare to
        self.xref = self.yref = [0,]

        # The index of the best value
        self.index = 0

        # The best values as a dictionary
        self.best_values = {}

        # The number of calculation files in the stack. This counter is
        # inremented each time calculation is run
        self.stack_count = 0

        # initialize all parameters with defaults values
        self.logger = logging.getLogger("RFACTOR")
        self.logger.info("Set default values =========================================")
        for p in (list(self.general_parameters)):
            p.set(p.default)
        self.logger.info("End of default values ======================================")
        #exit()


    def set_references(self, X, Y):
        self.xref = X
        self.yref = Y

    def run(self, *args, data=None, **kwargs):
        # Get the data object if provided
        #data = kwargs.pop('data', None)
        if data:
            assert isinstance(data, iodata.Data), ("Unsupported type for data"
                                                   "keyword.")
            self.iodata = data

        # Move to the working folder
        cwd = os.getcwd()
        os.chdir(self.folder)
        # Write the reference data
        np.savetxt('exp/experiment.txt', np.transpose([self.xref, self.yref]))

        # Write all the input calculation files
        # Number of input files
        noif = int(len(args)/2)
        for i in range(noif):
            X, Y = args[2*i], args[2*i+1]
            fname = os.path.join('calc',
                                 'calculation{:d}.txt'.format(self.stack_count))
            # And save to the working space
            np.savetxt(fname, np.transpose([X, Y]))
            self.stack_count += 1

        # Update the list of input calculation files
        self._params.calc_filename = []
        for i in range(self.stack_count):
            fname = os.path.join('calc', 'calculation{:d}.txt'.format(i))
            self._params.calc_filename.append(fname)

        # Write the input file
        self.io.write_input_file('comp_curves.dat')
        # And finally run
        _comp_curves.run()

        #######################################################################
        # Load the results
        #######################################################################
        self.rfc = []
        for i in range(self.stack_count):
            # Read results files and get the R-Factors, exp data and calc
            # data for file index 'i'
            rfc, exp_data, calc_data = self.io.load_results(index=i)
            # Update the list of R-Factors results
            self.rfc.append(rfc)
        # Update the variables values
        for i in range(noif):
            for k,v in kwargs.items():
                try:
                    vi = v[i]
                except (IndexError, TypeError) as err:
                    vi = v
                try:
                    self.variables[k].append(vi)
                except KeyError as err:
                    self.variables[k] = [vi,]

        #######################################################################
        # Analysis
        #######################################################################
        rfc = np.array(self.rfc)
        # Get the index of the minimum for each R-Factor (each column)
        all_min = np.argmin(rfc, axis=0)
        # Iterate over each run and get the number of R-Factors that are min
        # for this set
        all_counts = np.zeros(self.stack_count, dtype=int)
        for i in range(self.stack_count):
            all_counts[i] = len(np.where(all_min==i)[0])

        # The best set of variables (ie the run index) is the one with the
        # highest number of counts:
        self.index = np.argmax(all_counts)

        # Update the "best values" dict
        self.best_values = {k:self.variables[k][self.index] for k in
                            self.variables.keys()}

        # with np.printoptions(precision=6, linewidth=300, threshold=200):
            # print('rfc: ')
            # print(rfc)
            # print('all_min: ', all_min)
            # print('all_counts: ', all_counts)
            # print('variables: ', self.variables)

        #######################################################################
        # Store values
        #######################################################################
        # Three datasets will be created or existing ones will be reused if
        # any.
        dset_values_title = "CurveComparison Values"
        dset_results_title = "CurveComparison Results"
        dset_rfc_title = "CurveComparison R-Factors"

        # Create (or re-create) the datasets
        dset_values = self.iodata.add_dset(dset_values_title, overwrite=True)
        dset_values.add_columns(x=[], yref=[])
        view_values = dset_values.add_view("Comparison", xlabel="X", ylabel="Y",
                                           autoscale=True, overwrite=True)

        dset_results = self.iodata.add_dset(dset_results_title, overwrite=True)
        dset_results.add_columns(variable_set=list(range(self.stack_count)),
                                 counts=all_counts.copy())
        dset_results.add_columns(**self.variables)
        view_results = dset_results.add_view("R-Factor analysis",
                                             xlabel="variable set number",
                                             ylabel="counts",
                                             title=("Number of R-Factors "
                                                    "minima for each set of "
                                                    "variables"))

        dset_rfc = self.iodata.add_dset(dset_rfc_title, overwrite=True)
        dset_rfc.add_columns(rfactor_number=list(range(12)))
        view_rfc = dset_rfc.add_view("R-Factor results", xlabel="R-Factor #",
                                     ylabel="R-Factor value",
                                     title="", autoscale=True, marker="s")

        for i in range(self.stack_count):
            rfc, exp_data, calc_data = self.io.load_results(index=i)
            # Store the experimental data
            dset_values.x, dset_values.yref = exp_data.T
            # Append the calculated values
            ycalc = calc_data[:,1]
            dset_values.add_columns(**{"calc{:d}".format(i): ycalc})
            dset_rfc.add_columns(**{'variable_set{:d}'.format(i): rfc})

        # Plot the curves
        view_values.select('x', 'yref', legend='Reference values')
        title = ''
        for k,v in self.best_values.items():
            title += '{}={} '.format(k, v)
        view_values.select('x', "calc{:d}".format(self.index),
                           legend="Best calculated values")
        view_values.set_plot_options(title=title)

        view_results.select('counts')

        for i in range(self.stack_count):
            view_rfc.select('rfactor_number', 'variable_set{:d}'.format(i),
                            legend="variables set #{:d}".format(i))
        # Save the parameters
        for p in self.get_parameters():
            bundle = {'group': str(p.group),
                      'name': str(p.name),
                      'value': str(p.value),
                      'unit': '' if p.unit is None else str(p.unit)}
            dset_results.add_parameter(**bundle)

        # Move back to the initial folder
        os.chdir(cwd)
        # And return the Data object
        return self.iodata

    def get_parameters(self):
        """Get all the defined parameters in the calculator.

        :return: A list of all parameters objects.
        :rtype: List of :py:class:`parameters.Parameter`

        """
        _ = []
        for section in ('general', ):
            parameters = getattr(self, section + '_parameters')
            for p in parameters:
                _.append(p)
        return _



if __name__ == "__main__":
    pass
