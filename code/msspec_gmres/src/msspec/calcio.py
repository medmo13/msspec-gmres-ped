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
# Source file  : src/msspec/calcio.py
# Last modified: Mon, 27 Sep 2021 17:49:48 +0200
# Committed by : sylvain tricot <sylvain.tricot@univ-rennes1.fr>


"""
Module calcio
=============

"""

from datetime import datetime
import re
import numpy as np
import ase
from ase import Atom, Atoms
from ase.data import chemical_symbols
import os

from msspec.misc import UREG, LOGGER
from msspec.utils import get_atom_index


class PhagenIO(object):
    def __init__(self, phagen_parameters, malloc_parameters):
        self.parameters = phagen_parameters
        self.malloc_parameters = malloc_parameters

        self.tl = None
        self.nat = None
        self.nateqm = None
        self.ne = None
        self.ipot = None
        self.lmax_mode = None
        self.nlmax = None
        self.energies = None

    def write_input_file(self, filename='input.ms'): # noqa: C901
        # in input folder
        LOGGER.debug("Writing Phagen input file")
        atoms = self.parameters.atoms
        if not atoms.has('mt_radii'):
            for a in atoms:
                a.set('mt_radius', 0.)

        if not atoms.has('mt_radii_scale'):
            for a in atoms:
                a.set('mt_radius_scale', 1.)

        radii = atoms.get_array('mt_radii')
        radii_scale = atoms.get_array('mt_radii_scale')

        if np.all(radii == 0) and np.all(radii_scale == 1):
            self.parameters.norman = 'stdcrm'
        elif np.all(radii == 0) and np.any(radii_scale != 1):
            self.parameters.norman = 'scaled'
            radii = radii_scale
        else:
            self.parameters.norman = 'extrad'
            radii *= radii_scale

        parameters = []
        for parameter in self.parameters:
            name = parameter.name
            if name in ('ionicity', 'atoms'):
                # skip ionicity and atoms parameters as they are treated in
                # other sections
                continue
            value = parameter.value
            if isinstance(value, str) and not re.match(r'\..*\.', value):
                s = ' {}=\'{:s}\''.format(name, str(parameter))
            else:
                s = ' {}={:s}'.format(name, str(parameter))
            parameters.append(s)

        header = " &job\n"
        header += ',\n'.join(parameters) + "\n &end\n"
        header += (' c Auto-generated file on {}\n Computes the T-matrix and '
                   'radial matrix elements\n\n').format(
            datetime.ctime(datetime.now()))

        # cluster section
        try:
            absorber = atoms.absorber
        except AttributeError as err:
            print(err)
            absorber = 0

        cluster_section = ''
        absorber_line = ''
        line_format = '{:>10s}{:4d}{:11.6f}{:14.6f}{:14.6f}{:13.4f}\n'
        for iat, atom in enumerate(atoms):
            symbol = atom.symbol
            if symbol == 'X':
                symbol = 'ES'
            number = atom.number
            x, y, z = atom.position
            r = radii[iat]
            cluster_line = line_format.format(symbol, number, x, y, z, r)
            # store absober line
            if atom.index == absorber:
                absorber_line = cluster_line
            else:
                cluster_section += cluster_line

        cluster_section = absorber_line + cluster_section
        cluster_section += '{:10d}{:4d}{:10.0f}{:9.0f}{:9.0f}{:8.0f}\n'.format(
            -1, -1, 0, 0, 0, 0)

        # Ionicity section
        ionicity = self.parameters.ionicity
        ionicity_section = ''
        ionicity_format = '{:8d}{:9.2f}\n'
        numbers = np.unique(atoms.numbers)
        for Z in numbers:
            try:
                symbol = chemical_symbols[Z]
                charge = ionicity[symbol]
            except KeyError:
                charge = 0
            ionicity_section += ionicity_format.format(Z, charge)
        ionicity_section += '{:8d}\n'.format(-1)

        content = header + cluster_section + ionicity_section

        # Write the content to filename
        try:
            with open(filename, 'r') as fd:
                old_content = fd.read()
        except IOError:
            old_content = ''

        pat = re.compile(r' c .*\n')

        modified = False
        if pat.sub('', content) != pat.sub('', old_content):
            with open(filename, 'w') as fd:
                fd.write(content)
            modified = True

        return modified

    def write_include_file(self, filename='msxas3.inc'):
        # read the whole include file content
        with open(filename, 'r') as fd:
            content = fd.read()

        # backup the content in memory
        old_content = content

        # replace the content
        for p in self.malloc_parameters:
            content = re.sub(r'({:s}\s*=\s*)\d+'.format(p.name),
                             r'\g<1>{:d}'.format(p.value), content)

        # actually write to the file only if different from the previous file
        modified = False
        if content != old_content:
            with open(filename, 'w') as fd:
                fd.write(content)
            modified = True
        return modified

    def load_tl_file(self, filename='tmatrix.tl'):
        atom_data = []

        # load all the file in the string 'content'
        with open(filename, 'r') as fd:
            content = fd.read()

        # convert the file to a (nat x ne)  array
        #
        # first, split the content in a list for each atom
        pattern = re.compile(r'-+\s*ATOM.*-+')
        lines = pattern.split(content)

        # get the number of atoms (nat) and the number of energy points (ne)
        nat, ne, _, ipot, lmax_mode = list(map(int,
                                               content.split('\n')[0].split()))
        self.nat = nat
        self.ne = ne
        self.ipot = ipot
        self.lmax_mode = lmax_mode

        # extract atom data
        for line in lines:
            numbers_str = ''.join(line.strip().split('\n')).split()
            numbers = []
            for n in numbers_str:
                if not re.match(r'^\d+$', n):
                    numbers.append(float(n))
            if len(numbers) > 0:
                array = np.array(numbers).reshape((-1, 4))
                atom_data.append(array)

        # construct the data array
        data = []
        for i in range(nat):
            data.append([])
            for j in range(ne):
                data[i].append(atom_data[j * nat + i])

        self.tl = data
        return data

    def write_tl_file(self, filename='tmatrix.tl'):

        def get_lmaxs(ie):
            lmaxs = np.zeros(int(4 * np.ceil(self.nat / 4.)), dtype=int)
            for a in range(self.nat):
                lmaxs[a] = len(self.tl[a][ie]) - 1
            lmaxs = lmaxs.reshape((-1, 4))
            return lmaxs

        def get_energies(unit='eV'):
            emin = self.parameters.emin
            emax = self.parameters.emax
            delta = self.parameters.delta
            energies = np.arange(emin, emax, delta)
            if len(energies) == 0:
                energies = np.array([emin])
            if energies[-1] + delta / 2 < emax:
                energies = np.append(energies, energies[-1] + delta)

            # conversion in eV
            if unit == 'eV':
                energies = (energies * UREG.Ry).to('eV')

            return energies

        def custom_strfloat(f):
            mantissa, exponent = '{:.10E}'.format(f).split('E')
            m = format(float(mantissa) / 10, '.6f').replace('-0', '-')
            e = format(int(exponent) + 1, '+03d')
            return '     {}E{}'.format(m, e)

        with open(filename, 'w') as fd:
            fd.write('{:>9}{:>9}{:>9}{:>9}{:>9}\n'.format(self.nat, self.ne, 1,
                                                          self.ipot,
                                                          self.lmax_mode))

            nlmax = 0
            for ie in range(self.ne):
                # write all lmaxs for each energy set
                lmaxs = get_lmaxs(ie)
                nlmax = int(max(nlmax, np.max(lmaxs)))
                fmt1 = '{:>9}' * 4 + '\n'
                fmt2 = '{:12.4f}{:10.4f}'
                fmt3 = '{:18.6E}'
                fmt_line = '   {:9.4f} {:9.4f}{:18.6E}{:18.6E}\n'
                for _ in lmaxs:
                    fd.write(fmt1.format(*_))

                for ia in range(self.nat):
                    # write the atom header line
                    fd.write('{}ATOM{:>4}{}\n'.format('-' * 26 + '  ', ia + 1,
                                                      '  ' + '-' * 23))
                    for _a, _b, _c, _d in self.tl[ia][ie]:
                        # fd.write(fmt2.format(_a, _b))
                        # fd.write(custom_strfloat(_c))
                        # fd.write(custom_strfloat(_d))
                        #fd.write('\n')
                        fd.write(fmt_line.format(_a, _b, _c, _d))
            self.nlmax = nlmax
            self.energies = get_energies()

    def load_cluster_file(self, filename='cluster.clu'):
        data = np.loadtxt(filename, skiprows=1, usecols=(0, 2, 3, 4, 5, 6))
        atoms = self.parameters.atoms

        absorber_position = atoms[atoms.absorber].position
        positions = data[:, 2:5] + absorber_position

        proto_indices = []
        for atom in atoms:
            i = np.argmin(np.linalg.norm(positions - atom.position, axis=1))
            proto_index = int(data[i, -1])
            proto_indices.append(proto_index)
        atoms.set_array('proto_indices', np.array(proto_indices))
        self.nateqm = int(np.max([
            len(np.where(data[:, -1] == i)[0]) for i in range(1, self.nat + 1)
        ]))

    def load_prototypical_atoms(self, filename='molinpot3.out'):
        with open(filename, 'r') as fd:
            content = fd.readlines()

        atoms = self.parameters.atoms
        content = content[-len(atoms):]
        proto = Atoms()
        for line in content:
            datatxt = np.array(re.split(r'\s+', line.strip(' \n')))
            neq = int(datatxt[-1])
            if neq == 0:
                Z = int(datatxt[1])
                xyz = (datatxt[[2, 3, 4]].astype(float) * UREG.a0).to(
                       'angstrom').magnitude
                sym = chemical_symbols[Z]
                proto += Atom(sym, position=xyz)
        # now translate the cluster so that the protocluster is at the
        # same origin than the original cluster
        origin0 = atoms[atoms.absorber].position
        origin1 = proto[0].position
        proto.translate(origin0 - origin1)
        # for each atom in this prototypical cluster, find the corresponding
        # atom index in the original cluster
        indices = []
        for atom in proto:
            indices.append(get_atom_index(atoms, *atom.position))
        self.prototypical_atoms = atoms[indices]
        return self.prototypical_atoms

    def load_potential_file(self, filename='plot_vc.dat'):
        a_index = 0
        pot_data = []
        with open(filename, 'r') as fd:
            data = fd.readlines()

        for d in data:
            if d[1] == 'a':
                a_index += 1
                d = d.split()
                a = {'Symbol': d[1], 'distance': float(d[4]),
                     'coord': np.array([float(d[7]), float(d[8]),
                                        float(d[9])]),
                     'index': int(a_index), 'values': []}
                pot_data.append(a)
            else:
                pot_data[a_index - 1]['values'].append(
                    tuple(float(_) for _ in d.split()))

        # convert the values list to a numpy array
        for _pot_data in pot_data:
            v = _pot_data['values']
            _pot_data['values'] = np.asarray(v)

        return pot_data


class SpecIO(object):
    def __init__(self, parameters, malloc_parameters, phagenio):
        self.parameters = parameters
        self.malloc_parameters = malloc_parameters
        self.phagenio = phagenio

    def write_input_file(self, filename='spec.dat'): # noqa: C901
        def title(t, shift=4, width=79, center=True):
            if center:
                s = ('{}*{:^%ds}*\n' % (width - shift - 2)
                     ).format(' ' * shift, t)
            else:
                s = ('{}*{:%ds}*\n' % (width - shift - 2)
                     ).format(' ' * shift, t)
            return s

        def rule(tabs=(5, 10, 10, 10, 10), symbol='=', shift=4, width=79):
            s = ' ' * shift + '*' + symbol * (width - shift - 2) + '*\n'
            t = np.cumsum(tabs) + shift
            sep = list(s)
            for i in t:
                sep[i] = '+'
            return ''.join(sep)

        def fillstr(a, b, index, justify='left'):
            alist = list(a)

            if justify == 'left':
                offset = -len(b) + 1
            elif justify == 'center':
                offset = (-len(b) + 1) / 2
            elif justify == 'decimal':
                try:
                    offset = -(b.index('.') - 1)
                except ValueError:
                    offset = 0
            else:
                offset = 0

            for i, _ in enumerate(b):
                alist[int(index + offset + i)] = _
            return ''.join(alist)

        def create_line(legend='', index=49, dots=False):
            s = ' ' * 79 + '\n'
            if dots:
                s = fillstr(s, "..", 6, justify='right')
            s = fillstr(s, "*", 4)
            s = fillstr(s, "*", 78)
            s = fillstr(s, legend, index, justify='right')
            return s

        p = self.parameters

        content = rule(tabs=(), symbol='*')
        content += title('spec input file')
        content += rule(tabs=(), symbol='*')
        content += rule(tabs=(), symbol='=')

        content += title('CRYSTAL STRUCTURE :')
        content += rule()
        line = create_line("CRIST,CENTR,IBAS,NAT")
        line = fillstr(line, "CUB", 9, 'left')
        line = fillstr(line, "P", 19, 'left')
        line = fillstr(line, format(0, 'd'), 29, 'left')
        line = fillstr(line, str(p.get_parameter('extra_nat')), 39, 'left')
        content += line
        line = create_line("A,BSURA,CSURA,UNIT")
        line = fillstr(line, format(1., '.4f'), 9, 'decimal')
        line = fillstr(line, format(1., '.3f'), 19, 'decimal')
        line = fillstr(line, format(1., '.3f'), 29, 'decimal')
        line = fillstr(line, "ATU", 39, 'left')
        content += line
        line = create_line("ALPHAD,BETAD,GAMMAD")
        line = fillstr(line, format(90., '.2f'), 9, 'decimal')
        line = fillstr(line, format(90., '.2f'), 19, 'decimal')
        line = fillstr(line, format(90., '.2f'), 29, 'decimal')
        content += line
        line = create_line("H,K,I,L")
        line = fillstr(line, format(0, 'd'), 9, 'left')
        line = fillstr(line, format(0, 'd'), 19, 'left')
        line = fillstr(line, format(0, 'd'), 29, 'left')
        line = fillstr(line, format(1, 'd'), 39, 'left')
        content += line
        line = create_line("NIV,COUPUR,ITEST,IESURF")
        line = fillstr(line, format(8, 'd'), 9, 'left')
        line = fillstr(line, format(1.4, '.2f'), 19, 'decimal')
        line = fillstr(line, format(0, 'd'), 29, 'left')
        line = fillstr(line, format(1, 'd'), 39, 'left')
        content += line
        line = create_line("ATBAS,CHEM(NAT),NZAT(NAT)")
        line = fillstr(line, format(0., '.6f'), 9, 'decimal')
        line = fillstr(line, format(0., '.6f'), 19, 'decimal')
        line = fillstr(line, format(0., '.6f'), 29, 'decimal')
        line = fillstr(line, format("Mg", '>2s'), 39, 'right')
        line = fillstr(line, format(12, '>2d'), 43, 'right')
        content += line
        line = create_line("ATBAS,CHEM(NAT),NZAT(NAT)")
        line = fillstr(line, format(0., '.6f'), 9, 'decimal')
        line = fillstr(line, format(0.5, '.6f'), 19, 'decimal')
        line = fillstr(line, format(0., '.6f'), 29, 'decimal')
        line = fillstr(line, format("O", '>2s'), 39, 'right')
        line = fillstr(line, format(8, '>2d'), 43, 'right')
        content += line
        line = create_line("VECBAS")
        line = fillstr(line, format(1., '.6f'), 9, 'decimal')
        line = fillstr(line, format(0., '.6f'), 19, 'decimal')
        line = fillstr(line, format(0., '.6f'), 29, 'decimal')
        content += line
        line = create_line("VECBAS")
        line = fillstr(line, format(0., '.6f'), 9, 'decimal')
        line = fillstr(line, format(1., '.6f'), 19, 'decimal')
        line = fillstr(line, format(0., '.6f'), 29, 'decimal')
        content += line
        line = create_line("VECBAS")
        line = fillstr(line, format(0., '.6f'), 9, 'decimal')
        line = fillstr(line, format(0., '.6f'), 19, 'decimal')
        line = fillstr(line, format(1., '.6f'), 29, 'decimal')
        content += line
        line = create_line("IREL,NREL,PCREL(NREL)")
        line = fillstr(line, format(0, 'd'), 9, 'left')
        line = fillstr(line, format(0, 'd'), 19, 'left')
        line = fillstr(line, format(0., '.1f'), 29, 'decimal')
        line = fillstr(line, format(0., '.1f'), 39, 'decimal')
        content += line
        line = create_line("OMEGA1,OMEGA2,IADS")
        line = fillstr(line, format(28., '.2f'), 9, 'decimal')
        line = fillstr(line, format(0., '.2f'), 19, 'decimal')
        line = fillstr(line, format(1, 'd'), 29, 'left')
        content += line
        content += rule()

        content += title('TYPE OF CALCULATION :')
        content += rule()
        line = create_line("SPECTRO,ISPIN,IDICHR,IPOL")
        line = fillstr(line, str(p.calctype_spectro), 9, 'left')
        line = fillstr(line,
                       str(p.get_parameter('calctype_ispin')), 19, 'left')
        line = fillstr(line, str(p.calctype_idichr), 29, 'left')
        line = fillstr(line, str(p.calctype_ipol), 39, 'left')
        content += line
        line = create_line("I_AMP")
        line = fillstr(line, str(p.calctype_iamp), 9, 'left')
        content += line
        content += rule()

        content += title('PhD EXPERIMENTAL PARAMETERS :')
        content += rule()
        line = create_line("LI,S-O,INITL,I_SO")
        line = fillstr(line, str(p.ped_li), 9, 'left')
        line = fillstr(line, str(p.ped_so), 19, 'center')
        line = fillstr(line, str(p.ped_initl), 29, 'left')
        line = fillstr(line, str(p.ped_iso), 39, 'left')
        content += line
        line = create_line("IPHI,ITHETA,IE,IFTHET")
        line = fillstr(line, str(p.ped_iphi), 9, 'left')
        line = fillstr(line, str(p.ped_itheta), 19, 'left')
        line = fillstr(line, str(p.ped_ie), 29, 'left')
        line = fillstr(line, str(p.ped_ifthet), 39, 'left')
        content += line
        line = create_line("NPHI,NTHETA,NE,NFTHET")
        line = fillstr(line, str(p.ped_nphi), 9, 'left')
        line = fillstr(line, str(p.ped_ntheta), 19, 'left')
        line = fillstr(line, str(p.ped_ne), 29, 'left')
        line = fillstr(line, str(p.ped_nfthet), 39, 'left')
        content += line
        line = create_line("PHI0,THETA0,E0,R0")
        line = fillstr(line, str(p.get_parameter('ped_phi0')), 9, 'decimal')
        line = fillstr(line, str(p.get_parameter('ped_theta0')), 19, 'decimal')
        line = fillstr(line, str(p.get_parameter('ped_e0')), 29, 'decimal')
        line = fillstr(line, str(p.get_parameter('ped_r0')), 39, 'decimal')
        content += line
        line = create_line("PHI1,THETA1,E1,R1")
        line = fillstr(line, str(p.get_parameter('ped_phi1')), 9, 'decimal')
        line = fillstr(line, str(p.get_parameter('ped_theta1')), 19, 'decimal')
        line = fillstr(line, str(p.get_parameter('ped_e1')), 29, 'decimal')
        line = fillstr(line, str(p.get_parameter('ped_r1')), 39, 'decimal')
        content += line
        line = create_line("THLUM,PHILUM,ELUM")
        line = fillstr(line, str(p.get_parameter('ped_thlum')), 9, 'decimal')
        line = fillstr(line, str(p.get_parameter('ped_philum')), 19, 'decimal')
        line = fillstr(line, str(p.get_parameter('ped_elum')), 29, 'decimal')
        content += line
        line = create_line("IMOD,IMOY,ACCEPT,ICHKDIR")
        line = fillstr(line, str(p.get_parameter('ped_imod')), 9, 'decimal')
        line = fillstr(line, str(p.get_parameter('ped_imoy')), 19, 'decimal')
        line = fillstr(line, str(p.get_parameter('ped_accept')), 29, 'decimal')
        line = fillstr(line,
                       str(p.get_parameter('ped_ichkdir')), 39, 'decimal')
        content += line
        content += rule()

        content += title(' ' * 22 + 'LEED EXPERIMENTAL PARAMETERS :',
                         center=False)
        content += rule()
        line = create_line("IPHI,ITHETA,IE,IFTHET")
        line = fillstr(line, str(p.get_parameter('leed_iphi')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('leed_itheta')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('leed_ie')), 29, 'left')
        line = fillstr(line, str(p.get_parameter('leed_ifthet')), 39, 'left')
        content += line
        line = create_line("NPHI,NTHETA,NE,NFTHET")
        line = fillstr(line, str(p.get_parameter('leed_nphi')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('leed_ntheta')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('leed_ne')), 29, 'left')
        line = fillstr(line, str(p.get_parameter('leed_nfthet')), 39, 'left')
        content += line
        line = create_line("PHI0,THETA0,E0,R0")
        line = fillstr(line, str(p.get_parameter('leed_phi0')), 9, 'decimal')
        line = fillstr(line,
                       str(p.get_parameter('leed_theta0')), 19, 'decimal')
        line = fillstr(line, str(p.get_parameter('leed_e0')), 29, 'decimal')
        line = fillstr(line, str(p.get_parameter('leed_r0')), 39, 'decimal')
        content += line
        line = create_line("PHI1,THETA1,E1,R1")
        line = fillstr(line, str(p.get_parameter('leed_phi1')), 9, 'decimal')
        line = fillstr(line,
                       str(p.get_parameter('leed_theta1')), 19, 'decimal')
        line = fillstr(line, str(p.get_parameter('leed_e1')), 29, 'decimal')
        line = fillstr(line, str(p.get_parameter('leed_r1')), 39, 'decimal')
        content += line
        line = create_line("TH_INI,PHI_INI")
        line = fillstr(line, str(p.get_parameter('leed_thini')), 9, 'decimal')
        line = fillstr(line,
                       str(p.get_parameter('leed_phiini')), 19, 'decimal')
        content += line
        line = create_line("IMOD,IMOY,ACCEPT,ICHKDIR")
        line = fillstr(line, str(p.get_parameter('leed_imod')), 9, 'decimal')
        line = fillstr(line, str(p.get_parameter('leed_imoy')), 19, 'decimal')
        line = fillstr(line,
                       str(p.get_parameter('leed_accept')), 29, 'decimal')
        line = fillstr(line, str(p.get_parameter('leed_ichkdir')), 39,
                       'decimal')
        content += line
        content += rule()

        content += title(' ' * 21 + 'EXAFS EXPERIMENTAL PARAMETERS :',
                         center=False)
        content += rule()
        line = create_line("EDGE,INITL,THLUM,PHILUM")
        line = fillstr(line, str(p.get_parameter('exafs_edge')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('exafs_initl')), 19, 'left')
        line = fillstr(line,
                       str(p.get_parameter('exafs_thlum')), 29, 'decimal')
        line = fillstr(line, str(p.get_parameter('exafs_philum')), 39,
                       'decimal')
        content += line
        line = create_line("NE,EK_INI,EK_FIN,EPH_INI")
        line = fillstr(line, str(p.get_parameter('exafs_ne')), 9, 'left')
        line = fillstr(line,
                       str(p.get_parameter('exafs_ekini')), 19, 'decimal')
        line = fillstr(line,
                       str(p.get_parameter('exafs_ekfin')), 29, 'decimal')
        line = fillstr(line, str(p.get_parameter('exafs_ephini')), 39,
                       'decimal')
        content += line
        content += rule()

        content += title(' ' * 22 + 'AED EXPERIMENTAL PARAMETERS :',
                         center=False)
        content += rule()
        line = create_line("EDGE_C,EDGE_I,EDGE_A")
        line = fillstr(line, str(p.get_parameter('aed_edgec')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('aed_edgei')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('aed_edgea')), 29, 'left')
        content += line
        line = create_line("I_MULT,MULT")
        line = fillstr(line, str(p.get_parameter('aed_imult')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('aed_mult')), 19, 'center')
        content += line
        line = create_line("IPHI,ITHETA,IFTHET,I_INT")
        line = fillstr(line, str(p.get_parameter('aed_iphi')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('aed_itheta')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('aed_ifthet')), 29, 'left')
        line = fillstr(line, str(p.get_parameter('aed_iint')), 39, 'left')
        content += line
        line = create_line("NPHI,NTHETA,NFTHET")
        line = fillstr(line, str(p.get_parameter('aed_nphi')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('aed_ntheta')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('aed_nfthet')), 29, 'left')
        content += line
        line = create_line("PHI0,THETA0,R0")
        line = fillstr(line, str(p.get_parameter('aed_phi0')), 9, 'decimal')
        line = fillstr(line, str(p.get_parameter('aed_theta0')), 19, 'decimal')
        line = fillstr(line, str(p.get_parameter('aed_r0')), 29, 'decimal')
        content += line
        line = create_line("PHI1,THETA1,R1")
        line = fillstr(line, str(p.get_parameter('aed_phi1')), 9, 'decimal')
        line = fillstr(line, str(p.get_parameter('aed_theta1')), 19, 'decimal')
        line = fillstr(line, str(p.get_parameter('aed_r1')), 29, 'decimal')
        content += line
        line = create_line("IMOD,IMOY,ACCEPT,ICHKDIR")
        line = fillstr(line, str(p.get_parameter('aed_imod')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('aed_imoy')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('aed_accept')), 29, 'decimal')
        line = fillstr(line, str(p.get_parameter('aed_ichkdir')), 39, 'left')
        content += line
        content += rule()

        content += title(' ' * 19 + 'EIGENVALUE CALCULATION PARAMETERS :',
                         center=False)
        content += rule()
        line = create_line("NE,EK_INI,EK_FIN,I_DAMP")
        line = fillstr(line, str(p.get_parameter('eigval_ne')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('eigval_ekini')), 19,
                       'decimal')
        line = fillstr(line, str(p.get_parameter('eigval_ekfin')), 29,
                       'decimal')
        line = fillstr(line, str(p.get_parameter('eigval_idamp')), 39, 'left')
        content += line

        if p.calctype_spectro == "EIG":
            nlines = int(np.ceil(p.eigval_ne / 4.))
        else:
            nlines = 1
        table = np.full((nlines, 4), str(p.get_parameter('eigval_ispectrum_ne').value), dtype='U64')
        for i in range(nlines):
            line = create_line("I_SPECTRUM(NE)")
            for j, o in enumerate((9, 19, 29, 39)):
                line = fillstr(line, table[i, j], o, 'left')
            content += line

        line = create_line("I_PWM,METHOD,ACC,EXPO")
        line = fillstr(line, str(p.get_parameter('eigval_ipwm')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('eigval_method')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('eigval_acc')), 29, 'decimal')
        line = fillstr(line,
                       str(p.get_parameter('eigval_expo')), 39, 'decimal')
        content += line
        line = create_line("N_MAX,N_ITER,N_TABLE,SHIFT")
        line = fillstr(line, str(p.get_parameter('eigval_nmax')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('eigval_niter')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('eigval_ntable')), 29, 'left')
        line = fillstr(line, str(p.get_parameter('eigval_shift')), 39,
                       'decimal')
        content += line
        line = create_line("I_XN,I_VA,I_GN,I_WN")
        line = fillstr(line, str(p.get_parameter('eigval_ixn')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('eigval_iva')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('eigval_ign')), 29, 'left')
        line = fillstr(line, str(p.get_parameter('eigval_iwn')), 39, 'left')
        content += line
        line = create_line("L,ALPHA,BETA")
        line = fillstr(line, str(p.get_parameter('eigval_l')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('eigval_alpha')), 19,
                       'decimal')
        line = fillstr(line, str(p.get_parameter('eigval_beta')), 29,
                       'decimal')
        content += line
        content += rule()

        content += title(' ' * 24 + 'CALCULATION PARAMETERS :', center=False)
        content += rule()
        line = create_line("NO,NDIF,ISPHER,I_GR")
        line = fillstr(line, str(p.get_parameter('calc_no')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('calc_ndif')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('calc_ispher')), 29, 'left')
        line = fillstr(line, str(p.get_parameter('calc_igr')), 39, 'left')
        content += line
        line = create_line("I_REN,N_REN,REN_R,REN_I")
        line = fillstr(line, str(p.get_parameter('calc_iren')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('calc_nren')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('calc_renr')), 29, 'decimal')
        line = fillstr(line, str(p.get_parameter('calc_reni')), 39, 'decimal')
        content += line
        line = create_line("ISFLIP,IR_DIA,ITRTL,I_TEST")
        line = fillstr(line, str(p.get_parameter('calc_isflip')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('calc_irdia')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('calc_itrtl')), 29, 'left')
        line = fillstr(line, str(p.get_parameter('calc_itest')), 39, 'left')
        content += line
        line = create_line("NEMET,IEMET(NEMET)")
        line = fillstr(line, format(1, 'd'), 9, 'left')
        line = fillstr(line, format(1, 'd'), 19, 'left')
        line = fillstr(line, format(0, 'd'), 29, 'left')
        line = fillstr(line, format(0, 'd'), 39, 'left')
        content += line
        line = create_line("ISOM,NONVOL,NPATH,VINT")
        line = fillstr(line, str(p.get_parameter('calc_isom')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('calc_nonvol')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('calc_npath')), 29, 'left')
        line = fillstr(line, str(p.get_parameter('calc_vint')), 39, 'decimal')
        content += line
        line = create_line("IFWD,NTHOUT,I_NO,I_RA")
        line = fillstr(line, str(p.get_parameter('calc_ifwd')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('calc_nthout')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('calc_ino')), 29, 'left')
        line = fillstr(line, str(p.get_parameter('calc_ira')), 39, 'left')
        content += line

        nat = p.extra_nat
        nra_arr = np.ones((nat), dtype=int)
        thfwd_arr = np.ones((nat))
        path_filtering = p.extra_parameters['calculation'].get_parameter(
            'path_filtering').value
        if (path_filtering is not None and
                'backward_scattering' in path_filtering):
            ibwd_arr = np.ones((nat), dtype=int)
        else:
            ibwd_arr = np.zeros((nat), dtype=int)
        thbwd_arr = np.ones((nat))
        for at in p.extra_atoms:
            i = at.get('proto_index') - 1
            thfwd_arr[i] = at.get('forward_angle')
            thbwd_arr[i] = at.get('backward_angle')
            nra_arr[i] = at.get('RA_cut_off')
        for i in range(p.extra_nat):
            line = create_line("N_RA,THFWD,IBWD,THBWD(NAT)", dots=True)
            line = fillstr(line, format(nra_arr[i], 'd'), 9, 'left')
            line = fillstr(line, format(thfwd_arr[i], '.2f'), 19, 'decimal')
            line = fillstr(line, format(ibwd_arr[i], 'd'), 29, 'left')
            line = fillstr(line, format(thbwd_arr[i], '.2f'), 39, 'decimal')
            content += line

        line = create_line("IPW,NCUT,PCTINT,IPP")
        line = fillstr(line, str(p.get_parameter('calc_ipw')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('calc_ncut')), 19, 'left')
        line = fillstr(line,
                       str(p.get_parameter('calc_pctint')), 29, 'decimal')
        line = fillstr(line, str(p.get_parameter('calc_ipp')), 39, 'left')
        content += line
        line = create_line("ILENGTH,RLENGTH,UNLENGTH")
        line = fillstr(line, str(p.get_parameter('calc_ilength')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('calc_rlength')), 19,
                       'decimal')
        line = fillstr(line, str(p.get_parameter('calc_unlength')), 29, 'left')
        content += line
        line = create_line("IDWSPH,ISPEED,IATT,IPRINT")
        # Here, if 'vibrational_damping' is None, use the 'debye_waller'
        # approach and the debye model for the mean square displacements
        # and set the temeprature to 0K and the Debye temeparture to 500K.
        calc_param = p.extra_parameters['calculation']
        if calc_param.vibrational_damping is None:
            idwsph = format(0, 'd')
            idcm = format(2, 'd')
            temp = format(0., '.2f')
            td = format(500., '.2f')
            LOGGER.warning('Vibrational damping is disabled for this '
                           'calculation.')
        else:
            idwsph = str(p.get_parameter('calc_idwsph'))
            idcm = str(p.get_parameter('calc_idcm'))
            temp = str(p.get_parameter('calc_t'))
            td = str(p.get_parameter('calc_td'))
        ispeed = str(p.get_parameter('calc_ispeed'))
        line = fillstr(line, idwsph, 9, 'left')
        line = fillstr(line, ispeed, 19, 'left')
        line = fillstr(line, str(p.get_parameter('calc_iatt')), 29, 'left')
        line = fillstr(line, str(p.get_parameter('calc_iprint')), 39, 'left')
        content += line
        line = create_line("IDCM,TD,T,RSJ")
        line = fillstr(line, idcm, 9, 'left')
        line = fillstr(line, td, 19, 'decimal')
        line = fillstr(line, temp, 29, 'decimal')
        line = fillstr(line, str(p.get_parameter('calc_rsj')), 39, 'decimal')
        content += line
        line = create_line("ILPM,XLPM0")
        line = fillstr(line, str(p.get_parameter('calc_ilpm')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('calc_xlpm0')), 19, 'decimal')
        content += line

        nat = p.extra_nat
        nlines = int(np.ceil(nat / 4.))
        uj2_array = np.zeros((4 * nlines))
        # Now, for each atom in the cluster, get the mean_square_vibration and
        # store it in the index corresponding to the prototypical index
        for at in p.extra_atoms:
            i = at.get('proto_index') - 1
            msq_vib = at.get('mean_square_vibration')
            uj2_array[i] = msq_vib
        uj2_array = uj2_array.reshape((nlines, 4))
        for i in range(nlines):
            line = create_line("UJ2(NAT)", dots=True)
            for j, o in enumerate((9, 19, 29, 39)):
                line = fillstr(line, format(uj2_array[i, j], '.5f'), o,
                               'decimal')
            content += line
        content += rule()

        content += title(' ' * 17 + 'INPUT FILES (PHD, EXAFS, LEED, AED, '
                                    'APECS) :', center=False)
        content += rule(tabs=(), symbol='-')
        content += title(
            ' ' * 8 + 'NAME' + ' ' * 20 + 'UNIT' + ' ' * 16 + 'TYPE',
            center=False)
        content += rule(tabs=(5, 23, 7, 10))
        line = create_line("DATA FILE,UNIT")
        line = fillstr(line, str(p.get_parameter('input_data')), 9, 'right')
        line = fillstr(line, str(p.get_parameter('input_unit00')), 39, 'left')
        content += line
        line = create_line("PHASE SHIFTS/TL FILE,UNIT")
        line = fillstr(line, str(p.get_parameter('input_tl')), 9, 'right')
        line = fillstr(line, str(p.get_parameter('input_unit01')), 39, 'left')
        content += line
        line = create_line("RADIAL MATRIX ELTS FILE,UNIT")
        line = fillstr(line, str(p.get_parameter('input_rad')), 9, 'right')
        line = fillstr(line, str(p.get_parameter('input_unit02')), 39, 'left')
        content += line
        line = create_line("CLUSTER FILE,UNIT")
        line = fillstr(line, str(p.get_parameter('input_cluster')), 9, 'right')
        line = fillstr(line, str(p.get_parameter('input_unit03')), 39, 'left')
        content += line
        line = create_line("ADSORBATE FILE,UNIT")
        line = fillstr(line, str(p.get_parameter('input_adsorbate')), 9,
                       'right')
        line = fillstr(line, str(p.get_parameter('input_unit04')), 39, 'left')
        content += line
        line = create_line("K DIRECTIONS FILE,UNIT")
        line = fillstr(line, str(p.get_parameter('input_kdirs')), 9, 'right')
        line = fillstr(line, str(p.get_parameter('input_unit05')), 39, 'left')
        content += line
        content += rule(tabs=(5, 23, 7, 10))

        content += title(' ' * 21 + 'ADDITIONAL INPUT FILES (APECS) :',
                         center=False)
        content += title(' ' * 28 + '(AUGER ELECTRON)', center=False)
        content += rule(tabs=(), symbol='-')
        content += title(
            ' ' * 8 + 'NAME' + ' ' * 20 + 'UNIT' + ' ' * 16 + 'TYPE',
            center=False)
        content += rule(tabs=(5, 23, 7, 10))
        line = create_line("PHASE SHIFTS/TL FILE,UNIT")
        line = fillstr(line, str(p.get_parameter('input2_tl')), 9, 'right')
        line = fillstr(line, str(p.get_parameter('input2_unit06')), 39, 'left')
        content += line
        line = create_line("RADIAL MATRIX ELTS FILE,UNIT")
        line = fillstr(line, str(p.get_parameter('input2_rad')), 9, 'right')
        line = fillstr(line, str(p.get_parameter('input2_unit07')), 39, 'left')
        content += line
        line = create_line("K DIRECTIONS FILE,UNIT")
        line = fillstr(line, str(p.get_parameter('input2_kdirs')), 9, 'right')
        line = fillstr(line, str(p.get_parameter('input2_unit08')), 39, 'left')
        content += line
        content += rule(tabs=(5, 23, 7, 10))

        content += title(' ' * 29 + 'OUTPUT FILES :', center=False)
        content += rule(tabs=(), symbol='-')
        content += title(
            ' ' * 8 + 'NAME' + ' ' * 20 + 'UNIT' + ' ' * 16 + 'TYPE',
            center=False)
        content += rule(tabs=(5, 23, 7, 10))
        line = create_line("CONTROL FILE,UNIT")
        line = fillstr(line, str(p.get_parameter('output_log')), 9, 'right')
        line = fillstr(line, str(p.get_parameter('output_unit09')), 39, 'left')
        content += line
        line = create_line("RESULT FILE,UNIT")
        line = fillstr(line, str(p.get_parameter('output_res')), 9, 'right')
        line = fillstr(line, str(p.get_parameter('output_unit10')), 39, 'left')
        content += line
        line = create_line("SCATTERING FACTOR FILE,UNIT")
        line = fillstr(line, str(p.get_parameter('output_sf')), 9, 'right')
        line = fillstr(line, str(p.get_parameter('output_unit11')), 39, 'left')
        content += line
        line = create_line("AUGMENTED CLUSTER FILE,UNIT")
        line = fillstr(line,
                       str(p.get_parameter('output_augclus')), 9, 'right')
        line = fillstr(line, str(p.get_parameter('output_unit12')), 39, 'left')
        content += line
        content += rule(tabs=(5, 23, 7, 10))

        content += title(' ' * 26 + 'END OF THE DATA FILE', center=False)
        content += rule(tabs=())
        content += rule(tabs=(), symbol='*')

        try:
            with open(filename, 'r') as fd:
                old_content = fd.read()
        except IOError:
            old_content = ''

        modified = False
        if content != old_content:
            with open(filename, 'w') as fd:
                fd.write(content)
                LOGGER.debug("Writing Spec input file written in {}".format(filename))
            modified = True

        return modified

    def write_include_file(self, filename='spec.inc'):
        def get_li(level):
            orbitals = 'spdfghi'
            m = re.match(r'\d(?P<l>[%s])(\d/2)?' % orbitals, level)
            return orbitals.index(m.group('l'))

        requirements = {
            'NATP_M': self.phagenio.nat,
            'NATCLU_M': len(self.parameters.extra_atoms),
            'NAT_EQ_M': self.phagenio.nateqm,
            'N_CL_L_M': 0,
            'NE_M': self.phagenio.ne,
            'NL_M': self.phagenio.nlmax + 1,
            'LI_M': get_li(self.parameters.extra_level) + 1,
            'NEMET_M': 1,
            'NO_ST_M': self.parameters.calc_no,
        }

        # read the include file
        with open(filename, 'r') as fd:
            content = fd.read()

        # backup the content in memory
        old_content = content

        # for key in ('NATP_M', 'NATCLU_M', 'NE_M', 'NEMET_M', 'LI_M', 'NL_M',
        #             'NO_ST_M'):
        #     required = requirements[key]
        #     limit = self.malloc_parameters.get_parameter(key).value
        #     value = required if required > limit else limit
        #     content = re.sub(r'({:s}\s*=\s*)\d+'.format(key),
        #                      r'\g<1>{:d}'.format(value), content)

        for key in ('NAT_EQ_M', 'N_CL_N_M', 'NDIF_M', 'NSO_M', 'NTEMP_M',
                    'NODES_EX_M', 'NSPIN_M', 'NTH_M', 'NPH_M', 'NDIM_M',
                    'N_TILT_M', 'N_ORD_M', 'NPATH_M', 'NGR_M'):
            value = self.malloc_parameters.get_parameter(key).value
            content = re.sub(r'({:s}\s*=\s*)\d+'.format(key),
                             r'\g<1>{:d}'.format(value), content)

        for key, value in list(requirements.items()):
            content = re.sub(r'({:s}\s*=\s*)\d+'.format(key),
                             r'\g<1>{:d}'.format(value), content)

        modified = False
        if content != old_content:
            with open(filename, 'w') as fd:
                fd.write(content)
            modified = True
        return modified

    def write_kdirs_file(self, filename='kdirs.dat'):
        fwhm = 1.
        all_theta = self.parameters.extra_parameters['scan'].theta
        all_phi = self.parameters.extra_parameters['scan'].phi
        f = '{:7}{:4}{:6}\n'

        old_content = None
        try:
            with open(filename, 'r') as fd:
                old_content = fd.read()
        except IOError:
            pass

        content = ''
        content += f.format(2, 1, len(all_theta) * len(all_phi))
        content += f.format(1, len(all_phi), len(all_theta))
        for iphi, phi in enumerate(all_phi):
            for itheta, theta in enumerate(all_theta):
                s = '{:5}{:5}{:13.3f}{:11.3f}{:15e}\n'
                s = s.format(iphi + 1, itheta + 1, theta, phi, fwhm)
                content += s

        modified = False
        if content != old_content:
            with open(filename, 'w') as fd:
                fd.write(content)
            modified = True
        return modified

    def load_results(self, filename='results.dat'):
        rows2skip = {
            'PED': 27,
            'AED': 27,
            'EXAFS': 27,
            'LEED': 26,
            'EIG': 0
        }
        spectro = self.parameters.extra_parameters['global'].spectroscopy
        skip = rows2skip[spectro]

        data = np.loadtxt(filename, skiprows=skip, unpack=True)
        if len(data.shape) <= 1:
            data = data.reshape((2, -1))
        return data

    def load_facdif(self, filename='facdif1.dat'):
        data = np.loadtxt(filename, skiprows=1)
        return data

    def load_log(self, filename='spec.log'):
        pat = re.compile(r'ORDER\s+(\d+)\s+TOTAL NUMBER OF PATHS\s+:\s+(\d+)')
        with open(filename, 'r') as fd:
            content = fd.read()
        return pat.findall(content)


class CompCurveIO(object):
    def __init__(self, parameters):
        self.parameters = parameters

    def write_input_file(self, filename='comp_curve.dat'):
        def title(t, shift=4, width=78, center=True):
            if center:
                s = ('{}*{:^%ds}*\n' % (width - shift - 2)
                     ).format(' ' * shift, t)
            else:
                s = ('{}*{:%ds}*\n' % (width - shift - 2)
                     ).format(' ' * shift, t)
            return s

        def rule(tabs=(5, 10, 10, 10, 10), symbol='=', shift=4, width=78):
            s = ' ' * shift + '*' + symbol * (width - shift - 2) + '*\n'
            t = np.cumsum(tabs) + shift
            sep = list(s)
            for i in t:
                sep[i] = '+'
            return ''.join(sep)

        def fillstr(a, b, index, justify='left'):
            alist = list(a)

            if justify == 'left':
                offset = -len(b) + 1
            elif justify == 'center':
                offset = (-len(b) + 1) / 2
            elif justify == 'decimal':
                try:
                    offset = -(b.index('.') - 1)
                except ValueError:
                    offset = 0
            else:
                offset = 0

            for i, _ in enumerate(b):
                alist[int(index + offset + i)] = _
            return ''.join(alist)

        def create_line(legend='', index=49, dots=False):
            s = ' ' * 78 + '\n'
            if dots:
                s = fillstr(s, "..", 6, justify='right')
            s = fillstr(s, "*", 4)
            s = fillstr(s, "*", 77)
            s = fillstr(s, legend, index, justify='right')
            return s

        p = self.parameters

        content = rule(tabs=(), symbol='*')
        content += title('R-FACTOR, SIMILARITY INDEX, DISTANCE, GOODNESS OF FIT')
        content += title('KERNEL DISTANCE AND SHAPE ANALYSIS')
        content += rule(tabs=(), symbol='*')

        content += rule(tabs=(5,10,10,10,10))
        content += create_line('CALCULATION PARAMETERS : GENERAL', index=29)
        content += rule(tabs=(5,10,10,10,10))
        line = create_line("N_PAR,NORM,I_SCALE,I_NORM")
        line = fillstr(line, str(p.get_parameter('general_npar')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('general_norm')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('general_iscale')), 29, 'left')
        line = fillstr(line, str(p.get_parameter('general_inorm')), 39, 'left')
        content += line
        line = create_line("I_SYM,SYM,I_POSI")
        line = fillstr(line, str(p.get_parameter('general_isym')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('general_sym')), 19, 'decimal')
        line = fillstr(line, str(p.get_parameter('general_iposi')), 29, 'left')
        content += line
        line = create_line("I_DIST,I_CUR,I_SA,I_PRINT")
        line = fillstr(line, str(p.get_parameter('general_idist')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('general_icur')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('general_isa')), 29, 'left')
        line = fillstr(line, str(p.get_parameter('general_iprint')), 39, 'left')
        content += line

        content += rule(tabs=(5,10,10,10,10))
        content += create_line('CALCULATION PARAMETERS : WEIGHTS', index=29)
        content += rule(tabs=(5,10,10,10,10))
        line = create_line("I_WEIGHT,ALPHA,BETA,SIGMA")
        line = fillstr(line, str(p.get_parameter('weights_iweight')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('weights_alpha')), 19, 'decimal')
        line = fillstr(line, str(p.get_parameter('weights_beta')), 29, 'decimal')
        line = fillstr(line, str(p.get_parameter('weights_sigma')), 39, 'decimal')
        content += line
        line = create_line("I_SHIFT,MAXW")
        line = fillstr(line, str(p.get_parameter('weights_ishift')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('weights_maxw')), 19, 'decimal')
        content += line

        content += rule(tabs=(5,10,10,10,10))
        content += create_line('CALCULATION PARAMETERS : R-FACTORS', index=29)
        content += rule(tabs=(5,10,10,10,10))
        line = create_line("V_I")
        line = fillstr(line, str(p.get_parameter('rfc_vi')), 9, 'decimal')
        content += line

        content += rule(tabs=(5,10,10,10,10))
        content += create_line('CALCULATION PARAMETERS : SIMILARITY INDICES', index=29)
        content += rule(tabs=(5,10,10,10,10))
        line = create_line("ALHPAS,BETAS,N_BINS")
        line = fillstr(line, str(p.get_parameter('sim_alphas')), 9, 'decimal')
        line = fillstr(line, str(p.get_parameter('sim_betas')), 19, 'decimal')
        line = fillstr(line, str(p.get_parameter('sim_nbins')), 29, 'left')
        content += line

        content += rule(tabs=(5,10,10,10,10))
        content += create_line('CALCULATION PARAMETERS : DISTANCES', index=29)
        content += rule(tabs=(5,10,10,10,10))
        line = create_line("ALHPAD,I_BETA,L,SIGMAD")
        line = fillstr(line, str(p.get_parameter('dist_alphad')), 9, 'decimal')
        line = fillstr(line, str(p.get_parameter('dist_ibeta')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('dist_l')), 29, 'left')
        line = fillstr(line, str(p.get_parameter('dist_sigmad')), 39, 'decimal')
        content += line

        content += rule(tabs=(5,10,10,10,10))
        content += create_line('CALCULATION PARAMETERS : GOODNESS OF FIT', index=29)
        content += rule(tabs=(5,10,10,10,10))
        line = create_line("N_BING,ALPHAG")
        line = fillstr(line, str(p.get_parameter('gof_nbing')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('gof_alphag')), 19, 'decimal')
        content += line

        content += rule(tabs=(5,10,10,10,10))
        content += create_line('CALCULATION PARAMETERS : KERNEL DISTANCES', index=29)
        content += rule(tabs=(5,10,10,10,10))
        line = create_line("ALPHAK,L,SIGMAK")
        line = fillstr(line, str(p.get_parameter('kdist_alphak')), 9, 'decimal')
        line = fillstr(line, str(p.get_parameter('kdist_l')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('kdist_sigmak')), 29, 'decimal')
        content += line

        content += rule(tabs=(5,10,10,10,10))
        content += create_line('CALCULATION PARAMETERS : MOMENTS', index=29)
        content += rule(tabs=(5,10,10,10,10))
        line = create_line("N_GRID,N_MOM,BASIS")
        line = fillstr(line, str(p.get_parameter('mom_ngrid')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('mom_nmom')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('mom_basis')), 29, 'left')
        content += line
        line = create_line("I_ALG,MU,NU")
        line = fillstr(line, str(p.get_parameter('mom_ialg')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('mom_mu')), 19, 'decimal')
        line = fillstr(line, str(p.get_parameter('mom_nu')), 29, 'decimal')
        content += line

        content += rule(tabs=(5,10,10,10,10))
        content += create_line('CALCULATION PARAMETERS : CHORDS', index=29)
        content += rule(tabs=(5,10,10,10,10))
        line = create_line("I_CHORD,METHOD,VALUE,N_BINC")
        line = fillstr(line, str(p.get_parameter('chords_ichord')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('chords_method')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('chords_value')), 29, 'left')
        line = fillstr(line, str(p.get_parameter('chords_nbinc')), 39, 'left')
        content += line

        content += rule(tabs=(5,10,10,10,10))
        content += create_line('CALCULATION PARAMETERS : CHAIN CODES', index=29)
        content += rule(tabs=(5,10,10,10,10))
        line = create_line("N_CONNECT,SCALEC")
        line = fillstr(line, str(p.get_parameter('codes_nconnect')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('codes_scalec')), 19, 'decimal')
        content += line

        content += rule(tabs=(5,10,10,10,10))
        content += create_line('CALCULATION PARAMETERS : CONTOUR', index=29)
        content += rule(tabs=(5,10,10,10,10))
        line = create_line("NBIN,N_LEN,SH_AN,I_FOU")
        line = fillstr(line, str(p.get_parameter('cont_nbin')), 9, 'left')
        line = fillstr(line, str(p.get_parameter('cont_nlen')), 19, 'left')
        line = fillstr(line, str(p.get_parameter('cont_shan')), 29, 'left')
        line = fillstr(line, str(p.get_parameter('cont_ifou')), 39, 'left')
        content += line
        line = create_line("I_NORM")
        line = fillstr(line, str(p.get_parameter('cont_ifou')), 9, 'left')
        content += line

        content += rule(tabs=(5,10,10,10,10))
        content += create_line('EXPERIMENTAL INPUT FILE :', index=29)
        content += rule(tabs=(), symbol='-')
        line = create_line('')
        line = fillstr(line, 'NAME', 14, 'right')
        line = fillstr(line, 'TYPE', 58, 'right')
        content += line
        content += rule(tabs=(5,40))
        line = create_line('EXPERIMENT', index=49)
        line = fillstr(line, str(p.get_parameter('exp_filename')), 9, 'right')
        content += line

        content += rule(tabs=(5,40))
        content += create_line('CALCULATED INPUT FILE :', index=29)
        content += rule(tabs=(), symbol='-')
        line = create_line('')
        line = fillstr(line, 'NAME', 14, 'right')
        line = fillstr(line, 'PARAMETER 1', 49, 'right')
        line = fillstr(line, 'PARAMETER 2', 63, 'right')
        content += line
        content += rule(tabs=(5,40,7,7,7))
        calc_fnames = p.get_parameter('calc_filename').value
        for fname in calc_fnames:
            line = create_line('')
            line = fillstr(line, fname, 9, 'right')
            line = fillstr(line, str(p.get_parameter('calc_param1')), 56, 'decimal')
            line = fillstr(line, str(p.get_parameter('calc_param2')), 70, 'decimal')
            content += line

        content += rule(tabs=(5,40,7,7,7))
        content += create_line('END OF DATA FILE', index=31)
        content += rule(tabs=())
        content += rule(tabs=(),symbol='*')

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

    def load_results(self, index=0, prefix='rfc/experiment_rf'):
        data = []
        for i in range(1, 13):
            #data.append(np.loadtxt(prefix + f'{i:02d}' + '.txt')[-1])
            results = np.loadtxt(prefix + '{:02d}'.format(i) + '.txt')
            results = results.reshape((-1, 2))
            data.append(results[index,1])
        suffix = 'ren'
        exp = {'int': None, 'ren': None, 'chi': None, 'cdf': None}
        exp_ren = np.loadtxt(os.path.join('exp', 'div',
                                          'experiment_{}.txt'.format(suffix)))
        calc_ren = np.loadtxt(os.path.join('calc', 'div',
                                           'calculation{:d}_{}.txt'.format(index,suffix)))
        return data, exp_ren, calc_ren
