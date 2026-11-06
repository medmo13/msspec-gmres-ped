#!/usr/bin/env python
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
# Source file  : src/setup.py
# Last modified: mar. 07 avril 2020 17:01:42
# Committed by : "Sylvain Tricot <sylvain.tricot@univ-rennes1.fr>"

import sys
sys.path.insert(0, "msspec")
from setuptools import setup, find_packages
from version import __version__

with open('pip.freeze', 'r') as fd:
    REQUIREMENTS = fd.read().strip().split('\n')

if __name__ == "__main__":
    setup(name='msspec',
          version=__version__,
          include_package_data=True,
          packages=find_packages(include='msspec.*'),
          #setup_requires=SETUP_REQUIREMENTS,
          install_requires=REQUIREMENTS,

          author='Didier Sébilleau, Sylvain Tricot',
          author_email='sylvain.tricot@univ-rennes1.fr',
          maintainer='Sylvain Tricot',
          maintainer_email='sylvain.tricot@univ-rennes1.fr',
          url='https://msspec.cnrs.fr',
          description=('A multiple scattering package for sepectroscopies '
                       'using electrons to probe materials'),
          long_description="""MsSpec is a Fortran package to compute the
          cross-section of several spectroscopies involving one (or more)
          electron(s) as the probe. This package provides a python interface to
          control all the steps of the calculation.

          Available spectroscopies:
             * Photoelectron diffraction
             * Auger electron diffraction
             * Low energy electron diffraction
             * X-Ray absorption spectroscopy
             * Auger Photoelectron coincidence spectroscopy
             * Computation of the spectral radius""",
          download_url='https://msspec.cnrs.fr/downloads.html',
          # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
          classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Natural Language :: English',
            'Operating System :: Microsoft :: Windows :: Windows 10',
            'Operating System :: POSIX :: Linux',
            'Operating System :: MacOS :: MacOS X',
            'Programming Language :: Fortran',
            'Programming Language :: Python :: 3 :: Only',
            'Topic :: Scientific/Engineering :: Physics',
          ],
          keywords='spectroscopy atom electron photon multiple scattering',
          license='GPL',
          #entry_points={
          #    'console_scripts': ['msspec=msspec.cli:main']}
          )
