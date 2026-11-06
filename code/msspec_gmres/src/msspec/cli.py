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
# along with msspec.  If not, see <http://www.gnu.org/licenses/>.
#
# Source file  : src/msspec/cli.py
# Last modified: Mon, 27 Sep 2021 17:49:48 +0200
# Committed by : sylvain tricot <sylvain.tricot@univ-rennes1.fr>


import sys
import textwrap
import subprocess
from msspec.iodata import Data
from pip._internal.cli.main import main as pipmain


def print_info():
    s = textwrap.dedent("""
        This application is a take-away Python environment
        to run msspec.

        Below is a list of all packages bundled with this
        interpreter.
        """)
    print(s)
    pipmain(['list'])
    s = textwrap.dedent("""
        To run a script with this interpreter, please launch
        this application with your script as an argument.

        I let you with the interactive environment. Type in
        <CTRL>+D to exit.

        """)
    print(s)

def load_py(*args):
    subprocess.run(['python'] + list(args))

def load_hdf5(fname):
    data = Data.load(fname)
    data.view()


def main():
    argv = sys.argv
    argc = len(argv)
    if argc > 1:
        if argv[1].endswith('.hdf5'):
            load_hdf5(argv[1])
        else:
            load_py(*argv[1:])
    else:
        print_info()
        subprocess.run(['ipython', '--profile=msspec'])
        exit()
