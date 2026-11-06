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
# Source file  : src/msspec/version.py
# Last modified: Wed, 26 Oct 2022 17:15:24 +0200
# Committed by : Sylvain Tricot <sylvain.tricot@univ-rennes1.fr> 1666797324 +0200


import os

from importlib.metadata import version
import subprocess

# find the version number
# 1- Try to read it from the git info
# 2- If it fails, try to read it from the VERSION file
# 3- If it fails, try to read it from the distribution file

PKGNAME = 'msspec'

thisfile_path = os.path.abspath(__file__)
thisfile_dir = os.path.dirname(thisfile_path)

try:
    cmd = ["git describe|sed 's/-\([0-9]\+\)-.*/.dev\\1/g'"]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, 
            shell=True, cwd=thisfile_dir)
    __version__ = result.stdout.decode('utf-8').strip()
    if __version__ == "":
        raise
except Exception as err:
    try:
        versionfile = os.path.join(thisfile_dir, "./VERSION")
        with open(versionfile, "r") as fd:
            __version__ = fd.readline().strip()
    except Exception as err:
        try:
            __version__ = version(PKGNAME)
        except Exception as err:
            __version__ = "0.0.0"
