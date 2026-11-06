#!/usr/bin/env python
import glob
import logging
import os
import sys
from msspec.calculator import MSSPEC
from msspec.utils import get_atom_index
from msspec.utils import hemispherical_cluster
from msspec.utils import SPRKKRPotential
from sprkkr.calculator import SPRKKR

from ase.build import bulk



logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Create a copper cell
Cu = bulk('Cu')

# ########## SPRKKR part
if 'sprkkr' in sys.argv:
    # create a SPRKKR calculator
    calc = SPRKKR(label="Cu/Cu")

    # attach the atoms to the calculator object
    calc.set_atoms(Cu)

    # Here is how to modify input file
    # calc.input.control_section.set(DATASET="Fe", ADSI="SCF", POTFIL="Fe.pot")
    # calc.input.tau_section.set(nktab=250)
    # launch kkrscf
    calc.get_potential_energy()

    # calc2 = SPRKKR(label="Fe/Fe", task="phagen")
    # calc2.set_atoms(Fe)
    # calc2.input.control_section.set(POTFIL="Fe.pot_new")
    # calc2.phagen()

    #
    # EXPORT POTENTIAL FOR PHAGEN
    #
    #change task and command
    calc.set_command('PHAGEN')
    # Change output file
    calc.set_outfile('Cu_phagen.out')
    #to change task we need to replace input file tempate
    calc.set_inpfile("Cu_phagen.inp")
    calc.input.control_section.set(DATASET="PHAGEN", ADSI="PHAGEN")
    #set potetential file to converged potential
    conv_potfile=os.path.join(calc.potfile+'_new')

    calc.set_potfile(conv_potfile)
    #run given task
    calc.phagen()


# ######### MsSpec part
if 'msspec' in sys.argv:
    pot = SPRKKRPotential(Cu, "Cu/Cu.pot_new", *glob.glob("Cu/*PHAGEN.pot"))

    nplanes = 3
    cluster = hemispherical_cluster(Cu, planes=nplanes,
                                    emitter_plane=nplanes-1)
    cluster.absorber = get_atom_index(cluster, 0, 0, 0)

    calc = MSSPEC(folder="calc")
    calc.set_atoms(cluster)
    calc.tmatrix_parameters.potential = pot

    data = calc.get_theta_scan(level='2p3/2')
    data.view()

if len(sys.argv) <= 1:
    print("Please specify either 'sprkkr', 'msspec' keywords or both "
          "of them on the command line")
