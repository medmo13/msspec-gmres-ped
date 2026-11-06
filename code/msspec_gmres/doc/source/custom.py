# coding: utf8
# vim: set et ts=4 sw=4 sts mouse=a fdm=indent:

import os
import sys
import subprocess
import re
import numpy as np
from lxml import etree

import msspec
from msspec.calculator import MSSPEC
from msspec.misc import set_log_level
set_log_level("error")


def get_package_version():
    version, *dev_version = msspec.__version__.split('.post')
    return version, "".join(dev_version)

def generate_download_page():
    full_version = '.post'.join(get_package_version()).strip('.post')
    setupfile_path = '../../package/MsSpec-{}.setup'.format(full_version)
    docfile_path = '../../package/MsSpec-{}.pdf'.format(full_version)

    with open('downloads.rst', 'w') as fd:
        content="""        
##########################
Download and install notes
##########################

.. warning::
    The content of this page is outdated. An updated version will be soon available.

click :download:`here <{}>` to download the last version of MsSpec and 
:download:`here <{}>` for this website as a single pdf file


.. include:: install_notes.inc
""".format(setupfile_path, docfile_path)
        fd.write(content)

def modify_banner(inputfile='./title_template.svg', outputfile='./title.svg'):
    version, dev_version = get_package_version()
    xml = etree.parse(inputfile)
    old_content = etree.tostring(xml)
    svg = xml.getroot()
    elements = svg.findall(".//{http://www.w3.org/2000/svg}tspan")

    for element in elements:
        tid = element.get('id')
        if tid.startswith('msspec_version'):
            element.text = version
        if tid.startswith('dev_version'):
            if dev_version:
                element.text = "post release: " + dev_version
            else:
                element.text = ""

    new_content = etree.tostring(xml).decode('utf-8')
    if new_content != old_content:
        with open(outputfile, 'w') as fd:
            fd.write(new_content)

def generate_parameters(spectroscopy=None):
    def get_content(all_parameters):
        content = ""
        for p in all_parameters:
            content += ".. _{0}-{1}:\n\n".format(group.lower(), p.name)
            content += "{0}\n{1}\n\n".format(p.name, "-"*len(p.name))
            content += ".. admonition:: details\n\n"

            table  = ("\n.. csv-table::\n"
                      "\t:widths: 100, 100\n\n")
            table += "\t\"*Types*\", \"{}\"\n".format(', '.join([_.__name__ for _ in p.allowed_types]))
            table += "\t\"*Limits*\", \"{} <= value <= {}\"\n".format(str(p.low_limit), str(p.high_limit))
            table += "\t\"*Unit*\", \"{}\"\n".format(str(p.unit))
            if type(p.allowed_values) in (tuple, list, np.ndarray):
                table += "\t\"*Allowed values*\", \"{}\"\n".format(', '.join([str(_) for _ in p.allowed_values]))
            else:
                table += "\t\"*Allowed values*\", \"{}\"\n".format(str(p.allowed_values))
            table += "\t\"*Default*\", \"{}\"\n".format(str(p.default))
            table += "\n\n\n\n"

            content += table.replace('\n', '\n\t')
            content += "{}\n\n\n\n".format(p.docstring)
        return content

    content = ""
    c = MSSPEC(spectroscopy='PED' if spectroscopy==None else spectroscopy)
    c.shutdown()

    allp = {}
    for p in c.get_parameters():
        if p.group not in allp.keys():
            allp[p.group] = []
        if not(p.private):
            allp[p.group].append(p)


    common_groups = ("GlobalParameters", "MuffintinParameters", "TMatrixParameters", "SourceParameters",
                     "DetectorParameters", "ScanParameters", "CalculationParameters")

    if spectroscopy == None:
        groups = common_groups
        fn = 'spectroscopies/common_parameters.inc'
        for group in groups:
            title = "{}".format(group)
            content += "{1}\n{0}\n\n".format("=" * len(title), title)
            content += get_content(allp[group])
    else:
        group = spectroscopy + 'Parameters'
        fn = 'spectroscopies/{0}/{0}_parameters.inc'.format(spectroscopy.lower())
        title = "{}".format(group)
        content += "{1}\n{0}\n\n".format("=" * len(title), title)
        content += get_content(allp[group])

    with open(fn, 'w') as fd:
        fd.write(content)
