********************************
Photo-Electron Diffraction (PED)
********************************

Introduction
============
In PhotoElectron Diffraction, an incoming photon, with an energy in the X-ray range is absorbed by
an atom of the sample. A core electron of the absorbing atom is emitted and will eventually escape
from the sample after many scattering events toward a detector. This photo-electron is detected at
a given kinetic energy and for a given position (polar angle and azimutal angle) of the detector with
respect to the sample (see figure below). The distribution of electrons as a function of the sample 
polar or azimutal angles contains chemically resolved informations about the crystallography in the 
immediate proximity of the surface sample. This spectroscopy can also be done on Auger electrons. In 
this case, the technique is named Auger Electron Diffraction.

.. _ped_full_picture:

.. figure:: ../../full_picture.png
        :align: center
        :width: 80%

        The Full picture of a photoelectron diffraction process. a) The geometry of
        the experiment. b) A view of the multiple scattering process and c) Atomic
        energy level sketch of the normal and Auger photoemission process.


Quick reference
===============

To quickly start with MsSpec and Python, the easiest way is to read the :ref:`tutorials` section.
Here we summurize all the steps to perform a PED simulation:

        1. :ref:`Build your cluster <step1>` (thanks to the ASE python package)
        2. :ref:`Create a calculator <step2>` with :py:func:`MSSPEC`
        3. :ref:`Set the parameters <step3>` of your calculation.
        4. :ref:`Attach <step4>` your cluster to the calculator
        5. :ref:`Choose the absorber <step5>`
        6. :ref:`Compute <step6>` a scan
        7. :ref:`Plot <step7>` the results


.. _step1:

Build your cluster
------------------

Building a cluster means creating a list of atoms with their given positions in x, y, z coordinates.
It is easily done thanks to the `ase Python package <https://wiki.fysik.dtu.dk/ase/index.html>`_.

Because most of spectroscopies have a source and a detector in the same hemispherical space, a cluster
is often shaped as an half sphere. To create such atomic arrangements, special helper functions are provided
in the :py:mod:`utils` module.

For example to create an MgO cluster:

	.. literalinclude:: MgO.py
		:linenos:
		:lines: 1-15

.. only:: html

    will produce a cluster of 519 atoms like this:

    .. figure:: MgO.gif
        :align: center
        :width: 60%

        The shape of a typical (yet quite large) cluster used for a calculation (519 atoms of MgO).


.. only:: latex

    will produce a cluster of 519 atoms like this:

    .. figure:: MgO.png
        :align: center
        :width: 60%

        The shape of a typical (yet quite large) cluster used for a calculation (519 atoms of MgO).


.. _step2:

Create a calculator
-------------------

To create a claculator, you will use the :py:func:`calculator.MSSPEC` function. This function takes 4
keyword arguments:
        
        - **spectroscopy**, to specify the kind of spectroscopy. This is a string and 
          can be one of 'PED' for PhotoElectron Diffraction, 'AED' for Auger Electron
          Diffraction, 'APECS' for Auger PhotoElectron Coincidence Spectroscopy or
          'EXAFS' for Extended X-Ray Absorption Fine Structure.
        - **algorithm**, to choose between the matrix inversion method with the string
          'inversion' (best suited for lower kinetic energies < 100 eV), or the series
          expansion technique with 'expansion' or the correlation-expansion with
          'correlation'. The restarted GMRES solver, preconditioned by the exact
          tree factorization and with dense LU as the fallback, is selected with
          'gmres' (PED only). Its knobs are the ``linear_solver_tolerance``,
          ``linear_solver_max_iterations`` and ``linear_solver_restart`` parameters.
        - **polarization** to specify the light polarization. 'linear_qOz' or 'linear_xOy'
          for a linearly polarized light with the polarization vector in the :math:`(\vec{q}Oz)`
          or in the :math:`(xOy)` plane respectively. Finally choose 'circular' for circularly
          polarized light.
        - **folder**. Enter here the name of the folder used for temporary files.

The function returns a calculator object, so for example. To create a calculator for PhotoElectron 
Diffraction with the matrix inversion method:

.. code-block:: python

        calc = MSSPEC(spectroscopy = 'PED', algorithm = 'inversion')

To use the GMRES solver instead, with the same keyword:

.. code-block:: python

        calc = MSSPEC(spectroscopy = 'PED', algorithm = 'gmres')


.. _step3:

Set the parameters
------------------

A calculator has many parameters. They fall into 4 categories:

        * Muffin-tin parameters, to tweak the potential used for the phase shifts calculation
        * T-Matrix parameters, to control the T-matrix calculation
        * Calculation parameters, to tune the multiple scattering calculation: add atomic vibrations,
          add some filters to speed up the process, control the parameters of the series expansion method...
        * Spectroscopy dependent parameters. These parameters control -- for example -- the light source,
          the detector...

Each set of parameters is accessible through properties of the calculator object. For example, to tweak
the interstitial value of the Muffin Tin potential, use:

	>>> calc.muffintin_parameters.interstitial_potential = 12.1

To change the source energy, use:
	
	>>> calc.source_parameters.energy = 1253.0

All options are detailed in :ref:`this section <allparameters>`


.. _step4:

Attach your cluster
-------------------

Very easy! Juste use the :py:func:`set_atoms` function like this:

.. code-block:: python

        calc.set_atoms(cluster)

.. _step5:

Choose the absorber
-------------------

Set the absorber attribute of your cluster to the index of the atom you want it to be the absorber.
For example if the first atom of your cluster is the absorber

.. code-block:: python

        cluster.absorber = 0

The best way is to use a function to find the index based on the xyz coordinates of the atom. For
example to choose the closest atom of the origin:

.. code-block:: python
        
        cluster.absorber = get_atom_index(cluster, 0, 0, 0)

:py:func:`get_atom_index` is in the :py:mod:`utils` package so do not forget to import it. The
first argument is the cluster you will look for and the 3 next parameters are the x, y and z
coordinates.

.. _step6:

Compute
-------

You can compute 5 kinds of scans in PED spectroscopy:

        * A polar scan, with the :py:func:`get_theta_scan` method of the calculator object
        * An azimutal scan, with the :py:func:`MSSPEC.get_phi_scan` method
        * A stereographic scan, with the :py:func:`MSSPEC.get_theta_phi_scan` method
        * An energy scan, with the :py:func:`MSSPEC.get_energy_scan` method
        * The scattering factor, with the :py:func:`MSSPEC.get_scattering_factors` method

All these functions are used and detailed in the :ref:`tutorials <tutorials>`.

.. _step7:

Plot the results
----------------

Normally, the output of the previous functions is a :py:class:`iodata.Data` object. You can see
the results by typing:

>>> data = calc.get_theta_scan(...) # a polar scan for example
>>> data.view() # will popup a graphical window

