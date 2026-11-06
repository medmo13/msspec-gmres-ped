Playing with the scattering factor
----------------------------------

In this tutorial we will play with an important parameter of any multiple
scattering calculation: the *scattering factor*.  When a electron wave is
scattered by an atom, the electron trajectory is modified after this event. The
particle will most likely continue its trajectory in the same direction
(forward scattering), but, to a lesser extent and depending on the atom and on
the electron energy, the direction of the scattered electron can change. The
electron can even be backscattered.  The angular distribution of the electron
direction after the scattering event is the scattering factor.

In a paper published in 1998, T. Gerber *et al.* used the quite high
backscattering factor of Rhodium atoms to probe the distance of Oxygen atoms
adsorbed on a Rhodium surface. Some electrons coming from Oxygen atoms are
ejected toward the Rhodium surface. They are then backscattered and interfere
with the direct signal comming from Oxygen atoms (see the figure below).  They
demonstrated both experimentally and numerically with a sinle scattering
computation that this lead to a very accurate probe of adsorbed species that
can be sensitive to bond length changes of the order of :math:`\pm 0.02 \mathring{A}`.

.. figure:: RhO_fig0.png
        :align: center
        :width: 30%

        Interferences produced by the backscattering effect.


First, compute the scattering factor of both chemical species, Rh and O.

.. literalinclude:: sf.py
        :linenos:

Running the above script should produce this polar plot. You can see that for Rhodium,
the backscattering coefficient is still significant even at quite high kinetic energy
(here 723 eV).

.. figure:: RhO_fig1.png
        :align: center
        :width: 80%

        Polar representation of the scattering factor

Let an Oxygen atom being adsorbed at a distance :math:`z_0` of an fcc site of
the (111) Rh surface.  and compute the :math:`\theta-\phi` scan for different
values of :math:`z_0`. You can see on the stereographic projection 3 bright
circles representing fringes of constructive interference between the direct
O(1s) photoelectron wave and that backscattered by the Rhodium atoms. The
center of these annular shapes changes from bright to dark due to the variation
of the Oxygen atom height above the surface which changes the path difference.

.. image:: RhO_fig2a.png
        :align: center
        :height: 200px


.. only:: html

        .. figure:: RhO_fig2b.gif
                :align: center
                :width: 80%

                Stereographic projections of O(1s) emission at :math:`E_k` = 723 eV for an
                oxygen atom on top of a fcc site of 3 Rh atoms at various altitudes
                :math:`z_0`

.. only:: latex

        .. figure:: RhO_fig2b.png
                :align: center
                :width: 80%

                Stereographic projections of O(1s) emission at :math:`E_k` = 723 eV for an
                oxygen atom on top of a fcc site of 3 Rh atoms at various altitudes
                :math:`z_0`




Here is the script for the computation. (:download:`download <RhO.py>`)

.. literalinclude:: RhO.py
        :linenos:

.. .. note::
..         After runing this script, you will get 20 images in your folder. You can merge them in one animated gif image
..         like this:
.. 
..         .. code-block:: bash
.. 
..                 convert -delay 50 -loop 0 image*.png animation.gif

.. seealso::

        X-Ray Photoelectron Diffraction in the Backscattering Geometry: A Key to Adsorption Sites and Bond Lengths at Surfaces
                T. Gerber, J. Wider, E. Welti & J. Osterwalder, Phys. Rev. Lett. **81** (8) p1654 (1998)
                `[doi] <https://doi.org/10.1103/PhysRevLett.81.1654>`__

