From single scattering to multiple scattering effects
-----------------------------------------------------

In this tutorial we will observe one of the main difference between a single
and a multiple scatteirng process: the *defocusing effect*.

We will reproduce a calculation performed by M.-L Xu, J.J. Barton and M.A. Van Hove
in their paper of 1989 (see :ref:`below <refTuto2>` )
In the spirit of figure 3 of their paper, We create 3 atomic chains of Ni 
atoms (2, 3 and 5 atoms) tilted by 45° and we compare the intensity of the forward scattering 
peak as a function of the scattering order: for single scattering and for a multiple 
scattering at the 5th order.

Below is a commented version of this example. You can download it :download:`here
<Ni_chain.py>`

.. literalinclude:: Ni_chain.py
        :linenos:


The resulting ouput is reported below. We added a sketch of the atomic chains
on the left hand side of the figure. You can see that for the single scattering
computation (left column), the forward peak is growing with increasing the
number of atoms in the chain, while it is clearly damped when using the
multiple scattering approach. Electrons are focused by the second atom in the
chain and *over* focused again by the third atom leading to a diverging
trajectory for this electron which in turn lowers the signal.

.. figure:: Ni_fig1.png
        :align: center
        :width: 80%

        Figure 1. Polar scan of a Ni chain of 2-5 atoms for single and mutliple (5th order)
        scattering.


.. _refTuto2:

.. seealso::

        Electron scattering by atomic chains: Multiple-scattering effects
                M.-L Xu, J.J. Barton & M.A. Van Hove, Phys. Rev. B **39** (12) p8275 (1989)
                `[doi] <https://doi.org/10.1103/PhysRevB.39.8275>`__

