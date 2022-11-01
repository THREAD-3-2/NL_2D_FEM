============================
Documentation of `NL_2D_FEM`
============================

This Matlab code presents a novel method for
simulating the nonlinear dynamics of
highly flexible beam structures
based on a finite element discretization of
the geometrically exact beam model.
The finite element model,
currently limited to capturing in-plane motion,
is solved in the frequency domain
using a combination of harmonic balance and numerical continuation.
The solving scheme is automated with the solver MANLAB,
freely available to download at
`(the MANLAB website) <http://manlab.lma.cnrs-mrs.fr/>`_.
Several different computations are available to the user with this code,
notably classical forced response simulations,
but also computation of the nonlinear modes (i.e. the backbone curves). 

Contents
--------

.. toctree::
   installation
   walkthrough
   object
   theory
   gui
   useful
   inputs
   :maxdepth: 2
