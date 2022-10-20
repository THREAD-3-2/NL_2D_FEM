.. _theory:

Theory
======

Finite element model (FEM) of geometrically exact beam structures
-----------------------------------------------------------------

The principal advantage of the geometrically exact beam model
for the simulation of highly flexible structures is that
there is no truncation or linearization of the rotation terms,
meaning that the model is valid even at extreme amplitudes of cross-section rotation.
A continuum mechanics-based derivation of the geometrically exact beam model
based on the weak form / principle of virtual work is presented in detail in
`(O. Thomas, A. Sénéchal, J.-F. Deü, 2016) <https://link.springer.com/article/10.1007/s11071-016-2965-0>`_
and is not recopied here.
Instead, we focus primarily on the discretization of the model into finite elements.

The beam is discretized into :math:`N_e` elements of length :math:`L^e`,
where at each node :math:`i` of an element there are
three degrees of freedom :math:`\mathbf{u}_i = \left( u_i, w_i, \theta_i \right)`
corresponding to the axial displacement, transverse displacement and cross-section rotation, respectively.
The elements used in this code are two-node Timoshenko beam elements with linear shape functions.
The vector :math:`\mathbf{q}^e` contains all degrees of freedom at both nodes of the element:

.. math::
    :label: eqn_q_e
	
    \mathbf{q}^e = \begin{bmatrix}
        u_1 & w_1 & \theta_1 & u_2 & w_2 & \theta_2
    \end{bmatrix}^\intercal,

which are related to the interpolated displacements of the element
:math:`\mathbf{u}^e(x,t) = \left( u^e(x,t), w^e(x,t), \theta^e(x,t) \right)`
via the linear shape functions. 


Strain measure
~~~~~~~~~~~~~~

Based on these definitions, the continuous strains can be discretized into elemental strain values as:

.. math::
    :label: eqn_strains_disc
	
	\begin{align}
		&e^e = \left( 1 + \dfrac{u_2 - u_1}{L^e}\right)\cos{\theta^e} + \left(\dfrac{w_2 - w_1}{L^e}\right)\sin{\theta^e} - 1, \notag \\
		&\gamma^e = \left(\dfrac{w_2 - w_1}{L^e}\right)\cos{\theta^e} - \left( 1 + \dfrac{u_2 - u_1}{L^e}\right)\sin{\theta^e}, \\
		&\kappa^e = \dfrac{\theta_2 - \theta_1}{L^e}, \notag
	\end{align}
	
with :math:`e^e` and :math:`\gamma^e` related to the axial and transverse strains, respectively,
and :math:`\kappa^e` the curvature. 

Mass matrix
~~~~~~~~~~~

Next, discretization and evaluation of the different virtual work terms (recall the principle of virtual work: $\delta W_{\mathrm{inertial}} + \delta W_{\mathrm{internal}} = \delta W_{\mathrm{external}}$) leads to the definitions of the elementary vectors and matrices. Evaluation of the discretized inertial work component yields the expression for the elementary mass matrix:

.. math::
    :label: eqn_mass
	
    \mathbf{M}^e = \dfrac{\rho L^e}{6}\begin{bmatrix} 2A & 0 & 0 & A & 0 & 0 \\ 0 & 2A & 0 & 0 & A & 0\\ 0 & 0 & 2I & 0 & 0 & I \\ A & 0 & 0 & 2A & 0 & 0 \\ 0 & A & 0 & 0 & 2A & 0 \\ 0 & 0 & I & 0 & 0 & 2I\end{bmatrix},
	
with :math:`A` the cross-sectional area, :math:`\rho` the density, and :math:`I` the second moment of area. 


Elementary internal force vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The elementary internal force vector is likewise obtained by evaluating the discretized internal work component:

.. math::
    :label: eqn_internal_force_vec
	
    \mathbf{f}^e_\mathrm{int} = EA \bar{e} \begin{bmatrix} -\cos{\bar \theta} \\ -\sin{ \bar  \theta} \\ \bar \gamma \frac{L^e}{2} \\ \cos{\bar  \theta} \\ \sin{\bar  \theta} \\ \bar \gamma \frac{L^e}{2}\end{bmatrix} + kGA\bar \gamma \begin{bmatrix} \sin{\bar \theta} \\ -\cos{\bar \theta} \\  -\frac{L^e}{2}(1+\bar e) \\ -\sin{\bar \theta} \\ \cos{\bar \theta} \\ -\frac{L^e}{2}(1 + \bar e) \end{bmatrix} + EI \bar \kappa \begin{bmatrix} 0 \\ 0 \\ -1 \\ 0 \\ 0 \\ 1 \end{bmatrix},
	
where the :math:`\Bar{(\; )}` indicates a term evaluated at :math:`\Bar{\theta} = (\theta_1 + \theta_2)/2` and with :math:`E` the Young's modulus, :math:`G` the shear modulus and :math:`k` the shear correction factor. It is readily seen that the geometrical nonlinearities (the sine and cosine functions) appear in this internal force vector.


Equation of motion
~~~~~~~~~~~~~~~~~~

After assembly of the elemental matrices into their global counterparts according to traditional finite element procedures, the discretized nonlinear finite element model can be written in the general form:

.. math::
    :label: eqn_FEM_model

	\mathbf{M}\mathbf{\Ddot{q}} + \mathbf{C}\mathbf{\Dot{q}} + \mathbf{f}_{\mathrm{int}}(\mathbf{q}) = \mathbf{f}_{\mathrm{ext}},

with :math:`\mathbf{q}`, :math:`\mathbf{M}` and :math:`\mathbf{f}_\mathrm{int}`
the degree-of-freedom vector, the mass matrix and the internal force vector, respectively,
for the entire structure, all expressed in a global reference frame.
Note that a damping matrix :math:`\mathbf{D}` can be added into the model
depending on the desired simulation;
in the case of the forced response simulations (see below),
a linear Rayleigh-type damping is considered, proportional to the mass.
Additionally, the external force vector :math:`\mathbf{f}_\mathrm{ext}` is evaluated
according to any applied forces :math:`n(x,t)`, :math:`p(x,t)` and :math:`q(x,t)`,
representing, respectively, the distributed axial force, transverse force and bending moment.


Searching for periodic response
-------------------------------

Different simulations are possible depending on how Eq. :eq:`eqn_FEM_model` is written.
Here, two different nonlinear computations are proposed to the user,
namely the periodic response of the system under harmonic forcing (forced response)
and the periodic response of the free and undamped system (nonlinear modes).

Forced response
~~~~~~~~~~~~~~~

Consider first the periodic response of the system under harmonic forcing.
In this case, a harmonic force vector is applied for :math:`\mathbf{f}_\mathrm{ext}`:

.. math::
    :label: eqn_FEMforced

	\mathbf{M}\Ddot{\mathbf{q}}(t) + \mathbf{C}\Dot{\mathbf{q}}(t) + \mathbf{f}_\mathrm{int}[\mathbf{q}(t)] = \mathbf{F} \sin{\Omega t},
	
with :math:`\mathbf{F}` the vector of harmonic forcing amplitudes and
:math:`\Omega` the harmonic forcing frequency.
The amplitude of the response
(e.g. the maximum transverse displacement,
the amplitude of a single harmonic of a certain degree of freedom, etc.)
as a function of the response frequency yields the traditional forced response curves.


Nonlinear modes
~~~~~~~~~~~~~~~

Next we consider the second type of nonlinear computation,
that of the periodic response of the free and undamped system, equivalent to the nonlinear modes.
The nonlinear modes are graphically represented as the "backbone curve"
in the same frequency-amplitude plot as the forced response,
which is the curve representing the locus of all resonance points.
Within the nonlinear dynamics community, there are different definitions of the nonlinear modes,
aimed at extending the concept of the linear mode to the nonlinear regime.
For more insight on the definitions and usefulness of nonlinear modes,
see `(G. Kerschen, M. Peeters, et. al., 2009) <https://www.sciencedirect.com/science/article/pii/S0888327008001015>`_, `(C. Touzé, O. Thomas and A. Chaigne, 2004) <https://www.sciencedirect.com/science/article/pii/S0022460X03010083>`_ or `(S. Shaw, C. Pierre, 1991) <https://hal.archives-ouvertes.fr/hal-01310674>`_.

There are two methods proposed for computing the nonlinear modes.
The first represents the true physical definition of the nonlinear modes
and can be thought of as the reference backbone curve.
The second method makes use of a concept known as phase resonance to compute the same backbone curve.


Nonlinear normal modes computation
""""""""""""""""""""""""""""""""""

Recalling that the nonlinear modes can be thought of in a physical sense as
the periodic response of the free and undamped nonlinear system.
Therefore, Eq. :eq:`eqn_FEM_model` takes the form:

.. math::
    :label: eqn_FEMfree
	
    \mathbf{M}\Ddot{\mathbf{q}}(t)  + \mathbf{f}_\mathrm{int}[\mathbf{q}(t)] = \mathbf{0},

where it can be seen that the damping term :math:`\mathbf{D}\mathbf{\Dot{q}}`
and the forcing term :math:`\mathbf{f}_\mathrm{ext}` have been removed.
The backbone curve of a particular mode is traced when plotting the amplitude of
the response as a function of the oscillation frequency
(in fact, the backbone curve can be overlaid onto the forced response curves
to show the intersection with the points of nonlinear resonance, which is often done in the literature).


Phase resonance computation
"""""""""""""""""""""""""""

Another method for the computation of the nonlinear modes consists in
creating a phase resonance between the applied external forcing
and the displacement of the system :math:`\mathbf{u}`.
A detailed explanation is beyond the scope of this short summary,
but, in short, it has been shown that at a phase difference of :math:`\pi/2`
between the external forcing and the response :math:`\mathbf{u}`,
the external forcing term :math:`\mathbf{f}_\mathrm{ext}` exactly equals and,
therefore, cancels the damping term :math:`\mathbf{D}\mathbf{\Dot{q}}`,
leading to a simulation mathematically equivalent to Eq. :eq:`eqn_FEMfree`
(see e.g. `(M. Peeters, G. Kerschen and J. C. Golinval, 2011) <https://www.sciencedirect.com/science/article/pii/S0022460X10005559>`_).


Frequency domain resolution with MANLAB
---------------------------------------

The equation of motion for the simulation under consideration is solved in the frequency domain
by the solver MANLAB.
The computation is performed automatically without intervention on behalf of the user.
For this reason, only a high-level overview of MANLAB is presented here;
the user is encouraged to download the solver,
its documentation and some example test cases at
`(the MANLAB website) <http://manlab.lma.cnrs-mrs.fr/>`_.

The formalism of MANLAB requires that the system contain only
polynomial nonlinearities of quadratic order or less.
For this reason, an extra step known as the "quadratic recast"
of any non-quadratic nonlinearities
(in this case, notably the geometric nonlinearities :math:`\sin{\theta}`
and :math:`\cos{\theta}`) is required,
wherein additional variables are added into the system.
After this step, the FE equation of motion Eq. :eq:`eqn_FEM_model` takes the form of
a differential-algebraic system of equations (DAE).
This process is explained in detail in a recently-submitted paper,
`(M. Debeurre, A. Grolet, B. Cochelin and O. Thomas, submitted 2022) <https://hal.archives-ouvertes.fr/hal-03819580>`_.


Harmonic balance method
~~~~~~~~~~~~~~~~~~~~~~~

The unknowns of the system are solved in the frequency domain
using a combination of harmonic balance expansion and numerical continuation.
First, the harmonic balance method (HBM) is applied.
Each unknown :math:`x(t)` is assumed periodic
and is expanded in Fourier series up to harmonic :math:`H`:

.. math::
    :label: eqn_HBM
	
	x(t) = x_0 + \sum_{k=1}^H (x_{k}^c\cos{k\omega t} + x_{k}^s\sin{k\omega t}),
	
where :math:`x_0`, :math:`x_{k}^c` and  :math:`x_{k}^s` represent the Fourier coefficients
and :math:`\omega` the angular frequency.
Substituting the harmonic balance expansion Eq. :eq:`eqn_HBM`
into the quadratically-recast DAE equations governing the system,
the resulting system of (quadratic) algebraic equations can be written as:	

.. math::
    :label: eqn_residual

	\mathbf{R}(\mathbf{X},\omega,\lambda) = \mathbf{R}(\tilde{\mathbf{X}}) =  0,  

with :math:`\mathbf{X}` the vector containing the Fourier coefficients of all variables
and :math:`\lambda` a continuation parameter.
Eq. :eq:`eqn_residual` is the one solved by numerical continuation.


Asymptotic numerical continuation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The asymptotic numerical method (ANM) is used as the strategy for
numerical continuation of periodic solutions.
The details of this method can be found in many references
(see e.g.
`(B. Cochelin, N. Damil and M. Potier-Ferry, 1994) <https://onlinelibrary.wiley.com/doi/10.1002/nme.1620370706>`_ or
`(B. Cochelin and C. Vergez, 2009) <https://www.sciencedirect.com/science/article/pii/S0022460X09001217?via%3Dihub>`_).
Globally, the ANM seeks the solution in a power series expansion of
the unknowns :math:`\tilde{\mathbf{X}}` around :math:`a`, a pseudo arc-length parameter.
Since the process is entirely automated in MANLAB, a detailed investigation is left to the interested user.


