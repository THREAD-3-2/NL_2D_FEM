.. _walkthrough:

Complete walkthrough of a simple example
========================================

This example outlines the step-by-step procedure to compute the benchmark reference solution of
a clamped-clamped beam with a point periodic force,
found in `(A. Givois et al., 2019) <https://link.springer.com/article/10.1007/s11071-019-05021-6>`_.

Set the Matlab Path
-------------------

First run the :mod:`set_src_path` from the :mod:`NonLinearFEM` folder.


Create a new FE model
---------------------

Next, create a new folder.
Then, create a new script name :mod:`run_problem` that will contain all of the operations to be done.

In the script, clean the state of Matlab:

.. code-block::

	clear
	close all
	clc
	% dont forget to set the path (run set_src_path from the NonLinearFEM folder)
	
Then create a new FE model:

.. code-block::

	%% FEM model definition
	% create a new NL_2D_FEM model
	model = NL_2D_FEM; 
	

Input definition
----------------

The best way to define the inputs is to use a model data structure. However, we also describe how to build the model manually using the initialization methods of the :mod:`NL_2D_FEM` class.


Using a model data structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The model can be created by first building a model data structure. An example model is shown in the file :mod:`arthur_beam_problem`. Then the FE model is directly set up.

.. code-block::

	% load the model data structure
	mds = arthur_beam_problem();
	% set up the model using the data structure
	model = model.set_model_from_mds(mds);


Manual input definition
~~~~~~~~~~~~~~~~~~~~~~~

Geometry
""""""""

.. code-block::

	% Key nodes
	geom_node = [1 0 0 0;   % Point 1 coord (0,0)
				 2 L/pi 0 0; % Point 2 coord (L/pi,0) % for force 
				 3  L  0 0;];% Point 3 coord (L, 0)
	% Lines of the geometry
	geom_element = [1 1 2    % line 1 between point 1 and 2
					2 2 3];  % line 2 between point 2 and 3
	% Discretisation for each line
	discretisation = [15, 30]; % line 1 discretized with 15 elements
	% set geometry (optional if a mesh is directly provided)
	model = model.set_geom(geom_node, geom_element, discretisation);

Mesh
""""

.. code-block::
	
	% set compute a mesh 
	model = model.auto_mesh();
	
Properties
""""""""""

.. code-block::

	% Key nodes
	%% material inputs + adimentionalization
	b = 50; h = 1;  % constant relative to the beam section
	S = b*h;       % section area
	I = b*h^3/12;  % section second moment of area
	rho = 7.8e-9;  % material density
	E = 210e3;     % material Young Modulus
	poisson = 0.3;  % material Poisson ratio
	k = 1;           % shear area coeff
	alpha = 0.1;   % Raileight damping coeff (C = alpha M)
	% set properties
	model = model.set_prop(S, I, rho, E, poisson, k, alpha);
	
Boundary conditions
"""""""""""""""""""

.. code-block::

	%% Boundary condition input (Dirichlet condition only)
	mds.boundary.bc_node_list{1} = struct('node', 1,'dof', [1,2,3]); % node 1 is clamped
	mds.boundary.bc_node_list{2} = struct('node', 3,'dof', [1,2,3]); % node 3 is clamped
	% set boundary condition
	model = model.set_boundary(bc_node_list);


Visualization
"""""""""""""

.. code-block::
	
	%% Visualisation node Input (for results display)
	visu_node_list{1} = struct('node', 2 ,...
							'dof', [1 2]);
	% set visualized nodes
	model = model.set_visu(visu_node_list); 


Force definitions
"""""""""""""""""

.. code-block::

	% point periodic force
	periodic_ponctual_force_node_list{1} = struct('node', 2,'dof', [2],'amplitude', [0.1], 'harmonic', [1+1i] ); % complex amplitude f = re(amp) cos + im(amp) sin
	% dynamic loads
	model = model.set_periodic_loads('ponctual', periodic_ponctual_force_node_list);


Matrices and force vector initialization
----------------------------------------

.. code-block::

	% assemble mass matrix and force vector
	model = model.initialise_matrices_and_vector();


Static solution
---------------

.. code-block::

	% static solution
	[qs_full, res] = model.solve_static_problem();
	q0_full = model.vectors.null_vector;
	fig = figure(99); fig.Name='Static configuration'
	model.plot_deformed_mesh(q0_full, fig, '--k') % undef mesh
	model.plot_deformed_mesh(qs_full, fig, '-c') % def mesh
	% compute stresses
	[strain,stress] = model.strains_and_stress_at_gauss_point(qs_full);


Modal analysis
--------------

.. code-block::

	% modal analysis
	[shape, freq] = model.linear_modal_analysis(qs_full);
	% [shape, freq] = model.linear_modal_analysis();
	fig2 = figure(98); fig2.Name='Mode Shapes'
	model.plot_deformed_mesh(q0_full, fig2, '--k')
	model.plot_deformed_mesh(qs_full, fig2, '-c')
	model.plot_deformed_mesh(shape(:,1), fig2, '-r')
	model.plot_deformed_mesh(shape(:,2), fig2, '-g')
	model.plot_deformed_mesh(shape(:,3), fig2, '-b')


Linear forced analysis
----------------------

.. code-block::

	% linear analysis
	H = 1;
	target_mode = 1;
	Omega = linspace(freq(target_mode)*0.8, freq(target_mode)*1.2,500)*2*pi;
	[qp_full, bode] = model.linear_analysis(H, Omega, qs_full);
	figure
	subplot(2,1,1) % first harmonic amplitude; hold on
	plot(Omega, bode.amp_qp_full{1}(4,:)) % u node 2
	plot(Omega, bode.amp_qp_full{1}(5,:)) % v node 2
	xlabel('Omega'); ylabel('Amp H1')
	subplot(2,1,2) % first harmonic phase; hold on
	plot(Omega, bode.phase_qp_full{1}(5,:)) % phase u node 2
	plot(Omega, bode.phase_qp_full{1}(5,:)) % phase v node 2
	xlabel('Omega'); ylabel('Phase H1')


MANLAB analysis
---------------

MANLAB inputs
~~~~~~~~~~~~~

Define the MANLAB input data and initialize the MAN system:

.. code-block::

	%% MANLAB LAUNCHING SEQUENCE
	%% MANLAB INPUTS
	global U Section Diagram   % Global variables to export point from the diagram.
	H = 10;          % number of harmonics for the Fourier series   
	type = 'autonomous'; % type of system (can be 'forced' or 'autonomous')
	% type = 'forced'; % type of system (can be 'forced' or 'autonomous')
	target_mode = 1; % modeto be studied in NNM (autonomous) or in FRF (forced)
	angfreq = 'omega'; % 'omega' or constant value 
	%% Use the model to initialise MANLAB computation automatically
	% MANLAB structure of parameters for equation.m
	[nz, nz_aux, parameters] = model.set_MAN_parameters(H, type,  model, angfreq);
	% Construct MANLAB system (matlab object)
	sys = SystHBQ(nz,nz_aux,H,@equations_vector_NL_2D_FEM,@point_display,@global_display,parameters,type,'vectorial');


Nonlinear normal modes starting point
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Find the starting point for the nonlinear modes MANLAB computation:

.. code-block::

	omega0 = (freq(target_mode)*2*pi);
	lambda0 = 0;
	idx = sys.getcoord('cos',2 ,1); % dof to be imposed amplitude
	amp = 1e-5;    % imposed amplitude
	[Z0] = model.man_initial_point(H, omega0, qs_full, amp*shape(:,target_mode));
	U0 = sys.init_U0(Z0, omega0, lambda0);
	U0 = model.solve_MAN_system_at_fixed_amplitude(U0, idx, amp, sys); 


Forced response starting point
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Find the starting point for the forced response MANLAB computation:

.. code-block::

	omega0 = freq(target_mode)*2*pi*0.8;
	lambda0 = omega0; % continuation parameter initial value
	[qp_full, bode] = model.linear_analysis(H, omega0);    
	[Z0] = model.man_initial_point(H, omega0, qs_full, qp_full);
	U0 = sys.init_U0(Z0, omega0, lambda0);
	U0 = model.solve_MAN_system_at_fixed_frequency(U0, omega0, sys);


Display variables and call to MANLAB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Choose the display variables visualized during the computation and call MANLAB:

.. code-block::

	%%% Variable displayed in the projected bifurcation diagram.
	% To plot the coefficient of cos(h omega t) of variable number i with
	% respect to lambda you should write as follows:
	dispvars = [sys.getcoord('omega') sys.getcoord('cos',1,1);
				sys.getcoord('omega') sys.getcoord('sin',1,1);
				sys.getcoord('omega') sys.getcoord('cos',2,1);
				sys.getcoord('omega') sys.getcoord('sin',2,1)];
	%% Launch of Manlab with options
	Manlab('sys'       ,sys , ...
		'U0value'         ,U0, ...
		'order'           ,20, ...     % order of the series
		'ANMthreshold'    ,1e-10, ...   % threshold for the domain of validity of the series
		'Amax_max'        ,1e2, ...    % maximum value of the domain of validity of the series
		'NRthreshold'     ,1e-12, ...   % threshold for Newton-Raphson (NR) corrections
		'NRitemax'        ,50, ...     % Maximum number of iteration of NR algorithm
		'NRstart'         ,0, ...      % NR corrections for the user-defined starting point [on]/off
		'NRmethod'        ,0, ...      % NR corrections on/[off]
		'BifDetection'    ,1, ...      % Detection of bifurcation [on]/off
		'PointDisplay'    ,0, ...      % Point display [on]/off
		'GlobalDisplay'   ,0, ...      % Global display [on]/off
		'StabilityCheck'  ,0, ...      % Stability computation on/[off]
		'StabTol'         ,1e-6, ...   % Stability tolerance
		'displayvariables',dispvars);     % MANLAB run


Quick launch of a computation
-----------------------------

In what follows, all of the previous elementary functions have been used to provide a quick way to start a MANLAB computation.


Define the model
~~~~~~~~~~~~~~~~

.. code-block::

	clear
	close all
	clc
	%% Path of the SRC file.
	addpath(genpath('..\..\..\SRC'));
	addpath(genpath('..\..\NonlinearDyn_vectorialform'));
	%% FEM model definition
	% load the model data structure
	mds = toy_house_problem();
	% create a new NL_2D_FEM model
	model = NL_2D_FEM; 
	% set up the model using the data structure
	model = model.set_model_from_mds(mds);


Initialize the computation
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block::

	%% MANLAB LAUNCHING SEQUENCE
	%% MANLAB INPUTS
	global U Section Diagram   % Global variables to export point from the diagram in GUI
	H = 10;          % number of harmonics for the Fourier series   
	type = 'autonomous'; % type of system (can be 'forced' or 'autonomous')
	% type = 'forced'; % type of system (can be 'forced' or 'autonomous')
	target_mode = 1; % modeto be studied in NNM (autonomous) or in FRF (forced)
	angfreq = 'omega'; % 'omega' or constant value 
	% MANLAB structure of parameters for equation.m
	[nz, nz_aux, parameters] = model.set_MAN_parameters(H, type,  model, angfreq);
	% Construct MANLAB system (matlab object)
	sys = SystHBQ(nz,nz_aux,H,@equations_vector_NL_2D_FEM,@point_display,@global_display,parameters,type,'vectorial');
	%% compute static equilibrium, modal analysis and the MANLAB starting point
	[U0, omega0, lambda0] = model.initialise_MAN_computation(sys, type, target_mode);


Call to MANLAB
~~~~~~~~~~~~~~

Same as in the detailed version.
