.. _object:

Initial value problems 
=====================================

The finite element (FE) information and methods are defined in a Matlab object called :mod:`NL_2D_FEM`.
An instance of the FE model (i.e. a new FE model) can be created using the default constructor:
:mod:`model = NL_2D_FEM`.
This creates an object :mod:`model` containing the data of
a given problem and methods to perform the data analysis.

The attributes and methods of the :mod:`NL_2D_FEM`  class (denoted :mod:`model`) are given below:

Attributes
----------

1. :mod:`model.geom`: The geometry of the model (key points and lines of the geometry)

    * :mod:`model.geom.geom_node`: List of number and coordinates of the key points of the geometry. Each line is of the type (point number, x coordinate, y coordinate).
    * :mod:`model.geom.geom_element`: List of lines in the geometry connecting the nodes. Each line is of the type (line number, point 1, point 2).
    * :mod:`model.geom.discretisation`: List of number of elements in each line (list of integers, with size the number of lines of the geometry).

2. :mod:`model.prop`: Element cross-section and material properties (constant properties for all elements of the model)

    * :mod:`model.prop.area`: Area of the cross-section.
    * :mod:`model.prop.inertia`: Second moment of area of the cross-section.
    * :mod:`model.prop.shear_coeff_k`: Shear coefficient of the cross-section.
    * :mod:`model.prop.density`: Material density.
    * :mod:`model.prop.young_mod`: Material Young's modulus. 
    * :mod:`model.prop.shear_mod`: Material shear modulus.
    * :mod:`model.prop.poisson`: Material Poisson's ratio.
    * :mod:`model.prop.alpha`: Material damping coefficient (linear Rayleigh damping).

3. :mod:`model.mesh`: mesh of the FE model

    * :mod:`model.mesh.nodes`: List of number and coordinates of the nodes of the mesh. Each line is of the type (node number, x coordinate, y coordinate).
    * :mod:`model.mesh.connect`: Connectivity table, \textit{i.e.} list of elements of the mesh. Each line is of the type (element number, node 1, node 2).
    * :mod:`model.mesh.number_nodes`: Total number of nodes in the mesh.
    * :mod:`model.mesh.number_elements`: Total number of elements in the mesh.

4. :mod:`model.visu`: Visualization information for displaying results

    * :mod:`model.visu.visu_node_list`: List (cell array) of observed node and degree of freedom (dof). 

        ** :mod:`model.visu.visu_node_list{n}` is a structure with two fields: a) :mod:`mds.visu.visu_node_list{n}.node`: Node number, and b) :mod:`mds.visu.visu_node_list{n}.dof`: List of dof number(s).

    * :mod:`model.visu.dof`: List of indices of observed dof in the assembled vector. 

5. :mod:`model.boundary`: Boundary conditions information

    * :mod:`model.boundary.bc_node_list`: List (cell array) of boundary condition nodes and dof.

        ** :mod:`mds.boundary.bc_node_list{n}` is a structure with two fields: a) :mod:`mds.boundary.bc_node{n}.node`: Node number, and b) :mod:`mds.boundary.bc_node{n}.dof`: List of dof number(s).

    * :mod:`model.boundary.dof_list`: List of indices of all dof in the assembled vector. 
    * :mod:`model.boundary.prescribed_dof`: List of indices of all prescribed dof (displacement or rotation set to zero).
    * :mod:`model.boundary.active_dof`: List of indices of all non-boundary condition dof.

6. :mod:`model.loads`: Loads (static and periodic, point force and distributed loads)

:mod:`model.loads.static`: Static loads.

        * :mod:`model.loads.static.static_ponctual_force_node_list`: List (cell array) of forced node(s) and dof. 

            ** :mod:`model.loads.static.static_ponctual_force_node_list{n}` is a structure with three fields: a) :mod:`model.loads.static.static_ponctual_force_node_list{n}.node`: Node number, b) :mod:`model.loads.static.static_ponctual_force_node_list{n}.dof`: List of dof number(s), and c) :mod:`model.loads.static.static_ponctual_force_node_list{n}.amplitude`: List of force amplitudes (one amplitude for each dof in the previous field).

        * :mod:`model.loads.static.static_distributed_force_direction`: List of directions among (1) x-axis (2) y-axis  or (3) moment.
        * :mod:`model.loads.static.static_distributed_force_amplitude`: List of amplitudes (linear force density). 

:mod:`model.loads.periodic`: Periodic loads 

        * :mod:`model.loads.periodic.periodic_ponctual_force_node_list`: List (cell array) of forced node(s) and dof. 

            ** :mod:`model.loads.periodic.periodic_ponctual_force_node_list{n}` is a structure with three fields: a) :mod:`model.loads.periodic.periodicponctual_force_node_list{n}.node`: Node number, b) :mod:`model.loads.periodic.periodic_ponctual_force_node_list{n}.dof`: List of dof number(s), and c) :mod:`model.loads.periodic.periodic_ponctual_force_node_list{n}.amplitude`: List of force amplitudes (one amplitude for each dof in the previous field).

        * :mod:`model.loads.periodic.periodic_distributed_force_direction`: List of direction among (1) x-axis (2) y-axis  or (3) moment.
        * :mod:`model.loads.periodic.periodic_distributed_force_amplitude`: List of amplitudes (linear force density). 


7. :mod:`model.matrices`: Assembled matrices of the model

    * :mod:`model.matrices.mass`: Mass matrix.
    * :mod:`model.matrices.damping`: Damping matrix.
    * :mod:`model.matrices.stiffness_at_origin`: Linear stiffness matrix (around undeformed (reference) configuration).
    * :mod:`model.matrices.tangent_stiffness_at_qs`: Unused.


8. :mod:`model.vectors`: Assembled vectors of the model 

    * :mod:`model.vectors.null_vector`: Assembled vector of dof full of zeros.
    * :mod:`model.vectors.static_forces`: Assembled vector of static forces.
    * :mod:`model.vectors.periodic_forces`: Assembled vector (matrix) of periodic forces (each row corresponds to a harmonic).
    * :mod:`model.vectors.static_solution`: Unused.


9. :mod:`model.solver`: Solver information

    * :mod:`model.solver.type`: :mod:`'homemade'`: Default, solve NL system with a simple Newton method, or :mod:`'fsolve'`: solve NL system with Matlab :mod:`fsolve` from optimization toolbox). 
    * :mod:`model.solver.tol_res`: Residual tolerance.
    * :mod:`model.solver.tol_step`: Step tolerance.
    * :mod:`model.solver.n_iter_max`: Number of maximum iterations for Newton's method.



Methods
-------

Initialization methods
~~~~~~~~~~~~~~~~~~~~~~

1. :mod:`model.set_model_from_mds`: Set the model attributes from a model data structure :mod:`mds`.
2. :mod:`model.set_geom`: Set the geometry of the model (key points and lines of the geometry).
3. :mod:`model.set_prop`: Set the cross-section and material properties.
4. :mod:`model.set_boundary`: Set the boundary conditions.
5. :mod:`model.set_visu`: Set the visualization information for displaying result.
6. :mod:`model.set_static_loads`: Set the static loads (point force and distributed loads).
7. :mod:`model.set_periodic_loads`: Set the periodic loads (point force and distributed loads).
8. :mod:`model.set_mesh`: Set mesh from user inputs.


FEM methods
~~~~~~~~~~~

Mesh-related methods
""""""""""""""""""""

1. :mod:`model.auto_mesh`: Create a mesh from the geometry information.
2. :mod:`model.mesh_from_truss`: Used in :mod:`auto_mesh`.


Element-related methods
"""""""""""""""""""""""

1. :mod:`model.elementOrientation`: Compute element length and orientation angle.
2. :mod:`model.transformationMatrix`: Compute rotation matrix from an orientation angle.
3. :mod:`model.elementary_mass_matrix`: Elementary mass matrix of the Timoshenko beam model.
4. :mod:`model.elementary_tangent_stiffness_matrix`: Elementary tangent stiffness matrix from a	given (elementary) dof vector.
5. :mod:`model.elementary_internal_force_vector`: Elementary vector of nonlinear internal forces from a given (elementary) dof vector.


Assembly-related methods
""""""""""""""""""""""""

1. :mod:`model.assemble_constant_matrix`: Assemble a matrix depending on the input (\textit{e.g.} mass, stiffness).
2. :mod:`model.assemble_constant_vector`: Assemble a vector depending on the input (\textit{e.g.} [static, periodic, internal] forces).
3. :mod:`model.initialise_matrices_and_vector`: Assemble the mass and stiffness matrices and the static and periodic force vectors.
	
Mechanics-related methods
~~~~~~~~~~~~~~~~~~~~~~~~~

1. :mod:`model.strains_and_stress_at_gauss_point`: Return stains, stress and gradient information from a given dof vector.
2. :mod:`model.solve_static_problem`: Solve the static problem.
3. :mod:`model.linear_modal_analysis`: Linear modal analysis around a given configuration.
4. :mod:`model.buckling_analysis`: Simple buckling analysis.
5. :mod:`model.linear_analysis`: Linear forced response.


MAN-related methods
~~~~~~~~~~~~~~~~~~~

1. :mod:`model.set_MAN_parameters`: Set the :mod:`parameters` structure.
2. :mod:`model.initialise_MAN_computation`: Automatic initialization of a MAN computation. Performs a static, modal anlysis around a static configuration and computes a MAN starting point (:mod:`U0`) depending on the type of computation.
3. :mod:`model.solve_MAN_system_at_fixed_amplitude`: Solve the MAN system such that the first harmonic of a given dof is fixed (variable frequency); used in the initialization method.
4. :mod:`model.solve_MAN_system_at_fixed_frequency`: Solve the MAN system such that the frequency is fixed (variable amplitude); used in the initialization method.
5. :mod:`model.man_initial_point`: Returns a MAN point from a given static and periodic vector of dof.
6. :mod:`model.man_auxiliary_variable_vector`: Computes the MAN auxiliary variables vector from a given vector of dof.
7. :mod:`model.auxiliary_variable_expression`: Used in the MAN equation to define the auxiliary variables.


Equation (residual)-related methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. :mod:`model.static_residual`: Used to solve the static problem.
2. :mod:`model.newton_solver`: Homemade Newton solver used to solve the static problem or the MAN initial point (or use the Matlab function :mod:`fsolve` from the optimization toolbox).


Display-related methods
~~~~~~~~~~~~~~~~~~~~~~~

1. :mod:`model.plot_deformed_mesh`: Plot deformed structure from a given vector.
2. :mod:`model.plot_static_configuration`: Plot the static configuration.
3. :mod:`model.plot_mode_shape`: Plot the mode shapes from a given list. 
