.. _inputs:

Description of the inputs
=========================

Input by defining a model data base
-----------------------------------

The easiest way to define the data of a :mod:`NL_2D_FEM` object is to fill the field of a Matlab structure called the {\bf model data structure} (:mod:`mds`). An empty :mod:`mds` can be created as :mod:`mds = create_empty_mds()`. Otherwise, an existing :mod:`mds` can be simply copied and edited.

The :mod:`NL_2D_FEM` model :mod:`model` can be initialized automatically with the method :mod:`model = model.set_model_from_mds(mds)`.

Geometry input
~~~~~~~~~~~~~~

There are two ways to define the geometry:

	\begin{itemize}
		\item Provide a list of points and lines that will be discretized automatically.
		\item Provide a mesh directly by giving the coordinates of the nodes and the connectivity table.
	\end{itemize}
	
\paragraph{Global geometry input}

The field :mod:`geom` of the :mod:`mds` governs the global geometry input. There are three sub-fields to be defined:

	\begin{enumerate}
		\item :mod:`mds.geom.geom_node`: List of number and coordinates of the key points of the geometry. Each line is of the type `point number, x coordinate, y coordinate'.
		\item :mod:`mds.geom.geom_element`: List of lines of the geometry. Each line is of the type `line number, point 1, point 2'.
		\item :mod:`mds.geom.discretisation`: List of number of elements in each line (list of integers, with size the number of lines of the geometry).
	\end{enumerate}
	
If the geometry is inputted, the mesh is generated automatically.

\paragraph{User-defined mesh input}

Alternatively, the geometry input can be left empty, if the mesh is provided directly. The field :mod:`mesh` of the :mod:`mds` governs the mesh. There are two sub-field to be defined:

	\begin{enumerate}
		\item :mod:`mds.mesh.nodes`: List of number and coordinates of the nodes of the mesh. Each line is of the type `node number, x coordinate, y coordinate'.
		\item :mod:`mds.mesh.connect`: Connectivity table; list of elements of the mesh. Each line is of the type `element number, node 1, node 2'.
	\end{enumerate}
	
Properties input
~~~~~~~~~~~~~~~~

The field :mod:`prop` of the :mod:`mds` governs the cross-section profile and the material properties. The three sub-fields to be defined are:

	\begin{enumerate}
		\item :mod:`mds.prop.S`: Area of the cross-section (constant for all elements).
		\item :mod:`mds.prop.I`: Second moment of area of the cross-section (constant for all elements).
		\item :mod:`mds.prop.k`: Shear coefficient of the cross-section (constant for all elements).
		\item :mod:`mds.prop.rho`: Material density (constant for all elements).
		\item :mod:`mds.prop.E`: Material Young's modulus (constant for all elements).
		\item :mod:`mds.prop.poisson`: Material Poisson's ratio (constant for all elements).
		\item :mod:`mds.prop.alpha`: Material damping coefficient (linear Rayleigh damping, constant for all elements).
	\end{enumerate}

The shear modulus :math:`G` is computed automatically from the data.

Visualization input
~~~~~~~~~~~~~~~~~~~

The field :mod:`visu` of the :mod:`mds` governs the visualization input. There is only a single field, representative of a cell of structures:

	\begin{enumerate}
		\item :mod:`mds.prop.visu.visu_node_list`: List (cell array) of observed node(s) and dof. 
		\item :mod:`mds.prop.visu.visu_node_list{n}` is a structure with two fields:
		\begin{enumerate}
			\item :mod:`mds.prop.visu.visu_node_list{n}.node`: Node number.
			\item :mod:`mds.prop.visu.visu_node_list{n}.dof`: List of dof number(s).
		\end{enumerate}
	\end{enumerate}

Boundary conditions input
~~~~~~~~~~~~~~~~~~~~~~~~~

The field :mod:`boundary` of the :mod:`mds` governs the (Dirichlet) boundary conditions. There is only a single field, representative of a cell of structures:

	\begin{enumerate}
		\item :mod:`mds.boundary.bc_node_list`: List (cell array) of observed node(s) and dof. 
		\item :mod:`mds.boundary.bc_node_list{n}` is a structure with two fields:
		\begin{enumerate}
			\item :mod:`mds.boundary.bc_node{n}.node`: Node number(s).
			\item :mod:`mds.boundary.bc_node{n}.dof`: List of dof number(s).
		\end{enumerate}
	\end{enumerate}
	
Force definition
~~~~~~~~~~~~~~~~

\subsubsection*{Static force definition}

The field :mod:`force.static` of the :mod:`mds` governs the definition of static forces. 

\paragraph{Static point forces}

The first forces to consider are static point loads, which are defined as:

	\begin{enumerate}
		\item :mod:`mds.force.static.static_ponctual_force_node_list`: List (cell array) of forced node(s) and dof. 
		\item :mod:`mds.force.static.static_ponctual_force_node_list{n}` is a structure with three fields:
		\begin{enumerate}
			\item :mod:`mds.force.static.static_ponctual_force_node_list{n}.node`: Node number(s).
			\item :mod:`mds.force.static.static_ponctual_force_node_list{n}.dof`: List of dof number(s).
			\item :mod:`mds.force.static.static_ponctual_force_node_list{n}.amplitude`: List of force amplitudes (one amplitude for each dof in the previous field).
		\end{enumerate}
	\end{enumerate}

\paragraph{Static distributed forces}

Next, we consider static distributed forces, which are defined as:

	\begin{enumerate}
		\item :mod:`mds.force.static.static_distributed_force_direction`: List of directions, chosen from among (1) x-axis, (2) y-axis  or (3) moment.
		\item :mod:`mds.force.static.static_distributed_force_amplitude`: List of amplitudes (linear force density). 
	\end{enumerate}

\subsubsection*{Periodic force definition}

Periodic loadings can also be applied to the structure. The field :mod:`force.periodic` of the :mod:`mds` governs the definition of periodic forces. 

\paragraph{Periodic point forces}

As before, consider first periodic point loads, defined as:

	\begin{enumerate}
		\item :mod:`mds.force.periodic.periodic_ponctual_force_node_list`: List (cell array) of forced node(s) and dof. 
		\item :mod:`mds.force.periodic.periodic_ponctual_force_node_list{n}` is a structure with three fields:
		\begin{enumerate}
			\item :mod:`mds.force.periodic.periodic_ponctual_force_node_list{n}.node`: Node number(s).
			\item :mod:`mds.force.periodic.periodic_ponctual_force_node_list{n}.dof`: List of dof number(s).
			\item :mod:`mds.force.periodic.periodic_ponctual_force_node_list{n}.amplitude`: List of force amplitudes (one amplitude for each dof in the previous field). Each element of the list is a complex number :math:`z`.
			\item :mod:`mds.force.periodic.periodic_ponctual_force_node_list{n}.harmonic`: List of harmonics (one harmonic for each dof in the previous field). 
		\end{enumerate}
		For each forced dof, the force is written as: :math:`f(t) = real(z) \cos(h t) + \real(z) \sin(h t)`.
	\end{enumerate}
	
\paragraph{Periodic distributed forces}

Next, consider the definition of periodic distributed forces, defined as:

	\begin{enumerate}
		\item :mod:`mds.force.static.static_distributed_force_direction`: List of directions, chosen from among (1) x-axis, (2) y-axis  or (3) moment.
		\item :mod:`mds.force.static.static_distributed_force_amplitude`: List of amplitudes, complex numbers (linear force density). 
		\item :mod:`mds.force.static.static_distributed_force_harmonic`: List of harmonics. 
	\end{enumerate}