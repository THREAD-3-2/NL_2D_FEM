.. _useful:

Useful functions folder
=======================

The folder :mod:`USEFUL_FUNCTIONS` contains several functions used in the computation:

* :mod:`create_empty_mds`: Define an empty :mod:`mds`.
* :mod:`harmonicDecomposition`: Used in the point display method to visualize the harmonic content at a given point of the MANLAB diagram.
* :mod:`animate_solution`: Used in the point display method to animate the mesh at a given point of the MANLAB diagram.
* :mod:`man_residual_fixed_amplitude`: Residual used in the initialization.
* :mod:`man_residual_fixed_frequency`: Residual used in the initialization.
* :mod:`plot_deformed_mesh`: Same as in :mod:`model.plot_deformed_mesh`, for convenience.
	
In addition, a tool function :mod:`set_src_path_for_NL_2D_FEM` is provided to set the Matlab path.
It has to be run from the :mod:`NonLinearFEM` folder, before running any example or computations.
