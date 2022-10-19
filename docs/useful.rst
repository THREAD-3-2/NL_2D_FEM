.. _useful:

Useful functions folder
=======================

The folder :mod:`USEFUL_FUNCTIONS` contains several functions used in the computation:

	\begin{itemize}
		\item :mod:`create_empty_mds`: Define an empty :mod:`mds`.
		\item :mod:`animate_solution`: Used in the point display method to animate the mesh at a given point of the MANLAB diagram.
		\item :mod:`man_residual_fixed_amplitude`: Residual used in the initialization.
		\item :mod:`man_residual_fixed_frequency`: Residual used in the initialization.
		\item :mod:`plot_deformed_mesh`: Same as in :mod:`model.plot_deformed_mesh`, for convenience.
	\end{itemize}
	
In addition, a tool function :mod:`set_src_path_for_NL_2D_FEM` is provided to set the Matlab path. It has to be run from the :mod:`NonLinearFEM` folder, before running any example or computations.