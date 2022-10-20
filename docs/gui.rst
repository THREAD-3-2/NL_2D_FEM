.. _gui:

Code use
========

MANLAB Graphical User Interface (GUI)
-------------------------------------

When MANLAB is called, a Graphical User Interface (GUI) named "Manlab" appears along with the simulation windows. There are several different options available to the user. This section aims to explain the most important options that the user is most likely to make use of.

Nine different blocks house the GUI features. A basic description of their functionality is provided below:

* **Continuation**: The options housed in this block relate to the ANM continuation itself (if not specified, default values are used).
* **Correction**: There are two different types of Newton-Raphson (NR) corrections used to bring the continuation computation below the tolerance level indicated in "Threshold." The "Type" can be set to 0, 1 or 2 to implement no NR corrections, NR method 1 or NR method 2, respectively.
* **Bifurcation-Stability**: Checking the corresponding box in this block activates the bifurcation detector or stability analysis features of MANLAB.
* **Perturbation** This block enables the user to apply a perturbation to the system.
* **Display**: The options within this block allow the user to visualize the behavior of the system at a specific continuation point ("Point display" + "Select point") or for the entire computation ("Global display" + "Display diagram").
* **Current point**: The options in this block control the current point in the continuation; the current point can be chosen by the user, while other options allow the user to jump to a new point or reverse the calculation.
* **Diagram**: The "diagram" terminology in MANLAB refers to the Matlab structure that contains all of the computation information, including the computed solution branches, filled during ANM continuation. The options in this block allow the use to delete certain branches or parts of the computation, save the computation, or load a previously-saved computation.
* **Export**: The features contained in this block allow the user to export into the Matlab Workspace the computation at a single point ("Point"), an entire continuation branch ("Section") or the full computation ("Diagram").
	
Of the GUI features mentioned above, the following outlines in more detail those that the user is most likely to interact with:

* **Continuation** - Forward :math:`>>`: This represents the main button in the MANLAB GUI as the branches of solution are computed only when "Forward" is clicked.
* **Continuation** - Steps: This feature controls the number of solution branches computed by the ANM for each "Forward."
* **Current point** - Reverse tangent: This feature allows the user to reverse the direction of the continuation.
* **Current point** - Set: This feature sets the cursor at a (computed) continuation point selected by the user.
* **Current point** - Jump: This feature allows the user to indicate a location outside of calculated branches where they believe a solution exists. MANLAB uses Newton-Raphson corrections to converge to the nearest solution (if it exists). This feature is most often used when a point of internal resonance is encountered, in order to jump over the point of internal resonance.
	
As the branches of solution are computed by MANLAB, one or more figures (depending on if "Point display" and/or "Global display" are selected) show the evolution of the solution. After the number of branches indicated in "Steps" is reached, the solver (and all figures attached to the computation) will pause at their current location until "Forward" is clicked again. 
