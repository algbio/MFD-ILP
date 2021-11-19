# MFD-ILP
Fast and exact ILP-based solvers for the Minimum Flow Decomposition (MFD) problem, and variants of it.

The solvers are implemented using Python with the API for two different linear programming solvers: CPLEX and Gurobi.

# Requirements

Python:
  - itertools
  - more_itertools
  - math
  - os 
  - networkx 
  
 CPLEX Python API or Gurobi Python API (both version of the codes are available)
 Jupyter Notebook
 
 # Inputs
 For each solvers, an example of the inputs are available in "Example" folder. 
 
 # Different Formulations
 There are three different solvers available: the "Standard" files corresponds to the original and standard formulation; the "Inexact" files corresponds to the original formulations adjusted to incorporate inexact flow constraints and the "Subpath" files corresponds to the original formulations with the addition to the subpath constraints.
 
 # Running the solvers;
 In each solvers, in order to run each formulation, open the respective notebook and change the variable $path$ in the last cell to the folder where all the input files are. As reminder, all the input files are in Catfish format. See folder "Example" for sample of inputs.
