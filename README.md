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
 
 # Running the solvers
 In each solvers, in order to run each formulation, open the respective notebook and change the variable $path$ in the last cell to the folder where all the input files are. For the subpath constraints formulation, also change the $number_paths$ to the appropriated amount. The default value is 4.

As reminder, all the input files are in Catfish format. See folder "Example" for sample of inputs.

 # Outputs
 Each solvers outputs 2 files: the first file called "results_[CPLEX or Gurobi].txt" contains  the optimal number of $k$ flow paths and the runtime required to solve such instance, each instance is displayed in a single line; the second file called "results_[CPLEX or Gurobi]-details.txt" contains in each line the corresponding value of $w_k$ and the $k$ flow path associated with that solution, different instances are separated by "------------". 
 
