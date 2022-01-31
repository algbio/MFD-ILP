# ILP solvers for Minimum Flow Decomposition
This repository contains solvers for the Minimum Flow Decomposition (MFD) problem, and variants of it, described in the paper:

> Fernando H. C. Dias, Lucia Williams, Brendan Mumey, Alexandru I. Tomescu, **Fast, Flexible, and Exact Minimum Flow Decompositions via ILP**, To appear in the Proceedings of **RECOMB 2022** - 26th Annual International Conference on Research in Computational Molecular Biology, 2022, [**Full version**](https://arxiv.org/abs/2201.10923).

In the MFD problem, we are given a flow in a directed acyclic graph (DAG) with unique source *s* and unique sink *t*, and we need to decompose it into the minimum number of weighted paths (usually the weights are positive integers) from *s* to *t*, such that the weights of the paths sum up to the flow values, for every edge. In the image below, the flow is decomposable into 3 weighted paths, and this number is minimum.

![MFD Example](https://github.com/algbio/MFD-ILP/raw/main/mfd-example.png) 

For all MFD variants described in the paper above, we provide Jupyter notebooks (folder `notebooks`) implemented using Python with the API for two different linear programming solvers: CPLEX and Gurobi.

In addition, for standard MFD we provide a slightly more efficient stand-alone solver (folder `standalone`), which uses the Gurobi API.

# 1. Jupyter notebooks

## 1.1 Requirements

Python:
  - itertools
  - more_itertools
  - math
  - os 
  - networkx 
  
 CPLEX Python API or Gurobi Python API (both version of the codes are available)
 Jupyter Notebook
 
 ## 1.2 Inputs
 For each solver, an example of the inputs are available in "Example" folder. 
 
 ## 1.3 Different Formulations
 There are three different solvers available: the "Standard" files corresponds to the original and standard formulation; the "Inexact" files corresponds to the original formulations adjusted to incorporate inexact flow constraints and the "Subpath" files corresponds to the original formulations with the addition to the subpath constraints.
 
 ## 1.4 Running the solvers
 In each solvers, in order to run each formulation, open the respective notebook and change the variable $path$ in the last cell to the folder where all the input files are. For the subpath constraints formulation, also change the $number_paths$ to the appropriated amount. The default value is 4.

As reminder, all the input files are in Catfish format. See folder "Example" for sample of inputs.

 ## 1.5 Outputs
 Each solvers outputs 2 files: the first file called "results_[CPLEX or Gurobi].txt" contains  the optimal number of $k$ flow paths and the runtime required to solve such instance, each instance is displayed in a single line; the second file called "results_[CPLEX or Gurobi]-details.txt" contains in each line the corresponding value of $w_k$ and the $k$ flow path associated with that solution, different instances are separated by "------------". 

 # 2. Stand-alone solver for standard MFD

Run the solver as:

```
python3 mfd-solver-gurobi.py -i INPUT -o OUTPUT [-wt WEIGHTTYPE] [-t THREADS]

required arguments:
  -i INPUT, --input INPUT
                        Input filename
  -o OUTPUT, --output OUTPUT
                        Output filename

optional arguments:
  -wt WEIGHTTYPE, --weighttype WEIGHTTYPE
                        Type of path weights (default int+):
                           int+ (positive non-zero ints), 
                           float+ (positive non-zero floats).
  -t THREADS, --threads THREADS
                        Number of threads to use for the Gurobi solver; use 0 for all threads (default 0).
```

**NOTE**: Check `standalone/example.graph` for an example input graph. Note that, as opposed to the Jupyter notebooks, the stand-alone solver cannot read more than one graph from the input file. Encode only a single graph in the input file!

Example usage:

```
python3 standalone/mfd-solver-gurobi.py -i standalone/example.graph -o standalone/example.out
```
