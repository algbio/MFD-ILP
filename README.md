# Fast and exact ILP-based solvers for the Minimum Flow Decomposition problem, and variants of it

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
 For each solver, examples of the inputs are available in `Example` folder.
 
 ## 1.3 Different Formulations
There are three different solvers available: the "Standard" files corresponds to the original and standard formulation; the "Inexact" files corresponds to the original formulations adjusted to incorporate inexact flow constraints and the "Subpath" files corresponds to the original formulations with the addition to the subpath constraints.
 
 ## 1.4 Running the solvers
 For each solvers, in order to run each formulation, open the respective notebook and change the variable $path$ in the last cell to the folder where all the input files are. For the subpath constraints formulation, also change the $number_paths$ to the appropriated amount. The default value is 4.

As reminder, all the input files are in [Catfish](https://github.com/Kingsford-Group/catfishtest) format. See folder `Example` for sample of inputs.

 ## 1.5 Outputs
 Each solver outputs two files: the first file called `results_[CPLEX or Gurobi].txt` contains  the optimal number of $k$ flow paths and the runtime required to solve such instance, each instance is displayed in a single line; the second file called `results_[CPLEX or Gurobi]-details.txt` contains in each line the corresponding value of $w_k$ and the $k$ flow path associated with that solution, different instances are separated by "------------". 

# 2. Stand-alone solver for standard MFD

## 2.1 Running

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

**NOTE 1**: Check `standalone/example1.graph` and `standalone/example2.graph` for an example input graphs. Note that, as opposed to the Jupyter notebooks, the stand-alone solver cannot read more than one graph from the input file. Encode only a single graph in the input file!

**NOTE 2**: This graph format does not support parallel edges. If your graph has such edges, subdivide them (i.e. replace them with a path of two edges).

Example usage:

```
python3 standalone/mfd-solver-gurobi.py -i standalone/example2.graph -o standalone/example2.out
```

### 2.2 Example input / output for the standard formulation:

The flow in the above figure (left) can be encoded as (`6` is the number of nodes):

```
6
s a 6
s b 7
a b 2
a c 4
b c 9
c d 6
c t 7
d t 6
```

Its minimum flow decomposition in the figure (right) will be output as:

```
4 ['s', 'a', 'c', 'd', 't']
2 ['s', 'a', 'b', 'c', 'd', 't']
7 ['s', 'b', 'c', 't']
```

## 3 Stand-alone solver of MFD variants

### 3.1 For inexact formulation:

Run the solver as:

```
python3 mfd-inexact-solver-gurobi.py -i INPUT -o OUTPUT [-wt WEIGHTTYPE] [-t THREADS]

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

### 3.2 Example input / output for the inexact formulation:

The flow in the above figure (left) can be encoded as (`6` is the number of nodes):

```
6
s a 4 5
s b 6 8
a b 2 3
a c 3 5
b c 8 10
c d 5 7
c t 6 8
d t 5 7
```

Its minimum flow decomposition in the figure (right) will be output as:

```
2 ['s', 'a', 'b', 'c', 'd', 't']
3 ['s', 'a', 'c', 'd', 't']
6 ['s', 'b', 'c', 't']
```

### 3.3 For subpath formulation:

Run the solver as:

```
python3 mfd-subpath-solver-gurobi.py -i INPUT -s SUBPATHS -o OUTPUT [-wt WEIGHTTYPE] [-t THREADS]

required arguments:
  -i INPUT, --input INPUT
                        Input filename
  -s SUBPATHS, --input SUBPATHS
                        Subpaths filename
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

### 3.4 Example input / output for the subpath formulation:

The flow in the above figure (left) can be encoded as (`6` is the number of nodes):

```
6
s a 6
s b 7
a b 2
a c 4
b c 9
c d 6
c t 7
d t 6
```
With subpath files:

```
2
s a b 
a c
```

Its minimum flow decomposition in the figure (right) will be output as:

```
2 ['s', 'a', 'b', 'c', 'd', 't']
3 ['s', 'a', 'c', 'd', 't']
6 ['s', 'b', 'c', 't']
```

## 4 Installing Gurobi

Download the solver from [www.gurobi.com](www.gurobi.com), activate the (academic) license as instructed, and then install the Python API with:

```
pip3 install gurobipy
```
