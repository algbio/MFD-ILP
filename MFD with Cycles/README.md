# Minimum Flow Decomposition in Cyclic Graphs and its variants

This repository contains solvers for the Minimum Flow Decomposition (MFD) problem in cyclic graphs, and variants of it, described in the paper:

> Fernando H. C. Dias, Lucia Williams, Brendan Mumey, Alexandru I. Tomescu, **Minimal Flow Decomposition in Graphs with Cycles using Integer Linear Programming**, [**Full version**](https://arxiv.org/abs/2201.10923).

In the MFD problem, we are given a flow in a directed acyclic graph (DAG) with unique source *s* and unique sink *t*, and we need to decompose it into the minimum number of weighted elements (in this case paths or cycles, walks and trails) (usually the weights are positive integers) from *s* to *t*, such that the weights of the elements sum up to the flow values, for every edge. 

In the image below, the flow is decomposable into 4 weighted paths and cycles in a), into 3 trails in b) and into 2 walks in c). For eahc decomposition, the number of elements is minimum.

<img src="https://github.com/algbio/MFD-ILP/raw/main/fd_cycles.png" width="900" height="550">

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
There are three different solvers available: the "FDPC" files corresponds to the formulation where the elements are limited to paths or cycles; the "FDT" files corresponds to the original formulations adjusted to obtain only trails and the "FDW" files corresponds to the formulations that admit only walks.

 ## 1.4 Running the solvers
 For each solvers, in order to run each formulation, open the respective notebook and change the variable $path$ in the last cell to the folder where all the input files are.

As reminder, all the input files are in [Catfish](https://github.com/Kingsford-Group/catfishtest) format. See folder `Example` for sample of inputs.

 ## 1.5 Outputs
 Each solver outputs two files: the first file called `results_[CPLEX or Gurobi].txt` contains  the optimal number of $k$ flow paths and the runtime required to solve such instance, each instance is displayed in a single line; the second file called `results_[CPLEX or Gurobi]-details.txt` contains in each line the corresponding value of $w_k$ and the $k$ flow path associated with that solution, different instances are separated by "------------". 

# 2. Stand-alone solver for standard MFD

## 2.1 Running

Run the solver as:

```
python3  -i INPUT -o OUTPUT [-wt WEIGHTTYPE] [-t THREADS]

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

**NOTE 1**: 

**NOTE 2**: This graph format does not support parallel edges. If your graph has such edges, subdivide them (i.e. replace them with a path of two edges).

Example usage:

```
python3 
```

### 2.2 Example input / output for the standard formulation:

The flow in the above figure (left) can be encoded as (`6` is the number of nodes):

```

```

Its minimum flow decomposition in the figure (right) will be output as:

```

```

# 3. Dataset Creation
The dataset used to test these formulations are adapted from the available datasets in the literature. However, some alterations were necessary due to the specificities of our formulations. This section provides a complete detailed explanation of how such datasets were adapted and created. The dataset are named, according to their reference in the paper: **SRR020730 Salmon Adapted**, **Transportation Data** and ** *E. coli* Strains **.

## Dataset 1: SRR020730 Salmon Adapted
A specific sample from the [Catfish](https://github.com/Kingsford-Group/catfishtest) archive was selected for this dataset, and new cyclic graphs were created using the ground truth files. We took each ground truth for each graph and created a few copies of it. In each copy, we swapped two random nodes (except the source and sink). We built the graph using new and original ground truths, and that guarantees that at least one cycle is obtained.

For example, assume for a graph that the ground truth is:

```
4 ['s', 'a', 'c', 'd', 't']
2 ['s', 'a', 'b', 'c', 'd', 't']
7 ['s', 'b', 'c', 't']
```
By duplicating each ground truth and swapping each position randomly:
```
4 ['s', 'a', 'c', 'd', 't']
2 ['s', 'a', 'b', 'c', 'd', 't']
7 ['s', 'b', 'c', 't']
4 ['s', 'c', 'a', 'd', 't']
2 ['s', 'a', 'd', 'c', 'b', 't']
```
The resulting graph contains the original paths plus additional cycles.

This dataset is suitable for testing our FDPC formulation because it admits decompositions where the original paths are still obtained and new decompositions where cycles are obtained.

## Dataset 2: Transportation DAta

## Dataset 3: 


