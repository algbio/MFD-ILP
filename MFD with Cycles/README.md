# Minimum Flow Decomposition in Cyclic Graphs and its variants

This repository contains solvers for the Minimum Flow Decomposition (MFD) problem in cyclic graphs and variants of it described in the paper:

> Fernando H. C. Dias, Lucia Williams, Brendan Mumey, Alexandru I. Tomescu, **Minimal Flow Decomposition in Graphs with Cycles using Integer Linear Programming**, [**Full version**](https://arxiv.org/abs/2201.10923).

In the MFD problem, we are given a flow in a directed acyclic graph (DAG) with unique source *s* and unique sink *t*, and we need to decompose it into the minimum number of weighted elements (in this case, paths or cycles, walks and trails) (usually the weights are positive integers) from *s* to *t*, such that the weights of the elements sum up to the flow values, for every edge. 

In the image below, the flow is decomposable into 4 weighted paths and cycles in a), 3 trails in b), and 2 walks in c). For each decomposition, the number of elements is minimum.

<img src="https://github.com/algbio/MFD-ILP/raw/main/fd_cycles.png" width="900" height="550">

For all MFD variants described in the paper above, we provide Jupyter notebooks (folder `notebooks`) implemented using Python with the API for two different linear programming solvers: CPLEX and Gurobi.

In addition, for standard MFD, we provide a slightly more efficient standalone solver (folder `standalone`), which uses the Gurobi API.

# 1. Jupyter notebooks

## 1.1 Requirements

Python:
  - itertools
  - more_itertools
  - math
  - os 
  - networkx 
  
 CPLEX Python API or Gurobi Python API (both versions of the codes are available)
 Jupyter Notebook
 
 ## 1.2 Inputs
 For each solver, examples of the inputs are available in the `Example` folder.
 
 ## 1.3 Different Formulations
There are three different solvers available: the "FDPC" files correspond to the formulation where the elements are limited to paths or cycles; the "FDT" files correspond to the original formulations adjusted to obtain only trails, and the "FDW" files correspond to the formulations that admit only walks.

 ## 1.4 Running the solvers
 For each solver, to run each formulation, open the respective notebook and change the variable $path$ in the last cell to the folder where all the input files are.

As a reminder, all the input files are in [Catfish](https://github.com/Kingsford-Group/catfishtest) format. See folder `Example` for a sample of inputs.

 ## 1.5 Outputs
 Each solver outputs two files: the first file called `results_[CPLEX or Gurobi].txt` contains the optimal number of $k$ flow paths and the runtime required to solve such instance; each instance is displayed in a single line; the second file called `results_[CPLEX or Gurobi]-details.txt` contains in each line the corresponding value of $w_k$ and the $k$ flow path associated with that solution, different instances are separated by "------------". 

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
The dataset used to test these formulations are adapted from the available datasets in the literature. However, some alterations were necessary due to the specificities of our formulations. This section provides a complete detailed explanation of how such datasets were adapted and created. The dataset are named, according to their reference in the paper: **SRR020730 Salmon Adapted**, **Transportation Data** and **_E. coli_ Strains**.

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

## Dataset 2: Transportation Data

## Dataset 3: **_E. coli_ Strains**

We took 50 E. coli reference genomes from the dataset [3682 E. coli assemblies from NCBI](https://doi.org/10.5281/zenodo.6577996). We created a de Bruijn graph, a basic Bioinformatics data structure, from subsets of a number of *gt* genomes (*gt = 1...50*), as follows. We fix an integer *k=15*. A de Bruijn graph of order *k* built on a set of strings has a node for each k-mer, i.e. a string of length *k* in one of the strings (duplicated k-mers induce just a single node of the graph). We add an edge between two k-mers having a suffix-prefix overlap of (k-1) chracters if their merged (k+1)-mer appears in one of the strings. Since the genomes have repeated k-mers, then this graph will contain cycles.  Moreover, by construction, each of the *gt* genomes can be "traced" as a walk from its first k-mer to its last k-mer. We also add a global source *s* connected to all such start k-mers, and a global sink *t* connected from all such ending k-mers. 

To obtain flow values for the edges, we proceed as follows. First, for each of the *gt* genomes, we assign a weigth from log-normal distribution of mean -4 and variance 4, as in this [Bioinformatics package](http://alumni.cs.ucr.edu/~liw/rnaseqreadsimulator.html). Then, we initialize the flow values of all edge as *0*, and for each (k+1)-mer in a genome of weight *w*, we add up *w* to the flow vlue of the edge corresponding to that (k+1)-mer. Thus, each genome of weight *w* adds flow *w* to the *s*-*t* walk to which it corresponds in the graph. The resulting edge values are still a flow (statisfy flow conservation), and we can continue with the next genome.

Since E. coli genomes have more than 3 Million characters, and our *k* value is 15, we restrict the above procedure to non-overlapping windows of length *2000* inside each of the *gt* genomes. Moreover, to exclude trivial graphs, we don't consider graphs that have less than *50* simple cycles (computed with NetworkX function [simple_cycles](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.cycles.simple_cycles.html)). Note that the number of these simple cycles are not correlated with the number of cycles in a minimum flow decomposition, but are just indicative of the complexity of the graph.
