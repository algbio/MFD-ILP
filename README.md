# Fast and exact ILP-based solvers for the Minimum Flow Decomposition problem, and variants of it

This repository contains solvers for the Minimum Flow Decomposition (MFD) problem and variants of it. For MFD in directed acyclic graphs, the formulation is described in the paper:

> Fernando H. C. Dias, Lucia Williams, Brendan Mumey, Alexandru I. Tomescu, **Fast, Flexible, and Exact Minimum Flow Decompositions via ILP**, To appear in the Proceedings of **RECOMB 2022** - 26th Annual International Conference on Research in Computational Molecular Biology, 2022, [**Full version**](https://arxiv.org/abs/2201.10923).

In the MFD problem, we are given a flow in a directed acyclic graph (DAG) with unique source *s* and unique sink *t*, and we need to decompose it into the minimum number of weighted paths (usually the weights are positive integers) from *s* to *t*, such that the weights of the paths sum up to the flow values, for every edge. In the image below, the flow is decomposable into 3 weighted paths, and this number is minimum.

![MFD in DAGS Example](https://github.com/algbio/MFD-ILP/raw/main/mfd-example.png) 

For MFD in cyclic graphs, the formulation can be found at:

> Fernando H. C. Dias, Lucia Williams, Brendan Mumey, Alexandru I. Tomescu, ** Minimum Flow Decomposition in Graphs with Cycles via Integer Linear Programming[**Full version**](https://arxiv.org/abs/2201.10923).


In the image below, the flow is decomposable into 4 weighted paths and cycles in a), into 3 trails in b) and into 2 walks in c). For eahc decomposition, the number of elements is minimum.

<img src="https://github.com/algbio/MFD-ILP/raw/main/fd_cycles.png" width="900" height="550">


For all MFD variants, we provide Jupyter notebooks implemented using Python with the API for two different linear programming solvers: CPLEX and Gurobi.
For formulation involving DAGsa and cyclic graphs, details can be found in (folder `MFD in DAGs`) and (folder `MFD with Cycles`), respectively.

## Installing Gurobi

Download the solver from [www.gurobi.com](www.gurobi.com), activate the (academic) license as instructed, and then install the Python API with:

```
pip3 install gurobipy
```
