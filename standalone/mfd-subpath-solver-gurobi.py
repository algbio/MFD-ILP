import sys
import argparse
import gurobipy as gp
from gurobipy import GRB
import math
import os
from st_fd import *

## This lower bound assumes that there are no edges with flow zero.
## If there are, then it has to be changed to equal 
## the number of in-/out-neighbors with non-zero flow.

def decompose_flow(vertices, edges, out_neighbors, in_neighbors, source, sink, max_flow, K, threads, weighttype,subpaths,subpaths_weights):
        
    V = vertices
    E = set(edges.keys())
    data = dict()  
    W = max_flow
    
    try:
        #create extra sets
        T = [(i,j,k) for (i,j) in E for k in range(0,K)]
        SC = [k for k in range(0,K)]
        R = [(k,s) for k in range(0,K) for s in range(1,len(subpaths))]
        
        # Create a new model
        model = gp.Model("MFD")
        model.Params.LogToConsole = 0
        model.Params.Threads = threads

        print(f"INFO: Trying to decompose into {K} paths...")

        # Create variables
        x = model.addVars(T, vtype=GRB.BINARY, name="x")
        r = model.addVars(R,vtype=GRB.BINARY,name="r")

        if weighttype == 'int+':
            w = model.addVars(SC,vtype=GRB.INTEGER, name="w",lb=0)
            z = model.addVars(T, vtype=GRB.CONTINUOUS, name="z", lb=0)
        elif weighttype == 'float+':
            w = model.addVars(SC,vtype=GRB.CONTINUOUS, name="w",lb=0)
            z = model.addVars(T, vtype=GRB.CONTINUOUS, name="z", lb=0)
        else:
            sys.exit(f"ERROR: Unknown weight type {weighttype}")
    
        model.setObjective(GRB.MINIMIZE)
       
        # flow conservation
        for k in range(0,K):
            for i in V:
                if i == source:
                    model.addConstr(sum(x[i,j,k] for j in out_neighbors[i]) == 1)
                elif i == sink:
                    model.addConstr(sum(x[j,i,k] for j in in_neighbors[i]) == 1)
                else:
                    model.addConstr(sum(x[i,j,k] for j in out_neighbors[i]) - sum(x[j,i,k] for j in in_neighbors[i]) == 0)

        # flow explanation
        for (i,j) in E:
            model.addConstr(edges[(i,j)] == sum(z[i,j,k] for k in range(0,K)))

        # linearization
        for (i,j) in E:
            for k in range(0,K):
                model.addConstr(z[i,j,k] <= W*x[i,j,k])
                model.addConstr(w[k] - (1 - x[i,j,k])*W <= z[i,j,k])
                model.addConstr(z[i,j,k] <= w[k])
            
        # subpath contraints
        for k in range(0,K):
            for s in range(1,len(subpaths)):
                model.addConstr(sum(x[i,j,k] for (i,j) in subpaths[s]) >= len(subpaths[s])*r[k,s])
        
        model.addConstrs(sum(r[k,j] for k in range(0,K)) >= 1 for j in range(1,len(subpaths)))
        
        # objective function
        model.optimize()
        w_sol = [0]*len(range(0,K))
        x_sol = {}
    
        
        if model.status == GRB.OPTIMAL:
            data['message'] = 'solved'
            data['runtime'] = model.Runtime
            
            for v in model.getVars():
                if 'w' in v.VarName:
                    elements = v.VarName.replace('[',',').replace(']',',').split(',')
                    k = int(elements[1])
                    w_sol[k] = v.x                            
                
                if 'x' in v.VarName:          
                    elements = v.VarName.replace('[',',').replace(']',',').split(',')
                    i = elements[1]
                    j = elements[2]
                    k = int(elements[3])
                    x_sol[i,j,k] = v.x
                
            paths = extract_paths(x_sol, source, sink, out_neighbors, K)
            
            #print(x_sol,w_sol,paths)
            data['weights'] = w_sol
            data['paths'] = paths
        
        if model.status == GRB.INFEASIBLE:
            data['message'] = 'unsolved'
            
        
    except gp.GurobiError as e:
        sys.exit('Error code ' + str(e.errno) + ': ' + str(e))

    except AttributeError:
        sys.exit('Encountered an attribute error')
    
    return data

parser = argparse.ArgumentParser(
    description="""
    Decompose a network flow into a minimum number of weighted paths. 
    This script uses the Gurobi ILP solver.
    """,
    formatter_class=argparse.RawTextHelpFormatter
    )
parser.add_argument('-wt', '--weighttype', type=str, default='int+', help='Type of path weights (default int+):\n   int+ (positive non-zero ints), \n   float+ (positive non-zero floats).')
parser.add_argument('-t', '--threads', type=int, default=0, help='Number of threads to use for the Gurobi solver; use 0 for all threads (default 0).')

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i', '--input', type=str, help='Input filename', required=True)
requiredNamed.add_argument('-s', '--subpaths',type=str,help='Subpaths file', required=True)
requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)


args = parser.parse_args()

threads = args.threads
if threads == 0:
    threads = os.cpu_count()
print(f"INFO: Using {threads} threads for the Gurobi solver")

graph = read_input_standard(args.input, args.weighttype)
subpaths,subpaths_weights = read_input_subpath(args.subpaths)
edges = graph['edges']
vertices = graph['vertices']
in_neighbors = graph['in_neighbors']
out_neighbors = graph['out_neighbors']
source = graph['source']
sink = graph['sink']
max_flow = graph['max_flow']

lower_bound = lower_bound_degree(out_neighbors)
lower_bound = max(lower_bound, lower_bound_degree(in_neighbors))
lower_bound = max(lower_bound, lower_bound_distinct_flow_values(edges))

minK = lower_bound
maxK = len(edges)

w = None
paths = None
for K in range(minK, maxK + 1):
    data = decompose_flow(vertices, edges, out_neighbors, in_neighbors, source, sink, max_flow, K, threads, args.weighttype,subpaths,subpaths_weights)
    if data['message'] == "solved":
        w = data['weights']
        paths = data['paths']
        print(f"INFO: Found a decomposition into {K} paths")
        break

write_outout(args.output, paths, w, args.weighttype)