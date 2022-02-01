import sys
import argparse
import gurobipy as gp
from gurobipy import GRB
import math
import os

## This lower bound assumes that there are no edges with flow zero.
## If there are, then it has to be changed to equal 
## the number of in-/out-neighbors with non-zero flow.
def lower_bound_degree(neighbors):
    
    max_deg = 1
    for node, node_neighbors in neighbors.items():
        if len(node_neighbors) > max_deg:
            max_deg = len(node_neighbors)
    
    return max_deg

def lower_bound_distinct_flow_values(edges):
    
    flow_values = set(edges.values())
    return math.ceil(math.log2(len(flow_values)))

def get_extremity(neighbors, extremity_type):
    
    extremity = None
    for node, node_neighbors in neighbors.items():
        if len(node_neighbors) == 0:
            if extremity != None:
                sys.exit(f"ERROR: input graph has more than one {extremity_type}")
            else:
                extremity = node

    if extremity == None:
        sys.exit(f"ERROR: The input graph has no {extremity_type} (and thus it is not a DAG)")

    return extremity

def read_input(graphfile, weighttype):
    
    lines = open(graphfile,'r').readlines()
    
    num_nodes = 0
    edges = dict()
    out_neighbors = dict()
    in_neighbors = dict()
    max_flow = 0

    for line in lines:
        line = line.strip()
        if line[0] == '#' or line == '':
            continue
        elements = line.split(' ')
        if len(elements) == 1: # Number of nodes
            num_nodes = int(elements[0])
        elif len(elements) == 3: # True edge
            flow_value = float(elements[2])
            if flow_value < 0:
                sys.exit("ERROR: The input graph cannot have negative flow values")
            if weighttype.startswith('int'):
                if flow_value - int(flow_value) > 0:
                    sys.exit("ERROR: The input graph contains non-integer flow values, but the weights of the paths are required to be integers (--weighttype {weighttype})")
            edges[(elements[0], elements[1])] = flow_value
            max_flow = max(max_flow,flow_value)
            
            if elements[0] not in out_neighbors:
                out_neighbors[elements[0]] = []
            out_neighbors[elements[0]].append(elements[1])

            if elements[1] not in in_neighbors:
                in_neighbors[elements[1]] = []
            in_neighbors[elements[1]].append(elements[0])

            if elements[0] not in in_neighbors:
                in_neighbors[elements[0]] = []
            if elements[1] not in out_neighbors:
                out_neighbors[elements[1]] = []
        else:
            sys.exit("ERROR: input file contains an ill-formatted line")
    if num_nodes != len(in_neighbors):
        sys.exit(f"ERROR: expecting {num_nodes} nodes, the input graph has {len(in_neighbors)} nodes")

    source = get_extremity(in_neighbors, 'source')
    sink = get_extremity(out_neighbors, 'sink')

    return {'vertices': list(in_neighbors.keys()), 'edges' : edges, 'out_neighbors': out_neighbors, 'in_neighbors': in_neighbors, 'source': source, 'sink': sink, 'max_flow': max_flow}

def write_outout(outputfilename, paths, w, weighttype):
    
    outputfile = open(outputfilename, 'w')

    if paths == None:
        outputfile.write("ERROR: Did not find any decomposition")
    else:
        for k in range(0, len(paths)):
            if weighttype.startswith('int'):
                outputfile.write(f"{int(w[k])} {paths[k]}\n")
            else:
                outputfile.write(f"{w[k]} {paths[k]}\n")

    outputfile.close()

def extract_paths(x, source, sink, out_neighbors, K):
    
    paths = []
    for k in range(0,K):
        vertex = source
        path = [vertex]
        while vertex != sink:
            for out_neighbor in out_neighbors[vertex]:
                if x[vertex,out_neighbor,k] == 1:
                    vertex = out_neighbor
                    break
            path.append(vertex)
        paths.append(path)

    return paths

def decompose_flow(vertices, edges, out_neighbors, in_neighbors, source, sink, max_flow, K, threads, weighttype):
        
    V = vertices
    E = set(edges.keys())
    data = dict()  
    W = max_flow
    
    try:
        #create extra sets
        T = [(i,j,k) for (i,j) in E for k in range(0,K)]
        SC = [k for k in range(0,K)]
        
        # Create a new model
        model = gp.Model("MFD")
        model.Params.LogToConsole = 0
        model.Params.Threads = threads

        print(f"INFO: Trying to decompose into {K} paths...")

        # Create variables
        x = model.addVars(T, vtype=GRB.BINARY, name="x")

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
requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)


args = parser.parse_args()

threads = args.threads
if threads == 0:
    threads = os.cpu_count()
print(f"INFO: Using {threads} threads for the Gurobi solver")

graph = read_input(args.input, args.weighttype)
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
    data = decompose_flow(vertices, edges, out_neighbors, in_neighbors, source, sink, max_flow, K, threads, args.weighttype)
    if data['message'] == "solved":
        w = data['weights']
        paths = data['paths']
        print(f"INFO: Found a decomposition into {K} paths")
        break

write_outout(args.output, paths, w, args.weighttype)