import math
import sys
import more_itertools

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

def read_input_standard(graphfile, weighttype):
    
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

def read_input_subpath(graphfile):

    lines = open(graphfile,'r').read().split('\n')

    subpaths_number = 0
    subpaths = {}
    subpaths_weights = {}
    
    i = 0
    j = 1
    while(True):
        if i >= len(lines):
            break
        line = lines[i].split(" ")
        if line[0] == '#':
            i = i + 1
            continue
        if len(line) == 1: 
            subpaths_number = line[0]
        else:
            subpaths[j] = [(a,b) for a,b in more_itertools.pairwise(line[1:len(line)])]
            subpaths_weights[j] = line[0]
            j = j + 1
        i = i + 1

    return subpaths,subpaths_weights

def read_input_inexact(graphfile, weighttype):
    
    lines = open(graphfile,'r').readlines()
    
    num_nodes = 0
    edges = dict()
    fup = dict()
    fdown = dict()
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
        elif len(elements) == 4: # True edge
            flow_value = float(elements[2])
            flow_up = float(elements[3])
            if flow_value < 0:
                sys.exit("ERROR: The input graph cannot have negative flow values")
            if weighttype.startswith('int'):
                if flow_value - int(flow_value) > 0:
                    sys.exit("ERROR: The input graph contains non-integer flow values, but the weights of the paths are required to be integers (--weighttype {weighttype})")
            
            fup[(elements[0], elements[1])] = flow_value
            fdown[(elements[0], elements[1])] = flow_up

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

    return {'vertices': list(in_neighbors.keys()), 'fup' : fup,'fdown': fdown, 'out_neighbors': out_neighbors, 'in_neighbors': in_neighbors, 'source': source, 'sink': sink, 'max_flow': max_flow}

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