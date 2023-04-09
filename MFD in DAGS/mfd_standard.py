#!/usr/bin/env python
# coding: utf-8

import os
import sys
import argparse
import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from collections import deque
from bisect import bisect
from copy import deepcopy

def get_edge(raw_edge):

    parts = raw_edge.split()
    return int(parts[0]), int(parts[1]), float(parts[2])

def get_graph(raw_graph):

    graph = {
        'n': 0,
        'edges': list()
    }

    try:
        lines = raw_graph.split('\n')[1:]
        if not lines[-1]:
            lines = lines[:-1]
        graph['n'], graph['edges'] = int(lines[0]), [get_edge(raw_e) for raw_e in lines[1:]]

    finally:
        return graph


def read_input_graphs(graph_file):

    graphs_raw = open(graph_file, 'r').read().split('#')[1:]
    return [get_graph(raw_g) for raw_g in graphs_raw]

def read_input(graph_file):

    return read_input_graphs(graph_file)


def mfd_algorithm(data):

    data['message'] = 'unsolved'
    for i in range(2, len(data['graph'].edges) + 1):
        if fd_fixed_size(data, i)['message'] == 'solved':
            return data

    return data


def build_base_ilp_model(data, size):

    graph = data['graph']
    max_flow_value = data['max_flow_value']
    sources = data['sources']
    sinks = data['sinks']
    B = 2

    # create extra sets
    T = [(u, v, i, k) for (u, v, i) in graph.edges(keys=True) for k in range(size)]
    SC = list(range(size))

    # Create a new model
    model = gp.Model('MFD')
    model.setParam('LogToConsole', 0)
    model.setParam('Threads', threads)


    # Create variables
    x = model.addVars(T, vtype=GRB.BINARY, name='x')
    w = model.addVars(SC, vtype=GRB.INTEGER, name='w', lb=1)
    z = model.addVars(T, vtype=GRB.CONTINUOUS, name='z', lb=0)

    # flow conservation
    for k in range(size):
        for v in graph.nodes:
            if v in sources:
                model.addConstr(sum(x[v, w, i, k] for _, w, i in graph.out_edges(v, keys=True)) == 1)
            if v in sinks:
                model.addConstr(sum(x[u, v, i, k] for u, _, i in graph.in_edges(v, keys=True)) == 1)
            if v not in sources and v not in sinks:
                model.addConstr(sum(x[v, w, i, k] for _, w, i in graph.out_edges(v, keys=True)) - sum(x[u, v, i, k] for u, _, i in graph.in_edges(v, keys=True)) == 0)

    # flow superposition
    for (u, v, i, f) in graph.edges(keys=True, data='flow'):
        model.addConstr(f - sum(z[u, v, i, k] for k in range(size)) == 0)

    # linearization
    for (u, v, i) in graph.edges(keys=True):
        for k in range(size):
            model.addConstr(z[u, v, i, k] <= max_flow_value * x[u, v, i, k])
            model.addConstr(w[k] - (1 - x[u, v, i, k]) * max_flow_value <= z[u, v, i, k])
            model.addConstr(z[u, v, i, k] <= w[k])

    return model, x, w, z


def get_solution(model, data, size):

    data['weights'], data['solution'] = list(), list()

    if model.status == GRB.OPTIMAL:
        graph = data['graph']
        T = [(u, v, i, k) for (u, v, i) in graph.edges(keys=True) for k in range(size)]

        w_sol = [0] * len(range(size))
        paths = [list() for _ in range(size)]
        for k in range(size):
            w_sol[k] = round(model.getVarByName(f'w[{k}]').x)
        for (u, v, i, k) in T:
            if round(model.getVarByName(f'x[{u},{v},{i},{k}]').x) == 1:
                paths[k].append((u, v, i))
        for k in range(len(paths)):
            paths[k] = sorted(paths[k])

        data['weights'], data['solution'] = w_sol, paths

    return data


def update_status(data, model):

    if model.status == GRB.OPTIMAL:
        data['message'] = 'solved'
        data['runtime'] = model.Runtime

    if model.status == GRB.INFEASIBLE:
        data['message'] = 'unsolved'
        data['runtime'] = 0

    return data


def fd_fixed_size(data, size):

    # calculate a flow decomposition into size paths
    try:
        # Create a new model
        model, _, _, _ = build_base_ilp_model(data, size)

        # objective function
        model.optimize()

        data = update_status(data, model)
        data = get_solution(model, data, size)

    except gp.GurobiError as e:
        print(f'Error code {e.errno}: {e}', file=sys.stderr)

    except AttributeError:
        print('Encountered an attribute error', file=sys.stderr)

    return data

def output_paths(output,paths,weights):
    
    numberOfPaths = len(paths)

    for nP in range(0,numberOfPaths):
        nodes = set()
        for (i,j,k) in paths[nP]:
            nodes.add(i)
            nodes.add(j)
        
        output.write(str(weights[nP]))
        for i in nodes:
            output.write(' '.join(['',str(i)]))
        output.write('\n')

def compute_graph_metadata(graph):

    # creation of NetworkX Graph
    ngraph = nx.MultiDiGraph()
    ngraph.add_weighted_edges_from(graph['edges'], weight='flow')

    # calculating source, sinks
    sources = [x for x in ngraph.nodes if ngraph.in_degree(x) == 0]
    sinks = [x for x in ngraph.nodes if ngraph.out_degree(x) == 0]

    # definition of data
    return {
        'graph': ngraph,
        'sources': sources,
        'sinks': sinks,
        'max_flow_value': max(ngraph.edges(data='flow'), key=lambda e: e[-1])[-1] if len(ngraph.edges) > 0 else -1,
    }

def solve_instances(graphs,output_file, output_stats=False):

    output = open(output_file, 'w+')
    if output_stats:
        stats = open(f'{output_file}.stats', 'w+')

    for g, graph in enumerate(graphs):

        output.write(f'# graph {g}\n')
        if output_stats:
            stats.write(f'# graph {g}\n')

        if not graph['edges']:
            continue

        mfd = compute_graph_metadata(graph)

        if output_stats:
            is_unique_decomposition = True

        if len(mfd['graph'].edges) > 0:

            mfd = mfd_algorithm(mfd)
            paths,weights = mfd['solution'],mfd['weights']
            output_paths(output,paths,weights)


    output.close()
    if output_stats:
        stats.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='''
        Computes paths for Inexact Minimum Flow Decomposition.
        This script uses the Gurobi ILP solver.
        ''',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-t', '--threads', type=int, default=0,
                        help='Number of threads to use for the Gurobi solver; use 0 for all threads (default 0).')
    
 
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i', '--input', type=str, help='Input filename', required=True)
    requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)

    args = parser.parse_args()

    threads = args.threads
    if threads == 0:
        threads = os.cpu_count()
    print(f'INFO: Using {threads} threads for the Gurobi solver')

    solve_instances(read_input(args.input),args.output)
    print("Done")
