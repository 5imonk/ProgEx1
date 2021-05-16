# -*- coding: utf-8 -*-
"""
Created on Sun May 16 00:47:09 2021

@author: Simon Konzett
"""
from docplex.mp.model import Model
from docplex.mp.model_reader import ModelReader
from docplex.util.environment import get_environment
import os


def read(filename):
    with open(filename) as f:
        n_nodes = [int(x) for x in next(f).split()][0]
        n_edges = [int(x) for x in next(f).split()][0]
        edges_inp = []
        for line in f:  # read rest of lines
            edges_inp.append([int(x) for x in line.split()])

    # introduce nodes, edges and arcs
    nodes = range(0, n_nodes)
    nodes_m0 = range(1, n_nodes)
    arcs = []
    k_arcs = []
    edges = []
    c = []

    for e in edges_inp:
        arcs.append((e[1], e[2]))
        for i in nodes_m0:
            k_arcs.append((e[1], e[2], i))

        edges.append((e[1], e[2]))
        c.append(e[3])
    for e in edges_inp:
        arcs.append((e[2], e[1]))
        for i in nodes_m0:
            k_arcs.append((e[2], e[1], i))

        c.append(e[3])

    return c, n_nodes, n_edges, nodes, arcs, k_arcs, edges


def build_common_vars(mdl, nodes, arcs):
    node_vars = mdl.binary_var_list(nodes, name='x')
    arc_vars = mdl.binary_var_list(arcs, name='a')
    # edge_vars = mdl.binary_var_list(edges, name='e')
    return mdl, node_vars, arc_vars


def build_common_constraints(mdl, nodes, node_vars, arcs, arc_vars, k, n_edges):
    mdl.add_constraint(mdl.sum(node_vars[i]
                               for i in nodes) == k+1, ctname='k_nodes')
    mdl.add_constraint(mdl.sum(arc_vars[l] for l, a in enumerate(
        arcs)) == k, ctname='k_arcs')  # actually k-1 arcs

    mdl.add_constraint(node_vars[0] == 1, ctname='root_active')
    mdl.add_constraint(mdl.sum(arc_vars[l] for l, a in enumerate(
        arcs) if a[1] == 0) == 0, ctname='root_in')
    mdl.add_constraint(mdl.sum(arc_vars[l] for l, a in enumerate(
        arcs) if a[0] == 0) == 1, ctname='root_out')

    for l, a in enumerate(arcs):
        if l < n_edges:
            ctname = 'one_direction_' + str(a)
            mdl.add_constraint(
                arc_vars[l] + arc_vars[l+n_edges] <= node_vars[a[0]], ctname=ctname)
        else:
            ctname = 'one_direction_' + str(a)
            mdl.add_constraint(
                arc_vars[l] + arc_vars[l-n_edges] <= node_vars[a[0]], ctname=ctname)

    for j in nodes:
        if j > 0:
            ctname = 'one_predecessor_' + str(j)
            mdl.add_constraint(mdl.sum(arc_vars[l] for l, a in enumerate(
                arcs) if a[1] == j) == node_vars[j], ctname=ctname)

    return mdl


def build_MTZ_vars(mdl, nodes, k):
    u_vars = mdl.integer_var_list(nodes, lb=0, ub=k, name='u')

    return mdl, u_vars


def build_MTZ_constraints(mdl, nodes, node_vars, arcs, arc_vars, u_vars, k):
    mdl.add_constraint(u_vars[0] == 0, ctname='root_first')

    for i in nodes:
        if i > 0:
            ctname = 'u_lower_' + str(i)
            mdl.add_constraint(node_vars[i] <= u_vars[i], ctname=ctname)

    for l, a in enumerate(arcs):
        ctname = 'sequential_' + str(a)
        i = a[0]
        j = a[1]
        mdl.add_constraint(u_vars[i] + arc_vars[l] <=
                           u_vars[j] + k * (1 - arc_vars[l]), ctname=ctname)

    return mdl


def build_SCF_vars(mdl, arcs, k):
    #f_vars = mdl.integer_var_list(arcs, lb=0, ub=k, name='f')
    f_vars = mdl.integer_var_list(arcs, lb=0, name='f')
    return mdl, f_vars


def build_SCF_constraints(mdl, nodes, node_vars, arcs, arc_vars, f_vars, k):

    mdl.add_constraint(mdl.sum(f_vars[l] for l, a in enumerate(
        arcs) if a[0] == 0) == k, ctname='root_out_flow')

    for i in nodes:
        if i > 0:
            ctname = 'flow_' + str(i)
            mdl.add_constraint(mdl.sum(f_vars[l] for l, a in enumerate(arcs) if a[1] == i)
                               - mdl.sum(f_vars[l]
                                         for l, a in enumerate(arcs) if a[0] == i)
                               == node_vars[i], ctname=ctname)

    for l, a in enumerate(arcs):
        ctname = 'f_upper_' + str(a)
        mdl.add_constraint(f_vars[l] <= k*arc_vars[l], ctname=ctname)

    return mdl


def build_MCF_vars(mdl, k_arcs):
    #fk_vars = mdl.integer_var_list(k_arcs, lb=0, ub=1, name='fk')    
    fk_vars = mdl.integer_var_list(k_arcs, lb=0, name='fk')
    return mdl, fk_vars


def build_MCF_constraints(mdl, nodes, node_vars, arcs, arc_vars, k_arcs, fk_vars):

    for i in nodes:
        if i > 0:
            ctname = 'root_out_flow' + str(i)
            mdl.add_constraint(mdl.sum(fk_vars[m] for m, ka in enumerate(k_arcs)
                                       if ka[0] == 0 and ka[2] == i) == node_vars[i], ctname=ctname)

    for i in nodes:
        if i > 0:
            ctname = 'last_flow' + str(i)
            mdl.add_constraint(mdl.sum(fk_vars[m] for m, ka in enumerate(k_arcs)
                                        if ka[1] == i and ka[2] == i) == node_vars[i], ctname=ctname)

    for i in nodes:
        for j in nodes:
            if i > 0 and j > 0 and i != j:
                ctname = 'flow' + str((i, j))
                mdl.add_constraint(mdl.sum(fk_vars[m] for m, ka in enumerate(k_arcs) if ka[1] == i and ka[2] == j)
                                   - mdl.sum(fk_vars[m] for m, ka in enumerate(k_arcs) if ka[0] == i and ka[2] == j) == 0, ctname=ctname)

    m = 0
    for l, a in enumerate(arcs):
        for i in nodes:
            if i > 0:
                ctname = 'fk_upper_' + str((a[0], a[1], i))
                mdl.add_constraint(fk_vars[m] <= arc_vars[l], ctname=ctname)
                m += 1

    return mdl


def build_MCF2_constraints(mdl, nodes, node_vars, arcs, arc_vars, k_arcs, fk_vars):

    for i in nodes:
        if i > 0:
            ctname = 'root_out_flow' + str(i)
            mdl.add_constraint(mdl.sum(fk_vars[m] for m, ka in enumerate(k_arcs)
                                       if ka[0] == 0 and ka[2] == i) == node_vars[i], ctname=ctname)

    for i in nodes:
        if i > 0:
            ctname = 'last_flow' + str(i)
            mdl.add_constraint(mdl.sum(fk_vars[m] for m, ka in enumerate(k_arcs) if ka[1] == i and ka[2] == i) - mdl.sum(
                fk_vars[m] for m, ka in enumerate(k_arcs) if ka[0] == i and ka[2] == i) == node_vars[i], ctname=ctname)

    for i in nodes:
        for j in nodes:
            if i > 0 and j > 0 and i != j:
                ctname = 'flow' + str((i, j))
                mdl.add_constraint(mdl.sum(fk_vars[m] for m, ka in enumerate(k_arcs) if ka[1] == i and ka[2] == j)
                                   - mdl.sum(fk_vars[m] for m, ka in enumerate(k_arcs) if ka[0] == i and ka[2] == j) == 0, ctname=ctname)

    m = 0
    for l, a in enumerate(arcs):
        for i in nodes:
            if i > 0:
                ctname = 'fk_upper_' + str((a[0], a[1], i))
                mdl.add_constraint(fk_vars[m] <= arc_vars[l], ctname=ctname)
                m += 1

    return mdl


def write_log(mdl, logfile, k, f, pre_time):

    obj_val = mdl.objective_value
    details = mdl.get_solve_details()

    with open(logfile, 'a') as out_file:

        out_file.write(mdl.name + ', graph: ' + f + ', k=' + str(k) + '\n')
        out_file.write('status ' + details.status + '\n')
        out_file.write('pre_time ' + str(pre_time) + ' s\n')
        out_file.write('time ' + str(details.time) + ' s\n')
        out_file.write('problem ' + details.problem_type + '\n')
        out_file.write('gap ' + str(details.mip_relative_gap) + ' %\n')
        out_file.write('objective value: ' + str(obj_val) + os.linesep)

        out_file.close()


def write_in_files(mdl, run_log, k, f, pre_time):

    write_log(mdl, run_log, k, f, pre_time)
    write_log(mdl, 'solution_log.txt', k, f, pre_time)
    
    cwd = os.getcwd()

    with get_environment().get_output_stream(cwd+'\solutions\solution_' + mdl.name + '_' + f[:3]
                                             + '_k_' + str(k) + '.json') as fp:
        mdl.solution.export(fp, 'json')


def update_table(table, mdl, k, f, pre_time):
    # table = [['formulation', 'graph', 'k', 'objective_value', 'time', 'gap', 'status', 'pre_time'
    #       'n_vars', 'n_bin_vars', 'n_int_vars', 'n_cont_vars', 'n_constr', 'n_lin_constr',
    #       'n_eq_constr', 'n_le_constr', 'n_ge_constr']
    #     ]

    details = mdl.get_solve_details()
    stats = mdl.statistics
    table.append([mdl.name, f[:3], k, mdl.objective_value, details.time, details.mip_relative_gap, details.status, pre_time,
                  stats.number_of_variables, stats.number_of_binary_variables,
                  stats.number_of_integer_variables, stats.number_of_continuous_variables,
                  stats.number_of_constraints, stats.number_of_linear_constraints,
                  stats.number_of_eq_constraints, stats.number_of_le_constraints, stats.number_of_ge_constraints])
    return table


def initialize():
    filenames = []
    for i in range(1, 11):
        if i < 10:
            filenames.append('g0'+str(i)+'.dat')
        else:
            filenames.append('g'+str(i)+'.dat')

    run_log = 'run_log.txt'
    if os.path.exists(run_log):
        os.remove(run_log)
    else:
        print('The file does not exist')

    table = [['formulation', 'graph', 'k', 'objective_value', 'time', 'gap', 'status', 'pre_time',
              'n_vars', 'n_bin_vars', 'n_int_vars', 'n_cont_vars', 'n_constr', 'n_lin_constr',
              'n_eq_constr', 'n_le_constr', 'n_ge_constr']
             ]
    return filenames, run_log, table


def buildModel(modelname, k, f):
    c, n_nodes, n_edges, nodes, arcs, k_arcs, edges = read(f)
    
    mdl = Model(modelname)
    
    mdl, node_vars, arc_vars = build_common_vars(mdl, nodes, arcs)
    mdl = build_common_constraints(
        mdl, nodes, node_vars, arcs, arc_vars, k, n_edges)

    if modelname == 'MTZ':
        mdl, u_vars = build_MTZ_vars(mdl, nodes, k)
        mdl = build_MTZ_constraints(
            mdl, nodes, node_vars, arcs, arc_vars, u_vars, k)
    elif modelname == 'SCF':
        mdl, f_vars = build_SCF_vars(mdl, arcs, k)
        mdl = build_SCF_constraints(
            mdl, nodes, node_vars, arcs, arc_vars, f_vars, k)
    elif modelname == 'MCF':
        mdl, fk_vars = build_MCF_vars(mdl, k_arcs)
        mdl = build_MCF_constraints(
            mdl, nodes, node_vars, arcs, arc_vars, k_arcs, fk_vars)
    elif modelname == 'MCF2':
        mdl, fk_vars = build_MCF_vars(mdl, k_arcs)
        mdl = build_MCF2_constraints(
            mdl, nodes, node_vars, arcs, arc_vars, k_arcs, fk_vars)

    mdl.minimize(mdl.sum(arc_vars[l]*c[l] for l, a in enumerate(arcs)))  
    
    cwd = os.getcwd()
    if k == 2 and f == 'g01.dat':
        mdl.export_as_lp(cwd)
        
    mdl.export_as_lp(cwd + '\lp\lp_'+mdl.name+'k'+str(k)+str(f[:3])+'.lp')
    
    return mdl

def importModelfromLP(modelname, k, f):
    cwd = os.getcwd()
    lp_path = cwd + '\lp\lp_' + modelname + 'k' + str(k) + str(f[:3]) + '.lp'
    if os.path.isfile(lp_path):
        rd = ModelReader()
        mdl = rd.read_model(lp_path, model_name=modelname)    
    else:
        c, n_nodes, n_edges, nodes, arcs, k_arcs, edges = read(f)
        mdl = buildModel(modelname, k, f, nodes, arcs, k_arcs, n_edges, c)
        
    return mdl
