# -*- coding: utf-8 -*-
"""
Created on Sat May 15 21:34:28 2021

@author: Simon Konzett
"""

# from docplex.mp.model import Model
# from docplex.mp.model_reader import ModelReader
import pandas as pd
import time
# import os
# from kMST_Lib import read
# from kMST_Lib import build_common_constraints
# from kMST_Lib import build_common_vars
# from kMST_Lib import build_MTZ_constraints
# from kMST_Lib import build_MTZ_vars
# from kMST_Lib import build_SCF_constraints
# from kMST_Lib import build_SCF_vars
# from kMST_Lib import build_MCF_constraints
# from kMST_Lib import build_MCF_vars
# from kMST_Lib import build_MCF2_constraints
from kMST_Lib import write_in_files
from kMST_Lib import update_table
from kMST_Lib import initialize
from kMST_Lib import importModelfromLP
from kMST_Lib import buildModel

#modelnames = ['MTZ']
modelnames = ['MTZ', 'SCF', 'MCF', 'MCF2']
#modelnames = ['MCF', 'MCF2']

g_min = 1 #first graph
g_max = 4 #last graph

read_lp_model = False

filenames, run_log, table = initialize()

for i, f in enumerate(filenames):
    with open(f) as ff:
        n_nodes = [int(x) for x in next(ff).split()][0]   
                                        
    for modelname in modelnames:
        if i+1 <= g_max and i+1 >= g_min:
            # c, n_nodes, n_edges, nodes, arcs, k_arcs, edges = read(f)
            for sc in [0.2, 0.5]:                
                k = int(sc*(n_nodes-1))
                
                start_time = time.perf_counter()
                
                if read_lp_model: 
                    mdl = importModelfromLP(modelname, k, f)                   
                else:
                    mdl = buildModel(modelname, k, f)
                    
                    # mdl = Model(modelname)
                    # mdl.set_time_limit(720)
                    
                    # print('\n' + modelname + ', graph ' + f + ', k=' + str(k) + '\n')
        
                    # mdl, node_vars, arc_vars = build_common_vars(mdl, nodes, arcs)
                    # mdl = build_common_constraints(
                    #     mdl, nodes, node_vars, arcs, arc_vars, k, n_edges)
        
                    # if modelname == 'MTZ':
                    #     mdl, u_vars = build_MTZ_vars(mdl, nodes, k)
                    #     mdl = build_MTZ_constraints(
                    #         mdl, nodes, node_vars, arcs, arc_vars, u_vars, k)
                    # elif modelname == 'SCF':
                    #     mdl, f_vars = build_SCF_vars(mdl, arcs, k)
                    #     mdl = build_SCF_constraints(
                    #         mdl, nodes, node_vars, arcs, arc_vars, f_vars, k)
                    # elif modelname == 'MCF':
                    #     mdl, fk_vars = build_MCF_vars(mdl, k_arcs)
                    #     mdl = build_MCF_constraints(
                    #         mdl, nodes, node_vars, arcs, arc_vars, k_arcs, fk_vars)
                    # elif modelname == 'MCF2':
                    #     mdl, fk_vars = build_MCF_vars(mdl, k_arcs)
                    #     mdl = build_MCF2_constraints(
                    #         mdl, nodes, node_vars, arcs, arc_vars, k_arcs, fk_vars)
        
                    # mdl.minimize(mdl.sum(arc_vars[l]*c[l] for l, a in enumerate(arcs)))  
                
                mdl.set_time_limit(720)
                print('\n' + mdl.name + ', graph ' + f + ', k=' + str(k) + '\n') 
                
                end_time = time.perf_counter()
                
                pre_time = end_time - start_time                
                mdl.solve(log_output=True)  
                
                write_in_files(mdl, run_log, k, f, pre_time)
                table = update_table(table, mdl, k, f, pre_time)
                
                df = pd.DataFrame(table[1:], columns=table[0])
                df.to_csv('kMST.csv', sep =',', float_format="%.0f")
