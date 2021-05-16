# -*- coding: utf-8 -*-
"""
Created on Sat May 15 21:34:28 2021

@author: Simon Konzett
"""

import pandas as pd
import time
from kMST_Lib import write_in_files
from kMST_Lib import update_table
from kMST_Lib import initialize
from kMST_Lib import importModelfromLP
from kMST_Lib import buildModel

#modelnames = ['MTZ']
modelnames = ['MTZ', 'SCF', 'MCF', 'MCF2']
#modelnames = ['MCF', 'MCF2']

# graphs from to that are considered
g_min = 1 #first graph
g_max = 4 #last graph

# if the model constraints, variables are changed new lp models have to be written
# then turn on True
read_lp_model = True


filenames, run_log, table = initialize()

for i, f in enumerate(filenames):
    with open(f) as ff:
        n_nodes = [int(x) for x in next(ff).split()][0]   
                                        
    for modelname in modelnames:
        if i+1 <= g_max and i+1 >= g_min:            
            for sc in [0.2, 0.5]:                
                k = int(sc*(n_nodes-1))
                
                start_time = time.perf_counter()
                
                if read_lp_model: 
                    mdl = importModelfromLP(modelname, k, f)                   
                else:
                    mdl = buildModel(modelname, k, f)
                
                mdl.set_time_limit(720)
                print('\n' + mdl.name + ', graph ' + f + ', k=' + str(k) + '\n') 
                
                end_time = time.perf_counter()
                
                pre_time = end_time - start_time                
                mdl.solve(log_output=True)  
                
                write_in_files(mdl, run_log, k, f, pre_time)
                table = update_table(table, mdl, k, f, pre_time)
                
                df = pd.DataFrame(table[1:], columns=table[0])
                df.to_csv('kMST.csv', sep =',', float_format="%.0f")
