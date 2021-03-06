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
import warnings


#MTZ Miller Tucker Zemlin
#SCF Single Commodity Flow
#MCF weaker Multi Commmodity Formulation 
#MCF2 stronger Multi Commmodity Formulation 

modelnames = ['MTZ','SCF','MCF2']
#modelnames = ['MCF', 'MCF2']
#modelnames = ['MCF2']
#modelnames = ['MTZ', 'SCF']

# graphs from to that are considered
g_min = 1 #first graph
g_max = 5 #last graph

# set k as k = k_range |V|
k_range = [0.2, 0.5]

# if the model constraints, variables are changed new lp models have to be written
# then turn on True
# if no lp files are there formulation is made from scratch anyway
read_lp_model = True
if read_lp_model:
    warnings.warn("If read_lp_model Flag ist set TRUE changes made for variables and constraints are not yet reflected in the model as it is read from an LP file")    

filenames, run_log, table = initialize()

for i, f in enumerate(filenames):
    with open(f) as ff:
        n_nodes = [int(x) for x in next(ff).split()][0]   
                                        
    for modelname in modelnames:
        if i+1 <= g_max and i+1 >= g_min: 
            for sc in k_range:                
                k = int(sc*(n_nodes-1))
                
                start_time = time.perf_counter()
                
                if read_lp_model:                     
                    mdl = importModelfromLP(modelname, k, f)
                else:
                    mdl = buildModel(modelname, k, f)
                
                mdl.set_time_limit(900)
                print('\n' + mdl.name + ', graph ' + f + ', k=' + str(k) + '\n') 
                
                end_time = time.perf_counter()
                
                pre_time = end_time - start_time                
                mdl.solve(log_output=True)  
                
                write_in_files(mdl, run_log, k, f, pre_time)
                table = update_table(table, mdl, k, f, pre_time)
                
                df = pd.DataFrame(table[1:], columns=table[0])
                df.to_csv('kMST.csv', sep =',', float_format="%.4f")
