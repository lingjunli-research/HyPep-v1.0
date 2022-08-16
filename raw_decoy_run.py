# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 17:05:28 2022

@author: lawashburn
"""

import pickle

### SHS Module Programs
from input_data import database_sequence_data
import decoy_database as decoy
import SHS_algorithm
import identification_mode
import discovery_mode

### PARAMETER

       
#FDR algorithm choice, multiple choices, key: 1=Reverse, 2=Shuffle, 3=Random, 4=Hybrid
with open('fdr_algorithm.pkl', 'rb') as file_s:      
    # Call load method to deserialze
    fdr_selection = pickle.load(file_s)
    
#FDR %
with open('fdr.pkl', 'rb') as file_f:      
    # Call load method to deserialze
    fdr_percent = pickle.load(file_f)
    
#Selected sliding window size
with open('window_size.pkl', 'rb') as file_q:      
    # Call load method to deserialze
    window_size = pickle.load(file_q)
    



""" SHS MODULE """

### target and decoy databases in defafaultdict form

raw_decoy_run = SHS_algorithm.algorithm.main_algorithm(float(window_size), decoy.decoy_database.modified_decoy_database(int(fdr_selection)))

with open('raw_decoy_run.pkl', 'wb') as file_d:
    pickle.dump(raw_decoy_run, file_d)