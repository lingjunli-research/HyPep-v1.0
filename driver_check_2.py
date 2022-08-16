# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 16:41:36 2022

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
#discovery mode choice, boolean response, 1 = perform discovery mode, 0 = don't perform discovery mode
with open('discovery_mode.pkl', 'rb') as file_r:      
    # Call load method to deserialze
    to_discover = pickle.load(file_r)
       
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
    
#minimum alc
with open('min_alc.pkl', 'rb') as file_a:      
    # Call load method to deserialze
    malc_entry = pickle.load(file_a)

#minimum local confience
with open('discovery_alc.pkl', 'rb') as file_b:      
    # Call load method to deserialze
    disc_alc_entry_val_choice = pickle.load(file_b)

#motif file path
with open('motif_path.pkl', 'rb') as file_c:      
    # Call load method to deserialze
    motif_path = pickle.load(file_c)

#hypep FDR cutoff for discovery mode
with open('discovery_fdr_score.pkl', 'rb') as file_d:      
    # Call load method to deserialze
    hypep_score_cutoff_fdr = pickle.load(file_d)

if len(motif_path)>0: #motif decision
    motif_dec = 1
else:
    motif_dec = 0


""" SHS MODULE """

### target and decoy databases in defafaultdict form

raw_target_run = SHS_algorithm.algorithm.main_algorithm(float(window_size), database_sequence_data())
with open('raw_target_run.pkl', 'wb') as file_d:
    pickle.dump(raw_target_run, file_d)