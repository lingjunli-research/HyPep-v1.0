# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 17:45:04 2022

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


with open('raw_target_run.pkl', 'rb') as file_e:      
    # Call load method to deserialze
    raw_target_run = pickle.load(file_e)

with open('raw_decoy_run.pkl', 'rb') as file_f:      
    # Call load method to deserialze
    raw_decoy_run = pickle.load(file_f)
    
discovery_closest_to_user_FDR, discovery_FDR_filtered_score, discovery_FDR_filtered_target_list = identification_mode.FDR_Calc.FDR_filter(float(hypep_score_cutoff_fdr), 
                                                                                                                                          identification_mode.ID_modifications.precursor_sorted(raw_target_run),
                                                                                                                                          identification_mode.ID_modifications.precursor_sorted(raw_decoy_run))

with open('discovery_FDR_filtered_score.pkl', 'wb') as file_q:
    pickle.dump(discovery_FDR_filtered_score, file_q)
with open('discovery_closest_to_user_FDR.pkl', 'wb') as file_t:
    pickle.dump(discovery_closest_to_user_FDR, file_t)
with open('discovery_FDR_filtered_target_list.pkl', 'wb') as file_v:
    pickle.dump(discovery_FDR_filtered_target_list, file_v)
    
