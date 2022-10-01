# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 17:47:49 2022

@author: lawashburn
"""

import pickle
import pandas as pd

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
    raw_target_run_load = pickle.load(file_e)

raw_target_run = []

list_df = pd.DataFrame([sub.split(" ") for sub in raw_target_run_load])

list_df[3] = list_df[3].apply(pd.to_numeric,errors='coerce')
list_df_filtered = list_df

precursor_ID = list_df_filtered[0].values.tolist()
scan_ID = list_df_filtered[1].values.tolist()
db_ID = list_df_filtered[2].values.tolist()
score_ID = list_df_filtered[3].values.tolist()

top_range = len(list_df_filtered)

for a in range(0,top_range):
    
    nest = []
    
    precursor_ID_select = precursor_ID[a]
    scan_ID_select = scan_ID[a]
    db_ID_select = db_ID[a]
    score_ID_select = str(score_ID[a])
    
    nest.append(precursor_ID_select)
    nest.append(scan_ID_select)
    nest.append(db_ID_select)
    nest.append(score_ID_select)
    
    raw_target_run.append(nest)


with open('raw_decoy_run.pkl', 'rb') as file_f:      
    # Call load method to deserialze
    raw_decoy_run = pickle.load(file_f)
    
with open('discovery_FDR_filtered_score.pkl', 'rb') as file_g:      
    # Call load method to deserialze
    discovery_FDR_filtered_score = pickle.load(file_g)

    
with open('FDR_filtered_target_list.pkl', 'rb') as file_h:      
    # Call load method to deserialze
    FDR_filtered_target_list = pickle.load(file_h)

with open('discovery_FDR_filtered_target_list.pkl', 'rb') as file_i:      
    # Call load method to deserialze
    discovery_FDR_filtered_target_list = pickle.load(file_i)

if to_discover == 1:
    discovery_list = identification_mode.FDR_Calc.filter_by_score(raw_target_run, float(discovery_FDR_filtered_score))
    print(discovery_list)
    perfect_and_dup_precursor_removed_list = discovery_mode.precursor_sorted(discovery_list)
    ALC_filter_list = discovery_mode.ALC_filter(perfect_and_dup_precursor_removed_list, float(malc_entry))
    Local_Confidence_filter_list = discovery_mode.local_confidence_filter(ALC_filter_list,float(disc_alc_entry_val_choice))
    match_removed_list = identification_mode.ID_modifications.removing_matches(FDR_filtered_target_list,discovery_FDR_filtered_target_list)
    if motif_dec == 1:
        Motif_Filter_list = discovery_mode.motif_filter(FDR_filtered_target_list)
        
        with open('Motif_Filter_list.pkl', 'wb') as file_v:
            pickle.dump(Motif_Filter_list, file_v)