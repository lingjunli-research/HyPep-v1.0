# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 17:16:07 2022

@author: lawashburn
"""

import pickle

### SHS Module Programs
from input_data import database_sequence_data
import decoy_database as decoy
import SHS_algorithm
import identification_mode
import discovery_mode
import pandas as pd

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
    score_ID_select = score_ID[a]
    
    nest.append(precursor_ID_select)
    nest.append(scan_ID_select)
    nest.append(db_ID_select)
    nest.append(score_ID_select)
    
    raw_target_run.append(nest)

with open('raw_decoy_run.pkl', 'rb') as file_e:      
    # Call load method to deserialze
    raw_decoy_run_load = pickle.load(file_e)

raw_decoy_run = []

list_df = pd.DataFrame([sub.split(" ") for sub in raw_decoy_run_load])

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
    score_ID_select = score_ID[a]
    
    nest.append(precursor_ID_select)
    nest.append(scan_ID_select)
    nest.append(db_ID_select)
    nest.append(score_ID_select)
    
    raw_decoy_run.append(nest)


""" SHS MODULE """

### target and decoy databases in defafaultdict form

closest_to_user_FDR, FDR_filtered_score, FDR_filtered_target_list = identification_mode.FDR_Calc.FDR_filter(float(fdr_percent),
                                                                                                            identification_mode.ID_modifications.precursor_sorted(raw_target_run),
                                                                                                            identification_mode.ID_modifications.precursor_sorted(raw_decoy_run))

with open('FDR_filtered_score.pkl', 'wb') as file_d:
    pickle.dump(FDR_filtered_score, file_d)

with open('closest_to_user_FDR.pkl', 'wb') as file_e:
    pickle.dump(closest_to_user_FDR, file_e)
    
with open('FDR_filtered_target_list.pkl', 'wb') as file_f:
    pickle.dump(FDR_filtered_target_list, file_f)
    
dummy_v = 'dummy variable'