# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 12:37:45 2022

@author: lawashburn
"""

import os
import csv
import pandas as pd
import numpy as np
from datetime import datetime
now = datetime.now()

from user_input import top_FD_combined_out_file_name_path
from user_input import filtered_fragment_list_path
from user_input import target_list_out_results
from user_input import assign0_precursor_matches_out
from user_input import assign0_1_fragment_matches_out

spectra_import = top_FD_combined_out_file_name_path #path to directory containing spectra
fragment_list_import = filtered_fragment_list_path #path to directory with list of target fragment ions
target_list_import = target_list_out_results#path to directory with precursor list
working_directory = assign0_precursor_matches_out
final_dir = assign0_1_fragment_matches_out #path to directory for final, processed data


from user_input import sample_name
from user_input import rawconverter_formatted_out_file_name_path_txt
from user_input import precursor_error_cutoff
from user_input import fragment_error_cutoff
from user_input import precursor_charges
from user_input import fragment_charges

error_precursor = precursor_error_cutoff #+/- ppm, for precursor
error_fragment = fragment_error_cutoff #+/- Da, for fragment ion, charge state 1

h_mass = 1.00784

first_charge = 1

spectra_read= pd.read_csv(rawconverter_formatted_out_file_name_path_txt, sep=" ",skiprows=[0], names= ["m/z", "resolution", "charge", "intensity","MS2",'scan_number','precursor_charge','empty'])
spectra_value = pd.DataFrame()
spectra_value['Fragment m/z'] = spectra_read['m/z']
spectra_value['resolution'] = spectra_read['resolution']
spectra_value['charge'] = spectra_read['charge']
spectra_value['intensity'] = spectra_read['intensity']
spectra_value['Precursor'] = spectra_read['MS2']
spectra_value['Scan #'] = spectra_read['scan_number']
spectra_value['Precursor_Charge'] = spectra_read['precursor_charge']
spectra_value = spectra_value.drop(spectra_value[spectra_value['Precursor_Charge']==0].index) #remove charges equal to 0
spectra_value = spectra_value.drop(spectra_value[spectra_value['charge']!=0].index) #remove charges equal to 0
spectra_value = spectra_value.drop(spectra_value[spectra_value['intensity']<100].index) #remove intensity values less than 100

target_list = pd.read_csv(target_list_import)
sequence_targets = target_list['Target Sequence'].values.tolist()

spectra_value_first = spectra_value
spectra_value_first['charge'] = first_charge


precursor_matches = pd.DataFrame()

for a in sequence_targets:
    target_list_filtered = target_list[target_list['Target Sequence'] == a]
    
    for b in range(1,precursor_charges+1):
        precursor_target = target_list_filtered[str(b)].values.tolist()
        spectra_value_filtered = spectra_value[spectra_value['Precursor_Charge'] == b]
        
        for c in precursor_target:
            spectra_value_filtered['ppm error'] =  ((abs(spectra_value_filtered['Precursor'] - c))/c) * 1E6
            spectra_value_filtered = spectra_value_filtered[spectra_value_filtered['ppm error'] <= error_precursor]
            spectra_value_filtered['Possible sequence'] = a
            spectra_value_filtered['Theoretical precursor'] = c
            spectra_value_filtered['Theoretical precursor charge'] = b
            
            if len(spectra_value_filtered) > 0:
                precursor_matches = precursor_matches.append(spectra_value_filtered)
            else:
                pass
from user_input import assign0_precursor_matches_out
file_out_path = assign0_precursor_matches_out
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        precursor_matches.to_csv(filec,index=False)

updated_sequence_targets = precursor_matches['Possible sequence'].values.tolist()

sequence_no_dups = []
for o in updated_sequence_targets:
    if o not in sequence_no_dups:
        sequence_no_dups.append(o)

def get_file_names_with_strings(str_list):
    full_list = os.listdir(fragment_list_import)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

fragment_matches = pd.DataFrame()

for p in sequence_no_dups:
    d = first_charge
    file_query = 'Theoretical_b_y_fragment_list_' + p + '.csv'
    fragment_list = (get_file_names_with_strings([file_query]))
    for q in fragment_list:
        fragment_list_path = fragment_list_import + '\\' + q
        fragment_database = pd.read_csv(fragment_list_path)
        fragment_database_format = fragment_database.rename({'Mass (1)': '1','Mass (2)': '2','Mass (3)': '3','Mass (4)': '4' }, axis=1)
        seq_exp = precursor_matches[precursor_matches['Possible sequence'] == p]  
        theo_fragment = fragment_database_format[str(d)].values.tolist()
        seq_exp_filtered = seq_exp[seq_exp['charge'] == d]
        for u in theo_fragment:
                theo_frag_M = (u * d) - (h_mass * d)
                seq_exp_filtered['fragment M'] = (seq_exp_filtered['Fragment m/z'] * d) - (h_mass * d)
                seq_exp_filtered['theoretical fragment'] = u
                seq_exp_filtered['theoretical fragment M'] = theo_frag_M
                seq_exp_filtered['Da error'] = abs(seq_exp_filtered['fragment M'] - theo_frag_M)
                seq_exp_filter = seq_exp_filtered.sort_values(by='Da error')
                seq_exp_filter = seq_exp_filter[seq_exp_filter['Da error'] <= error_fragment]
                if len(seq_exp_filter) > 0:
                    fragment_matches = fragment_matches.append(seq_exp_filter)    
from  user_input import assign0_1_fragment_matches_out
file_out_path = assign0_1_fragment_matches_out
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        fragment_matches.to_csv(filec,index=False)
print('exporting fragment matches')
