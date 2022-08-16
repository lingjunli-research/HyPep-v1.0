# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 13:36:02 2022

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
from user_input import assign0_working_directory_path
from user_input import assign0_final_directory_path
from user_input import sample_name
from user_input import assign0_precursor_matches_out
from user_input import rawconverter_formatted_out_file_name_path_txt
spectra_import = rawconverter_formatted_out_file_name_path_txt #path to directory containing spectra
fragment_list_import = filtered_fragment_list_path #path to directory with list of target fragment ions
target_list_import = target_list_out_results#path to directory with inclusion lists
working_directory = assign0_working_directory_path
final_dir = assign0_final_directory_path #path to directory for final, processed data
from user_input import precursor_error_cutoff
from user_input import fragment_error_cutoff
from user_input import precursor_charges
from user_input import fragment_charges

precursor_charges = int(precursor_charges)
fragment_charges = int(fragment_charges)

precursor_charges = range(1,precursor_charges+1)
error_precursor = precursor_error_cutoff #+/- ppm, for precursor
error_fragment = fragment_error_cutoff #+/- Da, for fragment ion, charge state 1
h_mass = 1.00784

fragment_charges = range(1, fragment_charges+1)
first_charge = fragment_charges[0]
remaining_charges = fragment_charges[1:]

target_list = pd.read_csv(target_list_import)
sequence_targets = target_list['Target Sequence'].values.tolist()

spectra_read= pd.read_csv(spectra_import, sep=" ",skiprows=[0], names= ["m/z", "resolution", "charge", "intensity","MS2",'scan_number','precursor_charge','empty'])
spectra_value = pd.DataFrame()
spectra_value['Fragment m/z'] = spectra_read['m/z']
spectra_value['resolution'] = spectra_read['resolution']
spectra_value['charge'] = spectra_read['charge']
spectra_value['intensity'] = spectra_read['intensity']
spectra_value['Precursor'] = spectra_read['MS2']
spectra_value['Scan #'] = spectra_read['scan_number']
spectra_value['Precursor_Charge'] = spectra_read['precursor_charge']
spectra_value = spectra_value.drop(spectra_value[spectra_value['Precursor_Charge']==0].index) #remove charges equal to 0
spectra_value_zeroes = spectra_value.drop(spectra_value[spectra_value['charge']!=0].index) #remove charges not equal to 0
spectra_value_zeroes = spectra_value_zeroes.drop(spectra_value_zeroes[spectra_value_zeroes['intensity']<100].index) #remove intensity values less than 100


def get_file_names_with_strings(str_list):
    full_list = os.listdir(working_directory)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

#input info here for loop 1, charge1

all_fragment_matches = pd.DataFrame()
unavailable_spectra = pd.DataFrame()
all_used_fragments = pd.DataFrame()

for zz in fragment_charges:
    zs = str(zz)
    file_query = '_zero_reassign_fragment_matches.csv'
    fragment_list = (get_file_names_with_strings([file_query]))
    for a in fragment_list:
        fragment_list_path = working_directory + '\\' + a
        fragment_read = pd.read_csv(fragment_list_path,dtype=object)
        all_used_fragments = all_used_fragments.append(fragment_read)

    spectra_value_zeroes['Fragment m/z']=spectra_value_zeroes['Fragment m/z'].astype(float)
    spectra_value_zeroes['Scan #']=spectra_value_zeroes['Scan #'].astype(float)
    all_used_fragments['Fragment m/z']=all_used_fragments['Fragment m/z'].astype(float)
    all_used_fragments['Scan #']=all_used_fragments['Scan #'].astype(float)
    
    available_spectra = pd.merge(spectra_value_zeroes, all_used_fragments,how='left', left_on=['Fragment m/z', 'Scan #'],right_on=['Fragment m/z','Scan #'])
    available_spectra_filtered = available_spectra[available_spectra['Possible sequence'].isna()] 

    if len(available_spectra_filtered)>0:                                                                                                                                             
        #unavailable_spectra = unavailable_spectra.append(available_spectra_filtered)
        
        #true_available_spectra =  pd.merge(spectra_value_zeroes, unavailable_spectra,how='left', left_on=['Fragment m/z','Scan #'],right_on=['Fragment m/z','Scan #'])
        #true_available_spectra = true_available_spectra[true_available_spectra['Possible sequence'].isna()]
        
        true_working_spectra = available_spectra_filtered
        true_working_spectra = true_working_spectra.iloc[:,:7]
        true_working_spectra['charge'] = zz
        true_working_spectra = true_working_spectra.rename(columns={"resolution_x": "resolution", "intensity_x": "intensity", "Precursor_x": "Precursor", "Precursor_Charge_x": "Precursor_Charge"})
        true_working_spectra = true_working_spectra.drop(['charge_x'], axis=1)
        true_working_spectra2 = true_working_spectra.drop_duplicates()
    
        
        precursor_matches = pd.DataFrame()
        
        for a in sequence_targets:
            target_list_filtered = target_list[target_list['Target Sequence'] == a]
            
            for b in precursor_charges:
                precursor_target = target_list_filtered[str(b)].values.tolist() #change to not a string
                spectra_value_filtered = true_working_spectra2[true_working_spectra2['Precursor_Charge'] == b]
                
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
        file_name = sample_name + '_' + zs + '_zero_reassign_precursor_matches.csv'
        file_out_path = working_directory + '\\' + file_name
        with open(file_out_path,'w',newline='') as filec:
                writerc = csv.writer(filec)
                precursor_matches.to_csv(filec,index=False)
        
        #precursor_matches = pd.read_csv(assign0_precursor_matches_out)
        updated_sequence_targets = precursor_matches['Possible sequence'].values.tolist()
       
        sequence_no_dups = []
        for o in updated_sequence_targets:
            if o not in sequence_no_dups:
                sequence_no_dups.append(o)

        if len(sequence_no_dups)>0:        
            def get_fragment_file_names_with_strings(str_list):
                full_list = os.listdir(fragment_list_import)
                final_list = [nm for ps in str_list for nm in full_list if ps in nm]
                return final_list
            
            fragment_matches = pd.DataFrame()        
            #for d in remaining_charges:
            for d in fragment_charges:
                for p in sequence_no_dups:
                    file_query = 'Theoretical_b_y_fragment_list_' + p + '.csv'

                    fragment_list = (get_fragment_file_names_with_strings([file_query]))
                                       
                    for q in fragment_list:
                        fragment_list_path = fragment_list_import + '\\' + q
                        fragment_database = pd.read_csv(fragment_list_path)                        
                        fragment_database_format = fragment_database.rename({'1': '1','2': '2','3': '3','4': '4' }, axis=1)
                        seq_exp = precursor_matches[precursor_matches['Possible sequence'] == p]  
                        theo_fragment = fragment_database_format[str(d)].values.tolist() #change to not a string
                        seq_exp_filtered = seq_exp[seq_exp['charge'] == d]
                        for u in theo_fragment:
                                theo_frag_M = (u * d) - (h_mass * d)
                                seq_exp_filtered['fragment M'] = (seq_exp_filtered['Fragment m/z'] * d) - (h_mass * d)
                                seq_exp_filtered['theoretical fragment'] = u
                                seq_exp_filtered['theoretical fragment M'] = theo_frag_M
                                seq_exp_filtered['Da error'] = abs(seq_exp_filtered['fragment M'] - theo_frag_M)
                                #If doesn't work, clear the error restriction
                                seq_exp_filter = seq_exp_filtered.sort_values(by='Da error')
                                seq_exp_filter = seq_exp_filter[seq_exp_filter['Da error'] <= error_fragment]
                                if len(seq_exp_filter) > 0:
                                    fragment_matches = fragment_matches.append(seq_exp_filter)
                                else:
                                    pass                                
            else:
                pass
            
            file_name = sample_name + '_' + zs + '_zero_reassign_fragment_matches.csv'
            file_out_path = working_directory + '\\' + file_name
            with open(file_out_path,'w',newline='') as filec:
                    writerc = csv.writer(filec)
                    fragment_matches.to_csv(filec,index=False)
    
