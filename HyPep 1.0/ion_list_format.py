# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 10:43:41 2022

@author: lawashburn
"""

import csv
import pandas as pd
import os

from user_input import ion_list_directory_path
from user_input import target_list_out_results
from user_input import formatted_ion_list_path

ion_lists_path = ion_list_directory_path
target_sequence_list = target_list_out_results
final_directory = formatted_ion_list_path

target_list = pd.read_csv(target_sequence_list)
target_sequences = target_list['Target Sequence'].values.tolist()


def get_file_names_with_strings(str_list):
    full_list = os.listdir(ion_lists_path)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

ion_mod_full = pd.DataFrame()
for b in target_sequences:
    file_type2 = 'Theoretical_b_y_ion_'
    file_query2 = file_type2 + b + '.txt' #search query
    b_y_ion = (get_file_names_with_strings([file_query2]))
    b_y_ion_path = ion_lists_path + '\\' + b_y_ion[0] 
    prospector = pd.read_csv(b_y_ion_path,skiprows=[1]) #reads ion file
    prospector = prospector.dropna(axis=1, how='all')
    
    prospector_real = pd.DataFrame()
    prospector_real['mass'] = prospector['mass']
    prospector_real['ion'] = prospector['ion']
    prospector_real['formula'] = prospector['formula']
    prospector_real = prospector_real.dropna(how='all')
    prospector_real['ion name length'] = prospector_real['ion'].str.len()
    prospector_len_filter =  prospector_real[prospector_real['ion name length'] > 1] #removes all single AA entires
    prospector_b = prospector_real[prospector_real["ion"].str.contains('b')]
    prospector_y = prospector_real[prospector_real["ion"].str.contains('y')]
    
    ion_mod = pd.concat([prospector_b,prospector_y], axis=0,ignore_index=True)
    ion_mod = ion_mod.drop(columns=['ion name length'])
    
    out_path = final_directory + '\\Theoretical_b_y_ion_list_' + b + '.csv'
    with open(out_path,'w',newline='') as filec:
                            writerc = csv.writer(filec)
                            ion_mod.to_csv(filec,index=False)
    ion_mod_full = ion_mod_full.append(ion_mod)
