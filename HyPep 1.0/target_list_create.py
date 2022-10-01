# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 14:37:30 2022

@author: lawashburn
"""

import pandas as pd
import csv
import os
from user_input import ion_list_directory_path
from user_input import AMM_v_SHS_out_results
from user_input import sample_name
from datetime import datetime
from user_input import precursor_charges
from user_input import missing_fragment_files
now = datetime.now()

print('start time:', datetime.now())

ion_directory = ion_list_directory_path
working_directory = ion_list_directory_path
unique_amm_sequences = AMM_v_SHS_out_results


H_mass = 1.00784
charges = precursor_charges

def get_file_names_with_strings_list(str_list): #definition for finding a file containing a string in filename in specified directory
    full_list = os.listdir(ion_directory)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

amm_sequences = pd.read_csv(unique_amm_sequences)
target_sequences = amm_sequences['PTM'].values.tolist() #list of sequences

MH_list = []
ion_file_name = []
target_sequence = []

missing_seqs = []

for a in target_sequences: #for each sequence in the sequence list
    file_type = 'Theoretical_b_y_ion_'
    ion_query = file_type + a + '.txt' #search for ion list pertaining to the sequence
    ion_file = (get_file_names_with_strings_list([ion_query])) #search for the file based on query
    if len(ion_file) == 1: #ensures only one file matched
        for b in ion_file:
            ion_path = ion_directory + '\\' + b
            ion_read = pd.read_csv(ion_path, skiprows=[1]) #reads identified .csv
            ion_read_filter = ion_read.drop(ion_read[ion_read['ion']!= 'MH'].index) #keeps only the MH value from the ProteinProspector ion table
            MH = ion_read_filter['mass'].values.tolist() #selects the MH value
            for c in MH:
                MH_list.append(c)
                ion_file_name.append(b)
                target_sequence.append(a)
    if len(ion_file)>1: #throws an error if more than one ion file exists
        raise ValueError('More than one ion list located for:',a)
    if len(ion_file)==0: #throws an error if ion file does not exist
        missing_seqs.append(a)
    else:
        pass

missing_seqs2 = [value for value in missing_seqs if value != 'nan']

if len(missing_seqs2)>0:
    missing_sequence_export = pd.DataFrame()
    missing_sequence_export['Sequence'] = missing_seqs2
    with open(missing_fragment_files,'w',newline='') as file:
        writer = csv.writer(file)
        missing_sequence_export.to_csv(file,index=False)
        raise ValueError('see exported list of missing fragment lists')
else:
    pass
    
mh_table = pd.DataFrame() #creates a table of the identified MH values
mh_table['Target Sequence'] = target_sequence
mh_table['File Name'] = ion_file_name
mh_table['MH'] = MH_list
mh_table['M'] = mh_table['MH']-H_mass #calculates M from M+H
for b in charges:
    mh_table[str(b)] = (mh_table['M'] + (b*H_mass))/b #from M calculates +2 charge precursor

from user_input import target_list_out_results
result_name2 = target_list_out_results
with open(result_name2,'w',newline='') as file:
    writer = csv.writer(file)
    mh_table.to_csv(file,index=False)
