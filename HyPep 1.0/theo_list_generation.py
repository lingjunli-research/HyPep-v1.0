# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 13:49:40 2022

@author: lawashburn
"""

import pandas as pd
import csv
import os
from datetime import datetime
now = datetime.now()

print('start time:', datetime.now())

from user_input import formatted_ion_list_path
from user_input import filtered_fragment_list_path
from user_input import target_list_out_results
from user_input import sample_name

from user_input import fragment_charges
ion_directory = formatted_ion_list_path
working_directory = filtered_fragment_list_path
unique_amm_sequences = target_list_out_results


H_mass = 1.00784

def get_file_names_with_strings_list(str_list): #definition for finding a file containing a string in filename in specified directory
    full_list = os.listdir(ion_directory)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

amm_sequences = pd.read_csv(unique_amm_sequences)
target_sequences = amm_sequences['Target Sequence'].values.tolist() #list of sequences

MH_list = []
ion_file_name = []
target_sequence = []
theo_list_record = pd.DataFrame()
for a in target_sequences: #for each sequence in the sequence list
    file_type = 'Theoretical_b_y_ion_list_'
    ion_query = file_type + a + '.csv' #search for ion list pertaining to the sequence
    ion_file = (get_file_names_with_strings_list([ion_query])) #search for the file based on query 
    if len(ion_file)>1: #throws an error if more than one ion file exists
        raise ValueError('More than one ion list located for:',a)
    if len(ion_file)==0: #throws an error if ion file does not exist
        raise ValueError('No ion list located for:',a)
    if len(ion_file) == 1: #ensures only one file matched
        for b in ion_file:
            ion_path = ion_directory + '\\' + b
            ion_read = pd.read_csv(ion_path) #reads identified .csv
            ion_read['peptide'] = a
            
            theoretical_list = pd.DataFrame()
            theoretical_list['Peptide'] = ion_read['peptide']
            theoretical_list['1'] = ion_read['mass']
            
            for cc in fragment_charges:
                theoretical_list[str(cc)] = ((theoretical_list['1']-H_mass)+(cc*H_mass))/cc

            theoretical_list = theoretical_list.dropna()
            
            file_name = 'Theoretical_b_y_fragment_list_' + a + '.csv'
            file_path = working_directory + '\\' + file_name
            
            with open(file_path,'w',newline='') as file:
                writer = csv.writer(file)
                theoretical_list.to_csv(file,index=False)
            theo_list_record = theo_list_record.append(theoretical_list)