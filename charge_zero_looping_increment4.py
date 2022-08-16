# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 09:36:19 2022

@author: lawashburn
"""

import os
import csv
import pandas as pd
import numpy as np
from datetime import datetime
now = datetime.now()

spectra_import = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220407\Untarget_Raw_Files\Formatted\Brain_1_ms2_output_list.txt" #path to directory containing spectra
fragment_list_import = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220428\looping_zero_charge_overhaul\filtered_fragment_lists" #path to directory with list of target fragment ions
target_list_import = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220428\looping_zero_charge_overhaul\inclusion_lists\Brain_inclusion_list.csv"#path to directory with inclusion lists
working_directory = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\charge0_to_1\wd"
final_dir =r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\charge0_to_1\fd" #path to directory for final, processed data
spectra_dir = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\charge0_to_1\spectra_dir"
used_fragments_dir = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\charge0_to_1\used_fragments_dir"
#precursor_doc = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220428\looping_zero_charge_overhaul\zero_precursor_fragment_matching\working_directory\Untarget_Brain1_20ppmprecursor_matches.csv"
#charge_1 = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220425\zero_charge_test\Untarget_Brain1_20ppm_zero_1fragment_matches.csv"

sample_name = 'Untarget_Brain1_20ppm'

error_precursor = 20 #+/- ppm, for precursor
error_fragment = 0.02 #+/- Da, for fragment ion, charge state 1

precursor_charges = [1,2,3,4,5,6,7,8]
fragment_charges = [1,2,3,4]

h_mass = 1.00784

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
spectra_value['Scan #'] = spectra_value['Scan #'].round(0)
spectra_value['Fragment m/z'] = spectra_value['Fragment m/z'].round(4)
print(spectra_value)
spectra_value = spectra_value.drop(spectra_value[spectra_value['Precursor_Charge']==0].index) #remove charges equal to 0
spectra_value_zeroes = spectra_value.drop(spectra_value[spectra_value['charge']!=0].index) #remove charges equal to 0
spectra_value_zeroes = spectra_value_zeroes.drop(spectra_value_zeroes[spectra_value_zeroes['intensity']<100].index) #remove intensity values less than 100
print(spectra_value_zeroes)
def get_file_names_with_strings(str_list):
    full_list = os.listdir(working_directory)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

banned_scans = []
banned_fragments = []

all_used_fragments = pd.DataFrame()
for zz in remaining_charges:
    zs = str(zz)
    file_query = '_zero_reassign_fragment_matches.csv'
    fragment_list = (get_file_names_with_strings([file_query]))
    for a in fragment_list:
        fragment_list_path = working_directory + '\\' + a
        fragment_read = pd.read_csv(fragment_list_path)
        fragment_read = fragment_read.drop_duplicates()

        all_used_fragments = all_used_fragments.append(fragment_read)
        unavailable_scans = fragment_read['Scan #'].values.tolist()
        unavailable_fragments = fragment_read['Fragment m/z'].values.tolist()
        
        for b in unavailable_scans:
            banned_scans.append(b)
        
        for c in unavailable_fragments:
            banned_fragments.append(c)
    
    working_spectra = spectra_value_zeroes[~spectra_value_zeroes['Scan #'].isin(banned_scans)]
    working_spectra = spectra_value_zeroes[~spectra_value_zeroes['Scan #'].isin(banned_scans)]

