# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 11:54:58 2022

@author: lawashburn
"""

import pandas as pd
import csv
import os
from user_input import assign0_working_directory_path
zero_dir = assign0_working_directory_path
from user_input import non0_fragment_matches_out
normal_df = non0_fragment_matches_out
from user_input import combined_fragments_out
output_dir = combined_fragments_out
from user_input import sample_name

def get_file_names_with_strings(str_list):
    full_list = os.listdir(zero_dir)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

file_query = '_zero_reassign_fragment_matches.csv'
fragment_list = (get_file_names_with_strings([file_query]))

all_fragments = pd.DataFrame()

for a in fragment_list:
    a_path = zero_dir + '\\' + a
    a_df = pd.read_csv(a_path)
    
    all_fragments = all_fragments.append(a_df)

normal_df = pd.read_csv(normal_df)
all_fragments = all_fragments.append(normal_df)


from user_input import all_fragment_matches_report_out
file_out_path = all_fragment_matches_report_out
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        all_fragments.to_csv(filec,index=False)