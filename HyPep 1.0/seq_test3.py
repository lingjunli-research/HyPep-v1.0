# -*- coding: utf-8 -*-
"""
Created on Tue May 10 13:25:16 2022

@author: lawashburn
"""

import os
import csv
import pandas as pd
import numpy as np

import re
from datetime import datetime
now = datetime.now()

from user_input import formatted_ion_list_path
from user_input import all_fragment_matches_report_out
from user_input import sequence_working_directory_path
from user_input import sequence_final_directory_path
from user_input import sample_name
from user_input import precursor_error_cutoff
from user_input import fragment_error_cutoff
from user_input import fragment_charges
ion_list_import = formatted_ion_list_path #path to directory with b/y target ions
prelim_results = all_fragment_matches_report_out #path to output data from previous step

working_directory = sequence_working_directory_path #path to directory for all intermediate data
final_dir =sequence_final_directory_path#path to directory for final, processed data
#check_path = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\08_sequence_coverage_calculation\check"


error_marg = precursor_error_cutoff #+/- ppm, for precursor
error_ch1 = fragment_error_cutoff #+/- Da, for fragment

h_mass = 1.00784

precursor_fragment_matches = pd.read_csv(prelim_results)

sequence_unfilter = precursor_fragment_matches['Possible sequence'].values.tolist()



    
file_correct_path =  final_dir + '\\' + sample_name + '_all_coverage.csv'
file_correct = pd.read_csv(file_correct_path, on_bad_lines='skip')
file_filter = file_correct[file_correct['coverage'] != 'coverage']

from user_input import sequence_coverage_out_path
out_path = sequence_coverage_out_path
with open(out_path,'w',newline='') as filec:
                        writerc = csv.writer(filec)
                        file_filter.to_csv(filec,index=False)