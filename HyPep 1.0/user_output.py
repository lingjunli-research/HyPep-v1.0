# -*- coding: utf-8 -*-
"""
Created on Mon May  9 08:42:26 2022

@author: lawashburn
"""
from user_input import folder_generation_status
print('Step 1: Importing spectra')
from TopFD_rawConverter_Combined import precursor
print('Step 2: Preliminary AMM Matching')
from Prelim_AMM import prelim_AMM_out
print('Step 3: Comparing AMM and SHS')
from AMM_v_SHS import seq_overview
print('Step 4: Generating target list')
from target_list_create import mh_table
print('Step 5: Generating b and y target ion lists')
from ion_list_format import ion_mod_full
print('Step 6: Creating theoretical ion lists')
from theo_list_generation import theo_list_record
print('Step 7: Precursor matching')
from precursor_fragment_matching import precursor_matches
print('Step 8: Fragment matching')
from precursor_fragment_matching import fragment_matches
print('Step 9: Assigning 0 charges, pt 1')
from charge_zero_looping_increment import precursor_matches
print('Step 10: Assigning 0 charges, pt 2')
from charge_zero_looping_increment import fragment_matches
print('Step 11: Assigning 0 charges, pt 3')
from charge_zero_looping_increment3 import precursor_matches
print('Step 12: Assigning 0 charges, pt 4')
from charge_zero_looping_increment3 import fragment_matches
print('Step 13: Combining all matched fragments')
from fragment_combine import all_fragments
print('Step 14: Calculating sequence coverage')
from seq_test2 import file_filter
print('Step 15: Assigning scans to sequences')
from scan_loop_assign5 import results_report