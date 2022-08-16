# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 11:08:29 2022

@author: lawashburn
"""

import os
import pickle
import sys


print('SHS: Raw Target Database Construction')
from raw_target_run import raw_target_run
print('SHS: Raw Decoy Database Construction')
from raw_decoy_run import raw_decoy_run
print('SHS: Identification Mode')
from driver_identification_mode import dummy_v
print('SHS export')
from SHS_bridge_AMM import FDR_filtered_target_list_for_AMM
from user_input import folder_generation_status
print('Formatting .MS2 output')
from RawConverter_Formattingv2 import SIM_result
print('AMM: Importing spectra')
from TopFD_rawConverter_Combined import precursor
print('AMM: Preliminary AMM Matching')
from Prelim_AMM import prelim_AMM_out
print('AMM: Comparing AMM and SHS')
from AMM_v_SHS import seq_overview
print('AMM: Generating target list')
from target_list_create import mh_table
print('AMM: Generating b and y target ion lists')
from ion_list_format import ion_mod_full
print('AMM: Creating theoretical ion lists')
from theo_list_generation import theo_list_record
print('AMM: Precursor matching')
from precursor_fragment_matching import precursor_matches
print('AMM: Fragment matching')
from precursor_fragment_matching import fragment_matches
print('AMM: Assigning 0 charges, pt 1')
from charge_zero_looping_increment import precursor_matches
print('AMM: Assigning 0 charges, pt 2')
from charge_zero_looping_increment import fragment_matches
print('AMM: Assigning 0 charges, pt 3')
from charge_zero_looping_increment3 import precursor_matches
print('AMM: Assigning 0 charges, pt 4')
from charge_zero_looping_increment3 import fragment_matches
print('AMM: Combining all matched fragments')
from fragment_combine import all_fragments
print('AMM: Calculating sequence coverage')
from seq_test2 import file_filter
print('AMM: Assigning scans to sequences')
from scan_loop_assign5 import results_report
print('Finalizing results')
from final_report_generation import formatted_table
print('Analysis complete: results have been exported to the final results directory within selected output directory')