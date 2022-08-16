# -*- coding: utf-8 -*-
"""
Created on Tue May  3 13:28:20 2022

@author: lawashburn
"""

import os
import pickle
#import main
#from main import FDR_filtered_target_list_for_AMM
from hypep_import import out_path
from hypep_import import topfd_path
from hypep_import import rawconverter_path
from hypep_import import db_path
from hypep_import import pp_path
from hypep_import import sample_ID
from hypep_import import pro_mex
from hypep_import import p_err
from hypep_import import f_err
from hypep_import import max_p_charge
from hypep_import import max_f_charge
from hypep_import import loop


output_folder = out_path #folder in which all output directories will be generated

with open('SHS_ouput.pkl', 'rb') as file_h:      
    # Call load method to deserialze
    SHS_results_path = pickle.load(file_h)

SHS_results_path = SHS_results_path.rename(columns={"accession": "Accession Number", "Precursor_Sequence": "Sequence"})

SHS_results_path_for_AMM = r"SHS_pretend_results_empty.csv"

#SHS_results_path = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_GUI_linked_final_v2\input_files\SHS_pretend_results.csv"
#SHS_results_path = FDR_filtered_target_list_for_AMM ### a list not a file anymore
ion_list_directory_path = pp_path

sample_name = sample_ID
promex_cutoff = float(pro_mex)
precursor_error_cutoff = float(p_err)
fragment_error_cutoff = float(f_err)
precursor_charges = []
max_p_charge = int(max_p_charge)
max_f_charge = int(max_f_charge)
for a in range(1,max_p_charge+1):
    precursor_charges.append(a)

fragment_charges = []
for a in range(1,max_f_charge+1):
    fragment_charges.append(a)
number_loops = int(loop)

rawconverter_formatting_out_path = output_folder + '\\' + 'rawconverter_MS2_formatted'
topFD_rawconverter_combined_out_path = output_folder + '\\' + 'topFD_rawconverter_combined'
prelim_AMM_out_path = output_folder + '\\' + 'prelim_AMM_results'
AMM_v_SHS_out_path = output_folder + '\\' + 'AMM_v_SHS'
target_inclusion_list_path = output_folder + '\\' + 'inclusion_lists'
formatted_ion_list_path = output_folder + '\\' + 'ion_list_formatted'
formatted_ion_list_report_path = output_folder + '\\' + 'ion_list_formatted\\final_directory'
filtered_fragment_list_path = output_folder + '\\' + 'filtered_fragment_lists'
filtered_fragment_list_report_path = output_folder + '\\' + 'filtered_fragment_lists\\final_directory'
non0_precursor_fragment_matches_out = output_folder + '\\' + 'non0_precursor_fragment_matches'
assign0_working_directory_path = output_folder + '\\assign_charge_0\\working_directory'
assign0_final_directory_path = output_folder + '\\assign_charge_0\\final_directory'
combined_fragments_out = output_folder + '\\' + 'combined_fragment_matches'
sequence_working_directory_path = output_folder + '\\sequence_assign\\working_directory'
sequence_final_directory_path = output_folder + '\\sequence_assign\\final_directory'
looping_working_directory_path = output_folder + '\\scan_loop_assign\\working_directory'
looping_final_directory_path = output_folder + '\\scan_loop_assign\\working_directory'

results_export = output_folder + '\\final_results'

if not os.path.exists(rawconverter_formatting_out_path):
    os.makedirs(rawconverter_formatting_out_path)
if not os.path.exists(topFD_rawconverter_combined_out_path):
    os.makedirs(topFD_rawconverter_combined_out_path)
if not os.path.exists(prelim_AMM_out_path):
    os.makedirs(prelim_AMM_out_path)
if not os.path.exists(AMM_v_SHS_out_path):
    os.makedirs(AMM_v_SHS_out_path) 
if not os.path.exists(formatted_ion_list_path):
    os.makedirs(formatted_ion_list_path) 
if not os.path.exists(formatted_ion_list_report_path):
    os.makedirs(formatted_ion_list_report_path)     
if not os.path.exists(target_inclusion_list_path):
    os.makedirs(target_inclusion_list_path) 
if not os.path.exists(filtered_fragment_list_path):
    os.makedirs(filtered_fragment_list_path) 
if not os.path.exists(filtered_fragment_list_report_path):
    os.makedirs(filtered_fragment_list_report_path) 
if not os.path.exists(non0_precursor_fragment_matches_out):
    os.makedirs(non0_precursor_fragment_matches_out) 
if not os.path.exists(assign0_working_directory_path):
    os.makedirs(assign0_working_directory_path)
if not os.path.exists(assign0_final_directory_path):
    os.makedirs(assign0_final_directory_path)
if not os.path.exists(combined_fragments_out):
    os.makedirs(combined_fragments_out)
if not os.path.exists(sequence_working_directory_path):
    os.makedirs(sequence_working_directory_path)
if not os.path.exists(sequence_final_directory_path):
    os.makedirs(sequence_final_directory_path)
if not os.path.exists(looping_working_directory_path):
    os.makedirs(looping_working_directory_path)
if not os.path.exists(looping_final_directory_path):
    os.makedirs(looping_final_directory_path)
if not os.path.exists(results_export):
    os.makedirs(results_export)

rawconverter_formatted_out_file_name_path_txt = rawconverter_formatting_out_path + '\\' + sample_name + '_rawconverter_MS2_formatted.txt'
rawconverter_formatted_out_file_name_path_csv = rawconverter_formatting_out_path + '\\' + sample_name + '_rawconverter_MS2_formatted.csv'
top_FD_combined_out_file_name_path = topFD_rawconverter_combined_out_path + '\\' + sample_name + '_combined_spectrum.csv'
prelim_AMM_out_file_name_path = prelim_AMM_out_path + '\\' + sample_name + '_AMM_results.csv'
AMM_v_SHS_out_results = AMM_v_SHS_out_path + '\\' + sample_name + '_unique_AMM.csv'
missing_fragment_files = target_inclusion_list_path + '\\' + sample_name + '_missing_ion_list.csv'
target_list_out_results = target_inclusion_list_path + '\\' + sample_name + '_inclusion_list.csv'
formatted_ion_list_prefix = formatted_ion_list_path + '\\Theoretical_b_y_ion_list_'
formatted_ion_list_report = formatted_ion_list_report_path + '\\Theoretical_b_y_ion_list_report.csv'
formatted_fragment_list_prefix = filtered_fragment_list_path + '\\Theoretical_b_y_fragment_list_'
filtered_fragment_list_report = filtered_fragment_list_report_path + '\\Theoretical_b_y_fragment_list_report.csv'
non0_precursor_matches_out = non0_precursor_fragment_matches_out + '\\' + sample_name + 'precursor_matches.csv'
non0_fragment_matches_out = non0_precursor_fragment_matches_out + '\\' + sample_name + 'fragment_matches.csv'
assign0_precursor_matches_out = assign0_working_directory_path + '\\' + sample_name + 'precursor_matches.csv'
assign0_1_fragment_matches_out = assign0_working_directory_path + '\\' + sample_name + '_1_zero_reassign_fragment_matches.csv'
all_fragment_matches_report_out = combined_fragments_out + '\\' + sample_name + '_all_fragments.csv'
sequence_coverage_out_path = sequence_final_directory_path + '\\' + sample_name + '_all_coverage_formatted.csv'
results_out_path = results_export + '\\' + sample_name + '_combined_AMM_SHS_final_report.csv'
psm_matches_path = looping_final_directory_path + '\\' + sample_name + '_psm_matches_20220521.csv'

AMM_out_path = results_export + '\\' + sample_name + '_AMM_final_report.csv'
SHS_out_path = results_export + '\\' + sample_name + '_SHS_final_report.csv'

folder_generation_status = True