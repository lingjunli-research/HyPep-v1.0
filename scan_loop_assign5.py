# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 16:06:08 2022

@author: lawashburn
"""

import pandas as pd
import csv
from user_input import sequence_coverage_out_path
from user_input import SHS_results_path
from user_input import target_list_out_results
from user_input import looping_working_directory_path
from user_input import looping_final_directory_path
from user_input import sample_name
from user_input import number_loops
coverage_report_path = sequence_coverage_out_path
shs_priority_path = SHS_results_path
inclusion_list = target_list_out_results
final_directory = looping_final_directory_path
working_directory = looping_working_directory_path

#check_path = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\09_PSM_assign_looping\check"


loops = number_loops

unfiltered_inclusion_list = pd.read_csv(inclusion_list)
unfiltered_sequence_list = unfiltered_inclusion_list['Target Sequence'].values.tolist()

priority_inclusion_list = shs_priority_path
priority_sequence_list = priority_inclusion_list['Sequence'].values.tolist()

coverage_report = pd.read_csv(coverage_report_path)
coverage_report_filtered = coverage_report.drop_duplicates(subset=['coverage','sequence','sample','scan'])

results_report = pd.DataFrame()
banned_seqs = []
banned_scans = []
for b in range(0,loops): #shs priority
    coverage_seq_filter = coverage_report_filtered[~coverage_report_filtered['sequence'].isin(banned_seqs)]
    coverage_scan_filter = coverage_seq_filter[~coverage_seq_filter['scan'].isin(banned_scans)]
    current_seqs_dups = coverage_scan_filter['sequence'].values.tolist()
    
    current_seq = []
    for c in current_seqs_dups:
        if c not in current_seq:
            current_seq.append(c)
    
    search_list = []
    for z in priority_sequence_list:
        if z in current_seq:
            search_list.append(z)
        else:
            pass
    
    for y in current_seq:
        if y not in search_list:
            search_list.append(y)
    
    for e in search_list:
        coverage_seq_filter_update = coverage_scan_filter[~coverage_report_filtered['sequence'].isin(banned_seqs)]
        coverage_scan_filter_update = coverage_seq_filter_update[~coverage_seq_filter_update['scan'].isin(banned_scans)]
        coverage_filtered = coverage_scan_filter_update[coverage_scan_filter_update['sequence'] == e]

        if len(coverage_filtered)>1:
            scans_present = coverage_filtered['scan'].values.tolist()
            scan_counts = []
            for f in scans_present:
                count = scans_present.count(f)
                scan_counts.append(count)
            coverage_filtered['scan count'] = scan_counts
            coverage_filtered_sorted = coverage_filtered.sort_values(by=['coverage','scan count'], ascending=[False,True])
            matched_scan = coverage_filtered_sorted.drop_duplicates(subset=['sequence'],keep='first')
            results_report = results_report.append(matched_scan)
            
            sequence_remove = matched_scan['sequence'].iloc[0]
            banned_seqs.append(sequence_remove)
            
            scan_remove = matched_scan['scan'].iloc[0]
            banned_scans.append(scan_remove)
        if len(coverage_filtered) == 1:
            results_report = results_report.append(coverage_filtered)
            
            sequence_remove = coverage_filtered['sequence'].iloc[0]
            banned_seqs.append(sequence_remove)
            
            scan_remove = coverage_filtered['scan'].iloc[0]
            banned_scans.append(scan_remove)
        else:
            banned_seqs.append(e)

file_name =  sample_name + '_psm_matches_20220521.csv'
file_out_path = final_directory + '\\' +  file_name
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        results_report.to_csv(filec,index=False)