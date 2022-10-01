# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 15:48:58 2022

@author: lawashburn
"""

import pandas as pd
import csv

coverage_report_path = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\08_sequence_coverage_calculation\final_directory\Untarget_Brain1_20ppm_all_coverage_formatted.csv"
shs_priority_path = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\09_PSM_assign_looping\priority_list.csv"
inclusion_list = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\07inclusion_list\Brain_inclusion_list.csv"
final_directory = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\09_PSM_assign_looping\final_directory"
working_directory = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\09_PSM_assign_looping\working_directory"
check_path = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\09_PSM_assign_looping\check"
sample_name = 'Untarget_Brain_1'

loops = 5

unfiltered_inclusion_list = pd.read_csv(inclusion_list)
unfiltered_sequence_list = unfiltered_inclusion_list['Target Sequence'].values.tolist()

priority_inclusion_list = pd.read_csv(shs_priority_path)
priority_sequence_list = priority_inclusion_list['Priority sequence'].values.tolist()

for a in unfiltered_sequence_list:
    if a not in priority_sequence_list:
        priority_sequence_list.append(a)

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
    
    priority_seq_focused = []
    for d in priority_sequence_list:
        if d in current_seq:
            priority_seq_focused.append(d)

    for e in priority_seq_focused:
        coverage_seq_filter_update = coverage_scan_filter[~coverage_report_filtered['sequence'].isin(banned_seqs)]
        coverage_scan_filter_update = coverage_seq_filter_update[~coverage_seq_filter_update['scan'].isin(banned_scans)]
        coverage_filtered = coverage_scan_filter_update[coverage_scan_filter_update['sequence'] == e]

        file_name =  sample_name + str(b) + '_loops.csv'
        file_out_path = check_path + '\\' +  file_name
        with open(file_out_path,'w',newline='') as filec:
                writerc = csv.writer(filec)
                coverage_filtered.to_csv(filec,index=False)

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

    list_update = []
    for h in unfiltered_sequence_list:
        if h not in banned_seqs:
            list_update.append(h)
    priority_seq_focused = list_update

results_report = results_report.sort_values(by=['coverage'],ascending=False)


file_name =  sample_name + '_psm_matches.csv'
file_out_path = final_directory + '\\' +  file_name
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        results_report.to_csv(filec,index=False)
                