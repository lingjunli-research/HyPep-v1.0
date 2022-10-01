# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 08:24:12 2022

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


sequence = []
for a in sequence_unfilter:
    if a not in sequence:
        sequence.append(a)

rep1_possible_sequence = []
rep1_scan_no= []
rep1_theo_precursor = []
rep1_act_precursor = []
rep1_ppm_err = []
rep1_prec_act_charge = []
rep1_prec_theo_charge = []
rep1_act_frag_mz = []
rep1_act_frag_charge = []
rep1_intensity = []
rep1_resolution = []
rep1_frag_act_M = []
rep1_frag_theo_mz = []
rep1_frag_theo_M = []
rep1_frag_err = []
rep1_seq2 = []
rep1_ion = []
rep1_frag_mz = []
rep1_frag_M = []
rep1_frag_ch_mz = []
rep1_frag_ch =[]

for b in sequence:
    seq_out = []
    cov_out = []
    sample_out = []
    scan_out = []
    charge_out = []
    upper_only = []
    for c in b:
        remove_lower = lambda text: re.sub('[a-z]', '', text)
        c = remove_lower(c)
        c = ''.join(item for item in c if item.isalpha())
        upper_only.append(c)
    upper_only = [i for i in upper_only if i]
    #locates b/y ion table corresponding to the sequence
    def get_file_names_with_strings(str_list):
        full_list = os.listdir(ion_list_import)
        final_list = [nm for ps in str_list for nm in full_list if ps in nm]
        return final_list
    file_type2 = 'ion_list_'
    file_query2 = file_type2 + b + '.csv' #search query
    b_y_ion = (get_file_names_with_strings([file_query2]))
    b_y_ion_path = ion_list_import + '\\' + b_y_ion[0] 
    prospector = pd.read_csv(b_y_ion_path) #reads ion file
    
    precursor_fragment_matches_filter =  precursor_fragment_matches[precursor_fragment_matches['Possible sequence'] == b]
    #removes entries where duplicates across all lanes are present
    precursor_fragment_matches_filter = precursor_fragment_matches_filter.drop_duplicates() 
    #fragment_merge = precursor_fragment_matches_filter.merge(prospector,left_on='theoretical b/y ion',right_on='mass') #experimental and theoretical dataframes are merged on the basis of the theoretical fragment ion mass
    
    theoretical_fragments = pd.DataFrame()
    theoretical_fragments['Ion identity'] = prospector['ion']
    theoretical_fragments['fragment m/z'] = prospector['mass']
    theoretical_fragments['fragment M'] = (theoretical_fragments['fragment m/z']*1) - (h_mass * 1)
    
    for d in fragment_charges:
        theoretical_fragments[str(d)+' m/z'] = (theoretical_fragments['fragment M'] + (h_mass * d)) / d
    if len(fragment_charges) > 0:
        try:
            theoretical_fragments[['Ion identity', 'ion toss']] = theoretical_fragments['Ion identity'].str.split('-', expand=True)
        except Exception:
            pass
        try:
            theoretical_fragments[['Ion identity', 'ion toss2']] = theoretical_fragments['Ion identity'].str.split('+', expand=True)
        except Exception:
            pass
    theoretical_fragments['b status'] = np.where(theoretical_fragments['Ion identity'].str.contains('b'),True,False)
    theoretical_fragments['y status'] = np.where(theoretical_fragments['Ion identity'].str.contains('y'),True,False)
    theoretical_fragments['b status'] = theoretical_fragments['b status'].astype(int)
    theoretical_fragments['y status'] = theoretical_fragments['y status'].astype(int)
    theoretical_fragments['line status'] = theoretical_fragments['b status'] + theoretical_fragments['y status']

    theoretical_fragments_ion_filtered =  theoretical_fragments[theoretical_fragments['line status'] != 0]
    
    #restructuring table
    merge_prep_theoretical = pd.DataFrame()
    merge_prep_theoretical['Ion identity'] = theoretical_fragments_ion_filtered['Ion identity']
    merge_prep_theoretical['fragment m/z'] = theoretical_fragments_ion_filtered['fragment m/z']
    merge_prep_theoretical['fragment M'] = theoretical_fragments_ion_filtered['fragment M']
    merge_prep_theoretical['sequence'] = b

    for e in fragment_charges:
        merge_prep_theoretical[str(e)+' m/z'] = theoretical_fragments_ion_filtered[str(e)+' m/z']
    
    cm_seq = []
    cm_ion = []
    cm_frag_mz = []
    cm_fragM = []
    cm_ch_mz = []
    cm_ch = []
    
    for f in fragment_charges:
        charge_merge = pd.DataFrame()
        charge_merge['Sequence'] = merge_prep_theoretical['sequence']
        charge_merge['Ion identity'] = merge_prep_theoretical['Ion identity']
        charge_merge['fragment m/z'] = merge_prep_theoretical['fragment m/z']
        charge_merge['fragment M'] = merge_prep_theoretical['fragment M']
        charge_merge[str(f)+' m/z'] = merge_prep_theoretical[str(f)+' m/z']
        charge_merge['frag charge'] = f
        seq = charge_merge['Sequence'].values.tolist()
        for aa in seq:
            cm_seq.append(aa)
        
        ion = charge_merge['Ion identity'].values.tolist()
        for ab in ion:
            cm_ion.append(ab)
        
        mz = charge_merge['fragment m/z'].values.tolist()
        for ac in mz:
            cm_frag_mz.append(ac)
        
        M = charge_merge['fragment M'].values.tolist()
        for ad in M:
            cm_fragM.append(ad)
        
        ch_mz = charge_merge[str(f)+' m/z'].values.tolist()
        for ae in ch_mz:
            cm_ch_mz.append(ae)
            
        ch = charge_merge['frag charge'].values.tolist()
        for af in ch:
            cm_ch.append(af)

    theoretical_values = pd.DataFrame()
    theoretical_values['Sequence'] = cm_seq
    theoretical_values['Ion identity'] = cm_ion
    theoretical_values['fragment m/z'] = cm_frag_mz
    theoretical_values['fragment M'] = cm_fragM
    theoretical_values['frag charge m/z'] = cm_ch_mz
    theoretical_values['theoretical_charge'] = cm_ch
    sequence_merge = precursor_fragment_matches_filter.merge(theoretical_values,left_on='theoretical fragment',right_on='frag charge m/z') #experimental and theoretical dataframes are merged on the basis of the theoretical fragment ion mass
    sequence_merge =  sequence_merge[sequence_merge['charge'] == sequence_merge['theoretical_charge']]
    if len(sequence_merge) > 0:        
        file_name = sample_name + '_' + b + '_match_check.csv'
        file_path = working_directory + '\\' + file_name
        with open(file_path,'w',newline='') as filec:
                writerc = csv.writer(filec)
                sequence_merge.to_csv(filec,index=False) 
    
    pos_seq = sequence_merge['Possible sequence'].values.tolist()
    scan = sequence_merge['Scan #'].values.tolist()
    t_prec = sequence_merge['Theoretical precursor'].values.tolist()
    a_prec = sequence_merge['Precursor'].values.tolist()
    prec_err = sequence_merge['ppm error'].values.tolist()
    p_a_ch = sequence_merge['Precursor_Charge'].values.tolist()
    p_t_ch = sequence_merge['Theoretical precursor charge'].values.tolist()
    f_a_mz = sequence_merge['Fragment m/z'].values.tolist()
    f_a_ch = sequence_merge['charge'].values.tolist()
    inten = sequence_merge['intensity'].values.tolist()
    res = sequence_merge['resolution'].values.tolist()
    f_a_M = sequence_merge['fragment M_x'].values.tolist()
    f_t_mz = sequence_merge['theoretical fragment'].values.tolist()
    f_t_M = sequence_merge['theoretical fragment M'].values.tolist()
    f_err = sequence_merge['Da error'].values.tolist()
    seq3 = sequence_merge['Sequence'].values.tolist()
    ion = sequence_merge['Ion identity'].values.tolist()
    f_mz = sequence_merge['fragment m/z'].values.tolist()
    f_M = sequence_merge['fragment M_y'].values.tolist()
    f_ch = sequence_merge['frag charge m/z'].values.tolist()
    ch = sequence_merge['theoretical_charge'].values.tolist()
    
    for bb in pos_seq:
        rep1_possible_sequence.append(bb)
    for bc in scan:
        rep1_scan_no.append(bc)
    for bd in t_prec:
        rep1_theo_precursor.append(bd)
    for be in a_prec:
        rep1_act_precursor.append(be)
    for bf in prec_err:
        rep1_ppm_err.append(bf)
    for bg in p_a_ch:
        rep1_prec_act_charge.append(bg)
    for bh in p_t_ch:
        rep1_prec_theo_charge.append(bh)
    for bi in f_a_mz:
        rep1_act_frag_mz.append(bi)
    for bh in f_a_ch:
        rep1_act_frag_charge.append(bh)
    for bj in inten:
        rep1_intensity.append(bj)
    for bk in res:
        rep1_resolution.append(bk)
    for bl in f_a_M:
        rep1_frag_act_M.append(bl)
    for bm in f_t_mz:
        rep1_frag_theo_mz.append(bm)
    for bn in f_t_M:
        rep1_frag_theo_M.append(bn)
    for bo in f_err:
        rep1_frag_err.append(bo)
    for bp in seq3:
        rep1_seq2.append(bp)
    for bq in ion:
        rep1_ion.append(bq)
    for br in f_mz:
        rep1_frag_mz.append(br)
    for bs in f_M:
        rep1_frag_M.append(bs)
    for bt in f_ch:
        rep1_frag_ch_mz.append(bt)
    for bu in ch:
        rep1_frag_ch.append(bu)

total_merged = pd.DataFrame()
total_merged['Possible sequence'] = rep1_possible_sequence
total_merged['Scan #'] = rep1_scan_no
total_merged['Theoretical precursor'] = rep1_theo_precursor
total_merged['Actual precursor'] = rep1_act_precursor
total_merged['Precursor ppm error'] = rep1_ppm_err
total_merged['Precursor actual charge'] = rep1_prec_act_charge
total_merged['Precursor theoretical charge'] = rep1_prec_theo_charge
total_merged['Actual fragment m/z'] = rep1_act_frag_mz
total_merged['Actual fragment charge'] = rep1_act_frag_charge
total_merged['Fragment actual M'] = rep1_frag_act_M
total_merged['Fragment theoretical m/z (@ ch1)'] = rep1_frag_theo_mz
total_merged['Fragment theoretical M'] = rep1_frag_theo_M
total_merged['Fragment error'] = rep1_frag_err
total_merged['Ion identity'] = rep1_ion
total_merged['Fragment theoretical m/z'] = rep1_frag_ch_mz 
total_merged['Fragment theoretical charge'] = rep1_frag_ch
total_merged['Fragment intensity'] = rep1_intensity
total_merged['MS2 resolution'] = rep1_resolution 

scans_present = total_merged['Scan #'].values.tolist()
scans_sort = []
for zz in scans_present:
    if zz not in scans_sort:
        scans_sort.append(zz)

for za in scans_sort:
    scan_merged = total_merged[total_merged['Scan #'] == za]
    filtered_seq = scan_merged['Possible sequence'].values.tolist()
       
    seq_nodups = []
    for zb in filtered_seq:
        if zb not in seq_nodups:
            seq_nodups.append(zb)
    for zc in seq_nodups:
        seq_merged = scan_merged[scan_merged['Possible sequence'] == zc]
        seq_merged = seq_merged.drop_duplicates(subset=['Ion identity'])
        
        mid_out = working_directory + '\\' + sample_name + '_' + zc + '_' + str(za) + 'fragments.csv'
        with open(mid_out,'w',newline='') as filec:
                                    writerc = csv.writer(filec)
                                    seq_merged.to_csv(filec,index=False)
        remove_lower = lambda text: re.sub('[a-z]', '', text)
        zc_simple = remove_lower(zc)
        zc_simple = re.sub(r"[-()\"#/@;:<>{}`+=~|.!?,]", "", zc_simple)
        no_ions = len(zc_simple)-1
        seq_list1 = list(zc_simple)
        seq_list2 = list(zc_simple)
        seq_list3 = list(zc_simple)
        
        seq_merged['b status'] = np.where(seq_merged['Ion identity'].str.contains('b'),True,False)
        seq_merged['y status'] = np.where(seq_merged['Ion identity'].str.contains('y'),True,False)

        seq_merged_b = seq_merged[seq_merged['b status'] == True]
        seq_merged_y = seq_merged[seq_merged['y status'] == True]
        
        b_ions = seq_merged_b['Ion identity'].values.tolist()
        y_ions = seq_merged_y['Ion identity'].values.tolist()
        
        b_ions_no_dups = []
        for i in b_ions:
            if i not in b_ions_no_dups:
                b_ions_no_dups.append(i)

        y_ions_no_dups = []
        for j in y_ions:
            if j not in y_ions_no_dups:
                y_ions_no_dups.append(j)
        
        b_ion_loc = [k[1:] for k in b_ions_no_dups]
        for l in range(0,len(b_ion_loc)):
            b_ion_loc[l] = int(b_ion_loc[l])
        for m in b_ion_loc:
            seq_list1[m] = True
            
        y_ion_loc = [n[1:] for n in y_ions_no_dups]
        for o in range(0,len(y_ion_loc)):
            y_ion_loc[o] = int(y_ion_loc[o])
        for p in y_ion_loc:
            seq_list3[p] = True
            
        pos_ions_b = []
        pos_ions_y = []
        num_b = []
        
        for q in range(0,len(zc_simple)):
            pos_ions_b.append('b')
            pos_ions_y.append('y')
            num_b.append(q)
        
        num_y = num_b[::-1]
        #creates result table
        result = pd.DataFrame()
        result['b present'] = seq_list1
        result['b'] = pos_ions_b
        result['b loc'] = num_b
        result['sequence'] = seq_list2
        result['y'] = pos_ions_y
        result['y loc'] = num_y
        result['y present'] = seq_list3
        result.loc[(result['b present'] != True),'b present']=False #makes any not true values False
        result.loc[(result['y present'] != True),'y present']=False
        result['b loc'] = result['b loc'].map(str)
        result['b ions'] = [''.join(i) for i in zip(result["b"].map(str),result["b loc"])] #concatonates "b" and the number location
        result['y loc'] = result['y loc'].map(str)
        result['y ions'] = [''.join(i) for i in zip(result["y"].map(str),result["y loc"])]

        #formats result chart
        result_exp = pd.DataFrame()
        result_exp['b ion present'] = result['b present']
        result_exp['b ion'] = result['b ions']
        result_exp['sequence'] = result['sequence']
        result_exp['y ion'] = result['y ions']
        result_exp['y ion present'] = result['y present']

        res_out = working_directory + '\\' + sample_name + '_' + zc + '_' + str(za) + '_.txt'
        with open(res_out,'w',newline='') as filec:
                                    writerc = csv.writer(filec)
                                    result_exp.to_csv(filec,index=False)

        cov_result = result_exp
        cov_result['b ion present'] = cov_result['b ion present'].astype(int) #changes T/F to 0/1
        cov_result['y ion present'] = cov_result['y ion present'].astype(int)
        cov_result['Sum'] = abs(cov_result['y ion present'] + cov_result['b ion present']) #anything greater than 0 in the sum would indicate the ion is present
        total_col = cov_result['Sum']
        total = total_col[total_col != 0].count() #count all cells greater than zero
        coverage = (total / ((len(zc_simple))-1)) * 100 #percentage calculation
        coverage = round(coverage,4)
        coverage = str(coverage)
        
        seq_out.append(zc)
        cov_out.append(coverage)
        sample_out.append(sample_name)
        scan_out.append(za)
   
    summary = pd.DataFrame()
    summary['coverage'] = cov_out
    summary['sequence'] = seq_out
    summary['sample'] = sample_out
    summary['scan'] = scan_out

    out_path = final_dir + '\\' + sample_name + '_all_coverage.csv'
    with open(out_path,'a',newline='') as filec:
                            writerc = csv.writer(filec)
                            summary.to_csv(filec,index=False)

file_correct_path =  final_dir + '\\' + sample_name + '_all_coverage.csv'
file_correct = pd.read_csv(file_correct_path)
file_filter = file_correct[file_correct['coverage'] != 'coverage']

file_filter['coverage'] = file_filter['coverage'].astype(float)
file_filter_sorted = file_filter.sort_values(['coverage'])
from user_input import sequence_coverage_out_path
out_path = sequence_coverage_out_path
with open(out_path,'w',newline='') as filec:
                        writerc = csv.writer(filec)
                        file_filter_sorted.to_csv(filec,index=False)
