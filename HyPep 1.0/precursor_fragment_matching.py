# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 13:34:48 2022

@author: lawashburn
"""

import os
import csv
import pandas as pd
import numpy as np
from datetime import datetime
now = datetime.now()

from user_input import top_FD_combined_out_file_name_path
from user_input import filtered_fragment_list_path
from user_input import target_list_out_results
from user_input import non0_precursor_fragment_matches_out
from user_input import sample_name
from user_input import rawconverter_formatted_out_file_name_path_txt
spectra_import = rawconverter_formatted_out_file_name_path_txt #path to directory containing spectra
fragment_list_import = filtered_fragment_list_path #path to directory with list of target fragment ions
target_list_import = target_list_out_results#path to directory with inclusion lists
final_dir =  non0_precursor_fragment_matches_out #path to directory for final, processed data


from user_input import precursor_error_cutoff
from user_input import fragment_error_cutoff
from user_input import precursor_charges
from user_input import fragment_charges

error_precursor = precursor_error_cutoff #+/- ppm, for precursor
error_fragment = fragment_error_cutoff #+/- Da, for fragment ion, charge state 1


h_mass = 1.00784

print('loading files', datetime.now())

spectra_read= pd.read_csv(spectra_import, sep=" ",skiprows=[0], names= ["m/z", "resolution", "charge", "intensity","MS2",'scan_number','precursor_charge','empty'])
spectra_value = pd.DataFrame()
spectra_value['Fragment m/z'] = spectra_read['m/z']
spectra_value['resolution'] = spectra_read['resolution']
spectra_value['charge'] = spectra_read['charge']
spectra_value['intensity'] = spectra_read['intensity']
spectra_value['Precursor'] = spectra_read['MS2']
spectra_value['Scan #'] = spectra_read['scan_number']
spectra_value['Precursor_Charge'] = spectra_read['precursor_charge']
spectra_value = spectra_value.drop(spectra_value[spectra_value['Precursor_Charge']==0].index) #remove charges equal to 0
spectra_value = spectra_value.drop(spectra_value[spectra_value['charge']==0].index) #remove charges equal to 0
spectra_value = spectra_value.drop(spectra_value[spectra_value['intensity']<100].index) #remove intensity values less than 100
target_list = pd.read_csv(target_list_import)
sequence_targets = target_list['Target Sequence'].values.tolist()

MS1_sequence_storage = []
MS1_file_type_storage = []
MS1_scan_storage = []
MS1_theo_precursor_store = []
MS1_actual_precursor_store = []
MS1_precursor_theo_charge = []
MS1_precursor_actual_charge = []
MS2_fragment_charge = []
MS2_fragment_mz = []
MS1_intensity = []
MS1_resolution = []
error_storage = []

for a in sequence_targets:
    target_list_filtered = target_list[target_list['Target Sequence'] == a]
    
    for b in precursor_charges:
        precursor_target = target_list_filtered[str(b)].values.tolist()
        spectra_value_filtered = spectra_value[spectra_value['Precursor_Charge'] == b]
        
        for c in precursor_target:
            spectra_value_filtered['ppm error'] =  ((abs(spectra_value_filtered['Precursor'] - c))/c) * 1E6
            spectra_value_filtered = spectra_value_filtered[spectra_value_filtered['ppm error'] <= error_precursor]
            spectra_value_filtered['Possible sequence'] = a
            spectra_value_filtered['Theoretical precursor'] = c
            spectra_value_filtered['Theoretical precursor charge'] = b
            
            if len(spectra_value_filtered) > 0:
                
                sequence = spectra_value_filtered['Possible sequence'].values.tolist()
                scan = spectra_value_filtered['Scan #'].values.tolist()
                theo_precursor = spectra_value_filtered['Theoretical precursor'].values.tolist()
                actual_precursor = spectra_value_filtered['Precursor'].values.tolist()
                precursor_error = spectra_value_filtered['ppm error'].values.tolist()
                precursor_theretical_charge = spectra_value_filtered['Theoretical precursor charge'].values.tolist()
                actual_precursor_charge = spectra_value_filtered['Precursor_Charge'].values.tolist()
                fragment_charge = spectra_value_filtered['charge'].values.tolist()
                fragment_mz = spectra_value_filtered['Fragment m/z'].values.tolist()
                intensity = spectra_value_filtered['intensity'].values.tolist()
                resolution = spectra_value_filtered['resolution'].values.tolist()
                
                for d in sequence:
                    MS1_sequence_storage.append(d)
                for e in scan:
                    MS1_scan_storage.append(e)
                for f in theo_precursor:
                    MS1_theo_precursor_store.append(f)
                for g in actual_precursor:
                    MS1_actual_precursor_store.append(g)
                for h in precursor_error:
                    error_storage.append(h)
                for i in precursor_theretical_charge:
                    MS1_precursor_theo_charge.append(i)
                for j in actual_precursor_charge:
                    MS1_precursor_actual_charge.append(j)
                for k in fragment_charge:
                    MS2_fragment_charge.append(k)
                for l in fragment_mz:
                    MS2_fragment_mz.append(l)
                for m in intensity:
                    MS1_intensity.append(m)
                for n in resolution:
                    MS1_resolution.append(n)

precursor_matches = pd.DataFrame()
                    
precursor_matches['Possible Sequence'] = MS1_sequence_storage
precursor_matches['Scan #'] = MS1_scan_storage
precursor_matches['Theoretical precursor'] = MS1_theo_precursor_store
precursor_matches['Actual precursor'] = MS1_actual_precursor_store
precursor_matches['ppm error'] = error_storage
precursor_matches['Precursor theoretical charge'] = MS1_precursor_theo_charge
precursor_matches['Precursor actual charge'] = MS1_precursor_actual_charge
precursor_matches['Fragment m/z'] = MS2_fragment_mz
precursor_matches['Fragment ion charge'] = MS2_fragment_charge
precursor_matches['Fragment intensity'] = MS1_intensity
precursor_matches['Fragment resolution'] = MS1_resolution

from user_input import non0_precursor_matches_out
file_out_path = non0_precursor_matches_out
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        precursor_matches.to_csv(filec,index=False)

updated_sequence_targets = precursor_matches['Possible Sequence'].values.tolist()

sequence_no_dups = []
for o in updated_sequence_targets:
    if o not in sequence_no_dups:
        sequence_no_dups.append(o)

def get_file_names_with_strings(str_list):
    full_list = os.listdir(fragment_list_import)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

frag_possible_seqeuence_store = []
frag_scan_store = []
frag_theo_precursor_store = []
frag_act_precursor_store = []
frag_precursor_err_store = []
frag_precursor_act_charge_store = []
frag_precursor_theo_charge_store = []
frag_mz_store = []
frag_ion_charge_store = []
frag_intensity_store = []
frag_resoltuion_store = []
frag_act_M_store = []
frag_theo_mz_store = []
frag_theo_M_store = []
Da_err_store = []

for a in fragment_charges:
    for p in sequence_no_dups:
        file_query = 'Theoretical_b_y_fragment_list_' + p + '.csv'
        fragment_list = (get_file_names_with_strings([file_query]))
        for q in fragment_list:
            fragment_list_path = fragment_list_import + '\\' + q
            fragment_database = pd.read_csv(fragment_list_path)
            mass_str = str(a)
            #fragment_database_format = fragment_database.rename({'Mass (1)': '1','Mass (2)': '2','Mass (3)': '3','Mass (4)': '4','Mass (5)': '5','Mass (6)': '6' }, axis=1)
            fragment_database_format=fragment_database
            seq_exp = precursor_matches[precursor_matches['Possible Sequence'] == p]  
            theo_fragment = fragment_database_format[str(a)].values.tolist()
            seq_exp = seq_exp[seq_exp['Fragment ion charge'] == a]
            for u in theo_fragment:
                    theo_frag_M = (u * a) - (h_mass * a)
                    seq_exp['fragment M'] = (seq_exp['Fragment m/z'] * a) - (h_mass * a)
                    seq_exp['theoretical fragment'] = u
                    seq_exp['theoretical fragment M'] = theo_frag_M
                    seq_exp['Da error'] = abs(seq_exp['fragment M'] - theo_frag_M)
                    seq_exp_filter = seq_exp.sort_values(by='Da error')
                    seq_exp_filter = seq_exp_filter[seq_exp_filter['Da error'] <= error_fragment]
                    if len(seq_exp_filter) > 0:
                        sequence = seq_exp_filter['Possible Sequence'].values.tolist()
                        scan = seq_exp_filter['Scan #'].values.tolist()
                        theo_precursor = seq_exp_filter['Theoretical precursor'].values.tolist()
                        actual_precursor = seq_exp_filter['Actual precursor'].values.tolist()
                        precursor_error = seq_exp_filter['ppm error'].values.tolist()
                        precursor_theretical_charge = seq_exp_filter['Precursor theoretical charge'].values.tolist()
                        actual_precursor_charge = seq_exp_filter['Precursor actual charge'].values.tolist()
                        fragment_charge = seq_exp_filter['Fragment ion charge'].values.tolist()
                        fragment_mz = seq_exp_filter['Fragment m/z'].values.tolist()
                        intensity = seq_exp_filter['Fragment intensity'].values.tolist()
                        resolution = seq_exp_filter['Fragment resolution'].values.tolist()
                        frag_act_M = seq_exp_filter['fragment M'].values.tolist()
                        theo_frag = seq_exp_filter['theoretical fragment'].values.tolist()
                        theo_frag_M = seq_exp_filter['theoretical fragment M'].values.tolist()
                        Da_err = seq_exp_filter['Da error'].values.tolist()
                        
                        for v in sequence:
                            frag_possible_seqeuence_store.append(v)
                        for w in scan:
                            frag_scan_store.append(w)
                        for x in theo_precursor:
                            frag_theo_precursor_store.append(x)
                        for y in actual_precursor:
                            frag_act_precursor_store.append(y)
                        for z in precursor_error:
                            frag_precursor_err_store.append(z)
                        for aa in precursor_theretical_charge:
                            frag_precursor_theo_charge_store.append(aa)
                        for ab in actual_precursor_charge:
                            frag_precursor_act_charge_store.append(ab)
                        for ac in fragment_charge:
                            frag_ion_charge_store.append(ac)
                        for ad in fragment_mz:
                            frag_mz_store.append(ad)
                        for ae in intensity:
                            frag_intensity_store.append(ae)
                        for af in resolution:
                            frag_resoltuion_store.append(af)
                        for ag in frag_act_M:
                            frag_act_M_store.append(ag)
                        for ah in theo_frag:
                            frag_theo_mz_store.append(ah)
                        for ai in theo_frag_M:
                            frag_theo_M_store.append(ai)
                        for aj in Da_err:
                            Da_err_store.append(aj)
            
fragment_matches = pd.DataFrame()
fragment_matches['Fragment m/z']= frag_mz_store
fragment_matches['resolution']= frag_resoltuion_store
fragment_matches['charge']= frag_ion_charge_store
fragment_matches['intensity']= frag_intensity_store
fragment_matches['Precursor']= frag_act_precursor_store
fragment_matches['Scan #']= frag_scan_store
fragment_matches['Precursor_Charge']= frag_precursor_act_charge_store
fragment_matches['ppm error']= frag_precursor_err_store
fragment_matches['Possible sequence'] = frag_possible_seqeuence_store
fragment_matches['Theoretical precursor']= frag_theo_precursor_store
fragment_matches['Theoretical precursor charge']= frag_precursor_theo_charge_store
fragment_matches['fragment M']= frag_act_M_store
fragment_matches['theoretical fragment']= frag_theo_mz_store
fragment_matches['theoretical fragment M']= frag_theo_M_store
fragment_matches['Da error']= Da_err_store   

from user_input import non0_fragment_matches_out
file_out_path = non0_fragment_matches_out
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        fragment_matches.to_csv(filec,index=False)