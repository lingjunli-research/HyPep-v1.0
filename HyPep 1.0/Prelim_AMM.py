# -*- coding: utf-8 -*-
"""
Created on Mon May  2 09:34:42 2022

@author: lawashburn
"""

import pandas as pd
import csv
from datetime import datetime
now = datetime.now()

from datetime import datetime
now = datetime.now()
print('Preliminary accurate mass match in progress',datetime.now())

from datetime import datetime
now = datetime.now()
from user_input import db_path
from user_input import top_FD_combined_out_file_name_path
topFD_combined_out = pd.read_csv(top_FD_combined_out_file_name_path)

from user_input import prelim_AMM_out_path
from user_input import precursor_error_cutoff
from user_input import precursor_charges
from user_input import prelim_AMM_out_file_name_path


final_directory = prelim_AMM_out_path
error_cutoff = precursor_error_cutoff
h_mass = 1.00784
charges = precursor_charges

db = pd.read_csv(db_path)

if len(db)<1: #throws an error if more than one ion file exists
    raise ValueError('Database file is empty')

for a in charges:
    db[str(a)] = ((db['Monoisotopic Mass']-h_mass)+(h_mass * a))/a

sequences = db['Sequence'].values.tolist()
sequence_store = []
charge_store = []
target_mz_store = []
query_mz_store = []
ppm_err_store = []

query = topFD_combined_out

for b in sequences:
    db2 = db[db['Sequence'] == b]
    for c in charges:
        target_mz = db2[str(c)]
        query2 = query[query['Precursor actual charge'] == c]
        monomz = query2['Precursor actual m/z'].values.tolist()
        for d in target_mz:
            for e in monomz:
                ppm_err = ((abs(d-e))/d) * 1E6
                
                sequence_store.append(b)
                charge_store.append(c)
                target_mz_store.append(d)
                query_mz_store.append(e)
                ppm_err_store.append(ppm_err)
                
results_df = pd.DataFrame()
results_df['Sequence'] = sequence_store
results_df['Precursor actual charge'] = charge_store  
results_df['Theoretical precursor m/z'] = target_mz_store
results_df['Precursor actual m/z'] = query_mz_store
results_df['Precursor error (ppm)'] = ppm_err_store
results_df = results_df[results_df['Precursor error (ppm)'] <= error_cutoff]
results_df = results_df.sort_values(by='Precursor error (ppm)')
results_df = results_df.drop_duplicates(subset=['Sequence','Precursor actual charge'],keep='first').reset_index(drop=True)

complete_results_df = pd.merge(results_df, query,  how='left', left_on=['Precursor actual m/z','Precursor actual charge'], right_on = ['Precursor actual m/z','Precursor actual charge'])

complete_results_df = complete_results_df.sort_values(by='Precursor error (ppm)')
complete_results_df = complete_results_df.drop_duplicates(subset=['Sequence','Precursor actual charge'],keep='first').reset_index(drop=True)

complete_results_df_accession = complete_results_df.merge(db,left_on='Sequence',right_on='Sequence')
complete_results_df_accession  = complete_results_df_accession.iloc[: , :6]

file_path = prelim_AMM_out_file_name_path
with open(file_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        complete_results_df_accession.to_csv(filec,index=False)  
        


prelim_AMM_out = complete_results_df
