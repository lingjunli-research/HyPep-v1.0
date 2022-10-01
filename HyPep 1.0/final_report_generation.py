# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 07:46:44 2022

@author: lawashburn
"""


import pandas as pd
import csv

from user_input import results_out_path
from user_input import SHS_out_path
from user_input import SHS_results_path
from user_input import AMM_out_path
from user_input import psm_matches_path
from user_input import db_path

file_name =  SHS_out_path
file_out_path = SHS_out_path
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        SHS_results_path.to_csv(filec,index=False)


amm_results = pd.read_csv(psm_matches_path)
amm_results = amm_results.drop(columns=['scan count'],axis=1)

file_name =  AMM_out_path
file_out_path = AMM_out_path
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        amm_results.to_csv(filec,index=False)

database = pd.read_csv(db_path)

amm_w_accession = amm_results.merge(database,left_on='sequence',right_on='Sequence',how='left')

amm_w_accession['Accession Number'] = (amm_w_accession['Accession Number']).astype(int)
SHS_results_path['Accession Number'] = (SHS_results_path['Accession Number']).astype(int)

merged_results = amm_w_accession.merge(SHS_results_path,left_on='Accession Number',right_on='Accession Number',how='left')

formatted_table = pd.DataFrame()
formatted_table['Sample Name'] = merged_results['sample']
formatted_table['Sequence'] = merged_results['sequence']
formatted_table['AMM Coverage'] = merged_results['coverage']
formatted_table['SHS Score'] = merged_results['Score']
formatted_table['Scan # (AMM)'] = merged_results['scan']
formatted_table['% ALC with Scan # (SHS)'] = merged_results['Scan_Numbers[ALC%]']

file_name =  results_out_path
file_out_path = results_out_path
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        formatted_table.to_csv(filec,index=False)
        
#Exporting actual FDR
from driver_identification_mode import closest_to_user_FDR
from hypep_import import out_path
from hypep_import import fdr_percent

user_fdr = 'User-defined %FDR Threshold' + str(fdr_percent) + '%'
actual_fdr = 'Actual %FDR' + str(closest_to_user_FDR) + '%'

with open(out_path+'\\Actual %FDR.txt','w') as filec:
        filec.write(user_fdr)
        filec.write('\n')
        filec.write(actual_fdr)
