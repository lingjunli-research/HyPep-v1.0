# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 12:32:38 2022

@author: lawashburn
"""

import pandas as pd
import csv

no_ptm_list = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\09_PSM_assign_looping\SHS_priority_list.csv"
assecion_noPTM = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\09_PSM_assign_looping\NP_database_assession.csv"
assession_PTM = r"C:\Users\lawashburn\Desktop\ALC50_Mass_Search_Files\Crustacean_DB_Masses_sequences_updated20220407.csv"
final_directory = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\09_PSM_assign_looping"

no_ptm_list = pd.read_csv(no_ptm_list)
assecion_noPTM = pd.read_csv(assecion_noPTM)
assession_PTM = pd.read_csv(assession_PTM)

merge1 = assecion_noPTM.merge(no_ptm_list, left_on='Sequence', right_on='Sequence (PTMs removed)')
merge2 = merge1.merge(assession_PTM, left_on='Accession Number', right_on='accession')

priority_list = merge2['sequence']

file_name = 'priority_list.csv'
file_out_path = final_directory + '\\' + file_name
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        priority_list.to_csv(filec,index=False)