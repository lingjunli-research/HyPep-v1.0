# -*- coding: utf-8 -*-
"""
Created on Mon May  2 08:53:41 2022

@author: lawashburn
"""
#This script is fully functional (6/30/22)
import csv
import pandas as pd
from user_input import topfd_path
topFD_path = topfd_path
from user_input import rawconverter_path
from user_input import topFD_rawconverter_combined_out_path
final_directory = topFD_rawconverter_combined_out_path
from user_input import promex_cutoff
from user_input import sample_name
from RawConverter_Formattingv2 import SIM_result
from user_input import rawconverter_formatted_out_file_name_path_txt

topFD = pd.read_csv(topFD_path)
raw_converter = pd.read_csv(rawconverter_formatted_out_file_name_path_txt, sep=" ",skiprows=[0], 
                            names= ["m/z", "resolution", "charge", "intensity","MS2",'scan_number','precursor_charge','empty'])

topFD = topFD[topFD['promex_score'] >= promex_cutoff]

precursor_mz = []
precursor_z = []

topFD_mz = topFD['MonoMz'].values.tolist()
topFD_z = topFD['Charge'].values.tolist()

rawconv_mz = raw_converter['MS2'].values.tolist()
rawconv_z = raw_converter['precursor_charge'].values.tolist()

for a in topFD_mz:
    precursor_mz.append(a)

for b in rawconv_mz:
    precursor_mz.append(b)
    
for c in topFD_z:
    precursor_z.append(c)

for d in rawconv_z:
    precursor_z.append(d)

precursor = pd.DataFrame()
precursor['Precursor actual m/z'] = precursor_mz
precursor['Precursor actual charge'] = precursor_z

precursor = precursor.drop_duplicates()

from user_input import top_FD_combined_out_file_name_path
file_out_path = top_FD_combined_out_file_name_path
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        precursor.to_csv(filec,index=False)
                