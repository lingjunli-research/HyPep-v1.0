# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 14:05:16 2022

@author: lawashburn
"""

import pickle
import pandas as pd

#Simple file to pull all variables from GUI
#All values will be exported as strings


#with open('min_alc.pkl', 'rb') as file_a:      
#    # Call load method to deserialze
#    malc_entry = pickle.load(file_a)

#with open('discovery_alc.pkl', 'rb') as file_b:      
    # Call load method to deserialze
#    discovery_alc = pickle.load(file_b)
    
#with open('motif_path.pkl', 'rb') as file_c:      
    # Call load method to deserialze
#    motif_path = pickle.load(file_c)
    
#with open('discovery_fdr_score.pkl', 'rb') as file_d:      
#    # Call load method to deserialze
#    discovery_fdr_score = pickle.load(file_d)
    
with open('database_path.pkl', 'rb') as file_e:      
    # Call load method to deserialze
    db_path = pickle.load(file_e)
    
with open('fdr.pkl', 'rb') as file_f:      
    # Call load method to deserialze
    fdr_percent = pickle.load(file_f)

with open('max_precursor_z.pkl', 'rb') as file_g:      
    # Call load method to deserialze
    max_p_charge = pickle.load(file_g)

with open('max_fragment_z.pkl', 'rb') as file_h:      
    # Call load method to deserialze
    max_f_charge = pickle.load(file_h)
    
with open('precursor_error.pkl', 'rb') as file_g:      
    # Call load method to deserialze
    p_err = pickle.load(file_g)

with open('fragment_error.pkl', 'rb') as file_h:      
    # Call load method to deserialze
    f_err = pickle.load(file_h)
    
with open('promex_score.pkl', 'rb') as file_i:      
    # Call load method to deserialze
    pro_mex = pickle.load(file_i)
    
with open('loops.pkl', 'rb') as file_j:      
    # Call load method to deserialze
    loop = pickle.load(file_j)
    
with open('sample_name.pkl', 'rb') as file_k:      
    # Call load method to deserialze
    sample_ID = pickle.load(file_k)

with open('peaks_path.pkl', 'rb') as file_l:      
    # Call load method to deserialze
    peaks_path = pickle.load(file_l)

with open('top_fd_path.pkl', 'rb') as file_m:      
    # Call load method to deserialze
    topfd_path = pickle.load(file_m)

with open('RawConverter_path.pkl', 'rb') as file_n:      
    # Call load method to deserialze
    rawconverter_path = pickle.load(file_n)  

with open('theoretical_fragment_path.pkl', 'rb') as file_o:      
    # Call load method to deserialze
    pp_path = pickle.load(file_o)

with open('output_directory.pkl', 'rb') as file_p:      
    # Call load method to deserialze
    out_path = pickle.load(file_p)

with open('window_size.pkl', 'rb') as file_q:      
    # Call load method to deserialze
    window_size = pickle.load(file_q)  

with open('discovery_mode.pkl', 'rb') as file_r:      
    # Call load method to deserialze
    to_discover = pickle.load(file_r)

with open('fdr_algorithm.pkl', 'rb') as file_s:      
    # Call load method to deserialze
    fdr_selection = pickle.load(file_s)


#if len(motif_path)>0: #motif decision
#    motif_dec = 1
#else:
#    motif_dec = 0

db_ptm_mass_replace = pd.read_csv(db_path)
db_ptm_mass_replace = db_ptm_mass_replace.replace({'Amidated':0.9840,'Pyro-glu':-17.0265,'Sulfo':-17.0265,
                                                   'Oxidation':15.9949})    
