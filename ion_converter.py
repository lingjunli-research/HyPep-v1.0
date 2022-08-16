# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 14:33:40 2022

@author: lawashburn
"""
import csv
import pandas as pd
from amazon.ion import simpleion

ion_file_path = open(r"C:\Users\lawashburn\Documents\HyPep1.0\Theoretical_Fragment_Calculator\ADDMTEEAALQAAED.ions","rb")
data = ion_file_path.read()
iondata = simpleion.loads(data, single_value=False)   # Loading as ion data
print(iondata)