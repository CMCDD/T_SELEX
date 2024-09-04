#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 15:25:39 2023

@author: Kabelo Mokgopa
"""
import os
import pandas as pd

import Docking
from Docking import Mol_docking_calc

home_directory = os.path.expanduser( '~' )
sofware_path = 'software/HDOCKlite-v1.1'
directory_path = os.path.join(home_directory, sofware_path)
directory_path

os.chdir(home_directory)
df1 = pd.read_csv("Updatedaptamers_with_energy.csv")
r = '/home/s1800206/hdock_test/targets/hsa_miR_10b_5p.pdb'
l = '/home/s1800206/hdock_test/ligands/aptamers'
s = '/home/s1800206/software/HDOCKlite-v1.1'


p = Mol_docking_calc(data_frame= df1,MFE_column ='mfe_energiesRNAfold',receptor_name= "hsa_miR_10b_5p",receptor=r,ligands_directory=l,directory_path= s,Ap_folded=True)
print(p)