#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 14:40:34 2023

@author:Kabelo Mokgopa
"""

import os 

home_directory = os.path.expanduser( '~' )
home_directory

receptor_name = "hsa_miR_10b_5p"
os.chdir(f'{home_directory}/Docking/{receptor_name}/folded')
aptamers_directory = os.listdir(f'{home_directory}/Docking/{receptor_name}/folded')
aptamers_directory

for directory in aptamers_directory:
    present_directory = os.chdir(f'{home_directory}/Docking/{receptor_name}/folded/{directory}')
    current_path = f'{home_directory}/Docking/{receptor_name}/folded/{directory}'
    copy1 = 'cp {}/pymol_automated_model1_3r.pml {}'.format(home_directory, current_path)
    copy2 = 'cp {}/pymol_automated_model1_to_5_2r.pml {}'.format(home_directory, current_path)
    copy3 = 'cp {}/pymol_automated_model1_to_10_2r.pml {}'.format(home_directory, current_path)
    
    os.system(copy1)
    os.system(copy2)
    os.system(copy3)
    #os.system("pwd")
    os.system('pymol pymol_automated_model1_3r.pml')
    os.system('pymol pymol_automated_model1_to_5_2r.pml')
    os.system('pymol pymol_automated_model1_to_10_2r.pml')
    
