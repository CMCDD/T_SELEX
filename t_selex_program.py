#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 15:39:36 2024

@author: s1800206
"""

import pandas as pd
import os
import argparse
from T_SELEX_program import gen_aptamers, fold_and_composition, tertiary_structure, Mol_docking_calc, SSC, loops, mass, aptamerbase

def run_aptamer_analysis(gen_aptamers_seed, gen_aptamers_length, gen_aptamers_width, 
                         receptor_file, ligands_directory, software_path, 
                         output_csv, ssc_csv, aptamer_type):
    # Generate aptamers
    p = gen_aptamers(gen_aptamers_seed, gen_aptamers_length, gen_aptamers_width)
    print("Generated Aptamers:", p)
    
    # Fold and compose aptamers
    a = fold_and_composition(p)
    print("Folded and Composed Aptamers:", a)
    
    # Define paths
    home_directory = os.path.expanduser('~')
    directory_path = os.path.join(home_directory, software_path)
    
   
    os.chdir(home_directory)
    
    # Load DataFrame
    df1 = pd.read_csv(output_csv)
    
    # Perform molecular docking calculation
    p = Mol_docking_calc(data_frame=a, MFE_column='Minimum free Energy', receptor_name="6zsl",
                         receptor=receptor_file, ligands_directory=ligands_directory, 
                         directory_path=directory_path, Ap_folded=True)
    print("Aptamer:", a['Aptamer'], "MFE Structure:", a['MFE structure'])
    
    # Perform SSC calculation
    r = SSC(a, ssc_csv, 10)
    print(r) 
    
    # Get aptamer base
    l = aptamerbase(aptamer_type)
    print("Aptamer Base:", l)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run aptamer analysis.")
    
    parser.add_argument('--gen_aptamers_seed', type=int, required=True, help="Seed for generating aptamers.")
    parser.add_argument('--gen_aptamers_length', type=int, required=True, help="Length of aptamers to generate.")
    parser.add_argument('--gen_aptamers_width', type=int, required=True, help="Width of aptamers to generate.")
    parser.add_argument('--receptor_file', type=str, required=True, help="Path to receptor file.")
    parser.add_argument('--ligands_directory', type=str, required=True, help="Path to ligands directory.")
    parser.add_argument('--software_path', type=str, required=True, help="Path to the software directory.")
    parser.add_argument('--receptor_name', type=str, required=True, help="The name of the target macromolecule to help with storage of the results.") #avoid manual work of storing the results
    parser.add_argument('--output_csv', type=str, required=True, help="Path to the updated aptamers CSV file.")
    parser.add_argument('--ssc_csv', type=str, required=True, help="Path to the SSC CSV file.")
    parser.add_argument('--aptamer_type', type=str, required=True, help="Type of aptamers (e.g., RNA).")
    
    args = parser.parse_args()
    
    run_aptamer_analysis(args.gen_aptamers_seed, args.gen_aptamers_length, args.gen_aptamers_width,
                         args.receptor_file, args.ligands_directory, args.software_path,
                         args.output_csv, args.ssc_csv, args.aptamer_type)
