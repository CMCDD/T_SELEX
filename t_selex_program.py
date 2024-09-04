#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import os
import argparse
import T_SELEX_program
from T_SELEX_program import gen_aptamers, fold_and_composition, Mol_docking_calc, SSC, aptamerbase, tertiary_structure

def run_aptamer_analysis(gen_aptamers_seed=1, gen_aptamers_length='randomize', gen_aptamers_number=100, 
                         receptor_file=None, receptor_name=None, ligands_directory=None, software_path=None, 
                         output_csv="output.csv", ssc_csv=None, aptamer_type=None):
    # Print introductory information
    print('///////////////////////////////////////////////////////////////////////////////////////////////')
    print('///////////////////////////////////////////////////////////////////////////////////////////////')
    print('/                                                                                             /')
    print('/ This program is written by Kabelo under the supervision of prof K .Lobb and Dr T. Tshiwawa  /')
    print('/                                                                                             /')
    print('/ It is designed for aptamer analysis, including generation, folding, and virual screening.   /')
    print('/                                                                                             /')
    print('/           For questions or support, please contact Kabelo at racerpoet1@gmail.com           /')
    print('/                                                                                             /')
    print('///////////////////////////////////////////////////////////////////////////////////////////////')
    print('///////////////////////////////////////////////////////////////////////////////////////////////')
    

    
    print("\n.\\\\\\\\\\\\\\Setting Parameters:\\\\\\\\\\\\\ ")
    print(f'gen_aptamers_seed: {gen_aptamers_seed}')
    print(f'gen_aptamers_length: {gen_aptamers_length}')
    print(f'gen_aptamers_number: {gen_aptamers_number}')
    print(f'receptor_file: {receptor_file}')
    print(f'receptor_name: {receptor_name}')
    print(f'ligands_directory: {ligands_directory}')
    print(f'software_path: {software_path}')
    print(f'output_csv: {output_csv}')
    print(f'ssc_csv: {ssc_csv}')
    print(f'aptamer_type: {aptamer_type}')

    # Convert gen_aptamers_length to integer if it's not 'randomize'
    if gen_aptamers_length != 'randomize':
        try:
            gen_aptamers_length = int(gen_aptamers_length)
        except ValueError:
            raise ValueError("gen_aptamers_length must be 'randomize' or an integer.")
    
    # Generate aptamers
    p = gen_aptamers(gen_aptamers_seed, gen_aptamers_length, gen_aptamers_number)
    print("Generated Aptamers:")
    print(p)
    
    # Fold and compose aptamers
    a = fold_and_composition(p)
    print("Folded and Composed Aptamers:", a)
    
    # Generate tertiary structure
    k = tertiary_structure(a['Aptamer'], a['MFE structure'])
    
    # If receptor_file is provided, ensure required arguments are provided
    if receptor_file:
        if not all([ligands_directory, software_path, receptor_name]):
            raise ValueError("If --receptor_file is provided, then --ligands_directory, --software_path, and --receptor_name are also required.")
        
        
        home_directory = os.path.expanduser('~')
        directory_path = os.path.join(home_directory, software_path)
        
        os.chdir(home_directory)
        
        # Load DataFrame
        #df1 = pd.read_csv(output_csv)
        
        # Perform molecular docking calculation
        p = Mol_docking_calc(data_frame=a, MFE_column='Minimum free Energy', receptor_name=receptor_name,
                             receptor=receptor_file, ligands_directory=ligands_directory, 
                             directory_path=directory_path, Ap_folded=True)
        print("Aptamer:", a['Aptamer'], "MFE Structure:", a['MFE structure'])
    
        # Perform SSC calculation if SSC CSV is provided
        if ssc_csv:
            r = SSC(a, ssc_csv, 10)
            print(r)
    
    # Get aptamer base if argument provided
    if aptamer_type:
        l = aptamerbase(aptamer_type)
        print("Aptamer Base:", l)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run aptamer analysis.")
    
    parser.add_argument('--gen_aptamers_seed', type=int, default=1, help="Seed for generating aptamers.")
    parser.add_argument('--gen_aptamers_length', type=str, default='randomize', help="Length of aptamers to generate. Use 'randomize' for random length or an integer for a fixed length.")
    parser.add_argument('--gen_aptamers_number', type=int, default=100, help="Number of aptamers to generate.")
    parser.add_argument('--receptor_file', type=str, help="Path to receptor file. If provided, other arguments become required.")
    parser.add_argument('--ligands_directory', type=str, help="Path to ligands directory.")
    parser.add_argument('--software_path', type=str, help="Path to the software directory.")
    parser.add_argument('--receptor_name', type=str, help="The name of the target macromolecule to help with storage of the results.")
    parser.add_argument('--output_csv', type=str, default="output.csv", help="Path to the updated aptamers CSV file.")
    parser.add_argument('--ssc_csv', type=str, help="Path to the SSC CSV file. Optional.")
    parser.add_argument('--aptamer_type', type=str, help="Type of aptamers (e.g., RNA). Optional.")
    
    args = parser.parse_args()
    
    run_aptamer_analysis(args.gen_aptamers_seed, args.gen_aptamers_length, args.gen_aptamers_number,
                         args.receptor_file, args.receptor_name, args.ligands_directory, args.software_path,
                         args.output_csv, args.ssc_csv, args.aptamer_type)
