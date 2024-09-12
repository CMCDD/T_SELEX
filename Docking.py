#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 17:23:36 2023

@author: Kabelo Mokgopa
"""

import sys
import os 
import fnmatch
import csv
import pandas as pd
import tempfile





def sep(data_frame,target_column):
    
    nam = "aptamerd"
    
    folded =[]
    
    Unfolded =[]
    

    for a,b in zip(data_frame['Number'],data_frame[target_column]):
        
        if b ==0:
            
           lig_name = f'{nam}{a}.pdb'
           
           #print(lig_name)
           
           Unfolded.append(lig_name)
           
        if b !=0:
            
            lig_name= f'{nam}{a}.pdb'
            
            #print(lig_name)
            folded.append(lig_name)
            
        else:
            
            pass
        
    return folded,Unfolded







def aptamer_separator(data_frame,target_column):
    
    
    nam = "aptamerd"
    
    folded_best = []
    
    unfolded_best = []
    
    for a,b,c in zip(data_frame['Number'],data_frame['mfe_energiesRNAfold'],data_frame[target_column]):
        
        if b ==0 and c<-8:
            
            lig_name= f'{nam}{a}.pdb'
            
            
            
            unfolded_best.append(lig_name)
            
        elif b !=0 and c<-8:
            
            lig_name= f'{nam}{a}.pdb'
            folded_best.append(lig_name)
            
        else:
            
            pass
        
    return unfolded_best, folded_best





def Mol_docking_calc(data_frame,MFE_column,receptor_name,receptor,ligands_directory,directory_path,Ap_folded):
    my_files = os.listdir('/home/s1800206/hdock_test')

    #for files in my_files:
        #print(files)

    for file in my_files:
        if fnmatch.fnmatch(file,"ligand*"):
            pass
            #print(file)

    stream = os.popen('ls -la')
    output = stream.readlines()
    #print(output)



    # directory_path = '/home/s1800206/software/HDOCKlite-v1.1'
    # os.chdir(directory_path)
    # ligands = os.listdir("/home/s1800206/hdock_test/ligands")
    # receptor = "/home/s1800206/hdock_test/targets/target.pdb"

    # for ligand in ligands:
    #     if ligand.endswith(".pdb"):

    #         os.system('./hdock ' + ligand +' ' + receptor + ' -out Hdock.out')
    #         os.system('./createpl Hdock.out top100.pdb -nmax 100 -complex -models')


    directory_path = '/home/s1800206/software/HDOCKlite-v1.1'

    os.chdir(directory_path)

    # ligands_directory = '/home/s1800206/hdock_test/ligands'

    #receptor = '/home/s1800206/hdock_test/targets/hsa_miR_10b_5p.pdb'

    #receptor_name = "hsa_miR_10b_5p"

    # n= list(range(1, len(ligands_directory)))

    RSMD = []

    energy = []
    Ligand = []
    data = []
    #df = pd.read_csv("Updatedaptamers_with_energy.csv")




    folded ,unfolded = sep(data_frame,MFE_column)

    #print(unfolded_best)

    #ligands_directory = '/home/s1800206/hdock_test/ligands/aptamers'

    ligands = os.listdir(ligands_directory)

    lig=""
    if Ap_folded==True:
        
        for name1 in folded:

            for ligand in ligands:

                #print(ligand)

                if name1 == ligand:

                    #print(ligand)

                    if ligand.endswith('.pdb'):

                        name_ligand = ligand[:-4]

                        ligand_path = os.path.join(ligands_directory, ligand)


                        with tempfile.TemporaryDirectory() as tmpdirname:

                            os.system(f'cp {ligand_path} {tmpdirname}')

                            ligandfd = os.listdir(tmpdirname)

                            ligand_for_docking = os.path.join(tmpdirname,ligandfd[0])

                            command1 = './hdock {} {} -out Hdock.out'.format( receptor,ligand_for_docking)

                            print("RUNNING ",command1)

                            os.system(command1)

                            command2 = './createpl Hdock.out top100.pdb -nmax 100 -complex -models'

                            print("RUNNING ",command2)

                            os.system(command2)

                            parent_folder_name = "Docking"

                            folder_name = "receptor_name"

                            folded_aptamers = "folded"

                            #dirname = 'outputs_{}_{}'.format(name_ligand,receptor_name)
                            #if not os.path.exists(dirname):
                                #os.mkdir(dirname)
                            #out_path = os.path.join(directory_path,folder_name,dirname)

                            parent_folder_name = "Docking"

                            folder_name = "receptor_name"

                            folded_aptamers = "folded"

                            home_directory = os.path.expanduser( '~' )

                            folder_name = "receptor_name"

                            # Creating path for pasting the the models to overcome overwrite

                            path1 = os.path.join(home_directory,parent_folder_name)

                            if not os.path.exists(path1):

                                os.mkdir(path1)

                            path2 = os.path.join(home_directory,parent_folder_name,receptor_name)

                            if not os.path.exists(path2):

                                os.mkdir(path2)

                            path3 = os.path.join(home_directory,parent_folder_name,receptor_name, folded_aptamers)

                            if not os.path.exists(path3):

                                os.mkdir(path3)

                            outcome_path_for_docked_models = os.path.join(home_directory,parent_folder_name,receptor_name,folded_aptamers,name_ligand)

                            if not os.path.exists(outcome_path_for_docked_models):

                                os.mkdir(outcome_path_for_docked_models)


                            my_files = os.listdir('/home/s1800206/software/HDOCKlite-v1.1')

                            lig = name_ligand

                            Ligand.append(name_ligand)

                            e = None

                            rsmd = None

                            for file in my_files:

                                if fnmatch.fnmatch(file, "model*"):

                                    file_path = os.path.join(directory_path,file)

                                    cmd = 'cp {} {}'.format(file_path,outcome_path_for_docked_models)

                                    os.system(cmd)

                                if fnmatch.fnmatch(file,""):

                                    cmd2 = 'cp {} {}'.format(file_path,outcome_path_for_docked_models)

                                    os.system(cmd2)

                                if fnmatch.fnmatch(file,"Hdock.out"):

                                    file_path2 = os.path.join(directory_path,file)

                                    cmd3 = 'cp {} {}'.format(file_path2,outcome_path_for_docked_models)

                                    os.system(cmd3)

                            models_files = os.listdir(outcome_path_for_docked_models)

                            models_files

                            energy_models =[]

                            for model in models_files:

                                with open(model,"r") as f:

                                    lines = f.readlines()

                                    for line in lines:

                                        if "REMARK Score" in line:

                                            words =line.split()

                                            e = words[2]

                                            #print(f'{lig},{words[2]} kcal/mol')

                                            energy_models.append(float(words[2]))

                            FQ = round(sum(energy_models)/len(energy_models),2)
                            print(FQ)

                            ####
                    with open("model_1.pdb","r") as f:

                        lines = f.readlines()

                        for line in lines:

                            if "REMARK Score" in line:

                                words =line.split()

                                e = float(words[2])

                                print(f'{lig},{words[2]} kcal/mol')

                                energy.append(words[2])

                                # confidance score
                                import math

                                exponent = ((0.02) * (e + 150))

                                #print(exponent)

                                xpnt= math.exp(exponent)

                                #print(xpnt)

                                dinominator = 1+xpnt

                                con_score = 1/dinominator

                                cs = round(con_score *100,3)

                                print(cs)

                            if "REMARK RMSD" in line:

                                words = line.split()

                                rsmd = words[2]

                                RSMD.append(words[2])

                    data.append([lig, e, rsmd,FQ,cs ])


        headers = ["ligand","Energy","RMSD", "Fitness quality","Confidence Score (%)"]

        with open("Docking_files.csv","w",newline="") as csvfile:

            writer =csv.writer(csvfile)

            writer.writerow(headers)

            for row in data:

                writer.writerow(row)

        df= pd.read_csv("Docking_files.csv")

        df.sort_values(["Energy"], axis =0, ascending=True,inplace =True)

        df.to_csv('Docking_files.csv', index=False)

        df= pd.read_csv("Docking_files.csv")
        #os.system('cp')
        print("####################################################")
        print("#                                                  #")
        print("#       Please reference T-SELEX program           #")
        print("#                                                  #")
        print("####################################################")      

        print(df)
    else:
        
        for name1 in unfolded:

            for ligand in ligands:

                #print(ligand)

                if name1 == ligand:

                    #print(ligand)

                    if ligand.endswith('.pdb'):

                        name_ligand = ligand[:-4]

                        ligand_path = os.path.join(ligands_directory, ligand)


                        with tempfile.TemporaryDirectory() as tmpdirname:

                            os.system(f'cp {ligand_path} {tmpdirname}')

                            ligandfd = os.listdir(tmpdirname)

                            ligand_for_docking = os.path.join(tmpdirname,ligandfd[0])

                            command1 = './hdock {} {} -out Hdock.out'.format( receptor,ligand_for_docking)

                            print("RUNNING ",command1)

                            os.system(command1)

                            command2 = './createpl Hdock.out top100.pdb -nmax 100 -complex -models'

                            print("RUNNING ",command2)

                            os.system(command2)

                            parent_folder_name = "Docking"

                            folder_name = "receptor_name"

                            folded_aptamers = "folded"

                            #dirname = 'outputs_{}_{}'.format(name_ligand,receptor_name)
                            #if not os.path.exists(dirname):
                                #os.mkdir(dirname)
                            #out_path = os.path.join(directory_path,folder_name,dirname)

                            parent_folder_name = "Docking"

                            folder_name = "receptor_name"

                            folded_aptamers = "folded"

                            home_directory = os.path.expanduser( '~' )

                            folder_name = "receptor_name"

                            # Creating path for pasting the the models to overcome overwrite

                            path1 = os.path.join(home_directory,parent_folder_name)

                            if not os.path.exists(path1):

                                os.mkdir(path1)

                            path2 = os.path.join(home_directory,parent_folder_name,receptor_name)

                            if not os.path.exists(path2):

                                os.mkdir(path2)

                            path3 = os.path.join(home_directory,parent_folder_name,receptor_name, folded_aptamers)

                            if not os.path.exists(path3):

                                os.mkdir(path3)

                            outcome_path_for_docked_models = os.path.join(home_directory,parent_folder_name,receptor_name,folded_aptamers,name_ligand)

                            if not os.path.exists(outcome_path_for_docked_models):

                                os.mkdir(outcome_path_for_docked_models)


                            my_files = os.listdir('/home/s1800206/software/HDOCKlite-v1.1')

                            lig = name_ligand

                            Ligand.append(name_ligand)

                            e = None

                            rsmd = None

                            for file in my_files:

                                if fnmatch.fnmatch(file, "model*"):

                                    file_path = os.path.join(directory_path,file)

                                    cmd = 'cp {} {}'.format(file_path,outcome_path_for_docked_models)

                                    os.system(cmd)

                                if fnmatch.fnmatch(file,""):

                                    cmd2 = 'cp {} {}'.format(file_path,outcome_path_for_docked_models)

                                    os.system(cmd2)

                                if fnmatch.fnmatch(file,"Hdock.out"):

                                    file_path2 = os.path.join(directory_path,file)

                                    cmd3 = 'cp {} {}'.format(file_path2,outcome_path_for_docked_models)

                                    os.system(cmd3)

                            models_files = os.listdir(outcome_path_for_docked_models)

                            models_files

                            energy_models =[]

                            for model in models_files:

                                with open(model,"r") as f:

                                    lines = f.readlines()

                                    for line in lines:

                                        if "REMARK Score" in line:

                                            words =line.split()

                                            e = words[2]

                                            #print(f'{lig},{words[2]} kcal/mol')

                                            energy_models.append(float(words[2]))

                            FQ = round(sum(energy_models)/len(energy_models),2)
                            print(FQ)

                            ####
                    with open("model_1.pdb","r") as f:

                        lines = f.readlines()

                        for line in lines:

                            if "REMARK Score" in line:

                                words =line.split()

                                e = float(words[2])

                                print(f'{lig},{words[2]} kcal/mol')

                                energy.append(words[2])

                                # confidance score
                                import math

                                exponent = ((0.02) * (e + 150))

                                #print(exponent)

                                xpnt= math.exp(exponent)

                                #print(xpnt)

                                dinominator = 1+xpnt

                                con_score = 1/dinominator

                                cs = round(con_score *100,3)

                                print(cs)

                            if "REMARK RMSD" in line:

                                words = line.split()

                                rsmd = words[2]

                                RSMD.append(words[2])

                    data.append([lig, e, rsmd,FQ,cs ])


        headers = ["ligand","Energy","RMSD", "Fitness quality","Confidence Score (%)"]

        with open("Docking_files.csv","w",newline="") as csvfile:

            writer =csv.writer(csvfile)

            writer.writerow(headers)

            for row in data:

                writer.writerow(row)

        df= pd.read_csv("Docking_files.csv")

        df.sort_values(["Energy"], axis =0, ascending=True,inplace =True)

        df.to_csv('Docking_files.csv', index=False)

        df= pd.read_csv("Docking_files.csv")
        #os.system('cp')
        print("####################################################")
        print("#                                                  #")
        print("#       Please reference T-SELEX program           #")
        print("#                                                  #")
        print("####################################################")      

        print(df)
    return df         
