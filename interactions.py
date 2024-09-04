#!/usr/bin/env python3


def intarna(csv_file, aptamers_colname,target_seq,target_name,analysis):
    print(f"##################################################")
    print(f"#                                                #")
    print(f"#       Setting up Mutiscale RNA-RNA             #")
    print(f"#        interactions calculations              #")
    print(f'#        Please reference T-SELEX                #')
    print(f'#                                                #')
    print(f"##################################################")
    
    import numpy as np
    
    from fastcluster import linkage
    
    from scipy.spatial.distance import pdist
    
    from scipy.cluster.hierarchy import dendrogram
    
    import matplotlib.pyplot as plt
    
    import os
    
    import pandas as pd
    
    import matplotlib.pyplot as plt
    
    from collections import Counter
    
    from Bio import pairwise2
    
    from Bio.pairwise2 import format_alignment
    
    import RNA


    df = pd.read_csv(csv_file)
    
    aptamers = df[aptamers_colname].tolist()
    
    target_sequence = target_seq
    
    target_name = target_name  #MIMAT0003164
    
    set_built = "set -o noclobber"
    
    os.system(set_built)
    
    best_score = []
    
    folded_best_score = []
    
    unfolded_best_score = []
    
    good_score =[]
    
    folded_good_score = []
    
    unfolded_good_score = []
    
    better_score=[]
    
    folded_better_score = []
    
    unfolded_better_score = []
    
    folded_poor_score = []
    
    unfolded_poor_score = []
    
    both_best_and_good = []
    
    best_good_and_better = []
    
    No_interaction=[]
    
    folded_interactions = []
    
    Unfolded_interactions = []
    
    folded_interaction_energies = []
    
    unfolded_interaction_energies = []
    
    all_interacting_aptamers = []
    
    mfe_energiesRNAfold = []
    
    mfe_structureRNAfold = []
    
    folded_array = []
    
    unfolded_array =[]
    
    if analysis == False:

        for aptamer in aptamers:
            
            (mfe_structure,mfe_energy) = RNA.fold(aptamer)
            
            mfe_energiesRNAfold.append(mfe_energy)
            
            mfe_structureRNAfold.append(mfe_structure)

        df["MFE"] = mfe_energiesRNAfold
        
        df["MFE secondary structure"] =mfe_structureRNAfold
        
        df.to_csv("Updatedaptamers.csv")

        df1 = pd.read_csv("Updatedaptamers.csv")
        
        unfolded_aptamers = []
        
        folded_aptamers = []
        
        for k,i,j in zip(df1['Number'],df1['MFE'],df1['Aptamer']):
            
            if i==0:
                
                unfolded_aptamers.append([k,j])
                
            else:
                
                folded_aptamers.append([k,j])

        print(f'There are {len(folded_aptamers)} folded aptamers and {len(unfolded_aptamers)} Unfolded aptamers from the library. ')


        print("####################################################################")
        print("#                                                                  #")
        print("#                The results for folded Aptamers                   #")
        print("#                                                                  #")
        print("####################################################################") 


        for i in folded_aptamers:
            
            dirname = f"{target_name}_outputs"
            
            sub_dirname = "folded_outputs"
            
            if not os.path.exists(dirname):
                
                os.mkdir(dirname)
                
                print(f"Creating a directory named {target_name}_outputs....")

            if not os.path.exists(os.path.join(dirname, sub_dirname)):
                
                os.mkdir(os.path.join(dirname, sub_dirname))
                
                print(f"Creating a subdirectory named {sub_dirname}....")

            filename = f"intoutput_aptamer_{i[0]}_and_{target_name}.log"
            
            filepath = os.path.join(dirname, sub_dirname, filename)

            with open(filepath, "w") as f:
                
                f.write(f"################################################\n")
                
                f.write(f"#  these results are for folded aptamer {i[0]} #\n")
                
                f.write(f"#     against the target named {target_name}   #\n")
                
                f.write(f"#     Please remember to reference Intarna     #\n")
                
                f.write(f"#     and T-selex when writing report          #\n")
                
                f.write(f"################################################\n\n")

                cmd = f"IntaRNA -t {target_sequence} -q {i[1]} --outMode=D --seedBP=4"
                
                result = os.popen(cmd).read()

                f.write(result)

                f.write(f"################################################\n")
                
                f.write(f"#    the results below shows the pairwise      #\n")
                
                f.write(f"#     for the aptamer against the target       #\n")
                
                f.write(f"################################################\n\n")

                alignments = pairwise2.align.globalxx(target_sequence,aptamer)
                
                #test_alignments = pairwise2.align.localds( miR_21a,aptamer, blosum62 , -10, -1)
                
                for alignment in alignments:
                    
                    f.write(str(alignment))
                    
                    f.write(format_alignment(*alignment))


            f.close()
            
            with open(filepath,"r") as fin:
                
                lines=fin.readlines()
                
                for line in lines:
                    
                    if "interaction energy" in line:
                        
                        words=line.split()
                        
                        all_interacting_aptamers.append(aptamer)
                        
                        folded_interaction_energies.append(float(words[3]))
                        
                        if float(words[3]) < -8:
                            
                            folded_best_score.append(str(i[1]))
                            
                            best_score.append(str(i[1]))
                            
                            for base in i[1]:
                                
                                if base == 'U':
                                    
                                    print('\033[41m' + base, end="")
                                    
                                elif base == 'A':
                                    
                                    print('\033[44m' + base, end="")
                                    
                                elif base == 'G':
                                    
                                    print('\033[43m' + base, end="")
                                    
                                elif base == 'C':
                                    
                                    print('\033[42m' + base, end="")
                                    
                                else:
                                    
                                    print(base, end="")
                                    
                            print('\033[0m', words[2], words[3], words[4])
                            
                            folded_array.append([i[1], words[3]])

                            break

                        if -8 < float(words[3]) < -6:
                            
                            folded_good_score.append(str(i[1]))
                            
                            good_score.append(str(i[1]))
                            
                            for base in i[1]:
                                
                                if base == 'U':
                                    
                                    print('\033[41m' + base, end="")
                                    
                                elif base == 'A':
                                    
                                    print('\033[44m' + base, end="")
                                    
                                elif base == 'G':
                                    
                                    print('\033[43m' + base, end="")
                                    
                                elif base == 'C':
                                    
                                    print('\033[42m' + base, end="")
                                    
                                else:
                                    
                                    print(base, end="")
                                    
                            print('\033[0m',words[2], words[3],words[4])
                            
                            folded_array.append([i[1], words[3]])
                            
                            #print(aptamer, words[2], words[3])
                            break

                        if -6 < float(words[3]) < -4:  
                            
                            folded_better_score.append(str(i[1]))
                            
                            better_score.append(i[1])
                            
                            for base in i[1]:
                                
                                if base == 'U':
                                    
                                    print('\033[41m' + base, end="")
                                    
                                elif base == 'A':
                                    
                                    print('\033[44m' + base, end="")
                                    
                                elif base == 'G':
                                    
                                    print('\033[43m' + base, end="")
                                    
                                elif base == 'C':
                                    
                                    print('\033[42m' + base, end="")
                                    
                                else:
                                    
                                    print(base, end="")
                                    
                            print('\033[0m',words[2], words[3], words[4])
                            
                            folded_array.append([i[1], words[3]])
                            
                            #print(aptamer, words[2], words[3])
                            break

                        if -4 < float(words[3]) < -2:
                            
                            folded_poor_score.append(str(i[1]))
                            
                            for base in i[1]:
                                
                                if base == 'U':
                                    
                                    print('\033[41m' + base, end="")
                                    
                                elif base == 'A':
                                    
                                    print('\033[44m' + base, end="")
                                    
                                elif base == 'G':
                                    
                                    print('\033[43m' + base, end="")
                                    
                                elif base == 'C':
                                    
                                    print('\033[42m' + base, end="")
                                    
                                else:
                                    
                                    print(base, end="")
                                    
                            print('\033[0m',words[2], words[3],words[4])
                            
                            folded_array.append([i[1], words[3]])
                            
                            #print(aptamer,words[2], words[3])

                            break
                            
                        if float(words[3]) > -2:
                            
                            folded_array.append([i[1], words[3]])
                            
                    elif "no favorable interaction":
                        
                        No_interaction.append(aptamer)

                    else:
                        pass


        print("####################################################################")
        
        print("#                                                                  #")
        
        print("#              The results for unfolded Aptamers                   #")
        
        print("#                                                                  #")
        
        print("####################################################################") 

        #predicting the interactions for unfolded aptamres in the library
        
        for i in unfolded_aptamers:
            
            dirname = f"{target_name}_outputs"
            
            sub_dirname = "Unfolded_outputs"
            

            if not os.path.exists(dirname):
                
                os.mkdir(dirname)
                
                print(f"Creating a directory named {target_name}_outputs....")

            if not os.path.exists(os.path.join(dirname, sub_dirname)):
                
                os.mkdir(os.path.join(dirname, sub_dirname))
                
                print(f"Creating a subdirectory named {sub_dirname}....")
                

            filename = f"intoutput_aptamer_{i[0]}_and_{target_name}.log"
            
            filepath = os.path.join(dirname, sub_dirname, filename)
            

            with open(filepath, "w") as f:
                
                f.write(f"################################################\n")
                
                f.write(f"#    these results are for folded aptamer {i[0]} #\n")
                
                f.write(f"#     against the target named {target_name}   #\n")
                
                f.write(f"#     Please remember to reference Intarna     #\n")
                
                f.write(f"#     and T-selex when writing report          #\n")
                
                f.write(f"################################################\n\n")

                cmd = f"IntaRNA -t {target_sequence} -q {i[1]} --outMode=D --seedBP=4"
                
                result = os.popen(cmd).read()

                f.write(result)

                f.write(f"################################################\n")
                
                f.write(f"#    The results below shows the pairwise      #\n")
                
                f.write(f"#     for the aptamer against the target       #\n")
                
                f.write(f"################################################\n\n")

                alignments = pairwise2.align.globalxx(target_sequence,i[1])
                
                #test_alignments = pairwise2.align.localds( miR_21a,aptamer, blosum62 , -10, -1)
                
                for alignment in alignments:
                    
                    f.write(str(alignment))
                    
                    f.write(format_alignment(*alignment))
                    


            f.close()



            with open(filepath,"r") as fin:
                
                lines=fin.readlines()
                
                for line in lines:
                    
                    if "interaction energy" in line:
                        
                        words=line.split()
                        
                        all_interacting_aptamers.append(aptamer)
                        
                        unfolded_interaction_energies.append(float(words[3]))
                        
                        if float(words[3]) < -8:
                            
                            unfolded_best_score.append(str(i[1]))
                            
                            best_score.append(str(i[1]))
                            
                            for base in i[1]:
                                
                                if base == 'U':
                                    print('\033[41m' + base, end="")
                                    
                                elif base == 'A':
                                    
                                    print('\033[44m' + base, end="")
                                    
                                elif base == 'G':
                                    
                                    print('\033[43m' + base, end="")
                                    
                                elif base == 'C':
                                    
                                    print('\033[42m' + base, end="")
                                    
                                else:
                                    
                                    print(base, end="")
                                    
                            print('\033[0m', words[2], words[3], words[4])
                            
                            unfolded_array.append([i[1], words[3]])

                            break

                        if -8 < float(words[3]) < -6:
                            
                            unfolded_good_score.append(str(i[1]))
                            
                            good_score.append(str(i[1]))
                            
                            for base in i[1]:
                                
                                if base == 'U':
                                    print('\033[41m' + base, end="")
                                    
                                elif base == 'A':
                                    
                                    print('\033[44m' + base, end="")
                                    
                                elif base == 'G':
                                    
                                    print('\033[43m' + base, end="")
                                    
                                elif base == 'C':
                                    
                                    print('\033[42m' + base, end="")
                                    
                                else:
                                    
                                    print(base, end="")
                                    
                            print('\033[0m',words[2], words[3],words[4])
                            
                            unfolded_array.append([i[1], words[3]])
                            
                            #print(aptamer, words[2], words[3])
                            break

                        if -6 < float(words[3]) < -4: 
                            
                            unfolded_better_score.append(str(i[1]))
                            
                            better_score.append(str(i[1]))
                            
                            for base in i[1]:
                                
                                if base == 'U':
                                    
                                    print('\033[41m' + base, end="")
                                    
                                elif base == 'A':
                                    
                                    print('\033[44m' + base, end="")
                                    
                                elif base == 'G':
                                    
                                    print('\033[43m' + base, end="")
                                    
                                elif base == 'C':
                                    
                                    print('\033[42m' + base, end="")
                                    
                                else:
                                    
                                    print(base, end="")
                                    
                            print('\033[0m',words[2], words[3], words[4])
                            
                            unfolded_array.append([i[1], words[3]])
                            
                            #print(aptamer, words[2], words[3])
                            break

                        if -4 < float(words[3]) < -2:
                            
                            unfolded_poor_score.append(str(i[1]))
                            
                            for base in i[1]:
                                
                                if base == 'U':
                                    
                                    print('\033[41m' + base, end="")
                                    
                                elif base == 'A':
                                    
                                    print('\033[44m' + base, end="")
                                    
                                elif base == 'G':
                                    
                                    print('\033[43m' + base, end="")
                                    
                                elif base == 'C':
                                    
                                    print('\033[42m' + base, end="")
                                    
                                else:
                                    
                                    print(base, end="")
                                    
                            print('\033[0m',words[2], words[3],words[4])
                            
                            unfolded_array.append([i[1], words[3]])
                            
                            #print(aptamer,words[2], words[3])

                            break
                            
                        if float(words[3]) > -2:
                            
                            unfolded_array.append([i[1], words[3]])
                            
                    elif "no favorable interaction":
                        
                        No_interaction.append(aptamer)

                    else:
                        pass
                    
        
    else:
        
        
        for aptamer in aptamers:
            
            (mfe_structure,mfe_energy) = RNA.fold(aptamer)
            
            mfe_energiesRNAfold.append(mfe_energy)
            
            mfe_structureRNAfold.append(mfe_structure)

        df["MFE"] = mfe_energiesRNAfold
        
        df["MFE secondary structure"] =mfe_structureRNAfold
        
        df.to_csv("Updatedaptamers.csv")




        df1 = pd.read_csv("Updatedaptamers.csv")
        
        unfolded_aptamers = []
        
        folded_aptamers = []
        for k,i,j in zip(df1['Number'],df1['MFE'],df1['Aptamer']):
            
            if i==0:
                
                unfolded_aptamers.append([k,j])
                
            else:
                
                folded_aptamers.append([k,j])
                

        print(f'There are {len(folded_aptamers)} folded aptamers and {len(unfolded_aptamers)} Unfolded aptamers from the library. ')


        print("####################################################################")
        print("#                                                                  #")
        print("#                The results for folded Aptamers                   #")
        print("#                                                                  #")
        print("####################################################################") 

        for i in folded_aptamers:
            
            dirname = f"{target_name}_outputs"
            
            sub_dirname = "folded_outputs"

            if not os.path.exists(dirname):
                
                os.mkdir(dirname)
                
                print(f"Creating a directory named {target_name}_outputs....")
                

            if not os.path.exists(os.path.join(dirname, sub_dirname)):
                
                os.mkdir(os.path.join(dirname, sub_dirname))
                
                print(f"Creating a subdirectory named {sub_dirname}....")

            filename = f"intoutput_aptamer_{i[0]}_and_{target_name}.log"
            
            filepath = os.path.join(dirname, sub_dirname, filename)

            with open(filepath, "w") as f:
                
                f.write(f"################################################\n")
                
                f.write(f"#  these results are for folded aptamer {i[0]} #\n")
                
                f.write(f"#     against the target named {target_name}   #\n")
                
                f.write(f"#     Please remember to reference Intarna     #\n")
                
                f.write(f"#     and T-selex when writing report          #\n")
                
                f.write(f"################################################\n\n")

                cmd = f"IntaRNA -t {target_sequence} -q {i[1]} --outMode=D --seedBP=4"
                
                result = os.popen(cmd).read()

                f.write(result)

                f.write(f"################################################\n")
                
                f.write(f"#    the results below shows the pairwise      #\n")
                
                f.write(f"#     for the aptamer against the target       #\n")
                
                f.write(f"################################################\n\n")

                alignments = pairwise2.align.globalxx(target_sequence,aptamer)
                
                #test_alignments = pairwise2.align.localds( miR_21a,aptamer, blosum62 , -10, -1)
                
                for alignment in alignments:
                    
                    f.write(str(alignment))
                    
                    f.write(format_alignment(*alignment))


            f.close()



            with open(filepath,"r") as fin:
                
                lines=fin.readlines()
                
                for line in lines:
                    
                    if "interaction energy" in line:
                        
                        words=line.split()
                        
                        all_interacting_aptamers.append(aptamer)
                        
                        folded_interaction_energies.append(float(words[3]))
                        
                        if float(words[3]) < -8:
                            
                            folded_best_score.append(str(i[1]))
                            
                            best_score.append(str(i[1]))
                            
                            for base in i[1]:
                                
                                if base == 'U':
                                    
                                    print('\033[41m' + base, end="")
                                    
                                elif base == 'A':
                                    
                                    print('\033[44m' + base, end="")
                                    
                                elif base == 'G':
                                    
                                    print('\033[43m' + base, end="")
                                    
                                elif base == 'C':
                                    
                                    print('\033[42m' + base, end="")
                                    
                                else:
                                    
                                    print(base, end="")
                                    
                            print('\033[0m', words[2], words[3], words[4])
                            
                            folded_array.append([i[1], words[3]])

                            break

                        if -8 < float(words[3]) < -6:
                            
                            folded_good_score.append(str(i[1]))
                            
                            good_score.append(str(i[1]))
                            
                            for base in i[1]:
                                
                                if base == 'U':
                                    
                                    print('\033[41m' + base, end="")
                                    
                                elif base == 'A':
                                    
                                    print('\033[44m' + base, end="")
                                    
                                elif base == 'G':
                                    
                                    print('\033[43m' + base, end="")
                                    
                                elif base == 'C':
                                    
                                    print('\033[42m' + base, end="")
                                    
                                else:
                                    
                                    print(base, end="")
                                    
                            print('\033[0m',words[2], words[3],words[4])
                            
                            folded_array.append([i[1], words[3]])
                            
                            #print(aptamer, words[2], words[3])
                            break

                        if -6 < float(words[3]) < -4: 
                            
                            folded_better_score.append(str(i[1]))
                            
                            better_score.append(i[1])
                            
                            for base in i[1]:
                                
                                if base == 'U':
                                    print('\033[41m' + base, end="")
                                    
                                elif base == 'A':
                                    
                                    print('\033[44m' + base, end="")
                                    
                                elif base == 'G':
                                    
                                    print('\033[43m' + base, end="")
                                    
                                elif base == 'C':
                                    
                                    print('\033[42m' + base, end="")
                                    
                                else:
                                    
                                    print(base, end="")
                                    
                            print('\033[0m',words[2], words[3], words[4])
                            
                            folded_array.append([i[1], words[3]])
                            
                            #print(aptamer, words[2], words[3])
                            break

                        if -4 < float(words[3]) < -2:
                            
                            folded_poor_score.append(str(i[1]))
                            
                            for base in i[1]:
                                
                                if base == 'U':
                                    
                                    print('\033[41m' + base, end="")
                                    
                                elif base == 'A':
                                    
                                    print('\033[44m' + base, end="")
                                    
                                elif base == 'G':
                                    
                                    print('\033[43m' + base, end="")
                                    
                                elif base == 'C':
                                    
                                    print('\033[42m' + base, end="")
                                    
                                else:
                                    
                                    print(base, end="")
                                    
                            print('\033[0m',words[2], words[3],words[4])
                            
                            folded_array.append([i[1], words[3]])
                            
                            #print(aptamer,words[2], words[3])
                            break
                            
                        if float(words[3]) > -2:
                            
                            folded_array.append([i[1], words[3]])
                            
                    elif "no favorable interaction":
                        
                        No_interaction.append(aptamer)

                    else:
                        pass




        interaction_energies_folded = []
        
        for i in folded_array:
            
            interaction_energies_folded.append(float(i[1]))

        n = len(interaction_energies_folded)
        
        interaction_energies = np.array(interaction_energies_folded[:n])

        mean = np.mean(interaction_energies_folded)
        
        median = np.median(interaction_energies_folded)
        
        std_dev = np.std(interaction_energies_folded)
        
        variance = std_dev **2
        
        range_ = np.max(interaction_energies_folded) - np.min(interaction_energies_folded)
        
        min_val = np.min(interaction_energies_folded)
        
        max_val = np.max(interaction_energies_folded)
        
        print("N (size of dataset): ", len(interaction_energies_folded))
        
        print("Mean: ", mean, "kcal/mol")
        
        print("Median: ", median, "kcal/mol")
        
        print("Standard deviation): ", std_dev)
        
        print("Variance: ", variance)
        
        print("Range: ", range_, "kcal/mol")
        
        print("Minimum value: ", min_val, "kcal/mol")
        
        print("Maximum value: ", max_val, "kcal/mol")

        binary_seqs = np.array([list(seq) for seq in folded_best_score]).view(np.uint8)
        
        dist_matrix = pdist(binary_seqs, metric='hamming')

        # Build the tree
        linkage_matrix = linkage(dist_matrix, method='complete')
        
        fig, ax = plt.subplots(figsize=(30, 30))
        
        dendrogram(linkage_matrix, labels= folded_best_score, ax=ax)
        
        plt.savefig(f'best_score{target_name}.png')
        
        plt.show()

        binary_seqs1 = np.array([list(seq1) for seq1 in folded_good_score]).view(np.uint8)
        
        dist_matrix1 = pdist(binary_seqs1, metric='hamming')

        # Build the tree
        linkage_matrix1 = linkage(dist_matrix1, method='complete')
        
        fig, ax = plt.subplots(figsize=(50, 50))
        
        dendrogram(linkage_matrix1, labels= folded_good_score, ax=ax)
        
        plt.savefig(f'good_score_{target_name}.png')
        
        plt.show()


        binary_seqs2 = np.array([list(seq2) for seq2 in folded_better_score]).view(np.uint8)
        
        dist_matrix2 = pdist(binary_seqs2, metric='hamming')

        # Build the tree
        linkage_matrix2 = linkage(dist_matrix2, method='complete')
        
        fig, ax = plt.subplots(figsize=(50, 50))
        
        dendrogram(linkage_matrix2, labels= folded_better_score, ax=ax)
        
        plt.savefig(f'better_score_{target_name}.png')
        
        plt.show()



        print("####################################################################")
        print("#                                                                  #")
        print("#              The results for unfolded Aptamers                   #")
        print("#                                                                  #")
        print("####################################################################") 





        #predicting the interactions for unfolded aptamres in the library


        for i in unfolded_aptamers:
            
            dirname = f"{target_name}_outputs"
            
            sub_dirname = "Unfolded_outputs"

            if not os.path.exists(dirname):
                
                os.mkdir(dirname)
                
                print(f"Creating a directory named {target_name}_outputs....")
                

            if not os.path.exists(os.path.join(dirname, sub_dirname)):
                
                os.mkdir(os.path.join(dirname, sub_dirname))
                
                print(f"Creating a subdirectory named {sub_dirname}....")
                

            filename = f"intoutput_aptamer_{i[0]}_and_{target_name}.log"
            
            filepath = os.path.join(dirname, sub_dirname, filename)

            with open(filepath, "w") as f:
                
                f.write(f"################################################\n")
                f.write(f"#    these results are for folded aptamer {i[0]} #\n")
                f.write(f"#     against the target named {target_name}   #\n")
                f.write(f"#     Please remember to reference Intarna     #\n")
                f.write(f"#     and T-selex when writing report          #\n")
                f.write(f"################################################\n\n")

                cmd = f"IntaRNA -t {target_sequence} -q {i[1]} --outMode=D --seedBP=4"
                
                result = os.popen(cmd).read()

                f.write(result)

                f.write(f"################################################\n")
                f.write(f"#    The results below shows the pairwise      #\n")
                f.write(f"#     for the aptamer against the target       #\n")
                f.write(f"################################################\n\n")

                alignments = pairwise2.align.globalxx(target_sequence,i[1])
                
                for alignment in alignments:
                    
                    f.write(str(alignment))
                    
                    f.write(format_alignment(*alignment))


            f.close()



            with open(filepath,"r") as fin:
                
                lines=fin.readlines()
                
                for line in lines:
                    
                    if "interaction energy" in line:
                        
                        words=line.split()
                        
                        all_interacting_aptamers.append(aptamer)
                        
                        unfolded_interaction_energies.append(float(words[3]))
                        
                        if float(words[3]) < -8:
                            
                            unfolded_best_score.append(str(i[1]))
                            
                            best_score.append(str(i[1]))
                            
                            for base in i[1]:
                                
                                if base == 'U':
                                    
                                    print('\033[41m' + base, end="")
                                    
                                elif base == 'A':
                                    
                                    print('\033[44m' + base, end="")
                                    
                                elif base == 'G':
                                    
                                    print('\033[43m' + base, end="")
                                    
                                elif base == 'C':
                                    
                                    print('\033[42m' + base, end="")
                                    
                                else:
                                    
                                    print(base, end="")
                                    
                            print('\033[0m', words[2], words[3], words[4])
                            
                            unfolded_array.append([i[1], words[3]])

                            break

                        if -8 < float(words[3]) < -6:
                            
                            unfolded_good_score.append(str(i[1]))
                            
                            good_score.append(str(i[1]))
                            
                            for base in i[1]:
                                
                                if base == 'U':
                                    
                                    print('\033[41m' + base, end="")
                                    
                                elif base == 'A':
                                    
                                    print('\033[44m' + base, end="")
                                    
                                elif base == 'G':
                                    
                                    print('\033[43m' + base, end="")
                                    
                                elif base == 'C':
                                    
                                    print('\033[42m' + base, end="")
                                    
                                else:
                                    
                                    print(base, end="")
                                    
                            print('\033[0m',words[2], words[3],words[4])
                            
                            unfolded_array.append([i[1], words[3]])
                            
                            #print(aptamer, words[2], words[3])
                            break

                        if -6 < float(words[3]) < -4:  
                            
                            unfolded_better_score.append(str(i[1]))
                            
                            better_score.append(str(i[1]))
                            
                            for base in i[1]:
                                
                                if base == 'U':
                                    
                                    print('\033[41m' + base, end="")
                                    
                                elif base == 'A':
                                    
                                    print('\033[44m' + base, end="")
                                    
                                elif base == 'G':
                                    
                                    print('\033[43m' + base, end="")
                                    
                                elif base == 'C':
                                    
                                    print('\033[42m' + base, end="")
                                    
                                else:
                                    
                                    print(base, end="")
                                    
                            print('\033[0m',words[2], words[3], words[4])
                            
                            unfolded_array.append([i[1], words[3]])
                            
                            #print(aptamer, words[2], words[3])
                            break

                        if -4 < float(words[3]) < -2:
                            
                            unfolded_poor_score.append(str(i[1]))
                            
                            for base in i[1]:
                                
                                if base == 'U':
                                    
                                    print('\033[41m' + base, end="")
                                    
                                elif base == 'A':
                                    
                                    print('\033[44m' + base, end="")
                                    
                                elif base == 'G':
                                    
                                    print('\033[43m' + base, end="")
                                    
                                elif base == 'C':
                                    
                                    print('\033[42m' + base, end="")
                                    
                                else:
                                    
                                    print(base, end="")
                                    
                            print('\033[0m',words[2], words[3],words[4])
                            
                            unfolded_array.append([i[1], words[3]])
                            
                            #print(aptamer,words[2], words[3])

                            break
                            
                        if float(words[3]) > -2:
                            
                            unfolded_array.append([i[1], words[3]])
                            
                    elif "no favorable interaction":
                        
                        No_interaction.append(aptamer)

                    else:
                        pass



        interaction_energies_unfolded = []
        
        for i in unfolded_array:
            
            interaction_energies_unfolded.append(float(i[1]))

        n = len(interaction_energies_unfolded)
        
        interaction_energies = np.array(interaction_energies_unfolded[:n])

        mean1 = np.mean(interaction_energies_unfolded)
        
        median1 = np.median(interaction_energies_unfolded)
        
        std_dev1 = np.std(interaction_energies_unfolded)
        
        variance1 = std_dev1 **2
        
        range_1 = np.max(interaction_energies_unfolded) - np.min(interaction_energies_unfolded)
        
        min_val1 = np.min(interaction_energies_unfolded)
        
        max_val1 = np.max(interaction_energies_unfolded)
        
        print("N (size of dataset): ", len(interaction_energies_unfolded))
        
        print("Mean: ", mean1, "kcal/mol")
        
        print("Median: ", median1, "kcal/mol")
        
        print("Standard deviation): ", std_dev1)
        
        print("Variance: ", variance1)
        
        print("Range: ", range_1, "kcal/mol")
        
        print("Minimum value: ", min_val1, "kcal/mol")
        
        print("Maximum value: ", max_val1, "kcal/mol")

        binary_seqs = np.array([list(seq) for seq in unfolded_best_score]).view(np.uint8)
        
        dist_matrix = pdist(binary_seqs, metric='hamming')

        # Build the tree
        linkage_matrix = linkage(dist_matrix, method='complete')
        
        fig, ax = plt.subplots(figsize=(30, 30))
        
        dendrogram(linkage_matrix, labels= unfolded_best_score, ax=ax)
        
        plt.savefig(f'unfolded_best_score{target_name}.png')
        
        plt.show()



        binary_seqs1 = np.array([list(seq1) for seq1 in unfolded_good_score]).view(np.uint8)
        
        dist_matrix1 = pdist(binary_seqs1, metric='hamming')

        # Build the tree
        linkage_matrix1 = linkage(dist_matrix1, method='complete')
        
        fig, ax = plt.subplots(figsize=(50, 50))
        
        dendrogram(linkage_matrix1, labels= unfolded_good_score, ax=ax)
        
        plt.savefig(f'unfolded_good_score_{target_name}.png')
        
        plt.show()


        binary_seqs2 = np.array([list(seq2) for seq2 in unfolded_better_score]).view(np.uint8)
        
        dist_matrix2 = pdist(binary_seqs2, metric='hamming')

        # Build the tree
        linkage_matrix2 = linkage(dist_matrix2, method='complete')
        
        fig, ax = plt.subplots(figsize=(50, 50))
        
        dendrogram(linkage_matrix2, labels= unfolded_better_score, ax=ax)
        
        plt.savefig(f'unfolded_better_score_{target_name}.png')
        
        plt.show()






        # Recognising nucleotide patterns using N-gram method 
        number_pattern = 4
        
        no_top_pattern = 10

        nucleotide_pattern = []
        
        for aptamer in best_score:
            
          for i in range(len(aptamer)-number_pattern+1):
            
            nucleotide_pattern.append(aptamer[i:i+number_pattern])

        nucleotide_pattern_frequencies = dict(Counter(nucleotide_pattern))
        
        results1 = sorted(nucleotide_pattern_frequencies.items(), key=lambda x: x[1], reverse=True)[:no_top_pattern]
        
        print(len(best_score))
        
        print(f"The {no_top_pattern} most common {number_pattern}-grams in the words aptamers are: {results1}")



        nucleotide_pattern = []
        
        for aptamer in good_score:
            
          for i in range(len(aptamer)-number_pattern+1):
            
            nucleotide_pattern.append(aptamer[i:i+number_pattern])

        nucleotide_pattern_frequencies = dict(Counter(nucleotide_pattern))
        
        results2 = sorted(nucleotide_pattern_frequencies.items(), key=lambda x: x[1], reverse=True)[:no_top_pattern]
        
        print(len(good_score))
        
        print(f"The {no_top_pattern} most common {number_pattern}-grams in the aptemers are: {results2}")

        nucleotide_pattern = []
        
        for aptamer in better_score:
            
          for i in range(len(aptamer)-number_pattern+1):
            
            nucleotide_pattern.append(aptamer[i:i+number_pattern])

        nucleotide_pattern_frequencies = dict(Counter(nucleotide_pattern))
        
        results3 = sorted(nucleotide_pattern_frequencies.items(), key=lambda x: x[1], reverse=True)[:no_top_pattern]
        
        print(len(better_score))
        
        print(f"The {no_top_pattern} most common {number_pattern}-grams in aptamers are: {results3}")

        fig,ax = plt.subplots()

        x_values = [x[0] for x in results1]
        
        y_values = [x[1] for x in results1]
        
        plt.bar(x_values, y_values)
        
        plt.xlabel('pattern')
        
        plt.ylabel('Frequency')
        
        plt.title('pattern Frequencies for best score')
        
        plt.savefig(f'freq_bestscore_{target_name}.png')
        
        plt.show()


        x_values = [x[0] for x in results2]
        
        y_values = [x[1] for x in results2]
        
        plt.bar(x_values, y_values)
        
        plt.xlabel('pattern')
        
        plt.ylabel('Frequency')
        
        plt.title('parttern Frequencies for good score')
        
        plt.savefig(f'freq_goodscore_{target_name}.png')
        
        plt.show()


        x_values = [x[0] for x in results3]
        
        y_values = [x[1] for x in results3]
        
        plt.bar(x_values, y_values)
        
        plt.xlabel('pattern')
        
        plt.ylabel('Frequency')
        
        plt.title('parttern Frequencies for better score')
        
        plt.savefig(f'freq_betterscore{target_name}.png')
        
        plt.show()


                           #####################################################
                           #                                                   #
                           #   Box and wiska for ufolded and folded aptamers   #
                           #                                                   #
                           #####################################################


        import matplotlib.pyplot as plt

        data = [interaction_energies_unfolded, interaction_energies_folded]

        fig, ax = plt.subplots()

        ax.boxplot(data)

        ax.set_xticklabels(['unfolded', 'folded'])

        ax.set_ylabel('Interaction energy (kcal/mol)')

        ax.set_title('Box and Whiskers Plot')

        plt.show()


        import seaborn as sns
        import pandas as pd
        
        df = pd.DataFrame({'Data': interaction_energies_unfolded + interaction_energies_folded, 'Group': ['Unfolded'] * len(interaction_energies_unfolded) + ['Folded 2'] * len(interaction_energies_folded)})

        sns.kdeplot(data=df, x='Data', hue='Group', fill=False)

        plt.xlabel('Interaction energy (kcal/mol)')
        
        plt.ylabel('Density')
        
        plt.title('Density Plot for Data 1 and Data 2')

        plt.show()





                           #####################################################
                           #                                                   #
                           #      simple optimization of the best aptamer      #
                           #                                                   #
                           #####################################################


        # for folded aptamers 


        min_pair = folded_array[0]

        for pair in folded_array:
            
            pair[1] = float(pair[1])
            
            if pair[1] < min_pair[1]:
                
                min_pair = pair
                
                bst_folded_aptamer = min_pair[0]
                
        print(min_pair)
        
        print(bst_folded_aptamer)


        end_add = []
        
        start_add = []
        
        mid_insertion = []
        
        mid_substitution = []
        
        end_sub = []
        
        start_sub= []
        
        for i in results1:
            
            end_add.append(bst_folded_aptamer + i[0])
            
            start_add.append(i[0] + bst_folded_aptamer)
            
            firsthalf,remainder = divmod(len(bst_folded_aptamer),2)
            
            mid_insertion.append(bst_folded_aptamer[:firsthalf + remainder] + i [0] + bst_folded_aptamer[firsthalf + remainder:])
            position = len(bst_folded_aptamer)//2
            
            mid_substitution.append(bst_folded_aptamer[:position] + i[0] + bst_folded_aptamer[position+4:])
            #end_sub.append()



        for p,i,j,k,l in zip( results1,end_add,start_add,mid_insertion, mid_substitution):
            
            dirname = f"{target_name}_outputs"
            
            sub_dirname = "Folded_simpleOptimization_outputs"

            if not os.path.exists(dirname):
                
                os.mkdir(dirname)
                
                print(f"Creating a directory named {target_name}_outputs....")

            if not os.path.exists(os.path.join(dirname, sub_dirname)):
                
                os.mkdir(os.path.join(dirname, sub_dirname))
                
                print(f"Creating a subdirectory named {sub_dirname}....")

            filename = f"intoutput_aptamer_{p[0]}_and_{target_name}.log"
            
            filepath = os.path.join(dirname, sub_dirname, filename)

            with open(filepath, "w") as f:
                f.write(f"###################################################\n")
                f.write(f"#  These results are for folded optimized aptamer  #\n")
                f.write(f"#     against the target named {target_name}      #\n")
                f.write(f"#     Please remember to reference Intarna        #\n")
                f.write(f"#     and T-selex when writing report             #\n")
                f.write(f"###################################################\n\n")
                
                cmd1 = f"IntaRNA -t {target_sequence} -q {i} --outMode=D --seedBP=4"
                
                cmd2 = f"IntaRNA -t {target_sequence} -q {j} --outMode=D --seedBP=4"
                
                cmd3 = f"IntaRNA -t {target_sequence} -q {k} --outMode=D --seedBP=4"
                
                cmd4 = f"IntaRNA -t {target_sequence} -q {l} --outMode=D --seedBP=4"

                result1 = os.popen(cmd1).read()
                
                result2 = os.popen(cmd2).read()
                
                result3 = os.popen(cmd3).read()
                
                result4 = os.popen(cmd4).read()

                f.write(f"----------------------------------------------------------\n\n")
                
                f.write(f" End_adittion\n")
                
                f.write(f"----------------------------------------------------------\n")
                
                f.write(result1)
                
                f.write(f"----------------------------------------------------------\n\n")

                f.write(f" Start_adittion\n")
                
                f.write(f"----------------------------------------------------------\n")
                
                f.write(result2)
                
                f.write(f"----------------------------------------------------------\n\n")

                f.write(f" Middle_insertion\n")
                
                f.write(f"----------------------------------------------------------\n")
                
                f.write(result3)
                
                f.write(f"----------------------------------------------------------\n\n")

                f.write(f" Middle_substitution\n")
                
                f.write(f"----------------------------------------------------------\n")
                
                f.write(result4)
                
                f.write(f"----------------------------------------------------------\n\n")


            f.close()
            
            print(f"intoutput_aptamer_{p[0]}_and_{target_name}.log")
            
            with open(filepath,"r") as fin:
                
                lines=fin.readlines()
                
                for line in lines:
                    
                    if "interaction energy" in line:
                        
                        words=line.split()
                        
                        print(words[3], words[4])
                        
                        #print(words[3], words[4])
                        #print(words[3],words[4])
                        #print(words[3],words[4])

        ##  print(min(i[0]))


                                #####################################################
                                #                                                   #
                                #              Writing summarry log file            #
                                #                                                   #
                                #####################################################





        with open(f"Summary_file_of_{target_name}.log","w") as f:
            
            f.write(f"##################################################\n")
            f.write(f"#                                                #\n")
            f.write(f"#       Summary file for Aptamers against        #\n")
            f.write(f"#               Mi-RNA target                    #\n")
            f.write(f'#         Please reference T-SELEX               #\n')
            f.write(f'#                                                #\n')
            f.write(f"##################################################\n")
            
            f.write("\n")
            
            f.write(f"Target Name : {target_name} \n")
            
            f.write(f"Target Sequence : {target_sequence}\n")
            
            f.write(f"Number of queries : {len(aptamers)}\n")
            
            f.write(f"Number of query showed interaction : {len(interaction_energies)}\n")

            f.write("\n")
            
            f.write("Descriptive Statistical Analysis for folded aptamers\n")
            
            f.write("------------------------------------------------------\n\n")
            
            f.write(f"N (size of dataset):  {len(interaction_energies_folded)}\n")
            
            f.write(f"Mean: {mean} kcal/mol\n")
            
            f.write(f"Median:  {median} kcal/mol\n")
            
            f.write(f"Standard deviation):  {std_dev}\n")
            
            f.write(f"Variance: {variance}\n")
            
            f.write(f"Range: {range_} kcal/mol\n")
            
            f.write(f"Minimum value: {min_val} kcal/mol\n")
            
            f.write(f"Maximum value: {max_val} kcal/mol\n")
            
            f.write("\n")
            
            f.write("Descriptive Statistical Analysis for Unfolded aptamers\n")
            
            f.write("------------------------------------------------------\n\n")
            
            f.write(f"N (size of dataset): {len(interaction_energies_unfolded)}\n")
        
            f.write(f"Mean: {mean1} kcal/mol\n")
            
            f.write(f"Median: {median1} kcal/mol")
            
            f.write(f"Standard deviation): { std_dev1}\n")
            
            f.write(f"Variance:  {variance1}\n")
            
            f.write(f"Range: {range_1} kcal/mol\n")
            
            f.write(f"Minimum value: {min_val1} kcal/mol")
            
            f.write(f"Maximum value: {max_val1}kcal/mol")
            
            f.write("Patterns and Frequences\n")
            
            f.write("-----------------------\n")
            
            f.write("For or the aptamrs with best score:\n")
            
            f.write(f"{results1}\n")
            
            f.write("\n")
            
            f.write("For or the aptamrs with good score:\n")
            
            f.write(f"{results2}\n")
            
            f.write("\n")
            
            f.write("For or the aptamrs with better score:\n")
            
            f.write(f"{results3}\n")
            
        print(f"##################################################")
        print(f"#                                                #")
        print(f"#        Interactions calculations done         #")
        print(f"#        summury file file saved locally         #")
        print(f'#        Please reference T-SELEX                #')
        print(f'#                                                #')
        print(f"##################################################")
        
    return df1

