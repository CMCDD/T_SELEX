#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 15:25:39 2023

@author: Kabelo Mokgopa
"""

import tempfile
import os

def fastafile1(seq):
    
    try:
        
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fa") as temp_file:
            
            temp_file.write(seq.encode('utf-8'))
            
            temp_file.seek(0)
            
            data = temp_file.read()
            
            print("For the requested sequences (temporary file):")
            
            print(data.decode('utf-8'))
            
            name = temp_file.name 
            
            temp_file.close()
            
            return name 
        
    except Exception as e:
        
        print("Temporary file creation failed. Saving to a permanent file.")
        
        with open("seq2.fa", "w") as f:
            
            f.write(seq)
            
        print("For the requested sequences (permanent file):")
        
        print(seq)



def RNAcofold(seq1, seq2):
    
    seq = f'{seq1}&{seq2}'
    
    seq = seq.upper()
    
    name = fastafile1(seq) 
    
    print('=======')
    
    print('RESULTS')
    
    print('=======')   
    
    try:
        
        cmd = f'RNAcofold -a -d2 --noLP < {name} > seq1.out'
        
        os.system(cmd)
        
        with open("seq1.out", "r") as f:
            
            data = f.read()
            
            print(data)
            
    except Exception as e:
        
        cmd = f'RNAcofold -a -d2 --noLP < seq2.fa > seq2.out'
        
        os.system(cmd)
        
        with open("seq2.out", "r") as f:
            
            data = f.read()
            
            print(data)

# Example
#seq1 = "ACGAUCAGAGAUCAGAGCAUACGACAGCAG"

#seq2 = "ACGAAAAAAAGAGCAUACGACAGCAG"

#RNAcofold(seq1, seq2)


def fastafile(seq):
    
    try:
        
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fa") as temp_file:
            
            temp_file.write(seq.encode('utf-8'))
            
            temp_file.seek(0)
            
            data = temp_file.read()
            
            print("For the requested sequences (temporary file):")
            
            print(data.decode('utf-8'))
            
            seq1 = temp_file.name
            
            temp_file.close()
            
            return seq1
        
    except Exception as e:
        
        print("Temporary file creation failed. Saving to permanent file.")
        
        with open("seq2.fa", "w") as f:
            
            f.write(seq)
            
            print("For the requested sequences (permanent file):")
            
            print(seq)

def FullFold(seq):
    
    seq = seq.upper()
    
    seq1 = fastafile(seq)
    
    try:
      
        cmd = f'RNAfold -p -d2 --noLP < {seq1} > seq1.out'
        
        os.system(cmd)
        
        with open("seq1.out","r") as f:
            
            data = f.read()
            
            print(data)
    except Exception as e:
        cmd = f'RNAfold -p -d2 --noLP < seq1.fa > seq2.out'
        
        os.system(cmd)
        
        with open("seq2.out","r") as f:
            
            data = f.read()
            
            print(data)
            
            

def textfile1(seq,sec):
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fa") as temp_file:
            
            temp_file.write(seq.encode('utf-8') + b'\n')
            
            temp_file.write(sec.encode('utf-8') + b'\n')
            
            temp_file.seek(0)
            
            data = temp_file.read()
            
            print("For the requested sequences (temporary file)")
            
            print(data.decode('utf-8'))
            
            text = temp_file.name
            
            temp_file.close()
            
            return text
        
    except Exception as e:
        
        print("Temporary file creation failed. Saving to permanent file.")
        
        with open("text.txt", "w") as f:
            
            f.write(f'{seq} \n' )
            
            f.write(sec)
            
            print("For the requested sequences (permanent file):")
            
            print(seq)

            
def RNAeval(seq,sec):
    
    seq = seq.upper()
    
    sec = sec.upper()
    
    text = textfile1(seq,sec)
    
    try:
        
        print('=======')
        
        print('RESULTS')
        
        print('=======')
        
        cmd = f'RNAeval -v < {text}'
        
        os.system(cmd)
        
        #with open("seq3.out","r") as f:
            #data = f.read()
            #print(data)
    except Exception as e:
        
        text = 'text.txt'
        
        cmd = f'RNAeval -v < {text} '
        
        os.system(cmd)
        #





import tempfile
import os

def fastafile2(seq):
    
    try:
        
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fa") as temp_file:
            
            temp_file.write(seq.encode('utf-8'))
            
            temp_file.seek(0)
            
            data = temp_file.read()
            
            print("For the requested sequences (temporary file):")
            
            print(data.decode('utf-8'))
            
            name2 = temp_file.name 
            
            temp_file.close()
            
            return name2 
        
    except Exception as e:
        
        print("Temporary file creation failed. Saving to a permanent file.")
        
        with open("seq2.fa", "w") as f:
            
            f.write(seq)
            
        print("For the requested sequences (permanent file):")
        
        print(seq)

def RNAup(seq1, seq2):
    
    seq3 = f'{seq1}&{seq2}'
    
    seq3 = seq3.upper()
    
    name2 = fastafile2(seq3) 
    
    print('=======')
    
    print('RESULTS')
    
    print('=======')
    
    
    try:
        
        cmd = f"RNAup -d 2 --noLP --include_both -c 'S' < {name2} > RNAup.out"
        
        os.system(cmd)
        
        with open("RNAup.out", "r") as f:
            
            data = f.read()
            
            print(data)
            
        # Capture the relevant data from RNAup output, avoiding duplicates
        capture_data = False
        
        data_lines = []
        
        seen_lines = set()  
        
        with open("RNA_w25_u1.out", "r") as f:
            
            lines = f.readlines()
            
            for line in lines:
                
                if "pos      u4S        dG" in line:
                    
                    capture_data = True
                    
                elif "pos      u4S" in line:
                    
                    capture_data = False
                    
                if capture_data and line not in seen_lines:
                    
                    data_lines.append(line)
                    
                    seen_lines.add(line)
        
        
        import numpy as np
        
        import pandas as pd
        
        data5 = []
        
        for line in data_lines:
            
            words = line.split()
            
            p = [words[0],words[1],words[2]]
            
            data5.append(p)
            
        data6 = data5[1:]
        
        
        column_names = ["pos", "u4S", "dG"]
        
        df = pd.DataFrame(data6, columns=column_names)
        
        df.replace("NA", np.nan, inplace=True)
        
        df = df.astype(float) 
        
        import pandas as pd

       
      
        last_occurrence_index = None
        
        second_last_occurrence_index = None

        for index, row in df.iterrows():
            
            if row["pos"] == 1:
                
                second_last_occurrence_index = last_occurrence_index
                
                last_occurrence_index = index

        
        if second_last_occurrence_index is not None:
            
            df1 = df.loc[last_occurrence_index:]
            
            df2 = df.loc[second_last_occurrence_index:last_occurrence_index - 1]
            
        else:
            
            df1 = df  
            
            df2 = pd.DataFrame() 

       
        print("First DataFrame (last occurrence):")
        
        print(df1)
        
        print("\nSecond DataFrame (second-to-last occurrence):")
        
        print(df2)

     
        import matplotlib.pyplot as plt
        
        import numpy as np
        
        # Plotting the data
        x = df1["pos"].astype(int) 
        
        y1 = df1["u4S"].astype(float) 
        
        y2 = df1["dG"].astype(float) 

        plt.plot(x, y1, label="u4S")
        
        plt.plot(x, y2, label="dG")
        
        plt.xlabel("Position")
        
        plt.ylabel("ΔGi(kcal/mol)")
        
        plt.legend()
        
        plt.show()
        
        
        
        x = df2["pos"].astype(int)  
        
        y1 = df2["u4S"].astype(float)  
        
        y2 = df2["dG"].astype(float)
        

        plt.plot(x, y1, label="u4S")
        
        plt.plot(x, y2, label="dG")
        
        plt.xlabel("Position")
        
        plt.ylabel("ΔGi(kcal/mol)")
        
        plt.legend()
        
        plt.show()
        
        return df
        
    except Exception as e:
        print(e)





def fastafile3(seq):
    
    import matplotlib.pyplot as plt
    
    import tempfile
    
    import subprocess
    
    import os
    
    import pandas as pd
    
    import csv
    
    try:
        
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fa") as temp_file:
            
            temp_file.write(seq.encode('utf-8'))
            
            temp_file.seek(0)
            
            data = temp_file.read()
            
            print("For the requested sequences (temporary file):")
            
            print(data.decode('utf-8'))
            
            name2 = temp_file.name 
            
            temp_file.close()
            
            return name2 
        
    except Exception as e:
        print("Temporary file creation failed. Saving to a permanent file.")
        with open("seq2.fa", "w") as f:
            f.write(seq)
            name3 = f.name
            return name3
        print("For the requested sequences (permanent file):")
        print(seq)
           
            
def barrier_calculations(seq):
    
    import os
    
    seq = seq.upper()
    
    sequences = fastafile3(seq)
    
    try:
      
        cmd1 = f'RNAsubopt -d2 --noLP -s -e 30.3 < {sequences} > RNAsubopt.out'
        
        os.system(cmd1)
        
        cmd2 = f'barriers --rates -G RNA-noLP --max 50 --minh 0.1 < RNAsubopt.out > barriers.out'
        
        os.system(cmd2)
        
        with open("barriers.out","r") as f:
            
            data = f.read()
            
            print(data)
            
    except Exception as e:
        
        cmd = f'RNAfold -p -d2 --noLP < seq2.fa > seq2.out'
        
        os.system(cmd)
        
        with open("seq2.out","r") as f:
            
            data = f.read()
            
            print(data)
            
    lines = data.strip().split('\n')

    header = lines[0]
    
    rows = lines[1:]
    
    import csv
    
    csv_file = "output.csv"

    with open(csv_file, mode='w', newline='') as file:
        
        writer = csv.writer(file)

        writer.writerow(["Sequence", "SEc Str", "Free Energy", "label", "Energy barrier"])

        for row in rows:
            
            parts = row.split()
            
            sequence = header
            
            pattern = parts[1]
            
            value = parts[2]
            
            number = parts[3]
            
            another_value = parts[4]
            
            writer.writerow([sequence, pattern, value, number, another_value])

    print(f"CSV file '{csv_file}' has been created.")
    
    import pandas as pd
    
    df = pd.read_csv('output.csv')
    
    #print(df)
    return df

def barriers(seq):
    
    p= barrier_calculations(seq)
    
    try:
        
        import subprocess

        file_path = 'tree.ps'  

        subprocess.call(['gv', file_path])
    
    except Exception as e:
        
        import matplotlib.pyplot as plt

        eps_file = 'tree.ps'
        
        fig, ax = plt.subplots(figsize=(10, 10)) 
        
        img = plt.imread(eps_file)
        
        ax.imshow(img)

        ax.axis('off')

        dpi = 600 

        # Save the figure as a high-resolution PDF
        output_file = 'output.pdf'
        
        plt.savefig(output_file, dpi=dpi, format='pdf', bbox_inches='tight', pad_inches=0)

        plt.show()



#######################################################################################


def fastafile11(seq):
    
    try:
        
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fa") as temp_file:
            
            temp_file.write(seq.encode('utf-8'))
            
            temp_file.seek(0)
            
            data = temp_file.read()
            
            print("For the requested sequences (temporary file):")
            
            print(data.decode('utf-8'))
            
            name = temp_file.name 
            
            temp_file.close()
            
            return name 
        
    except Exception as e:
        
        print("Temporary file creation failed. Saving to a permanent file.")
        
        with open("seq2.fa", "w") as f:
            
            f.write(seq)
            
        print("For the requested sequences (permanent file):")
        
        print(seq)

def RNAcofold22(seq1, seq2):
    
    seq = f'{seq1}&{seq2}'
    
    seq = seq.upper()
    
    name = fastafile11(seq) 
    
    print('=======')
    
    print('RESULTS')
    
    print('=======')   
    
    try:
        
        cmd = f'RNAcofold -a -d2 --noLP < {name} > seq1.out'
        
        os.system(cmd)
        
        with open("seq1.out", "r") as f:
            
            data = f.read()
            
            print(data)
            
    except Exception as e:
        
        cmd = f'RNAcofold -a -d2 --noLP < seq2.fa > seq2.out'
        
        os.system(cmd)
        
        with open("seq2.out", "r") as f:
            
            data = f.read()
            
            print(data)

