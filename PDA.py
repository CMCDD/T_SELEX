#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 12:46:13 2023

@author: s1800206
"""



def DMBA(df1,label1, df2,label2,df3,label3, df4,label4,df5,label5, df6,label6):
    
    label1 = str(label1)
    
    # docking multiscale behavior analysis
    
    import pandas as pd
    
    import matplotlib.pyplot as plt
   
    
    if df3 is None and df4 is None and df5 is None and df6 is None and label3 is None and label4 is None and label5 is None and label6 is None:
        
        ligands = df1['ligand'].tolist()

        x = list(range(len(ligands)))
        
        y1 = df1['Energy']
        y2 = df2['Energy']
        
        plt.scatter(x, y1, label=label1, color='red', marker='*', s =1)
        
        plt.scatter(x, y2, label=label2, color='green', marker='^', s =1)
        
        plt.xlabel('best>least docking aptamers')
        
        plt.ylabel('Docking score')
        
        plt.title('Scatter Plot Example')

        plt.legend()

        plt.show()

        print('red:'+label1)
        
        print('green:'+label2)
        
    elif df4 is None and df5 is None and df6 is None and label4 is None and label5 is None and label6 is None:
        
        ligands = df1['ligand'].tolist()

        x = list(range(len(ligands)))
        
        y3 = df3['Energy']
        
        y1 = df1['Energy']
        
        y2 = df2['Energy']
        
        plt.scatter(x, y3, label=label3, color='blue', marker='o', s =1)
        
        plt.scatter(x, y1, label=label1, color='red', marker='*', s =1)
        
        plt.scatter(x, y2, label=label2, color='green', marker='^', s =1)
        
        plt.xlabel('best>least docking aptamers')
        
        plt.ylabel('Docking score')
        
        plt.title('Scatter Plot Example')
        
        print('red:' +label1)
        
        print('green:'+label2)
        
        print('blue:'+label3)
        
        
    elif df5 is None and df6 is None and label5 is None and label6 is None:
        
        ligands = df1['ligand'].tolist()

        x = list(range(len(ligands)))
        
        y3 = df3['Energy']
        
        y1 = df1['Energy']
        
        y2 = df2['Energy']
        
        y4 = df4['Energy']
        
        plt.scatter(x, y3, label=label3, color='blue', marker='o', s =1)
        
        plt.scatter(x, y1, label=label1, color='red', marker='*', s =1)
        
        plt.scatter(x, y2, label=label2, color='green', marker='^', s =1)
        
        plt.scatter(x, y4, label= label4, color='yellow', marker='^', s =1)
        
        plt.xlabel('best>least docking aptamers')
        
        plt.ylabel('Docking score')
        
        plt.title('Scatter Plot Example')
        
        print('red:'+label1)
        
        print('green:'+label2)
        
        print('blue:'+label3)
        
        print('yellow:'+label4)
        
        
    elif df6 is None and label6 is None:
        
        ligands = df1['ligand'].tolist()
       
        x = list(range(len(ligands)))
        
        y3 = df3['Energy']
        
        y1 = df1['Energy']
        
        y2 = df2['Energy']
        
        y4 = df4['Energy']
        
        y5 = df5['Energy']
        
        plt.scatter(x, y3, label=label3, color='blue', marker='o', s =1)
        
        plt.scatter(x, y1, label=label1, color='red', marker='*', s =1)
        
        plt.scatter(x, y2, label=label2, color='green', marker='^', s =1)
        
        plt.scatter(x, y4, label= label4, color='yellow', marker='^', s =1)
        
        plt.scatter(x, y5, label= label5, color='purple', marker='*', s =1)
        
        plt.xlabel('best>least docking aptamers')
        
        plt.ylabel('Docking score')
        
        plt.title('Scatter Plot Example')
        print('red:'+label1)
        print('green:'+label2)
        print('blue:'+label3)
        print('yellow:'+label4)
        print('purple:'+label5)
        
        
    else:
        ligands = df1['ligand'].tolist()
       
        x = list(range(len(ligands)))
        
        y3 = df3['Energy']
        
        y1 = df1['Energy']
        
        y2 = df2['Energy']
        
        y4 = df4['Energy']
        
        y5 = df5['Energy']
        
        y6 = df6['Energy']
        
        plt.scatter(x, y3, label=label3, color='blue', marker='o', s =1)
        
        plt.scatter(x, y1, label=label1, color='red', marker='*', s =1)
        
        plt.scatter(x, y2, label=label2, color='green', marker='^', s =1)
        
        plt.scatter(x, y4, label= label4, color='yellow', marker='^', s =1)
        
        plt.scatter(x, y5, label= label5, color='purple', marker='*', s =1)
        
        plt.scatter(x, y6, label= label6, color='orange', marker='o', s =1)
        
        plt.xlabel('best>least docking aptamers')
        
        plt.ylabel('Docking score')
        
        plt.title('Scatter Plot Example')
        
        print('red:'+label1)
        
        print('green:'+label2)
        
        print('blue:'+label3)
        
        print('yellow:'+label4)
        
        print('purple:'+label5)
        
        print('orange:'+label6)


def BMA(df, receptor_name):
    
    import matplotlib.pyplot as plt
    
    import numpy as np
    
    import sys
    
    import os 
    
    import fnmatch
    
    import csv
    
    import pandas as pd
    
    import os 
    lig_list =[]

    for i in df['ligand']:

        if  len(lig_list) <10:
            lig_list.append(i)
        else:
            break
    
    ligands = df['ligand'].tolist()
    
    home_directory = os.path.expanduser( '~' )

    
    os.chdir(f'{home_directory}/Docking/{receptor_name}/folded')
    
    aptamers_directory = os.listdir(f'{home_directory}/Docking/{receptor_name}/folded')
    
    aptamers_directory
    
    x_list = list(range(1,101))
    
    overal_list = []
    
    energy_list = []
    
    rmsd_list =[]
    
    for directory in ligands[:10]:
        
        present_directory = os.chdir(f'{home_directory}/Docking/{receptor_name}/folded/{directory}')
        
        modelfiles = os.listdir(present_directory)

        for x_value in x_list:

            for file in modelfiles:
                
                name = f'model_{x_value}.pdb'
                
                if name == file:
                    
                    with open (file,'r') as f:
                        
                        lines = f.readlines()
                        
                        for line in lines:
                            
                            if "REMARK Score" in line:
                                
                                words=line.split()
                                
                                energy = float(str(words[2]))
                                
                                
                                energy_list.append(float(energy))


                            if "REMARK RMSD" in line:
                                
                                words=line.split()
                                
                                rmsd = float(str(words[2]))
                                
                                rmsd_list.append(float(rmsd))



    energy_list
    
    print(rmsd_list)

    num_sublists = 10
    
    sublist_size = len(energy_list) // num_sublists

    num_subl = 10
    
    sublsize = len(rmsd_list) // num_subl

    sublists = []
    
    subls =[]
    
    for i in range(0, len(energy_list), sublist_size):
        
        sublist = energy_list[i:i + sublist_size]
        
        sublists.append(sublist)


    for r in range(0, len(rmsd_list), sublsize):
        
        sublist2 = rmsd_list[r:r + sublsize]
        
        subls.append(sublist2)


   
    for sublist in sublists:
        
        pass

    import matplotlib.pyplot as plt

    for sublist in sublists:
        
        plt.plot(x_list, sublist, linestyle='-')

    plt.xlabel('Models')
    
    plt.ylabel(' Docking Score ')
    
    plt.legend(ligands[:10])

    plt.show()


    for sublist in subls:
        
        plt.plot(x_list, sublist, linestyle='-')
        
    plt.xlabel('Models')
    
    plt.ylabel(' RMSD(Å)')
    
    plt.legend(ligands[:10])

    plt.show()

    variances = [np.var(y) for y in subls]

    plot_labels = ligands[:10]

    sorted_variances, sorted_labels = zip(*sorted(zip(variances, plot_labels), reverse=True))

    plt.figure(figsize=(10, 6))
    
    plt.barh(sorted_labels, sorted_variances, color='skyblue')
    
    plt.xlabel('Variance')
    
    plt.title('Line Graph Variance Plot (Ranked by Fluctuation)')
    
    plt.show()

    print("Ranked Line Graphs by Fluctuation:")
    
    for i, label in enumerate(sorted_labels):
        
        print(f"{i+1}. {label}: Variance = {sorted_variances[i]}")
    
    fig = plt.figure()
    
    ax = fig.add_subplot(111, projection='3d')

   
    x = df['Energy']
    
    y = df['Fitness quality'] 
    
    z = df['Confidence Score (%)']  



    x1 = df['Energy'] 
    
    y1 = df['Fitness quality']
    
    z1 = df['Confidence Score (%)']
 
    ax.scatter(x, y, z, c='r', marker='o')  
    
    ax.set_xlabel('Docking score')
    
    ax.set_ylabel('Models quality')
    
    ax.set_zlabel('confidence score (%)')

    
    plt.show()


import math

def PDCA(df, receptor_name):
    
    import matplotlib.pyplot as plt
    
    import numpy as np
    
    import sys
    
    import os 
    
    import fnmatch
    
    import csv
    
    import pandas as pd
    
    import os 
    
    stdev_list = []

    x_list = list(range(1,101))

    overal_list = []

    energy_list = []

    rmsd_list =[]
    ligands = df['ligand'].tolist()
    
    home_directory = os.path.expanduser( '~' )

    for directory in ligands:
        
        directory
        
        present_directory = os.chdir(f'{home_directory}/Docking/{receptor_name}/folded/{directory}')
        
        modelfiles = os.listdir(present_directory)

        for x_value in x_list:

            for file in modelfiles:
                
                name = f'model_{x_value}.pdb'
                
                if name == file:
                    
                    with open (file,'r') as f:
                        
                        lines = f.readlines()
                        
                        for line in lines:
                            
                            if "REMARK Score" in line:
                                
                                words=line.split()
                                
                                energy = float(str(words[2]))
                                
                                energy_list.append(float(energy))

    energy_list

    num_sublists = len(ligands)
    
    sublist_size = len(energy_list) // num_sublists

    sublists = []

    subls =[]

    for i in range(0, len(energy_list), sublist_size):

        sublist = energy_list[i:i + sublist_size]

        sublists.append(sublist)

    def variance(data):

        n = len(data)

        mean = sum(data) / n

        return sum((x - mean) ** 2 for x in data) / (n - 1)

    def standard_deviation(data):

        var = variance(data)

        std_dev = math.sqrt(var)

        return std_dev

    Zscore_list = []

    for i in sublists:

        stdev = standard_deviation(i)

        FQ = round(sum(i)/len(i),2)

        Zscore = (i[0] - FQ)/stdev

        Zscore_list.append(Zscore)

    df['Zcore'] = Zscore_list

    x = list(range(len(ligands)))

    y = df['Zcore']

    plt.scatter(x, y, label= None, color='purple', marker='o', s =4)

    plt.xlabel('best>least docking aptamers')

    plt.ylabel('Zscore')

    plt.title('Scatter Plot Example')

    plt.legend()

    plt.show()

    x =df['Energy']

    y = df['Zcore']

    plt.scatter(x, y, label= None, color='green', marker='o', s =4)

    plt.xlabel('Docking Score ')

    plt.ylabel('Zscore')

    plt.legend()

    plt.show()

    import statsmodels.api as sm

    import statsmodels.formula.api as smf

    def get_r2_statsmodels(x, y, k=1):

        xpoly = np.column_stack([x**i for i in range(k+1)])

        return sm.OLS(y, xpoly).fit().rsquared

    def get_r2_statsmodels_formula(x, y, k=1):

        formula = 'y ~ 1 + ' + ' + '.join('I(x**{})'.format(i) for i in range(1, k+1))

        data = {'x': x, 'y': y}

        return smf.ols(formula, data).fit().rsquared

    get_r2_statsmodels(x,y)

    X = df['Energy']
    
    Y = df['Fitness quality']
    
    X_mean = np.mean(X)
    
    Y_mean = np.mean(Y)

    num = 0
    
    den = 0
    
    for i in range(len(X)):
        
        num += (X[i] - X_mean)*(Y[i] - Y_mean)
        
        den += (X[i] - X_mean)**2
        
    m = num / den
    
    c = Y_mean - m*X_mean

    print(f"y = {round(m,2)}x {round(c,2)}")

    corr_matrix = np.corrcoef(X, Y)
    
    corr = corr_matrix[0,1]
    
    R_sq = corr**2

    print(R_sq)

    plt.scatter(X, Y, label='Fitiness Qualitiy', color='red', marker='o', s =5)
    
    plt.xlabel('Docking Score')
    
    plt.ylabel('Fitness Quality')
    
    Y_pred = m*X + c

    plt.plot([min(X), max(X)], [min(Y_pred), max(Y_pred)], color='green') # predicted
    
    plt.show()
    
    x = list(range(len(ligands)))
    
    y = df['Energy']
    
    y2 = df['Fitness quality']
    
    plt.scatter(x, y, label='Docking Score', color='blue', marker='o', s =5)
    
    plt.scatter(x, y2, label='Fitiness Qualitiy', color='red', marker='o', s =5)
    
    plt.xlabel('best>least docking aptamers')
    
    plt.ylabel('Docking Score')
    
    plt.title('Scatter Plot Example')

    plt.legend()

    plt.show()

