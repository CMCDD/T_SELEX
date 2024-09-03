# T_SELEX


<div style="text-align: center;">
    <img src="pictures/logo.png" alt="T_SELEX" style="width: 90%; height: auto;">
</div>


## Introduction
T_SELEX is a Python package designed for generating RNA aptamer libraries and predicting their secondary and tertiary structures, RNA-RNA interactions, and performing macromolecular docking. The package is primarily designed to run on a Linux machine and offers a structured workflow for high-throughput RNA aptamer analysis. 

## Prerequisites

- **Operating System:** Linux
- **Python Version:** Python 3.x
- **Required Python Packages:** `pandas`, `matplotlib`, `biopython`
- **External Dependencies:**
  - Selenium
  - ViennaRNA
  - Hdock
  - RNAComposer
  - IntaRNA

## Installation

### Internal Dependencies

To install the required internal dependencies, execute the following command in the terminal from the home directory:

Step 1: Create a Conda environment named t_selex

```bash
conda create --name t_selex python=3.9
```
Step 2: Activate the Conda environment

```bash
conda activate t_selex
```
Step 3: Navigate to the site-packages directory


```bash
cd anaconda3/envs/t_selex/lib/python3.9/site-packages/
```

Step 4: Clone the software repository from GitHub
```bash
git clone https://github.com/CMCDD/T_SELEX.git
```
Step 5: Navigate into the cloned repository directory
```bash
cd T_SELEX
```
Step 6: Run the Python script to install dependencies
```bash
python install_dependencies.py
```

### Hdock extension
- First dowload the stand alone software from the [HDock webserver](http://hdock.phys.hust.edu.cn/) 
- Secondly create the directory called sofware and paste the HDOCKlite-v1.1 directory there.




## Usage(Program)
To run the T_SELEX as program program with the given parameters, you can use the following command:


```bash
t_selex_program.py --gen_aptamers_seed 1 --gen_aptamers_length 33 --gen_aptamers_width 10 --receptor_file /home/s1800206/Downloads/6zsl.pdb --receptor_name "6zsl" --ligands_directory /home/s1800206/hdock_test --software_path software/HDOCKlite-v1.1

```

###### Arguments
- `--gen_aptamers_seed` : Sets the seed for generating aptamers to ensure reproducibility.
- `--gen_aptamers_length` : Specifies the length of the aptamers to be generated.
- `--gen_aptamers_number` : Defines the numbers of the aptamers.
- `--receptor_file` : Provides the path to the receptor file, which in this case is the PDB file 6zsl.pdb.
- `--receptor_name` : Provides assist with managing the output files.
- `--ligands_directory` : Points to the directory where the ligand pdb files are located.
- `--software_path' : Indicates the path to the docking software, in this case, HDOCKlite-v1.1.


## Usage(API)
The T_SELEX package is structured to follow a step-by-step workflow:

#### Step 1: RNA Aptamer Library Generation
##### Option 1: Generate Random Sequences
Use the gen_aptamers() function from the to generate RNA sequences. This function generates random sequences using the Base Randomization Algorithm (BRA) which put the bases of the sequnces randomly.
```python

import pandas as pd
import T_SELEX_program
from T_SELEX_program import gen_aptamers,fold_and_composition,tertiary_structure,Mol_docking_calc

# This will generate 100 RNA sequnces with the random lengths ranging from 16  to 60
RNA = gen_aptamers(seed =1, length="randomize", aptamers_num =100)

```

###### Arguments

- `seed`: Initial value of the random number generator. If `None`, a random seed is used.
- `length`: Length of the sequences. Can be a fixed value or "randomize" for lengths between 16 to 60.
- `aptamers_num`: Number of sequences to generate.

##### Option 2: Download RNA Sequences from Aptamerbase
Use the aptamerbase() function from to dowload the aptamer sequences dataset.

```python

# for downloading DNA aptamers
dna_aptamers = aptamerbase("DNA")

# for downloading DNA aptamers
dna_aptamers = aptamerbase("RNA")

```

#### Step 2: RNA Folding and Secondary Structures

Employ the fold_and_composition() function to accurately predict the Minimum Free Energy (MFE) secondary structures. This step involves folding the RNA sequences generated in Step 1 and predicting their secondary structures based on the calculated energy.
```python
import pandas as pd
import T_SELEX_program
from T_SELEX_program import gen_aptamers,fold_and_composition,tertiary_structure,Mol_docking_calc

# This will generate 100 RNA sequnces with the random lengths ranging from 16  to 60
RNA = gen_aptamers(seed =1, length="randomize", aptamers_num =100)
print(RNA)
df = fold_and_composition(RNA)

```
#### Step 3: Tertiary Structures Prediction

Utilize the tertiary_structure() function to forecast the 3D structures of the sequences. In this step, the secondary structures obtained in Step 2 are used as input to predict the tertiary structures of RNA molecules.

```python
import pandas as pd
import T_SELEX_program
from T_SELEX_program import gen_aptamers,fold_and_composition,tertiary_structure,Mol_docking_calc

# This will generate 100 RNA sequnces with the random lengths ranging from 16  to 60
RNA = gen_aptamers(seed =1, length="randomize", aptamers_num =100)
print(RNA)

df = fold_and_composition(RNA)

#Generating tertiary strucures
tertiary structures = tertiary_structure(df["Aptamer"],df["MFE structure"] )

```

#### Step 4: RNA-RNA Interactions Prediction(This is optinal steps when you have the RNA molecule as the target molecules.)

In oerder to predict RNA-RNA interactions at large scale using sequences of the intarna() function is available. This function predicts the interactions between RNA molecules, providing insights into potential binding partners and complex formation.
```python
import pandas as pd
import T_SELEX_program
from T_SELEX_program import intarna

target_seq = "CCAGAGGUUGUAACGUUGUCUAUAUAUACCCUGUAGAACCGAAUUUGUGUGGUAUCCGUAUAGUCACAGAUUCGAUUCUAGGGGAAUAUAUGGUCGAUGCAAAAACUUCA"
p= intarna('aptamers.csv', 'Aptamer',target_seq,"Pre_miR10b",True)
print(p)
# 'aptamers.csv' the csvfile that contain library of sequences and 'Aptamer' is the column name


```

###### Arguments
- `aptamers.csv` is the CSV file containing a library of RNA sequences.
- `Aptamer` is the column name in the CSV file that contains the sequences to be analyzed.
- `target_seq` is the RNA sequence for which interactions are predicted.
- `Pre_miR10b` specifies the name of the target sequence.
- `True` indicates additional settings or options for the prediction analysis and the saving of the results.

#### Step 5: Macromolecular Docking

Leverage the Mol_docking_calc() function provided by the Docking.py script for proficient macromolecular docking. This step involves docking the RNA molecules with their target macromolecules to study their binding affinities and structural conformations
```python
import os
import pandas as pd
import T_SELEX_program
from T_SELEX_program import gen_aptamers,fold_and_composition,tertiary_structure,Mol_docking_calc,
import os
aptmers_sequences = gen_aptamers(1,33,10)
print(p)
df = fold_and_composition(aptamers_sequences)
tertiary structures = tertiary_structure(df['Aptamer'],df['MFE structure'] )

home_directory = os.path.expanduser( '~' )
sofware_path = 'software/HDOCKlite-v1.1'
r = '/home/s1800206/Downloads/6zsl.pdb'
l = '/home/s1800206/hdock_test'
s = '/home/s1800206/software/HDOCKlite-v1.1'
virtual_screening = Mol_docking_calc(data_frame= a,MFE_column ='Minimum free Energy',receptor_name= "6zsl",receptor=r,ligands_directory=l,directory_path= s,Ap_folded=True)

```

###### Arguments
- `data_frame`: DataFrame containing aptamer sequences and their MFE structures.
- `MFE_column`: Column name in the DataFrame that holds the MFE values.
- `receptor_name`: Name of the receptor file.
- `receptor`: Path to the receptor PDB file.
- `ligands_directory`: Directory containing ligand files.
- `directory_path`: Path to the HDOCK software directory.
- `Ap_folded`: Boolean indicating wheater to take folded or non folded aptamers.


#### Step 6: Docking Results Analysis

Conduct a comprehensive analysis of docking results utilizing the generated CSV file. This encompasses regression analysis, Z-score calculations, and data visualization through graphical plots.




## Example

This script showcases the use of the T_SELEX_program module for aptamer design and analysis, specifically focusing on interactions with human bovine serum proteins. It begins by generating a set of aptamers with specified parameters, such as the number of aptamers and their length. Following this, the script predicts the secondary structures and folding compositions of these aptamers and calculates their tertiary structures based on the predicted MFE (Minimum Free Energy) structures. The script then sets up the necessary file paths for the receptor proteins human bovine serum proteins and ligand directories required for molecular docking calculations. It performs docking simulations to evaluate how well the aptamers bind to different human bovine serum receptors, saving the results in a designated directory. After completing the docking calculations, the script conducts post-docking analysis using various functions from the T_SELEX_program module. This includes comparing the docking results across multiple receptors to assess the binding interactions and evaluating the quality of these interactions through additional analytical methods. Overall, the workflow integrates aptamer design, docking simulations, and comprehensive analysis to study aptamer-receptor interactions with human bovine serum proteins.

```python
import pandas as pd
import T_SELEX_program
from T_SELEX_program import gen_aptamers, fold_and_composition, tertiary_structure, Mol_docking_calc, SSC, loops, mass, aptamerbase
import os
p = gen_aptamers(1, 33, 10)
print(p)
a = fold_and_composition(p)
print(a)
k = tertiary_structure(a['Aptamer'],a['MFE structure'] )


home_directory = os.path.expanduser('~')
sofware_path = 'software/HDOCKlite-v1.1'
directory_path = os.path.join(home_directory, sofware_path)
directory_path

os.chdir(home_directory)
df1 = pd.read_csv("Updatedaptamers_with_energy.csv")
r = '/home/s1800206/Downloads/1ao6_clean.pdb'
r1 = '/home/s1800206/Downloads/1bm0_clean.pdb'
r2 = '/home/s1800206/Downloads/4l9k_clean.pdb'
l = '/home/s1800206/Downloads'
s = '/home/s1800206/software/HDOCKlite-v1.1'


p = Mol_docking_calc(data_frame=a, MFE_column='Minimum free Energy', receptor_name="1ao6",
                     receptor=r, ligands_directory=l, directory_path=s, Ap_folded=True)

p1 = Mol_docking_calc(data_frame=a, MFE_column='Minimum free Energy', receptor_name="1bm0",
                     receptor=r1, ligands_directory=l, directory_path=s, Ap_folded=True)

p2 = p = Mol_docking_calc(data_frame=a, MFE_column='Minimum free Energy', receptor_name="4l9k",
                     receptor=r2, ligands_directory=l, directory_path=s, Ap_folded=True)


#post dockong analysis
from T_SELEX_program import DMBA,BMA,PDCA


res = DMBA(df1=p,label1="1ao6", df2=p1,label2="1bm0",df3=p2,label3 ="4l9k")
resl1 = BMA(df1=p,label1="1ao6")
resl2 = BMA(df1=p1,label1="1bm0")
resl3 = BMA(df1=p2,label1="4l9k")



reslp1 = PDCA(df1=p,label1="1ao6")
reslp2 = PDCA(df1=p1,label1="1bm0")
reslp3 = PDCA(df1=p2,label1="4l9k")

```


### results

```bash






['UACCAUCUGGCCCGGGAGUAUGCUCAUAGUAAU', 'UCUCAAUUCAGGGGGGUCUCGGGUGUCUGGGAU', 'CCUGGUCGUAAUAACUCCUUCCGCGCCAGUGCU',
 'AUAACGUCUAAGCCUCCAAAUACCCAGCGGCCU', 'ACACAUCAUAAGACUGGAUGAUGCCUCACCUUC', 'ACAAGAAACUAGGACUGCGCCCUAUACGGGAAC',
 'GCGAGCACCGGAUUCGCAAUUGGGAUUUUGUCU', 'GUAUCGUACGACGAUAGGCACUGAUUUGUGUUC', 'CGCUAUCACAGCCCGUUUCCUAACAUCGUCCUC',
'CUGAUAAGCGAAUAUGUUAUUUCCGAAGAAAGG']
csv file name data.csv is created and save with all data
   Number  ... Minimum free Energy
0       1  ...                -1.7
1       2  ...                -8.2
2       3  ...                -6.5
3       4  ...                -2.2
4       5  ...                -4.9
5       6  ...                -3.6
6       7  ...                -5.1
7       8  ...                -8.0
8       9  ...                -1.2
9      10  ...                -0.1

[10 rows x 8 columns]
Processing aptamer 1 with secondary structure method ..((....)).....(((....)))........
Navigated to RNAcomposer.
Submitted aptamer 1 for folding.
Results for aptamer 1 saved to output_aptamer1.txt.
PDB file for aptamer 1 downloaded.
Processing aptamer 2 with secondary structure method ((((......)))).(((((((....)))))))
Navigated to RNAcomposer.
Submitted aptamer 2 for folding.
Results for aptamer 2 saved to output_aptamer2.txt.
PDB file for aptamer 2 downloaded.
Processing aptamer 3 with secondary structure method .((((.(((.............)))))))....
Navigated to RNAcomposer.
Submitted aptamer 3 for folding.
Results for aptamer 3 saved to output_aptamer3.txt.
...........
......
....
...
..
.
####################################################
#                                                  #
#       Please reference T-SELEX program           #
#                                                  #
####################################################
      ligand  Scores  RMSD  Fitness quality  Confidence Score (%)
0  aptamerd1 -352.40  66.02          -263.69                98.284
1  aptamerd7 -341.53  51.74          -245.67                97.876
2  aptamerd9 -339.23  54.46          -249.63                97.779
3  aptamerd2 -320.27  72.24          -248.79                96.787
4  aptamerd8 -319.85  61.13          -244.87                96.761
5  aptamerd6 -311.51  73.91          -244.47                96.196
6  aptamerd3 -298.08  67.27          -236.41                95.081
7  aptamerd4 -297.31  80.35          -251.79                95.008
8  aptamerd5 -285.65  75.85          -234.41                93.779
..............
.............
..........
.....
...
####################################################
#                                                  #
#       Please reference T-SELEX program           #
#                                                  #
####################################################
     ligand  Energy   RMSD  Fitness quality  Confidence Score (%)
0  aptamerd1 -360.74  56.36          -263.76                98.544
1  aptamerd7 -343.07  34.17          -248.75                97.939
2  aptamerd9 -330.47  54.76          -254.31                97.365
3  aptamerd2 -326.52  72.30          -251.16                97.154
4  aptamerd6 -321.52  26.70          -244.65                96.864
5  aptamerd8 -314.52  79.86          -249.88                96.410
6  aptamerd4 -313.51  85.65          -252.20                96.339
7  aptamerd5 -310.66  51.59          -239.80                96.133
8  aptamerd3 -282.08  90.95          -239.60                93.349






```
