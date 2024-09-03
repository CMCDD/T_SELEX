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

