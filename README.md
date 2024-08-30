# T_SELEX
T_SELEX
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

```bash
pip install -r requirements.txt
```
/
### External Dependencies
Detailed instructions for installing the external dependencies can be found in the external_dependencies_installation_guide.txt located in the home directory

## Usage
The T_SELEX package is structured to follow a step-by-step workflow:

#### Step 1: RNA Aptamer Library Generation
##### Option 1: Generate Random Sequences
Use the gen_aptamers() function from the secondary.py script. This function generates random sequences using the Base Randomization Algorithm (BRA).
```python
def gen_aptamers(seed, length, aptamers_num):
    # Function implementation
```

###### Arguments

- `seed`: Initial value of the random number generator. If `None`, a random seed is used.
- `length`: Length of the sequences. Can be a fixed value or "randomize" for lengths between 16 to 60.
- `aptamers_num`: Number of sequences to generate.

##### Option 2: Download RNA Sequences from Aptamerbase
Use the aptamerbase() function from the secondary.py script.

```python
def aptamerbase(n_type):
    # Function implementation
```

#### Step 2: RNA Folding and Secondary Structures

Employ the fold_and_composition() function to accurately predict the Minimum Free Energy (MFE) secondary structures. This step involves folding the RNA sequences generated in Step 1 and predicting their secondary structures based on the calculated energy.

#### Step 3: Tertiary Structures Prediction

Utilize the tertiary_structure() function to forecast the 3D structures of the sequences. In this step, the secondary structures obtained in Step 2 are used as input to predict the tertiary structures of RNA molecules.

#### Step 4: RNA-RNA Interactions Prediction

Harness the predictive capability of the intarna() function available in the interactions.py script. This function predicts the interactions between RNA molecules, providing insights into potential binding partners and complex formation.


#### Step 5: Macromolecular Docking

Leverage the Mol_docking_calc() function provided by the Docking.py script for proficient macromolecular docking. This step involves docking the RNA molecules with their target macromolecules to study their binding affinities and structural conformations
```python
import os
import pandas as pd

from T_SELEX import Docking
from Docking import Mol_docking_calc

home_directory = os.path.expanduser( '~' )
sofware_path = 'software/HDOCKlite-v1.1'
directory_path = os.path.join(home_directory, sofware_path)

df1 = pd.read_csv("Updatedaptamers_with_energy.csv")
r = '/home/s1800206/hdock_test/targets/hsa_miR_10b_5p.pdb'
l = '/home/s1800206/hdock_test/ligands/aptamers'
s = '/home/s1800206/software/HDOCKlite-v1.1'


p = Mol_docking_calc(data_frame= df1,MFE_column ='mfe_energiesRNAfold',receptor_name= "hsa_miR_10b_5p",receptor=r,ligands_directory=l,directory_path= s,Ap_folded=True)
print(p)

```
#### Step 6: Docking Results Analysis

Conduct a comprehensive analysis of docking results utilizing the generated CSV file. This encompasses regression analysis, Z-score calculations, and data visualization through graphical plots.
