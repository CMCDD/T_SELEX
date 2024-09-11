import time
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.service import Service 
service = Service() 
import selenium
print(selenium.__version__)
# this is important
options = webdriver.ChromeOptions() 
options.add_argument("--headless=new")
driver = webdriver.Chrome(service=service, options=options)
import numpy as np
from fastcluster import linkage
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
import math
import time
import random
import csv
from IPython.display import FileLink
from seqfold import fold, dg, dg_cache, dot_bracket
import RNA
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import pandas as pd
import subprocess

# RNA aptamer generator function

def gen_aptamers(seed,length, aptamers_num):
    
    if seed == None:
        if length == "randomize":

            aptamers = set()
            
            while len(aptamers)< aptamers_num:

                p= list(range(16,60))

                l = random.choices(p)

                aptamer = ''.join(random.choices('ACUG', k=l[0]))   
                
                if aptamer not in aptamers:

                    aptamers.add(aptamer)
        else:

            aptamers = set()

            aptamer = ''.join(random.choices('ACUG', k=length))

            while len(aptamers)< aptamers_num:

                aptamer = ''.join(random.choices('ACUG', k=length))   
                
                if aptamer not in aptamers:

                    aptamers.add(aptamer)
    else:

        random.seed(seed)

        if length == "randomize":

            aptamers = set()
            
            while len(aptamers)< aptamers_num:

                p= list(range(16,60))

                l = random.choices(p)

                aptamer = ''.join(random.choices('ACUG', k=l[0]))   
                
                if aptamer not in aptamers:

                    aptamers.add(aptamer)
        else:

            aptamers = set()

            aptamer = ''.join(random.choices('ACUG', k=length))

            while len(aptamers)< aptamers_num:

                aptamer = ''.join(random.choices('ACUG', k=length))   
                
                if aptamer not in aptamers:

                    aptamers.add(aptamer)
                                                        
    return list(aptamers)

def fold_and_composition(aptamers_list):

  import RNA

  data =[]

  for i, aptamer in enumerate(aptamers_list):

      G = 0

      A = 0

      C = 0

      U = 0

      for char in aptamer:

          if char == "G":

              G += 1

          if char == "A":

              A += 1

          if char == "C":

              C += 1

          if char == "U":

              U += 1

      (mfe_structure,mfe_energy) = RNA.fold(aptamer)



      data.append([i+1, aptamer, U, G, A, C, mfe_structure,mfe_energy ])

  # Define the column headers
  headers = ['Number', 'Aptamer', 'Us', 'Gs', 'As', 'Cs','MFE structure','Minimum free Energy']

  import csv

  with open('data.csv', 'w', newline='') as csvfile:

      writer = csv.writer(csvfile)

      writer.writerow(headers)

      for row in data:

          writer.writerow(row)

  import pandas as pd

  df = pd.read_csv('data.csv')
  print("csv file name data.csv is created and save with all data")
  
  return df


def tertiary_structure(aptamer_list, secondary_structure):
   import os
   from selenium import webdriver
   from selenium.webdriver.chrome.service import Service
   from selenium.webdriver.chrome.options import Options
   from selenium.webdriver.common.by import By
   from selenium.webdriver.common.keys import Keys
   from selenium.webdriver.support.ui import WebDriverWait
   from selenium.webdriver.support import expected_conditions as EC
   import time
   import pandas as pd
   import os
    chrome_options = Options()
    chrome_options.add_argument("--headless")  
    chrome_options.add_argument("--no-sandbox")
    chrome_options.add_argument("--disable-dev-shm-usage")

    # Set up the driver
    driver = webdriver.Chrome(options=chrome_options)

    def process_aptamer(name, aptamer, sec):
        try:
            print(f"Processing aptamer {name} with secondary structure method {sec}")

            driver.get("https://rnacomposer.cs.put.poznan.pl/")
            print("Navigated to RNAcomposer.")

            # Wait for the page to load
            #time.sleep(20)

            search_box = driver.find_element(By.NAME, "content")
            search_box.clear()
            search_box.send_keys(f'>aptamerd{name}')
            search_box.send_keys(Keys.RETURN)
            search_box.send_keys(aptamer)
            search_box.send_keys(Keys.RETURN)
            search_box.send_keys(sec)
            search_box.send_keys(Keys.RETURN)

            compose = driver.find_element(By.NAME, "send")
            compose.click()
            print(f"Submitted aptamer {name} for folding.")

            # Wait for the process to complete
            #time.sleep(25)

            results = driver.find_elements(By.CLASS_NAME, "task-log")
            output_file_name = f"output_aptamer{name}.txt"
            with open(output_file_name, "w") as file:
                for result in results:
                    file.write(result.text + "\n")
            print(f"Results for aptamer {name} saved to {output_file_name}.")

            pdb_link_text = "aptamerd" if sec == "IPknot" else "aptamer"
            pdb = WebDriverWait(driver, 200).until(
                EC.presence_of_element_located((By.PARTIAL_LINK_TEXT, pdb_link_text))
            )
            pdb.click()
            print(f"PDB file for aptamer {name} downloaded.")

            # Record downloaded aptamer
            with open('downloaded_aptamers.txt', 'a') as downloaded_file:
                downloaded_file.write(f"{name}\n")

            # Record last successfully processed aptamer
            with open('last_processed_aptamer.txt', 'w') as last_processed_file:
                last_processed_file.write(f"{name}\n")

        except Exception as e:
            print(f"An error occurred while processing aptamer {name}: {e}")
            # Record the last aptamer attempted but failed
            with open('last_failed_aptamer.txt', 'w') as last_failed_file:
                last_failed_file.write(f"{name}\n")

    def resume_processing(start_index):
        nonlocal secondary_structure
        if isinstance(secondary_structure, pd.Series):
            secondary_structure = secondary_structure.tolist()

        downloaded_aptamers = set()
        if os.path.exists('downloaded_aptamers.txt'):
            with open('downloaded_aptamers.txt', 'r') as downloaded_file:
                downloaded_aptamers = set(line.strip() for line in downloaded_file)

        last_failed_aptamer = 1
        if os.path.exists('last_failed_aptamer.txt'):
            with open('last_failed_aptamer.txt', 'r') as last_failed_file:
                try:
                    last_failed_aptamer = int(last_failed_file.read().strip())
                except ValueError:
                    last_failed_aptamer = 1

        if isinstance(secondary_structure, list):
            for name, (aptamer, sec) in enumerate(zip(aptamer_list, secondary_structure), start=1):
                if str(name) in downloaded_aptamers:
                    continue
                if name < last_failed_aptamer:
                    continue
                process_aptamer(name, aptamer, sec)
        elif secondary_structure in ["default", "CF", "Context", "CONTRA", "IPknot", "RNAstructure"]:
            for name, aptamer in enumerate(aptamer_list, start=1):
                if str(name) in downloaded_aptamers:
                    continue
                if name < last_failed_aptamer:
                    continue
                process_aptamer(name, aptamer, secondary_structure)
        else:
            print("Unsupported secondary structure method.")

    # Determine the starting index
    start_index = 1
    if os.path.exists('last_processed_aptamer.txt'):
        with open('last_processed_aptamer.txt', 'r') as last_processed_file:
            try:
                start_index = int(last_processed_file.read().strip()) + 1
            except ValueError:
                start_index = 1

    resume_processing(start_index)

    # Close the driver
    driver.quit()
    print("Driver closed.")

    finished_text = "Please note that the calculations are done..."
    return finished_text






                                    ############################################################################
                                    #                                                                          #
                                    #                      Neighrest naigbour model for RNA                    #
                                    #                               New model                                  #
                                    #                                                                          #
                                    ############################################################################




new_model = {

    'GC/CG': {'ΔG°37': -3.46, 'ΔH°': -16.52, "ΔS°":-42.13},

    'CC/GG': {'ΔG°37': -3.28, 'ΔH°': -13.94, "ΔS°":-34.41},

    'GA/CU': {'ΔG°37': -2.42, 'ΔH°': -13.75, "ΔS°":-36.53},

    'CG/GC': {'ΔG°37': -2.33, 'ΔH°': -9.61, "ΔS°":-23.46},

    'AC/UG': {'ΔG°37': -2.25, 'ΔH°': -11.98, "ΔS°":-31.37},

    'CA/GU': {'ΔG°37': -2.07, 'ΔH°': -10.47, "ΔS°":-27.08},

    'AG/UC': {'ΔG°37': -2.01, 'ΔH°': -9.34, "ΔS°":-23.66},

    'UA/AU': {'ΔG°37': -1.29, 'ΔH°': -9.16, "ΔS°":-25.40},

    'AU/UA': {'ΔG°37': -1.09, 'ΔH°': -8.91, "ΔS°":-25.22},

    'AA/UU': {'ΔG°37': -0.94, 'ΔH°': -7.44, "ΔS°":-20.90},

    "Initiation": {"ΔG°37": 4.10,"ΔH°": 4.66,"ΔS°": 1.78},

    "Symmetry": {"ΔG°37":0.43,"ΔS°":-1.38},

    "AU End on AU": {"ΔG°37": 0.22,"ΔH°": 4.36,"ΔS°": 13.35},

    "AU End on CG": {"ΔG°37": 0.44, "ΔH°": 3.17,"ΔS°": 8.79},

    'GC/UG': {'ΔG°37': -2.23, 'ΔH°': -14.73, "ΔS°":-40.32},

    'CU/GG': {'ΔG°37': -1.93, 'ΔH°': -9.26, "ΔS°":-23.64},

    'GG/CU': {'ΔG°37': -1.80, 'ΔH°': -12.41, "ΔS°":-34.23},

    'CG/GU': {'ΔG°37': -1.05, 'ΔH°': -5.64, "ΔS°":-14.83},

    'AU/UG': {'ΔG°37': -0.76, 'ΔH°': -9.23, "ΔS°":-27.32},

    'GA/UU': {'ΔG°37': -0.60, 'ΔH°': -10.58, "ΔS°":-32.19},

    'UG/GU': {'ΔG°37': -0.38, 'ΔH°': -8.76, "ΔS°":-27.04},

    'UA/GU': {'ΔG°37': -0.22, 'ΔH°': -2.72, "ΔS°":-8.08},

    'GG/UU': {'ΔG°37': -0.20, 'ΔH°': -9.06, "ΔS°":-28.57},

    'GU/UG': {'ΔG°37': -0.19, 'ΔH°': -7.66, "ΔS°":-24.11},

    'AG/UU': {'ΔG°37': 0.02, 'ΔH°': -5.10, "ΔS°":-16.53},

    'GGUC/CUGG': {'ΔG°37': -3.80, 'ΔH°': -32.49, "ΔS°":-92.57},

    'AU End on GU': {'ΔG°37': -0.71, 'ΔH°': 5.16,"ΔS°": 18.96},

    'GU End on CG': {'ΔG°37': 0.13, 'ΔH°': 3.91, "ΔS°": 12.17},

    'GU End on AU': {'ΔG°37': -0.31, 'ΔH°': 3.65, "ΔS°":-12.78},

    'GU End on GU': {'ΔG°37': -0.74, 'ΔH°': 6.23, "ΔS°":22.47}

}



def thermodynamics_properties(rna_sequence):

    gibb_start =0

    enthalpy_start =0

    entropy_start =0

    d = '/'

    b =[]

    c = []

    gibbs_pairing_list = []

    enthalpy_pairing_list = []

    entropy_pairing_list = []

    self_complement_strant = rna_sequence[::-1]

    features = []

    pairings = []

    complement_pairs = []

    RNA_bases = ['A', 'U', 'G', 'C']

    #constant

    R = 0.001987203611

    ########________________##############

    for i in range(len(rna_sequence) - 1):

        if rna_sequence[i] in RNA_bases and rna_sequence[i+1] in RNA_bases:

            pairing = rna_sequence[i] + rna_sequence[i+1]

            if pairing in ['AU', 'AG', 'AC', 'UA', 'UG', 'UC', 'GA', 'GU', 'GC', 'CA', 'CU', 'CG','UU','CC','GG','AA']:

                pairings.append(pairing + d)

    for i in range(len(self_complement_strant) - 1):

        if self_complement_strant[i] in RNA_bases and self_complement_strant[i+1] in RNA_bases:

            pairing = self_complement_strant[i] + self_complement_strant[i+1]

            if pairing in ['AU', 'AG', 'AC', 'UA', 'UG', 'UC', 'GA', 'GU', 'GC', 'CA', 'CU', 'CG','UU','CC','GG','AA']:

                complement_pairs.append(pairing)

    for i,j in zip(pairings,complement_pairs):

        features.append(i+j)


    middle_index_features = len(features)/2

    middle_index_features = round(float(middle_index_features)) -1

    #if middle_index_features________

    #print(middle_index_features)

    for i in features:

        for j in new_model:

            if j ==i and i == features[middle_index_features]:

                #print(i)
                for p in new_model[j]:

                    if p == "ΔG°37":

                        #print(p + ": " + str(new_model[j][p]))

                        gibbs_pairing_list.append(new_model[j][p])

                    if p == 'ΔH°':

                        #print(p + ": " + str(new_model[j][p]))

                        enthalpy_pairing_list.append(new_model[j][p])

                    if p == 'ΔS°':

                        #print(p + ": " + str(new_model[j][p]))

                        entropy_pairing_list.append(new_model[j][p])


            elif j ==i and i!= features[middle_index_features]:

                for p in new_model[j]:

                    if p == "ΔG°37":

                        #print(p + ": " + str(new_model[j][p]))

                        gibbs_pairing_list.extend([new_model[j][p],new_model[j][p]])

                    if p == 'ΔH°':

                        #print(p + ": " + str(new_model[j][p]))

                        enthalpy_pairing_list.extend([new_model[j][p], new_model[j][p]])

                    if p == 'ΔS°':

                        #print(p + ": " + str(new_model[j][p]))

                        entropy_pairing_list.extend([new_model[j][p], new_model[j][p]])




    #for i in gibbs_pairing_list:

        #b.extend([i,i])

   #b.remove(b[0])

    #print(enthalpy_pairing_list)

    #print(sum(gibbs_pairing_list))

    #print(sum(enthalpy_pairing_list))

    #print(sum(entropy_pairing_list))

    end_terminalAU_on_AU = [ "AU/UA","AA/UU","UU/AA","AU/UA"]

    end_terminalAU_on_CG = [ "GU/CA","GA/CU","CU/GA","CA/GU"]

    end_terminalAU_on_GU = [ "GA/UU","GU/UA","UA/GU","UU/GA"]

    end_terminalGU_on_CG = [ "GA/CU","GU/CA","CA/GU","CU/GA"]

    end_terminalGU_on_AU = [ "AG/UU","AU/UG","UG/AU","GU/AG"]

    end_terminalGU_on_GU = [ "GG/UU","GU/UG","UG/GU","UU/GG"]

    for i in end_terminalAU_on_AU:

        if features[-1] ==i:

            gibb_start = new_model["AU End on AU"]["ΔG°37"]

            enthalpy_start = new_model["AU End on AU"]["ΔH°"]

            entropy_start = new_model["AU End on AU"]["ΔS°"]

        else:

            pass

    for j in end_terminalAU_on_CG:

        if features[-1] ==j:

            gibb_start = new_model["AU End on CG"]["ΔG°37"]

            enthalpy_start = new_model["AU End on CG"]["ΔH°"]

            entropy_start = new_model["AU End on CG"]["ΔS°"]
        else:

            pass

    for p in end_terminalAU_on_GU:

        if features[-1] ==p:

            gibb_start = new_model["AU End on GU"]["ΔG°37"]

            enthalpy_start = new_model["AU End on GU"]["ΔH°"]

            entropy_start = new_model["AU End on GU"]["ΔS°"]

        else:

            pass

    for p in end_terminalGU_on_CG:

        if features[-1] ==p:

            gibb_start = new_model["GU End on CG"]["ΔG°37"]

            enthalpy_start = new_model["GU End on CG"]["ΔH°"]

            entropy_start = new_model["GU End on CG"]["ΔS°"]

        else:

            pass

    for p in end_terminalGU_on_AU:

        if features[-1] ==p:

            gibb_start = new_model["GU End on AU"]["ΔG°37"]

            enthalpy_start = new_model["GU End on AU"]["ΔH°"]

            entropy_start = new_model["GU End on AU"]["ΔS°"]

        else:

            pass

    for p in end_terminalGU_on_GU:

        if features[-1] ==p:

            gibb_start = new_model["GU End on GU"]["ΔG°37"]

            enthalpy_start = new_model["AU End on GU"]["ΔH°"]

            entropy_start = new_model["AU End on GU"]["ΔS°"]

        else:

            pass

    #print(enthalpy_start)

    #print(gibb_start)

    total_gibbs_energy = new_model["Symmetry"]["ΔG°37"] + new_model["Initiation"]["ΔG°37"] + sum(gibbs_pairing_list) +2* gibb_start

    total_enthalpy_energy = new_model["Initiation"]["ΔH°"] + sum(enthalpy_pairing_list) + 2* enthalpy_start

    total_entropy =  new_model["Symmetry"]["ΔS°"] + new_model["Initiation"]["ΔS°"] + sum(entropy_pairing_list) + 2* entropy_start

    Tm = total_enthalpy_energy / ((total_enthalpy_energy-total_gibbs_energy)/310.15+R*math.log(0.0001/1))-273.15

    Tm_adjust = (total_enthalpy_energy / ((total_enthalpy_energy-total_gibbs_energy)/310.15+R*math.log(0.0001/1))-273.15) +16.6*math.log10(0.05)

    return features, round(total_gibbs_energy,2), round(total_enthalpy_energy,2), round(total_entropy,2), round(Tm,2)



# Mass calculation for each sequence

def mass(sequence):

    G = 0

    A = 0

    C = 0

    U = 0

    for char in sequence:

        if char == "G":

            G += 1

        if char == "A":

            A += 1

        if char == "C":

            C += 1

        if char == "U":

            U += 1

    mw = (A * 329.2) + (U * 306.2) + (C * 305.2) + (G * 345.2) + 159

    units = "g/mol"

    return  mw





                                    ############################################################################
                                    #                                                                          #
                                    #                     algorithm for downloading and cleaning               #
                                    #                             aptamaer base from git                       #
                                    #                                                                          #
                                    ############################################################################
def aptamerbase(n_type):

  import seaborn as sns

  import matplotlib.pyplot as plt

  from statsmodels.distributions.empirical_distribution import ECDF

  import numpy as np

  # this data base is downloaded from aptmerbase projects
  # the article to reference when using this baseset
  #(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3308162/)
  # this functiong preprocess the data, by removing empty columns


  df = pd.read_csv(url_data)

  nan_indices_column_A = df[df['sequence'].isna()].index

  nan_indices_column_A

  drop_list =[]

  for i in nan_indices_column_A:

    #print(i)

    drop_list.append(i)

  df = df.drop(df.index[drop_list])


  if n_type == 'RNA':

    df['sequence'] = df['sequence'].astype(str)

    df = df[~df['sequence'].str.contains('T|t', na=False)]

    df = df.reset_index(drop=True)

    df['sequence'] = df['sequence'].str.upper()


  elif n_type == 'DNA':

    df['sequence'] = df['sequence'].astype(str)

    df = df[~df['sequence'].str.contains('U|u', na=False)]

    df = df.reset_index(drop=True)

    df['sequence'] = df['sequence'].str.upper()


  else:
      error1 =print("You entered wrong n_type")

      error2 =print("for RNA aptamers .... aptamerbase(n_type ='RNA')")

      error3 =print("for DNA aptamers .... aptamerbase(n_type ='DNA')")

      df = {error1,
            error2,
            error3}

  #count plot for leagth analysis

  lengths_seq = []

  sequences = df['sequence'].astype(str).tolist()

  for sequence in sequences:

    length_of_sequence = len(sequence)

    lengths_seq.append(length_of_sequence)

  df['sequence_length'] = lengths_seq

  

  hist_seql, edges_seql = np.histogram(df['sequence_length'], bins=np.arange(min(df['sequence_length']), max(df['sequence_length']) + 2))

  plt.figure(figsize=(12, 8))

  plt.bar(edges_seql[:-1], hist_seql, width=0.4, align='center', color='red')

  plt.xlabel('sequence_length')

  plt.ylabel('Count')

  plt.show()

  # Plot the cumulative EDF
  data = df['sequence_length']

  ecdf = ECDF(data)

  # Plot the cumulative EDF

  plt.step(ecdf.x, ecdf.y, label='Cumulative EDF')

  plt.xlabel('sequence_length')

  plt.ylabel('Cumulative Probability')

  plt.title('Cumulative Empirical Distribution Function')

  plt.legend()

  plt.show()


  # for plots with mean and meadian

  plt.step(ecdf.x, ecdf.y, label='Cumulative EDF', linewidth=2, color='blue')

  # Adding lines for mean and median

  mean_value = data.mean()

  median_value = data.median()

  print(f'Mean : {mean_value}')

  print(f'Meadian: {median_value}')

  plt.axvline(mean_value, color='red', linestyle='dotted', label='Mean')

  plt.axvline(median_value, color='green', linestyle='dotted', label='Median')

  plt.xlabel('sequence_length')

  plt.ylabel('Cumulative Probability')

  plt.legend()

  plt.show()

  return df





                                    ############################################################################
                                    #                                                                          #
                                    #                     Comparation functions  for aptamers                  #
                                    #                             input should be a list                       #
                                    #                                                                          #
                                    ############################################################################

def find_similar_aptamers(aptamers, aptamers_tags, n):

    similar_aptamers = []

    highest_scores = {tag: 0 for tag in aptamers_tags}

    for i in range(len(aptamers)):

        for j in range(i+1, len(aptamers)):

            count = 0

            similar_positions = []

            for k in range(len(aptamers[i])):

                if aptamers[i][k] == aptamers[j][k]:

                    count += 1

                    similar_positions.append(k)

            if count >= n:

                score_i = highest_scores[aptamers_tags[i]]

                score_j = highest_scores[aptamers_tags[j]]

                if count > score_i and count > score_j:

                    highest_scores[aptamers_tags[i]] = count

                    highest_scores[aptamers_tags[j]] = count

                    similar_aptamers.append([aptamers_tags[i], aptamers_tags[j], similar_positions, count])

    return similar_aptamers





def SSC(input_data, output_csv, min_score):
    def find_similar_aptamers(aptamers, aptamers_tags, n):
        similar_aptamers = []
        highest_scores = {tag: 0 for tag in aptamers_tags}

        for i in range(len(aptamers)):
            for j in range(i + 1, len(aptamers)):
                count = 0
                similar_positions = []

                for k in range(len(aptamers[i])):
                    if aptamers[i][k] == aptamers[j][k]:
                        count += 1
                        similar_positions.append(k)

                if count >= n:
                    score_i = highest_scores[aptamers_tags[i]]
                    score_j = highest_scores[aptamers_tags[j]]

                    if count > score_i and count > score_j:
                        highest_scores[aptamers_tags[i]] = count
                        highest_scores[aptamers_tags[j]] = count
                        similar_aptamers.append([aptamers_tags[i], aptamers_tags[j], similar_positions, count])

        return similar_aptamers

    
    df = input_data
    df["Aptamer-tag"] = df['Number'].astype(str) + "-" + df['Aptamer'].astype(str)

    aptamers = df['Aptamer'].tolist()
    aptamers_tags = df['Aptamer-tag'].tolist()

    
    header = ['Aptamer1', 'Aptamer2', 'Position', 'Score', 'CombinedAptamers']

    with open(output_csv, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)

        existing_combined_aptamers = set()

        for score in range(len(aptamers[0]) - 1, min_score - 1, -1):
            similar_aptamers = find_similar_aptamers(aptamers, aptamers_tags, score)

            for aptamer_list in similar_aptamers:
                combined_aptamer = aptamer_list[0] + " and " + aptamer_list[1]

                if combined_aptamer in existing_combined_aptamers:
                    continue

                existing_combined_aptamers.add(combined_aptamer)

                row = [aptamer_list[0], aptamer_list[1], '|'.join(map(str, aptamer_list[2])), aptamer_list[3], combined_aptamer]
                writer.writerow(row)
    return df

# Example usage
"""SSC('aptamers.csv', 'SimilarAptamers.csv', 10)"""

#df = pd.read_csv('aptamers.csv')

#df["Aptamer-tag"] = df['Number'].astype(str) +"-"+ df['Aptamer'].astype(str)


#aptamers = df['Aptamer'].tolist()

#aptamers_tags = df['Aptamer-tag'].tolist()





#header = ['Aptamer1', 'Aptamer2', 'Position', 'Score', 'CombinedAptamers']

#with open("SimilarAptamers.csv", "w", newline='') as f:

#    writer = csv.writer(f)

#    writer.writerow(header)
#
 #   existing_combined_aptamers = set()

#    for score in range(len(aptamers[0])-1, 9, -1):

 #       similar_aptamers = find_similar_aptamers(aptamers, aptamers_tags, score)

#        for aptamer_list in similar_aptamers:

#           combined_aptamer = aptamer_list[0] +" and "+ aptamer_list[1]

#            if combined_aptamer in existing_combined_aptamers:

#                continue

 #           existing_combined_aptamers.add(combined_aptamer)

#            row = [aptamer_list[0], aptamer_list[1], '|'.join(map(str, aptamer_list[2])), aptamer_list[3], combined_aptamer]

#            writer.writerow(row)

#           #print(f"Aptamer tags {aptamer_list[0]} and {aptamer_list[1]} have {aptamer_list[2]} matches at positions {aptamer_list[3]} (score {score})")




def process_aptamers_for_plot(csv_file, excel_file):
    sys.setrecursionlimit(10000)


	# Polygenetic tree






	#bin_seqs = np.array([list(aptamers) for seq in aptamers]).view(np.uint8)

	# Calculate pairwise distances between sequences
	#dist_matrix = pdist(bin_seqs, metric='hamming')

	#linkage_matrix = linkage(dist_matrix, method='complete')
	#fig, ax = plt.subplots(figsize=(100, 100))
	#dendrogram(linkage_matrix, labels=aptamers, ax=ax)


	#plt.show()


	# Converting csv file to excel file so i remove the duplicates

    # Convert CSV file to Excel file
    csvfile = pd.read_csv(csv_file)
    
    with pd.ExcelWriter(excel_file) as writer:
    
        csvfile.to_excel(writer, index=False)
        
    print(f"Excel file is generated named {excel_file} (this was converted from CSV file)")

    df = pd.read_excel(excel_file)

    # Fill NaN values and create pivot tables
    df['Aptamer1'] = df['Aptamer1'].fillna('NULL')
    
    df2 = df.pivot_table(columns=['Aptamer1'], aggfunc='size')
    
    print(df2)

    df['Aptamer2'] = df['Aptamer2'].fillna('NULL')
    
    df3 = df.pivot_table(columns=['Aptamer2'], aggfunc='size')
    
    print(df3)

    # Prepare list of aptamers for plotting
    Aptamers1 = df['Aptamer1'].tolist()
    
    Aptamers2 = df['Aptamer2'].tolist()
    
    new_list_for_plot = set()

    for tag1, tag2 in zip(Aptamers1, Aptamers2):
    
        if tag1 not in new_list_for_plot:
        
            new_list_for_plot.add(tag1)
            
        if tag2 not in new_list_for_plot:
        
            new_list_for_plot.add(tag2)

    print(len(new_list_for_plot))
    
    if len(new_list_for_plot) != len(set(new_list_for_plot)):
    
        print("The list contains repeated values.")
        
    else:
    
        print("The list does not contain repeated values.")

    
    dictionary = {}
    
    aptamers = list(new_list_for_plot)  
    
    for tag in new_list_for_plot:
    
        if aptamers:
        
            dictionary[tag] = aptamers.pop(0)
            
        else:
        
            dictionary[tag] = None

    
    aptamers_for_plot = []
    
    for tag in new_list_for_plot:
    
        aptamer = dictionary.get(tag)
        
        if aptamer is not None:
        
            aptamers_for_plot.append(aptamer)

    print(aptamers_for_plot)
    
    binary_seqs = np.array([list(aptamer) for aptamer in aptamers_for_plot]).view(np.uint8)
    
    dist_matrix = pdist(binary_seqs, metric='hamming')

    
    linkage_matrix = linkage(dist_matrix, method='complete')
    
    fig, ax = plt.subplots(figsize=(50, 50))
    
    dendrogram(linkage_matrix, labels=aptamers_for_plot, ax=ax)

    plt.show()







def get_mfe_structure(sequence):
    
    result = subprocess.run(
        ['RNAfold', '--noPS'],
        input=sequence.encode(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    
    output = result.stdout.decode()
    
    mfe_structure = output.split('\n')[1].strip()
    
    return mfe_structure

def extract_loops_and_pairs(sequence, structure):
    
    loops = []
    pairs = []
    current_region = []
    inside_loop = False
    inside_pair = False
    
    for char, base in zip(structure, sequence):
    
        if char == '.':
        
            if inside_pair:
            
                # Closing pair region
                pairs.append(''.join(current_region))
                
                current_region = []
                
                inside_pair = False
                
            if not inside_loop:
            
                inside_loop = True
                
                current_region = [base]
                
            else:
            
                current_region.append(base)
                
        elif char in '()':
        
            if inside_loop:
            
                # Closing loop region
                loops.append(''.join(current_region))
                
                current_region = []
                
                inside_loop = False
                
            if not inside_pair:
            
                inside_pair = True
                
                current_region = [base]
                
            else:
            
                current_region.append(base)
                
    if inside_loop:
    
        loops.append(''.join(current_region))
        
    if inside_pair:
    
        pairs.append(''.join(current_region))
        
    return loops, pairs

def count_base_frequencies(sequences):
    
    base_counts = pd.DataFrame(index=['A', 'U', 'G', 'C'], columns=['Loop Count', 'Pair Count'])
    
    base_counts = base_counts.fillna(0)
    
    overall_counts = pd.DataFrame(index=['A', 'U', 'G', 'C'], columns=['Loop Count', 'Pair Count', 'Overall Count'])
    
    overall_counts = overall_counts.fillna(0)
    
    for sequence in sequences:
    
        loops, pairs = extract_loops_and_pairs(sequence, get_mfe_structure(sequence))
        
        # Count loop bases
        for loop in loops:
            for base in base_counts.index:
                base_counts.loc[base, 'Loop Count'] += loop.count(base)
        
        # Count paired bases
        for pair in pairs:
            for base in base_counts.index:
                base_counts.loc[base, 'Pair Count'] += pair.count(base)
    
    overall_counts['Loop Count'] = base_counts['Loop Count']
    
    overall_counts['Pair Count'] = base_counts['Pair Count']
    
    overall_counts['Overall Count'] = overall_counts['Loop Count'] + overall_counts['Pair Count']
    
    return base_counts, overall_counts

def analyze_sequences(sequences):
    
    results = []
    base_counts = pd.DataFrame(index=['A', 'U', 'G', 'C'], columns=['Loop Count', 'Pair Count'])
    base_counts = base_counts.fillna(0)
    
    for sequence in sequences:
        mfe_structure = get_mfe_structure(sequence)
        loops, pairs = extract_loops_and_pairs(sequence, mfe_structure)
        counts, overall_counts = count_base_frequencies([sequence])
        
        results.append({
            'Sequence': sequence,
            'Loops': loops,
            'Pairs': pairs,
            'Base Frequencies': counts
        })
        
        
        base_counts += counts

    
    overall_counts = count_base_frequencies(sequences)[1]

    return results, base_counts, overall_counts

def write_log(file_path, results, base_counts, overall_counts):
    
    with open(file_path, 'w') as file:
    
        for result in results:
        
            file.write(f"Sequence: {result['Sequence']}\n")
            
            file.write("Extracted Loops:\n")
            
            for loop in result['Loops']:
            
                file.write(f"{loop}\n")
                
            file.write("Extracted Pairs:\n")
            
            for pair in result['Pairs']:
            
                file.write(f"{pair}\n")
                
            file.write("Base Frequencies in Loops and Pairs:\n")
            
            file.write(f"{result['Base Frequencies']}\n")
            
            file.write("\n" + "-"*40 + "\n")
        
        file.write("Overall Base Frequencies:\n")
        
        file.write(f"{base_counts}\n")
        
        file.write("Aggregated Base Frequencies:\n")
        
        file.write(f"{overall_counts}\n")

def loops(sequences, log_file):

    results, base_counts, overall_counts = analyze_sequences(sequences)
    
    # Print the results to the console
    for result in results:
        print(f"Sequence: {result['Sequence']}")
        print("Extracted Loops:")
        for loop in result['Loops']:
            print(loop)
        print("Extracted Pairs:")
        for pair in result['Pairs']:
            print(pair)
        print("Base Frequencies in Loops and Pairs:")
        print(result['Base Frequencies'])
        print("\n" + "-"*40 + "\n")
    
   
    write_log(log_file, results, base_counts, overall_counts)




#loops(sequences, log_file)

