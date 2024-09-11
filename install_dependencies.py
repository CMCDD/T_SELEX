import subprocess
import sys

def run_command(command):
    
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")

def install_conda_packages():
    
    conda_commands = [
        "conda install -y bioconda::viennarna",
        "conda install -y bioconda::intarna",
        "conda install -y conda-forge::biopython",
        "conda install -y conda-forge::scipy",
        "conda install numpy"
    ]
    
    for command in conda_commands:
        run_command(command)

def install_pip_packages():
    
    pip_commands = [
        "pip install selenium",
        "pip install fastcluster",
        "pip install pandas",
        "pip install numpy",
        "pip install matplotlib",
        "pip install seqfold",
        "pip install seaborn",
        "pip install viennarna"
    ]
    
    for command in pip_commands:
        run_command(command)

if __name__ == "__main__":
    print("Starting Conda installations...")
    install_conda_packages()
    
    print("Starting Pip installations...")
    install_pip_packages()
    
    print("All dependencies have been installed.")
