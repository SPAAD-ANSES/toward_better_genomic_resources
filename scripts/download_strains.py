# =============================================================================
# Genomic Data Downloader from NCBI
# =============================================================================
# This script automates the process of dowloading genomic assemblies from NCBI using
# a list of accession numbers provided in a `.tsv` file (e.g., from NCBI Pathogen Detection).
#
# Features:
# - Automatically downloads genome files via the NCBI Datasets CLI.
# - Skips downloads for assemblies that already exist locally.
# - Handles extraction and renaming of downloaded files.
# - Cleans up temporary directories after each download.
#
# Requirements:
# - Conda environment named 'genbank_request' with 'ncbi-datasets-cli' version 15.30.0.
#   Activate the environment before running:
#
#       conda activate genbank_request
# - The conda environment is available in the 'genbank_request.yaml' file.
#
# Notes:
# - This script is designed to be robust and idempotent, meaning it can be run
#   multiple times without redownloading existing data.
# - Make sure your '.tsv' input file contains a column with valid NCBI assembly accessions.


# === Standard Libraries ===
import pandas as pd
import os, subprocess, shutil, glob
import time 



# === Script ===
# Load df 
df_path = "Vibrio-parahaemolyticus_isolates.tsv" # df path from NCBI pathogen detection metadata
df = pd.read_csv(df_path, sep='\t', low_memory=False).fillna('')
df.columns

# Create assembly list from df 
assembly_list = df["Assembly"].unique().tolist()

# Delete empty elments in assembly_list 
assembly_list = [value for value in assembly_list if value != '']

# Count the elements in assembly_list 
count_assemblies = len(assembly_list)

print('Number of assemblies:', count_assemblies)
print('')

not_file = []
empty_file = []

base_path = "../Data/Vibrio-parahaemolyticus/" # output data

## TEST 
# assembly_list = ['GCA_034023795.1']

for accession in assembly_list: 
    assembly_path = os.path.join(base_path, accession + '.fa')
    # Check if the file exist 
    if not os.path.exists(assembly_path): 
        not_file.append(accession) 
        
        try: 
            # download assemblies 
            print(f"Dowloading {accession}...")
            print('')
            subprocess.run(["datasets", "download", "genome", "accession", accession], check=True)
            # Extract the downloaded zip file 
            shutil.unpack_archive('ncbi_dataset.zip')
            # Find the .fna file corresponding to the assembly
            fna_file = glob.glob(f"ncbi_dataset/data/{accession}/*.fna")[0]
            
            # Copy the .fna file to the destination directory with a .fa extension
            shutil.copy(fna_file, assembly_path)
            print(f"Downloaded and copied {fna_file} to {assembly_path}")
            print('')
            
            # Delete the directory 
            shutil.rmtree(os.path.dirname(fna_file))
            print(f"Deleted the directory {os.path.dirname(fna_file)}")
            print('')
            
            time.sleep(2)
            
        except: 
            pass
             
    else: 
        # Check if the file is empty 
        if os.path.getsize(assembly_path) == 0: 
            empty_file.append(accession)
        else:
            continue


print("Number of files that do not exist:", len(not_file))
print("Number of empty files:", len(empty_file))