# =============================================================================
# Duplicated Assembly Detection
# =============================================================================
# This script identifies potentially duplicated assemblies based on identical
# metadata and highly similar genomic profiles.

# Steps:
# 1. Filter by identical metadata:
#    - Identify rows that are completely identical across *selected* metadata fields
#      (e.g., source, location, outbreak ID, etc.).
#    - Assemblies with identical metadata are flagged as potential duplicates.
#
# 2. Compare MLST (MultiLocus Sequence Typing) profiles:
#    - For metadata-duplicated rows, check if they share the same Sequence Type (ST).
#    - Matching STs suggest the assemblies may represent the same or very closely
#      related strains.
#
# 3. Evaluate Mash distances:
#    - Use Mash results to assess genomic distances based on k-mer profiles.
#    - Assemblies with distances < 0.001 are considered highly similar
#      (i.e., ≥99.9% identity), suggesting near-identical genomes.
#
# 4. Compare core genome MLST (cgMLST) profiles:
#    - Assemblies with cgMLST allele differences < 1 are considered likely duplicates.
#
# 5. Retain a single representative assembly:
#    - If multiple assemblies still meet all duplication criteria, a single representative 
#     is selected.
#


# === Standard Libraries ===
# Python version: 3.9.13
import os 
import random

# === Data Manipulation ===
import pandas as pd     #Pandas version: 1.4.4
import numpy as np      #Numpy version: 1.21.5


# =============================================================================
# Paths and Input Files
# =============================================================================
# Define base directories 


# === Base directory ===
base_dir = "../duplicated_analysis"

# === Metadata files ===
base_dir_metadata = "../_metadata"
path_metadata_LM = f'{base_dir_metadata}/quality_control_data/LM_metadata_quality.csv'
path_metadata_VP = f'{base_dir_metadata}/quality_control_data/VP_metadata_quality.csv'

# === ST files ===
path_ST_LM = f'{base_dir}/data_mash/LM_ST_CC_Lineage_updated.csv'
path_ST_VP = f'{base_dir}/data_mash/VP_ST_updated.csv'

# === Mash results  ===
mash_dir_LM = "../mash_Listeria/mash_dist_fasta"
mash_dir_VP = "../mash_Vibrio/mash_dist_fasta"
 
# === cgMLST results ===
cgMLST_dir_LM = f'{base_dir}/data_mash/matrix_alleles_cgmlst_LM.tsv'
cgMLST_dir_VP = f'{base_dir}/data_mash/matrix_alleles_cgmlst_VP.tsv'


# === List of metadata columns for deduplication  ===
## These columns provide important information about the context of sampling and genotyping
columns_to_keep = [
    'TaxID',  # Taxonomic identifier for the organism
    'PFGE secondary enzyme pattern',  # PFGE pattern based on secondary enzyme
    'PFGE primary enzyme pattern',  # PFGE pattern based on primary enzyme
    'Host disease',  # Disease associated with the host
    'Outbreak',  # Information about whether the sample is linked to an outbreak
    'Lat/Lon',  # Geographic location (latitude/longitude) where the sample was collected
    'Collection date',  # Date when the sample was collected
    'Collected by',  # Person or organization that collected the sample
    'Stress genotypes',  # Genotypes related to stress resistance
    'Host',  # Host organism from which the sample was taken
    'AST phenotypes',  # Antimicrobial susceptibility testing (AST) results
    'Virulence genotypes',  # Genotypes related to virulence factors
    'Scientific name',  # Scientific name of the organism
    'K-mer group',  # Group based on k-mer analysis
    'Source type',  # Type of source (e.g., food, environment, clinical)
    'IFSAC category',  # IFSAC (Interagency Food Safety Analytics Collaboration) food category
    'AMR genotypes core',  # Core antimicrobial resistance genotypes
    'Serovar',  # Serological classification of the organism
    'Location',  # Specific location or region where the sample was obtained
    'Isolation source',  # Source from which the organism was isolated (e.g., food, water, host)
    'Isolation type',  # Type of isolation (e.g., environment, clinical, food)
    'Food origin',  # Origin of the food sample, if applicable
    'SNP cluster',  # Cluster based on SNP (single nucleotide polymorphism) analysis
    'Min-same',  # Minimum SNP distance within the same cluster
    'Min-diff',  # Minimum SNP distance between different clusters
    'AMR genotypes'  # Antimicrobial resistance genotypes
]


# === Functions  ===
def metadata_type_dfs(df, df_st, group_columns):
    """
    Separates a DataFrame into multiple DataFrames grouped metadata information, and then by ST 

    Parameters:
    df: the metadata df to be grouped and separated.
    df_st: df containing ST information 
    group_columns (list): columns to group by for metadata-based information

    Returns:
    tuple: A tuple containing:
        - separated_dfs (dict): df grouped by metadata.
        - separated_dfs_unique (dict): df grouped metadata.
        - assembly_lists (dict): List of assemblies. key: description of the metadata. Value: list of assemblies.
        - list_assemblies (list): List of lists of assemblies.
        - potential_duplicates (DataFrame): Rows with completely identical metadata (potential duplicates).
    """
    ## Replace NaN, None, or empty values with an empty string 
    df_metadata = df.fillna('')
    df_st = df_st.fillna('')
    
    ## Merge metadata df with the ST df 
    df = pd.merge(df_metadata, df_st, on='Assembly', how='left', indicator=True)
    
    ## Check merge 
    ckeck_merge = df['_merge'].value_counts().to_frame()
    print("Check merge of metadata and ST df: %s"%(ckeck_merge['count'].sum()==len(df_metadata)))
    print(ckeck_merge)
    print('')
    
    # Identify duplicates before filtering the df to ensure full metadata comparison
    potential_duplicates = df_metadata[df_metadata.duplicated(subset=group_columns, keep=False)]
    
    # Count total duplicated rows
    duplicate_count = potential_duplicates.shape[0]
    print(f"Number of possible duplicated rows - Only metadata: {duplicate_count}")
    
    ## Filter the df to keep only the columns in group_columns that exist in the df
    df_metadata_filtered = df_metadata[[col for col in group_columns if col in df_metadata.columns]]
    print("Shape of the filtered Listeria DataFrame:", df_metadata_filtered.shape)
    
    # Initialize dictionaries and lists to store grouped data and results
    separated_dfs = {}
    separated_dfs_unique = {}
    assembly_lists = {}
    list_assemblies = []
    
    # Group by the specified metadata columns
    metadata_groups = df.groupby(group_columns)
    
    for group, groups_df in metadata_groups: 
        
        group_str = "_".join(map(str,group)).replace(' ', '_')
        metadata_key = f'{group_str}'
        
        # Further group by ST within each metadata group 
        st_groups = groups_df.groupby('ST')
        
        for st, st_df in st_groups: 
            key = f'{metadata_key}_ST_{st}'
        
            # If the group contains more than one row, store it in separated_dfs
            if len(st_df) > 1:
                separated_dfs[key] = st_df
                assemblies = st_df['Assembly'].tolist()
                assembly_lists[key] = assemblies
                list_assemblies.append(assemblies)
                
            else:
                # If there's only one row in the group, store it as unique
                separated_dfs_unique[key] = st_df
    
    print("")
    print("Number of dfs to continue with the distance comparison:", len(separated_dfs.keys()))
    print("Number of unique dfs (genomes) after filtering by metadata and ST:", len(separated_dfs_unique.keys()))
    
    flattened_list = [item for sublist in list_assemblies for item in sublist]
    print("Number of total assemblies (to define matrix dimension):", len(flattened_list))

    return separated_dfs, separated_dfs_unique, assembly_lists, list_assemblies



def analyze_assemblies(assembly_lists, mash_dir, patho):
    """
    Analyzes assembly distances by comparing pairs of assemblies and checking their distances from Mash distance files.
    
    Parameters:
    assembly_lists (dict): Dictionary where keys are group identifiers (description) and values are lists of assembly names.
    mash_dir (str): Directory where Mash distance files are stored.
    
    Returns:
        df containing results of comparisons with distances < 0.001.
    """
    data = []  # List to store results for comparisons
    total_comparisons = 0  # Counter for total number of comparisons made
    
    data_sup = [] 
    
    # Iterate over each key-value pair in assembly_lists
    for key, assemblies in assembly_lists.items():
        # Compare each pair of assemblies within the list
        for i in range(len(assemblies)):
            for j in range(i + 1, len(assemblies)):
                assembly1 = assemblies[i]
                assembly2 = assemblies[j]
                
                total_comparisons += 1  
                
                # Construct the filename for the Mash distance file for the first assembly
                file_name = f"{assembly1}.dist"
                file_path = os.path.join(mash_dir, file_name)
                
                # Check if the distance file exists
                if os.path.exists(file_path): 
                    # Check if the file is empty
                    if os.path.getsize(file_path) == 0: 
                        print(file_path, 'is empty')
                    else: 
                        # Open and read the distance file
                        with open(file_path, 'r') as file: 
                            # Look for the line corresponding to the second assembly
                            for line in file: 
                                if patho == 'Listeria': 
                                    row_lm = '../Data/'  ## Fasta directory 
                                    row_name = row_lm
                                elif patho == 'Vibrio' : 
                                    row_vp = '../Data/'
                                    row_name = row_vp
                                else: 
                                    pass 
                                
                                if line.startswith(f"{row_name}{assembly2}.fa"): 
                                    parts = line.strip().split('\t')  # Split the line into parts
                                    value = float(parts[1])  # Extract the distance value
                                    
                                    # If the distance is less than 0.001, record the result
                                    if value < 0.001:
                                        result = { 
                                            'Assembly_1': assembly1,
                                            'Assembly_2': assembly2,
                                            'value': value, 
                                            'Description': key                
                                        }
                                        data.append(result)
                                    else: 
                                        result_ = { 
                                            'Assembly_1': assembly1,
                                            'Assembly_2': assembly2,
                                            'value': value, 
                                            'Description': key                
                                        }
                                        data_sup.append(result_)
                                
                else: 
                    print(file_path, 'does not exist')
    
    # Convert the list of results to a DataFrame
    output_df = pd.DataFrame(data)
    # Print summary statistics
    # print('Total comparisons:', total_comparisons)
    print('Number of comparisons to keep (assemblies < 0.001 mash distance):', len(output_df))
    
    # Df with assemblies comp > 0.001
    output_df_ = pd.DataFrame(data_sup)
    print('Number of comparisons to keep (assemblies > 0.001 mash distance):', len(output_df_))
    
    return output_df, output_df_

def compare_profiles(cgMLST_df, mash_results_df): 
    """
    function to compare allelic profiles of assemblies and determine if they are identical or not. 
    We drop the columns in the df where there are non-allelic information. 
    
    Parameters: 
    cgMLST_df: df contaning the allelic profiles 
    mash_results_df: df contaning assemblies with mash results (duplicated values)
    
    Retunrs: 
    df : df contaning with a new column "same_allelic_profile" where TRUE if all allelic profiles have the same values, 
            or FALSE if there if at least 1 difference between the profiles  
    
    """
    mash_results_df = mash_results_df.copy()
    results = []
    
    # Iterate over each row in mash_results_df
    for index, row in mash_results_df.iterrows():
        assembly_1 = row['Assembly_1']
        assembly_2 = row['Assembly_2']
        
        # Extract corresponding rows (profiles) from cgMLST_df
        profile_1 = cgMLST_df[cgMLST_df['Assembly'] == assembly_1]
        profile_2 = cgMLST_df[cgMLST_df['Assembly'] == assembly_2]

        # Check if both profiles exist in cgMLST_df
        if not profile_1.empty and not profile_2.empty:
            
            # Concatenate both profiles into a temporary DataFrame
            tmp_df = pd.concat([profile_1, profile_2], axis=0)
            
            # Drop columns that contain NaN in either profile
            tmp_df_clean = tmp_df.dropna(axis=1)
            
            # Compare remaining values in the cleaned df
            values_equal = tmp_df_clean.iloc[0, 1:].equals(tmp_df_clean.iloc[1, 1:])
            
            # Add the comparison result (True/False) to the list
            results.append(values_equal)
        else:
            # If either profile is missing, append False
            results.append(False)
            
        if (index + 1) % 1000 == 0: 
            print(f"Processed {index + 1} rows...")

            
    # Add the results to mash_results_df as a new column 'same_allelic_profile'      
    mash_results_df['same_allelic_profile'] = results
        
            
    return mash_results_df

def assemblies_to_drop(df):
    """
    Identifies assemblies to drop based on clustering of assembly pairs. Assemblies that belong to the same cluster are grouped together,
    and for each cluster, only one assembly is retained while others are marked for removal.
    
    Parameters:
    df : df containing pairs of assemblies with their respective distances.

    Returns:
    list: List of assemblies to drop, where each assembly is not the first one in its cluster.
    """
    # Initialize the dictionary with the first pair of assemblies as the first cluster
    cluster = {0: set([df.iloc[0]['Assembly_1'], df.iloc[0]['Assembly_2']])}
    
    # Iterate over the df starting from the second row
    for index, row in df.iloc[1:].iterrows():
        added = False  # Flag to check if the current pair is added to any existing cluster
        
        # Check if the current assembly pair belongs to any existing cluster
        for key, value in cluster.items():
            if row['Assembly_1'] in value or row['Assembly_2'] in value:
                # Add both assemblies to the existing cluster
                value.update([row['Assembly_1'], row['Assembly_2']])
                added = True  # Set flag to True as the pair is added to a cluster
                break
        
        # If the pair doesn't belong to any cluster, create a new cluster
        if not added:
            cluster[len(cluster)] = set([row['Assembly_1'], row['Assembly_2']])
    
    # Initialize a list to collect assemblies to drop
    out_seq = []
    
    # Iterate over all clusters
    for assemblies in cluster.values():
        assemblies_list = list(assemblies)  # Convert the set of assemblies to a list
        random.shuffle(assemblies_list)  # Shuffle the list to randomize the order
        # Add all assemblies except the first one (to retain one assembly per cluster)
        out_seq.extend(assemblies_list[1:])
    
    return out_seq


# === Listeria deduplication  ===
# Read metadata of Listeria 
df_metadata_LM = pd.read_csv(path_metadata_LM, sep=',', low_memory=False)
df_metadata_LM = df_metadata_LM.drop('_merge', axis=1)
print(df_metadata_LM.columns)
print("Shape Listeria metadata dataframe", df_metadata_LM.shape)
print('')

## Read ST df of Listeria 
df_st_LM = pd.read_csv(path_ST_LM)
## Select only assemblies of interest 
df_st_LM = df_st_LM[df_st_LM['Assembly'].isin(df_metadata_LM['Assembly'])]
df_st_LM['ST'] = df_st_LM['ST'].astype("Int64").astype(str)
print("Shape Listeria ST dataframe", df_st_LM.shape)

## Identify duplicated metadata 
separated_dfs, separated_dfs_unique, assembly_lists, list_assemblies = metadata_type_dfs(df_metadata_LM, df_st_LM, columns_to_keep)

## Deduplication using MASH dataframe 
## Calcul of duplicated values
LM_duplicated_assemblies_df, LM_mash_distace_assemblies_df= analyze_assemblies(assembly_lists, mash_dir_LM, 'Listeria')


## Read cgMLST profiles
cgMSLT_lm = pd.read_csv(cgMLST_dir_LM, sep='\t', low_memory=False)
cgMSLT_lm.dtypes

for col in cgMSLT_lm.columns[1:]:
    cgMSLT_lm[col] = cgMSLT_lm[col].astype('Int64')

print("Shape Listeria cgMLST dataframe:", cgMSLT_lm.shape)

## Unique assemblies in duplicated_assemblies df 
unique_assemblies_mash = pd.unique(LM_duplicated_assemblies_df[['Assembly_1', 'Assembly_2']].values.ravel('K'))
print("Number of unique assemblies after filtering using mash distance:", len(unique_assemblies_mash))

## Filter the cgMLST df based on the Assemblies of interest: unique_assemblies_mash
cgMSLT_lm_filtered = cgMSLT_lm[cgMSLT_lm['Assembly'].isin(unique_assemblies_mash)]
print("Shape Listeria filtered cgMLST dataframe:", cgMSLT_lm_filtered.shape)

## Compare allelic profiles to determine duplicated values 
LM_duplicated_assemblies_df_ = compare_profiles(cgMSLT_lm_filtered, LM_duplicated_assemblies_df)
LM_duplicated_assemblies_df_.head(8)

## Number of True and False values after compare allelic profiles of assemblies
count_allelic_profiles = LM_duplicated_assemblies_df_['same_allelic_profile'].value_counts().reset_index()
print("Number of True and False values after compare allelic profiles of assemblies")
count_allelic_profiles

## Drop rows where 'same_allelic_profile' is False
LM_duplicated_assemblies_df__ = LM_duplicated_assemblies_df_[LM_duplicated_assemblies_df_['same_allelic_profile'] == True]
print('Number of comparisons to keep after cgMLST filtering:', len(LM_duplicated_assemblies_df__))

## Unique assemblies in duplicated_assemblies df 
unique_assemblies_cgMLST = pd.unique(LM_duplicated_assemblies_df__[['Assembly_1', 'Assembly_2']].values.ravel('K'))
print("Number of unique assemblies after filtering using cgMLST profile:", len(unique_assemblies_cgMLST))

## Dropping assemblies to create the new df 
## List of assemblies to drop: source attribution 
df_grouped = LM_duplicated_assemblies_df__.groupby('Description')
list_assemblies_to_drop_LM = []

for name, group in df_grouped: 
    result = assemblies_to_drop(group)
    list_assemblies_to_drop_LM.extend(result)

print("Number of assemblies to drop from the df: ", len(list_assemblies_to_drop_LM))

## Save the new dataframe without duplicated data 
df_LM_cleaned = df_metadata_LM[~df_metadata_LM['Assembly'].isin(list_assemblies_to_drop_LM)]
print("Number of assemblies after deduplication - LM : ", len(df_LM_cleaned))

df_LM_cleaned.to_csv(f'{base_dir_metadata}/deduplicated_data/LM_metadata_deduplicated.csv', index=False)

df_LM_cleaned_only_assemblies = df_LM_cleaned['Assembly']
print("Number of assemblies after deduplication - LM : ", len(df_LM_cleaned_only_assemblies))
df_LM_cleaned_only_assemblies.to_csv(f'{base_dir_metadata}/deduplicated_data/LM_assemblies_deduplicated.csv', index=False)

print("Shape Listeria ST dataframe", df_LM_cleaned.shape)

# === Vibrio deduplication  ===
## Read metadata of Listeria 
df_metadata_VP = pd.read_csv(path_metadata_VP, sep=',', low_memory=False)
df_metadata_VP = df_metadata_VP.drop('_merge', axis=1)
print(df_metadata_VP.columns)
print("Shape Vibrio metadata dataframe", df_metadata_VP.shape)
print('')

## Read ST df of Listeria 
df_st_VP = pd.read_csv(path_ST_VP)
## Select only assemblies of interest 
df_st_VP = df_st_VP[df_st_VP['Assembly'].isin(df_metadata_VP['Assembly'])]
df_st_VP['ST'] = df_st_VP['ST'].astype("Int64").astype(str)
print("Shape Vibrio ST dataframe", df_st_VP.shape)

## Identify duplicated metadata 
separated_dfs, separated_dfs_unique, assembly_lists, list_assemblies = metadata_type_dfs(df_metadata_VP, df_st_VP, columns_to_keep)

## Deduplication using MASH dataframe 
## Calcul of duplicated values
VP_duplicated_assemblies_df, VP_mash_distace_assemblies_df = analyze_assemblies(assembly_lists, mash_dir_VP, 'Vibrio')

## Deduplication using cgMLST profiles dataframe 
## Read cgMLST profiles
cgMSLT_vp = pd.read_csv(cgMLST_dir_VP, sep='\t', low_memory=False)
cgMSLT_vp.dtypes

for col in cgMSLT_vp.columns[1:]:
    cgMSLT_vp[col] = cgMSLT_vp[col].astype('Int64')

print("Shape Listeria cgMLST dataframe:", cgMSLT_vp.shape)

## Unique assemblies in duplicated_assemblies df 
unique_assemblies_mash = pd.unique(VP_duplicated_assemblies_df[['Assembly_1', 'Assembly_2']].values.ravel('K'))
print("Number of unique assemblies after filtering using mash distance:", len(unique_assemblies_mash))

## Filter the cgMLST df based on the Assemblies of interest: unique_assemblies_mash
cgMSLT_vp_filtered = cgMSLT_vp[cgMSLT_vp['Assembly'].isin(unique_assemblies_mash)]
print("Shape Vibrio filtered cgMLST dataframe:", cgMSLT_vp_filtered.shape)

## Compare allelic profiles to determine duplicated values 
VP_duplicated_assemblies_df_ = compare_profiles(cgMSLT_vp_filtered, VP_duplicated_assemblies_df)
VP_duplicated_assemblies_df_.head(8)

## Number of True and False values after compare allelic profiles of assemblies
count_allelic_profiles = VP_duplicated_assemblies_df_['same_allelic_profile'].value_counts().reset_index()
print("Number of True and False values after compare allelic profiles of assemblies")
count_allelic_profiles

## Drop rows where 'same_allelic_profile' is False
VP_duplicated_assemblies_df__ = VP_duplicated_assemblies_df_[VP_duplicated_assemblies_df_['same_allelic_profile'] == True]
print('Number of comparisons to keep after cgMLST filtering:', len(VP_duplicated_assemblies_df__))

## Unique assemblies in duplicated_assemblies df 
unique_assemblies_cgMLST = pd.unique(VP_duplicated_assemblies_df__[['Assembly_1', 'Assembly_2']].values.ravel('K'))
print("Number of unique assemblies after filtering using cgMLST profile:", len(unique_assemblies_cgMLST))

## Dropping assemblies to create the new df 
## List of assemblies to drop: source attribution 
df_grouped = VP_duplicated_assemblies_df__.groupby('Description')
list_assemblies_to_drop_VP = []

for name, group in df_grouped: 
    result = assemblies_to_drop(group)
    list_assemblies_to_drop_VP.extend(result)

print("Number of assemblies to drop from the df: ", len(list_assemblies_to_drop_VP))

## Save the new dataframe without duplicated data 
df_VP_cleaned = df_metadata_VP[~df_metadata_VP['Assembly'].isin(list_assemblies_to_drop_VP)]
print("Number of assemblies after deduplication - VP : ", len(df_VP_cleaned))

df_VP_cleaned.to_csv(f'{base_dir_metadata}/deduplicated_data/VP_metadata_deduplicated.csv', index=False)

df_VP_cleaned_only_assemblies = df_VP_cleaned['Assembly']
print("Number of assemblies after deduplication - VP : ", len(df_VP_cleaned_only_assemblies))
df_VP_cleaned_only_assemblies.to_csv(f'{base_dir_metadata}/deduplicated_data/VP_assemblies_deduplicated.csv', index=False)

print("Shape Vibrio ST dataframe", df_VP_cleaned.shape)
