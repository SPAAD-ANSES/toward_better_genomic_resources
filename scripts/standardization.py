# =============================================================================
# Normalization and Standardization of Source, Geolocation, and Temporal Data
# =============================================================================
# This script performs normalization and standardization of the 'Location' and
# 'Lat/Lon' columns using geocoding. It enriches the dataset by extracting and
# standardizing geographical components such as:
#   - Continent
#   - ISO (3-letter country code)
#   - ISO 2 (2-letter country code)
#   - Country
#   - State
#   - Region
#   - County
#   - City
#   - Latitude
#   - Longitude
#
# Additionally, it merges manually curated FoodEx2 codes into the dataset,
# aligning food names or categories with the dataset.
#
# Finally, the script standardizes the 'Date' column to extract:
#   - Full standardized date (YYYY-MM-DD)
#   - Year 
#   - Only_Year 


# Requirements:
# - Conda environment named 'standardization'.
#   Activate the environment before running:
#
#       conda activate genbank_request
# - The conda environment is available in the 'standardization.yaml' file.



# === Standard Libraries ===
import sys     
import re      
import os 
from datetime import datetime     
import gensim
from gensim.utils import simple_preprocess

# === Data Manipulation ===
import pandas as pd     
import numpy as np      

# === Geolocation & Country Info ===
from geopy.geocoders import Nominatim              # Geocoding (convert address to coordinates)
from geopy.exc import GeocoderTimedOut             # Handle timeouts from geopy
from countryinfo import CountryInfo                # Country metadata
import importlib.metadata
import pycountry                                   # ISO country codes and names

# === Usage explanation ===
if len(sys.argv) < 2:
    print("Usage: python standardization.py <organism> [<df_path>]")
    print("Example 1: python standardization.py Listeria")
    print("Example 2: python standardization.py Vibrio C:/data/vibrio.tsv")
    print("ERROR: You must provide both the organism and the data file path.")
    sys.exit(1)
    
# === Get arguments ===
organims = sys.argv[1] # Listeria Vibrio 

organims_path = ""
foodex_path = ""
if organims == "Listeria": 
    df_path = f"{organims_path}/Listeria_monocytogenes_isolates.tsv"
    foodex_df = f"{foodex_path}/Listeria_to_foodex2.csv"
if organims == "Vibrio": 
    df_path = f"{organims_path}/Vibrio-parahaemolyticus_isolates.tsv"
    foodex_df = f"{foodex_path}/Vibrio_to_foodex.csv"
        
# === Validate file path ===
if not os.path.exists(df_path):
    print(f"ERROR: The file path does not exist: {df_path}")
    sys.exit(1)
    
# === Confirmation output ===
print(f"Organism set to: {organims}")
print(f"Data file path: {df_path}")


# === Functions === 
def split_colon(location):
    """
    Function to separate the colon (:) and add a space between the country name 
    and the supplementary information.
    
    Parameters: 
    - location: String from the Location column representing the country.
    
    returns: 
    if the location have two parts: first part: second part 
    if the location do not have two parts: location 
    """
    # Check if the input is a string
    if isinstance(location, str):
        # Split the string using colon (:) as a delimiter
        parts = location.split(':')
        
        # Check if there are exactly two parts after splitting
        if len(parts) == 2:
            # Return the first part stripped of leading and trailing whitespace,
            # followed by a colon and a space, and then the second part stripped
            # of leading and trailing whitespace
            return parts[0].strip() + ': ' + parts[1].strip()
        else:
            # Return the original location string if there are not exactly two parts
            return location 


def divide_dataframe_into_segments(df, segment_size=10):
    """
    Function to divide a DataFrame into segments of specified size.

    Parameters:
    - df: The DataFrame to be divided.
    - segment_size: The size of each segment. Default is 10.

    Returns:
    A list of DataFrames, each representing a segment of the original DataFrame.
    """
    num_segments = len(df) // segment_size + 1
    segmented_dfs = []

    for i in range(num_segments):
        start_index = i * segment_size
        end_index = (i + 1) * segment_size
        segment_df = df.iloc[start_index:end_index]
        segmented_dfs.append(segment_df)

    return segmented_dfs  

# Create a geolocator
geolocator = Nominatim(user_agent=f"{organims}", timeout=10)

def get_coordinates(place): 
    
    """
    Function to get the coordinates from country or city name.
    
    Parameters: 
    - place: A string representing the name of the country or city for which coordinates are requested.
    
    Returns:
    A list containing two elements: latitude and longitude coordinates.
    If the location cannot be obtained, returns [NaN, NaN].
    """
    coordinates = []
    location = geolocator.geocode(place, language='en')
    if location: 
        latitude = location.latitude
        longitude = location.longitude
        coordinates.extend([latitude, longitude])
    else: 
        coordinates.extend([np.nan, np.nan])
        print(f"Unable to obtain location for: {place}")
    return coordinates

def get_coordinates_with_retry(location, max_retries=3):
    
    """
    Function to retrieve the coordinates from a country or city name after multiple (default: 3) retries.
    This is implemented to handle connection problems or timeouts during geocoding.
    
    Parameters: 
    - location: A string representing the name of the country or city for which coordinates are requested.
    - max_retries = 3 : An integer specifying the maximum number of retry attempts. Default is 3.
    
    Returns:
    A list containing two elements: latitude and longitude coordinates.
    If the location cannot be obtained even after multiple retries, raises an Exception.
    """
    retries = 0
    while retries < max_retries:
        try:
            return get_coordinates(location)
        except GeocoderTimedOut:
            retries += 1
    raise Exception("Geocoding failed after multiple retries")

def get_lat_lon(location):
    """
    Function to find the latitude and longitude by removing the last word from the location.

    Parameters:
    - location: The location from which latitude and longitude are to be extracted.

    Returns:
    A tuple containing latitude and longitude coordinates.
    """
    # Define a regular expression pattern for special characters
    words = re.findall(r'\w+', location)

    for i in range(len(words), 0, -1):
        # Construct partial location by removing the last word
        partial_location = ' '.join(words[:-1]).strip()
        print(f'Partial location of "{location}" is: "{partial_location}"')
        try:
            # Attempt to get coordinates for the partial location
            return get_coordinates_with_retry(partial_location)
        except Exception as e:
            print(f"Error for location {partial_location}: {str(e)}")

def get_lat_lon4(location):
    """
    Function to find the latitude and longitude by considering only the words before ',' in the location.

    Parameters:
    - location: The location from which latitude and longitude are to be extracted.

    Returns:
    A tuple containing latitude and longitude coordinates.
    """
    # Split the location by ':'
    words = location.split(',')
    for i in range(len(words), 0, -1):
        # Construct partial location by joining the words before ':'
        partial_location = ' '.join(words[:-1]).strip()
        print(f'Partial location of "{location}" is: "{partial_location}"')
        try:
            # Attempt to get coordinates for the partial location
            return get_coordinates_with_retry(partial_location)
        except Exception as e:
            print(f"Error for location {partial_location}: {str(e)}")

def get_lat_lon5(location):
    """
    Function to find the latitude and longitude by considering only the first word of each sentence before ',' in the location.

    Parameters:
    - location: The location from which latitude and longitude are to be extracted.

    Returns:
    A tuple containing latitude and longitude coordinates.
    """
    # Split the location by ','
    split_country = location.split(',')
    for country in split_country:
        # Extract the first word of the sentence
        partial_location = country.split()[0]
        print(f'Partial location of "{location}" is: "{partial_location}"')
        try:
            # Attempt to get coordinates for the first word of the sentence
            return get_coordinates_with_retry(partial_location)
        except Exception as e:
            print(f"Error for location {partial_location}: {str(e)}")

def reverse_geocode(lat, lon):
    """
    Function to perform reverse geocoding for a given pair of latitude and longitude coordinates
    and obtain the country information.

    Parameters:
    - lat: A float representing the latitude coordinate.
    - lon: A float representing the longitude coordinate.

    Returns:
    A tuple containing the geographical information derived from the reverse geocoding process.
    The tuple consists of the following elements:
    - city: The name of the city.
    - county: The name of the county.
    - region: The name of the region.
    - state: The name of the state.
    - country: The name of the country.

    If reverse geocoding fails for any reason, returns (None, None) indicating no information found.
    """
    
    try: 
        location = geolocator.reverse(f"{lat}, {lon}", language='en')
        city = location.raw.get('address', {}).get('city', '')
        county = location.raw.get('address', {}).get('county', '')
        region = location.raw.get('address',{}).get('region', '')
        state = location.raw.get('address', {}).get('state', '')
        country = location.raw.get('address', {}).get('country', '')
        #iso_alpha2 = location.raw.get('address', {}).get('ISO3166-2-lvl4', '')
        #iso = location.raw.get('address', {}).get('country_code', '')
        return city, county, region, state, country
    except Exception as e: 
        return None, None
    

def get_iso(country): 
    
    """
    Function to retrieve ISO names of two and three characters for a given country name.

    Parameters:
    - country: A string representing the name of the country.

    Returns:
    A tuple containing ISO names of two and three characters.
    The tuple consists of the following elements:
    - iso_code: The ISO name of two characters.
    - iso_code2: The ISO name of three characters.

    If ISO names cannot be obtained for any reason, returns (None, None) indicating no information found.
    """
    
    if not country: 
        return None 
    try: 
        country_info = CountryInfo(country)
        if country_info:
            iso_code = country_info.iso(2)
            iso_code2 = country_info.iso(3)
            return iso_code, iso_code2
        else: 
            return None 
    except Exception as e: 
        return None

def get_iso_from_Lat_Lon(lat, lon): 
    """
    Function to perform the reverse geocoding for a given pair of Latitude and longitude coordinates
    and obtain the ISO code.

    Parameters:
    - lat: A float representing the latitude coordinate.
    - lon: A float representing the longitude coordinate.
    
    Return: 
    The ISO name of two characteres. 
    
    If the reverse geocoding fails for anty reason, retunr None indicating no information found. 
    """
    try: 
        location = geolocator.reverse(f"{lat}, {lon}", language='en')
        iso = location.raw.get('address', {}).get('country_code', '')
        return iso 
    except Exception as e: 
        return None   
        

def get_continent(country):
    
    """
    Function to retrieve the continent based on the country code.

    Parameters:
    - country: A string representing the ISO country code (two or three characters).

    Returns:
    A string representing the name of the continent to which the country belongs.

    If the country code is not provided or if the continent information cannot be obtained, returns None.
    """
    
    if not country:
        return None
    try: 
        country_info = CountryInfo(country)
        if country_info:
            continent = country_info.region()
            return continent
        else:
            return None
    except Exception as e: 
        return None
 
 
def parse_coordinates(coordenada):
    """ 
    Parse latitude and longitude from a coordinate string
    To change the N,S,W,E coordinates to '+' and '-' 
    
    Parameters:
    - coordenada: A string representing the coordinates.
    
    Returns:
    A tuple containing latitude and longitude as strings.
    """
    # Check if the input is a float value
    if isinstance(coordenada, float): 
        return float('nan'), float('nan')  # Return NaN values if the input is a float
    
    # Split the coordinate string into parts
    partes = coordenada.split()
    
    # Extract latitude and handle negative directions
    latitud = f"{'-' if partes[1].lower() == 's' else ''}{partes[0]}".rstrip()
    
    # Extract longitude and handle negative directions
    longitud = f"{'-' if partes[3].lower() == 'w' else ''}{partes[2]} ".rstrip()
    
    return latitud, longitud


def all_elements_start_with_digit(data):
    """
    Check if all elements in the given list start with a digit.

    Args:
    - data: A list of strings to be checked.

    Returns:
    - True if all elements in the list start with a digit, else False.
    """
    for element in data:
        if not (isinstance(element, str) and re.search(r'^\d', element)):
            return False
    return True



def try_parsing_date(text):
    """
    Function to parse dates from text.

    Parameters:
    - text: The text containing date information.

    Returns:
    A datetime object representing the parsed date if successful, otherwise returns NaN.
    """
    if pd.isna(text):
        return np.nan
    if not isinstance(text, str):
        text = str(text)
        
    if re.search('.*[0-9]{4}\/.*[0-9]{4}', text):
        text = text.split('/')[0]
    for fmt in ('%b-%Y', '%Y', '%d-%b-%Y', '%Y-%m-%d', '%Y-%m'):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            pass
    return np.nan


## === Pipeline ===
## Columns to keep for analysis 
columns_to_keep = ["Assembly", # Assembly ID 
                   'Isolation type', 'Isolation source', 'Host', 'Source type', 'Host disease', # Source 
                   'Location', 'Lat/Lon', # Location 
                   'Collection date'] # Date 


df= pd.read_csv(df_path, sep='\t', low_memory=False)
df = df[columns_to_keep].copy()
print(f"Number of {organims} isolates available on NCBI Pathogens (March 2024): %s"%(len(df)))
print('')

# Count the number of empty values in Assambly column 
empty_values_count_assambly = df.Assembly.isnull().sum()
print("Number of empty values in Assambly column:", empty_values_count_assambly)

## Delete rows with no values in a Assambly column
df = df.dropna(subset=['Assembly'])
print("Dataframe new dimention", df.shape)
print(f"Number of {organims} isolates available after deleting empty values in assembly column: %s"%(len(df)))
print('')

## Columns description
print("Number of unique values in Location column:", df.Location.nunique())
print("Number of unique values in Collection date column:", df['Collection date'].nunique())
matrice_initial = len(df[['Isolation type', 'Isolation source','Host', 'Source type', 'Host disease']].drop_duplicates())
print("Number of unique values describing the source columns: %s"%(matrice_initial))
print('')

## Foodex2 dataframe 

foodex2 = pd.read_csv(foodex_df, sep=';').fillna('')
print("Total number of descriptive metadata in the matrix annotated in FoodEx2 (post preprocess): %s" % (len(foodex2)))

foodex2_ = foodex2[foodex2['Code'] != 'ROOT']
print("Number of descriptive metadata in the matrix annotated in FoodEx2: %s" % (len(foodex2_)))
print("and not annotated: %s"%(len(foodex2) - len(foodex2_)))

foodex2__ = foodex2[foodex2['Ambiguity'] == False]
print("Number of descriptive metadata in the matrix annotated in FoodEx2 without ambiguity: %s" % (len(foodex2__)))
print("and number of annotations considered ambiguous: %s" % (len(foodex2) - len(foodex2__)))

## convert to lowercase the values of Location column 
df.loc[:,'Location_'] = df.loc[:,'Location'].str.lower()
df.loc[:,'Lat/Lon_'] = df.loc[:,'Lat/Lon'].str.lower()

## Rename NaN values 
NaN_val = {
    'nan' : 'NaN', 
    'not collected' : 'NaN', 
    'not available: to be reported later' : 'NaN', 
    '-' : 'NaN', 
    'other' : 'NaN', 
    '' : 'NaN', 
    'none' : 'NaN', 
    'not available' : 'NaN',
    'not provided' : 'NaN',
    ' ' :  'NaN'
}
df.loc[:,'Location_'] = df.loc[:,'Location_'].astype(str).replace(NaN_val)
df.loc[:,'Location_'] = df.loc[:,'Location_'].replace('NaN', np.nan)#.fillna('')

df.loc[:,'Lat/Lon_'] = df.loc[:,'Lat/Lon_'].astype(str).replace(NaN_val)
df.loc[:,'Lat/Lon_'] = df.loc[:,'Lat/Lon_'].replace('NaN', np.nan)#.fillna('')

print("Total of NaN values in Location column:", df['Location_'].isnull().sum())
print("Total of NaN values in Lat/Lon column:", df['Lat/Lon_'].isnull().sum())

## Apply homogenization to the 'Location' column and convert the values of the Location column to uppercase
df.loc[:,'Location_'] = df.loc[:,'Location_'].apply(split_colon).str.upper()

print("")
df = df.copy()
## Rename in the Location column the old country names 
old_country = ['USSR', 'Czechoslovakia', 'Yugoslavia']
## Check if any old country names are present in the 'Location' column
for country in old_country:
    if any(df['Location_'].str.contains(country, case=False, na=False)):
        print(f"The DataFrame contains '{country}' in the 'Location' column.")
    else: 
        pass 
## RENAME old country name 
# if organims == "Listeria": 
df.loc[~df['Location_'].isna() & df['Location_'].str.contains('USSR'), 'Location_'] = 'RUSSIA'

## Check if there are any ocean, sea, or gulf names present in the 'Location' column
if any(df['Location_'].str.contains('OCEAN|SEA|GULF', case=False, na=False)):
    mask = ~df['Location_'].str.contains(':', na=False)
    mask &= df['Location_'].str.contains('OCEAN|SEA|GULF', case=False, na=False)  
    print(f"The DataFrame contains '{df.loc[mask, 'Location_'].tolist()}' in the 'Location' column")
    # Replace the identified rows with NaN
    df.loc[mask, 'Location_'] = np.nan
   
print("")

## Check if the checked words are present in the 'Location' column
mask_Location = df[df['Location_'].str.contains(r'^(?!.*:).*?(?:OCEAN|SEA|GULF)', case=False, na=False)]
mask_Location = mask_Location['Location_'].to_list()

all_old_country = '|'.join(old_country)
mask_old_country = df[df['Location_'].str.contains(all_old_country, case=False, na=False)]
mask_old_country = mask_old_country['Location_'].to_list()

all_mask = mask_old_country + mask_Location
all_mask = '|'.join(all_mask)

if not all_mask: 
    print(f"Oceans and/or old country names still in the dataframe: '{False}'")   
elif df['Location_'].str.contains(all_mask, case=False, na=False).any(): 
    print(f"Oceans and/or old country names still in the dataframe: '{True}'")

print("")

df.loc[~df['Location_'].isna() & df['Location_'].str.contains('USSR'), 'Location'] = 'RUSSIA'

## Check if there are any ocean, sea, or gulf names present in the 'Location' column
if any(df['Location_'].str.contains('OCEAN|SEA|GULF', case=False, na=False)):
    mask = ~df['Location_'].str.contains(':', na=False)
    mask &= df['Location_'].str.contains('OCEAN|SEA|GULF', case=False, na=False)  
    print(f"The DataFrame contains '{df.loc[mask, 'Location_'].tolist()}' in the 'Location' column")
    # Replace the identified rows with NaN
    df.loc[mask, 'Location_'] = np.nan
   
print("")

## Preprossess the dataframe before including the location information 
# Change None to NaN values is location column 
df = df.fillna(value=np.nan)

# Select rows with not empty values in 'Location' 
df_location_ = df[df['Location_'].notnull()]  ##54260 values 

# Drop duplicates in 'Location', 'Lat/Lon'
df_location__ = df_location_.drop_duplicates(subset=['Location_', 'Lat/Lon_'])

# Replace ":" for a white space in the 'Location' column
df_location__ = df_location__.copy()
df_location__.loc[:,'Location_'] = df_location__.loc[:,'Location_'].str.replace(':', ',')

# Unique values in location column
print("Number of rows in the dataframe after dropping empty values in the Location column:", len(df_location_))
print("Number of unique values in the Location column:", df_location__.Location.nunique())
print("Number of unique values after dropping duplicates:", len(df_location__))
print("Dataframe shape:", df_location__.shape)

## Include the Latitude and Logitude column in the dataframe 
## get_coordinates_with_retry will include the Latitude and Longitude columns  
df_location__ = df_location__.copy()
# Divide the DataFrame into segments
segments = divide_dataframe_into_segments(df_location__, segment_size=25)

df_location__.loc[:,'Latitude'] = np.nan  # Initialize Latitude column
df_location__.loc[:,'Longitude'] = np.nan  # Initialize Longitude column

df_location__ = df_location__.copy()

# Add latitude and longitude information to each segment
for segment_df in segments:    
    for index, row in segment_df.iterrows():
        lat, lon = get_coordinates_with_retry(row['Location_'])
        df_location__.loc[index, 'Latitude'] = lat
        df_location__.loc[index, 'Longitude'] = lon


## Select rows with not empty values in 'Location' but Latitude and Longitude columns are empty
df_location___ = df_location__[df_location__['Location_'].notnull() & df_location__['Latitude'].isnull()]

# Unique values in location column
print("Number of rows in the dataframe with empty values in the Latitude column:", len(df_location___))
print("Number of unique values in the Location column with empty values in the Latitude column:", df_location___.Location.nunique())

df_location___ = df_location___.copy()

if len(df_location___) != 0:
    # Divide the DataFrame into segments
    segments = divide_dataframe_into_segments(df_location___, segment_size=10)

    # Add latitude and longitude information to each segment
    for segment_df in segments:    
        for index, row in segment_df.iterrows():
            lat, lon = get_lat_lon(row['Location_'])
            df_location__.loc[index, 'Latitude'] = lat
            df_location__.loc[index, 'Longitude'] = lon
else: 
    pass 

## Select rows with not empty values in 'Location' but Latitude and Longitude columns still empty 
df_location____ = df_location__[df_location__['Location_'].notnull() & df_location__['Latitude'].isnull()]

# Unique values in location column
print("Number of rows in the dataframe with empty values in the Latitude column:", len(df_location____))
print("Number of unique values in the Location column with empty values in the Latitude column:", df_location____.Location.nunique())

## Include the Latitude and Longitude in where the Latitude column is empty 
## get_lat_lon3 to keep the first 2 words in the Location column and
## then search in geolocator the correspondent Latitude and Logitude column 

df_location____ = df_location____.copy()

if len(df_location____) != 0:
    # Divide the DataFrame into segments
    segments = divide_dataframe_into_segments(df_location____, segment_size=10)

    # Add latitude and longitude information to each segment
    for segment_df in segments:    
        for index, row in segment_df.iterrows():
            lat, lon = get_lat_lon5(row['Location_'])
            df_location__.loc[index, 'Latitude'] = lat
            df_location__.loc[index, 'Longitude'] = lon
else: 
    pass 


## Select rows with not empty values in 'Location' but Latitude and Longitude columns still empty 
df_location_____ = df_location__[df_location__['Location_'].notnull() & df_location__['Latitude'].isnull()]

# Unique values in location column
print("Number of rows in the dataframe with empty values in the Latitude column:", len(df_location_____))
print("Number of unique values in the Location column with empty values in the Latitude column:", df_location_____.Location.nunique())

## Include the Latitude and Longitude in where the Latitude column is empty 
## get_lat_lon5 to keep the first word in the Location column and
## then search in geolocator the correspondent Latitude and Logitude column 

df_location_____ = df_location_____.copy()

if len(df_location_____) != 0:
    # Divide the DataFrame into segments
    segments = divide_dataframe_into_segments(df_location_____, segment_size=10)

    # Add latitude and longitude information to each segment
    for segment_df in segments:    
        for index, row in segment_df.iterrows():
            lat, lon = get_lat_lon5(row['Location_'])
            df_location__.loc[index, 'Latitude'] = lat
            df_location__.loc[index, 'Longitude'] = lon
else: 
    pass 

## Include 'City', 'County', 'Region', 'State', 'Country' 'ISO, and 'ISO_2' columns to the dataframe 
## Use of the Latitude and Longitude columns and the function reverse_geocode to find the City, County, Region, State, Country name 
## Use the function get_iso to include the columns ISO (two letters code) and ISO_2 (three letters code) to the dataframe 
## Use of the Country column to get the Continent name using the function get_continent 

df_location__ = df_location__.copy()
result_columns = ['City', 'County', 'Region', 'State', 'Country']
df_location__[result_columns] = df_location__.apply(lambda row: pd.Series(reverse_geocode(row['Latitude'], row['Longitude'])), axis=1)
result_columns_iso = ['ISO', 'ISO_2']
df_location__[result_columns_iso] = df_location__['Country'].apply(get_iso).apply(pd.Series)
df_location__.loc[:,'Continent'] = df_location__.loc[:,'Country'].apply(get_continent)
df_location__.head()


## Check information 
df_location__ = df_location__.copy()
## Filter the dataframe to keep the rows where the Lat/Lon column is not empty 
## Create a Lat and Lon columns using the parse_coordinates function, 
## use the get_iso_from_Lat_Lon function to obtain the ISO column (ISO_comp) 
## Compare the results obtained using the ISO column and the ISO_comp column 

# Filter Lat/Lon column (not empty rows)
df_location__comp = df_location__[df_location__['Lat/Lon_'].notnull()]

df_location__comp = df_location__comp.copy() 
# Parse the Lat/Lon column 
coordinates_df = df_location__comp['Lat/Lon_'].apply(parse_coordinates).apply(pd.Series)
df_location__comp[['Lat', 'Lon']] = coordinates_df

# Add the ISO column to compare 
df_location__comp["ISO_comp"] = df_location__comp.apply(lambda row: pd.Series(get_iso_from_Lat_Lon(row['Lat'], row['Lon'])), axis=1)

print("Number of rows having information in the Lat/Lon column:", len(df_location__comp))

## Compare the ISO and the ISO_comp columns 
df_location__comp.loc[:,'ISO_comp'] = df_location__comp.loc[:,'ISO_comp'].str.upper()
df_location__comp['Comparison'] = df_location__comp['ISO'] == df_location__comp['ISO_comp']

## Check if all values in the 'Comparison' column are True
if df_location__comp['Comparison'].all():
    print("Country names found from the 'Lat/Lon' column and/or the 'Latitude' and 'Longitude' columns (value calculated from the 'Location' column), correspond to the same country:", True)
else: 
    print("Country names found from the 'Lat/Lon' column and/or the 'Latitude' and 'Longitude' columns (value calculated from the 'Location' column), correspond to the same country:", False)
    # Print the row where the comparaison is False 
    print("")
    false_rows = df_location__comp[~df_location__comp['Comparison']]
    print(false_rows[['Location', 'Country', 'Lat', 'Lon', 'Latitude', 'Longitude', 'ISO', 'ISO_comp']])
    
print("")
# Count the occurrences of True and False values in the 'Comparison' column
count_true_false_ = df_location__comp['Comparison'].value_counts().reset_index()
print("Number of occurences in the Comparaison column: ")
print(count_true_false_)

## Add ANTARCTICA information
df_location__ = df_location__.fillna(value=np.nan)
# Check that the "Country", "ISO", and "ISO_2" columns are not empty
if df_location__['Country'].notnull().all() and df_location__['ISO'].notnull().all() and df_location__['ISO_2'].notnull().all():
    print("The Country and ISO columns are not empty:", True)
else: 
    print("Country and/or ISO columns contain empty values: Check Antarctica rows")     
    
    # Print rows with empty values in the 'Country' column
    empty = df_location__[df_location__["Country"].isna()]
    
    # Check if 'Location' starts with 'ANTARCTICA' for filtered rows
    antarctica_empty = empty[empty['Location_'].str.startswith('ANTARCTICA', na=False)]

    # Check if any rows meet the condition
    if not antarctica_empty.empty:
        print("")
        print('Adding country information to Antarctica rows:')
        df_location__.loc[index, 'Country'] = "ANTARCTICA"  
        df_location__.loc[index, 'Continent'] = "ANTARCTICA"  
        df_location__.loc[index, 'ISO'] = 'AQ'
        df_location__.loc[index, 'ISO_2'] = 'ATA'
        

df_location__ = df_location__.fillna(value=np.nan)
# Check that the "Country", "ISO", and "ISO_2" columns are not empty
if df_location__['Country'].notnull().all() and df_location__['ISO'].notnull().all() and df_location__['ISO_2'].notnull().all():
    print("The Country and ISO columns are not empty:", True)
else: 
    print("Country and/or ISO columns contain empty values")
    
    # Drop rows where both 'Location' and 'Country' are NaN (No information found using Lat/Lon column)
    print("Location and Country information were not available: Drop rows ")
    df_location__ = df_location__.dropna(subset=['Location_', 'Country'], how='all')
    
    # Print rows with empty values in the 'Country' column
    empty = df_location__[df_location__["Country"].isna()]
    
    if len(empty) != 0: 
        print("Location with no information country: ", empty["Location_"]) 
        print("")
        print("Add Location using only country name in the Location column")
        segments = divide_dataframe_into_segments(empty, segment_size=10)
    
        # Add latitude and longitude information to each segment
        for segment_df in segments:    
            for index in segment_df.index:
                lat, lon = get_lat_lon5(segment_df.loc[index, 'Location_2'])
                df_location__.loc[index, 'Latitude'] = lat
                df_location__.loc[index, 'Longitude'] = lon
                
                # Fill empty values with country information
                result_columns = ['City', 'County', 'Region', 'State', 'Country']
                df_location__.loc[index, result_columns] = reverse_geocode(lat, lon) 
                iso_data = get_iso(df_location__.loc[index, 'Country'])  
                df_location__.loc[index, ['ISO', 'ISO_2']] = iso_data
                df_location__.loc[index, 'Continent'] = get_continent(df_location__.loc[index, 'Country'])  
                   
    else: 
        empty_ = df_location__[df_location__["ISO"].isna()]
        print(empty_)
        
        if len(empty_) != 0: 
            print("Location with no information country (ISO): ", empty_["Location_"]) 
            print("")
            print("Update country information using ISO code")
            segments = divide_dataframe_into_segments(empty_, segment_size=10)
            
            # Add latitude and longitude information to each segment
            for segment_df in segments:    
                for index in segment_df.index:
                    country_name = get_iso(segment_df.loc[index, 'Location_2'])
                    df_location__.loc[index, ['ISO', 'ISO_2']] = country_name
                    df_location__.loc[index, 'Continent'] = get_continent(df_location__.loc[index, 'ISO'])
   
        else: 
            pass
                            
print("")            
print("Dataframe shape - after adding the City, County, Region, State, Country, ISO, and continent information:", df_location__.shape)

# Check again that the "Country", "ISO", and "ISO_2" columns are not empty after updates
if df_location__['Country'].notnull().all() and df_location__['ISO'].notnull().all() and df_location__['ISO_2'].notnull().all():
    print("The Country and ISO columns are not empty after updates in country information:", True)
else: 
    print("Country and/or ISO columns still contain empty values after updates")
    
    # Drop rows with no values in the 'Country' column
    print("Drop rows with no information in the country column")
    df_location__.dropna(subset=['Country'], inplace=True)

    print("Location with no information country after dropping rows:", df_location__[df_location__["Country"].isna()]["Location_"])

print("")    
print("Dataframe shape - after adding the City, County, Region, State, Country, ISO, and continent information:", df_location__.shape)

## Check the geolocator package. 
## Verify that the country in the Location column and the country in the Country column are the same
## Create a new dataframe to compare the country name between the Location and the Country column (add Location_2 column)

# Dataframe with colon (:) in the Location column 
df_location__check = df_location_.copy()
# List of column names to keep
selected_columns = ['Location_', 'Lat/Lon_']
# Create a new dataframe with only the selected columns
df_location__check = df_location__check[selected_columns]

# Drp NaN values from the Location column 
df_location__check = df_location__check.dropna(subset="Location_")

# Create Location_2 column to keep only the countries names (delete elements after :) in df_location__check dataframe 
df_location__check['Location_2'] = df_location__check['Location_'].astype(str).apply(lambda x: x.split(':')[0].strip()).to_frame()

# Replace ":" for a white space in the 'Location' column
df_location__check.loc[:,'Location_'] = df_location__check.loc[:,'Location_'].str.replace(':', ',')

# Merge the df_location__check and df_location__ based on the 'Location' column
df_location__c_merged = pd.merge(df_location__check, df_location__, on=['Location_', 'Lat/Lon_'], how='inner', indicator=True)

# Count the _merge values in the '_merge' column
if (df_location__c_merged['_merge'].astype(str) == 'both').all():
    print("Check the merge:", True)
else: 
    print("Check the merge:", False)
count_merge = df_location__c_merged['_merge'].value_counts().reset_index()
print("Number of occurences in the _merge column: ")
print(count_merge)

print('')

# Create a mapping dictionary to handle 'USA' and 'United States' , 'KOREA' and 'SOUTH KOREA', 'VIET NAM' and  'VIETNAM', 'REUNION' and 'GUADELOUPE' to France ... 
mapping = {'USA': 'UNITED STATES',
           'PUERTO RICO' : 'UNITED STATES', 
           'KOREA' : 'SOUTH KOREA', 
           'VIET NAM': 'VIETNAM',
           'REUNION': 'FRANCE', 
           'GUADELOUPE' : 'FRANCE', 
           'HONG KONG': 'CHINA', 
           'CZECH REPUBLIC' : 'CZECHIA', 
           'GAMBIA': 'THE GAMBIA', 
           'MUMBAI' : 'INDIA',
           'FRENCH GUIANA': 'FRANCE'}
# Replace 'USA' with 'United States' in 'Location_2'
df_location__c_merged['Location_2'] = df_location__c_merged['Location_2'].replace(mapping)  
    
## Compare the Location_2 and the Country columns 
df_location__c_merged['Country'] = df_location__c_merged['Country'].str.upper()
df_location__c_merged['Comparison'] = df_location__c_merged['Location_2'] == df_location__c_merged['Country']

## Check if all values in the 'Comparison' column are True
if df_location__comp['Comparison'].all():
    print("The information in 'Location_2' column and 'Country' column are the same:", True)
else: 
    print("The information in 'Location_2' column and 'Country' column are the same:", False)
    # Print the row where the comparaison is False 
    print("")
    false_rows = df_location__c_merged[~df_location__c_merged['Comparison']]
    print("Different rows: ")
    print(false_rows[['Location', 'Location_','Location_2', 'Country', 'Latitude', 'Longitude']])
    
print("")
# Count the occurrences of True and False values in the 'Comparison' column
count_true_false_ = df_location__c_merged['Comparison'].value_counts().reset_index()
print("Number of occurences in the Comparaison column: ")
print(count_true_false_)

## Check if all values in the 'Comparison' column are True
if df_location__c_merged['Comparison'].all():
    print("The information in 'Location_2' column and 'Country' column are the same:", True)
else: 
    print("The information in 'Location_2' column and 'Country' column are the same:", False)
    # Print the rows where the comparison is False 
    print("")
    false_rows = df_location__c_merged[~df_location__c_merged['Comparison']]
    
    # Update the Location information columns in tha false rows 
    # Iterate through false_rows and update coordinates 
    for index, row in false_rows.iterrows(): 
        location = row['Location_']
        if organims == "Vibrio": 
            lat, lon = get_lat_lon5(location)
        elif organims == "Listeria": 
            lat, lon = get_lat_lon4(location)
            
        df_location__c_merged.at[index, 'Latitude'] = lat
        df_location__c_merged.at[index, 'Longitude'] = lon
        
        city, county, region, state, country = reverse_geocode(lat, lon)
        df_location__c_merged.at[index, 'City'] = city
        df_location__c_merged.at[index, 'County'] = county
        df_location__c_merged.at[index, 'Region'] = region
        df_location__c_merged.at[index, 'State'] = state
        df_location__c_merged.at[index, 'Country'] = country
    
    
        iso = get_iso_from_Lat_Lon(lat,lon)
        df_location__c_merged.at[index, 'ISO'] = iso.upper()  
        
        iso_2 = get_iso(iso)
        df_location__c_merged.at[index, 'ISO_2'] = iso_2[1]
        
        continent = get_continent(iso)
        df_location__c_merged.at[index, 'Continent'] = continent
        
    print("")
    print("Updated DataFrame using the country name")

## Drop Location_2 columns, _merge and Comparaison column and rename datafame 
columns_to_drop = ['Location_2', '_merge', 'Comparison']
df_location__c_merged.drop(columns=columns_to_drop, inplace=True)
df_location__c = df_location__c_merged.copy()
print(df_location__c.columns)

## Re-Check the geolocator package. 
## Verify that the country in the Location column and the country in the Country column are the same
## Create a new dataframe to compare the country name between the Location and the Country column (add Location_2 column)

# Create Location_2 column to keep only the countries names (delete elements after ,) in df_location__check dataframe 
df_location__c['Location_2'] = df_location__c['Location_'].astype(str).apply(lambda x: x.split(',')[0].strip()).to_frame()

# Create a mapping dictionary to handle 'USA' and 'United States' , 'KOREA' and 'SOUTH KOREA', 'VIET NAM' and  'VIETNAM', 'REUNION' and 'GUADELOUPE' to France ... 
mapping = {'USA': 'UNITED STATES',
           'PUERTO RICO' : 'UNITED STATES', 
           'KOREA' : 'SOUTH KOREA', 
           'VIET NAM': 'VIETNAM',
           'REUNION': 'FRANCE', 
           'GUADELOUPE' : 'FRANCE', 
           'HONG KONG': 'CHINA', 
           'CZECH REPUBLIC' : 'CZECHIA', 
           'GAMBIA': 'THE GAMBIA', 
           'MUMBAI' : 'INDIA',
           'FRENCH GUIANA': 'FRANCE'}
# Replace 'USA' with 'United States' in 'Location_2'
df_location__c['Location_2'] = df_location__c['Location_2'].replace(mapping)  
    
## Compare the Location_2 and the Country columns 
df_location__c['Country'] = df_location__c['Country'].str.upper()
df_location__c['Comparison'] = df_location__c['Location_2'] == df_location__c['Country']

## Check if all values in the 'Comparison' column are True
if df_location__c['Comparison'].all():
    print("The information in 'Location_2' column and 'Country' column are the same:", True)
else: 
    print("The information in 'Location_2' column and 'Country' column are the same:", False)
    # Print the row where the comparaison is False 
    print("")
    false_rows = df_location__c[~df_location__c['Comparison']]
    print("Different rows: ")
    print(false_rows[['Location', 'Location_','Location_2', 'Country']])
    print(false_rows)
    
print("")
# Count the occurrences of True and False values in the 'Comparison' column
count_true_false_ = df_location__c['Comparison'].value_counts().reset_index()
print("Number of occurences in the Comparaison column: ")
print(count_true_false_)


## Include location into the df using Lat/Lon columns 
## Include Latitude and Logitude column by the Lat/Lon column to add the Location information 
# Select rows with empty values in 'Location' but Lat/Lon column is not empty
df_location_emp = df[df['Location_'].isnull() & df['Lat/Lon_'].notnull()]

## Preposses the dartaframe before include the information in the Location column 
# Verify if the element from the Lat/Lon column starts with a digits
lat_list = df_location_emp['Lat/Lon_'].unique()

# Drop duplicates in 'Location', 'Lat/Lon', ''
df_location_emp_ = df_location_emp.drop_duplicates(subset=['Lat/Lon_'])

print("All element in the lat/lon columns with empty values in the Location column start with a digit:", all_elements_start_with_digit(lat_list))
print("Number of rows in the Lat/Lon column with empty values in the Location column:", len(df_location_emp))
print("Number of unique values after dropping duplicates in the Lat/Lon column:", len(df_location_emp_))
print("Dataframe shape:", df_location_emp_.shape)


## Add the Latitude and Logitude columns 
## Parse the coordinates from the Lat/Lon column using the parse_coordinates function 
df_location_emp_ = df_location_emp_.copy() 

if len(df_location_emp) != 0:
    coordinates_df = df_location_emp_['Lat/Lon_'].apply(parse_coordinates).apply(pd.Series)
    df_location_emp_[['Latitude', 'Longitude']] = coordinates_df

    ## Verify if the Latitude and Longirude columns had been included in the dataframe 
    print("Verify the new shape of the dataframe:", df_location_emp_.shape)


    ## check that the "Latitude" and "Longitude" columns are not empty 
    if df_location_emp_['Latitude'].notnull().all() and df_location_emp_['Longitude'].notnull().all():
        print("The Latitude and Longitude columns are not empty:", True)
    else: 
        print("Latitude and/or Longitude columns contain empty values")
        
    ## Include 'City', 'County', 'Region', 'State', 'Country' 'ISO, and 'ISO_2' columns to the dataframe 
    ## Use of the Latitude and Longitude columns and the function reverse_geocode to find the City, County, Region, State, Country name 
    ## Use the function get_iso to include the columns ISO (two letters code) and ISO_2 (three letters code) to the dataframe 
    ## Use of the Country column to get the Continent name using the function get_continent 

    df_location_emp_ = df_location_emp_.copy()

    result_columns = ['City', 'County', 'Region', 'State', 'Country']
    df_location_emp_[result_columns] = df_location_emp_.apply(lambda row: pd.Series(reverse_geocode(row['Latitude'], row['Longitude'])), axis=1)
    result_columns_iso = ['ISO', 'ISO_2']
    df_location_emp_[result_columns_iso] = df_location_emp_['Country'].apply(get_iso).apply(pd.Series)
    df_location_emp_.loc[:,'Continent'] = df_location_emp_.loc[:,'Country'].apply(get_continent)
    print(df_location_emp_.head())

    ## Verify if the Latitude and Longirude columns had been included in the dataframe 
    print("Verify the new shape of the dataframe:", df_location_emp_.shape)

    df_location_emp_['Country'] = df_location_emp_['Country'].str.upper()

    # Check that the "Country", "ISO", and "ISO_2" columns are not empty
    if df_location_emp_['Country'].notnull().all() and df_location_emp_['ISO'].notnull().all() and df_location_emp_['ISO_2'].notnull().all():
        print("The Country and ISO columns are not empty:", True)
    else: 
        print("Country and/or ISO columns contain empty values")
        
        # Drop rows where both 'Location' and 'Country' are NaN (No information found using Lat/Lon column)
        print("Location and Country information were not available: Drop rows ")
        df_location__concat = df_location_emp_.dropna(subset=['Location_', 'Country'], how='all')


    ## Concatenate the df_location_emp_ dataframe and the df_location__ dataframe 
    df_location__concat = pd.concat([df_location__c, df_location_emp_], ignore_index=True)
    
else: 
    df_location_emp_ = pd.DataFrame() 
    df_location__concat = df_location__c.copy()
    

columns_to_keep = ['Location',
       'Lat/Lon', 'Latitude', 'Longitude', 'City', 'County',
       'Region', 'State', 'Country', 'ISO', 'ISO_2', 'Continent']

location_df = df_location__concat[columns_to_keep]
location_df = location_df.drop_duplicates().fillna('')
print("Number of rows in the location dataframe:" ,len(location_df))
print(location_df.head())

"""
Regroupe tous les termes dedescriptions dans un même champs
Puis pivot table pour voir si il y a des duplicats: c'est à dire des lignes pourlesquelles, si on supprime IFSAC les annotations 
foodex2 deviennent différentes pour une même description de source
"""
foodex2['description'] = foodex2[['Isolation type', 'Isolation source','Host', 'Source type', 'Host disease', 'Ambiguity']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
foodex2['count'] = [1]*len(foodex2)

foodex2_pivot = foodex2[['description', 'Code', 'count']].pivot_table(index='description', columns='Code', values='count').fillna(0)
foodex2_pivot_ = foodex2_pivot.sum(axis = 1).to_frame()
foodex2_pivot_[foodex2_pivot_[0]>1]

## Foodex2 
# Replace empty strings with NaN
df_source = df.copy().fillna('')

df_source.loc[:, 'Isolation source'] = df_source.loc[:, 'Isolation source'].apply(lambda x: ' '.join(simple_preprocess(x)))
df_source.loc[:, 'Isolation type'] = df_source.loc[:, 'Isolation type'].apply(lambda x: ' '.join(simple_preprocess(x)))
df_source.loc[:, 'Host'] = df_source.loc[:, 'Host'].apply(lambda x: ' '.join(simple_preprocess(x)))
df_source.loc[:, 'Host disease'] = df_source.loc[:, 'Host disease'].apply(lambda x: ' '.join(simple_preprocess(x)))
df_source.loc[:, 'Source type'] = df_source.loc[:, 'Source type'].apply(lambda x: ' '.join(simple_preprocess(x)))

#remove IFSAC Codification
foodex2__ = foodex2__.copy()
foodex2_reduce = foodex2__[['Isolation type', 'Isolation source', 'Host', 'Source type',
       'Host disease', 'Code', 'Foodex2']].drop_duplicates()
df_foodex2 = pd.merge(df_source, foodex2_reduce, on=['Isolation type', 'Isolation source','Host', 'Source type', 'Host disease'], how='left', indicator='_foodex')
print("Merged DataFrame has the same number of rows as the original: %s"%(len(df_foodex2) == len(df)))

df_foodex2 = df_foodex2[['Assembly', 'Isolation type', 'Isolation source', 'Host', 'Source type',
       'Host disease', 'Location', 'Lat/Lon', 'Collection date', 'Code', 'Foodex2', '_foodex']]

print(df_foodex2.tail(2))
print(df_foodex2.shape)

df_foodex2['Foodex2'] = df_foodex2['Foodex2'].replace('NA', np.nan).replace('?', np.nan).replace('ROOT', np.nan)
df_foodex2['Code'] = df_foodex2['Code'].replace('NA', np.nan).replace('?', np.nan).replace('0ROOT', np.nan)
print('')

df_foodex2 = df_foodex2.astype(str).fillna('')
merge_stats = df_foodex2['_foodex'].value_counts().to_frame(name='count')
print("Check merge between the database and FoodEx2 annotations: %s"%(merge_stats['count'].sum() == len(df)))
print('')
print("Details of the merge between the database and FoodEx2 annotations")
print(merge_stats)

df_foodex2 = df_foodex2[df_foodex2['_foodex'] == "both"]
print('')

df_foodex2 = df_foodex2[df_foodex2['Foodex2'] != "nan"]
print(df_foodex2.shape)

## convert to lowercase the values of Location column 
df_foodex2_ = df_foodex2.copy()
df_foodex2_.loc[:,'Location'] = df_foodex2_.loc[:,'Location'].str.lower()
df_foodex2_.loc[:,'Lat/Lon'] = df_foodex2_.loc[:,'Lat/Lon'].str.lower()

location_df.loc[:,'Location'] = location_df.loc[:,'Location'].str.lower()
location_df.loc[:,'Lat/Lon'] = location_df.loc[:,'Lat/Lon'].str.lower()

## Rename NaN values 
NaN_val = {
    'nan' : 'NaN', 
    'not collected' : 'NaN', 
    'not available: to be reported later' : 'NaN', 
    '-' : 'NaN', 
    'other' : 'NaN', 
    '' : 'NaN', 
    'none' : 'NaN', 
    'not available' : 'NaN',
    'not provided' : 'NaN',
    ' ' :  'NaN'
}
df_foodex2_.loc[:,'Location'] = df_foodex2_.loc[:,'Location'].astype(str).replace(NaN_val)
df_foodex2_.loc[:,'Location'] = df_foodex2_.loc[:,'Location'].replace('NaN', np.nan).fillna('')
location_df.loc[:,'Location'] = location_df.loc[:,'Location'].astype(str).replace(NaN_val)
location_df.loc[:,'Location'] = location_df.loc[:,'Location'].replace('NaN', np.nan).fillna('')

df_foodex2_.loc[:,'Lat/Lon'] = df_foodex2_.loc[:,'Lat/Lon'].astype(str).replace(NaN_val)
df_foodex2_.loc[:,'Lat/Lon'] = df_foodex2_.loc[:,'Lat/Lon'].replace('NaN', np.nan).fillna('')
location_df.loc[:,'Lat/Lon'] = location_df.loc[:,'Lat/Lon'].astype(str).replace(NaN_val)
location_df.loc[:,'Lat/Lon'] = location_df.loc[:,'Lat/Lon'].replace('NaN', np.nan).fillna('')

## Apply homogenization to the 'Location' column and convert the values of the Location column to uppercase
df_foodex2_.loc[:,'Location'] = df_foodex2_.loc[:,'Location'].apply(split_colon).str.upper()
location_df.loc[:,'Location'] = location_df.loc[:,'Location'].apply(split_colon).str.upper()

# for col in ['Location', 'Lat/Lon']:
#     location_df[col] = location_df[col].astype(str).str.strip().str.lower()
#     df_foodex2[col] = df_foodex2[col].astype(str).str.strip().str.lower()
    
columnas_clave = ['Location', 'Lat/Lon']

df_merged = pd.merge(df_foodex2_, location_df, on=columnas_clave, how='left', indicator='_location')
df_merged = df_merged.fillna('')
merge_stats = df_merged['_location'].value_counts().to_frame()
print('Détail du merge de la base de données et des annotations foodex2')
print(merge_stats)

### Parse data df 
df_date = df_merged.copy()
# df_date = df_date[['Assembly', 'Collection date']]

df_date.loc[:,'Collection date'] = df_date.loc[:,'Collection date'].astype(str).replace(NaN_val)
df_date.loc[:,'Collection date'] = df_date.loc[:,'Collection date'].replace('NaN', np.nan)#.fillna('')
print("Total of NaN values in Collection date column:", df_date['Collection date'].isnull().sum())

df_date['Date'] = df_date['Collection date'].apply(try_parsing_date)

# Crea la columna "year" a partir de la columna de fecha de colección
# df_foodex2_location['Year'] = df_foodex2_location['Date'].apply(lambda x: x.year if not pd.isnull(x) else np.nan)
df_date['Year'] = df_date['Date'].dt.year

# Verificamos si los valores en la columna "Collection date" son numericos y tienen 4 dígitos
# 4 digitos = only year disponible en la columna 
# only_year column = True if only the year is available in the original dataframe, else false 
df_date['only_year'] = df_date['Collection date'].astype(str).str.isdigit() & (df_date['Collection date'].astype(str).str.len() == 4)

print(df_date.shape)
df_date.tail()

final_df = df_date.copy().fillna('')
final_df.replace('', np.nan, inplace=True)
print("Number of isolates with information in the Foodex2 column:", final_df.Foodex2.value_counts().sum())
print("Number of isolates with information in the Country column:", final_df.Country.value_counts().sum())
print("Number of isolates with information in the Location column:", final_df.Date.value_counts().sum())

df= pd.read_csv(df_path, sep='\t', low_memory=False)
columns_to_keep = ["Assembly", # Assembly ID 
                   'Isolation type', 'Isolation source', 'Host', 'Source type', 'Host disease', # Source 
                   'Location', 'Lat/Lon', # Location 
                   'Collection date'] # Date 
df = df[columns_to_keep].copy()
print(f"Number of {organims} isolates available on NCBI Pathogens (March 2024): %s"%(len(df)))
print('')

final_df_ = final_df.drop(columns=['Isolation type', 'Isolation source', 'Host', 'Source type', 'Host disease',
                   'Location', 'Lat/Lon', 
                   'Collection date', '_foodex', '_location'])

final_df__ = pd.merge(df,  final_df_, on="Assembly", indicator=True)
final_df__ = final_df__.fillna('')
merge_stats = final_df__['_merge'].value_counts().to_frame()
print('Détail du merge de la base de données et des annotations foodex2')
print(merge_stats)

final_df__ = final_df__.drop(columns="_merge")

final_df__ = final_df__.copy().fillna('')
final_df__.replace('', np.nan, inplace=True)
print("Number of isolates with information in the Foodex2 column:", final_df__.Foodex2.value_counts().sum())
print("Number of isolates with information in the Country column:", final_df__.Country.value_counts().sum())
print("Number of isolates with information in the Location column:", final_df__.Date.value_counts().sum())
print(final_df__.head())

## Only column assembly 
if organims == "Vibrio": 
    name = "VP"
if organims == "Listeria": 
    name = "LM"
    
final_df__[['Assembly']].to_csv(f'assemblies_values_standardization_{name}.csv', index=False)
## All dataframe 
final_df__.to_csv(f'{name}_standardization.csv', index=False)