### Scripts

Python scripts used in the data processing pipeline:

- `download_genomes.py`: Downloads genome assemblies from NCBI repository. 
- `deduplication.py`: Identifies and removes duplicate genome records.
- `standardization.py`: Standardizes metadata: location, source, and date columns creation. 

#### `env/`

Contains the conda environment specification:

- `genbak_request.yaml`: conda environment file for download_genomes.py
- `standardization.yaml`: conda environment file for standardization.py