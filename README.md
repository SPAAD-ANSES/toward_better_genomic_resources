# Toward better genomic resources: Data and Scripts

This repository contains the data, scripts, and environments specifications used in the study _"Toward better public genomic resources:  two standardized datasets for L. monocytogenes and V. parahaemolyticus"_.

## Repository structure 
```text
.
├── data/
│ ├── Listeria_monocytogenes_data.csv
│ ├── Listeria_monocytogenes_data.json
│ ├── Vibrio_parahaemolyticus_data.csv
│ └── Vibrio_parahaemolyticus_data.json
├── scripts/
│ ├── download_genomes.py
│ ├── deduplication.py
│ └── standardization.py
│ ├──── env/
|   ├── genbak_request.yaml
|   └── standardization.yaml
├── tables/
│ ├── Listeria_to_foodex2.csv
│ ├── Vibrio_to_foodex.csv
│ ├── foodex2.txt
│ └── new_foodex2.txt
├── figures/
│ ├── Listeria/
│ └── Vibrio/

```

### `data/`
Cointains the final processed datasets for 
- ***Listeria monocytogenes***: In CSV and JSON formats 
- ***Vibrio parahaemolyticus***: In CSV and JSON formats 

These datasets include 
- genome metadata: 
    - Raw metadata columns used for standardization.
    - standardized source, location, and date descriptions. 
- Quality metrics of the genomes. 
- Sequence Type information of each genome.  

### `scripts/`

Python scripts used in the data processing pipeline:

- `download_genomes.py`: Downloads genome assemblies from NCBI repository. 
- `deduplication.py`: Identifies and removes duplicate genome records.
- `standardization.py`: Standardizes metadata: location, source, and date columns creation. 

#### `scripts/env/`

Contains the conda environment specification:

- `genbak_request.yaml`: conda environment file for download_genomes.py
- `standardization.yaml`: conda environment file for standardization.py

### `tables/`

Mapping and ontology tables for source standardization:
- Listeria_to_foodex2.csv: is the correspondence table between text-free description of the source from the NCBI portal and the foodEx2 ontology for ***Listeria monocytogenes***.
- Vibrio_to_foodex.csv: correspondence table between text-free description of the source from the NCBI portal and the foodEx2 ontology for ***Vibrio parahaemolyticus***.
- foodex2.txt: contains complete foodEx2 ontology used for source strandardization.
- new_foodex2.csv: Extended/custom ontology terms proposed during this work.

### `figures/`
Contains all figures related to: 
- Quality control (e.g., Total lenght, N50, GC%, etc).
- Metadata exploration (e.g., metadata distribution, source distribution, etc.).

In a HTML or PNG  format for direct use or viasualization.

## Genome Quality analysis 

The following bioinformatics tools were used to assess genome quality and completeness. 

### Quast (v5.2.0)

Quast was used to compute genome assembly quality metrics: 

```bash
quast -o quast_output genome_assembly.fasta
```

### BUSCOs (v5.7.0)
BUSCOs was used to assess genome completeness using the followig databases: 
- **Listeria monocytogenes**: bacillales_odb10
- **Vibrio parahaemolyticus**: vibrionales_odb10

```bash
busco -i genome_assembly.fasta -l bacillales_odb10 -o busco_output -m genome -f --download_path ./busco_db --offline \
&& cp busco_output/short_summary*.json busco_summary.json
```

### Kraken2 (v2.1.2)

Kraken2 was used to check taxonomy identity of each assembly (detect possible misclassified samples): 

```bash
kraken2 --db /path/to/minikraken2_v2_8GB_201904_UPDATE --report /path/to/outptu_report --output /path/to/output /path/to/assembly
```

### Bakta (v1.9.3) - CDS counting 
In this study, Bakta was used  specifically to extract the number of Coding DNA Sequence (CDS) per genome.

- Database: BAKTA-DB-v5.1

```bash
bakta genome_assembly.fasta \
      --db /path/to/BAKTA-DB-v5.1 \
      --output bakta_output \
      --threads 4 \
      --force
```

## Genome Typing analysis 

The following bioinformatics tool was used to assess genome typing.

### MSLT (v.2.23.0)

To obtain the MLST (Multi-Locus Sequence Typing) profile for each genome, the following command was used:

```bash
mlst --threads 4 --nopath --novel novel_alleles.txt --label GCA_000315135.1 --blastdb /path/to/blastdb genome_assembly.fasta > mlst_result.txt
```

### cgMLST (chewbbaca v3.0)

To obtain the cgMLST (core-genome Multi-Locus Sequence Typing) profile for each genome, the following command was used:

```bash
chewBBACA.py AlleleCall -i /path/to/assemblies/  -g /path/to/schema/ -o /path/to/results/ --minimum-length 144 --cpu 30 --gl /path/to/schema/listGenes.txt
```

cgMLST schema were obtained from:
- Moura, A. et al. 2016 for *Listeria moncytogenes*
- Gonzalez-Escalona, N. et al. 2017 for *Vibrio parahaemolyticus*

These schemas were provided in FASTA fromat and prepared using the following command:

```bash
chewBBACA.py PrepExternalSchema -g /path/to/fasta_files_of_schema/ -o /path/to/schema/ 
```
## Genomic distance analysis

### Mash (v2.3)
Pairwise Mash distances between genomes assemblies were computed by first generating sketches of all assemblies, followed by distance calculation based on these sketches results:

```bash
mash sketch -o /path/to/sketch_output/ /path/to/assembly/
```
```bash
mash dist -t path/to/sketch_output/target path/to/sketch_output/queries > /path/to/dist_output/
```