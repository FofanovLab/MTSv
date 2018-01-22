# Setup Environment on Monsoon
### Load the Python3 version of Anaconda
```
$ module load anaconda/3.latest
```

### Create Conda environment with Python3 and required packages.
The environment only needs to be created once.
```
$ conda create --name biopy3 python=3.5.2 biopython pandas ete3
```
The name `biopy3` can be changed to whatever name you want.  

### Activate Conda Environment
```
$ source activate biopy3
```

### Deactivate Conda Environment
```
$ source deactivate
```

If you already have a Python3.5.2 environment you can just use conda to install the required packages.
```
$ source activate ENV_NAME
$ conda install -c etetoolkit ete3 ete_toolchain
$ conda install pandas
$ conda install biopython
```

# MTSv Summary
The `MTSv_summary.py` script summarizes the number of hits per taxon. The total number of reads mapped, the number of unique mapped reads, and the number of signature hits per taxon per sample (samples are in the same order as the FASTQ files that were passed to `MTSv-readprep`).

### Input
**Required Positional Arguments**  
`project_name`: Provide a prefix that will be used to name all output files.  
`collapse_file`: Path to MTSv-collapse output file  
**Optional Arguments**  
`--out_path`: Directory to write output (Default: ./)  

### Output

An Excel file <`out_path`/`program_name`_summary.xlsx> in the following format:  


| TaxID | Division | Sci. Name          | Sample | Total Hits | Unique Hits | Signature Hits |
|-------|----------|--------------------|--------|------------|-------------|----------------|
| 1392  | Bacteria | Bacillus anthracis | 0      | 12         | 2           | 2              |
| 1392  | Bacteria | Bacillus anthracis | 1      | 13         | 1           | 0              |
| 1392  | Bacteria | Bacillus anthracis | 2      | 15         | 2           | 10             |
| 1396  | Bacteria | Bacillus cereus    | 0      | 10         | 1           | 0              |
| 1396  | Bacteria | Bacillus cereus    | 1      | 13         | 1           | 0              |
| 1396  | Bacteria | Bacillus cereus    | 2      | 15         | 1           | 0              |

These results correspond to the following input file:
```
R1_10_13_15:1392,1396
R2_2_0_10:1392
```

### Usage
```
python MTSv_summary.py --help
usage: MTSv Summary [-h] [-o OUT_PATH] PROJECT_NAME COLLAPSE_FILE

Summarize number of hits for each taxa, including signature hits.

positional arguments:
  PROJECT_NAME          Project name and output file prefix
  COLLAPSE_FILE         Path to MTSv-collapse output file

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_PATH, --out_path OUT_PATH
                        Output directory (default: ./)
```
### Example Slurm Script
Change nauid to your nauid and modify `out_path` to test and run script  

```
#!/bin/bash
#SBATCH --job-name=MTSv-summary
#SBATCH --output=/scratch/nauid/output.txt 

module load anaconda/3.latest
source activate biopy3

python -u MTSv_summary.py test \
/scratch/tf362/vedro/merge/merged_results.txt \
--out_path /scratch/nauid/path/to/output/ \
```

# MTSv Extract
The `MTSv_extract.py` script extracts all unique read sequences that aligned to a provided taxid or species name.

### Input
**Required Positional Arguments**  
`project_name`: provide prefix that will be used to name all output files  
`collapse_file`: Path to MTSv-collapse output file  
`signature_file`: Path to MTSv-inform output file  
`reads_fasta`: Path to FASTA file produced from MTSv-readprep  

**Required Mutually Exclusive Arguments**  
`taxid`: the taxid to extract  
`species`: the species name to extract

**Optional Arguments**  
`--out_path`: Directory to write output (Default: ./)

### Output
`PROJECT_NAME_TAXID_all.fasta`: Contains all sequence reads that aligned to taxid.  
`PROJECT_NAME_TAXID_signature.fasta`: Contains sequence reads that only aligned to this taxid and no other taxids.

### Usage
```
python MTSv_extract.py --help
usage: MTSv Extract [-h] [-o OUT_PATH] (-t TAXID | -s SPECIES)
                    PROJECT_NAME COLLAPSE_FILE SIGNATURE_FILE READS_FASTA

Extracts all sequences, including signature hits, that aligned to a given
taxid or species name.

positional arguments:
  PROJECT_NAME          Project name and output file prefix
  COLLAPSE_FILE         Path to MTSv-collapse output file
  SIGNATURE_FILE        Path to MTSv-inform output file
  READS_FASTA           Path to FASTA file from MTSv-readprep

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_PATH, --out_path OUT_PATH
                        Output directory (default: ./)
  -t TAXID, --taxid TAXID
                        Extract sequences by taxid (default: None)
  -s SPECIES, --species SPECIES
                        Extract sequences by species name (default: None)

```
### Example Slurm Script
Change nauid to your nauid and modify `out_path` to test and run script
```
#!/bin/bash
#SBATCH --job-name=MTSv-extract
#SBATCH --output=/scratch/nauid/output.txt      
module load anaconda/3.latest
source activate biopy3

python -u MTSv_extract.py test \
/scratch/tf362/vedro/merge/merged_results.txt \
/scratch/tf362/vedro/inform/informative.txt \
/scratch/tf362/vedro/readprep/seg50_minqual15-qualthresh-3.fasta \
--out_path /scratch/nauid/path/to/output/ \
--taxid 9606 # human
```
